# cython: language_level=2
# Cython wrapper for the pynest functions.
from libcpp.string cimport string
from libcpp.vector cimport vector
cimport cpython # we need this for the unicode UTF-8 conversion

import signal

class NESTError(Exception):
    def __init__(self, msg) :
        Exception.__init__(self, msg)

# This imports the C++ class wrappers
cimport classes

cdef public class PyToken[object PyToken, type PyTokenType]:
     """
     Python wrapper of SLI's Token class.
     """
     cdef classes.Token *thisptr

     def __dealloc__(self):
         if self.thisptr:
            del self.thisptr

include "datamanager.pyx"

sli_func = None
broadcast = None


def cynest_signal_handler(signal,frame):
    raise KeyboardInterrupt()


cdef class NESTEngine:
    cdef classes.NESTEngine *thisptr
    cdef SLIDataContainer sli_container
    cdef bint _protected

    def __cinit__(self):
        self.thisptr= new classes.NESTEngine()
        self.sli_container = SLIDataContainer()
        self.sli_container.initialize(self.thisptr)
        self._protected = True
        
    def __dealloc__(self):
        print ("CyNEST says good bye.")
        del self.thisptr
       
    def run_protected (self):
        if self._protected:
            return True
        else:
            return False

    def init(self, argv, modulepath):
        """
        Startup the NEST engine.
        argv: list of strings with NEST's command line arguments
        modulepath: the path to CyNEST's startup file
        The function will return true if startup succeded and
        false otherwise.
        This function is part of the low-level API.
        """
        if len(argv[0]) == 0:
            argv[0]='cynest'

        argv_bytes= [ str.encode('UTF-8') for str in argv]
        cdef bytes modulepath_bytes=modulepath.encode('UTF-8')
        result= self.thisptr.init(argv_bytes, modulepath_bytes)

        if result:
           signal.signal(signal.SIGINT, cynest_signal_handler)
        return result



    def add_command(self, cmd):
        self.sli_container.add_command(cmd.encode('UTF-8'))

    def get_pytoken(self, cmd):
        return self.sli_container.get_pytoken(cmd.encode('UTF-8'))

    def generate_arg_pytoken(self, cmd):
        return self.sli_container.generate_arg_pytoken(cmd.encode('UTF-8'))


    def register_cython_model(self, model):
        self.thisptr.register_cython_model(model.encode('UTF-8'))



    def push(self, value):
        """
        Push a Python object onto NEST's operand stack.
        value can be almost any Python object.
        Not pushable are dictionaries with non-string keys.
        This function is part of the low-level API.
        """

        return self.thisptr.push(value)
   
    def pop(self):
        """
        Pop the top value from NEST's operand stack and return it as
        Python object.
        If the stack is empty, a NESTError exception is raised.
        This function is part of the low-level API.
        """
        return self.thisptr.pop()

    def run(self, command):
        """
        Execute a SLI command string.
        SLI is the native language of NEST. This function takes a string which is then parsed and executed by
        NEST. The string may contain an arbitrary sequence of NEST commands.
        If errors occur during execution, they are reported back to Python as NESTError exections.
        This function is part of the low-level API.
        """
        cdef bytes command_bytes=command.encode('UTF-8')

        result = self.sli_container.run(command_bytes)
        self._protected = False

        if result is invalid_cmd:
            print (NESTError("Cannot generate PyToken for the following command: " + command + "\nThe command will be executed in standard mode."))
            result = self.thisptr.run(command_bytes)
            self._protected = True
        elif result is composed_protected_cmd:
            result = self.thisptr.run(command_bytes)
            self._protected = True
        elif result is composed_unprotected_cmd:
            result = self.thisptr.run(command_bytes)
            self._protected = False

        signal.signal(signal.SIGINT, cynest_signal_handler)

        return result

    cpdef run_pytoken(self, PyToken command):
        """
        Execute a PyDatum object.
        This function is part of the low-level API.
        """
        return self.thisptr.run_token(command.thisptr[0])

    def pop_pytoken(self):
         t= PyToken()
         t.thisptr= self.thisptr.pop_token()
         return t

    cpdef push_pytoken(self, PyToken obj):
          """
          Push a token to NEST's operand stack.
          This function is part of the low-level API.
          """
          if obj.thisptr:
            return self.thisptr.push_token(obj.thisptr[0])
          else:
            raise NESTError("Cannot push empty PyToken.")

    def push_connections(self, connectome):
        """
        Push a list of dictionaries with connection information to NEST. Each dictionary is converted to a NEST ConnectionDatum
        which encapsulates the five-tuple defining a connection.
        This function is part of the low-level API.
        """
        return self.thisptr.push_connections(connectome)
    
    def check_engine(self):
        """
        Returns true if the object is properly initialized and false otherwise.
        """
        return self.thisptr.check_engine()

    def convergent_connect(self, pre, post, weight, delay, model):
        self.add_command('ConvergentConnect')
        cdef PyToken cmd = self.get_pytoken('ConvergentConnect')
        cdef PyToken m = self.generate_arg_pytoken('/'+model)

        if weight == None and delay == None:
            for d in post :
                self.push(pre)
                self.push(d)
                self.push_pytoken(m)
                self.run_pytoken(cmd) # always run in unprotected mode

        elif weight != None and delay != None:
            weight = broadcast(weight, len(pre), (float,), "weight")
            if len(weight) != len(pre):
                raise NESTError("weight must be a float, or sequence of floats of length 1 or len(pre)")
            delay = broadcast(delay, len(pre), (float,), "delay")
            if len(delay) != len(pre):
                raise NESTError("delay must be a float, or sequence of floats of length 1 or len(pre)")
        
            for d in post:
                self.push(pre)
                self.push(d)
                self.push(weight)
                self.push(delay)
                self.push_pytoken(m)
                self.run_pytoken(cmd) # always run in unprotected mode

        else:
            raise NESTError("Both 'weight' and 'delay' have to be given.")

    def divergent_connect(self, pre, post, weight, delay, model):
        self.add_command('DivergentConnect')
        cdef PyToken cmd = self.get_pytoken('DivergentConnect')
        cdef PyToken m = self.generate_arg_pytoken('/' + model)

        if weight == None and delay == None:
            for s in pre :
                self.push(s)
                self.push(post)
                self.push_pytoken(m)
                self.run_pytoken(cmd) # always run in unprotected mode

        elif weight != None and delay != None:
            weight = broadcast(weight, len(post), (float,), "weight")
            if len(weight) != len(post):
                raise NESTError("weight must be a float, or sequence of floats of length 1 or len(post)")
            delay = broadcast(delay, len(post), (float,), "delay")
            if len(delay) != len(post):
                raise NESTError("delay must be a float, or sequence of floats of length 1 or len(post)")

            for s in pre :
                self.push(s)
                self.push(post)
                self.push(weight)
                self.push(delay)
                self.push_pytoken(m)
                self.run_pytoken(cmd) # always run in unprotected mode

        else:
            raise NESTError("Both 'weight' and 'delay' have to be given.")


    def data_connect1(self, list pre, list params, model):
        self.add_command('DataConnect_i_dict_s')
        cdef PyToken cmd1 = self.get_pytoken('DataConnect_i_dict_s')
        cdef PyToken m = self.generate_arg_pytoken('/' + model)

        for s,p in zip(pre,params):
            self.push(s)
            self.push(p)
            self.push_pytoken(m)
            self.run_pytoken(cmd1) # sure it is not run_protected



    def data_connect2(self, list pre, list params, model):
        self.add_command('Connect_i_D_i')
        cdef PyToken cmd2 = self.get_pytoken('Connect_i_D_i')
        self.run('synapsedict') #sure unprotected
        self.run('/'+ model+ ' get') 

        cdef int model_id = self.pop()
        cdef dict params_dict = {}
        cdef int i = 0
        cdef int j = 0
        cdef int length = min(len(pre), len(params))

        for j from 0 <= j < length by 1:
            for i from 0 <= i < len(params[j]['target']) by 1:
                for key in params[j]:
                    params_dict[key] = params[j][key][i]

                self.thisptr.push(pre[j])
                self.thisptr.push(params_dict)
                self.thisptr.push(model_id)
                self.thisptr.run_token(cmd2.thisptr[0])



    def random_convergent_connect(self, pre, post, n, weight, delay, model, options):
        # store current options, set desired options
        old_options = None
        error = False
        if options:
            old_options = sli_func('GetOptions', '/RandomConvergentConnect', litconv=True)
            del old_options['DefaultOptions'] # in the way when restoring
            sli_func('SetOptions', '/RandomConvergentConnect', options,
                 litconv=True)

        if weight == None and delay == None:
            sli_func(
                '/m Set /n Set /pre Set { pre exch n m RandomConvergentConnect } forall',
                post, pre, n, '/'+model, litconv=True)
    
        elif weight != None and delay != None:
            weight = broadcast(weight, n, (float,), "weight")
            if len(weight) != n:
                raise NESTError("weight must be a float, or sequence of floats of length 1 or n")
            delay = broadcast(delay, n, (float,), "delay")
            if len(delay) != n:
                raise NESTError("delay must be a float, or sequence of floats of length 1 or n")

            sli_func(
                '/m Set /d Set /w Set /n Set /pre Set { pre exch n w d m RandomConvergentConnect } forall',
                post, pre, n, weight, delay, '/'+model, litconv=True)
    
        else:
            error = True

        # restore old options
        if old_options:
            sli_func('SetOptions', '/RandomConvergentConnect', old_options, litconv=True)

        if error:
            raise NESTError("Both 'weight' and 'delay' have to be given.")


    def random_divergent_connect(self, pre, post, n, weight, delay, model, options):
        # store current options, set desired options
        old_options = None
        error = False
        if options:
            old_options = sli_func('GetOptions', '/RandomDivergentConnect', litconv=True)
            del old_options['DefaultOptions'] # in the way when restoring
            sli_func('SetOptions', '/RandomDivergentConnect', options, litconv=True)

        if weight == None and delay == None:
            sli_func(
                '/m Set /n Set /post Set { n post m RandomDivergentConnect } forall',
                pre, post, n, '/'+model, litconv=True)

        elif weight != None and delay != None:
            weight = broadcast(weight, n, (float,), "weight")
            if len(weight) != n:
                raise NESTError("weight must be a float, or sequence of floats of length 1 or n")
            delay = broadcast(delay, n, (float,), "delay")
            if len(delay) != n:
                raise NESTError("delay must be a float, or sequence of floats of length 1 or n")

            sli_func(
                '/m Set /d Set /w Set /n Set /post Set { n post w d m RandomDivergentConnect } forall',
                pre, post, n, weight, delay, '/'+model, litconv=True)

        else:
            error = True

        # restore old options
        if old_options:
            sli_func('SetOptions', '/RandomDivergentConnect', old_options, litconv=True)
    
        if error:
            raise NESTError("Both 'weight' and 'delay' have to be given.")
    


cdef class NameDatum:
     """
     Python wrapper of SLI's NameDatum class. 
     This class is used to store SLI commands and variables.
     """
     cdef classes.NameDatum *thisptr
     def __cinit__(self, name):
         cdef string name_b=name.encode('UTF-8')
         self.thisptr = new classes.NameDatum(name_b)

     def __dealloc__(self):
         del self.thisptr

     def as_PyToken(self):
         """
         Create a new Token with the NameDatum.
         """
         t=PyToken()
         cdef classes.NameDatum *name_ptr=new classes.NameDatum(self.thisptr[0])
         t.thisptr= new classes.Token(<classes.Datum *>name_ptr)
         return t

     def str(self):
         """
         Return a string representation of the NameDatum.
         """
         return self.thisptr.toString()


cdef public object Token_to_PyObject(classes.Token *arg):
     """
     Convert a Datum pointer to a Python object.
     This function is exposed to C/C++ and used by the DatumToPythonConverter to
     encapsulate arbitrary Tokens in PyToken objects.
     """
     dat=PyToken()
     dat.thisptr=arg
     return dat



# Additional modules
include "modules/objects.pyx"
include "modules/graphics/graphics.pyx"
