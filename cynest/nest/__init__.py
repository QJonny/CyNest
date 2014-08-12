"""
Initializer of CyNEST.
"""

import sys, os, atexit

# The following is a workaround to make MPI-enabled NEST import
# properly. The basic problem is that the shared object pynestkernel
# dynamically opens other libraries that open other libraries...
try:
    try:
        import dl
    except:
        import DLFCN as dl
    sys.setdlopenflags(dl.RTLD_NOW|dl.RTLD_GLOBAL)
except:
    # this is a hack for Python 2.6 on Mac, where RTDL_NOW is nowhere to
    # be found. See ticket #397
    import ctypes
    sys.setdlopenflags(ctypes.RTLD_GLOBAL) 

import cynest.hl_api
import cynest.kernel as _kernel

hl_api.nest = _kernel

hl_api.schedulerObj = _kernel.Scheduler()
hl_api.timeObj = _kernel.Time(_kernel.ms(0.0))
hl_api.ticObj = _kernel.tic(0)
hl_api.stepObj = _kernel.step(0)
hl_api.msObj = _kernel.ms(0.0)
hl_api.ms_stampObj = _kernel.ms_stamp(0.0)

_kernel.engine= _kernel.NESTEngine()


sli_push = _kernel.engine.push
hl_api.sps = sli_push
sps = sli_push

sli_pop = _kernel.engine.pop
hl_api.spp = sli_pop
spp = sli_pop

hl_api.dtc1 = _kernel.engine.data_connect1
hl_api.dtc2 = _kernel.engine.data_connect2
hl_api.cvc = _kernel.engine.convergent_connect
hl_api.dvc = _kernel.engine.divergent_connect
hl_api.rcc = _kernel.engine.random_convergent_connect
hl_api.rdc = _kernel.engine.random_divergent_connect
hl_api.reg = _kernel.engine.register_cython_model
   
def sli_run(*args):
    raise NESTError("CyNEST is not initialized properly. Please call init() first.")

def sli_func(s, *args, **kwargs):
    """This function is a convenience function for executing the 
       sequence sli_push(args); sli_run(s); y=sli_pop(). It takes
       an arbitrary number of arguments and may have multiple
       return values. The number of return values is determined by
       the SLI function that was called.

       Keyword arguments:
       namespace - string: The sli code is executed in the given SLI namespace.
       litconv   - bool  : Convert string args beginning with / to literals.
       
       Examples:
         r,q = sli_func('dup rollu add',2,3)
         r   = sli_func('add',2,3)
         r   = sli_func('add pop',2,3)
         l   = sli_func('CreateLayer', {...}, namespace='topology')
         opt = sli_func('GetOptions', '/RandomConvergentConnect', litconv=True)
    """

    # check for namespace
    slifun = 'sli_func'  # version not converting to literals
    if 'namespace' in kwargs:
        s = kwargs['namespace'] + ' using ' + s + ' endusing'
    elif 'litconv' in kwargs:
        if kwargs['litconv']:
            slifun = 'sli_func_litconv'
    elif len(kwargs) > 0:
        hl_api.NESTError("'namespace' and 'litconv' are the only valid keyword arguments.")
    
    sli_push(args)       # push array of arguments on SLI stack
    sli_push(s)          # push command string
    sli_run(slifun)      # SLI support code to execute s on args
    r=sli_pop()          # return value is an array

    if len(r) == 1:        # 1 return value is no tuple
        return r[0]
 
    if len(r) != 0:   
       return tuple(r)   # convert array to tuple

kernel_sr = _kernel.engine.run
hl_api.sr = sli_run
sr = sli_run
hl_api.sli_func = sli_func
_kernel.sli_func = sli_func
_kernel.broadcast = hl_api.broadcast

initialized = False

def catching_sr(cmd):
    """
    Send a command string to the NEST kernel to be executed.
    catching_sr is a wrapper of the kernel_sr to raise errors as Python errors.
    """

    kernel_sr('{ '+cmd+' } runprotected')
    if _kernel.engine.run_protected() and not sli_pop():
        errorname = sli_pop()
        message = sli_pop()
        commandname = sli_pop()
        raise hl_api.NESTError(errorname + ' in ' + commandname + message)


def catch_errors(catchErrors = True):
    """Switch between the catching and non-catching versions or sr"""

    global sr, sli_run

    if catchErrors:
        sli_run = catching_sr
    else:
        sli_run = kernel_sr

    sr = sli_run
    hl_api.sr = sli_run


def init(argv) :
    """Initialize. argv is passed to the NEST kernel."""

    global initialized

    if initialized:
        raise hl_api.NESTError("NEST already initialized.")
        return

    quiet = False
    if argv.count("--quiet") :
        quiet = True
        argv.remove("--quiet")

    initialized |= _kernel.engine.init(argv, __path__[0])

    if initialized :

        if not quiet :
            kernel_sr("pywelcome")
        catch_errors(True)

        # Dirty hack to get models completion in iPython shell
        try:
            __IPYTHON__
        except NameError:
            pass
        else:
            try:
                import keyword
                keyword.kwlist += hl_api.Models()
            except ImportError:
                pass

def test ():
    """ Runs a battery of unit tests on PyNEST """
    import cynest.tests
    import unittest

    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(cynest.tests.suite())


# graphics module
_kernel.setGraphicsParameters(hl_api, os.path.dirname(os.path.realpath(__file__)))

# wrapper function around setGraphicsSimulator (since global Simulator has to be manually changed (no transitivity))
def SetGraphicsSimulator(s):
    global Simulate, hl_api
    _kernel.setGraphicsSimulator(s)
    Simulate = hl_api.Simulate
    
initNetworkVisualizer = _kernel.initNetworkVisualizer
# end of graphics module


if not 'DELAY_PYNEST_INIT' in os.environ:
    init(sys.argv)

    from cynest.hl_api import *
