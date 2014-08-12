import re

cdef bytes invalid_cmd = "<+*-_Invalid Command_-*+>".encode('UTF-8')
cdef bytes composed_protected_cmd = "<+*-_Composed Protected Command_-*+>".encode('UTF-8')
cdef bytes composed_unprotected_cmd = "<+*-_Composed Unprotected Command_-*+>".encode('UTF-8')

cdef class SLIDataContainer:
    cdef classes.NESTEngine *nest_engine
    cdef commands

    def __cinit__(self):
        self.commands = {}

    cdef initialize(self, classes.NESTEngine *nest):
        self.nest_engine = nest

    cdef PyToken generate_func_pytoken(self, string cmd):
        cdef PyToken t = PyToken()

        self.nest_engine.run('/'.encode("UTF-8") + cmd + ' lookup'.encode("UTF-8"))

        result = self.nest_engine.pop()

        if result is True:
            t.thisptr = self.nest_engine.pop_token()

        return t

    cdef PyToken generate_arg_pytoken(self, string cmd):
        cdef PyToken t = PyToken()

        self.nest_engine.run("/ajvehwlksjdbjds ".encode("UTF-8") + cmd + " def".encode("UTF-8"))
        self.nest_engine.run("/ajvehwlksjdbjds lookup".encode("UTF-8"))

        if self.nest_engine.pop():
            t.thisptr= self.nest_engine.pop_token()

        return t
        

    cdef bint is_command(self, string cmd):
        for key in self.commands.keys():
            if cmd.compare(key) == 0:
                return True

        return False

    cdef bint add_command(self, string cmd):
        cdef PyToken token

        if not self.is_command(cmd):
            token = self.generate_func_pytoken(cmd)

            if token.thisptr:
                self.commands[cmd] = token
                return True
            else:
                return False
        else:
            return True


    cdef PyToken get_pytoken(self, string cmd):
        if self.is_command(cmd):
            return self.commands[cmd]
        else:
            return None

    cdef run(self, string cmd):
        cdef bytes command
        cdef PyToken t
        cdef composed_cmd

        m = re.search('^{ (.+?) } runprotected$'.encode("UTF-8"), cmd)
        if hasattr(m, 'group'):
            command = m.group(1)
            composed_cmd = composed_protected_cmd
        else:
            command = cmd
            composed_cmd = composed_unprotected_cmd

        if re.match('^[^ /]+$'.encode("UTF-8"), command):
            if self.add_command(command):
                t = self.commands[command]
                return self.nest_engine.run_token(t.thisptr[0])
            else:
                return invalid_cmd
        else:
            return composed_cmd


