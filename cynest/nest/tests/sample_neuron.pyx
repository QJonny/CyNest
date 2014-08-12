include "/home/jonny/Programs/Nest/include/Neuron.pyx"



cdef class sample_neuron(Neuron):
    cdef double V_m
    cdef double V_th
    cdef double E_L
    cdef int param

    def __cinit__(self):
        self.V_m = 0.0
        self.V_th = 0.0
        self.E_L = 0.0
        self.param = 12

    cpdef calibrate(self):
        pass


    cpdef update(self):
        self.param = self.param + 1

    cpdef getStatus(self):
        cdef dict d = {}

        d["V_m"] = self.V_m
        d["V_th"] = self.V_th
        d["E_L"] = self.E_L
        d["param"] = self.param
        return d

    cpdef setStatus(self, d):
        if "V_m" in d:
            self.V_m = d["V_m"]
        if "V_th" in d:
            self.V_th = d["V_th"]
        if "E_L" in d:
            self.E_L = d["E_L"]
        if "param" in d:
            self.param = d["param"]

