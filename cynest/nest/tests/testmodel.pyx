include "/home/jonny/Programs/Nest/include/Neuron.pyx"

import math

cdef class C:
    cdef double p1
    cdef double p2

    def __cinit__(self):
        self.p1 = 1.0
        self.p2 = 250.0


cdef class testmodel(Neuron):
    cdef C params

    def __cinit__(self):
        self.params = C()

    cpdef calibrate(self):
        pass


    cpdef update(self):
        self.params.p1 = math.exp(self.ex_spikes + self.in_spikes) + self.params.p1
        if self.params.p1 >= self.params.p2:
            self.spike = 1 # True
            self.params.p1 = 0.0
        else:
            self.spike = 0 # False

    cpdef getStatus(self):
        cdef dict d = {}

        d["p1"] = self.params.p1
        d["p2"] = self.params.p2
        return d

    cpdef setStatus(self, d):
        if d.has_key("p1"):
            self.params.p1 = d["p1"]
        if d.has_key("p2"):
            self.params.p2 = d["p2"]

