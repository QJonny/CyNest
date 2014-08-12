include "/home/jonny/Programs/Nest/include/Neuron.pyx"

import math
import sys

cdef class cython_iaf_psc_delta_c_members(Neuron):
    cdef double tau_m
    cdef double C_m
    cdef double t_ref
    cdef double E_L
    cdef double I_e_
    cdef double V_th
    cdef double V_min_
    cdef double V_reset
    cdef bint with_refr_input_
    cdef double y0_
    cdef double y3_
    cdef double r_
    cdef double refr_spikes_buffer_
    cdef double P30_
    cdef double P33_
    cdef int RefractoryCounts_
    cdef double ms_resolution

    def __cinit__(self):
        self.tau_m   = 10.0  # ms
        self.C_m     = 250.0 # pF
        self.t_ref   =  2.0 # ms
        self.E_L     = -70.0 # mV
        self.I_e_     = 0.0 # pA
        self.V_th    = -55.0 - self.E_L # mV, rel to U0_
        self.V_min_   = -sys.float_info.max # relative U0_-55.0-U0_),  // mV, rel to U0_
        self.V_reset = -70.0 - self.E_L
        self.with_refr_input_ = False
        self.y0_      = 0.0
        self.y3_      = 0.0
        self.r_       = 0
        self.refr_spikes_buffer_ = 0.0
        self.P30_ = 0.0
        self.P33_ = 0.0
        self.RefractoryCounts_ = 0

    cpdef calibrate(self):
        self.ms_resolution = Time_get_resolution().get_ms()
        self.P33_ = math.exp(-self.ms_resolution/self.tau_m)
        self.P30_ = 1/self.C_m*(1-self.P33_)*self.tau_m
        self.RefractoryCounts_ = Time(ms(self.t_ref)).get_steps()


    cpdef update(self):
        if self.r_ == 0:
            # neuron not refractory
            self.y3_ = self.P30_*(self.y0_ + self.I_e_) + self.P33_*self.y3_ + (self.ex_spikes + self.in_spikes)
        
            # if we have accumulated spikes from refractory period, 
            # add and reset accumulator
            if self.with_refr_input_ and self.refr_spikes_buffer_ != 0.0:
                self.y3_ += self.refr_spikes_buffer_
                self.refr_spikes_buffer_ = 0.0
      
            # lower bound of membrane potential
            if self.y3_ < self.V_min_:
                self.y3_ = self.V_min_	 
        
        else: # neuron is absolute refractory
            # read spikes from buffer and accumulate them, discounting
            # for decay until end of refractory period
            if self.with_refr_input_:
                self.refr_spikes_buffer_ += (self.ex_spikes + self.in_spikes) * math.exp(-self.r_ * self.ms_resolution / self.tau_m)
            else:
                self.ex_spikes = 0  # clear buffer entry, ignore spike
                self.in_spikes = 0

            self.r_ -= 1
        
        self.y0_ = self.currents
        
        # threshold crossing
        if self.y3_ >= self.V_th:
            self.r_ = self.RefractoryCounts_
            self.y3_ = self.V_reset
            self.spike = 1 # True
        else:
            self.spike = 0 # False


    cpdef getStatus(self):
        cdef dict d = {}
        d["tau_m"] = self.tau_m
        d["C_m"] = self.C_m
        d["t_ref"] = self.t_ref
        d["E_L"] = self.E_L
        d["I_e_"] = self.I_e_
        d["V_th"] = self.V_th
        d["V_min_"] = self.V_min_
        d["V_reset"] = self.V_reset
        d["with_refr_input_"] = self.with_refr_input_
        d["y0_"] = self.y0_
        d["y3_"] = self.y3_
        d["r_"] = self.r_
        d["refr_spikes_buffer_"] = self.refr_spikes_buffer_
        d["P30_"] = self.P30_
        d["P33_"] = self.P33_
        d["RefractoryCounts_"] = self.RefractoryCounts_
        return d

    cpdef setStatus(self, d):
        if self.tau_m in d:
            self.tau_m = d["tau_m"]
        if self.C_m in d:
            self.C_m = d["C_m"]
        if self.t_ref in d:
            self.t_ref = d["t_ref"]
        if self.E_L in d:
            self.E_L = d["E_L"]
        if self.I_e_ in d:
            self.I_e_ = d["I_e_"]
        if self.V_th in d:
            self.V_th = d["V_th"]
        if self.V_min_ in d:
            self.V_min_ = d["V_min_"]
        if self.V_reset in d:
            self.V_reset = d["V_reset"]
        if self.with_refr_input_ in d:
            self.with_refr_input_ = d["with_refr_input_"]
        if self.y0_ in d:
            self.y0_ = d["y0_"]
        if self.y3_ in d:
            self.y3_ = d["y3_"]
        if self.r_ in d:
            self.r_ = d["r_"]
        if self.refr_spikes_buffer_ in d:
            self.refr_spikes_buffer_ = d["refr_spikes_buffer_"]
        if self.P30_ in d:
            self.P30_ = d["P30_"]
        if self.P33_ in d:
            self.P33_ = d["P33_"]
        if self.RefractoryCounts_ in d:
            self.RefractoryCounts_ = d["RefractoryCounts_"]

