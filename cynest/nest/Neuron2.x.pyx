# cython: language_level=2

cdef class Neuron:
    # Standard Parameters
    cdef double currents
    cdef double in_spikes
    cdef double ex_spikes
    cdef long t_lag
    cdef long spike
    cdef double current_value
    
    def __cinit__(self):
        pass

    cpdef calibrate(self):
        pass
        
    cpdef update(self):
        pass
        
    cpdef getStatus(self):
        return {}
        
    cpdef setStatus(self, params):
        pass
       
    # We convert the address of the variables into long in order to extract the pointers

    cpdef getPCurrents(self):
        return <long>(&(self.currents))
        
    cpdef getPIn_Spikes(self):
        return <long>(&(self.in_spikes))
        
    cpdef getPEx_Spikes(self):
        return <long>(&(self.ex_spikes))
       
    cpdef getPT_Lag(self):
        return <long>(&(self.t_lag))

    cpdef getPSpike(self):
        return <long>(&(self.spike))

    cpdef getPCurrent_Value(self):
        return <long>(&(self.current_value))

# This class contains the totality of the imported objects
# needed for accessing classes on the project side
cdef class ObjectManager:
    cdef object schedulerObj
    cdef object timeObj
    cdef object ticObj
    cdef object stepObj
    cdef object msObj
    cdef object ms_stampObj

    def setScheduler(self, obj):
        self.schedulerObj = obj

    def setTime(self, obj):
        self.timeObj = obj

    def setTic(self, obj):
        self.ticObj = obj

    def setStep(self, obj):
        self.stepObj = obj

    def setMs(self, obj):
        self.msObj = obj

    def setMs_stamp(self, obj):
        self.ms_stampObj = obj


cdef ObjectManager objectManager = ObjectManager()


# These methods are called from the project in order to
# fill the different imported objects

def setScheduler(obj):
    objectManager.setScheduler(obj)

def setTime(obj):
    objectManager.setTime(obj)

def setTic(obj):
    objectManager.setTic(obj)

def setStep(obj):
    objectManager.setStep(obj)

def setMs(obj):
    objectManager.setMs(obj)

def setMs_stamp(obj):
    objectManager.setMs_stamp(obj)



# version function (returns 2 or 3, depending on the python release)
def getVersion():
   return 2


# When creating a new objects (of the ones below), their internal 
# real object must be cloned in order to create another instance.
# Thus we call the create method of the imported object

cdef class Unit:
    cdef object ob



cdef class tic(Unit):
    def __cinit__(self, t):
        self.ob = objectManager.ticObj.create(t)



cdef class step(Unit):
    def __cinit__(self, t):
        self.ob = objectManager.stepObj.create(t)



cdef class ms(Unit):
    def __cinit__(self, t):
        self.ob = objectManager.msObj.create(t)



cdef class ms_stamp(Unit):
    def __cinit__(self, t):
        self.ob = objectManager.ms_stampObj.create(t)



cdef class Time:
    cdef object time_ob

    def __cinit__(self, Unit t):
        self.time_ob = objectManager.timeObj.create(t.ob)

    def get_tics(self):
        return self.time_ob.get_tics()

    def get_steps(self):
        return self.time_ob.get_steps()

    def get_ms(self):
        return self.time_ob.get_ms()

    def set_to_zero(self):
        self.time_ob.set_to_zero()

    def calibrate(self):
        self.time_ob.calibrate()

    def advance(self):
        self.time_ob.advance()

    def is_grid_time(self):
        return self.time_ob.is_grid_time()

    def is_neg_inf(self):
        return self.time_ob.is_neg_inf()

    def is_pos_inf(self):
        return self.time_ob.is_pos_inf()

    def is_finite(self):
        return self.time_ob.is_finite()

    def is_step(self):
        return self.time_ob.is_step()

    def succ(self):
        return self.time_ob.succ()

    def pred(self):
        return self.time_ob.pred()


# These are static methods (on the project side)
# Note that the imported objects itself are accessed
# instead of copies of them. This proves the static type
# of these methods.

cpdef Time_get_resolution():
    return objectManager.timeObj.get_resolution()

cpdef Time_set_resolution(d):
    objectManager.timeObj.set_resolution(d)

cpdef Time_reset_resolution():
    objectManager.timeObj.reset_resolution()

cpdef Time_resolution_is_default():
    return objectManager.timeObj.resolution_is_default()

cpdef Time_get_ms_per_tic():
    return objectManager.timeObj.get_ms_per_tic()

cpdef Time_get_tics_per_ms():
    return objectManager.timeObj.get_tics_per_ms()

cpdef Time_get_tics_per_step():
    return objectManager.timeObj.get_tics_per_steps()

cpdef Time_get_old_tics_per_step():
    return objectManager.timeObj.get_old_tics_per_step()

cpdef Time_get_tics_per_step_default():
    return objectManager.timeObj.get_tics_per_step_default()

cpdef Time_min():
    return objectManager.timeObj.min()

cpdef Time_max():
    return objectManager.timeObj.max()

cpdef Time_pos_inf():
    return objectManager.timeObj.pos_inf()

cpdef Time_neg_inf():
    return objectManager.timeObj.neg_inf()




cpdef Scheduler_get_modulo(v):
    return objectManager.schedulerObj.get_modulo(v)

cpdef Scheduler_get_slice_modulo(v):
    return objectManager.schedulerObj.get_slice_modulo(v)

cpdef Scheduler_get_min_delay():
    return objectManager.schedulerObj.get_min_delay()

cpdef Scheduler_get_max_delay():
    return objectManager.schedulerObj.get_max_delay()

