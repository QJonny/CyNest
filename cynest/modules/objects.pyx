

# In order to avoid the problem of nested classes,
# not supported by Cython, a UnitManager class has
# been created, which can embed on the C++ side
# objects of ms, tic, step and ms_stamp type (classes.Time._ _ _).
# The code is the following:
# 1: tic, 2: step, 3: ms, 4: ms_stamp
# The usefullness is that a classes.Time objects can be created
# based on the internal state of UnitManager

cdef class Unit:
    cdef classes.UnitManager *thisptr

    def __dealloc__(self):
        del self.thisptr

cdef class tic(Unit):
    def __cinit__(self, t):
        self.thisptr= new classes.UnitManager(1, <long>t)

    def create(self, t):
        return tic(t)


cdef class step(Unit):
    def __cinit__(self, t):
        self.thisptr= new classes.UnitManager(2, <long>t)
       
    def create(self, t):
        return step(t)


cdef class ms(Unit):
    def __cinit__(self, t):
        self.thisptr= new classes.UnitManager(3, <double>t)

    def create(self, t):
        return ms(t)


cdef class ms_stamp(Unit):
    def __cinit__(self, t):
        self.thisptr= new classes.UnitManager(4, <double>t)

    def create(self, t):
        return ms_stamp(t)




cdef class Time:
    cdef classes.Time thisptr

    def __cinit__(self, Unit t):
        # This function generates a classes.Time object
        # based of the internal state of UnitManager,
        # which is in turn based on the real type of t
        self.thisptr = t.thisptr[0].generateTime()

    # Creates a Time object from a classes.Time one
    # (somehow ugly)
    cdef Time getTime(self, classes.Time t):
        cdef Time tm = Time(ms(0.0))
        tm.thisptr = t
        return tm

    # This method creates a new instance of the object.
    # It is actually a way for creating new instances
    # from an object instead of a class definition.
    # That's useful in the imported pyx files
    # since they have no access to the class definitions
    # but they can handle general objects
    def create(self, t):
        return Time(t)

    def get_tics(self):
        return self.thisptr.get_tics()

    def get_steps(self):
        return self.thisptr.get_steps()

    def get_ms(self):
        return self.thisptr.get_ms()

    def set_to_zero(self):
        self.thisptr.set_to_zero()

    def calibrate(self):
        self.thisptr.calibrate()

    def advance(self):
        self.thisptr.advance()

    def is_grid_time(self):
        return self.thisptr.is_grid_time()

    def is_neg_inf(self):
        return self.thisptr.is_neg_inf()

    def is_pos_inf(self):
        return self.thisptr.is_pos_inf()

    def is_finite(self):
        return self.thisptr.is_finite()

    def is_step(self):
        return self.thisptr.is_step()

    def succ(self):
        return self.getTime(self.thisptr.succ())

    def pred(self):
        return self.getTime(self.thisptr.pred())


# static methods

    def get_resolution(self):
        return self.getTime(classes.get_resolution())
        
    def set_resolution(self, d):
        classes.set_resolution(d)

    def reset_resolution(self):
        classes.reset_resolution()

    def resolution_is_default(self):
        return classes.resolution_is_default()

    def get_ms_per_tic(self):
        return classes.get_ms_per_tic()

    def get_tics_per_ms(self):
        return classes.get_tics_per_ms()

    def get_tics_per_step(self):
        return classes.get_tics_per_steps()

    def get_old_tics_per_step(self):
        return classes.get_old_tics_per_step()

    def get_tics_per_step_default(self):
        return classes.get_tics_per_step_default()

    def min(self):
        return self.getTime(classes.min())

    def max(self):
        return self.getTime(classes.max())

    def pos_inf(self):
        return self.getTime(classes.pos_inf())

    def neg_inf(self):
        return self.getTime(classes.neg_inf())



cdef class Scheduler:
    def get_modulo(self, v):
        return classes.get_modulo(v)

    def get_slice_modulo(self, v):
        return classes.get_slice_modulo(v)

    def get_min_delay(self):
        return classes.get_min_delay()

    def get_max_delay(self):
        return classes.get_max_delay()
