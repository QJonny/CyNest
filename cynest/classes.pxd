# Cython wrapper for the pynest kernel functions.
# This is a cython 'header' file that is read by other modules.
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "psignal.h":
  cdef int SLIsignalflag
  cdef void SLISignalHandler(int)
  
cdef extern from "datum.h":
  cdef cppclass Datum:
        pass

cdef extern from "namedatum.h":
    cdef cppclass NameDatum:
        NameDatum(string)
        NameDatum(NameDatum)
        string toString()

cdef extern from "token.h":
    cdef cppclass Token:
        Token(Token)
        Token(Datum *)

cdef extern from "cynestkernel.h":
    cdef cppclass NESTEngine:
        NESTEngine()
        bint init(vector[string] argv, string)
        bint push(object)
        bint push_token(Token)
        bint push_connections(object)
        bint run(string)
        bint run_token(Token)
        object pop()
        Token* pop_token()
        Datum* PyObject_as_Datum(object)
        bint check_engine()
        void register_cython_model(string)


cdef extern from "object_manager.h":
    cdef cppclass Time:
        long get_tics()
        long get_steps()
        double get_ms()
        void set_to_zero()
        void calibrate()
        void advance()
        Time succ()
        Time pred()
        bint is_grid_time()
        bint is_neg_inf()
        bint is_pos_inf()
        bint is_finite()
        bint is_step()


    cdef cppclass UnitManager:
        UnitManager(int, long)
        UnitManager(int, double)
        Time generateTime()

cdef extern from "object_manager.h" namespace "nest::Time":
    Time get_resolution()
    void set_resolution(double)
    void reset_resolution()
    bint resolution_is_default()
    double get_ms_per_tic()
    double get_tics_per_ms()
    long get_tics_per_step()
    long get_old_tics_per_step()
    long get_tics_per_step_default()
    Time min()
    Time max()
    Time pos_inf()
    Time neg_inf()





cdef extern from "object_manager.h" namespace "nest::Scheduler":
    unsigned int get_modulo(unsigned int)
    unsigned int get_slice_modulo(unsigned int)
    unsigned int get_min_delay()
    unsigned int get_max_delay()
