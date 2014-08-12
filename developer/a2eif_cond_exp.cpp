/*
 *  a2eif_cond_exp.cpp
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2010 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#include "a2eif_cond_exp.h"
#include "nest_names.h"

#ifdef HAVE_GSL_1_11

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include "universal_data_logger_impl.h"
#include <limits>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <cstdio>

/* ---------------------------------------------------------------- 
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap<nest::a2eif_cond_exp> nest::a2eif_cond_exp::recordablesMap_;

namespace nest
{
  /*
   * template specialization must be placed in namespace
   *
   * Override the create() method with one call to RecordablesMap::insert_() 
   * for each quantity to be recorded.
   */
  template <>
  void RecordablesMap<a2eif_cond_exp>::create()
  {
    // use standard names whereever you can for consistency!
    insert_(names::V_m, 
	    &a2eif_cond_exp::get_y_elem_<a2eif_cond_exp::State_::V_M>);
    insert_(names::g_ex, 
	    &a2eif_cond_exp::get_y_elem_<a2eif_cond_exp::State_::G_EXC>);
    insert_(names::g_in, 
	    &a2eif_cond_exp::get_y_elem_<a2eif_cond_exp::State_::G_INH>);
    insert_("w1", 
	    &a2eif_cond_exp::get_y_elem_<a2eif_cond_exp::State_::W1>);
    insert_("w2", 
	    &a2eif_cond_exp::get_y_elem_<a2eif_cond_exp::State_::W2>);
  }
}

extern "C"
int nest::a2eif_cond_exp_dynamics (double, const double y[], double f[], void* pnode)
{
  // a shorthand
  typedef nest::a2eif_cond_exp::State_ S;

  // get access to node so we can almost work as in a member function
  assert(pnode);
  const nest::a2eif_cond_exp& node =  *(reinterpret_cast<nest::a2eif_cond_exp*>(pnode));

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[]. 
  
  // The following code is verbose for the sake of clarity. We assume that a
  // good compiler will optimize the verbosity away ...

  // shorthand for state variables
  const double_t& V     = y[S::V_M  ];
  const double_t& g_ex  = y[S::G_EXC];
  const double_t& g_in  = y[S::G_INH];
  const double_t& w1    = y[S::W1   ];
  const double_t& w2    = y[S::W2   ];

  const double_t I_syn_exc = g_ex * (V - node.P_.E_ex);
  const double_t I_syn_inh = g_in * (V - node.P_.E_in);
  const double_t I_spike   = node.P_.Delta_T * std::exp((V - node.P_.V_th) / node.P_.Delta_T);

  // Membrane potential dV/dt
  f[S::V_M  ] = ( -node.P_.g_L *( (V-node.P_.E_L) - I_spike) 
	       - I_syn_exc - I_syn_inh - w1 - w2 + node.P_.I_e + node.B_.I_stim_) / node.P_.C_m;
  // Synaptic conductances dg_ex/dt, dg_in/dt
  f[S::G_EXC ] = -g_ex/ node.P_.tau_syn_ex;
  f[S::G_INH ] = -g_in/ node.P_.tau_syn_in;

  // Adaptation currents dw1/dt, dw2/dt.
  f[S::W1   ] = ( node.P_.a1 * (V - node.P_.E_L) - w1 ) / node.P_.tau_w1;
  f[S::W2   ] = ( node.P_.a2 * (V - node.P_.E_L) - w2 ) / node.P_.tau_w2;

  return GSL_SUCCESS;
}

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
nest::a2eif_cond_exp::Parameters_::Parameters_()
  : V_peak_    (   0.0    ),  // mV
    V_reset_   ( -50.0    ),  // mV
    t_ref_     (   3.0    ),  // ms
    g_L        (  14.0    ),  // nS
    C_m        ( 260.0    ),  // pF
    E_ex       (   0.0    ),  // mV
    E_in       ( -85.0    ),  // mV
    E_L        ( -69.0    ),  // mV
    Delta_T    (   2.4    ),  // mV
    tau_w1     (   5.5    ),  // ms
    a1         (   0.0    ),  // nS
    b1         (1250.0    ),  // pA
    tau_w2     ( 100.0    ),  // ms
    a2         (   0.0    ),  // nS
    b2         ( 110.0    ),  // pA
    V_th       ( -47.0    ),  // mV
    tau_syn_ex (   0.2    ),  // ms
    tau_syn_in (   2.0    ),  // ms
    I_e        (   0.0    )   // pA
{
  recordablesMap_.create();
}

nest::a2eif_cond_exp::State_::State_(const Parameters_ &p)
  : r_(0)
{
  y_[0] = p.E_L;  // initialize to reversal potential
  for ( size_t i = 1; i < STATE_VEC_SIZE; ++i )
    y_[i] = 0;
}

nest::a2eif_cond_exp::State_::State_(const State_ &s)
  : r_(s.r_)
{
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
    y_[i] = s.y_[i];
}

nest::a2eif_cond_exp::State_& nest::a2eif_cond_exp::State_::operator=(const State_ &s)
{
  assert(this != &s);  // would be bad logical error in program
  
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
    y_[i] = s.y_[i];
  r_ = s.r_;
  return *this;
}

/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::a2eif_cond_exp::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d,names::C_m,        C_m);
  def<double>(d,names::V_th,       V_th);
  def<double>(d,names::t_ref,      t_ref_);
  def<double>(d,names::g_L,        g_L);
  def<double>(d,names::E_L,        E_L); 
  def<double>(d,names::V_reset,    V_reset_);
  def<double>(d,names::E_ex,       E_ex);
  def<double>(d,names::E_in,       E_in);
  def<double>(d,names::tau_syn_ex, tau_syn_ex);
  def<double>(d,names::tau_syn_in, tau_syn_in);
  def<double>(d,"a1",              a1);
  def<double>(d,"b1",              b1);
  def<double>(d,"tau_w1",          tau_w1);
  def<double>(d,"a2",              a2);
  def<double>(d,"b2",              b2);
  def<double>(d,"tau_w2",          tau_w2);
  def<double>(d,names::Delta_T,    Delta_T);
  def<double>(d,names::I_e,        I_e);
  def<double>(d,names::V_peak,     V_peak_);
}

void nest::a2eif_cond_exp::Parameters_::set(const DictionaryDatum &d)
{
  updateValue<double>(d,names::V_th,       V_th);
  updateValue<double>(d,names::V_peak,     V_peak_);
  updateValue<double>(d,names::t_ref,      t_ref_);
  updateValue<double>(d,names::E_L,        E_L);
  updateValue<double>(d,names::V_reset,    V_reset_);
  updateValue<double>(d,names::E_ex,       E_ex);
  updateValue<double>(d,names::E_in,       E_in);
  updateValue<double>(d,names::C_m,        C_m);
  updateValue<double>(d,names::g_L,        g_L);
  updateValue<double>(d,names::tau_syn_ex, tau_syn_ex);
  updateValue<double>(d,names::tau_syn_in, tau_syn_in);
  updateValue<double>(d,"a1",              a1);
  updateValue<double>(d,"b1",              b1);
  updateValue<double>(d,"tau_w1",          tau_w1);
  updateValue<double>(d,"a2",              a2);
  updateValue<double>(d,"b2",              b2);
  updateValue<double>(d,"tau_w2",          tau_w2);
  updateValue<double>(d,names::Delta_T,    Delta_T);
  updateValue<double>(d,names::I_e,        I_e);

  //  Restrictions for physical model consistency
  if ( V_reset_ >= V_peak_ ) 
    throw BadProperty("Reset potential must be smaller than spike cut-off threshold.");     

  if ( C_m <= 0 )
    throw BadProperty("Capacitance must be strictly positive.");
    
  if ( t_ref_ < 0 )
    throw BadProperty("Refractory time cannot be negative.");
      
  if ( tau_syn_ex <= 0 || tau_syn_in <= 0 || tau_w1 <= 0 || tau_w2 <= 0 )
    throw BadProperty("All time constants must be strictly positive.");
}

void nest::a2eif_cond_exp::State_::get(DictionaryDatum &d) const
{
  def<double>(d,names::V_m,  y_[V_M]);
  def<double>(d,names::g_ex, y_[G_EXC]);
  def<double>(d,names::g_in, y_[G_INH]);
  def<double>(d,"w1",        y_[W1]);
  def<double>(d,"w2",        y_[W2]);
}

void nest::a2eif_cond_exp::State_::set(const DictionaryDatum &d, const Parameters_ &)
{
  updateValue<double>(d,names::V_m,  y_[V_M]);
  updateValue<double>(d,names::g_ex, y_[G_EXC]);
  updateValue<double>(d,names::g_in, y_[G_INH]);
  updateValue<double>(d,"w1",        y_[W1]);
  updateValue<double>(d,"w2",        y_[W2]);

  if ( y_[G_EXC] < 0 || y_[G_INH] < 0 )
    throw BadProperty("Conductances must not be negative.");
}

nest::a2eif_cond_exp::Buffers_::Buffers_(a2eif_cond_exp &n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

nest::a2eif_cond_exp::Buffers_::Buffers_(const Buffers_ &, a2eif_cond_exp &n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ---------------------------------------------------------------- 
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

nest::a2eif_cond_exp::a2eif_cond_exp()
  : Archiving_Node(), 
    P_(), 
    S_(P_),
    B_(*this)
{
}

nest::a2eif_cond_exp::a2eif_cond_exp(const a2eif_cond_exp &n)
  : Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{
}

nest::a2eif_cond_exp::~a2eif_cond_exp()
{
  // GSL structs may not have been allocated, so we need to protect destruction
  if ( B_.s_ ) gsl_odeiv_step_free(B_.s_);
  if ( B_.c_ ) gsl_odeiv_control_free(B_.c_);
  if ( B_.e_ ) gsl_odeiv_evolve_free(B_.e_);
}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::a2eif_cond_exp::init_state_(const Node &proto)
{
  const a2eif_cond_exp &pr = downcast<a2eif_cond_exp>(proto);
  S_ = pr.S_;
}

void nest::a2eif_cond_exp::init_buffers_()
{
  B_.spike_exc_.clear();       // includes resize
  B_.spike_inh_.clear();       // includes resize
  B_.currents_.clear();        // includes resize
  Archiving_Node::clear_history();

  B_.logger_.reset();

  B_.step_ = Time::get_resolution().get_ms();

  // We must integrate this model with high-precision to obtain decent results
  B_.IntegrationStep_ = std::min(0.01, B_.step_);

  static const gsl_odeiv_step_type* T1 = gsl_odeiv_step_rkf45;
  
  if ( B_.s_ == 0 )
    B_.s_ = gsl_odeiv_step_alloc (T1, State_::STATE_VEC_SIZE);
  else 
    gsl_odeiv_step_reset(B_.s_);

  if ( B_.c_ == 0 )  
    B_.c_ = gsl_odeiv_control_yp_new (1e-6,1e-6);
  else
    gsl_odeiv_control_init(B_.c_, 1e-6, 1e-6, 0.0, 1.0);
    
  if ( B_.e_ == 0 )  
    B_.e_ = gsl_odeiv_evolve_alloc(State_::STATE_VEC_SIZE);
  else 
    gsl_odeiv_evolve_reset(B_.e_);
  
  B_.sys_.function  = a2eif_cond_exp_dynamics; 
  B_.sys_.jacobian  = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params    = reinterpret_cast<void*>(this);

  B_.I_stim_ = 0.0;
}

void nest::a2eif_cond_exp::calibrate()
{
  B_.logger_.init();  // ensures initialization in case mm connected after Simulate
  V_.RefractoryCounts_ = Time(Time::ms(P_.t_ref_)).get_steps();
  assert(V_.RefractoryCounts_ >= 0);  // since t_ref_ >= 0, this can only fail in error
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void nest::a2eif_cond_exp::update(const Time &origin, const long_t from, const long_t to)
{
  assert ( to >= 0 && (delay) from < Scheduler::get_min_delay() );
  assert ( from < to );
  assert ( State_::V_M == 0 );

  for ( long_t lag = from; lag < to; ++lag )
  {
    double t = 0.0;

    if ( S_.r_ > 0 )
      --S_.r_;

    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t
    while ( t < B_.step_ )
    {
      const int status = gsl_odeiv_evolve_apply(B_.e_, B_.c_, B_.s_, 
						&B_.sys_,             // system of ODE
						&t,                   // from t
						B_.step_,             // to t <= step
						&B_.IntegrationStep_, // integration step size
						S_.y_);               // neuronal state

      if ( status != GSL_SUCCESS )
        throw GSLSolverFailure(get_name(), status);

      // spikes are handled inside the while-loop
      // due to spike-driven adaptation
      if ( S_.r_ > 0 ) // refractoriness => clamp potential
        S_.y_[State_::V_M] = P_.V_reset_;
      else if ( S_.y_[State_::V_M] >= P_.V_peak_ ) // no refractoriness
      {
	S_.y_[State_::V_M]  = P_.V_reset_;
	S_.y_[State_::W1]   += P_.b1; // spike-driven adaptation
	S_.y_[State_::W2]   += P_.b2; // spike-driven adaptation
	S_.r_               = V_.RefractoryCounts_;
	
	set_spiketime(Time::step(origin.get_steps() + lag + 1));
	SpikeEvent se;
	network()->send(*this, se, lag);
      }
    }

    // add incoming spikes
    S_.y_[State_::G_EXC] += B_.spike_exc_.get_value(lag);
    S_.y_[State_::G_INH] += B_.spike_inh_.get_value(lag);
    // set new input current
    B_.I_stim_ = B_.currents_.get_value(lag);
    // log state data
    B_.logger_.record_data(origin.get_steps() + lag);
  }
}
  
void nest::a2eif_cond_exp::handle(SpikeEvent &e)
{
  assert ( e.get_delay() > 0 );

  if(e.get_weight() > 0.0)
    B_.spike_exc_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			    e.get_weight() * e.get_multiplicity());
  else
    B_.spike_inh_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			    -e.get_weight() * e.get_multiplicity());  // keep conductances positive 
}

void nest::a2eif_cond_exp::handle(CurrentEvent &e)
{
  assert ( e.get_delay() > 0 );

  const double_t c=e.get_current();
  const double_t w=e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), w*c);
}

void nest::a2eif_cond_exp::handle(DataLoggingRequest &e)
{
  B_.logger_.handle(e);
}

#endif //HAVE_GSL_1_11
