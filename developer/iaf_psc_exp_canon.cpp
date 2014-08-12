/*
 *  iaf_psc_exp_canon.cpp
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2007-2008 by
 *  The NEST Initiative 
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#include "iaf_psc_exp_canon.h"

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "analog_data_logger_impl.h"

#include <limits>

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
nest::iaf_psc_exp_canon::Parameters_::Parameters_()
  : tau_m_  ( 10.0     ), // [ms]
    tau_ex_ (  2.0     ), // [ms]
    tau_in_ (  1.0     ), // [ms]
    c_m_    (250.0     ), // [pF]
    t_ref_  (  2.0     ), // [ms]
    E_L_    (-70.0     ), // [mV]
    I_e_    (  0.0     ), // [pA]
    U_th_   (-55.0-E_L_), // [mV], rel. to E_L_
    U_min_  (-std::numeric_limits<double_t>::infinity()), // [mV]
    U_reset_(-70.0-E_L_), // [mV], rel. to E_L_
    Interpol_(iaf_psc_exp_canon::LINEAR)
{}

nest::iaf_psc_exp_canon::State_::State_()
  : y1_ex_(0.0),
    y1_in_(0.0),
    y2_(0.0),  
    is_refractory_(false),
    last_spike_step_(-1),
    last_spike_offset_(0.0)
{}

/* ---------------------------------------------------------------- 
 * Paramater and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::iaf_psc_exp_canon::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d, names::tau_m, tau_m_);
  def<double>(d, names::tau_syn_ex, tau_ex_);
  def<double>(d, names::tau_syn_in, tau_in_);
  def<double>(d, names::C_m, c_m_);
  def<double>(d, names::t_ref, t_ref_);
  def<double>(d, names::E_L, E_L_);
  def<double>(d, names::I_e, I_e_);
  def<double>(d, names::V_th, U_th_+E_L_);
  def<double>(d, names::V_min, U_min_+E_L_);
  def<double>(d, names::V_reset, U_reset_+E_L_);
  def<long>(d, names::Interpol_Order, Interpol_);
}

void nest::iaf_psc_exp_canon::Parameters_::set(const DictionaryDatum& d)
{
  updateValue<double>(d, names::tau_m, tau_m_);
  updateValue<double>(d, names::tau_syn_ex, tau_ex_);
  updateValue<double>(d, names::tau_syn_in, tau_in_);
  updateValue<double>(d, names::C_m, c_m_);
  updateValue<double>(d, names::t_ref, t_ref_);
  updateValue<double>(d, names::E_L, E_L_);
  updateValue<double>(d, names::I_e, I_e_);

  if ( updateValue<double>(d, names::V_th, U_th_) ) 
    U_th_ -= E_L_;

  if ( updateValue<double>(d, names::V_min, U_min_) ) 
    U_min_ -= E_L_;

  if ( updateValue<double>(d, names::V_reset, U_reset_) ) 
    U_reset_ -= E_L_;
  
  long_t tmp;
  if ( updateValue<long_t>(d, names::Interpol_Order, tmp) )
    if ( NO_INTERPOL <= tmp && tmp < END_INTERP_ORDER )
      Interpol_ = static_cast<interpOrder>(tmp);
    else
      throw BadProperty("Interpolation order must be 0, 1, 2 or 3.");

  if ( tau_m_ <= 0 || tau_ex_ <= 0 ||  tau_in_ <= 0 )
    throw BadProperty("All time constants must be strictly positive.");
    
  if ( c_m_ <= 0 )
    throw BadProperty("Capacitance must be strictly positive.");
    
  if ( t_ref_ < 0 )
    throw BadProperty("Refractory time must not be negative.");

  if ( U_reset_ >= U_th_ )
    throw BadProperty("Reset potential must be smaller than threshold.");

  if ( U_reset_ < U_min_ )
    throw BadProperty("Reset potential must be greater equal minimum potential.");
}

void nest::iaf_psc_exp_canon::State_::get(DictionaryDatum &d, 
                                            Parameters_ const &p) const
{
  def<double>(d, names::V_m, y2_ + p.E_L_); // Membrane potential
  def<double>(d, names::t_spike, Time(Time::step(last_spike_step_)).get_ms());
  def<double>(d, names::offset, last_spike_offset_);
  def<bool>(d, names::is_refractory, is_refractory_);
}

void nest::iaf_psc_exp_canon::State_::set(DictionaryDatum const &d,
                                            Parameters_ const &p)
{
  if ( updateValue<double>(d, names::V_m, y2_) )
    y2_ -= p.E_L_;
    
  if ( y2_ < p.U_min_ )
    throw BadProperty("Membrane potential cannot be smaller than minimum potential.");
}

/* ---------------------------------------------------------------- 
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::iaf_psc_exp_canon::iaf_psc_exp_canon()
  : Node(), 
    P_(), 
    S_()
{}

nest::iaf_psc_exp_canon::iaf_psc_exp_canon(iaf_psc_exp_canon const &n)
  : Node(n), 
    P_(n.P_), 
    S_(n.S_)
{}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::iaf_psc_exp_canon::init_state_(Node const &proto)
{
  iaf_psc_exp_canon const &pr = downcast<iaf_psc_exp_canon>(proto);
  S_ = pr.S_;
}

void nest::iaf_psc_exp_canon::init_buffers_()
{
  B_.events_.resize();
  B_.events_.clear(); 
  B_.potentials_.clear_data(); // includes resize
}

void nest::iaf_psc_exp_canon::calibrate()
{
  V_.h_ms_ = Time::get_resolution().get_ms();

  // pre-compute matrix for full time step
  V_.expm1_tau_m_  = numerics::expm1(-V_.h_ms_/P_.tau_m_);
  V_.expm1_tau_ex_ = numerics::expm1(-V_.h_ms_/P_.tau_ex_); 
  V_.expm1_tau_in_ = numerics::expm1(-V_.h_ms_/P_.tau_in_);
  V_.P20_    = -P_.tau_m_/P_.c_m_ * V_.expm1_tau_m_;
  V_.P21_ex_ =  P_.tau_m_/P_.c_m_ * 1/P_.tau_m_ / (1/P_.tau_m_ - 1/P_.tau_ex_)
                                  * (V_.expm1_tau_ex_-V_.expm1_tau_m_);
  V_.P21_in_ =  P_.tau_m_/P_.c_m_ * 1/P_.tau_m_ / (1/P_.tau_m_ - 1/P_.tau_in_)
                                  * (V_.expm1_tau_in_-V_.expm1_tau_m_);

  // t_ref_ is the refractory period in ms
  // refractory_steps_ is the duration of the refractory period in whole
  // steps, rounded down
  V_.refractory_steps_ = Time(Time::ms(P_.t_ref_)).get_steps();
  assert ( V_.refractory_steps_ >= 0 ); // since t_ref_ >= 0, this can only fail in error
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void nest::iaf_psc_exp_canon::update(Time const &origin,
                                       long_t const from, long_t const to)
{
  assert ( to >= 0 );
  assert ( static_cast<delay>(from) < Scheduler::get_min_delay() );
  assert ( from < to );

  // at start of slice, tell input queue to prepare for delivery
  if ( from == 0 )
    B_.events_.prepare_delivery();

  // Neurons may have been initialized to superthreshold potentials.
  // We need to check for this here and issue spikes at the beginning of
  // the interval.
  if ( S_.y2_ >= P_.U_th_ )
    emit_instant_spike_(origin, from,
                          V_.h_ms_*(1-std::numeric_limits<double_t>::epsilon()));

  for ( long_t lag = from; lag < to; ++lag )
  {
    // time at start of update step
    long_t const T = origin.get_steps() + lag;

    // if neuron returns from refractoriness during this step, place
    // pseudo-event in queue to mark end of refractory period
    if ( S_.is_refractory_ && (T + 1 - S_.last_spike_step_ == V_.refractory_steps_) )
      B_.events_.add_refractory(T, S_.last_spike_offset_);

    // save state at beginning of interval for spike-time computation
    V_.y1_ex_before_ = S_.y1_ex_;
    V_.y1_in_before_ = S_.y1_in_;
    V_.y2_before_    = S_.y2_;

    // get first event
    double_t ev_offset;
    double_t ev_weight;
    bool end_of_refract;

    if ( !B_.events_.get_next_spike(T, ev_offset, ev_weight, end_of_refract) )
    {
    // No incoming spikes, handle with fixed propagator matrix.
    // Handling this case separately improves performance significantly
    // if there are many steps without input spikes.

      if ( !S_.is_refractory_ )
      {
        // update membrane potential
        S_.y2_ = V_.P20_*P_.I_e_ + V_.P21_ex_*S_.y1_ex_
                                 + V_.P21_in_*S_.y1_in_
                                 + V_.expm1_tau_m_*S_.y2_ + S_.y2_;

        // lower bound of membrane potential
        S_.y2_ = ( S_.y2_ < P_.U_min_ ? P_.U_min_ : S_.y2_ ); 
      }

      // update synaptic current
      S_.y1_ex_ = V_.expm1_tau_ex_*S_.y1_ex_ + S_.y1_ex_;
      S_.y1_in_ = V_.expm1_tau_in_*S_.y1_in_ + S_.y1_in_;

      // The following must not be moved before the y1_ update,
      // since the spike-time computation within emit_spike_ depends
      // on all state variables having their values at the end of the
      // interval.
      if ( S_.y2_ >= P_.U_th_ )
        emit_spike_(origin, lag, 0, V_.h_ms_);
    }
    else
    {
      // We only get here if there is at least one event,
      // which has been read above. We can therefore use
      // a do-while loop.

      // Time within step is measured by offsets, which are h at the beginning
      // and 0 at the end of the step.
      double_t last_offset = V_.h_ms_; // start of step

      do
      {
        // time is measured backward: inverse order in difference
        double_t const ministep = last_offset-ev_offset;

        propagate_(ministep);

        // check for threshold crossing during ministep
        // this must be done before adding the input, since
        // interpolation requires continuity
        if ( S_.y2_ >= P_.U_th_ )
          emit_spike_(origin, lag, V_.h_ms_-last_offset, ministep);

        // handle event
        if ( end_of_refract )
          S_.is_refractory_ = false; // return from refractoriness
        else
        {
          if ( ev_weight >= 0.0 )
            S_.y1_ex_ += ev_weight; // spike input ex.
	  else
            S_.y1_in_ += ev_weight; // spike input in.
	}

        // store state
        V_.y1_ex_before_ = S_.y1_ex_;
        V_.y1_in_before_ = S_.y1_in_;
        V_.y2_before_    = S_.y2_;
        last_offset      = ev_offset;
      }
      while ( B_.events_.get_next_spike(T, ev_offset, ev_weight,
                                          end_of_refract) );

      // no events remaining, plain update step across remainder 
      // of interval
      if ( last_offset > 0 ) // not at end of step, do remainder
      {
	propagate_(last_offset);
        if ( S_.y2_ >= P_.U_th_ )
          emit_spike_(origin, lag, V_.h_ms_-last_offset, last_offset);
      }
    } // else

    // voltage logging
    B_.potentials_.record_data(origin.get_steps()+lag, S_.y2_+P_.E_L_);
  }
}

// function handles exact spike times
void nest::iaf_psc_exp_canon::handle(SpikeEvent &e)
{
  assert ( e.get_delay() > 0 );

  // We need to compute the absolute time stamp of the delivery time
  // of the spike, since spikes might spend longer than min_delay_
  // in the queue.  The time is computed according to Time Memo, Rule 3.
  long_t const Tdeliver = e.get_stamp().get_steps() + e.get_delay() - 1;

  B_.events_.add_spike(e.get_rel_delivery_steps(network()->get_slice_origin()),
                         Tdeliver, e.get_offset(), e.get_weight() * e.get_multiplicity());
}

void nest::iaf_psc_exp_canon::handle(PotentialRequest &e)
{
  B_.potentials_.handle(*this, e);
}

/* ---------------------------------------------------------------- 
 * Auxiliary functions
 * ---------------------------------------------------------------- */

inline
nest::Time nest::iaf_psc_exp_canon::get_spiketime() const
{
  return Time::step(S_.last_spike_step_);
}

inline
void nest::iaf_psc_exp_canon::set_spiketime(Time const &now)
{
  S_.last_spike_step_ = now.get_steps();
}

void nest::iaf_psc_exp_canon::propagate_(double_t const dt)
{
  double_t const ps_e_TauEx = numerics::expm1(-dt/P_.tau_ex_); // needed in any case 
  double_t const ps_e_TauIn = numerics::expm1(-dt/P_.tau_in_); // needed in any case 

  // y2_ remains unchanged at 0.0 while neuron is refractory
  if ( !S_.is_refractory_ )
  {
    double_t const ps_e_Tau  = numerics::expm1(-dt/P_.tau_m_);
    double_t const ps_P20    = -P_.tau_m_/P_.c_m_ * ps_e_Tau;
    double_t const ps_P21_ex =  P_.tau_m_/P_.c_m_ * 1/P_.tau_m_ / (1/P_.tau_m_ - 1/P_.tau_ex_)
                                                  * (ps_e_TauEx-ps_e_Tau);
    double_t const ps_P21_in =  P_.tau_m_/P_.c_m_ * 1/P_.tau_m_ / (1/P_.tau_m_ - 1/P_.tau_in_)
                                                  * (ps_e_TauIn-ps_e_Tau);

    S_.y2_ = ps_P20*P_.I_e_ + ps_P21_ex*S_.y1_ex_
                            + ps_P21_in*S_.y1_in_
                            + ps_e_Tau*S_.y2_ + S_.y2_;

    // lower bound of membrane potential
    S_.y2_ = ( S_.y2_ < P_.U_min_ ? P_.U_min_ : S_.y2_); 
  }

  // now the synaptic components
  S_.y1_ex_ = ps_e_TauEx*S_.y1_ex_ + S_.y1_ex_;
  S_.y1_in_ = ps_e_TauIn*S_.y1_in_ + S_.y1_in_;

  return;
}

void nest::iaf_psc_exp_canon::emit_spike_(Time const &origin, long_t const lag, 
                                            double_t const t0,  double_t const dt)
{
  // we know that the potential is subthreshold at t0, super at t0+dt

  // compute spike time relative to beginning of step
  double_t const spike_offset = V_.h_ms_ - (t0+thresh_find_(dt));
  set_spiketime(Time::step(origin.get_steps() + lag + 1));
  S_.last_spike_offset_ = spike_offset;

  // reset neuron and make it refractory
  S_.y2_ = P_.U_reset_;
  S_.is_refractory_ = true; 

  // send spike
  SpikeEvent se;
  se.set_offset(S_.last_spike_offset_);
  network()->send(*this, se, lag);

  return;
}

void nest::iaf_psc_exp_canon::emit_instant_spike_(Time const &origin, long_t const lag,
                                                    double_t const spike_offset) 
{
  assert ( S_.y2_ >= P_.U_th_ ); // ensure we are superthreshold

  // set stamp and offset for spike
  set_spiketime(Time::step(origin.get_steps() + lag + 1));
  S_.last_spike_offset_ = spike_offset;

  // reset neuron and make it refractory
  S_.y2_ = P_.U_reset_;
  S_.is_refractory_ = true;

  // send spike
  SpikeEvent se;
  se.set_offset(S_.last_spike_offset_);
  network()->send(*this, se, lag);

  return;
}

/* ---------------------------------------------------------------- 
 * Functions for threshold crossing interpolation
 * ---------------------------------------------------------------- */

// switch
inline
nest::double_t nest::iaf_psc_exp_canon::thresh_find_(double_t const dt) const
{
  switch ( P_.Interpol_ )
  {
    case NO_INTERPOL: return dt;
    case LINEAR     : return thresh_find1_(dt);
    case QUADRATIC  : return thresh_find2_(dt);
    case CUBIC      : return thresh_find3_(dt);
    default:
      throw BadProperty("Invalid interpolation order in iaf_psc_exp_canon.");
  }
  return 0;
}

// linear interpolation
nest::double_t nest::iaf_psc_exp_canon::thresh_find1_(double_t const dt) const
{
  double_t tau = (P_.U_th_-V_.y2_before_) * dt / (S_.y2_-V_.y2_before_);

  return tau;
}

// quadratic interpolation  
nest::double_t nest::iaf_psc_exp_canon::thresh_find2_(double_t const dt) const
{
  double_t const h_sq = dt*dt;

  // f' at left border used as third condition
  double_t const derivative = -V_.y2_before_/P_.tau_m_
                                + (P_.I_e_ + V_.y1_ex_before_ + V_.y1_in_before_) / P_.c_m_; 

  double_t const a = (-V_.y2_before_/h_sq) + (S_.y2_/h_sq) - (derivative/dt);
  double_t const b = derivative;
  double_t const c = V_.y2_before_;

  double_t const sqr_ = std::sqrt(b*b - 4*a*c + 4*a*P_.U_th_);
  double_t const tau1 = (-b+sqr_) / (2*a);
  double_t const tau2 = (-b-sqr_) / (2*a);

  if ( tau1 >= 0 )
    return tau1;
  else if ( tau2 >= 0 )
    return tau2;
  else
    return thresh_find1_(dt);
}

// cubic interpolation
nest::double_t nest::iaf_psc_exp_canon::thresh_find3_(double_t const dt) const
{
  double_t const h_ms_ = dt;
  double_t const h_sq  = h_ms_*h_ms_;
  double_t const h_cb  = h_sq*h_ms_;

  // f' at left border used as third condition
  double_t const deriv_t1 = -V_.y2_before_/P_.tau_m_
                              + (P_.I_e_ + V_.y1_ex_before_ + V_.y1_in_before_) / P_.c_m_;

  // f' at right border used as fourth condition
  double_t const deriv_t2 = -S_.y2_/P_.tau_m_
                              + (P_.I_e_ + S_.y1_ex_ + S_.y1_in_) / P_.c_m_;

  double_t const w3_ =  (2*V_.y2_before_ / h_cb) - (2*S_.y2_ / h_cb)
                         + (deriv_t1 / h_sq) + (deriv_t2 / h_sq) ;
  double_t const w2_ = -(3*V_.y2_before_ / h_sq) + (3*S_.y2_ / h_sq)
                         - (2*deriv_t1 / h_ms_) - (deriv_t2 / h_ms_);
  double_t const w1_ = deriv_t1;
  double_t const w0_ = V_.y2_before_;

  // normal form: x^3 + r*x^2 + s*x + t with coefficients: r, s, t
  double_t const r    = w2_/w3_;
  double_t const s    = w1_/w3_;
  double_t const t    = (w0_-P_.U_th_) / w3_;
  double_t const r_sq = r*r;

  // substitution y = x + r/3:  y^3 + p*y + q == 0 
  double_t const p = -r_sq/3 + s;
  double_t const q = 2*(r_sq*r) / 27 - r*s / 3 + t;

  // discriminante
  double_t const D = std::pow((p/3), 3) + std::pow((q/2), 2);

  double_t tau1;
  double_t tau2;
  double_t tau3;

  if ( D < 0 )
  {
    double_t const roh = std::sqrt(-(p*p*p) / 27);
    double_t const phi = std::acos(-q / (2*roh));
    double_t const a   = 2*std::pow(roh, (1.0/3.0));

    tau1 = (a*std::cos(phi/3)) - r/3;
    tau2 = (a*std::cos(phi/3 + 2*numerics::pi/3)) - r/3;
    tau3 = (a*std::cos(phi/3 + 4*numerics::pi/3)) - r/3;
  }
  else
  {
    double_t const sgnq = ( q >= 0 ? 1 : -1 );
    double_t const u    = -sgnq*std::pow(std::fabs(q)/2.0 + std::sqrt(D), 1.0/3.0);
    double_t const v    = -p/(3*u);

    tau1= (u+v) - r/3;
    if ( tau1 >= 0 )
      return tau1;
    else
      return thresh_find2_(dt);
  }

  // set tau to the smallest root above 0
  double tau = ( tau1 >= 0 ? tau1 : 2*h_ms_ );

  if (( tau2 >=0 ) && ( tau2 < tau ))
    tau = tau2;
  if (( tau3 >=0 ) && ( tau3 < tau ))
    tau = tau3;
  return ( tau <= h_ms_ ? tau : thresh_find2_(dt) );
}
