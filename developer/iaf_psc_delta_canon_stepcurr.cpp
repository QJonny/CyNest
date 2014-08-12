/*
 *  iaf_psc_delta_canon_stepcurr.cpp
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2004-2006 by
 *  The NEST Initiative 
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

/* iaf_psc_delta_canon_stepcurr is a neuron where the potential jumps on each spike arrival. */

#include "config.h"

#include "exceptions.h"
#include "iaf_psc_delta_canon_stepcurr.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"  
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include "analog_data_logger_impl.h"

#include <limits>

namespace nest
{

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
iaf_psc_delta_canon_stepcurr::Parameters_::Parameters_()
  : tau_m_  ( 10.0    ),  // ms
    c_m_    (250.0    ),  // pF
    t_ref_  (  2.0    ),  // ms
    E_L_    (-70.0    ),  // mV
    I_e_    (  0.0    ),  // pA
    U_th_   (-55.0-E_L_), // mV, rel to E_L_
    U_min_  (-std::numeric_limits<double_t>::max()),  // mV
    U_reset_(-70.0-E_L_)  // mV, rel to E_L_
{}

iaf_psc_delta_canon_stepcurr::State_::State_()
  : U_(0.0), //  or U_ = U_reset_;
    last_spike_step_(-1),
    last_spike_offset_(0.0),
    is_refractory_(false),
    with_refr_input_(false)
{}

/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void iaf_psc_delta_canon_stepcurr::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d, names::E_L, E_L_);
  def<double>(d, names::I_e, I_e_);
  def<double>(d, names::V_th, U_th_+E_L_);
  def<double>(d, names::V_min, U_min_+E_L_);
  def<double>(d, names::V_reset, U_reset_+E_L_);
  def<double>(d, names::C_m, c_m_);
  def<double>(d, names::tau_m, tau_m_);
  def<double>(d, names::t_ref, t_ref_);
}

void iaf_psc_delta_canon_stepcurr::Parameters_::set(const DictionaryDatum& d)
{
  updateValue<double>(d, names::tau_m, tau_m_);
  updateValue<double>(d, names::C_m, c_m_);
  updateValue<double>(d, names::t_ref, t_ref_);
  updateValue<double>(d, names::E_L, E_L_);
  updateValue<double>(d, names::I_e, I_e_);

  if (updateValue<double>(d, names::V_th, U_th_)) 
    U_th_ -= E_L_;

  if (updateValue<double>(d, names::V_min, U_min_)) 
    U_min_ -= E_L_;

  if (updateValue<double>(d, names::V_reset, U_reset_)) 
    U_reset_ -= E_L_;
  
  if ( U_reset_ >= U_th_ )
    throw BadProperty("Reset potential must be smaller than threshold.");

  if ( U_reset_ < U_min_ )
    throw BadProperty("Reset potential must be greater equal minimum potential.");
    
  if ( c_m_ <= 0 )
    throw BadProperty("Capacitance must be strictly positive.");
    
  if ( t_ref_ < 0 )
    throw BadProperty("Refractory time must not be negative.");
    
  if ( tau_m_ <= 0 )
    throw BadProperty("All time constants must be strictly positive.");
}

void iaf_psc_delta_canon_stepcurr::State_::get(DictionaryDatum &d, 
                                            const Parameters_& p) const
{
  def<double>(d, names::V_m, U_ + p.E_L_); // Membrane potential
  def<double>(d, names::t_spike, Time(Time::step(last_spike_step_)).get_ms());
  def<double>(d, names::offset, last_spike_offset_);
  def<bool>(d, names::is_refractory, is_refractory_);
  def<bool>(d, names::refractory_input, with_refr_input_);
}

void iaf_psc_delta_canon_stepcurr::State_::set(const DictionaryDatum& d, const Parameters_& p)
{
  if ( updateValue<double>(d, names::V_m, U_) )
    U_ -= p.E_L_;
}

/* ---------------------------------------------------------------- 
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

iaf_psc_delta_canon_stepcurr::iaf_psc_delta_canon_stepcurr()
  : Node(), 
    P_(), 
    S_()
{}

iaf_psc_delta_canon_stepcurr::iaf_psc_delta_canon_stepcurr(const iaf_psc_delta_canon_stepcurr& n)
  : Node(n), 
    P_(n.P_), 
    S_(n.S_)
{}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void iaf_psc_delta_canon_stepcurr::init_state_(const Node& proto)
{
  const iaf_psc_delta_canon_stepcurr& pr = downcast<iaf_psc_delta_canon_stepcurr>(proto);
  S_ = pr.S_;
}

void iaf_psc_delta_canon_stepcurr::init_buffers_()
{
  B_.events_.resize();
  B_.events_.clear();
  B_.currents_.clear();        // includes resize
  B_.potentials_.clear_data(); // includes resize
}

void iaf_psc_delta_canon_stepcurr::calibrate()
{
  V_.h_ms_ = Time::get_resolution().get_ms(); 

  V_.exp_t_ = std::exp(-V_.h_ms_/P_.tau_m_);
  V_.expm1_t_ = numerics::expm1(-V_.h_ms_/P_.tau_m_);
  V_.v_inf_ = P_.I_e_ * P_.tau_m_ / P_.c_m_;
  V_.I_contrib_ = -V_.v_inf_ * V_.expm1_t_;

  // t_ref_ is the refractory period in ms
  // refractory_steps_ is the duration of the refractory period in whole
  // steps, rounded down
  V_.refractory_steps_ = Time(Time::ms(P_.t_ref_)).get_steps();
  assert(V_.refractory_steps_ >= 0);  // since t_ref_ >= 0, this can only fail in error
}

void iaf_psc_delta_canon_stepcurr::update(Time const & origin, 
				 const long_t from, const long_t to)
{
  assert(to >= 0);
  assert(static_cast<delay>(from) < Scheduler::get_min_delay());
  assert(from < to);

  // at start of slice, tell input queue to prepare for delivery
  if ( from == 0 )
    B_.events_.prepare_delivery();

  /*
    The psc_delta neuron can fire only 
    1. precisely upon spike arrival
    2. in between spike arrivals when threshold is reached due
       to the DC current.
    3. if the membrane potential is superthreshold at the BEGINNING
       of a slice due to initialization.
    In case 1, the spike time is known immediately.
    In case 2, the spike time can be found by solving the membrane
    equation.
    In case 3, the spike time is defined to be at from+epsilon.
    Thus, we can take arbitrary time steps within (from, to]. 
    Since slice_ring_buffer's delivery mechanism is built on time slices,
    we still need to step through the individual time steps to check
    for events.
    In typical network simulations and for typical step sizes, the probability
    of steps without input spikes is small. So the outer loop is over steps.
   */
  
  // check for super-threshold at beginning
  if ( S_.U_ >= P_.U_th_ )
    emit_instant_spike_(origin, from, 
			V_.h_ms_*(1-std::numeric_limits<double_t>::epsilon()));

  for ( long_t lag = from ; lag < to ; ++lag )
  {
    // time at start of update step
    const long_t T = origin.get_steps() + lag;

    // Time within step is measured by offsets, which are h at the beginning
    // and 0 at the end of the step.
    double_t t = V_.h_ms_;
    
    // place pseudo-event in queue to mark end of refractory period
    if (S_.is_refractory_ && (T + 1 - S_.last_spike_step_ == V_.refractory_steps_))
      B_.events_.add_refractory(T, S_.last_spike_offset_);

    // get first event
    double_t ev_offset;
    double_t ev_weight;
    bool     end_of_refract;

    if ( !B_.events_.get_next_spike(T, ev_offset, ev_weight, end_of_refract) )
    { // No incoming spikes, handle with fixed propagator matrix.
      // Handling this case separately improves performance significantly
      // if there are many steps without input spikes.

      // update membrane potential
      if ( !S_.is_refractory_ )
      {
	/* The following way of updating U_ is numerically more precise
	   than the more natural approach

	      U_ = exp_t_ * U_ + I_contrib_;

           particularly when U_ * exp_t_ is close to -I_contrib_.
	*/

	/*
	  this is the infinite time value of v, if the stepwise constant current
	  would go on forever
	*/
	double_t I_contrib_t = - S_.it_ * P_.tau_m_ / P_.c_m_ * V_.expm1_t_;
 
	S_.U_ = V_.I_contrib_ + I_contrib_t + V_.expm1_t_ * S_.U_ + S_.U_;

	S_.U_ = S_.U_ < P_.U_min_ ? P_.U_min_ : S_.U_;  // lower bound on potential
	if ( S_.U_ >= P_.U_th_ )
	  emit_spike_(origin, lag, 0);  // offset is zero at end of step

	// We exploit here that the refractory period must be at least
	// one time step long. So even if the spike had happened at the
	// very beginning of the step, the neuron would remain refractory
        // for the rest of the step. We can thus ignore the time reset
	// issued by emit_spike_().
      }
      // nothing to do if neuron is refractory
    }
    else
    {
      // We only get here if there is at least one event, 
      // which has been read above.  We can therefore use 
      // a do-while loop.
 
      do {

	if ( S_.is_refractory_ )
	{
	  // move time to time of event
	  t = ev_offset;

	  // normal spikes need to be accumulated
	  if ( !end_of_refract )
	  {
	    if ( S_.with_refr_input_ )
	      V_.refr_spikes_buffer_ += ev_weight 
		* std::exp( - (   ( S_.last_spike_step_ - T - 1 ) * V_.h_ms_ 
				  - ( S_.last_spike_offset_ - ev_offset )
				  + P_.t_ref_ 
			      )
			    / P_.tau_m_
		  );
	  }
	  else 
	  {
	    // return from refractoriness---apply buffered spikes
	    S_.is_refractory_ = false;

	    if ( S_.with_refr_input_ )
	    {
	      S_.U_ += V_.refr_spikes_buffer_;
	      V_.refr_spikes_buffer_ = 0.0;
	    }

	    // check if buffered spikes cause new spike
	    if ( S_.U_ >= P_.U_th_ )
	      emit_instant_spike_(origin, lag, t);
	  }

	  // nothing more to do in this loop iteration
	  continue;
	}

	// we get here only if the neuron is not refractory
	// update neuron to time of event
	// time is measured backward: inverse order in difference
	propagate_(t-ev_offset);
	t = ev_offset;

	// Check whether we have passed the threshold. If yes, emit a
	// spike at the precise location of the crossing.
        // Time within the step need not be reset to the precise time
	// of the spike, since the neuron will be refractory for the
	// remainder of the step.
	// The event cannot be a return-from-refractoriness event,
	// since that violates the assumption that the neuron was not
	// refractory. We can thus ignore the event, since the neuron
        // is refractory after the spike and ignores all input.
	if ( S_.U_ >= P_.U_th_ ) 
	{
	  emit_spike_(origin, lag, t);
	  continue;
	}

	// neuron is not refractory: add input spike, check for output
	// spike
	S_.U_ += ev_weight;
	if ( S_.U_ >= P_.U_th_ )
	  emit_instant_spike_(origin, lag, t);

      } while ( B_.events_.get_next_spike(T, ev_offset, ev_weight, end_of_refract) );

      // no events remaining, plain update step across remainder 
      // of interval
      if ( !S_.is_refractory_ && t > 0 )  // not at end of step, do remainder
      {
	propagate_(t);
	if ( S_.U_ >= P_.U_th_ ) 
	  emit_spike_(origin, lag, 0);
      }

    }  // else

    // set new input current
    S_.it_ = B_.currents_.get_value(lag);

    // voltage logging
    B_.potentials_.record_data(origin.get_steps()+lag, S_.U_ + P_.E_L_);
  }  

}                           

void iaf_psc_delta_canon_stepcurr::propagate_(const double_t dt)
{
  assert(!S_.is_refractory_);  // should not be called if neuron is
			    // refractory

  // see comment on regular update above
  const double expm1_dt = numerics::expm1(-dt/P_.tau_m_);
  const double_t v_inf = V_.v_inf_ + S_.it_ * P_.tau_m_ / P_.c_m_;
  S_.U_ = - v_inf * expm1_dt + S_.U_ * expm1_dt + S_.U_;

  return;
}
                     
void iaf_psc_delta_canon_stepcurr::emit_spike_(Time const &origin, const long_t lag,
					    const double_t offset_U) 
{
  assert(S_.U_ >= P_.U_th_);  // ensure we are superthreshold
 
  // compute time since threhold crossing
  double_t v_inf = V_.v_inf_ + S_.it_ * P_.tau_m_ / P_.c_m_;
  double_t dt = - P_.tau_m_ * std::log((v_inf - S_.U_) / (v_inf - P_.U_th_));

  // set stamp and offset for spike
  set_spiketime(Time::step(origin.get_steps() + lag + 1));
  S_.last_spike_offset_ = offset_U + dt;

  // reset neuron and make it refractory
  S_.U_ = P_.U_reset_;
  S_.is_refractory_ = true; 

  // send spike
  SpikeEvent se;
  se.set_offset(S_.last_spike_offset_);
  network()->send(*this, se, lag);

  return;
}

void iaf_psc_delta_canon_stepcurr::emit_instant_spike_(Time const &origin, const long_t lag,
						    const double_t spike_offs) 
{
  assert(S_.U_ >= P_.U_th_);  // ensure we are superthreshold
 
  // set stamp and offset for spike
  set_spiketime(Time::step(origin.get_steps() + lag + 1));
  S_.last_spike_offset_ = spike_offs;

  // reset neuron and make it refractory
  S_.U_ = P_.U_reset_;
  S_.is_refractory_ = true; 

  // send spike
  SpikeEvent se;
  se.set_offset(S_.last_spike_offset_);
  network()->send(*this, se, lag);

  return;
}

void iaf_psc_delta_canon_stepcurr::handle(SpikeEvent & e)
{
  assert(e.get_delay() > 0);

  /* We need to compute the absolute time stamp of the delivery time
     of the spike, since spikes might spend longer than min_delay_
     in the queue.  The time is computed according to Time Memo, Rule 3.
  */
  const long_t Tdeliver = e.get_stamp().get_steps() + e.get_delay() - 1;
  B_.events_.add_spike(e.get_rel_delivery_steps(network()->get_slice_origin()), 
		    Tdeliver, e.get_offset(), e.get_weight() * e.get_multiplicity());
}

void iaf_psc_delta_canon_stepcurr::handle(CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  const double_t c=e.get_current();
  const double_t w=e.get_weight();

  // add stepwise constant current; MH 2009-10-14
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
		      w*c);
}

void iaf_psc_delta_canon_stepcurr::handle(PotentialRequest& e)
{
  B_.potentials_.handle(*this, e);
}

void iaf_psc_delta_canon_stepcurr::set_spiketime(Time const & now)
{
  S_.last_spike_step_ = now.get_steps();
}
 
} // namespace

