/*
 *  theta_neuron_ps.cpp
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2009 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#include "theta_neuron_ps.h"

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include "universal_data_logger_impl.h"
#include <limits>

#include <iomanip>
#include <iostream>
#include <cstdio>

/* ---------------------------------------------------------------- 
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap<nest::theta_neuron_ps> nest::theta_neuron_ps::recordablesMap_;

namespace nest // template specialization must be placed in namespace
{
  // Override the create() method with one call to
  // RecordablesMap::insert_() for each quantity to be recorded.
  template <>
  void RecordablesMap<theta_neuron_ps>::create()
  {
    // Use standard names whereever you can for consistency!
    insert_(names::theta,
	    &theta_neuron_ps::get_theta_);
  }
}

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
nest::theta_neuron_ps::Parameters_::Parameters_()
  : I_e(-1e-4)
{
  recordablesMap_.create();
}

nest::theta_neuron_ps::State_::State_(const Parameters_ &)
  : theta(0.0)
{
}

nest::theta_neuron_ps::State_::State_(const State_ &s)
  : theta(s.theta)
{
}

nest::theta_neuron_ps::Buffers_::Buffers_(theta_neuron_ps &n)
  : logger_(n)
{
  // The other member variables are left uninitialised or are
  // automatically initialised by their default constructor.
}

nest::theta_neuron_ps::Buffers_::Buffers_(const Buffers_ &, theta_neuron_ps &n)
  : logger_(n)
{
  // The other member variables are left uninitialised or are
  // automatically initialised by their default constructor.
}

/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::theta_neuron_ps::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d, names::I_e, I_e);
}

void nest::theta_neuron_ps::Parameters_::set(const DictionaryDatum& d)
{
  updateValue<double>(d, names::I_e, I_e);
}

void nest::theta_neuron_ps::State_::get(DictionaryDatum &d) const
{
  def<double>(d, names::theta, theta);
}

void nest::theta_neuron_ps::State_::set(const DictionaryDatum& d, const Parameters_&)
{
  if ( updateValue<double>(d, names::theta, theta) &&
       ( theta > std::acos(-1.0) || theta < -std::acos(-1.0) ) )
    throw BadProperty("Initial neuron phase must be in (-pi, pi).");
}

/* ---------------------------------------------------------------- 
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

nest::theta_neuron_ps::theta_neuron_ps()
  : Archiving_Node(), 
    P_(), 
    S_(P_),
    B_(*this)
{
}

nest::theta_neuron_ps::theta_neuron_ps(const theta_neuron_ps& n)
  : Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{
}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::theta_neuron_ps::init_state_(const Node &proto)
{
  const theta_neuron_ps &pr = downcast<theta_neuron_ps>(proto);

  S_ = pr.S_;
}

void nest::theta_neuron_ps::init_buffers_()
{
  Archiving_Node::clear_history();

  B_.events_.resize();
  B_.events_.clear();
  B_.currents_.clear(); // includes resize

  B_.logger_.reset();
}

void nest::theta_neuron_ps::calibrate()
{
  B_.logger_.init();  // ensures initialization in case mm connected after Simulate
  V_.h_ms = Time::get_resolution().get_ms();
  V_.I_0  = P_.I_e;
  if ( V_.I_0 == 0.0 )
    throw BadProperty("Total external input must not be zero.");
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void nest::theta_neuron_ps::update(Time const &origin, long_t const from, long_t const to)
{
  assert ( to >= 0 && (delay) from < Scheduler::get_min_delay() );
  assert ( from < to );

  // at start of slice, tell input queue to prepare for delivery
  if ( from == 0 )
    B_.events_.prepare_delivery();

  for ( long_t lag = from; lag < to; ++lag )
  {
    // time at start of update step
    const long_t T = origin.get_steps() + lag;

    // propagation step equals simulation resolution in the absence
    // of input spikes
    double_t last_event_offset = V_.h_ms;

    // event info
    double_t ev_offset;
    double_t ev_weight;
    bool     end_of_refract;

    while ( B_.events_.get_next_spike(T, ev_offset, ev_weight, end_of_refract) )
    { // there are incoming spikes
      propagate_(origin, lag, last_event_offset, last_event_offset-ev_offset, ev_weight);
      last_event_offset = ev_offset;
    }
    propagate_(origin, lag, last_event_offset, last_event_offset, 0.0);

    // log state data
    B_.logger_.record_data(T);

    // set new input current
    V_.I_0 = P_.I_e+B_.currents_.get_value(lag);
    if ( V_.I_0 == 0.0 )
      throw BadProperty("Total external input must not be zero.");
  }
}

void nest::theta_neuron_ps::propagate_(Time const &origin, long_t const lag, double_t const last_offset, double_t const dt, double_t const w)
{
  assert ( V_.I_0 != 0.0 );

  if ( dt == 0.0 )
    return;

  assert ( dt > 0.0 );

  // theta at the left interval border
  double_t const theta_before = S_.theta;

  double_t const beta  = std::sqrt(std::abs(V_.I_0));
  double_t const c1_jk = std::tan(theta_before/2.0) / beta;

  // propagate theta (theta_minus)
  // and integrate spike input (theta_plus)
  if ( V_.I_0 < 0.0 )
  {
    double_t const c2_jk = std::tanh(-beta*dt);

    S_.theta = 2.0*std::atan(w + beta*(c1_jk + c2_jk) / (1.0 + c1_jk*c2_jk));
  }
  else
  {
    double_t const c2_jk = std::tan(beta*dt);

    S_.theta = 2.0*std::atan(w + beta*(c1_jk + c2_jk) / (1.0 - c1_jk*c2_jk));
  }

  // We do not need the exact spike time, so we apply the usual
  // retrospective spike test and compute the exact spike time in case
  // of a positive outcome. As we cannot check for theta > pi due to
  // the cyclic behaviour, we assume that a spike occured within the
  // interval if (approx.) theta > pi/2 at the left interval border
  // and theta < -pi/2 at the right one.
  if ( theta_before > 1.5 && S_.theta < -1.5 )
  {
    double_t t_spike;

    if ( V_.I_0 < 0.0 )
      t_spike = 1.0/beta * atanh(1.0/c1_jk);
    else
      t_spike = 1.0/beta * std::atan(1.0/c1_jk);

    assert ( t_spike >= 0 );

    // threshold crossing caused by the incoming spike
    if ( t_spike > dt )
       t_spike = dt;

    set_spiketime(Time::step(origin.get_steps()+lag+1));

    SpikeEvent se;

    se.set_offset(last_offset-t_spike);
    network()->send(*this, se, lag);
  }
}

void nest::theta_neuron_ps::handle(SpikeEvent &e)
{
  assert ( e.get_delay() > 0 );

  // We need to compute the absolute time stamp of the delivery time
  // of the spike, since spikes might spend longer than min_delay_
  // in the queue.  The time is computed according to Time Memo, Rule 3.
  const long_t Tdeliver = e.get_stamp().get_steps() + e.get_delay() - 1;

  B_.events_.add_spike(e.get_rel_delivery_steps(network()->get_slice_origin()), 
		       Tdeliver, e.get_offset(), e.get_weight() * e.get_multiplicity());
}

void nest::theta_neuron_ps::handle(CurrentEvent &e)
{
  assert ( e.get_delay() > 0 );

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
			 e.get_weight() * e.get_current());
}

void nest::theta_neuron_ps::handle(DataLoggingRequest &e)
{
  B_.logger_.handle(e);
}

