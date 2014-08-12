/*
 *  theta_neuron.cpp
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

#include "theta_neuron.h"

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

nest::RecordablesMap<nest::theta_neuron> nest::theta_neuron::recordablesMap_;

namespace nest // template specialization must be placed in namespace
{
  // Override the create() method with one call to
  // RecordablesMap::insert_() for each quantity to be recorded.
  template <>
  void RecordablesMap<theta_neuron>::create()
  {
    // Use standard names whereever you can for consistency!
    insert_(names::theta,
	    &theta_neuron::get_theta_);
  }
}

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
nest::theta_neuron::Parameters_::Parameters_()
  : I_e(-1e-4)
{
  recordablesMap_.create();
}

nest::theta_neuron::State_::State_(const Parameters_ &)
  : theta(0.0)
{
}

nest::theta_neuron::State_::State_(const State_ &s)
  : theta(s.theta)
{
}

nest::theta_neuron::State_& nest::theta_neuron::State_::operator=(const State_ &s)
{
  if ( this == &s ) // avoid assignment to self
    return *this;

  theta = s.theta;
  return *this;
}

nest::theta_neuron::Buffers_::Buffers_(theta_neuron &n)
  : logger_(n)
{
  // The other member variables are left uninitialised or are
  // automatically initialised by their default constructor.
}

nest::theta_neuron::Buffers_::Buffers_(const Buffers_ &, theta_neuron &n)
  : logger_(n)
{
  // The other member variables are left uninitialised or are
  // automatically initialised by their default constructor.
}

/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::theta_neuron::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d, names::I_e, I_e);
}

void nest::theta_neuron::Parameters_::set(const DictionaryDatum &d)
{
  updateValue<double>(d, names::I_e, I_e);
}

void nest::theta_neuron::State_::get(DictionaryDatum &d) const
{
  def<double>(d, names::theta, theta);
}

void nest::theta_neuron::State_::set(const DictionaryDatum &d, const Parameters_ &)
{
  if ( updateValue<double>(d, names::theta, theta) &&
       ( theta > std::acos(-1.0) || theta < -std::acos(-1.0) ) )
    throw BadProperty("Initial neuron phase must be in (-pi, pi).");
}

/* ---------------------------------------------------------------- 
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

nest::theta_neuron::theta_neuron()
  : Archiving_Node(), 
    P_(), 
    S_(P_),
    B_(*this)
{
}

nest::theta_neuron::theta_neuron(const theta_neuron &n)
  : Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{
}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::theta_neuron::init_state_(const Node &proto)
{
  const theta_neuron &pr = downcast<theta_neuron>(proto);

  S_ = pr.S_;
}

void nest::theta_neuron::init_buffers_()
{
  Archiving_Node::clear_history();

  B_.spikes_.clear();   // includes resize
  B_.currents_.clear(); // includes resize

  B_.logger_.reset();
}

void nest::theta_neuron::calibrate()
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

void nest::theta_neuron::update(Time const &origin, const long_t from, const long_t to)
{
  assert ( to >= 0 && (delay) from < Scheduler::get_min_delay() );
  assert ( from < to );

  assert ( V_.I_0 != 0.0 );

  for ( long_t lag = from; lag < to; ++lag )
  {
    // theta at the left interval border
    double_t const theta_before = S_.theta;

    double_t const beta  = std::sqrt(std::abs(V_.I_0));
    double_t const c1_jk = std::tan(theta_before/2.0) / beta;

    // spike input
    double_t const w = B_.spikes_.get_value(lag);

    // propagate theta (theta_minus)
    // and integrate spike input (theta_plus)
    if ( V_.I_0 < 0.0 )
    {
      double_t const c2_jk = std::tanh(-beta*V_.h_ms);

      S_.theta = 2.0*std::atan(w + beta*(c1_jk + c2_jk) / (1.0 + c1_jk*c2_jk));
    }
    else
    {
      double_t const c2_jk = std::tan(beta*V_.h_ms);

      S_.theta = 2.0*std::atan(w + beta*(c1_jk + c2_jk) / (1.0 - c1_jk*c2_jk));
    }

    // We do not need the exact spike time, so we apply the usual
    // retrospective spike test. As we cannot check for theta > pi due
    // to the cyclic behaviour, we assume that a spike occured within
    // the interval if (approx.) theta > pi/2 at the left interval
    // border and theta < -pi/2 at the right one.
    if ( theta_before > 1.5 && S_.theta < -1.5 )
    {
      set_spiketime(Time::step(origin.get_steps()+lag+1));

      SpikeEvent se;

      network()->send(*this, se, lag);
    }

    // log state data
    B_.logger_.record_data(origin.get_steps()+lag);

    // set new input current
    V_.I_0 = P_.I_e+B_.currents_.get_value(lag);
    if ( V_.I_0 == 0.0 )
      throw BadProperty("Total external input must not be zero.");
  }
}

void nest::theta_neuron::handle(SpikeEvent &e)
{
  assert ( e.get_delay() > 0 );

  B_.spikes_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
		       e.get_weight() * e.get_multiplicity() );
}

void nest::theta_neuron::handle(CurrentEvent &e)
{
  assert ( e.get_delay() > 0 );

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
			 e.get_weight() * e.get_current());
}

void nest::theta_neuron::handle(DataLoggingRequest &e)
{
  B_.logger_.handle(e);
}

