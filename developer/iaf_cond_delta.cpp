/*
 *  iaf_cond_delta.cpp
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2005 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

/* iaf_cond_delta is a neuron with conductance synapses where the potential jumps 
 * on each spike arrival. */

#include "iaf_cond_delta.h"

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include <limits>

#include "analog_data_logger_impl.h"
#include "event.h"

#include <iomanip>
#include <iostream>
#include <cstdio>

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
nest::iaf_cond_delta::Parameters_::Parameters_()
  : E_L_       (-70.0     ),  // mV
    V_th_      (-55.0-E_L_),  // mV, relative to E_L
    V_reset_   (-70.0-E_L_),  // mV, relative to E_L
    E_ex_      (  0.0-E_L_),  // mV, relative to E_L
    E_in_      (-80.0-E_L_),  // mV, relative to E_L
    t_ref_     (  2.0     ),  // ms
    g_L_       ( 16.6667  ),  // nS, note, tau_m=C_m/g_L
    C_m_       (250.0     )  // pF
{}

nest::iaf_cond_delta::State_::State_()
  : V_m_   (0.0),
    r_    (0)
{}

/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::iaf_cond_delta::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d,names::E_L,          E_L_); 

  def<double>(d,names::V_th,         V_th_+E_L_);
  def<double>(d,names::V_reset,      V_reset_+E_L_);
  def<double>(d,names::E_ex,         E_ex_+E_L_);
  def<double>(d,names::E_in,         E_in_+E_L_);

  def<double>(d,names::t_ref,        t_ref_);
  def<double>(d,names::g_L,          g_L_);
  def<double>(d,names::C_m,          C_m_);
}

void nest::iaf_cond_delta::Parameters_::set(const DictionaryDatum& d)
{
  // allow setting the membrane potential
  updateValue<double>(d,names::E_L,     E_L_);

  if (updateValue<double>(d, names::V_th, V_th_)) 
    V_th_ -= E_L_;
  if(updateValue<double>(d, names::V_reset, V_reset_))
    V_reset_ -= E_L_;
  if (updateValue<double>(d, names::E_ex, E_ex_))
    E_ex_ -= E_L_;
  if (updateValue<double>(d, names::E_in, E_in_))
    E_in_ -= E_L_;
    
  updateValue<double>(d,names::t_ref,   t_ref_); 
  updateValue<double>(d,names::g_L,     g_L_);
  updateValue<double>(d,names::C_m,     C_m_);

  if ( V_reset_ >= V_th_ )
    throw BadProperty("Reset potential must be smaller than threshold.");   
  if ( C_m_ <= 0 )
    throw BadProperty("Capacitance must be strictly positive.");    
  if ( g_L_ <= 0 )
    throw BadProperty("Leak conductance must be strictly positive.");    
  if ( t_ref_ < 0 )
    throw BadProperty("Refractory time cannot be negative.");
      
}

void nest::iaf_cond_delta::State_::get(DictionaryDatum &d, const Parameters_& p) const
{
  def<double>(d, names::V_m, V_m_ + p.E_L_); // Membrane potential
}

void nest::iaf_cond_delta::State_::set(const DictionaryDatum& d, const Parameters_& p)
{
  if (updateValue<double>(d, names::V_m, V_m_))
    V_m_ -= p.E_L_;
}

/* ---------------------------------------------------------------- 
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

nest::iaf_cond_delta::iaf_cond_delta()
  : Archiving_Node(), 
    P_(), 
    S_()
{}

nest::iaf_cond_delta::iaf_cond_delta(const iaf_cond_delta& n)
  : Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_)
{}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::iaf_cond_delta::init_state_(const Node& proto)
{
  const iaf_cond_delta& pr = downcast<iaf_cond_delta>(proto);
  S_ = pr.S_;
}

void nest::iaf_cond_delta::init_buffers_()
{
  B_.spikes_ex_.clear();       // includes resize
  B_.spikes_in_.clear();       // includes resize
  B_.potentials_.clear_data(); // includes resize
  Archiving_Node::clear_history();
}

void nest::iaf_cond_delta::calibrate()
{
  const double h = Time::get_resolution().get_ms(); 
 
  V_.P_L_ = std::exp(-h*P_.g_L_/P_.C_m_);

  V_.RefractoryCounts_ = Time(Time::ms(P_.t_ref_)).get_steps();
  assert(V_.RefractoryCounts_ >= 0);  // since t_ref_ >= 0, this can only fail in error
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void nest::iaf_cond_delta::update(Time const & origin, const long_t from, const long_t to)
{
   
  assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
  assert(from < to);

  for ( long_t lag = from ; lag < to ; ++lag )
  {
    
    if ( S_.r_ == 0 )
    {
      // neuron not refractory
      S_.V_m_ = V_.P_L_*S_.V_m_;
      S_.V_m_ += (P_.E_ex_ - S_.V_m_)*B_.spikes_ex_.get_value(lag) + (P_.E_in_ - S_.V_m_)*B_.spikes_in_.get_value(lag);

    }
    else // neuron is absolute refractory
    {
      B_.spikes_ex_.get_value(lag);  // clear buffer entry, ignore spike
      B_.spikes_in_.get_value(lag);  // clear buffer entry, ignore spike
      --S_.r_;
    }
   
    // threshold crossing
    if (S_.V_m_ >= P_.V_th_)
    {
      S_.r_ = V_.RefractoryCounts_;
      S_.V_m_ = P_.V_reset_;
        
      // EX: must compute spike time
      set_spiketime(Time::step(origin.get_steps()+lag+1));
        
      SpikeEvent se;
      network()->send(*this, se, lag);
    }

    // voltage logging
    B_.potentials_.record_data(origin.get_steps()+lag, S_.V_m_ + P_.E_L_);
  }
}

void nest::iaf_cond_delta::handle(SpikeEvent & e)
{
  assert(e.get_delay() > 0);

  if(e.get_weight() > 0.0)
    B_.spikes_ex_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			 e.get_weight() * e.get_multiplicity() );
  else
    B_.spikes_in_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			 -e.get_weight() * e.get_multiplicity() );  // note: conductance is positive for both exc and inh synapses but we use the sign of spike event weight to distinguish which kind of synapse it is
}

void nest::iaf_cond_delta::handle(PotentialRequest& e)
{
  B_.potentials_.handle(*this, e);
}

