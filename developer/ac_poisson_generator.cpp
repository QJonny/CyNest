/*
 *  ac_poisson_generator.cpp
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

#include "exceptions.h"
#include "ac_poisson_generator.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "arraydatum.h"
#include "dictutils.h"
#include "numerics.h"
#include "analog_data_logger_impl.h"

#include <cmath>
#include <limits>

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameter
 * ---------------------------------------------------------------- */
    
nest::ac_poisson_generator::Parameters_::Parameters_()
  : om_(),
    phi_(),
    dc_(0),
    ac_()
{}

/*---------------------------
 * Copy constructor EN
 *--------------------------*/
nest::ac_poisson_generator::Parameters_::Parameters_(const Parameters_& p )
  : om_(p.om_),
    phi_(p.phi_),
    dc_(p.dc_),
    ac_(p.ac_)
{}

nest::ac_poisson_generator::Parameters_& 
nest::ac_poisson_generator::Parameters_::operator=(const Parameters_& p)
{
  dc_ = p.dc_;

  om_.resize(p.om_.size());
  phi_.resize(p.phi_.size());
  ac_.resize(p.ac_.size());

  for(nest::int_t i = 0; i < om_.size(); ++i)
    {
      om_[i] = p.om_[i];
    }

  for(nest::int_t i = 0; i < phi_.size(); ++i)
    {
      phi_[i] = p.phi_[i];
    }

  for(nest::int_t i = 0; i < ac_.size(); ++i)
    {
      ac_[i] = p.ac_[i];
    }

  return *this;
}

nest::ac_poisson_generator::State_::State_()
  : last_spike_(Time::step(-1))
{}


// The following initialisations existed in the old
// version of ac_poisson_generator, but are not 
// initialised in the new version following the 
// model scheme introduced in r7458
//   N_osc_(om_.size()),  
//   poisson_dev_(dc_)


/* ---------------------------------------------------------------- 
 * Parameter extraction and manipulation functions
 * ---------------------------------------------------------------- */

void nest::ac_poisson_generator::Parameters_::get(DictionaryDatum &d) const
{
  (*d)["DC"] = dc_ * 1000.0; 

  //  std::valarray<double_t> freq();
  store_array_(d, "Freq", om_ / ( 2.0 * numerics::pi / 1000.0)); // in Hz
  store_array_(d, "Phi", phi_);
  store_array_(d, "AC", ac_ * 1000.0);  // convert to s^-1

  // warn about mismatch between oscillator arrays
  if ( ! ( om_.size() == phi_.size() && om_.size() == ac_.size() ) )
    network()->message(SLIInterpreter::M_INFO,
		       "ac_poisson_generator::get_status",
	      "Amplitude, frequency and phase arrays have different sizes.\n"
	      "Ensure they have equal numbers of elements before starting\n"
              "a simulation!");
}  

void nest::ac_poisson_generator::State_::get(DictionaryDatum &d) const
{
  (*d)["last_spike"] = last_spike_.get_ms();
}  

void nest::ac_poisson_generator::Parameters_::set(const DictionaryDatum& d)
{
  if ( updateValue<double>(d, "DC", dc_) )
    dc_ /= 1000.0;           // scale to ms^-1
  
  if ( extract_array_(d, "Freq", om_) )
    om_ *= 2.0 * numerics::pi / 1000.0;

  extract_array_(d, "Phi", phi_);

  if ( extract_array_(d, "AC", ac_) )
    {
      ac_ /= 1000.0;

    }

  // warn about mismatch between oscillator arrays
  if ( ! ( om_.size() == phi_.size() && om_.size() == ac_.size() ) )
    network()->message(SLIInterpreter::M_WARNING,
		       "ac_poisson_generator::get_status",
	      "Amplitude, frequency and phase arrays have different sizes.\n"
	      "Ensure they have equal numbers of elements before starting\n"
              "a simulation!");
}

bool nest::ac_poisson_generator::
Parameters_::extract_array_(const DictionaryDatum &d, 
			    const std::string& dname,
			    std::valarray<double>& data) const
{
  if ( d->known(dname) )
  {
    ArrayDatum *ad = dynamic_cast<ArrayDatum *>((*d)[dname].datum());
    if ( ad == 0 )
      throw BadProperty();

    const size_t nd = ad->size();
    data.resize(nd);
    for ( size_t n = 0 ; n < nd ; ++n )
      {
	data[n] = getValue<double>((*ad)[n]);
      }

    return true;
  }
  else
    return false;
}


void nest::ac_poisson_generator::
Parameters_::store_array_(DictionaryDatum &d, 
			  const std::string& dname,
			  std::valarray<double> data) const
{
  // copy data to vector; valarrays have no iterators, se we must loop
  std::vector<double> tmp(data.size());
  for ( size_t n = 0 ; n < data.size() ; ++n )
    {
    tmp[n] = data[n];
    }

  // store in dictionary
  (*d)[dname] = new ArrayDatum(tmp);

  return;
} 

/* ---------------------------------------------------------------- 
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::ac_poisson_generator::ac_poisson_generator()
  : Node(),
    device_(),
    P_(),
    S_(),
    B_()
{
}

nest::ac_poisson_generator::ac_poisson_generator(const ac_poisson_generator&n)
  : Node(n),
    device_(n.device_),
    P_(n.P_),
    S_(n.S_),
    B_(n.B_)
{
}

// ----------------------------------------------------

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::ac_poisson_generator::init_state_(const Node& proto)
{ 
  const ac_poisson_generator& pr = 
    downcast<ac_poisson_generator>(proto);
  
  device_.init_state(pr.device_);
  S_ = pr.S_;
}

void nest::ac_poisson_generator::init_buffers_()
{ 
  device_.init_buffers();
  B_.rates_.clear_data();

  B_.N_osc_ = P_.om_.size();
  if ( ! ( P_.ac_.size() == B_.N_osc_ && P_.phi_.size() == B_.N_osc_ ) )
  {
    network()->message(SLIInterpreter::M_ERROR,
		       "ac_poisson_generator::calibrate",
		       "Frequency, AC amplitude, and phase arrays "
		       "must have equal size!");
    throw DimensionMismatch();
  }

  /* state vector:
     - dc component is constant and not included
     - rate is accumulated directly in update
     - vector contains two elements for each oscillator
  */
  B_.state_oscillators_.resize(2*B_.N_osc_, 0);

  // remaining elements, see Rotter & Diesmann, 3.1.3
  // phi = pi/2 is stable, since tan(phi)^2 appears in denominator only
  // assumes IEEE arithmetic with 1/Inf = 0
  // sign is handled below when assigning to S_
  const std::valarray<double_t> q = 
    P_.ac_ / std::sqrt(1.0 + std::pow(std::tan(P_.phi_), 2.0));

  // valarray slices do not appear to be reliable under Tru64 CXX,
  // so one better uses for loops.
  for ( size_t n = 0 ; n < B_.N_osc_ ; ++n )
    {
    B_.state_oscillators_[2 * n    ] = 
      std::cos(P_.phi_[n]) >= 0 ? q[n] : -q[n];
    B_.state_oscillators_[2 * n + 1] = std::sin(P_.phi_[n]) >= 0 
      ?  std::sqrt(std::pow(P_.ac_[n],2.0)-std::pow(q[n],2.0))
      : -std::sqrt(std::pow(P_.ac_[n],2.0)-std::pow(q[n],2.0));
  }

}

void nest::ac_poisson_generator::calibrate()
{
  // This check is placed here so that it only gets executed on actual
  // instances, not on Model instances during a ResetKernel.
  // TODO: Make work for parallel simulation!
  if ( network()->get_num_threads() > 1 || network()->get_num_processes() > 1 )
    throw BadProperty("ac_poisson_generator is presently not suitable"
                      " for parallel simulation.");

  device_.calibrate();

  // time resolution
  const double h = Time::get_resolution().get_ms(); 

  V_.sins_.resize(B_.N_osc_);
  V_.coss_.resize(B_.N_osc_);
  V_.sins_ = std::sin(h * P_.om_);         // block elements
  V_.coss_ = std::cos(h * P_.om_);

  B_.rates_.calibrate();

  return;
}


void nest::ac_poisson_generator::update(Time const& origin, 
			 const long_t from, const long_t to)
{
  assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
  assert(from < to);

  const long_t start = origin.get_steps();

  // time resolution
  const double h = Time::get_resolution().get_ms(); 

  // random number generator
  librandom::RngPtr rng = net_->get_rng(get_thread());

  // We iterate the dynamics even when the device is turned off,
  // but do not issue spikes while it is off. In this way, the 
  // oscillators always have the right phase.  This is quite 
  // time-consuming, so it should be done only if the device is
  // on most of the time.

  for ( long_t lag = from ; lag < to ; ++lag )
  {
    // update oscillator blocks, accumulate rate as sum of DC and N_osc_ AC elements
    // rate is instantaneous sum of state
    double r = P_.dc_;

    for ( unsigned int n = 0 ; n < B_.N_osc_ ; ++n )
    {
      const unsigned int offs = 2 * n;  // index of first block elem 
      const double_t new1 = V_.coss_[n] * B_.state_oscillators_[offs] - V_.sins_[n] * B_.state_oscillators_[offs+1];
    
      B_.state_oscillators_[offs+1] = 
	V_.sins_[n] * B_.state_oscillators_[offs] + 
	V_.coss_[n] * B_.state_oscillators_[offs+1];
      B_.state_oscillators_[offs]    = new1;

      r += B_.state_oscillators_[offs+1];
    }

    // store rate in Hz
    B_.rates_.record_data(origin.get_steps()+lag, r * 1000);

    // create spikes
    if ( r > 0 && device_.is_active(Time::step(start+lag)) )
    {
      V_.poisson_dev_.set_lambda(r * h);
      ulong_t n_spikes = V_.poisson_dev_.uldev(rng);
      while ( n_spikes-- )
      {
      	SpikeEvent se;
        network()->send(*this, se, lag);
        S_.last_spike_ = Time::step(origin.get_steps()+lag+1);
      }
    }

  }
}                       

void nest::ac_poisson_generator::handle(PotentialRequest& e)
{
  B_.rates_.handle(*this, e);
}




