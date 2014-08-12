/*
 *  inh_gamma_generator.cpp
 *
 *  Based loosely on the poisson_generator
 *
 *  AUTHOR: Eilif Muller
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

#include "inh_gamma_generator.h"

#ifdef HAVE_GSL

#include "network.h"
#include "dict.h"
#include "doubledatum.h"
#include "integerdatum.h"
#include "arraydatum.h"
#include "dictutils.h"
#include "exceptions.h"
#include <gsl/gsl_sf.h>


namespace nest {

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

nest::inh_gamma_generator::Parameters_::Parameters_()
  : a_(), b_(), tbins_()

{
  tbins_.clear();
  tbins_.push_back(0.0);
  a_.clear();
  a_.push_back(1.0);
  b_.clear();
  b_.push_back(1.0);

}

nest::inh_gamma_generator::Parameters_::Parameters_(const Parameters_& op)
  : a_(op.a_), b_(op.b_), tbins_(op.tbins_)
{

  //synapse_context_ = n.synapse_context_;

  /*
  // clear tspike_
  size_t i;

  for (i=0;i<tspike_.size();i++) {
    tspike_[i]=0.0;
    //tlast_[i]=-1000.0;
    tlast_[i]=0.0;
  }

  bin_position_ = 0;
  */

}

nest::inh_gamma_generator::State_::State_()
  : tspike_(), tlast_(), bin_position_(0)
{

  /*
  rmax_ = 1.0;

  tspike_.clear();
  tlast_.clear();
  bin_position_ = 0;
  */

}

nest::inh_gamma_generator::State_::State_(const State_& op)
  : tspike_(op.tspike_), tlast_(op.tlast_), bin_position_(op.bin_position_)
{

  /*
  rmax_ = 1.0;

  tspike_.clear();
  tlast_.clear();
  bin_position_ = 0;
  */

}



/* ---------------------------------------------------------------- 
 * Parameter extraction and manipulation functions
 * ---------------------------------------------------------------- */


void nest::inh_gamma_generator::Parameters_::get(DictionaryDatum &d) const
{
  // we need to convert to vector of doubles first
  //std::vector<double_t> d_times(spike_times_.size());
  //for (size_t i = 0; i < d_times.size(); ++i)
  //  d_times[i] = spike_times_[i].get_ms();

  (*d)["tbins"] = DoubleVectorDatum(new std::vector<double_t>(tbins_));
  (*d)["a"] = DoubleVectorDatum(new std::vector<double_t>(a_));
  (*d)["b"] = DoubleVectorDatum(new std::vector<double_t>(b_));

}

void nest::inh_gamma_generator::State_::get(DictionaryDatum &d) const
{

  (*d)["tspike"] = DoubleVectorDatum(new std::vector<double_t>(tspike_));
  (*d)["tlast"] = DoubleVectorDatum(new std::vector<double_t>(tlast_));

}


 
void nest::inh_gamma_generator::Parameters_::set(const DictionaryDatum& d, State_& s)
{


  // start/stop properties

  Name tbins_name("tbins");

  if(d->known(tbins_name))
    {
      tbins_ = getValue<std::vector<double_t> >(d->lookup(tbins_name));

      // ensure that tbins are sorted
      if ( tbins_.size()>1 )
	{ 
	  std::vector<double_t>::const_iterator prev = tbins_.begin();
	  for ( std::vector<double_t>::const_iterator next = prev + 1;
		next != tbins_.end() ; ++next, ++prev )
	    if ( *prev > *next )
	      // parent call get_status handles revert in case of error
	      throw BadProperty("tbins must be sorted in non-descending order.");
	}


    }
  
  /*
  Name tspike_name("tspike");

  if(d->known(tspike_name))
    {
      tspike_ = getValue<std::vector<double_t> >(d->lookup(tspike_name));
    }

  Name tlast_name("tlast");

  if(d->known(tlast_name))
    {
      tlast_ = getValue<std::vector<double_t> >(d->lookup(tlast_name));
    }
  */

  Name a_name("a");

  if(d->known(a_name))
    {
      a_ = getValue<std::vector<double_t> >(d->lookup(a_name));
    }

  Name b_name("b");

  if(d->known(b_name))
    {
      b_ = getValue<std::vector<double_t> >(d->lookup(b_name));

      // setting rmax with calibrate in parent call "set_status" ...
    }

  if (b_.size()!=a_.size() || b_.size()!=tbins_.size()) {
    // parent call get_status handles revert in case of error
      throw BadProperty("Size of a, b, tbins must be the same.");
  }
  

}

/* ---------------------------------------------------------------- 
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::inh_gamma_generator::inh_gamma_generator()
  : Node(), device_(),P_(),S_(), exprand_dev_(), num_connections_(0), rmax_(0.0)
{

  
  calibrate();

  //synapse_context_ = -1;

}

nest::inh_gamma_generator::inh_gamma_generator(inh_gamma_generator const &n)
  : Node(n), device_(n.device_),
    P_(n.P_), S_(n.S_), exprand_dev_(n.exprand_dev_),
    num_connections_(n.num_connections_), rmax_(n.rmax_)
{
  calibrate();

}

void nest::inh_gamma_generator::init_state_(const Node& proto)
{ 
  const inh_gamma_generator& pr = downcast<inh_gamma_generator>(proto);

  device_.init_state(pr.device_);
  S_ = pr.S_;

  //network()->message(SLIInterpreter::M_WARNING, "inh_gamma_generator::init_state_", 
  //"Dynamic state initialization not fully implemented, be careful!");
}


void nest::inh_gamma_generator::init_buffers_()
{ 
  device_.init_buffers();

}


void nest::inh_gamma_generator::calibrate()
{

  // Auto calculation of S_.rmax_ from P_.b_

  // rmax is upper bound of hazard, and rate of surrogate Poisson process
  // for gamma for a given b, it converges to monotonically increasing to 1/b
  // so we must loop over b and find max of 1/b  -> that is rmax

  rmax_ = 0.0;
  for (size_t i=0; i<P_.b_.size(); i++) {
    if (1.0/P_.b_[i]>rmax_) rmax_ = 1.0/P_.b_[i];
  }
  rmax_ *= 1.02;  // 2% margin for safety.

  device_.calibrate();
}

//
// Time Evolution Operator
//
void nest::inh_gamma_generator::update(Time const & sliceT0, const long_t from, const long_t to)
{

  assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
  assert(from < to);

  if ( rmax_ <= 0.0 ) {
    // no spikes to be generated
    return;
  }

  const Time tstop  = sliceT0 + Time::step(to);

  Time spike_time;

  double_t t_bin_stop;
  double_t dt = Time::get_resolution().get_ms();
  //Time time;
  librandom::RngPtr rng = net_->get_rng(get_thread());
  double c;
  double r;

  // need to loop over because tbins_ specifies
  // resolution on which a and b should vary
  // (inhomogeneous process)

  // loop through tbins_ histogram until tstop
  // (bin_position_++ is at bottom of loop)
  while (P_.tbins_[S_.bin_position_] <= tstop.get_ms() ) {
    
    // get right side of tbin if there is one and <tstop, or get tstop
    if (S_.bin_position_+1 < P_.tbins_.size() && P_.tbins_[S_.bin_position_+1]<tstop.get_ms() )
      t_bin_stop = P_.tbins_[S_.bin_position_+1];
    else 
      t_bin_stop = tstop.get_ms();
    
    // loop over all realizations
    for (size_t p=0;p<S_.tspike_.size();p++) {

      // spike of surrogate process still in the presen tbin
      while (S_.tspike_[p]<=t_bin_stop) {

	// calc the gamma hazard
	// b is in seconds so convert to milliseconds
	calc_h(&c,S_.tspike_[p]-S_.tlast_[p],P_.a_[S_.bin_position_],P_.b_[S_.bin_position_]*1000.0,
	       dt/10.0);

	// now get rate in Hz
	c*=1000.0;

	// check if we should really emit a spike (thinning step)

	// uniform random # on [0,rmax_)
	r = rng->drand()*rmax_;

	if (r<c) {
	// spike!
	S_.tlast_[p] = S_.tspike_[p];

	spike_time = Time(Time::ms( S_.tspike_[p]));

	// send spike if device is active at spike time
	if (device_.is_active(spike_time)) {
	    SpikeEvent se;

	    se.set_offset(Time(Time::step(spike_time.get_steps())).get_ms() - S_.tspike_[p]);
	    // we need to subtract one from stamp which is added again in send()
	    long_t lag = Time(spike_time - sliceT0).get_steps() - 1;

	    se.set_port(p);
	    network()->send(*this, se, lag);
	}

	}

	// get next "possible" spike time
	S_.tspike_[p] += (1000.0/rmax_)*exprand_dev_(rng);  //rmax_ is in Hz
	//tspike_[p] += exprand_dev_(rng)*10.0;  //rmax_ is in Hz
	//tspike_[p] += 1.0;


      } // while tspike_[p]<time



    } // for p


    // go to next bin position if its before tstop
    if (S_.bin_position_+1 < P_.tbins_.size() && P_.tbins_[S_.bin_position_+1]<tstop.get_ms() )
      S_.bin_position_++;
    else
      // no other bin position, then we are done.
      break;

  } // while tbins

}


void nest::inh_gamma_generator::calc_h(double h[], double t,double a, double b, double dt)
{

  double Hpre, Hpost;

  if (t < dt) {
    h[0] = 0.0;
    return;
  }

  // The following is PyGSL python code
  // which yields the correct result
  /*
  def calc_h(x,a,b):

    import pygsl
    from pygsl.sf import gamma_inc_Q

    dt = 1e-4

    Hpre = -log(gamma_inc_Q(5.0,(0.1-dt)/b)[0])[0]

    Hpost = -log(gamma_inc_Q(5.0,(0.1+dt)/b)[0])[0]

    return 0.5*(Hpost-Hpre)/dt
  */

  // This code in C should look like:

  gsl_sf_result result;
  gsl_sf_gamma_inc_Q_e (a, (t-dt)/b, &result);
  Hpre = -log(result.val);
  gsl_sf_gamma_inc_Q_e (a, (t+dt)/b, &result);
  Hpost = -log(result.val);

  h[0] = 0.5*(Hpost-Hpre)/dt;

}


}



#endif //HAVE_GSL
