/*
 *  inh_gamma_generator_ng.cpp
 *
 *  Based loosely on the poisson_generator
 *
 *  AUTHOR: Eilif Muller
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2007 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#include "inh_gamma_generator_ng.h"

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

nest::inh_gamma_generator_ng::inh_gamma_generator_ng()
  : Node(), StimulatingDevice<SpikeEvent>(),
          exprand_dev_()
{
  calibrate();

  tbins_.clear();
  tbins_.push_back(0.0);
  a_.clear();
  a_.push_back(1.0);
  b_.clear();
  b_.push_back(1.0);
  rmax_ = 1.0;

  tspike_.clear();
  tlast_.clear();
  bin_position_ = 0;
  synapse_context_ = -1;

}

nest::inh_gamma_generator_ng::inh_gamma_generator_ng(inh_gamma_generator_ng const &pg)
  : Node(pg), StimulatingDevice<SpikeEvent>(pg),
          exprand_dev_(pg.exprand_dev_)
{
  calibrate();

  tbins_ = pg.tbins_;
  a_ = pg.a_;
  b_ = pg.b_;
  rmax_ = pg.rmax_;
  bin_position_ = pg.bin_position_;
  num_connections_ = pg.num_connections_;
  tspike_ = pg.tspike_;
  tlast_ = pg.tlast_;
  synapse_context_ = pg.synapse_context_;

  // clear tspike_
  size_t i;

  for (i=0;i<tspike_.size();i++) {
    tspike_[i]=0.0;
    //tlast_[i]=-1000.0;
    tlast_[i]=0.0;
  }

  bin_position_ = 0;
}

void nest::inh_gamma_generator_ng::init_parameters_(Node const* model)
{
  inh_gamma_generator_ng const* proto = dynamic_cast<inh_gamma_generator_ng const*>(model);
  assert(proto != 0);
  
  StimulatingDevice<SpikeEvent>::init_parameters(*proto);
  network()->warning("inh_gamma_generator_ng::init_parameters_",
    "Parameter initialization not fully implemented, be careful!");
}

void nest::inh_gamma_generator_ng::init_dynamic_state_(Node const* model)
{ 
  inh_gamma_generator_ng const* proto = dynamic_cast<inh_gamma_generator_ng const*>(model);
  assert(proto != 0);
  
  StimulatingDevice<SpikeEvent>::init_dynamic_state(*proto);
  network()->warning("inh_gamma_generator_ng::init_dynamic_state_",
    "Dynamic state initialization not fully implemented, be careful!");
}

void nest::inh_gamma_generator_ng::calibrate()
{
  StimulatingDevice<SpikeEvent>::calibrate();

  // rate_ is in Hz, dt in ms, so we have to convert
  // from s to ms
  //poisson_dev_.set_lambda(Time::get_resolution().get_ms() * rate_*1e-3);
}

//
// Time Evolution Operator
//
void nest::inh_gamma_generator_ng::update(Time const & T, const long_t from, const long_t to)
{

  assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
  assert(from < to);

  assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
  assert(from < to);

  double_t time;
  librandom::RngPtr rng = net_->get_rng(get_thread());
  double c;
  double r;

  for (long_t lag = from; lag < to; lag++)
    {

      // get time in ms.
      time = (T + Time(Time::step(lag))).get_ms();

      if ( !is_active(T + Time(Time::step(lag))) || rmax_ <= 0.0 )
        {
          // no spikes to be generated
          return;
        }

      // find position in histogram for this time step
      while (bin_position_+1 < tbins_.size() && tbins_[bin_position_+1] <= time) 
	bin_position_++;


      for (size_t p=0;p<tspike_.size();p++) {

	// spike still in the future?
	while (tspike_[p]<=time) {

	  // get hazard

	  // b is in seconds so convert to milliseconds

	  calc_h(&c,tspike_[p]-tlast_[p],a_[bin_position_],b_[bin_position_]*1000.0,
		 Time::get_resolution().get_ms() / 10.0);

	  // now get rate in Hz
	  c*=1000.0;

	  // check if we should really emit a spike (thinning step)

	  // uniform random # on [0,rmax_)
	  r = rng->drand()*rmax_;

	  if (r<c) {
	    // spike!
	    tlast_[p] = tspike_[p];

	    SpikeEvent se;
	    se.set_port(p);
	    network()->send(*this, se, lag);

	  }

	  // get next "possible" spike time
	  tspike_[p] += 1000.0/rmax_*exprand_dev_(rng);  //rmax_ is in Hz

	} // while tspike_[p]<time



      } // for p



    } // for lag

}

void nest::inh_gamma_generator_ng::get_status(DictionaryDatum &d) const
{
  // start/stop properties
  StimulatingDevice<SpikeEvent>::get_status(d);

  (*d)["tbins"] = DoubleVectorDatum(new std::vector<double>(tbins_));
  (*d)["a"] = DoubleVectorDatum(new std::vector<double>(a_));
  (*d)["b"] = DoubleVectorDatum(new std::vector<double>(b_));
  (*d)["rmax"] = new DoubleDatum(rmax_);
  (*d)["tspike"] = DoubleVectorDatum(new std::vector<double>(tspike_));
  (*d)["tlast"] = DoubleVectorDatum(new std::vector<double>(tlast_));
  (*d)["synapse_type_id"] = new IntegerDatum(synapse_context_);
}

void nest::inh_gamma_generator_ng::set_status(const DictionaryDatum &d)
{

  // start/stop properties

  double tmp;
  if ( updateValue<double>(d, "rmax", tmp) )
  {
      if ( tmp < 0.0 )
        throw BadProperty();

      rmax_ = tmp;

      calibrate();
  }

  Name tbins_name("tbins");

  if(d->known(tbins_name))
    {
      tbins_ = = getValue<std::vector<double> >(d->lookup(tbins_name));
    }
  
  Name tspike_name("tspike");

  if(d->known(tspike_name))
    {
      tspike_ = getValue<std::vector<double> >(d->lookup(tspike_name));
    }

  Name tlast_name("tlast");

  if(d->known(tlast_name))
    {
      tlast_ = getValue<std::vector<double> >(d->lookup(tlast_name));
    }

  Name a_name("a");

  if(d->known(a_name))
    {
      a_ = getValue<std::vector<double> >(d->lookup(a_name));
    }

  Name b_name("b");

  if(d->known(b_name))
    {
      b_ = getValue<std::vector<double> >(d->lookup(b_name));
    }


}

void nest::inh_gamma_generator_ng::calc_h(double h[], double t,double a, double b, double dt)
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
