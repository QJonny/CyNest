/*
 *  pseudoneuron.cpp
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2004-2005 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  This file is CONFIDENTIAL and PROPRIETARY.
 */

#include "pseudoneuron.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "arraydatum.h"
#include "dictutils.h"

nest::pseudoneuron::pseudoneuron()
  : Device(),
    exp_dev_(),
    rate_(0.0),
    h_(network()->get_resolution().get_ms()),
    overflow_(),
    spikes_rcvd_(0)
{
}

//!< @bug is it correct to copy exp_dev_?
nest::pseudoneuron::pseudoneuron(const pseudoneuron& n)
  : Device(n),
    exp_dev_(),  // has no member vars
    rate_(n.rate_),
    h_(n.h_),
    overflow_(n.overflow_),
    spikes_rcvd_(n.spikes_rcvd_)
{
}

void nest::pseudoneuron::init()
{
}

void nest::pseudoneuron::calibrate()
{
  // Set up overflow_ buffer here so that it gets set up properly
  // even when connecting a spike_detector using "sd pg Connect".
  // New elements are filled with 0, existing elements left untouched.
  overflow_.resize(targets_.size(), 0.0);
}

//
// Time Evolution Operator
//
void nest::pseudoneuron::update(thread thrd, Time const & origin,
					const long_t from, const long_t to)
{
  assert(to-from <= Scheduler::get_min_delay()); // kludge til cleanup

  // check if anything to do 
  if ( rate_ <= 0 )
    return;

  // get start time for slice
  const long_t start = origin.get_steps();

  // obtain rng
  librandom::RngPtr rng = network()->get_rng(thrd);

  // We now go through all targets and, for each send all spikes 
  // within T_limit_ms.  
  for (size_t prt = 0; prt < targets_.size(); ++prt)
  {
    Connection& tgt = targets_[prt];
    if (!tgt.is_valid())
      continue;
    
    // check if overflow has been set already, if not, set one
    // @todo this initialization should be moved to calibrate, once
    // that is thread-aware.
    assert(prt < overflow_.size());
    if (!overflow_[prt])
      overflow_[prt] = 1000.0 * exp_dev_(rng) / rate_;

    // We need to know in which time step we are, so we can check
    // if the device is active and have the timestamp for dispatch.
    // Computing both from a double_t Poisson spike time would be
    // to cumbersome, presumably.
    for (long_t step = start+from; step < start+to; ++step)
    {
      if (is_active(Time::step(step)))
      {
        // Compute time-limit in ms: 
        // We fire all spikes up to but excluding T_limit_ms.
        const double_t T_limit_ms = Time(Time::step(step+1)).get_ms();
        
        // now generate spikes
        while (overflow_[prt] <= T_limit_ms)
        {
          // We could use an EventInstance<SpikeEvent> here, as it is only
          // a local variable and thus not using the pool allocator.  But
          // since we have made SpikeEvent a concrete class, we may use it
          // directly.
          SpikeEvent se;
          se.set_sender(*this);
          se.set_offset(T_limit_ms - overflow_[prt]);
          network()->send(*this, se, step-start);
          
          // draw new spike time
          overflow_[prt] += 1000.0 * exp_dev_(rng)/rate_;
        }
      }
    }
  }
}

void nest::pseudoneuron::event_hook(thread t, DSSpikeEvent& e)
{

}

void nest::pseudoneuron::get_status(DictionaryDatum &d) const
{
  // start/stop properties
  Device::get_status(d);

  (*d)["rate"] = new DoubleDatum(rate_);
  (*d)["spike_rcvd"] = new IntegerDatum(spikes_rcvd_);

}

void nest::pseudoneuron::set_status(const DictionaryDatum &d)
{

  // start/stop properties
  Device::set_status(d);
  
  double tmp;
  if (updateValue<double>(d, "rate", tmp)) 
  {
    if (tmp < 0.0)
      throw BadProperty();
    rate_ = tmp;
  }
}

void nest::pseudoneuron::handle(SpikeEvent& e)
{
  assert(e.get_delay() > 0);
  const long_t Tdeliver = e.get_rel_delivery_steps(network()->get_slice_origin());
  assert(Tdeliver >= 0);

  ++spikes_rcvd_;
}

nest::port nest::pseudoneuron::connect(Node& r)
{
  EventInstance<SpikeEvent> se;
  return net_->connect(*this, r, se);
}
