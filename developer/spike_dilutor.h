/*
 *  spike_dilutor.h
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2011 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#ifndef SPIKE_DILUTOR_H
#define SPIKE_DILUTOR_H

#include "nest.h"
#include "event.h"
#include "node.h"
#include "stimulating_device.h"
#include "connection.h"
#include "ring_buffer.h"

namespace nest
{

/* BeginDocumentation
Name: spike_dilutor - repeats incoming spikes with a certain probability.

Description:
  The device repeats incoming spikes with a certain probability.
  Targets will receive diffenrent spike trains.

Remarks:
  In parallel simulations, a copy of the device is present on each process
  and spikes are collected only from local sources.

Parameters:
   The following parameters appear in the element's status dictionary:
   p_copy double - Copy probability

Sends: SpikeEvent

Author: Adapted from mip_generator by Kunkel, Oct 2011
SeeAlso: mip_generator
*/

  class spike_dilutor : public Node
  {

  public:

    spike_dilutor();
    spike_dilutor(const spike_dilutor & rhs);

    bool has_proxies() const {return false;}
    bool local_receiver() const {return true;}

    using Node::connect_sender;
    using Node::handle;
    using Node::event_hook;

    port check_connection(Connection &, port);
    port connect_sender(SpikeEvent &, port);
    void handle(SpikeEvent &);

    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &) ;

  private:

    void init_state_(const Node &);
    void init_buffers_();
    void calibrate();

    void update(Time const &, const long_t, const long_t);

    void event_hook(DSSpikeEvent&);

    // ------------------------------------------------------------

    /**
     * Store independent parameters of the model.
     */
    struct Parameters_ {
      double_t p_copy_; //!< copy probability for each incoming spike

      Parameters_();  //!< Sets default parameter values
      Parameters_(const Parameters_ &);

      void get(DictionaryDatum &) const;  //!< Store current values in dictionary
      void set(const DictionaryDatum &);  //!< Set values from dicitonary
    };

    struct Buffers_ {
      RingBuffer n_spikes_;
    };

    // ------------------------------------------------------------

    StimulatingDevice<SpikeEvent> device_;
    Parameters_ P_;
    Buffers_ B_;
  };

  inline
  port spike_dilutor::check_connection(Connection & c, port receptor_type)
  {
    DSSpikeEvent e;
    e.set_sender(*this);
    c.check_event(e);
    return c.get_target()->connect_sender(e, receptor_type);
  }

  inline
  port spike_dilutor::connect_sender(SpikeEvent &, port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return 0;
  }

  inline
  void spike_dilutor::get_status(DictionaryDatum & d) const
  {
    P_.get(d);
    device_.get_status(d);
  }

  inline
  void spike_dilutor::set_status(const DictionaryDatum & d)
  {
    Parameters_ ptmp = P_;  // temporary copy in case of errors
    ptmp.set(d);            // throws if BadProperty

    // We now know that ptmp is consistent. We do not write it back
    // to P_ before we are also sure that the properties to be set
    // in the parent class are internally consistent.
    device_.set_status(d);

    // if we get here, temporaries contain consistent set of properties
    P_ = ptmp;
  }

} // namespace nest

#endif
