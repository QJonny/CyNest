/*
 *  lossy_connection.h
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

#ifndef LOSSY_CONNECTION_H
#define LOSSY_CONNECTION_H

/* BeginDocumentation
  Name: lossy_synapse - Static synapse which transmits spike with a certain probability.

  Description:
   lossy_synapse is a static synapse that transmits spikes with probability p_transmit.

  Parameters:
   p_transmit   double - Spike transmission probability

  Transmits: SpikeEvent

  FirstVersion: Nov 2011
  Author: Susanne Kunkel
  SeeAlso: synapsedict,  static_synapse
*/

#include "connection_het_wd.h"
#include "nestmodule.h"

namespace nest
{
  class LossyConnection : public ConnectionHetWD
  {

    public:

    /**
     * Default Constructor.
     * Sets default values for all parameters. Needed by GenericConnectorModel.
     */
    LossyConnection();

    /**
     * Copy constructor.
     * Needs to be defined properly in order for GenericConnector to work.
     */
    LossyConnection(const LossyConnection &);

    /**
     * Default Destructor.
     */
    ~LossyConnection() {}

    /**
     * Get all properties of this connection and put them into a dictionary.
     */
    void get_status(DictionaryDatum & d) const;

    /**
     * Set properties of this connection from the values given in dictionary.
     */
    void set_status(const DictionaryDatum & d, ConnectorModel & cm);

    /**
     * Set properties of this connection from position p in the properties
     * array given in dictionary.
     */
    void set_status(const DictionaryDatum & d, index p, ConnectorModel & cm);

    /**
     * Create new empty arrays for the properties of this connection in the given
     * dictionary. It is assumed that they are not existing before.
     */
    void initialize_property_arrays(DictionaryDatum & d) const;

    /**
     * Append properties of this connection to the given dictionary. If the
     * dictionary is empty, new arrays are created first.
     */
    void append_properties(DictionaryDatum & d) const;

    /**
     * Send an event to the receiver of this connection.
     * \param e The event to send
     * \param t_lastspike Point in time of last spike sent
     * \param cp common properties of all synapses (empty)
     */
    void send(Event & e, double_t, const CommonSynapseProperties &);

    // overloaded for all supported event types
    using Connection::check_event;
    void check_event(SpikeEvent &) {}
    void check_event(DSSpikeEvent &) {}

  private:

    // data members of each connection
    double_t p_transmit_;
  };

  inline
  void LossyConnection::send(Event & e, double_t, const CommonSynapseProperties &)
  {
    SpikeEvent e_spike = static_cast<SpikeEvent &>(e);

    librandom::RngPtr rng = NestModule::get_network().get_rng(target_->get_thread());
    ulong_t n_spikes_in = e_spike.get_multiplicity();
    ulong_t n_spikes_out = 0;

    for ( ulong_t n = 0; n < n_spikes_in; n++ )
      if ( rng->drand() < p_transmit_ )
	n_spikes_out++;

    if ( n_spikes_out > 0 )
    {
      e_spike.set_multiplicity(n_spikes_out);
      e.set_receiver(*target_);
      e.set_weight(weight_);
      e.set_delay(delay_);
      e.set_rport(rport_);
      e();
    }

    e_spike.set_multiplicity(n_spikes_in);  // probably not needed
  }

} // of namespace nest

#endif // of #ifndef LOSSY_CONNECTION_H
