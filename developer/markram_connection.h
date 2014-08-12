/*
 *  markram_connection.h
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2004 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#ifndef MARKRAM_CONNECTION_H
#define MARKRAM_CONNECTION_H

#include "connection_het_wd.h"

/* BeginDocumentation
  Name: markram_synapse - Synapse type with phenomenological depression and facilitation dynamics

  Description:

   Implemented is the synapse model described in Eq (2) and (3) in [1]
   or En (6) in [2].  Whereby Eq (2) in [1] seems to have an error in
   the subscript of u_{n+1}.  It should be u_{n}.

   The markram_synapse is a simplification of the tsodyks_synapse
   (synaptic time course is neglected).
   

  References:
   [1] Markram, Wang, Tsodyks (1998) Differential Signaling via the same axon 
       of neocortical pyramidal neurons.  PNAS, vol 95, pp. 5323-5328.

   [2] D. Sussillo, T. Toyoizumi, and W. Maass. Self-tuning of neural circuits through
   short-term synaptic plasticity. Journal of Neurophysiology, 97:4079-4095, 2007.

  FirstVersion: Sept 2007
  Author: Moritz Helias, Eilif Muller
  SeeAlso: synapsedict, tsodyks_synapse, stdp_synapse, static_synapse
*/

/**
 * Class representing a synapse with Markram short term plasticity.
 * A suitale Connector containing these connections can be obtained from the template GenericConnector.
 */

namespace nest {

  //class CommonProperties;

class MarkramConnection : public ConnectionHetWD
{
 public:

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  MarkramConnection();

  /**
   * Default Destructor.
   */
  virtual ~MarkramConnection() {}

  // Import overloaded virtual function set to local scope. 
  using Connection::check_event;

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  virtual void get_status(DictionaryDatum & d) const;
  
  /**
   * Set properties of this connection from the values given in dictionary.
   */
  virtual void set_status(const DictionaryDatum & d, ConnectorModel &cm);

  /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  virtual void set_status(const DictionaryDatum & d, index p, ConnectorModel &cm);

  /**
   * Create new empty arrays for the properties of this connection in the given
   * dictionary. It is assumed that they are not existing before.
   */
  void initialize_property_arrays(DictionaryDatum & d) const;

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  virtual void append_properties(DictionaryDatum & d) const;

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param t_lastspike Point in time of last spike sent.
   * \param cp Common properties to all synapses (empty).
   */
  void send(Event& e, double_t t_lastspike, const CommonSynapseProperties &cp);

  // overloaded for all supported event types
  void check_event(SpikeEvent&) {}
 
 private:
  double_t tau_fac_;   //!< [ms] time constant for fascilitation
  double_t tau_rec_;   //!< [ms] time constant for recovery
  double_t U_;         //!< asymptotic value of probability of release
  double_t R_;         //!< amount of resources in recovered state
  double_t u_;         //!< actual probability of release
};


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
inline
void MarkramConnection::send(Event& e, double_t t_lastspike, const CommonSynapseProperties &)
{
  double_t h = e.get_stamp().get_ms() - t_lastspike;

  // t_lastspike_ = 0 initially

  // TODO: use expm1 here instead, where applicable
  // propagation t_lastspike -> t_spike

  // iterative form of:
  // dU/dt = -(u-U)/tau_fac
  u_ = (u_-U_)*std::exp(-h/tau_fac_) + U_;
  
  // iterative form of:
  // dR/dt = -(R-1.0)/tau_rec
  R_ = (R_-1.0)*std::exp(-h/tau_rec_) + 1.0;

  // postsynaptic step caused by incoming spike
  double_t delta_y_tsp = u_*R_;

  // delta functions for u,R
  R_ -= R_*u_;
  u_ += U_*(1.0-u_);

  // send the spike to the target
  e.set_receiver(*target_);
  e.set_weight( weight_ * delta_y_tsp );
  e.set_delay( delay_ );
  e.set_rport( rport_ );
  e();
}
 
} // namespace

#endif // MARKRAM_CONNECTION_H
