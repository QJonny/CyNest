/*
 *  static_connection_mult0.h
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2004-2010 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */


/* BeginDocumentation
  Name: static_synapse_mult0 - Synapse type for debugging.

  Description:
  This version is static_synapse is for debugging only. It sets
  the multiplicity of the event to 0 before sending it, and thus
  may uncover interferences in the spike delivery due to the
  re-use of the Event. See #426.

  FirstVersion: May 2010
  Author: Hans Ekkehard Plesser, Jochen Martin Eppler, Moritz Helias

  Transmits:
  SpikeEvent
  
  Remarks: Refactored for new connection system design, March 2007

  SeeAlso: synapsedict, tsodyks_synapse, stdp_synapse
*/

#ifndef STATICCONNECTION_MULT0_H
#define STATICCONNECTION_MULT0_H

#include "connection_het_wd.h"

namespace nest
{

/**
 * Class representing a static connection. A static connection has the properties weight, delay and receiver port.
 * This class also serves as the base class for dynamic synapses (like TsodyksConnection, STDPConnection).
 * A suitale Connector containing these connections can be obtained from the template GenericConnector.
 */
class StaticConnectionMult0 : public ConnectionHetWD
{

 public:

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  StaticConnectionMult0() : ConnectionHetWD() {}

  /**
   * Default Destructor.
   */
  ~StaticConnectionMult0() {}

  // overloaded for all supported event types
  using Connection::check_event;
  void check_event(SpikeEvent&) {}

  void send(Event& e, double_t, const CommonSynapseProperties &);
};

inline
void StaticConnectionMult0::send(Event& e, double_t, const CommonSynapseProperties &)
{
  SpikeEvent* se = dynamic_cast<SpikeEvent*>(&e);
  assert(se);

  se->set_multiplicity(0);
  se->set_weight(weight_);
  se->set_delay(delay_);
  se->set_receiver(*target_);
  se->set_rport(rport_);
  (*se)();
}

} // namespace

#endif /* #ifndef STATICCONNECTION_MULT0_H */
