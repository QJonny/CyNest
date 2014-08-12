/*
 *  maturing_connector.h
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

/*
 * first version
 * author: Moritz Helias
 * date: october 2006
 */

/* BeginDocumentation
 * Name: connector_maturating - Connector to create maturating synapses.
 * A newly created synapse remains in a premature state not producing
 * postsynaptic potentials.
 * It measures the correlation between pre- and postsynaptic spike trains
 * and maturates after a given time, if this correlation exceeds a threshold.
 * Otherwise, the synapse dies.
 *
 * Parameters (via connector_stdp_properties) :
 * t_d          time until decision will take place
 * w            weight of a mature synapse
 * 
 */

#ifndef MATURING_CONNECTOR_H
#define MATURING_CONNECTOR_H

#include "generic_connector_model.h"
#include "generic_connector.h"
#include "maturing_connection.h"
#include <cmath>


namespace nest
{


class MaturingConnector :
  public GenericConnectorBase< MaturingConnection,
                               MaturingCommonProperties,
                               GenericConnectorModelBase< MaturingConnection, MaturingCommonProperties, MaturingConnector >
                             >
  {

  // this is the ConnectorModel specialization which belongs to this kind of Connector specialization
  typedef GenericConnectorModelBase< MaturingConnection, MaturingCommonProperties, MaturingConnector > GCMT;

 public:

  /**
   * Default constructor.
   * \param cm ConnectorModel, which created this Connector.
   */
  MaturingConnector(GCMT &cm);

  /**
   * Register a new connection at the sender side.
   */ 
  void register_connection(Node&, Node&);

  /**
   * Register a new connection at the sender side.
   * Use given weight and delay.
   */ 
  void register_connection(Node&, Node&, double_t, double_t);

  /**
   * Register a new connection at the sender side. 
   * Use given dictionary for parameters.
   */ 
  void register_connection(Node&, Node&, DictionaryDatum&);

  /**
   * Send an event to this connector.
   * The connector will propagate it to all its targets.
   */
  void send(Event& e);
 
 private:

  // activation state of the NMDA receptors
  double_t h_;

  // recovery of synaptic resources
  double_t x_;
  
};

} // namespace

#endif /* #ifndef MATURING_CONNECTOR_H */

