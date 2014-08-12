/*
 *  stdp_connection.h
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

#ifndef SCHEMMEL_S1_CONNECTION_H
#define SCHEMMEL_S1_CONNECTION_H

/* BeginDocumentation

  Name: schemmel_s1_synapse - Synapse type for
spike-timing dependent plasticity, and short-term depression and
faciliation dynamics which emulates the FACETS stage 1 VLSI Hardware
synapse described in [1,2].  

Note: Spike pairing is nearest neighbor as defined in [2].

NOTE: The weight dependent lookup is presently unimplemented pending further
details thereof.

  Description: schemmel_s1_synapse is a connector to create synapses
  with spike time dependent plasticity and short-term depression and
  facilitation dynamics (as defined in [1,2]). 
  
  Parameters:

  STDP:
   tau_plus      time constant of STDP window, potentiation 
   tau_minus     time constant of STDP window, depression
   Note: Tau_minus in the target neuron is not used
   as a nearest neighbor spike pairing scheme is implemented.

   A_plus       potentiation factor
   A_minus      stdp depression factor

  Dynamic synapses (short-term depression+facilitation):

   mode (int)   0=depressing, 1 = facilitating, else disabled
   tau_rec      recovery time constant
   U_SE         usable synaptic efficacy
   Note: A_SE in [2] is represented here by the connection weight
   Synaptic dynamics do not alter the connection weight,
   but do alter the weight of the event transmitted to the target
   each pre-synaptic spike.

  References: 

[1] J. Schemmel et al. (2006). Implementing Synaptic 
    Plasticity in a VLSI Spiking Neural Network Model.  Proc. of the 
    2006 International Joint Conference on Neural Networks (IJCNN
    2006), 1--6, IEEE Press.  

[2] J. Schemmel et al. (2007).  Modeling Synaptic Plasticity within 
    Networks of Highly Accelerated I&F Neurons.  Proc. of the International 
    Symposium on Circuits and Systems 2007 (ISCAS2007), in press. 


  FirstVersion: Nov 2007
  Author: Eilif Muller, modified from stdp_synapse.
  SeeAlso: synapsedict, stdp_synapse, tsodyks_synapse, static_synapse

*/

#include "connection_het_wd.h"
#include "archiving_node.h"
#include "generic_connector.h"
#include <cmath>

namespace nest
{
  //class CommonProperties;

  class SchemmelS1Connection : public ConnectionHetWD {

  public:
  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  SchemmelS1Connection();
  
  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  SchemmelS1Connection(const SchemmelS1Connection &);

  /**
   * Default Destructor.
   */
  ~SchemmelS1Connection() {}

  // Import overloaded virtual function set to local scope. 
  using Connection::check_event;
  
  /*
   * This function calls check_connection on the sender and checks if the receiver
   * accepts the event type and receptor type requested by the sender.
   * Node::check_connection() will either confirm the receiver port by returning
   * true or false if the connection should be ignored.
   * We have to override the base class' implementation, since for STDP
   * connections we have to call register_stdp_connection on the target neuron
   * to inform the Archiver to collect spikes for this connection.
   *
   * \param s The source node
   * \param r The target node
   * \param receptor_type The ID of the requested receptor type
   * \param t_lastspike last spike emitted by presynaptic neuron (in ms)
   */
  void check_connection(Node & s, Node & r, port receptor_type, double_t t_lastspike);

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status(DictionaryDatum & d) const;
  
  /**
   * Set properties of this connection from the values given in dictionary.
   */
  void set_status(const DictionaryDatum & d, ConnectorModel &cm);

  /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void set_status(const DictionaryDatum & d, index p, ConnectorModel &cm);

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
   * \param t_lastspike Point in time of last spike sent.
   * \param cp common properties of all synapses (empty).
   */
  void send(Event& e, double_t t_lastspike, const CommonSynapseProperties &cp);

  // overloaded for all supported event types
  void check_event(SpikeEvent&) {}
  
 private:

  inline double_t facilitate_(double_t w, double_t kplus);

  inline double_t depress_(double_t w, double_t kminus);

  // data members of each connection

  //STDP

  double_t tau_plus_;
  double_t tau_minus_;
  double_t A_plus_;
  double_t A_minus_;

  //Short-term dynamics

  double_t tau_rec_;
  double_t U_SE_;
  double_t I_; // inactivate synaptic efficacy (state variable)
  long_t mode_;

  };


inline
double_t SchemmelS1Connection::facilitate_(double_t w, double_t kplus)
{
  // Note: Weight update lookup table should be implemented here

  // An additive potentiation rule is implemented for now.

  return w+kplus;


}

inline 
double_t SchemmelS1Connection::depress_(double_t w, double_t kminus)
{
  // Note: Weight update lookup table should be implemented here

  // An multiplicative depression rule is implemented for now.

  double_t new_w = w - ( w * kminus);
  return new_w > 0.0 ? new_w : 0.0;

}


inline 
  void SchemmelS1Connection::check_connection(Node & s, Node & r, port receptor_type, double_t t_lastspike)
{
  ConnectionHetWD::check_connection(s, r, receptor_type, t_lastspike);
  r.register_stdp_connection(t_lastspike - Time(Time::step(delay_)).get_ms());
}

/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
inline
void SchemmelS1Connection::send(Event& e, double_t t_lastspike, const CommonSynapseProperties &)
{
  // Schemmel Stage 1 STDP and short-term depressing/facilitation dynamics
  
  double_t t_spike = e.get_stamp().get_ms();
  // t_lastspike_ = 0 initially
  double_t dendritic_delay = Time(Time::step(delay_)).get_ms(); 

  //get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque<histentry>::iterator start;
  std::deque<histentry>::iterator finish;    
  target_->get_history(t_lastspike - dendritic_delay, 
		       t_spike - dendritic_delay, &start, &finish);

  //facilitation+depression due to post-synaptic spikes since last
  //pre-synaptic spike (process new spikes since last spike)

  double_t minus_prepost_dt,minus_postpre_dt;
  while (start != finish)
  {
    minus_prepost_dt = t_lastspike - (start->t_ + dendritic_delay);
    minus_postpre_dt =  (start->t_ + dendritic_delay) - t_spike;

    ++start;
    if (minus_prepost_dt < 0) {
      weight_ = facilitate_(weight_, A_plus_ * std::exp(minus_prepost_dt / tau_plus_));
    }

    if (minus_postpre_dt < 0) {
      weight_ = depress_(weight_, A_minus_ * std::exp(minus_postpre_dt / tau_minus_));
    }


  }


  // short-term synaptic dynamics:


  double_t h = t_spike - t_lastspike;

  I_ = I_*std::exp(-h/tau_rec_); 

  double_t ds_weight_;

  if (mode_==0) {
    // depression
    ds_weight_ = 1.0 - I_ > 0.0 ? 1.0-I_ : 0.0;
  }
  else if (mode_==1) {
    // facilitation
    ds_weight_ = I_;
  }
  else ds_weight_ = 1.0; // disable short term dynamics

  // inactive synaptic efficacy jump

  I_ += U_SE_*(1.0-I_);


  e.set_receiver(*target_);
  e.set_weight(weight_*ds_weight_);
  e.set_delay(delay_);
  e.set_rport(rport_);
  e();

}

} // of namespace nest

#endif // of #ifndef SCHEMMEL_S1_CONNECTION_H
