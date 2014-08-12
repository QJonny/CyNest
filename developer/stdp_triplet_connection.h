/*
 *  stdp_triplet_connection.h
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
 *  Author: Eilif Muller (modified from stdp_synapse)
 *
 */

#ifndef STDP_TRIPLET_CONNECTION_H
#define STDP_TRIPLET_CONNECTION_H

/* BeginDocumentation
  Name: stdp_triplet_synapse - Synapse type for spike-timing dependent plasticity accounting for spike triplets as described in [1].

  Description:
   stdp_triplet_synapse is a connector to create synapses with spike time 
   dependent plasticity accounting for spike triplets (as defined in [1]). 
   
   Here, a multiplicative weight dependence is added (in contrast to [1]) 
   to depression resulting in a stable weight distribution.
  
  Parameters:

   tau_plus     time constant of STDP window, potentiation 
                (tau_minus defined in post-synaptic neuron)

   tau_x        time constant of triplet potentiation

   tau_y        time constant of triplet depression

   A_2p         weight of pair potentiation rule

   A_2m         weight of pair depression rule

   A_3p         weight of triplet potentiation rule

   A_3m         weight of triplet depression rule


  References:

   [1] J.-P. Pfister & W. Gerstner (2006) Triplets of Spikes in a Model 
   of Spike Timing-Dependent Plasticity.  The Journal of Neuroscience 
   26(38):9673-9682; doi:10.1523/JNEUROSCI.1425-06.2006



  FirstVersion: Nov 2007
  Author: Moritz Helias, Abigail Morrison, Eilif Muller
  SeeAlso: synapsedict, stdp_synapse, tsodyks_synapse, static_synapse
*/

#include "static_connection.h"
#include "archiving_node.h"
#include "generic_connector.h"
#include <cmath>

namespace nest
{
  //class CommonProperties;

  class STDPTripletConnection : public StaticConnection {

  public:
  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  STDPTripletConnection();
  
  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  STDPTripletConnection(const STDPTripletConnection &);

  /**
   * Default Destructor.
   */
  ~STDPTripletConnection() {}

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
   * \param t_lastspike last spike emitted by presynaptic neuron
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
  // this needs to be commented out on Kei MH 11-06-21
  using Connection::check_event;
  void check_event(SpikeEvent&) {}
  
 private:

  inline double_t facilitate_(double_t w, double_t kplus, double_t ky);

  inline double_t depress_(double_t w, double_t kminus, double_t kx);

  // data members of each connection
  double_t tau_plus_;
  double_t tau_x_;

  double_t A_2p_;
  double_t A_2m_;

  double_t A_3p_;
  double_t A_3m_;

  double_t Kplus_;
  double_t Kx_;

  };


inline
double_t STDPTripletConnection::facilitate_(double_t w, double_t kplus, double_t ky)
{
  return w+(A_2p_+A_3p_*ky)*kplus;

}

inline 
double_t STDPTripletConnection::depress_(double_t w, double_t kminus, double_t kx)
{
  //double new_w = w - w * kminus*(A_2m_ + A_3m_*kx); 
  double new_w = w - kminus*(A_2m_ + A_3m_*kx); 
  return new_w > 0.0 ? new_w : 0.0;
}


inline 
  void STDPTripletConnection::check_connection(Node & s, Node & r, port receptor_type, double_t t_lastspike)
{
  Connection::check_connection(s, r, receptor_type, t_lastspike);
  r.register_stdp_connection(t_lastspike - Time(Time::step(delay_)).get_ms());
}

/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
inline
void STDPTripletConnection::send(Event& e, double_t t_lastspike, const CommonSynapseProperties &)
{
  // synapse STDP depressing/facilitation dynamics
  
  double_t t_spike = e.get_stamp().get_ms();
  // t_lastspike_ = 0 initially
  double_t dendritic_delay = Time(Time::step(delay_)).get_ms(); 
    
  //get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque<histentry>::iterator start;
  std::deque<histentry>::iterator finish;    
  target_->get_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay, 
			       &start, &finish);
  //facilitation due to post-synaptic spikes since last pre-synaptic spike
  double_t minus_dt;
  double_t ky;
  while (start != finish)
  {      
    // EM: post-synaptic spike is delay by dendritic_delay so that
    // it is effectively late by that much at the synapse.
    minus_dt = t_lastspike - (start->t_ + dendritic_delay);
    // EM: dendritic delay?  Accounted for in the time minus_dt
    // no need to modify ky below to account for it.
    // EM: subtract 1.0 yields the triplet_Kminus value just prior to
    // the post synaptic spike, implementing the t-epsilon in 
    // Pfister et al, 2006
    ky = start->triplet_Kminus_-1.0;
    //ky = target_->get_triplet_K_value(start)-1.0;
    //ky = target_->get_triplet_K_value(start->t_-dendritic_delay);
    ++start;
    if (minus_dt == 0)
      continue;
    weight_ = facilitate_(weight_, Kplus_ * std::exp(minus_dt / tau_plus_),ky);
  }

  //depression due to new pre-synaptic spike
  Kx_ *= std::exp((t_lastspike - t_spike) / tau_x_);

  // dendritic delay means we must look back in time by that amount
  // for determining the K value, because the K value must propagate
  // out to the synapse
  weight_ = depress_(weight_, target_->get_K_value(t_spike - dendritic_delay),Kx_);

  Kx_ += 1.0;

  e.set_receiver(*target_);
  e.set_weight(weight_);
  e.set_delay(delay_);
  e.set_rport(rport_);
  e();

  Kplus_ = Kplus_ * std::exp((t_lastspike - t_spike) / tau_plus_) + 1.0;

}

} // of namespace nest

#endif // of #ifndef STDP_TRIPLET_CONNECTION_H
