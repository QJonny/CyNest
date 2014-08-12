/*
 *  stdp_connection_hom.h
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

#ifndef STDP_CONNECTION_RSNN_SPIKEPAIRING_HOM_H
#define STDP_CONNECTION_RSNN_SPIKEPAIRING_HOM_H

/* BeginDocumentation
  Name: stdp_rsnn_spikepairing_synapse_hom - Synapse type for spike-timing dependent
   plasticity using homogeneous parameters, i.e. all synapses have the same parameters.

  Description:
   stdp_synapse is a connector to create synapses with spike time 
   dependent plasticity (as defined in [1]). Here the weight dependence
   exponent can be set separately for potentiation and depression.
   The reduced symmetric nearest-neighbor spike pairing scheme [5]
   as implemented in the FACETS (BrainScaleS) hardware systems [6] is used
   instead of an all-to-all spike pairing scheme.
   Consequently tau_minus_ is calculated within this synapse and not
   at the neuron site via Kplus_ like in stdp_connection_hom.
   
  Examples:
   multiplicative STDP [2]  mu_plus = mu_minus = 1.0
   additive STDP       [3]  mu_plus = mu_minus = 0.0
   Guetig STDP         [1]  mu_plus = mu_minus = [0.0,1.0]
   van Rossum STDP     [4]  mu_plus = 0.0 mu_minus = 1.0 

  Parameters:
   tau_plus        double - Time constant of STDP window, potentiation in ms 
   tau_minus_stdp  double - Time constant of STDP window, depression in ms 
   lambda          double - Step size
   alpha           double - Asymmetry parameter (scales depressing increments as alpha*lambda)
   mu_plus         double - Weight dependence exponent, potentiation
   mu_minus        double - Weight dependence exponent, depression
   Wmax            double - Maximum allowed weight

  Transmits: SpikeEvent
   
  References:
   [1] Guetig et al. (2003) Learning Input Correlations through Nonlinear
       Temporally Asymmetric Hebbian Plasticity. Journal of Neuroscience

   [2] Rubin, J., Lee, D. and Sompolinsky, H. (2001). Equilibrium
       properties of temporally asymmetric Hebbian plasticity, PRL
       86,364--367

   [3] Song, S., Miller, K. D. and Abbott, L. F. (2000). Competitive
       Hebbian learning through spike-timing-dependent synaptic
       plasticity, Nature Neuroscience 3:9,919--926

   [4] van Rossum, M. C. W., Bi, G-Q and Turrigiano, G. G. (2000). 
       Stable Hebbian learning from spike timing-dependent
       plasticity, Journal of Neuroscience, 20:23,8812--8821

   [5] Morrison, A., Diesmann, M., and Gerstner, W. (2008).
       Phenomenological models of synaptic plasticity based on
       spike-timing, Biol. Cybern., 98,459--478

   [6] Schemmel, J., Gruebl, A., Meier, K., and Mueller, E. (2006).
       Implementing synaptic plasticity in a VLSI spiking neural
       network model, In Proceedings of the 2006 International
       Joint Conference on Neural Networks, pp.1--6, IEEE Press

  FirstVersion: March 2006
  Author: Moritz Helias, Abigail Morrison, Thomas Pfeil
  SeeAlso: synapsedict, tsodyks_synapse, static_synapse
*/

#include "connection_het_wd.h"
#include "archiving_node.h"
#include <cmath>

namespace nest
{

  /**
   * Class containing the common properties for all synapses of type STDPRSNNSpikePairingConnectionHom.
   */
  class STDPRSNNSpikePairingHomCommonProperties : public CommonSynapseProperties
    {
      friend class STDPRSNNSpikePairingConnectionHom;

    public:

      /**
       * Default constructor.
       * Sets all property values to defaults.
       */
      STDPRSNNSpikePairingHomCommonProperties();

      /**
       * Get all properties and put them into a dictionary.
       */
      void get_status(DictionaryDatum & d) const;

      /**
       * Set properties from the values given in dictionary.
       */
      void set_status(const DictionaryDatum & d, ConnectorModel& cm);

      // overloaded for all supported event types
      void check_event(SpikeEvent&) {}

    private:

      // data members common to all connections
      double_t tau_plus_;
      double_t tau_minus_;
      double_t lambda_;
      double_t alpha_;
      double_t mu_plus_;
      double_t mu_minus_;
      double_t Wmax_;
    };



  /**
   * Class representing an STDP connection with homogeneous parameters, i.e. parameters are the same for all synapses.
   */
  class STDPRSNNSpikePairingConnectionHom : public ConnectionHetWD
  {

  public:
  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  STDPRSNNSpikePairingConnectionHom();

  /**
   * Copy constructor from a property object.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  STDPRSNNSpikePairingConnectionHom(const STDPRSNNSpikePairingConnectionHom &);

  /**
   * Default Destructor.
   */
  virtual ~STDPRSNNSpikePairingConnectionHom() {}

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
   * \param t_lastspike last spike produced by presynaptic neuron (in ms)
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
   */
  void send(Event& e, double_t t_lastspike, const STDPRSNNSpikePairingHomCommonProperties &);

  // overloaded for all supported event types
  using Connection::check_event;
  void check_event(SpikeEvent&) {}

 private:

  double_t facilitate_(double_t w, double_t kplus, const STDPRSNNSpikePairingHomCommonProperties &cp);
  double_t depress_(double_t w, double_t kminus, const STDPRSNNSpikePairingHomCommonProperties &cp);

  };


inline
double_t STDPRSNNSpikePairingConnectionHom::facilitate_(double_t w, double_t kplus, const STDPRSNNSpikePairingHomCommonProperties &cp)
{
  double_t norm_w = (w / cp.Wmax_) + (cp.lambda_ * std::pow(1.0 - (w/cp.Wmax_), cp.mu_plus_) * kplus);
  return norm_w < 1.0 ? norm_w * cp.Wmax_ : cp.Wmax_;
}

inline
double_t STDPRSNNSpikePairingConnectionHom::depress_(double_t w, double_t kminus, const STDPRSNNSpikePairingHomCommonProperties &cp)
{
  double_t norm_w = (w / cp.Wmax_) - (cp.alpha_ * cp.lambda_ * std::pow(w/cp.Wmax_, cp.mu_minus_) * kminus);
  return norm_w > 0.0 ? norm_w * cp.Wmax_ : 0.0;
}


inline
  void STDPRSNNSpikePairingConnectionHom::check_connection(Node & s, Node & r, port receptor_type, double_t t_lastspike)
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
void STDPRSNNSpikePairingConnectionHom::send(Event& e, double_t t_lastspike, const STDPRSNNSpikePairingHomCommonProperties &cp)
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
  double_t minus_dt = 0;
  double_t plus_dt = 0;

  if(start != finish) //take only first postspike after last prespike
  {
    minus_dt = t_lastspike - (start->t_ + dendritic_delay);
  }

  if(start != finish) //take only last postspike before current spike
  {
    --finish;
    plus_dt = (finish->t_ + dendritic_delay) - t_spike;
  }

  if(minus_dt != 0){
    weight_ = facilitate_(weight_, std::exp(minus_dt / cp.tau_plus_), cp);
  }

  if(plus_dt != 0){
    weight_ = depress_(weight_, std::exp(plus_dt / cp.tau_minus_), cp);
  }

  e.set_receiver(*target_);
  e.set_weight(weight_);
  e.set_delay(delay_);
  e.set_rport(rport_);
  e();

  }
} // of namespace nest

#endif // of #ifndef STDP_CONNECTION_HOM_H
