/*
 *  tsodyks2_stdp_doublet_connection.h
 *
 *  This file is part of NEST
 *
 *  Edited by Giuseppe Chindemi
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

#ifndef TSODYKS2_STDP_DOUBLET_CONNECTION_H
#define TSODYKS2_STDP_DOUBLET_CONNECTION_H

#include "connection_het_wd.h"

/* BeginDocumentation
  Name: tsodyks2_stdp_doublet_synapse - Synapse type with short term plasticity and basic 
   spike timing dependent plasticity (STDP).

  Description:
   This synapse model implements synaptic short-term depression and short-term facilitation
   according to [1] and [2]. It solves Eq (2) from [1] and modulates U according to eq. (2) of [2].

   This connection merely scales the synaptic weight, based on the spike history and the 
   parameters of the kinetic model. Thus, it is suitable for all types of synaptic dynamics,
   that is current or conductance based.

   The parameter A_se from the publications is represented by the
   synaptic weight. The variable x in the synapse properties is the
   factor that scales the synaptic weight.

   Parameters: 
     The following parameters can be set in the status dictionary:
     U          double - probability of release increment (U1) [0,1], default=0.5
     u          double - Maximum probability of release (U_se) [0,1], default=0.5
     x          double - current scaling factor of the weight, default=U 
     tau_rec    double - time constant for depression in ms, default=800 ms
     tau_rec    double - time constant for facilitation in ms, default=0 (off)

  Notes:

     Under identical conditions, the tsodyks2_synapse produces
     slightly lower peak amplitudes than the tsodyks_synapse. However,
     the qualitative behavior is identical. The script
     test_tsodyks2_synapse.py in the examples compares the two synapse
     models.


  References:
   [1] Tsodyks, M. V., & Markram, H. (1997). The neural code between neocortical pyramidal neurons 
       depends on neurotransmitter release probability. PNAS, 94(2), 719-23.
   [2] Fuhrmann, G., Segev, I., Markram, H., & Tsodyks, M. V. (2002). Coding of temporal 
       information by activity-dependent synapses. Journal of neurophysiology, 87(1), 140-8.

  Transmits: SpikeEvent
       
  FirstVersion: October 2011
  Author: Marc-Oliver Gewaltig, based on tsodyks_synapse by Moritz Helias
  SeeAlso: tsodyks_synapse, synapsedict, stdp_synapse, static_synapse
*/


/**
 * Class representing a synapse with Tsodyks short term plasticity and basic spike timing dependent 
 * plasticity (STDP), based on the iterative formula. A suitable Connector containing these 
 * connections can be obtained from the template GenericConnector.
 */

namespace nest {

class Tsodyks2STDPDoubletConnection : public ConnectionHetWD
{
 public:

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  Tsodyks2STDPDoubletConnection();

  /**
   * Default Destructor.
   */
  ~Tsodyks2STDPDoubletConnection() {}

  /*
   * STDP ???
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
   * \param cp Common properties to all synapses (empty).
   */
  void send(Event& e, double_t t_lastspike, const CommonSynapseProperties &cp);

  // overloaded for all supported event types
  using Connection::check_event;
  void check_event(SpikeEvent&) {}
  
 private:
  double_t U_;       //!< unit increment of a facilitating synapse (U1)
  double_t u_;       //!< probability of release (Use)
  double_t x_;       //!< current fraction of the synaptic weight 
  double_t tau_rec_; //!< [ms] time constant for recovery
  double_t tau_fac_; //!< [ms] time constant for facilitation

  // STDP ---------------------------------------------------------------------
  double_t potentiate_(double_t w, double_t kplus);
  double_t depress_(double_t w, double_t kminus);

  // data members of each connection
  double_t tau_plus_;
  double_t lambda_;
  double_t alpha_;
  double_t mu_plus_;
  double_t mu_minus_;  
  double_t Wmax_;
  double_t Kplus_;
  //---------------------------------------------------------------------------
};


//*****************************************************************************
//* STDP BEGIN ****************************************************************
//*****************************************************************************
/**
 * Potentiation.
 * TODO Verify if the formula is compliant with the one used in the detailed column
 * TODO Define and set parameters
 */
inline
double_t Tsodyks2STDPDoubletConnection::potentiate_(double_t w, double_t kplus)
{
  double_t norm_w = (w / Wmax_) + (lambda_ * std::pow(1.0 - (w/Wmax_), mu_plus_) * kplus);
  return norm_w < 1.0 ? norm_w * Wmax_ : Wmax_;
}


/**
 * Depression.
 * TODO Verify if the formula is compliant with the one used in the detailed column
 * TODO Define and set parameters
 */
inline 
double_t Tsodyks2STDPDoubletConnection::depress_(double_t w, double_t kminus)
{
  double_t norm_w = (w / Wmax_) - (alpha_ * lambda_ * std::pow(w/Wmax_, mu_minus_) * kminus);
  return norm_w > 0.0 ? norm_w * Wmax_ : 0.0;
}


/**
 * ???
 * TODO What is it? How does it work?
 */
inline 
void Tsodyks2STDPDoubletConnection::check_connection(Node & s, Node & r, port receptor_type, double_t t_lastspike)
{
  ConnectionHetWD::check_connection(s, r, receptor_type, t_lastspike);

  // For a new synapse, t_lastspike contains the point in time of the last spike.
  // So we initially read the history(t_last_spike - dendritic_delay, ...,  T_spike-dendritic_delay]
  // which increases the access counter for these entries.
  // At registration, all entries' access counters of history[0, ..., t_last_spike - dendritic_delay] will be 
  // incremented by the following call to Archiving_Node::register_stdp_connection().
  // See bug #218 for details.
  r.register_stdp_connection(t_lastspike - Time(Time::step(delay_)).get_ms());
}
//*****************************************************************************
//* STDP END ******************************************************************
//*****************************************************************************

/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
inline
void Tsodyks2STDPDoubletConnection::send(Event& e, double_t t_lastspike, const CommonSynapseProperties &)
{
  double_t h = e.get_stamp().get_ms() - t_lastspike;  
  double_t f = std::exp(-h/tau_rec_);
  double_t u_decay = (tau_fac_ < 1.0e-10) ? 0.0 : std::exp(-h/tau_fac_);

  x_= x_*(1.0-u_)*f + u_*(1.0-f); // Eq. 2 from reference [1]
  u_ *= u_decay; 
  u_+= U_*(1.0-u_); // for tau_fac=0 and u_=0, this will render u_==U_

  // send the spike to the target
//  e.set_receiver(*target_);
//  e.set_weight( x_*weight_ );
//  e.set_delay( delay_ );
//  e.set_rport( rport_ );
//  e();

  // STDP-----------------------------------------------------------------------
  // TODO
  // synapse STDP depressing/potentiation dynamics

  double_t t_spike = e.get_stamp().get_ms();
  // t_lastspike_ = 0 initially
  double_t dendritic_delay = Time(Time::step(delay_)).get_ms();

  //get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque<histentry>::iterator start;
  std::deque<histentry>::iterator finish;

  // For a new synapse, t_lastspike contains the point in time of the last spike.
  // So we initially read the history(t_last_spike - dendritic_delay, ...,  T_spike-dendritic_delay]
  // which increases the access counter for these entries.
  // At registration, all entries' access counters of history[0, ..., t_last_spike - dendritic_delay] have been 
  // incremented by Archiving_Node::register_stdp_connection(). See bug #218 for details.
  target_->get_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay,
                         &start, &finish);
  //facilitation due to post-synaptic spikes since last pre-synaptic spike
  double_t minus_dt;
  while (start != finish)
  {
    minus_dt = t_lastspike - (start->t_ + dendritic_delay);
    ++start;
    if (minus_dt == 0)
      continue;
    weight_ = potentiate_(weight_, Kplus_ * std::exp(minus_dt / tau_plus_));
  }

  //depression due to new pre-synaptic spike
  weight_ = depress_(weight_, target_->get_K_value(t_spike - dendritic_delay));

  e.set_receiver(*target_);
  e.set_weight(weight_*x_); // TODO Check if reasonable, have a look to the detailed model.
  e.set_delay(delay_);
  e.set_rport(rport_);
  e();

 
  Kplus_ = Kplus_ * std::exp((t_lastspike - t_spike) / tau_plus_) + 1.0;
  // --------------------------------------------------------------------------
}
 
} // namespace

#endif // TSODYKS2_STDP_DOUBLET_CONNECTION_H
