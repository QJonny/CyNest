/*
 *  annealing_connection.h
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

#ifndef ANNEALING_CONNECTION_H
#define ANNEALING_CONNECTION_H

#include "connection_het_wd.h"
#include "normal_randomdev.h"
#include <iostream>

/* BeginDocumentation
  Name: annealing_synapse - Synapse type with short term plasticity.

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
     U          double - probability of release increment [0,1], default=0.5
     u          double - Maximum probability of release [0,1], default=0.5
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
 * Class representing a synapse with Tsodyks short term plasticity, based on the iterative formula
 * A suitable Connector containing these connections can be obtained from the template GenericConnector.
 */

namespace nest {
  /**
   * Class containing the common properties for all synapses of type ConnectionHomWD.
   */
  class AnnealingCommon : public CommonSynapseProperties
  {
    friend class AnnealingConnection;
    
  public:
    
    /**
     * Default constructor.
     * Sets all property values to defaults.
     */
    AnnealingCommon();
    
    /**
     * Get all properties and put them into a dictionary.
     */
    void get_status(DictionaryDatum & d) const;
    
    /**
     * Set properties from the values given in dictionary.
     */
    void set_status(const DictionaryDatum & d, ConnectorModel& cm);
    
    // data members common to all connections
    bool     with_noise;   //!< True if synapses are perturbed
    bool     update_means; //!< True if perturbed synapse parameters become new means.
    long_t   epoch;        //!< learning epoch counter
    int_t    mode;         //!< Noise mode for weight update
    double_t A_upper;
    double_t A_lower;
    double_t A_std;
    double_t U_upper;
    double_t U_lower;
    double_t U_std;
    double_t D_upper;
    double_t D_lower;
    double_t D_std;
    double_t F_upper;
    double_t F_lower;
    double_t F_std;
    librandom::NormalRandomDev normal_dev;  //!< random deviate generator
  };


class AnnealingConnection : public ConnectionHetWD
{
 public:

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  AnnealingConnection();
  /**
   * Copy constructor to propagate common properties.
   */
  AnnealingConnection(const AnnealingConnection&);

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
  void send(Event& e, double_t t_lastspike, AnnealingCommon &cp);

  // overloaded for all supported event types
  using Connection::check_event;
  void check_event(SpikeEvent&) {}
  
 private:
  double_t U_;       //!< unit increment of a facilitating synapse (U)
  double_t u_;       //!< dynamic value of probability of release
  double_t x_;       //!< current fraction of the synaptic weight 
  double_t tau_rec_; //!< [ms] time constant for recovery from depression (D)
  double_t tau_fac_; //!< [ms] time constant for facilitation (F)
  double_t A_;       //!< Noisy weight; we use the weight as the mean.

  // Variables for synapse annealing
  double_t U_mean_;
  double_t D_mean_;
  double_t F_mean_;
  long_t   epoch_; //!< Current learning epoch of the synapse

};


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 * \param cp Common properties object, containing the annealing parameters.
 */
inline
void AnnealingConnection::send(Event& e, double_t t_lastspike, AnnealingCommon &cp)
{
  double_t h = e.get_stamp().get_ms() - t_lastspike;  

  // has the synapse already received the noise or update_means signal?
  //std::cerr << "send: " << (size_t) &cp << '\n';
  //std::cerr << "epoch: "<< cp.epoch << '\n';

  if(epoch_ <= cp.epoch)
    {
      epoch_=cp.epoch; // mark synapse as updated
      //std::cerr << "epoch =" << epoch_ << std::endl;

      // If trial was successful, we update the state variables to new values. 
      if(cp.update_means)
	{
	  weight_=A_;
	  U_mean_=U_;
	  D_mean_= tau_rec_;
	  F_mean_=tau_fac_;
	}

      if(cp.with_noise)
	{ 
	  Network *net=Node::network();
	  int vp=target_->get_vp();
	  librandom::RngPtr rng=net->get_rng(vp);
	  
	  A_ = weight_ * (1. + cp.A_std* cp.normal_dev(rng));
	  U_ = U_mean_ *(1. + cp.U_std* cp.normal_dev(rng));
	  tau_rec_ = D_mean_* (1. + cp.D_std * cp.normal_dev(rng));
	  tau_fac_ = F_mean_*(1. + cp.F_std * cp.normal_dev(rng));

	  switch (cp.mode)
	    {
	    case -1: // Just allow negative changes
 	      A_ = std::min(A_, weight_);
	      break;
	    case 1:
	      A_ = std::max(A_, weight_);
	      break;
	    default:
	      break;
	    }

	  if(A_ < cp.A_lower)
	    A_=cp.A_lower;
	  if(A_> cp.A_upper)
	    A_=cp.A_upper;
	  
	  if(U_ < cp.U_lower)
	    U_=cp.U_lower;
	  else if (U_>cp.U_upper)
	    U_=cp.U_upper;
	  
	  if(tau_rec_< cp.D_lower)
	    tau_rec_=cp.D_lower;
	  else if (tau_rec_>cp.D_upper)
	    tau_rec_=cp.D_upper;
	  
	  if(F_mean_>0.0) // otherwise the synapse cannot facilitate
	    {
	      if(tau_fac_<cp.F_lower)
		tau_fac_=cp.F_lower;
	      else if (tau_fac_>cp.F_upper)
		tau_fac_=cp.F_upper;
	    }
	} else
	{
	  // Just assign the mean values
	  A_=weight_;
	  U_=U_mean_;
	  tau_rec_=D_mean_;
	  tau_fac_=F_mean_;
	}
    }

  const double_t f = std::exp(-h/tau_rec_);
  const double_t u_decay = (tau_fac_ < 1.0e-10) ? 0.0 : std::exp(-h/tau_fac_);
  
  x_= x_*(1.0-u_)*f + u_*(1.0-f); // Eq. 2 from reference [1]
  u_ *= u_decay; 
  u_+= U_*(1.0-u_); // for tau_fac=0 or u_=0, this will render u_==U_

  // send the spike to the target
  e.set_receiver(*target_);
  e.set_weight( x_*A_ );
  e.set_delay( delay_ );
  e.set_rport( rport_ );
  e();
}
 
} // namespace

#endif // ANNEALING_CONNECTION_H
