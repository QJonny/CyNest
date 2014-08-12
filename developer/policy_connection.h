/*
 *  policy_connection.h
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

#ifndef POLICY_CONNECTION_H
#define POLICY_CONNECTION_H

/* BeginDocumentation
  Name: policy_synapse - Synapse type for dopamine modulated synaptic plasticity. The synapse implements the policy of actor-critic TD learning.   

  Description:
   policy_synapse is a connector to create synapses with dopamine modulated plasticity;  the plasticity rule depends on pre- and postsynaptic spikes and on an additional volume signal (dopamine) coming from a non-local neuron pool; the plasticity rule is based on the policy update of actor-critic TD learning.
   

  
  Parameters:
   vt_      volume_transmitter from which the synapse receives information about volume signal modulating the plasticity

 
  FirstVersion: August 2008
  Author: Wiebke Potjans
  
*/

#include "connection_het_wd.h"
#include "archiving_node.h"
#include <cmath>
#include "volume_transmitter.h"

namespace nest
{

  /**
   * Class containing the common properties for all synapses of type POLICYConnection.
   */
  class PolicyCommonProperties : public CommonSynapseProperties
    {
      friend class PolicyConnection;

    public:

      /**
       * Default constructor.
       * Sets all property values to defaults.
       */
      PolicyCommonProperties();
   
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

      Node* get_node();

 
    private:
      volume_transmitter * vt_;
      double_t tau_d_; 
      double_t tau_pre_;
      double_t tau_epsilon_;
      double_t tau_a_;
      double_t beta_;
      double_t weight_min_;
      double_t weight_max_;
      double_t dopa_base_;
    };



  /**
   * Class representing an POLICY connection with homogeneous parameters, i.e. parameters are the same for all synapses.
   */
  class PolicyConnection : public ConnectionHetWD
  {

  public:
  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  PolicyConnection();
  
  /**
   * Copy constructor from a property object.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  PolicyConnection(const PolicyConnection &);

  /**
   * Default Destructor.
   */
  virtual ~PolicyConnection() {}

  // Import overloaded virtual function set to local scope. 
  using Connection::check_event;

  /*
   * This function calls check_connection on the sender and checks if the receiver
   * accepts the event type and receptor type requested by the sender.
   * Node::check_connection() will either confirm the receiver port by returning
   * true or false if the connection should be ignored.
   * We have to override the base class' implementation, since for POLICY
   * connections we have to call register_policy_connection on the target neuron
   * to inform the Archiver to collect spikes for this connection.
   *
   * \param s The source node
   * \param r The target node
   * \param receptor_type The ID of the requested receptor type
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
  void send(Event& e, double_t t_lastspike, const PolicyCommonProperties &);


  // overloaded for all supported event types
  void check_event(SpikeEvent&) {}

  
  void trigger_update_weight(const vector<spikecounter> &dopa_spikes, const PolicyCommonProperties &cp);
  void update_weight(double_t this_update, const PolicyCommonProperties &cp);

  private:
   // data members of each connection
  double_t last_update_;
  double_t last_dopa_spike_;
  double_t dopa_trace_;
  double_t pre_k_;
  double_t last_spike_; 
 
  };



  inline
    void PolicyConnection::trigger_update_weight(const vector<spikecounter> &dopa_spikes, const PolicyCommonProperties &cp)
  {
    //get spike history of postsynaptic neuron in range (t1,t2]
    std::deque<histentry>::iterator start;
    std::deque<histentry>::iterator finish;
    target_->get_history(last_update_ , dopa_spikes.back().spike_time_, &start, &finish);
    for(uint_t i=0; i < dopa_spikes.size(); i++) //update from dopa spike to dopa spike
      {
	if(dopa_spikes.at(i).spike_time_ >last_update_)
	  {
	    double_t this_dopa_spike = dopa_spikes.at(i).spike_time_;
	    
	    while((start->t_ < this_dopa_spike)&&(start != finish)) //update from post spike to post spike
	      {
		update_weight(start->t_, cp);
		++start;
	      }
	    
	    update_weight(this_dopa_spike, cp);
	    
	    dopa_trace_ = dopa_trace_*std::exp((last_dopa_spike_-this_dopa_spike)/cp.tau_d_)+dopa_spikes.at(i).multiplicity_/cp.tau_d_; 
	    last_dopa_spike_ = this_dopa_spike;
	  }
      }
  }
  
 

  inline
    void PolicyConnection::update_weight(double_t this_update, const PolicyCommonProperties &cp)
    {
     
      double_t post_trace = target_->get_K_value(this_update);
      
      double_t exp_pre_t1 = std::exp((last_spike_ - last_update_)/cp.tau_pre_);
      double_t exp_pre_t2 = std::exp((last_spike_ - this_update)/cp.tau_pre_);
      double_t exp_dopa_t1 = std::exp((last_dopa_spike_ - last_update_)/cp.tau_d_);
      double_t exp_dopa_t2 = std::exp((last_dopa_spike_ - this_update)/cp.tau_d_);
      double_t exp_epsilon_t1 = std::exp((last_spike_ - last_update_)/cp.tau_epsilon_);
      double_t exp_epsilon_t2 = std::exp((last_spike_ - this_update)/cp.tau_epsilon_); 

      double_t exp_post_t1 = std::exp((this_update - last_update_)/cp.tau_a_);
      double_t exp_post_t2 = std::exp((this_update - this_update)/cp.tau_a_);

      
      weight_ = weight_ + (pre_k_*(dopa_trace_-cp.dopa_base_)*post_trace*(1.0/(1.0/cp.tau_pre_+1.0/cp.tau_epsilon_+1.0/cp.tau_d_+1.0/cp.tau_a_)*(exp_pre_t2*exp_epsilon_t2*exp_dopa_t2*exp_post_t2 - exp_pre_t1*exp_epsilon_t1*exp_dopa_t1*exp_post_t1) -  1.0/(1.0/cp.tau_d_+1.0/cp.tau_pre_+1.0/cp.tau_a_)*(exp_pre_t2*exp_dopa_t2*exp_post_t2 - exp_pre_t1*exp_dopa_t1*exp_post_t1)))*cp.beta_;
      
    
      /*
       std::ofstream f("/home/potjans/td_dopa/2_states/policy/weight_policy.txt",std::ios::app);
       f<<this_update<<"\t"<<dopa_trace_<<"\t"<<cp.dopa_base_<<std::endl;
       f.close();
      */

       if(weight_ < cp.weight_min_) weight_ = cp.weight_min_;
       if(weight_ > cp.weight_max_) weight_ = cp.weight_max_;
       
       
       last_update_ = this_update;
  }


inline 
  void PolicyConnection::check_connection(Node & s, Node & r, port receptor_type, double_t t_lastspike)
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
void PolicyConnection::send(Event& e, double_t t_lastspike, const PolicyCommonProperties &cp)
{
 
  double_t t_spike = e.get_stamp().get_ms();
  
  last_spike_ = t_lastspike; //we need last spike time for weight update 
  vector<spikecounter> dopa_spikes = cp.vt_->deliver_spikes();
   //get spike history of postsynaptic neuron in range (t1,t2]
  std::deque<histentry>::iterator start;
  std::deque<histentry>::iterator finish;
  target_->get_history(last_update_ , t_spike , &start, &finish);

  for (uint_t i=0; i < dopa_spikes.size(); i++)
    { 
      if((dopa_spikes.at(i).spike_time_ > last_update_)&&(dopa_spikes.at(i).spike_time_ <= t_spike)) 
	{
	  double_t this_dopa_spike = dopa_spikes.at(i).spike_time_;

	  while((start != finish)&&(start->t_ < this_dopa_spike))
	    {
	      update_weight(start->t_, cp);
	      ++start;
	    }
	  update_weight(this_dopa_spike, cp);
	  dopa_trace_ = dopa_trace_*std::exp((last_dopa_spike_-this_dopa_spike)/cp.tau_d_)+dopa_spikes.at(i).multiplicity_/cp.tau_d_;
	  last_dopa_spike_ = this_dopa_spike;
	}
    }


  while((start != finish)&&(start->t_ < t_spike))
    {
      update_weight(start->t_, cp);
      ++start;
    }
  update_weight(t_spike, cp);


  
 
  //update for presynaptic activity trace after weight update, as weight update require information about activity trace at last spike 
    pre_k_ = pre_k_*std::exp((t_lastspike - t_spike) / cp.tau_pre_) + 1.0/cp.tau_pre_;
    last_spike_ = t_spike;
 
  e.set_receiver(*target_);
  e.set_weight(weight_);
  e.set_delay(delay_);
  e.set_rport(rport_);
  e();

  
  
}

} // of namespace nest

#endif // of #ifndef POLICY_CONNECTION_H
