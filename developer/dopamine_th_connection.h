/*
 *  dopamine_th_connection.h
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

#ifndef DOPAMINE_TH_CONNECTION_H
#define DOPAMINE_TH_CONNECTION_H

/* BeginDocumentation
  Name: dopamine_th_synapse - Synapse type for dopamine modulated plasticity. The synapse implements the value function of actor-critic TD learning. 

  Description:
   dopamine_th_synapse is a connector to create synapses with dopamine modulated plasticity;
   the plasticity rule depends on pre- and postsynaptic spikes and on an additional volume signal (dopamine) coming from a non-local neuron pool; synapses are potentially active  if the presynaptic activity trace crosses a threshold, limit0; when the presynaptic activity trace decreases below a threshold, limit2, the potentially active synapses become active, as long as the postsynaptic activity trace decreases below a third threshold, limit1; this plasticity condition is described in Potjans W., Morrison A., Diesmann M.  Neural Computation 21, 301-339 (2009);
   the volume signal is made available to the synapse via the volume transmitter which collects all spikes from the neuron pool producing the volume signal; the dopamine synapse calculates based on this information a volume signal (low-pass filtered population rate); the plasticity rule is based on the value function update of actor-critic TD learning.
   
  
  
  Parameters:
   vt_      volume_transmitter from which the synapse receives information about volume signal modulating the plasticity


  FirstVersion: August 2009
  Author: Wiebke Potjans, Abigail Morrison
  
*/

#include "connection_het_wd.h"
#include "archiving_node.h"  
#include <cmath>
#include "volume_transmitter.h"

namespace nest
{

  /**
   * Class containing the common properties for all synapses of type DOPAMINE_TH_Connection.
   */
  class DOPAMINE_TH_CommonProperties : public CommonSynapseProperties
    {
      friend class DOPAMINE_TH_Connection;

    public:

      /**
       * Default constructor.
       * Sets all property values to defaults.
       */
      DOPAMINE_TH_CommonProperties();
   
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
      double_t tau_post_; 
      double_t tau_epsilon_; 
      double_t tau_pre_;
      double_t td_alpha_;
      double_t gamma_;
      double_t dopa_base_;
      double_t c_;
      double_t weight_min_;
      double_t weight_max_;
      double_t limit0_;
      double_t limit1_;
      double_t limit2_;
      
     
    };



  /**
   * Class representing an DOPAMINE_TH_ connection with homogeneous parameters, i.e. parameters are the same for all synapses.
   */
  class DOPAMINE_TH_Connection : public ConnectionHetWD
  {

  public:
  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  DOPAMINE_TH_Connection();
  
  /**
   * Copy constructor from a property object.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  DOPAMINE_TH_Connection(const DOPAMINE_TH_Connection &);

  /**
   * Default Destructor.
   */
  virtual ~DOPAMINE_TH_Connection() {}

  // Import overloaded virtual function set to local scope. 
  using Connection::check_event;

  /*
   * This function calls check_connection on the sender and checks if the receiver
   * accepts the event type and receptor type requested by the sender.
   * Node::check_connection() will either confirm the receiver port by returning
   * true or false if the connection should be ignored.
   * We have to override the base class' implementation, since for DOPAMINE
   * connections we have to call register_dopamine_connection on the target neuron
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
  void send(Event& e, double_t t_lastspike, const DOPAMINE_TH_CommonProperties &);


  // overloaded for all supported event types
  void check_event(SpikeEvent&) {}

  
  void trigger_update_weight(const vector<spikecounter> &dopa_spikes, const DOPAMINE_TH_CommonProperties &cp);
  void update_weight(double_t this_update, const DOPAMINE_TH_CommonProperties &cp);

  
  

  private:
 
 // data members of each connection
  double_t last_update_;
  double_t last_dopa_spike_;
  double_t dopa_trace_;
  double_t pre_k_;
  double_t last_spike_; 
  char pre_state_; 
 
  };




  inline
    void DOPAMINE_TH_Connection::trigger_update_weight(const vector<spikecounter> &dopa_spikes, const DOPAMINE_TH_CommonProperties &cp)
  {
    //get spike history of postsynaptic neuron in range (t1,t2]
    std::deque<histentry>::iterator start;
    std::deque<histentry>::iterator finish;
    target_->get_history(last_update_ , dopa_spikes.back().spike_time_, &start, &finish);
    for (uint_t i=0; i<dopa_spikes.size(); i++)
      {
	if(dopa_spikes[i].spike_time_ > last_update_) 
	  {
	    double_t this_dopa_spike = dopa_spikes[i].spike_time_;
	    while((start->t_ <  this_dopa_spike)&&(start != finish))
	      {
		update_weight(start->t_, cp);
		++start;
	      }
	    update_weight(this_dopa_spike, cp);
	    dopa_trace_ = dopa_trace_*std::exp((last_dopa_spike_-this_dopa_spike)/cp.tau_d_)+dopa_spikes[i].multiplicity_/cp.tau_d_;
	    last_dopa_spike_ = this_dopa_spike;
	  }
      }
  }
  

  inline
    void DOPAMINE_TH_Connection::update_weight(double_t this_update, const DOPAMINE_TH_CommonProperties &cp)
  {
 bool plastic = false;
      double_t pre_k_now = pre_k_*std::exp((last_spike_ - this_update) / cp.tau_pre_);

    //pre_state==pot_active
    if(pre_state_=='p'){
      if(pre_k_now <= cp.limit2_)//synapse is inactive
	{
	  plastic=true;
	}
      if((pre_k_now>cp.limit2_)&&(pre_k_now<cp.limit1_))//synapse is active
	{
	  plastic=true;
	}
    }

  //pre_state==active
    if(pre_state_=='a')
      {
	if(pre_k_now<=cp.limit2_){
	  plastic=true;
	}
	if((pre_k_now > cp.limit2_)&&(pre_k_now<cp.limit1_))
	  {
	    plastic=true;
	  }
      }

 //the actual state is defined 
  if(pre_state_=='i'){
    if(pre_k_now<=cp.limit2_) pre_state_ = 'i';
    if((pre_k_now>cp.limit2_)&&(pre_k_now<cp.limit0_)) pre_state_ = 'i';
    if(pre_k_now>=cp.limit0_) pre_state_ =  'p';
  }

  if(pre_state_=='p'){
    if(pre_k_now<=cp.limit2_) pre_state_ = 'i';
    if((pre_k_now>cp.limit2_)&&(pre_k_now<cp.limit1_)) pre_state_ = 'a';
    if(pre_k_now>=cp.limit1_) pre_state_= 'p';
  }

  if(pre_state_=='a'){
    if(pre_k_now<=cp.limit2_) pre_state_= 'i';
    if((pre_k_now>cp.limit2_)&&(pre_k_now<cp.limit1_)) pre_state_= 'a';
    if(pre_k_now>=cp.limit1_) pre_state_= 'p';
  }

 
      


  if(plastic)
    {
	double_t post_trace = target_ -> get_K_value(this_update);
	double_t exp_pre_t1 = std::exp((last_spike_ - last_update_)/cp.tau_pre_);
	double_t exp_pre_t2 = std::exp((last_spike_ - this_update)/cp.tau_pre_);
	double_t exp_dopa_t1 = std::exp((last_dopa_spike_ - last_update_)/cp.tau_d_);
	double_t exp_dopa_t2 = std::exp((last_dopa_spike_ - this_update)/cp.tau_d_);
	double_t exp_post_t1 = std::exp((this_update - last_update_)/cp.tau_post_);
	double_t exp_post_t2 = std::exp((this_update - this_update)/cp.tau_post_);
	double_t exp_epsilon_t1 = std::exp((last_spike_ - last_update_)/cp.tau_epsilon_); 
	double_t exp_epsilon_t2 = std::exp((last_spike_ - this_update)/cp.tau_epsilon_); 

	/*
	  weight_ = weight_ + (pre_k_*(dopa_trace_ - cp.dopa_base_)*(1.0/(1.0/cp.tau_d_+1.0/cp.tau_pre_)*(-exp_pre_t2*exp_dopa_t2 + exp_pre_t1*exp_dopa_t1))
	  + (cp.gamma_-1.0)*cp.c_*pre_k_*post_trace/cp.tau_post_/(1.0/cp.tau_pre_+1.0/cp.tau_epsilon_+1.0/cp.tau_post_)*(-exp_pre_t2*exp_epsilon_t2*exp_post_t2 + exp_pre_t1*exp_epsilon_t1*exp_post_t1))*cp.td_alpha_;
	*/
	
	weight_=weight_+cp.td_alpha_*pre_k_*
	  ((dopa_trace_ - cp.dopa_base_)*
	   (1.0/(1.0/cp.tau_pre_+1.0/cp.tau_epsilon_+1.0/cp.tau_d_)*
	    (exp_pre_t2*exp_epsilon_t2*exp_dopa_t2 - exp_pre_t1*exp_epsilon_t1*exp_dopa_t1)
	    - 1.0/(1.0/cp.tau_d_+1.0/cp.tau_pre_)*
	    (exp_pre_t2*exp_dopa_t2 - exp_pre_t1*exp_dopa_t1))
	   
	   
	   - (1.0-cp.gamma_)*cp.c_*post_trace/cp.tau_post_*
	   (1.0/(1.0/cp.tau_pre_+1.0/cp.tau_epsilon_+1.0/cp.tau_post_)*
	    (exp_pre_t2*exp_epsilon_t2*exp_post_t2 - exp_pre_t1*exp_epsilon_t1*exp_post_t1) 
	    - 1.0/(1.0/cp.tau_post_+1.0/cp.tau_pre_)*
	    (exp_pre_t2*exp_post_t2 - exp_pre_t1*exp_post_t1)));
    
    
    
    /*
    std::ofstream f("/home/potjans/td_dopa/baseline/baseline-30.txt",std::ios::app);
    f<<this_update<<"\t"<<dopa_trace_<<std::endl;
    f.close();
    */

  

    /*
    
    weight_=weight_;
    
    std::ostringstream ss;
    ss << "/home/potjans/td_dopa/baseline/baseline-" << weight_ << ".txt";
    std::ofstream f(ss.str().c_str(),std::ios::app);

    
    f<<this_update<<"\t"<<dopa_trace_<<std::endl; 
    f.close();
    */

  
	if(weight_ < cp.weight_min_) weight_ = cp.weight_min_;
	if(weight_ > cp.weight_max_) weight_ = cp.weight_max_;
	  }
	last_update_ = this_update;
      
  }



inline 
  void DOPAMINE_TH_Connection::check_connection(Node & s, Node & r, port receptor_type, double_t t_lastspike)
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
void DOPAMINE_TH_Connection::send(Event& e, double_t t_lastspike, const DOPAMINE_TH_CommonProperties &cp)
{
  double_t t_spike = e.get_stamp().get_ms();
  last_spike_ = t_lastspike; //we need last spike time for weight update 
  const vector<spikecounter>& dopa_spikes = cp.vt_->deliver_spikes();
  //get spike history of postsynaptic neuron in range (t1,t2]
  std::deque<histentry>::iterator start;
  std::deque<histentry>::iterator finish;
  target_->get_history(last_update_ , t_spike , &start, &finish);
  for (uint_t i=0; i<dopa_spikes.size(); i++)
    {
      if((dopa_spikes[i].spike_time_ > last_update_)&&(dopa_spikes[i].spike_time_ <= t_spike)) 
	{
	  double_t this_dopa_spike = dopa_spikes[i].spike_time_;
	 
	  while((start != finish)&&(start->t_ < this_dopa_spike))
	    {
	      update_weight(start->t_, cp);
	      ++start;
	    }
	 
	  update_weight(this_dopa_spike,  cp);
	  dopa_trace_ = dopa_trace_*std::exp((last_dopa_spike_-this_dopa_spike)/cp.tau_d_)+dopa_spikes[i].multiplicity_/cp.tau_d_;
	  last_dopa_spike_ = this_dopa_spike;
	}
    }
  
  while((start != finish)&&(start->t_ <  t_spike))
    {
      //attention, pre == post not considered !!!
      update_weight(start->t_,  cp);
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

#endif // of #ifndef DOPAMINE_TH_CONNECTION_H
