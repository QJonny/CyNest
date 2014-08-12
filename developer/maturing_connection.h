/*
 *  maturing_connection.h
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
 * date: april 2007
 */


#ifndef MATURING_CONNECTION_H
#define MATURING_CONNECTION_H

// use this for additional diagnostic evaluations
//#define DEBUG_MATURATION

#include "connection_het_wd.h"
#include "archiving_node.h"


namespace nest {

/**
 * Class storing properties which are common to all maturing connections.
 */
class MaturingCommonProperties : public CommonSynapseProperties {

  // grant access to private members for MaturingConnection  
friend class MaturingConnection;
friend class MaturingConnector;

 public:
  /**
   * Default constructor.
   * Sets all property values to defaults.
   */
  MaturingCommonProperties();
   
  /**
   * Get all properties and put them into a dictionary.
   */
  void get_status(DictionaryDatum & d) const;
  
  /**
   * Set properties from the values given in dictionary.
   */
  void set_status(const DictionaryDatum & d, ConnectorModel& cm);
 
  /**
   * Increment counter of died synapses.
   */
  void inc_synapse_died();

  /**
   * Increment counter of matured synapses.
   */
  void inc_synapse_matured();

  // overloaded for all supported event types
  void check_event(SpikeEvent&) {}

 private:
  double_t theta_b_;   //* lowest threshold on calcium influx, nothing happens below
  double_t theta_l_;   //* lower threshold on calcium influx, below negative counting
  double_t theta_h_;   //* higher threshold on calcium influx, above positive counting
  double_t theta_m0_;  //* threshold for maturation
  double_t theta_m1_;  //* threshold for dieing
  double_t tau_nmda_;  //* time constant of NMDA receptor deactivation
  double_t t_d0_;      //* time, after which decision premature -> mature/died is being made
  double_t t_d1_;      //* time, after which mature -> died is being made
  double_t tau_rec_;   //* recovery time of synaptic efficacy. 0.0 -> no depression
  double_t u_;         //* probability of quantal release for synaptic depression
  bool ignore_multapses_; //* whether multiple connections should be ignored upon connect
  long n_died_; //* number of synapses died
  long n_matured_; //* number of matured synapses
};

inline 
void MaturingCommonProperties::inc_synapse_died()
{
  n_died_++;
}


inline
void MaturingCommonProperties::inc_synapse_matured()
{
  n_matured_++;
}



/**
 *   Class storing information for one maturing connection.
 */
class MaturingConnection : public ConnectionHetWD
{

  public:
  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  MaturingConnection();
  
  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  MaturingConnection(const MaturingConnection &);

  /**
   * Default Destructor.
   */
  ~MaturingConnection() {}
  
  /*
   * This function calls check_connection on the sender and checks if the receiver
   * accepts the event type and receptor type requested by the sender.
   * Node::check_connection() will either confirm the receiver port by returning
   * true or false if the connection should be ignored.
   * We have to override the base class' implementation, since for maturing
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
   * Dummy implementation. Needed to be able to use GenericConnectorBase.
   */
  inline
  void send(Event&, double_t, MaturingCommonProperties &)
  { }

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param t_lastspike Point in time of last spike sent.
   * \param cp common properties of all synapses.
   * \param h activation state of NMDA receptors
   * \param x synaptic resources in readily releasable pool
   * \return true, if synapse has died
   */
  inline
  bool send(Event& e, double_t t_lastspike, MaturingCommonProperties &cp, double_t h, double_t x);

 private: 
  
  // dynamic variables stored for each connection
  
  bool mature_; // state of synapse {premature, mature}
  long m_;      // correlation random walk variable
  long n_post_; // count of postsynaptic spikes since birth of synapse
  double_t t0_; // time when synapse was created, -1, if mature synapse
#ifdef DEBUG_MATURATION
  long n_plus_; //* number of positive jumps
  long n_minus_; //* number of negative jumps
  long n_pre_; //* number of presynaptic spikes
#endif

};

inline 
  void MaturingConnection::check_connection(Node & s, Node & r, port receptor_type, double_t t_lastspike)
{
  ConnectionHetWD::check_connection(s, r, receptor_type, t_lastspike);
  r.register_stdp_connection(t_lastspike - Time(Time::step(delay_)).get_ms());
}

inline
bool MaturingConnection::send(Event& e, double_t t_lastspike, MaturingCommonProperties &cp, double_t h0, double_t x)
{
  //std::cout << "MaturingConnection::send" << std::endl;

#ifdef DEBUG_MATURATION
	n_pre_++;
#endif

  
  bool bDeleted = false;

  double_t t_spike = e.get_stamp().get_ms();

  //get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque<histentry>::iterator start;
  std::deque<histentry>::iterator finish;    

  double_t dendritic_delay = Time(Time::step(delay_)).get_ms();

  target_->get_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay,
		       &start, &finish);

  double_t td = cp.t_d0_;
  double_t theta_m = cp.theta_m0_;

  if (mature_) // synapse is mature
    {
      td = cp.t_d1_;
      theta_m = cp.theta_m1_;

      // so send spike actually
      e.set_receiver(*target_);
      e.set_weight(weight_ * x * cp.u_ );
      e.set_delay(delay_);
      e.set_rport(rport_);
      e();
    }
   
  // go through all postsynaptic spikes
  for ( ; start != finish; ++start )
    {
      // propagate NMDA activation state to point in time t = t_post + dendritic delay,
      // when the postsynaptic spike arrives at the synapse	
      double_t h = h0 * std::exp( -(start->t_ + dendritic_delay - t_lastspike) / cp.tau_nmda_ );

      if (h > cp.theta_b_ && h <= cp.theta_l_) {
#ifdef DEBUG_MATURATION
	n_minus_++;
#endif
	m_--;
      }
      else
	if (h > cp.theta_h_) {
#ifdef DEBUG_MATURATION
	  n_plus_++;
#endif
	  m_++;
	}

      //f_syn << m << ' ' << start->t_ + dendritic_delay << ' ' << h << '\n';
      
      n_post_++;
    }

  if (t_spike - t0_ > td)
    {      
      //std::cout << "MaturingConnection::send: reached checkpoint, m = " << m_ << std::endl;

      if ( m_ > theta_m /* std::sqrt(static_cast<double>(d->n_post)) */ ) // synapse is mature
	{	
	  if (!mature_)
	    cp.inc_synapse_matured(); // count matured synapses

	  mature_ = true; // mark as mature	
	  t0_ = t_spike;
	  m_ = 0;
	  n_post_ = 0;
	}        
      else 
	{
	  //std::cout << "synapse dies" << std::endl;
	  // unregister this neuron at the postsynaptic neuron
	  target_->unregister_stdp_connection(t_spike - dendritic_delay);

	  // report this synapse as deleted
	  bDeleted = true;
	
	  cp.inc_synapse_died();    // count died synapses
	}
    }

  // return, whether this synapse has died
  return bDeleted;
}




} // namespace nest

#endif // #ifndef MATURING_CONNECTION_H
