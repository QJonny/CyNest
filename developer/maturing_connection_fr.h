/*
 *  maturing_connection_fr.h
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

/* BeginDocumentation
  Name: maturing_synapse_fr - Synapse type for maturing and for pruning connections simulating a finite reservoir of
  CaMKII particcles, which are used as a correlation measure. 

  Description:
   maturing_synapse_fr measures the correlation between the presynaptic and the postsynaptic spike train. This is done
   by simulating the CaMKII dynamics in a spine, which is triggered by Ca influx. We assume the total amount of CaMKII
   to be finite.

  FirstVersion: October 2007
  Author: Moritz Helias

  SeeAlso: synapsedict, maturing_synapse
*/
 

#ifndef MATURING_CONNECTION_FR_H
#define MATURING_CONNECTION_FR_H

// use this for additional diagnostic evaluations
//#define DEBUG_MATURATION

#include "connection_het_wd.h"
#include "archiving_node.h"
#include "binomial_randomdev.h"

namespace nest {

/**
 * Class storing properties which are common to all maturing connections.
 */
class MaturingCommonPropertiesFr : public CommonSynapseProperties {

  // grant access to private members for MaturingConnectionFr  
friend class MaturingConnectionFr;
friend class MaturingConnectorFr;

 public:
  /**
   * Default constructor.
   * Sets all property values to defaults.
   */
  MaturingCommonPropertiesFr();
   
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
 inline 
 void inc_synapse_died();

  /**
   * Increment counter of matured synapses.
   */
 inline 
 void inc_synapse_matured();


 private:
  double_t theta_b_;   //* lowest threshold on calcium influx, nothing happens below
  double_t theta_l_;   //* lower threshold on calcium influx, below negative counting
  double_t theta_h_;   //* higher threshold on calcium influx, above positive counting
  long_t theta_m0_;  //* threshold for maturation, maturation, if x > theta_m0_
  long_t theta_m1_;  //* threshold for dieing, death, if x < theta_m1_
  double_t tau_nmda_;  //* time constant of NMDA receptor deactivation
  long_t N_camkii_;  //* total number of CaMKII particles
  double_t p_;         //* probability of a particle to be activated, given a high Ca influx
  double_t q_;         //* probability of a particle to be deactivated, given a low Ca influx
  double_t tau_rec_;   //* recovery time of synaptic efficacy. 0.0 -> no depression
  double_t u_;         //* probability of quantal release for synaptic depression
  bool ignore_multapses_; //* whether multiple connections should be ignored upon connect
  long n_died_; //* number of synapses died
  long n_matured_; //* number of matured synapses
  librandom::BinomialRandomDev binomial_dev_;  //!< random deviate generator for number of activated/deactivated particles
};

inline 
void MaturingCommonPropertiesFr::inc_synapse_died()
{
  n_died_++;
}


inline
void MaturingCommonPropertiesFr::inc_synapse_matured()
{
  n_matured_++;
}


/**
 *   Class storing information for one maturing connection.
 */
class MaturingConnectionFr : public ConnectionHetWD
{

  public:
  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  MaturingConnectionFr();
  
  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  MaturingConnectionFr(const MaturingConnectionFr &);

  /**
   * Default Destructor.
   */
  ~MaturingConnectionFr() {}
  
  // Import overloaded virtual function set to local scope. 
  using Connection::check_event;

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
   * \param t_lastspike last spike produced by presynaptic neuron
   */
  void check_connection(Node & s, Node & r, port receptor_type, double_t t_lastspike);

  /**
   * Overload for supported event types.
   */
  inline
  void check_event(SpikeEvent&)
  { }

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
  void send(Event&, double_t, MaturingCommonPropertiesFr &)
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
  bool send(Event& e, double_t t_lastspike, MaturingCommonPropertiesFr &cp, nest::Network &net, double_t h, double_t x);

 private: 
  
  // dynamic variables stored for each connection
  
  bool mature_; // state of synapse {premature, mature}
  long m_;      // correlation random walk variable
  long n_post_; // count of postsynaptic spikes since birth of synapse

#ifdef DEBUG_MATURATION
  long n_plus_; //* number of positive jumps
  long n_minus_; //* number of negative jumps
  long n_pre_; //* number of presynaptic spikes
  bool b_last_event_; //* indicates whether last event was plus of minus
#endif

};

inline 
void MaturingConnectionFr::check_connection(Node & s, Node & r, port receptor_type, double_t t_lastspike)
{
  ConnectionHetWD::check_connection(s, r, receptor_type, t_lastspike);
  r.register_stdp_connection(t_lastspike - Time(Time::step(delay_)).get_ms());
}

inline
bool MaturingConnectionFr::send(Event& e, double_t t_lastspike, MaturingCommonPropertiesFr &cp, nest::Network &net, double_t h0, double_t x)
{
  //std::cout << "MaturingConnectionFr::send" << std::endl;

#ifdef DEBUG_MATURATION
	n_pre_++;
#endif

  
  bool bDeleted = false;

  double_t t_spike = e.get_stamp().get_ms();

  //get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque<histentry>::iterator start;
  std::deque<histentry>::iterator finish;    

  double_t dendritic_delay = Time(Time::step(delay_)).get_ms();

  // for a new synapse, t_lastspike == 0.0, so
  // we initially read the whole history until including
  // t_spike-dendritic_delay, increasing the access counter for all
  // these entries MH 08-04-20
  target_->get_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay,
		       &start, &finish);

  if (mature_) // synapse is mature
    {
      // so send spike actually
      e.set_receiver(*target_);
      e.set_weight(weight_ * x * cp.u_ );
      e.set_delay(delay_);
      e.set_rport(rport_);
      e();
    }

  //std::cout.precision(12);
  //std::cout << "presynaptic spike at " << t_spike << " h = " << h0 << '\n';
  
  // do not take into account double presynaptic spikes
  if (t_lastspike == t_spike) return false;

  double_t t_last_post = 0.0;

  // go through all postsynaptic spikes
  for ( ; start != finish && !bDeleted; ++start )
    {
      // reject double spikes, since this it is unphysologic that they might both arrive back in the dendrite
      if (t_last_post == start->t_) continue;

      //std::cout << "postsynaptic spike at " << start->t_ << '\n';
      // propagate NMDA activation state to point in time t = t_post + dendritic delay,
      // when the postsynaptic spike arrives at the synapse	
      double_t h = h0 * std::exp( -(start->t_ + dendritic_delay - t_lastspike) / cp.tau_nmda_ );
      
      // draw a random number and make step stochastically
      //librandom::RngPtr rng = net.get_rng(target_->get_thread());
      //double_t r = rng->drand();
      //if (r < 0.249+0.249) { // do something at all
      //double r1 = rng->drand();
      
      //double p_pl = 0.44;
      //if (b_last_event_) // was a plus event
      //p_pl = 0.56;
      //else
      //	p_pl = 0.44;

      //      if (r1 < 1.0-p_pl) {
      if (h > cp.theta_b_ && h <= cp.theta_l_) {      
#ifdef DEBUG_MATURATION
	  //std::cout << "down, h =" << h << '\n';
	n_minus_++;
	b_last_event_ = false;
#endif

	cp.binomial_dev_.set_p_n(cp.q_, m_);
	librandom::RngPtr rng = net.get_rng(target_->get_thread());
	unsigned long m_minus = cp.binomial_dev_.uldev(rng);   //!< draw integer, threaded

	m_ -= m_minus;
      }
      else
	//{

	if (h > cp.theta_h_)
	  {
#ifdef DEBUG_MATURATION
	    //std::cout << "up, h=" << h << '\n';
	    n_plus_++;
	    b_last_event_ = true;
#endif

	    cp.binomial_dev_.set_p_n(cp.p_, cp.N_camkii_-m_);
	    librandom::RngPtr rng = net.get_rng(target_->get_thread());
	    unsigned long m_plus = cp.binomial_dev_.uldev(rng);   //!< draw integer, threaded

	    m_ += m_plus;
	  }
      //}
      //f_syn << m << ' ' << start->t_ + dendritic_delay << ' ' << h << '\n';
      
      n_post_++;
    

      if ( !mature_ && m_ > cp.theta_m0_ ) // synapse is mature
        {	
	  cp.inc_synapse_matured(); // count matured synapses

	  mature_ = true; // mark as mature	
	  //n_post_ = 0;
        }

      if ( mature_ && m_ < cp.theta_m1_ )
        {
	  //std::cout << "synapse dies" << std::endl;
	  // unregister this neuron at the postsynaptic neuron
	  // tell the archiver until where we have read the history
	  target_->unregister_stdp_connection(t_spike - dendritic_delay);

	  // report this synapse as deleted
	  bDeleted = true;
	
	  cp.inc_synapse_died();    // count died synapses
        }

      t_last_post = start->t_;
    }

  // return, whether this synapse has died
  return bDeleted;
}


} // namespace nest

#endif // #ifndef MATURING_CONNECTION_FR_H
