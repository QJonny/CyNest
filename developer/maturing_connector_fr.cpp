/*
 *  maturing_connector.cpp
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


#include "network.h"
#include <exception>
#include "nestmodule.h"
#include "maturing_connection_fr.h"
#include "maturing_connector_fr.h"

namespace nest
{

  MaturingConnectorFr::MaturingConnectorFr(GCMT &cm) :
    GenericConnectorBase< MaturingConnectionFr,
                          MaturingCommonPropertiesFr,
                          GCMT >(cm),
    h_(0.0),
    x_(1.0)
  { }


/**
 * Register a new connection at the sender side.
 */ 
void MaturingConnectorFr::register_connection(Node& s, Node& r)
{
  // only connect, if we do not ignore multapses
  // or if there is no connection to this target
  if ( !connector_model_.get_common_properties().ignore_multapses_ ||
       connection_exists(r) == invalid_port_ )
    {
      GenericConnectorBase< MaturingConnectionFr, MaturingCommonPropertiesFr, GCMT >::register_connection(s, r);
    } 
}

/**
 * Register a new connection at the sender side.
 * Use given weight and delay.
 */ 
void MaturingConnectorFr::register_connection(Node& s, Node& r, double_t w, double_t d)
{
  // only connect, if we do not ignore multapses
  // or if there is no connection to this target
  if ( !connector_model_.get_common_properties().ignore_multapses_ ||
       connection_exists(r) == invalid_port_ )
    {
      GenericConnectorBase< MaturingConnectionFr, MaturingCommonPropertiesFr, GCMT >::register_connection(s, r, w, d);
    } 
}

/**
 * Register a new connection at the sender side. 
 * Use given dictionary for parameters.
 */ 
void MaturingConnectorFr::register_connection(Node& s, Node& r, DictionaryDatum& d)
{
  // only connect, if we do not ignore multapses
  // or if there is no connection to this target
  if ( !connector_model_.get_common_properties().ignore_multapses_ ||
       connection_exists(r) == invalid_port_ )
    {
      GenericConnectorBase< MaturingConnectionFr, MaturingCommonPropertiesFr, GCMT >::register_connection(s, r, d);
    } 
}


void MaturingConnectorFr::send(Event& e)
{

  double_t t_spike = e.get_stamp().get_ms();

  //std::ofstream f_syn("syn_nest.dat", std::ios_base::app | std::ios_base::out);
 
  if (t_spike != t_lastspike_) // if spike is the same as before, do not send it again
    {
      // calculate recovery from synaptic depression      
      double_t tau_rec = connector_model_.get_common_properties().tau_rec_;
      if (tau_rec > 0.0)
	x_ = 1.0 + (x_ - 1.0)*std::exp(-(t_spike-t_lastspike_)/tau_rec);	      

      size_t i=0;
      for (ConnIter it = connections_.begin(); it != connections_.end(); i++)
	{
	  e.set_port(i);

	  bool synapse_died = it->send(e, t_lastspike_, connector_model_.get_common_properties(), connector_model_.network(), h_, x_);

	  if (synapse_died)
	    {
	      // delete this synapse
	      // deletion changes the port numbers,
	      // such that they stay in ascending order
	      it = connections_.erase(it);
	      //std::cout << "deleted synapse. remaining = " << connections_.size() << '\n';

	    }
	  else
	    ++it;
	}

      // propagate activation state of NMDA receptors to this spike time
      // if t_pre == t_post + dendritic_delay,
      // presynaptic spike will affect h AFTER postsynaptic spike
      h_ = h_ * std::exp(-(t_spike - t_lastspike_) / connector_model_.get_common_properties().tau_nmda_) + 1.0;

      // calculate synaptic depression      
      if (tau_rec > 0.0)
	x_ -= connector_model_.get_common_properties().u_ * x_;	

    }

  t_lastspike_ = t_spike;

}



} // namespace
