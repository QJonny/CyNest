
/*
 *  stdp_connection.cpp
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
#include "dictdatum.h"
#include "connector_model.h"
#include "common_synapse_properties.h"
#include "policy_connection.h"
#include "event.h"
#include "nestmodule.h"

namespace nest
{
  //
  // Implementation of class PolicyCommonProperties.
  //

  PolicyCommonProperties::PolicyCommonProperties() :
    CommonSynapseProperties(),
    vt_(0),
    tau_d_(100.0),
    tau_pre_(300.0),
    tau_epsilon_(1000.0),
    tau_a_(300.0),  
    beta_(0.7),
    weight_min_(30.0),
    weight_max_(1000.0),
    dopa_base_(0.197907222305)
    
  { }

  void PolicyCommonProperties::get_status(DictionaryDatum & d) const
  {
    CommonSynapseProperties::get_status(d);
    if(vt_!= 0)  
      def<long_t>(d, "vt", vt_->get_gid());
    else def<long_t>(d, "vt", -1);
    def<double_t>(d, "tau_d", tau_d_);
    def<double_t>(d, "tau_pre", tau_pre_);
    def<double_t>(d, "tau_epsilon", tau_epsilon_);
    def<double_t>(d, "tau_a", tau_a_);
    def<double_t>(d, "beta", beta_);
    def<double_t>(d, "weight_min", weight_min_);
    def<double_t>(d, "weight_max", weight_max_);
    def<double_t>(d, "dopa_base", dopa_base_);
   
    
  }
  
  void PolicyCommonProperties::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    CommonSynapseProperties::set_status(d, cm);

    

    long_t vtgid;
    if(updateValue<long_t>(d, "vt", vtgid))
      {
	vt_ = dynamic_cast<volume_transmitter *>(NestModule::get_network().get_node(vtgid));

	if(vt_==0)
	  throw BadProperty("Dopamine source must be volume transmitter");
      }
    updateValue<double_t>(d, "tau_d", tau_d_);
    updateValue<double_t>(d, "tau_pre", tau_pre_);
    updateValue<double_t>(d, "tau_epsilon", tau_epsilon_);
    updateValue<double_t>(d, "tau_a", tau_a_);
    updateValue<double_t>(d, "beta", beta_);
    updateValue<double_t>(d, "weight_min", weight_min_);
    updateValue<double_t>(d, "weight_max", weight_max_);
    updateValue<double_t>(d, "dopa_base", dopa_base_);

  }


  Node* PolicyCommonProperties::get_node()
   {
     if(vt_==0)
       throw BadProperty("No volume transmitter has been assigned to the dopamine synapse.");
     else
       return vt_;
   }

  
  //
  // Implementation of class PolicyConnection.
  //

   PolicyConnection::PolicyConnection() : 
     last_update_(0),
     last_dopa_spike_(0),
     dopa_trace_(0.15365),
     pre_k_(0.0),
     last_spike_(0.0)
     
    {
    }

  PolicyConnection::PolicyConnection(const PolicyConnection &rhs) :
    ConnectionHetWD(rhs)
  {
     
     last_update_ = rhs.last_update_;
     last_dopa_spike_ = rhs.last_dopa_spike_;
     dopa_trace_ = rhs.dopa_trace_;
     pre_k_ = rhs.pre_k_;
     last_spike_ = rhs.last_spike_;
    

  }

  void PolicyConnection::get_status(DictionaryDatum & d) const
   {

    // base class properties, different for individual synapse
   ConnectionHetWD::get_status(d);
 // own properties, different for individual synapse
   
    
   }
  
  void PolicyConnection::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    // base class properties
    ConnectionHetWD::set_status(d, cm);

//     if (d->known("tau_plus") ||
// 	d->known("lambd") ||
// 	d->known("alpha") ||
// 	d->known("mu_plus") ||
// 	d->known("mu_minus") ||
// 	d->known("Wmax") )
//       {
// 	cm.network().error("DOPAMINEConnection::set_status()", "you are trying to set common properties via an individual synapse.");
//       }

    
  }

   /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void PolicyConnection::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
 
    ConnectionHetWD::set_status(d, p, cm);
    
  }
     
  //}

   void PolicyConnection::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);
    
  }
  
  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void PolicyConnection::append_properties(DictionaryDatum & d) const
  {
  ConnectionHetWD::append_properties(d);
  
  }

} // of namespace nest
