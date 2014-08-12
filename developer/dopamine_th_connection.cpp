
/*
 *  dopamine_th_connection.cpp
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2004 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.

 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#include "network.h"
#include "dictdatum.h"
#include "connector_model.h"
#include "common_synapse_properties.h"
#include "dopamine_th_connection.h"
#include "event.h"
#include "nestmodule.h"

namespace nest
{
  //
  // Implementation of class DOPAMINE_TH_CommonProperties.
  //

  DOPAMINE_TH_CommonProperties::DOPAMINE_TH_CommonProperties() :
    CommonSynapseProperties(),
    vt_(0),
    tau_d_(100.0),
    tau_post_(250.0),
    tau_epsilon_(1000.0),
    tau_pre_(300.0),
    td_alpha_(1.0),
    gamma_(0.9),
    dopa_base_(0.197907222305),
    c_(38.75),
    weight_min_(0.0),
    weight_max_(500.0),
    limit0_(0.0036),
    limit1_(0.0031),
    limit2_(0.001)
        
  { }

  void DOPAMINE_TH_CommonProperties::get_status(DictionaryDatum & d) const
  {
    CommonSynapseProperties::get_status(d);
    if(vt_!= 0)
      def<long_t>(d, "vt", vt_->get_gid());
    else def<long_t>(d, "vt", -1);
    def<double_t>(d, "tau_d", tau_d_);
    def<double_t>(d, "tau_post", tau_post_);
    def<double_t>(d, "tau_epsilon", tau_epsilon_);
    def<double_t>(d, "tau_pre", tau_pre_);
    def<double_t>(d, "td_alpha", td_alpha_);
    def<double_t>(d, "gamma", gamma_);
    def<double_t>(d, "dopa_base", dopa_base_);
    def<double_t>(d, "c", c_);
    def<double_t>(d, "weight_min", weight_min_);
    def<double_t>(d, "weight_max", weight_max_);
    def<double_t>(d, "limit0", limit0_);
    def<double_t>(d, "limit1", limit1_);
    def<double_t>(d, "limit2", limit2_);
    
  }
  
  void DOPAMINE_TH_CommonProperties::set_status(const DictionaryDatum & d, ConnectorModel &cm)
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
    updateValue<double_t>(d, "tau_post", tau_post_);
    updateValue<double_t>(d, "tau_epsilon", tau_epsilon_);
    updateValue<double_t>(d, "tau_pre", tau_pre_);
    updateValue<double_t>(d, "td_alpha", td_alpha_);
    updateValue<double_t>(d, "gamma", gamma_);
    updateValue<double_t>(d, "dopa_base", dopa_base_);
    updateValue<double_t>(d, "c", c_);
    updateValue<double_t>(d, "weight_min", weight_min_);
    updateValue<double_t>(d, "weight_max", weight_max_);
    updateValue<double_t>(d, "limit0", limit0_);
    updateValue<double_t>(d, "limit1", limit1_);
    updateValue<double_t>(d, "limit2", limit2_);
}


  Node* DOPAMINE_TH_CommonProperties::get_node()
   {
     if(vt_==0)
       throw BadProperty("No volume transmitter has been assigned to the dopamine synapse.");
     else
       return vt_;
   }

  //
  // Implementation of class STDPConnectionHom.
  //

   DOPAMINE_TH_Connection::DOPAMINE_TH_Connection() :
     last_update_(0),
     last_dopa_spike_(0),
     dopa_trace_(0.15365),
     pre_k_(0.0),
     last_spike_(0.0),
     pre_state_('i')
     
   {
   }

  DOPAMINE_TH_Connection::DOPAMINE_TH_Connection(const DOPAMINE_TH_Connection &rhs) :
    ConnectionHetWD(rhs)
  {
     
     last_update_ = rhs.last_update_;
     last_dopa_spike_ = rhs.last_dopa_spike_;
     dopa_trace_ = rhs.dopa_trace_;
     pre_k_ = rhs.pre_k_;
     last_spike_ = rhs.last_spike_;
     pre_state_ = rhs.pre_state_; 
  }

  void DOPAMINE_TH_Connection::get_status(DictionaryDatum & d) const
   {

    // base class properties, different for individual synapse
   ConnectionHetWD::get_status(d);
     
   // own properties, different for individual synapse
   //  def<double_t>(d, "td_alpha", td_alpha_);
   // def<double_t>(d, "tau_post", tau_post_);
   //def<double_t>(d, "gamma", gamma_);
   //def<double_t>(d, "dopa_base", dopa_base_);
   }
  
  void DOPAMINE_TH_Connection::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    // base class properties
    ConnectionHetWD::set_status(d, cm);
     
    // updateValue<double_t>(d, "td_alpha", td_alpha_);
    //updateValue<double_t>(d, "tau_post", tau_post_);
    //updateValue<double_t>(d, "gamma", gamma_);
    //updateValue<double_t>(d, "dopa_base", dopa_base_);    
  }

   /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void DOPAMINE_TH_Connection::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);
    //   set_property<double_t>(d, "td_alphas", p, td_alpha_);
    //set_property<double_t>(d, "tau_posts", p, tau_post_);
    //set_property<double_t>(d, "gammas", p, gamma_);
    //set_property<double_t>(d, "dopa_bases", p, dopa_base_);
  }
     
  //}

   void DOPAMINE_TH_Connection::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);
    //    initialize_property(d, "td_alphas");
    //initialize_property(d, "tau_posts");
    //initialize_property(d, "gammas");
    //initialize_property(d, "dopa_bases");
  }
  
  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void DOPAMINE_TH_Connection::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);
    // append_property<double_t>(d, "td_alphas", td_alpha_);
    //append_property<double_t>(d, "tau_posts", tau_post_);
    //append_property<double_t>(d, "gammas", gamma_);
    //append_property<double_t>(d, "dopa_bases", dopa_base_);
  }

} // of namespace nest
