
/*
 *  maturing_connection.cpp
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

#include "maturing_connection.h"
#include "network.h"
#include "nestmodule.h"

namespace nest {

  ////////////////////////////////////////////////////
  // Implementation of class MaturingCommonProperties
  ////////////////////////////////////////////////////

  MaturingCommonProperties::MaturingCommonProperties() :
    CommonSynapseProperties(),
    theta_b_(0.1),
    theta_l_(0.45),
    theta_h_(0.8),
    theta_m0_(1.0),
    theta_m1_(1.0),
    tau_nmda_(35.0),
    t_d0_(30000.0),
    t_d1_(30000.0),
    tau_rec_(0.0),
    u_(1.0),
    ignore_multapses_(false),
    n_died_(0),
    n_matured_(0)
  { }

  void MaturingCommonProperties::get_status(DictionaryDatum & d) const
  {
    CommonSynapseProperties::get_status(d);

    // properties common to all synapses
    def<double_t>(d, "t_d0", t_d0_);
    def<double_t>(d, "t_d1", t_d1_);
    def<double_t>(d, "tau_nmda", tau_nmda_);
    def<double_t>(d, "theta_b", theta_b_);
    def<double_t>(d, "theta_l", theta_l_);
    def<double_t>(d, "theta_h", theta_h_);
    def<double_t>(d, "theta_m0", theta_m0_);
    def<double_t>(d, "theta_m1", theta_m1_);
    def<double_t>(d, "tau_rec", tau_rec_);
    def<double_t>(d, "u", u_);
    def<bool>(d, "ignore_multapses", ignore_multapses_);
    def<long>(d, "n_died", n_died_);
    def<long>(d, "n_matured", n_matured_);
  }

  void MaturingCommonProperties::set_status(const DictionaryDatum & d, ConnectorModel& cm)
  {  
    CommonSynapseProperties::set_status(d, cm);
    updateValue<double_t>(d, "t_d0", t_d0_);
    updateValue<double_t>(d, "t_d1", t_d1_);
    updateValue<double_t>(d, "tau_nmda", tau_nmda_);
    updateValue<double_t>(d, "theta_b", theta_b_);
    updateValue<double_t>(d, "theta_l", theta_l_);
    updateValue<double_t>(d, "theta_h", theta_h_);
    updateValue<double_t>(d, "theta_m0", theta_m0_);  
    updateValue<double_t>(d, "theta_m1", theta_m1_);
    updateValue<double_t>(d, "tau_rec", tau_rec_);
    updateValue<double_t>(d, "u", u_);
    updateValue<bool>(d, "ignore_multapses", ignore_multapses_);
    updateValue<long>(d, "n_died", n_died_);
    updateValue<long>(d, "n_matured", n_matured_);
  }


  //////////////////////////////////////////////
  // Implementation of class MaturingConnection
  //////////////////////////////////////////////

  MaturingConnection::MaturingConnection() :
    ConnectionHetWD(),
    mature_(false),
    m_(0),
    n_post_(0),
    t0_(0.0)
  {  
#ifdef DEBUG_MATURATION
    n_plus_ = 0;
    n_minus_ = 0;
    n_pre_ = 0;
#endif
}

  
  MaturingConnection::MaturingConnection(const MaturingConnection &mc) :
    ConnectionHetWD(mc),
    mature_(mc.mature_),
    m_(mc.m_),
    n_post_(mc.n_post_)
  {
    // set time of birth for this synapse to current network time
    t0_ = NestModule::get_network().get_slice_origin().get_ms();
#ifdef DEBUG_MATURATION
    n_plus_ = 0;
    n_minus_ = 0;
    n_pre_ = 0;
#endif
  }
 

  void MaturingConnection::get_status(DictionaryDatum & d) const
  {
    ConnectionHetWD::get_status(d);

    def<bool>(d, "mature", mature_);
    def<long>(d, "m", m_);
    def<double_t>(d, "t0", t0_);
    def<long>(d, "n_post", n_post_);
#ifdef DEBUG_MATURATION
    def<long>(d, "n_plus", n_plus_);
    def<long>(d, "n_minus", n_minus_);
    def<long>(d, "n_pre", n_pre_);
#endif
  }
  
  void MaturingConnection::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, cm);

    updateValue<bool>(d, "mature", mature_);
    updateValue<long>(d, "m", m_);
    updateValue<double_t>(d, "t0", t0_);  
    updateValue<long>(d, "n_post", n_post_);
  }
  
  void MaturingConnection::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);

    set_property<bool>(d, "matures", p, mature_);
    set_property<long>(d, "ms", p, m_);
    set_property<double_t>(d, "t0s", p, t0_);
    set_property<long>(d, "n_posts", p, n_post_);
  }

  void MaturingConnection::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);

    initialize_property_array(d, "matures");
    initialize_property_array(d, "ms");
    initialize_property_array(d, "t0s");
    initialize_property_array(d, "n_posts");
#ifdef DEBUG_MATURATION
    initialize_property_array(d, "n_pluss");
    initialize_property_array(d, "n_minuss");
    initialize_property_array(d, "n_pres");
#endif
  }

  void MaturingConnection::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);

    append_property<bool>(d, "matures", mature_);
    append_property<long>(d, "ms", m_);
    append_property<double_t>(d, "t0s", t0_);
    append_property<long>(d, "n_posts", n_post_);
#ifdef DEBUG_MATURATION
    append_property<long>(d, "n_pluss", n_plus_);
    append_property<long>(d, "n_minuss", n_minus_);
    append_property<long>(d, "n_pres", n_pre_);
#endif
  }


} // namespace nest
