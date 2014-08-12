
/*
 *  maturing_connection_fr.cpp
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
 * date: October 2007
 */

#include "maturing_connection_fr.h"
#include "network.h"
#include "nestmodule.h"

namespace nest {

  ////////////////////////////////////////////////////
  // Implementation of class MaturingCommonProperties
  ////////////////////////////////////////////////////

  MaturingCommonPropertiesFr::MaturingCommonPropertiesFr() :
    CommonSynapseProperties(),
    theta_b_(0.1),
    theta_l_(0.45),
    theta_h_(0.8),
    theta_m0_(60),
    theta_m1_(40),
    tau_nmda_(35.0),
    N_camkii_(100),
    p_(0.01),
    q_(0.01),
    tau_rec_(0.0),
    u_(1.0),
    ignore_multapses_(false),
    n_died_(0),
    n_matured_(0),
    binomial_dev_(0.5, 10)
  { }

  void MaturingCommonPropertiesFr::get_status(DictionaryDatum & d) const
  {
    CommonSynapseProperties::get_status(d);

    // properties common to all synapses
    def<double_t>(d, "theta_b", theta_b_);
    def<double_t>(d, "theta_l", theta_l_);
    def<double_t>(d, "theta_h", theta_h_);
    def<long_t>(d, "theta_m0", theta_m0_);
    def<long_t>(d, "theta_m1", theta_m1_);
    def<double_t>(d, "tau_nmda", tau_nmda_);
    def<long_t>(d, "N_camkii", N_camkii_);
    def<double_t>(d, "p", p_);
    def<double_t>(d, "q", q_);
    def<double_t>(d, "tau_rec", tau_rec_);
    def<double_t>(d, "u", u_);
    def<bool>(d, "ignore_multapses", ignore_multapses_);
    def<long>(d, "n_died", n_died_);
    def<long>(d, "n_matured", n_matured_);
  }

  void MaturingCommonPropertiesFr::set_status(const DictionaryDatum & d, ConnectorModel& cm)
  {  
    CommonSynapseProperties::set_status(d, cm);
    updateValue<double_t>(d, "theta_b", theta_b_);
    updateValue<double_t>(d, "theta_l", theta_l_);
    updateValue<double_t>(d, "theta_h", theta_h_);
    updateValue<long_t>(d, "theta_m0", theta_m0_);  
    updateValue<long_t>(d, "theta_m1", theta_m1_);
    updateValue<double_t>(d, "tau_nmda", tau_nmda_);
    updateValue<long_t>(d, "N_camkii", N_camkii_);
    updateValue<double_t>(d, "p", p_);
    updateValue<double_t>(d, "q", q_);
    updateValue<double_t>(d, "tau_rec", tau_rec_);
    updateValue<double_t>(d, "u", u_);
    updateValue<bool>(d, "ignore_multapses", ignore_multapses_);
    updateValue<long>(d, "n_died", n_died_);
    updateValue<long>(d, "n_matured", n_matured_);
  }


  ////////////////////////////////////////////////
  // Implementation of class MaturingConnectionFr
  ////////////////////////////////////////////////

  MaturingConnectionFr::MaturingConnectionFr() :
    ConnectionHetWD(),
    mature_(false),
    m_(0),
    n_post_(0)
  {  
#ifdef DEBUG_MATURATION
    n_plus_ = 0;
    n_minus_ = 0;
    n_pre_ = 0;
    b_last_event_ = false;
#endif
}

  
  MaturingConnectionFr::MaturingConnectionFr(const MaturingConnectionFr &mc) :
    ConnectionHetWD(mc),
    mature_(mc.mature_),
    m_(mc.m_),
    n_post_(mc.n_post_)
  {
    // set time of birth for this synapse to current network time
    //t0_ = NestModule::get_network().get_slice_origin().get_ms();
#ifdef DEBUG_MATURATION
    n_plus_ = 0;
    n_minus_ = 0;
    n_pre_ = 0;
    b_last_event_ = false;
#endif
  }
 

  void MaturingConnectionFr::get_status(DictionaryDatum & d) const
  {
    ConnectionHetWD::get_status(d);

    def<bool>(d, "mature", mature_);
    def<long>(d, "m", m_);
    def<long>(d, "n_post", n_post_);
#ifdef DEBUG_MATURATION
    def<long>(d, "n_plus", n_plus_);
    def<long>(d, "n_minus", n_minus_);
    def<long>(d, "n_pre", n_pre_);
    def<bool>(d, "last_event", b_last_event_);
#endif
  }
  
  void MaturingConnectionFr::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, cm);

    updateValue<bool>(d, "mature", mature_);
    updateValue<long>(d, "m", m_);
    updateValue<long>(d, "n_post", n_post_);
#ifdef DEBUG_MATURATION
    updateValue<long>(d, "n_plus", n_plus_);
    updateValue<long>(d, "n_minus", n_minus_);
    updateValue<long>(d, "n_pre", n_pre_);
    updateValue<bool>(d, "last_event", b_last_event_);

#endif

  }
  
  void MaturingConnectionFr::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);

    set_property<bool>(d, "matures", p, mature_);
    set_property<long>(d, "ms", p, m_);
    set_property<long>(d, "n_posts", p, n_post_);
#ifdef DEBUG_MATURATION
    set_property<long>(d, "n_plus", p, n_plus_);
    set_property<long>(d, "n_minus", p, n_minus_);
    set_property<long>(d, "n_pre", p, n_pre_);
#endif
 }

  void MaturingConnectionFr::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);

    initialize_property_array(d, "matures");
    initialize_property_array(d, "ms");
    initialize_property_array(d, "n_posts");
#ifdef DEBUG_MATURATION
    initialize_property_array(d, "n_pluss");
    initialize_property_array(d, "n_minuss");
    initialize_property_array(d, "n_pres");
#endif
  }

  void MaturingConnectionFr::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);

    append_property<bool>(d, "matures", mature_);
    append_property<long>(d, "ms", m_);
    append_property<long>(d, "n_posts", n_post_);
#ifdef DEBUG_MATURATION
    append_property<long>(d, "n_pluss", n_plus_);
    append_property<long>(d, "n_minuss", n_minus_);
    append_property<long>(d, "n_pres", n_pre_);
#endif
  }


} // namespace nest
