/*
 *  tsodyks_connection.cpp
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

#include "annealing_connection.h"
#include "network.h"
#include "connection.h"
#include "connector_model.h"
#include "nest_names.h"
#include "dictutils.h"

namespace nest
{
  
  AnnealingCommon::AnnealingCommon()
    :
    CommonSynapseProperties(),
    with_noise(true),
    update_means(false),
    epoch(0),
    mode(0),
    A_upper(1000.),//!< Rather arbitrary large number
    A_lower(0.01),
    A_std(0.01),
    U_upper(0.95),
    U_lower(0.05),
    U_std(0.01),
    D_upper(2000.0),
    D_lower(1.0),
    D_std(0.01),
    F_upper(3000.0),
    F_lower(1.0),
    F_std(0.01)
  {}
  
  
  void AnnealingCommon::get_status(DictionaryDatum & d) const
  {
    CommonSynapseProperties::get_status(d);
    
    def<bool>(d, names::with_noise, with_noise);
    def<bool>(d, "update_means", update_means);
    def<long>(d, "training_epoch", epoch);
    def<long>(d, "mode", mode);
    def<double_t>(d, names::A_upper, A_upper);
    def<double_t>(d, names::A_lower, A_lower);
    def<double_t>(d, names::A_std, A_std);
    def<double_t>(d, names::U_upper, U_upper);
    def<double_t>(d, names::U_lower, U_lower);
    def<double_t>(d, names::U_std, U_std);
    def<double_t>(d, names::D_upper, D_upper);
    def<double_t>(d, names::D_lower, D_lower);
    def<double_t>(d, names::D_std, D_std);
    def<double_t>(d, names::F_upper, F_upper);
    def<double_t>(d, names::F_lower, F_lower);
    def<double_t>(d, names::F_std, F_std); 
  }
  
  void AnnealingCommon::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    CommonSynapseProperties::set_status(d, cm);
    
    updateValue<bool>(d, names::with_noise, with_noise);
    updateValue<bool>(d, "update_means", update_means);
    updateValue<long>(d, "training_epoch", epoch);
    updateValue<long>(d, "mode", mode);
    updateValue<double_t>(d, names::A_upper, A_upper);
    updateValue<double_t>(d, names::A_lower, A_lower);
    updateValue<double_t>(d, names::A_std, A_std);
    updateValue<double_t>(d, names::U_upper, U_upper);
    updateValue<double_t>(d, names::U_lower, U_lower);
    updateValue<double_t>(d, names::U_std, U_std);
    updateValue<double_t>(d, names::D_upper, D_upper);
    updateValue<double_t>(d, names::D_lower, D_lower);
    updateValue<double_t>(d, names::D_std, D_std);
    updateValue<double_t>(d, names::F_upper, F_upper);
    updateValue<double_t>(d, names::F_lower, F_lower);
    updateValue<double_t>(d, names::F_std, F_std);   
  }
  
  AnnealingConnection::AnnealingConnection() :
    ConnectionHetWD(),
    U_(0.5),
    u_(U_),
    x_(U_),
    tau_rec_(800.0),
    tau_fac_(10.0),
    A_(weight_),
    U_mean_(U_),
    D_mean_(tau_rec_),
    F_mean_(tau_fac_),
    epoch_(0)
  {
  }

  AnnealingConnection::AnnealingConnection(const AnnealingConnection &rhs) :
    ConnectionHetWD(rhs),
    U_(rhs.U_),
    u_(rhs.u_),
    x_(rhs.x_),
    tau_rec_(rhs.tau_rec_),
    tau_fac_(rhs.tau_fac_),
    A_(rhs.A_),
    U_mean_(rhs.U_mean_),
    D_mean_(rhs.D_mean_),
    F_mean_(rhs.F_mean_),
    epoch_(rhs.epoch_)
  {
  }


  void AnnealingConnection::get_status(DictionaryDatum & d) const
  {
    ConnectionHetWD::get_status(d);
    def<double_t>(d, names::dU, U_);
    def<double_t>(d, names::u, u_);
    def<double_t>(d, names::tau_rec, tau_rec_);
    def<double_t>(d, names::tau_fac, tau_fac_);
    def<double_t>(d, names::x, x_);
    def<double_t>(d, names::A, A_); 
    def<double_t>(d, names::U_mean, U_mean_); 
    def<double_t>(d, names::D_mean, D_mean_); 
    def<double_t>(d, names::F_mean, F_mean_);
    def<long_t>(d, names::epoch, epoch_); 
  }
  
  void AnnealingConnection::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {

    ConnectionHetWD::set_status(d, cm);
    
    updateValue<double_t>(d, names::dU, U_);
    updateValue<double_t>(d, names::u, u_);
    updateValue<double_t>(d, names::tau_rec, tau_rec_);
    updateValue<double_t>(d, names::tau_fac, tau_fac_);
    updateValue<double_t>(d, names::x, x_);
    updateValue<double_t>(d, names::weight, weight_); 
    updateValue<double_t>(d, names::A, A_); 
    updateValue<double_t>(d, names::U_mean, U_mean_); 
    updateValue<double_t>(d, names::D_mean, D_mean_); 
    updateValue<double_t>(d, names::F_mean, F_mean_);
    updateValue<long_t>(d, names::epoch, epoch_); 
  }

  /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void AnnealingConnection::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);

    set_property<double_t>(d, names::dUs, p, U_);
    set_property<double_t>(d, names::us, p, u_);
    set_property<double_t>(d, names::xs, p, x_);
    set_property<double_t>(d, names::tau_recs, p, tau_rec_);
    set_property<double_t>(d, names::tau_facs, p, tau_fac_);
    set_property<double_t>(d, names::A,p, A_); 
    set_property<double_t>(d, names::U_mean, p, U_mean_); 
    set_property<double_t>(d, names::D_mean, p, D_mean_); 
    set_property<double_t>(d, names::F_mean, p, F_mean_);
    set_property<long_t>(d, names::epoch, p, epoch_); 
 }

  void AnnealingConnection::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);

    initialize_property_array(d, names::dUs); 
    initialize_property_array(d, names::us); 
    initialize_property_array(d, names::tau_recs);  
    initialize_property_array(d, names::tau_facs);  
    initialize_property_array(d, names::xs); 
    initialize_property_array(d, names::A); 
    initialize_property_array(d, names::U_mean); 
    initialize_property_array(d, names::D_mean); 
    initialize_property_array(d, names::F_mean);
    initialize_property_array(d, names::epoch); 
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void AnnealingConnection::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);

    append_property<double_t>(d, names::dUs, U_); 
    append_property<double_t>(d, names::us, u_); 
    append_property<double_t>(d, names::tau_recs, tau_rec_);  
    append_property<double_t>(d, names::tau_facs, tau_fac_);  
    append_property<double_t>(d, names::xs, x_); 
    append_property<double_t>(d, names::A, A_); 
    append_property<double_t>(d, names::U_mean, U_mean_); 
    append_property<double_t>(d, names::D_mean, D_mean_); 
    append_property<double_t>(d, names::F_mean, F_mean_);
    append_property<long_t>(d, names::epoch, epoch_); 
  }

} // of namespace nest
