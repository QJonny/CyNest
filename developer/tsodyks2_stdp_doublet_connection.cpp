/*
 *  tsodyks_connection.cpp
 *
 *  This file is part of NEST
 *
 *  Edited by Giuseppe Chindemi
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

#include "tsodyks2_stdp_doublet_connection.h"
#include "network.h"
#include "connector_model.h"
#include "nest_names.h"

namespace nest
{

  Tsodyks2STDPDoubletConnection::Tsodyks2STDPDoubletConnection() :
    ConnectionHetWD(),
    U_(0.5),
    u_(U_),
    x_(U_),
    tau_rec_(800.0),
    tau_fac_(0.0),

    // STDP
    tau_plus_(20.0),
    lambda_(0.01),
    alpha_(1.0),
    mu_plus_(1.0),
    mu_minus_(1.0),
    Wmax_(100.0),
    Kplus_(0.0)
  {}

  void Tsodyks2STDPDoubletConnection::get_status(DictionaryDatum & d) const
  {
    ConnectionHetWD::get_status(d);

    def<double_t>(d, names::dU, U_);
    def<double_t>(d, names::u, u_);
    def<double_t>(d, names::tau_rec, tau_rec_);
    def<double_t>(d, names::tau_fac, tau_fac_);
    def<double_t>(d, names::x, x_);

    // STDP
    // TODO Add partameters to names.h
    def<double_t>(d, "tau_plus", tau_plus_);
    def<double_t>(d, "lambda", lambda_);
    def<double_t>(d, "alpha", alpha_);
    def<double_t>(d, "mu_plus", mu_plus_);
    def<double_t>(d, "mu_minus", mu_minus_);
    def<double_t>(d, "Wmax", Wmax_);
  }
  
  void Tsodyks2STDPDoubletConnection::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, cm);
    
    updateValue<double_t>(d, names::dU, U_);
    updateValue<double_t>(d, names::u, u_);
    updateValue<double_t>(d, names::tau_rec, tau_rec_);
    updateValue<double_t>(d, names::tau_fac, tau_fac_);
    updateValue<double_t>(d, names::x, x_);

    // STDP
    updateValue<double_t>(d, "tau_plus", tau_plus_);
    updateValue<double_t>(d, "lambda", lambda_);
    updateValue<double_t>(d, "alpha", alpha_);
    updateValue<double_t>(d, "mu_plus", mu_plus_);
    updateValue<double_t>(d, "mu_minus", mu_minus_);
    updateValue<double_t>(d, "Wmax", Wmax_);
  }

  /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void Tsodyks2STDPDoubletConnection::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);

    set_property<double_t>(d, names::dUs, p, U_);
    set_property<double_t>(d, names::us, p, u_);
    set_property<double_t>(d, names::xs, p, x_);
    set_property<double_t>(d, names::tau_recs, p, tau_rec_);
    set_property<double_t>(d, names::tau_facs, p, tau_fac_);

    // STDP
    set_property<double_t>(d, "tau_pluss", p, tau_plus_);
    set_property<double_t>(d, "lambdas", p, lambda_);
    set_property<double_t>(d, "alphas", p, alpha_);
    set_property<double_t>(d, "mu_pluss", p, mu_plus_);
    set_property<double_t>(d, "mu_minuss", p, mu_minus_);
    set_property<double_t>(d, "Wmaxs", p, Wmax_);
  }

  void Tsodyks2STDPDoubletConnection::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);

    initialize_property_array(d, names::dUs); 
    initialize_property_array(d, names::us); 
    initialize_property_array(d, names::tau_recs);  
    initialize_property_array(d, names::tau_facs);  
    initialize_property_array(d, names::xs);

    // STDP
    initialize_property_array(d, "tau_pluss"); 
    initialize_property_array(d, "lambdas"); 
    initialize_property_array(d, "alphas"); 
    initialize_property_array(d, "mu_pluss"); 
    initialize_property_array(d, "mu_minuss");
    initialize_property_array(d, "Wmaxs");
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void Tsodyks2STDPDoubletConnection::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);

    append_property<double_t>(d, names::dUs, U_); 
    append_property<double_t>(d, names::us, u_); 
    append_property<double_t>(d, names::tau_recs, tau_rec_);  
    append_property<double_t>(d, names::tau_facs, tau_fac_);  
    append_property<double_t>(d, names::xs, x_); 

    // STDP
    append_property<double_t>(d, "tau_pluss", tau_plus_); 
    append_property<double_t>(d, "lambdas", lambda_); 
    append_property<double_t>(d, "alphas", alpha_); 
    append_property<double_t>(d, "mu_pluss", mu_plus_); 
    append_property<double_t>(d, "mu_minuss", mu_minus_);
    append_property<double_t>(d, "Wmaxs", Wmax_);
  }

} // of namespace nest
