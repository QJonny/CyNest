
/*
 *  stdp_triplet_connection.cpp
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
 *  Author: Eilif Muller (modified from stdp_synapse)
 *
 */
#include "network.h"
#include "dictdatum.h"
#include "connector_model.h"
#include "common_synapse_properties.h"
#include "stdp_triplet_connection.h"
#include "event.h"


namespace nest
{

  STDPTripletConnection::STDPTripletConnection() :
    StaticConnection(),
    tau_plus_(20.0),
    tau_x_(700.0),
    A_2p_(0.01),
    A_2m_(0.01),
    A_3p_(0.01),
    A_3m_(0.0),
    Kplus_(0.0),
    Kx_(0.0)
  { }


  STDPTripletConnection::STDPTripletConnection(const STDPTripletConnection &rhs) :
    StaticConnection(rhs)
  {
    tau_plus_ = rhs.tau_plus_;
    tau_x_ = rhs.tau_x_;
    A_2p_ = rhs.A_2p_;
    A_2m_ = rhs.A_2m_;
    A_3p_ = rhs.A_3p_;
    A_3m_ = rhs.A_3m_;
    Kplus_ = rhs.Kplus_;
    Kx_ = rhs.Kx_;
  }

  void STDPTripletConnection::get_status(DictionaryDatum & d) const
  {
    StaticConnection::get_status(d);

    def<double_t>(d, "tau_plus", tau_plus_);
    def<double_t>(d, "tau_x", tau_x_);
    def<double_t>(d, "A_2p", A_2p_);
    def<double_t>(d, "A_2m", A_2m_);
    def<double_t>(d, "A_3p", A_3p_);
    def<double_t>(d, "A_3m", A_3m_);
    def<double_t>(d, "Kplus", Kplus_);
    def<double_t>(d, "Kx", Kx_);
  }
  
  void STDPTripletConnection::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    StaticConnection::set_status(d, cm);

    updateValue<double_t>(d, "tau_plus", tau_plus_);
    updateValue<double_t>(d, "tau_x", tau_x_);
    updateValue<double_t>(d, "A_2p", A_2p_);
    updateValue<double_t>(d, "A_2m", A_2m_);
    updateValue<double_t>(d, "A_3p", A_3p_);
    updateValue<double_t>(d, "A_3m", A_3m_);
    updateValue<double_t>(d, "Kplus", Kplus_);    
    updateValue<double_t>(d, "Kx", Kx_);    
  }

   /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void STDPTripletConnection::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
    StaticConnection::set_status(d, p, cm);

    set_property<double_t>(d, "tau_pluss", p, tau_plus_);
    set_property<double_t>(d, "tau_xs", p, tau_x_);
    set_property<double_t>(d, "A_2ps", p, A_2p_);
    set_property<double_t>(d, "A_2ms", p, A_2m_);
    set_property<double_t>(d, "A_3ps", p, A_3p_);
    set_property<double_t>(d, "A_3ms", p, A_3m_);

    set_property<double_t>(d, "Kpluss", p, Kplus_);
    set_property<double_t>(d, "Kxs", p, Kx_);
  }


  void STDPTripletConnection::initialize_property_arrays(DictionaryDatum & d) const
  {
    StaticConnection::initialize_property_arrays(d);

    initialize_property_array(d, "tau_pluss");
    initialize_property_array(d, "tau_xs");
    initialize_property_array(d, "A_2ps");
    initialize_property_array(d, "A_2ms");

    initialize_property_array(d, "A_3ps");
    initialize_property_array(d, "A_3ms");
    initialize_property_array(d, "Kpluss");
    initialize_property_array(d, "Kxs");

  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void STDPTripletConnection::append_properties(DictionaryDatum & d) const
  {
    StaticConnection::append_properties(d);

    append_property<double_t>(d, "tau_pluss", tau_plus_);
    append_property<double_t>(d, "tau_xs", tau_x_);
    append_property<double_t>(d, "A_2ps", A_2p_);
    append_property<double_t>(d, "A_2ms", A_2m_);
    append_property<double_t>(d, "A_3ps", A_3p_);
    append_property<double_t>(d, "A_3ms", A_3m_);

    append_property<double_t>(d, "Kpluss", Kplus_);
    append_property<double_t>(d, "Kxs", Kx_);

  }

} // of namespace nest
