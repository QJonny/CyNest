/*
 *  markram_connection.cpp
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

#include "markram_connection.h"
#include "network.h"
#include "connector_model.h"

namespace nest
{

  MarkramConnection::MarkramConnection() :
    ConnectionHetWD(),
    tau_fac_(0.0),
    tau_rec_(800.0),
    U_(0.5),
    R_(1.0),
    u_(0.5)
  { }

  void MarkramConnection::get_status(DictionaryDatum & d) const
  {
    ConnectionHetWD::get_status(d);

    def<double_t>(d, "U", U_);
    def<double_t>(d, "tau_rec", tau_rec_);
    def<double_t>(d, "tau_fac", tau_fac_);
    def<double_t>(d, "R", R_);
    def<double_t>(d, "u", u_);
  }
  
  void MarkramConnection::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, cm);

    updateValue<double_t>(d, "U", U_);
    updateValue<double_t>(d, "tau_rec", tau_rec_);
    updateValue<double_t>(d, "tau_fac", tau_fac_);
    updateValue<double_t>(d, "R", R_);
    updateValue<double_t>(d, "u", u_);

    if (R_ > 1.0 || R_ < 0.0) {
       cm.network().message(SLIInterpreter::M_ERROR, "MarkramConnection::set_status()", "require 0<=R<=1.0");
    }

    if (u_ > 1.0 || u_ < 0.0) {
       cm.network().message(SLIInterpreter::M_ERROR, "MarkramConnection::set_status()", "require 0<=u<=1.0");
    }
    if (U_ > 1.0 || U_ < 0.0) {
       cm.network().message(SLIInterpreter::M_ERROR, "MarkramConnection::set_status()", "require 0<=U<=1.0");
    }


  }

  /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void MarkramConnection::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);

    set_property<double_t>(d, "Us", p, U_);
    set_property<double_t>(d, "tau_facs", p, tau_fac_);
    set_property<double_t>(d, "tau_recs", p, tau_rec_);
    set_property<double_t>(d, "Rs", p, R_);
    set_property<double_t>(d, "us", p, u_);
  }

  void MarkramConnection::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);

    initialize_property_array(d, "Us"); 
    initialize_property_array(d, "tau_facs"); 
    initialize_property_array(d, "tau_recs"); 
    initialize_property_array(d, "Rs");
    initialize_property_array(d, "us");
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void MarkramConnection::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);

    append_property<double_t>(d, "Us", U_); 
    append_property<double_t>(d, "tau_facs", tau_fac_);
    append_property<double_t>(d, "tau_recs", tau_rec_);  
    append_property<double_t>(d, "Rs", R_); 
    append_property<double_t>(d, "us", u_);
  }

} // of namespace nest
