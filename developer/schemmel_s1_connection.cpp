
/*
 *  schemmel_s1_connection.cpp
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
 *  A modification of stdp_synapse by Eilif Muller
 *
 */
#include "network.h"
#include "dictdatum.h"
#include "connector_model.h"
#include "common_synapse_properties.h"
#include "schemmel_s1_connection.h"
#include "event.h"


namespace nest
{

  SchemmelS1Connection::SchemmelS1Connection() :
    ConnectionHetWD(),
    tau_plus_(20.0),
    tau_minus_(20.0),
    A_plus_(1.0),
    A_minus_(1.0),
    tau_rec_(800.0),
    U_SE_(0.5),
    I_(0.0),
    mode_(0)
  { }


  SchemmelS1Connection::SchemmelS1Connection(const SchemmelS1Connection &rhs) :
    ConnectionHetWD(rhs)
  {
    tau_plus_ = rhs.tau_plus_;
    tau_minus_ = rhs.tau_minus_;
    A_plus_ = rhs.A_plus_;
    A_minus_ = rhs.A_minus_;
    tau_rec_ = rhs.tau_rec_;
    U_SE_ = rhs.U_SE_;
    I_ = rhs.I_;
    mode_ = rhs.mode_;
  }

  void SchemmelS1Connection::get_status(DictionaryDatum & d) const
  {
    ConnectionHetWD::get_status(d);

    def<double_t>(d, "tau_plus", tau_plus_);
    def<double_t>(d, "tau_minus", tau_minus_);
    def<double_t>(d, "A_plus", A_plus_);
    def<double_t>(d, "A_minus", A_minus_);

    def<double_t>(d, "tau_rec", tau_rec_);
    def<double_t>(d, "U_SE", U_SE_);
    def<double_t>(d, "I", I_);
    def<long_t>(d, "mode", mode_);

  }
  
  void SchemmelS1Connection::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, cm);

    updateValue<double_t>(d, "tau_plus", tau_plus_);
    updateValue<double_t>(d, "tau_minus", tau_minus_);
    updateValue<double_t>(d, "A_plus", A_plus_);
    updateValue<double_t>(d, "A_minus", A_minus_);

    updateValue<double_t>(d, "tau_rec", tau_rec_);
    updateValue<double_t>(d, "U_SE", U_SE_);
    updateValue<double_t>(d, "I", I_);
    updateValue<long_t>(d, "mode", mode_);


    if (I_ > 1.0 || I_ < 0.0) {
       cm.network().message(SLIInterpreter::M_ERROR, "MarkramConnection::set_status()", "require 0<=I<=1.0");
    }

    if (U_SE_ > 1.0 || U_SE_ < 0.0) {
       cm.network().message(SLIInterpreter::M_ERROR, "MarkramConnection::set_status()", "require 0<=U_SE<=1.0");
    }


  }

   /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void SchemmelS1Connection::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);

    set_property<double_t>(d, "tau_pluss", p, tau_plus_);
    set_property<double_t>(d, "tau_minuss", p, tau_minus_);
    set_property<double_t>(d, "A_pluss", p, A_plus_);
    set_property<double_t>(d, "A_minuss", p, A_minus_);

    set_property<double_t>(d, "tau_recs", p, tau_rec_);
    set_property<double_t>(d, "U_SEs", p, U_SE_);
    set_property<double_t>(d, "Is", p, I_);
    set_property<long_t>(d, "modes", p, mode_);

  }


  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void SchemmelS1Connection::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);

    initialize_property_array(d, "tau_pluss");
    initialize_property_array(d, "tau_minuss");
    initialize_property_array(d, "A_pluss");
    initialize_property_array(d, "A_minuss");

    initialize_property_array(d, "tau_recs");
    initialize_property_array(d, "U_SEs");
    initialize_property_array(d, "Is");
    initialize_property_array(d, "modes");

  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void SchemmelS1Connection::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);

    append_property<double_t>(d, "tau_pluss", tau_plus_);
    append_property<double_t>(d, "tau_minuss", tau_minus_);
    append_property<double_t>(d, "A_pluss", A_plus_);
    append_property<double_t>(d, "A_minuss", A_minus_);

    append_property<double_t>(d, "tau_recs", tau_rec_);
    append_property<double_t>(d, "U_SEs", U_SE_);
    append_property<double_t>(d, "Is", I_);
    append_property<long_t>(d, "modes", mode_);

  }

} // of namespace nest
