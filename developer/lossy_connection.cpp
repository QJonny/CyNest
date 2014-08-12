/*
 *  lossy_connection.cpp
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2011 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

//#include "network.h"
#include "dictdatum.h"
#include "connector_model.h"
#include "lossy_connection.h"
#include "event.h"

namespace nest
{

  LossyConnection::LossyConnection() :
    ConnectionHetWD(),
    p_transmit_(1.0)
  {}

  LossyConnection::LossyConnection(const LossyConnection & rhs) :
    ConnectionHetWD(rhs)
  {
    p_transmit_ = rhs.p_transmit_;
  }

  void LossyConnection::get_status(DictionaryDatum & d) const
  {
    ConnectionHetWD::get_status(d);
    def<double_t>(d, "p_transmit", p_transmit_);
  }

  void LossyConnection::set_status(const DictionaryDatum & d, ConnectorModel & cm)
  {
    ConnectionHetWD::set_status(d, cm);
    updateValue<double_t>(d, "p_transmit", p_transmit_);

    if ( p_transmit_ < 0 || p_transmit_ > 1 )
      throw BadProperty("Spike transmission probability must be in [0, 1].");
  }

   /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */
  void LossyConnection::set_status(const DictionaryDatum & d, index p, ConnectorModel & cm)
  {
    ConnectionHetWD::set_status(d, p, cm);
    set_property<double_t>(d, "p_transmits", p, p_transmit_);

    if ( p_transmit_ < 0 || p_transmit_ > 1 )
      throw BadProperty("Spike transmission probability must be in [0, 1].");
  }

  void LossyConnection::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);
    initialize_property_array(d, "p_transmits");
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void LossyConnection::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);
    append_property<double_t>(d, "p_transmits", p_transmit_);

    if ( p_transmit_ < 0 || p_transmit_ > 1 )
      throw BadProperty("Spike transmission probability must be in [0, 1].");
  }

} // of namespace nest
