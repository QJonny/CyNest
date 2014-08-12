
/*
 *  stdp_connection_hom.cpp
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
#include "stdp_connection_rsnn_spikepairing_hom.h"
#include "event.h"

namespace nest
{
  //
  // Implementation of class STDPRSNNSpikePairingHomCommonProperties.
  //

  STDPRSNNSpikePairingHomCommonProperties::STDPRSNNSpikePairingHomCommonProperties() :
    CommonSynapseProperties(),
    tau_plus_(20.0),
    tau_minus_(20.0),
    lambda_(0.01),
    alpha_(1.0),
    mu_plus_(1.0),
    mu_minus_(1.0),
    Wmax_(100.0)
  { }

  void STDPRSNNSpikePairingHomCommonProperties::get_status(DictionaryDatum & d) const
  {
    CommonSynapseProperties::get_status(d);

    def<double_t>(d, "tau_plus", tau_plus_);
    def<double_t>(d, "tau_minus", tau_minus_);
    def<double_t>(d, "lambda", lambda_);
    def<double_t>(d, "alpha", alpha_);
    def<double_t>(d, "mu_plus", mu_plus_);
    def<double_t>(d, "mu_minus", mu_minus_);
    def<double_t>(d, "Wmax", Wmax_);
  }

  void STDPRSNNSpikePairingHomCommonProperties::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    CommonSynapseProperties::set_status(d, cm);

    updateValue<double_t>(d, "tau_plus", tau_plus_);
    updateValue<double_t>(d, "tau_minus", tau_minus_);
    updateValue<double_t>(d, "lambda", lambda_);
    updateValue<double_t>(d, "alpha", alpha_);
    updateValue<double_t>(d, "mu_plus", mu_plus_);
    updateValue<double_t>(d, "mu_minus", mu_minus_);
    updateValue<double_t>(d, "Wmax", Wmax_);
  }


  //
  // Implementation of class STDPRSNNSpikePairingConnectionHom.
  //

  STDPRSNNSpikePairingConnectionHom::STDPRSNNSpikePairingConnectionHom(){ }

  STDPRSNNSpikePairingConnectionHom::STDPRSNNSpikePairingConnectionHom(const STDPRSNNSpikePairingConnectionHom &rhs) :
    ConnectionHetWD(rhs)
  { }

  void STDPRSNNSpikePairingConnectionHom::get_status(DictionaryDatum & d) const
  {
    // base class properties, different for individual synapse
    ConnectionHetWD::get_status(d);
  }

  void STDPRSNNSpikePairingConnectionHom::set_status(const DictionaryDatum & d, ConnectorModel &cm)
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
// 	cm.network().error("STDPRSNNSpikePairingConnectionHom::set_status()", "you are trying to set common properties via an individual synapse.");
//       }
  }

   /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */
  void STDPRSNNSpikePairingConnectionHom::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);

     if (d->known("tau_pluss") ||
         d->known("lambds") ||
         d->known("alphas") ||
         d->known("mu_pluss") ||
         d->known("mu_minuss") ||
         d->known("Wmaxs") )
     {
       cm.network().message(SLIInterpreter::M_ERROR, "STDPRSNNSpikePairingConnectionHom::set_status()", "you are trying to set common properties via an individual synapse.");
     }
  }

  void STDPRSNNSpikePairingConnectionHom::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void STDPRSNNSpikePairingConnectionHom::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);
  }

} // of namespace nest
