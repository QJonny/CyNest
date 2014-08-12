/*
 *  inh_gamma_generator_ng.h
 *
 *  Based loosely on the poisson_generator
 *
 *  AUTHOR: Eilif Muller
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2007 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#ifndef INH_GAMMA_GENERATOR_NG_H
#define INH_GAMMA_GENERATOR_NG_H

#include "config.h"

#ifdef HAVE_GSL

/****************************************/
/* class inh_gamma_generator_ng         */
/*                  Vers. 1.0       EBM */
/*                  Implementation: EBM */
/****************************************/

// Based on:

/****************************************/
/* class poisson_generator              */
/*                  Vers. 1.0       hep */
/*                  Implementation: hep */
/****************************************/

#include "nest.h"
#include "event.h"
#include "node.h"
#include "stimulating_device.h"
#include "exp_randomdev.h"
#include "connection.h"

namespace nest
{
//! class inh_gamma_generator_ng                           
/*! Class inh_gamma_generator_ng simulates an inhomogeneous
gamma renewal process.
*/


/* BeginDocumentation

  Name: inh_gamma_generator_ng - simulates an inhomogeneous gamma renewal process, with indepedent spike trains for each target.
  
  Description: The inh_gamma_generator_ng (next generation) simulates
  an inhomogeneous gamma renewal process, i.e. interspike intervals
  are gamma distributed. Unlike inh_gamma_generator, it provides
  indepedent spike trains for each of its targets.

  It will generate a _unique_ spike train for each of it's targets. If
  you do not want this behavior and need the same spike train for all
  targets, you use the inh_gamma_generator.

  Note:
  inh_gamma_generator_ng requires synapse types with a _S suffix for
  "selective".  For example: static_synapse_S,
  'stdp_pl_synapse_hom_S'.  If non-selective synapses are used, the
  behaviour is undefined.

  inh_gamma_generator_ng supports only targets of uniform synapse type.
  It will throw an exception if it is used incorrectly.  Multiple 
  generators should be used if multiple synapse types are required.
    
Parameters:
   The following parameters appear in the element's status dictionary:

   rmax     - max firing rate, i.e. max(1/b) (double, var)
   b        - time histogram of parameter b of gamma distribution in units of seconds (vector<double>, var)
   a        - time histogram of parameter a of gamma distribution (dimensionless) (vector<double>,var)
   tbins    - time bins of time histogram of a,b in units of milliseconds (vector<double>,var)

Remarks:

  The spike-train generation algorithm is based on thinning algortithm
  for a general hazard function in: Devroye, L "Non-Uniform Random
  Variate Generation", Springer Verlag, New York, 1986

  Gamma distribution parameters a,b are as in Matlab or GSL (GNU
  Scientific Library)

  FIRING RATE rmax should not approach 1/dt ... the NEST simulation
  resolution

  See: E.Muller et. al. "Spike-Frequency Adapting Neural Ensembles:
  Beyond Mean-Adaptation and Renewal Theories", Neural Computation
  19(11).

SeeAlso: inh_gamma_generator, poisson_generator_ps, Device, parrot_neuron
Author: Eilif Muller
References: GNU Scientific Library, http://www.gnu.org/software/gsl
*/


  class inh_gamma_generator_ng : public Node, public StimulatingDevice<SpikeEvent>
  {

  public:
	 
    /** 
     * The generator is threaded, so the RNG to use is determined
     * at run-time, depending on thread.
     */
    inh_gamma_generator_ng();
    inh_gamma_generator_ng(inh_gamma_generator_ng const&);

    bool has_proxies() const {return false;}


    port check_connection(Connection&, port);
    
    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &) ;

  protected:
    
    void init_parameters_(Node const*);
    void init_dynamic_state_(Node const*);
    void calibrate();
    void update(Time const &, const long_t, const long_t);

    void calc_h(double h[], double t, double a, double b, double dt);

    
  private: 
    librandom::ExpRandomDev exprand_dev_;  //!< random deviate generator

    std::vector<double>  a_;   //!< process parameter a dimensionless - time histogram - sorted with tbins
    std::vector<double>  b_;   //!< process parameter in seconds - time histogram - sorted with tbins
    std::vector<double>  tbins_;  //!< bins for a,b time histogram - sorted in ascending order in milliseconds
    double rmax_; // !< max 1/b in experiment... sets rate of thinned poisson process in Hz
    std::vector<double>  tspike_;   //!< next spike of poisson process to be thinned for each target in milliseconds
    std::vector<double>  tlast_;    //!< last spike of gamma process in milliseconds
    size_t bin_position_; //!< stores current time bin position in time histogram of a,b

    size_t num_connections_;

    long synapse_context_; //!< inh_gamma_generator_ng supports sending events via one synapse type.  If multiple synapse types are required, multiple nodes should be used, one per synapse type.

    
  };

  inline
  port inh_gamma_generator_ng::check_connection(Connection& c, port receptor_type)
  {
    // support only one synapse type for targets
    if (synapse_context_!=-1) {
      if (synapse_context_!=net_->get_synapse_context())
        throw IllegalConnection();
    }
    else 
      synapse_context_=net_->get_synapse_context();
    
    SpikeEvent e;
    e.set_sender(*this);
    c.check_event(e);

    // set the port to the number of connections
    e.set_port(tspike_.size());

    // add one more to the tspike_ records list
    tspike_.push_back(0.0);
    //tlast_.push_back(-1000.0);
    tlast_.push_back(0.0);

    return c.get_target()->connect_sender(e, receptor_type);
  }
  
} // namespace nest


#endif //HAVE_GSL
#endif //INH_GAMMA_GENERATOR_NG_H
