/*
 *  inh_gamma_generator.h
 *
 *  Based loosely on the spike_generator
 *
 *  AUTHOR: Eilif Muller
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2010 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#ifndef INH_GAMMA_GENERATOR_H
#define INH_GAMMA_GENERATOR_H

#include "config.h"

#ifdef HAVE_GSL

/****************************************/
/* class inh_gamma_generator            */
/*                  Vers. 2.0       EBM */
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
//! class inh_gamma_generator
/*! Class inh_gamma_generator simulates an inhomogeneous
gamma renewal process.
*/


/* BeginDocumentation

  Name: inh_gamma_generator - simulates an inhomogeneous gamma renewal process, with indepedent spike trains for each target.
  
  Description: The inh_gamma_generator (next generation) simulates
  an inhomogeneous gamma renewal process, i.e. interspike intervals
  are gamma distributed. Unlike inh_gamma_generator, it provides
  indepedent spike trains for each of its targets.

  It will generate a _unique_ spike train for each of it's targets. If
  you do not want this behavior and need the same spike train for all
  targets, you use the inh_gamma_generator.

  Note:
  inh_gamma_generator requires synapse types with a _S suffix for
  "selective".  For example: static_synapse_S,
  'stdp_pl_synapse_hom_S'.  If non-selective synapses are used, the
  behaviour is undefined.

  inh_gamma_generator supports only targets of uniform synapse type.
  Presently, it cannot throw an exception if it is used incorrectly in
  this regard.  Multiple generators should be used if multiple synapse
  types are required.
    
Parameters:
   The following parameters appear in the element's status dictionary:

   b        - time histogram of parameter b of gamma distribution in units of seconds (vector<double>, var)
   a        - time histogram of parameter a of gamma distribution (dimensionless) (vector<double>,var)
   tbins    - time bins of time histogram of a,b in units of milliseconds (vector<double>,var)

Remarks:

  The spike-train generation algorithm is based on thinning algortithm
  for a general hazard function in: Devroye, L "Non-Uniform Random
  Variate Generation", Springer Verlag, New York, 1986

  Gamma distribution parameters a,b are as in Matlab or GSL (GNU
  Scientific Library)

  FIRING RATE rmax of the surrogate Poisson process (rmax = max(1/b))
  should not approach 1/dt, the NEST simulation resolution.

  See: E.Muller et. al. "Spike-Frequency Adapting Neural Ensembles:
  Beyond Mean-Adaptation and Renewal Theories", Neural Computation
  19(11).

SeeAlso: poisson_generator_ps, Device, parrot_neuron
Author: Eilif Muller
References: GNU Scientific Library, http://www.gnu.org/software/gsl
*/


  class inh_gamma_generator : public Node
  {

  public:
	 
    /** 
     * The generator is threaded, so the RNG to use is determined
     * at run-time, depending on thread.
     */
    inh_gamma_generator();
    inh_gamma_generator(const inh_gamma_generator&);

    bool has_proxies() const {return false;}


    port check_connection(Connection&, port);
    
    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &) ;

  private:
    
    void init_state_(const Node&);
    void init_buffers_();
    void calibrate();
    void update(Time const &, const long_t, const long_t);

    void calc_h(double h[], double t, double a, double b, double dt);

    // ------------------------------------------------------------

    struct State_ {

      std::vector<double_t>  tspike_;   //!< next spike of poisson process to be thinned for each target in milliseconds
      std::vector<double_t>  tlast_;    //!< last spike of gamma process in milliseconds
      size_t bin_position_; //!< stores current time bin position in time histogram of a,b
      //long synapse_context_; //!< inh_gamma_generator supports sending events via one synapse type.  If multiple synapse types are required, multiple nodes should be used, one per synapse type.

      State_();  //!< Sets default state value
      State_(const State_&);  //!< Sets default state value


      void get(DictionaryDatum&) const;  //!< Store current values in dictionary

    };

    // ------------------------------------------------------------
    
    struct Parameters_ {

      std::vector<double_t>  a_;   //!< process parameter a dimensionless - time histogram - sorted with tbins
      std::vector<double_t>  b_;   //!< process parameter in seconds - time histogram - sorted with tbins
      std::vector<double_t>  tbins_;  //!< bins for a,b time histogram - sorted in ascending order in milliseconds

      Parameters_();  //!< Sets default parameter values
      Parameters_(const Parameters_&);  //!< Recalibrate all times

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary
      
      /**
       * Set values from dicitonary.
       * @note State is passed so that the position can be reset if the 
       *       spike_times_ or spike_weights_ vector has been filled with
       *       new data, or if the origin was reset.
       */
      void set(const DictionaryDatum&, State_&);  
    };
        
    // ------------------------------------------------------------

    StimulatingDevice<SpikeEvent> device_;
    Parameters_ P_;
    State_      S_;
    librandom::ExpRandomDev exprand_dev_;  //!< random deviate generator
    size_t num_connections_;
    double_t rmax_; // !< max 1/b in experiment... sets rate of thinned poisson process in Hz

    
  };

  inline
  port inh_gamma_generator::check_connection(Connection& c, port receptor_type)
  {
    // support only one synapse type for targets
    //if (synapse_context_!=-1) {
    //  if (synapse_context_!=net_->get_synapse_context())
    //    throw IllegalConnection();
    //}
    //else 
    //  synapse_context_=net_->get_synapse_context();
    
    SpikeEvent e;
    e.set_sender(*this);
    c.check_event(e);

    // set the port to the number of connections
    e.set_port(S_.tspike_.size());

    // add one more to the tspike_ records list
    S_.tspike_.push_back(0.0);
    //tlast_.push_back(-1000.0);
    S_.tlast_.push_back(0.0);

    return c.get_target()->connect_sender(e, receptor_type);
  }

inline
void inh_gamma_generator::get_status(DictionaryDatum &d) const
{
  P_.get(d);

  // for debugging
  S_.get(d);
  (*d)["rmax"] = new DoubleDatum(rmax_);


  device_.get_status(d);
}

inline
void inh_gamma_generator::set_status(const DictionaryDatum &d)
{
  Parameters_ ptmp = P_;  // temporary copy in case of errors
  ptmp.set(d, S_);        // throws if BadProperty

  // We now know that ptmp is consistent. We do not write it back
  // to P_ before we are also sure that the properties to be set
  // in the parent class are internally consistent.
  device_.set_status(d);

  // if we get here, temporary contains consistent set of properties
  P_ = ptmp;

  // setup rmax 
  calibrate();
}


  
} // namespace nest


#endif //HAVE_GSL
#endif //INH_GAMMA_GENERATOR_H
