/*
 *  iaf_cond_delta.h
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2005 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#ifndef IAF_COND_DELTA_H
#define IAF_COND_DELTA_H

#include "config.h"

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "connection.h"

#include "analog_data_logger.h"

/* BeginDocumentation

   Name: iaf_cond_delta - Leaky integrate-and-fire neuron model with delta post-synaptic conductances

   Description:

   iaf_cond_delta is an implementation of a spiking neuron using IAF dynamics with
   conductance-based synapses. Incoming spike events induce a post-synaptic change 
   of conductance modelled by a delta function, and hence the potential jumps discretely
   on each spike arrival by an amount that depends on the membrane potential and the
   reversal potential.

   The delta function is normalised such that the weight w of an event is time-integrated 
   conductance divided by capacitance, a dimensionless quantity (\int g/C delta(t-t_ev)dt),
   multiplied by -1 in the case of an inhibitory input. (The sign of the weight is purely
   to distinguish which type of event it is: w<0 iff event is inhibitory)

   The conductance event induces a membrane potential displacement of 
   dV = w*(E_ex - V) for excitatory input
   dV = -w*(E_in - V) for inhibitory input

   The threshold crossing is followed by an absolute refractory period
   during which the membrane potential is clamped to the resting potential.
   Spikes arriving while the neuron is refractory are discarded.

   Parameters: 
   The following parameters can be set in the status dictionary.

   V_m        double - Membrane potential in mV 
   E_L        double - Leak reversal potential in mV.
   C_m        double - Capacity of the membrane in pF
   t_ref      double - Duration of refractory period in ms. 
   V_th       double - Spike threshold in mV.
   V_reset    double - Reset potential of the membrane in mV.
   E_ex       double - Excitatory reversal potential in mV.
   E_in       double - Inhibitory reversal potential in mV.
   g_L        double - Leak conductance in nS;

   Sends: SpikeEvent

   Receives: SpikeEvent, PotentialRequest

   Author: Chris Trengove

   SeeAlso: iaf_psc_delta, iaf_psc_exp, iaf_cond_exp
*/

namespace nest
{
  class Network;
  
  class iaf_cond_delta : public Archiving_Node
  {
    
  public:        
    
    typedef Node base;
    
    iaf_cond_delta();
    iaf_cond_delta(const iaf_cond_delta&);

    /**
     * Import sets of overloaded virtual functions.
     * We need to explicitly include sets of overloaded
     * virtual functions into the current scope.
     * According to the SUN C++ FAQ, this is the correct
     * way of doing things, although all other compilers
     * happily live without.
     */
    using Node::connect_sender;
    using Node::handle;

    port check_connection(Connection&, port);
    
    void handle(SpikeEvent &);
    void handle(PotentialRequest &);
    
    port connect_sender(SpikeEvent &, port);
    port connect_sender(PotentialRequest &, port);
    
    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);
    
  private:
    void init_state_(const Node& proto);
    void init_buffers_();
    void calibrate();
    void update(Time const &, const long_t, const long_t);

  public:

    // ---------------------------------------------------------------- 

    /** 
     * Independent parameters of the model. 
     * These parameters must be passed to the iteration function that
     * is passed to the GSL ODE solvers. Since the iteration function
     * is a C++ function with C linkage, the parameters can be stored
     * in a C++ struct with member functions, as long as we just pass
     * it by void* from C++ to C++ function. The struct must be public,
     * though, since the iteration function is a function with C-linkage,
     * whence it cannot be a member function of iaf_cond_delta.
     * @note One could achieve proper encapsulation by an extra level
     *       of indirection: Define the iteration function as a member
     *       function, plus an additional wrapper function with C linkage.
     *       Then pass a struct containing a pointer to the node and a
     *       pointer-to-member-function to the iteration function as void*
     *       to the wrapper function. The wrapper function can then invoke
     *       the iteration function on the node (Stroustrup, p 418). But
     *       this appears to involved, and the extra indirections cost.
     */
    struct Parameters_ {
      double_t E_L_;        //!< Leak reversal Potential (aka resting potential) in mV
      double_t V_th_;       //!< Threshold Potential in mV
      double_t V_reset_;    //!< Reset Potential in mV
      double_t E_ex_;       //!< Excitatory reversal Potential in mV
      double_t E_in_;       //!< Inhibitory reversal Potential in mV
      double_t t_ref_;      //!< Refractory period in ms
      double_t g_L_;        //!< Leak Conductance in nS; note, tau_m=C_m/g_L
      double_t C_m_;        //!< Membrane Capacitance in pF
  
      Parameters_();  //!< Sets default parameter values

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary
      void set(const DictionaryDatum&);  //!< Set values from dicitonary
    };
    
  private:
    // ---------------------------------------------------------------- 

    /**
     * State variables of the model.
     */
    struct State_ {
      double_t V_m_;  //!< membrane potential
      int_t    r_;     //!< number of refractory steps remaining

      State_();  //!< Default initialization

      void get(DictionaryDatum&, const Parameters_&) const;
      void set(const DictionaryDatum&, const Parameters_&);
    };    

    // ---------------------------------------------------------------- 

    /**
     * Buffers of the model.
     */
    struct Buffers_ {
      /** buffers and sums up incoming exc and inh spikes */
      RingBuffer spikes_ex_;
      RingBuffer spikes_in_;

      /**
       * Buffer for membrane potential.
       */
      AnalogDataLogger<PotentialRequest>           potentials_;
    };

    // ---------------------------------------------------------------- 

    /**
     * Internal variables of the model.
     */
    struct Variables_ { 

      double_t P_L_; // propagator of the passive leak
      int_t    RefractoryCounts_;
    };

    // ---------------------------------------------------------------- 

    Parameters_ P_;
    State_      S_;
    Variables_  V_;
    Buffers_    B_;
  };

  
  inline
  port iaf_cond_delta::check_connection(Connection& c, port receptor_type)
  {
    SpikeEvent e;
    e.set_sender(*this);
    c.check_event(e);
    return c.get_target()->connect_sender(e, receptor_type);
  }

  inline
  port iaf_cond_delta::connect_sender(SpikeEvent&, port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
 
  inline
  port iaf_cond_delta::connect_sender(PotentialRequest& pr, port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    B_.potentials_.connect_logging_device(pr);
    return 0;
  }

  inline
  void iaf_cond_delta::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    S_.get(d, P_);
    Archiving_Node::get_status(d);
  }

  inline
  void iaf_cond_delta::set_status(const DictionaryDatum &d)
  {
    Parameters_ ptmp = P_;  // temporary copy in case of errors
    ptmp.set(d);                       // throws if BadProperty
    State_      stmp = S_;  // temporary copy in case of errors
    stmp.set(d, ptmp);                 // throws if BadProperty

    // We now know that (ptmp, stmp) are consistent. We do not 
    // write them back to (P_, S_) before we are also sure that 
    // the properties to be set in the parent class are internally 
    // consistent.
    Archiving_Node::set_status(d);

    // if we get here, temporaries contain consistent set of properties
    P_ = ptmp;
    S_ = stmp;
  }
  
} // namespace

#endif //IAF_COND_DELTA_H
