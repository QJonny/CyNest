/*
 *  theta_neuron_ps.h
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2009 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#ifndef THETA_NEURON_PS_H
#define THETA_NEURON_PS_H

#include "config.h"

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "../precise/slice_ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"
#include "recordables_map.h"

/* BeginDocumentation

Name: theta_neuron_ps - Ermentrout-Kopell canonical model handling precise spike times

Description:
  Propagation of theta and computation of exact spike times according
  to [1] (scaling parameter alpha omitted) but adapted to the
  time-driven canonical scheme for precise spike timing [2]. Note that
  the total external input I_0 (including I_e) must not be 0. A
  current can be injected, e.g. if I_e is set to a value different
  from 0 and a dc_generator with an amplitude different from -I_e is
  connected to the neuron.
  Theta can be recorded by the multimeter.

Parameters: 
  The following parameters can be set in the status dictionary.

  I_e      double - Constant external input (dimensionless)
  theta    double - Neuron phase (dimensionless)

References:
  [1] McKennoch S, Voegtlin T, and Bushnell L (2008) Spike-Timing
      Error Backpropagation in Theta Neuron Networks. Neural
      Comput. Aug 4. [Epub ahead of print]
  [2] Morrison A, Straube S, Plesser HE, and Diesmann M (2006) Exact
      Subthreshold Integration with Continuous Spike Times in Discrete
      Time Neural Network Simulations. Neural Comput 19(1):47-79.

Sends: SpikeEvent

Receives: SpikeEvent
          CurrentEvent

Author: Kunkel
*/

namespace nest
{
  class theta_neuron_ps : public Archiving_Node
  {
    
    // Boilerplate function declarations -----------------------------

  public:
    
    theta_neuron_ps();
    theta_neuron_ps(const theta_neuron_ps&);

    /*
     * Import all overloaded virtual functions that we
     * override in this class. For background information, 
     * see http://www.gotw.ca/gotw/005.htm.
     */

    using Node::connect_sender;
    using Node::handle;

    port check_connection(Connection &, port);

    port connect_sender(SpikeEvent &, port);
    port connect_sender(CurrentEvent &, port);
    port connect_sender(DataLoggingRequest &, port);
    
    void handle(SpikeEvent &);
    void handle(CurrentEvent &);
    void handle(DataLoggingRequest &);
        
    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);

    bool is_off_grid() const
    {
      return true;
    }

  private:
    void init_state_(const Node &proto);
    void init_buffers_();
    void calibrate();
    void update(Time const &, const long_t, const long_t);
    void propagate_(Time const &, long_t const, double_t const, double_t const, double_t const);

    // END Boilerplate function declarations -------------------------

    // Friends -------------------------------------------------------

    // The next two classes need to be friends to access the State_ class/member
    friend class RecordablesMap<theta_neuron_ps>;
    friend class UniversalDataLogger<theta_neuron_ps>;

  private:

    // Parameters class ----------------------------------------------

    //! Model parameters
    struct Parameters_ {
      double_t I_e; //!< Constant external input
  
      Parameters_(); //!< Set default parameter values

      void get(DictionaryDatum &) const; //!< Store current values in dictionary
      void set(const DictionaryDatum &); //!< Set values from dicitonary
    };
    
    // State variables class -----------------------------------------

    //! State variables of the model
    struct State_ {
      double_t theta; //!< Neuron phase

      State_(const Parameters_ &); //!< Default initialization
      State_(const State_ &);

      void get(DictionaryDatum&) const; //!< Store current values in dictionary

      /**
       * Set state from values in dictionary.
       * Requires Parameters_ as argument to, e.g. check bounds.
       */
      void set(const DictionaryDatum &, const Parameters_ &);
    };    

    // Buffers class -------------------------------------------------

    /**
     * Buffers of the model.
     * Buffers are on par with state variables in terms of persistence,
     * i.e., initalized only upon first Simulate call after ResetKernel
     * or ResetNetwork, but are implementation details hidden from the user.
     */
    struct Buffers_ {
      Buffers_(theta_neuron_ps &); //!< Sets buffer pointers to 0
      Buffers_(const Buffers_ &, theta_neuron_ps &); //!< Sets buffer pointers to 0

      //! Logger for all analog data
      UniversalDataLogger<theta_neuron_ps> logger_;

      /**
       * Buffers incoming spikes (precise spike times)
       * @note Handles also pseudo-events marking return from refractoriness.
       */
      SliceRingBuffer events_;

      //! Buffers incoming currents
      RingBuffer currents_;
    };

    // Variables class ----------------------------------------------

    /**
     * Internal variables of the model.
     * Variables are re-initialized upon each call to Simulate.
     */
    struct Variables_ {
      double_t h_ms; //!< Simulation resolution in [ms]
      double_t I_0;  //!< Total external input including I_e
    };

    // Access functions for UniversalDataLogger ----------------------

    //! Read out neuron phase, used by UniversalDataLogger
    double_t get_theta_() const { return S_.theta; }

    // Data members --------------------------------------------------

    // keep the order of these lines, seems to give best performance
    Parameters_ P_;
    State_      S_;
    Variables_  V_;
    Buffers_    B_;

    //! Mapping of recordables names to access functions
    static RecordablesMap<theta_neuron_ps> recordablesMap_;
  };


  // Boilerplate inline function definitions -------------------------

  inline
  port theta_neuron_ps::check_connection(Connection &c, port receptor_type)
  {
    SpikeEvent e;
    e.set_sender(*this);
    c.check_event(e);
    return c.get_target()->connect_sender(e, receptor_type);
  }

  inline
  port theta_neuron_ps::connect_sender(SpikeEvent &, port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
 
  inline
  port theta_neuron_ps::connect_sender(CurrentEvent &, port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
 
  inline
  port theta_neuron_ps::connect_sender(DataLoggingRequest& dlr,
				       port receptor_type) 
  { 
    if (receptor_type != 0) 
      throw UnknownReceptorType(receptor_type, get_name()); 
    return B_.logger_.connect_logging_device(dlr, recordablesMap_); 
  } 

  inline
  void theta_neuron_ps::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    S_.get(d);
    Archiving_Node::get_status(d);

    (*d)[names::recordables] = recordablesMap_.get_list();
  }

  inline
  void theta_neuron_ps::set_status(const DictionaryDatum &d)
  {
    Parameters_ ptmp = P_; // temporary copy in case of errors
    ptmp.set(d);           // throws if BadProperty
    State_      stmp = S_; // temporary copy in case of errors
    stmp.set(d, ptmp);     // throws if BadProperty

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

#endif //THETA_NEURON_PS_H
