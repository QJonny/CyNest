/*
 *  a2eif_cond_exp_HW.h
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

#ifndef A2EIF_COND_EXP_HW_H
#define A2EIF_COND_EXP_HW_H

#include "config.h"

#ifdef HAVE_GSL_1_11

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"
#include "recordables_map.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

/*----------------------------*/
/*!!! IN DEVELOPMENT STATE !!!*/
/*----------------------------*/

/* BeginDocumentation
Name: a2eif_cond_exp_HW - Conductance based exponential integrate-and-fire neuron model according to Brette and Gerstner (2005)* with two adaptation terms** and hardware restrictions/simplifications.

Description:
a2eif_cond_exp_HW is the adaptive exponential integrate and fire neuron with two adaptation terms for short-term
adaptation as well as long-term adaptation. This model includes hardware restrictions of the FACETS system. Some parameters are fixed because of simplification.
Synaptic conductances are modelled as exponential functions.

This implementation uses the embedded 4th order Runge-Kutta-Fehlberg solver with adaptive stepsize to integrate
the differential equation.

The membrane potential is given by the following differential equation:
C dV/dt= -g_L(V-E_L) +g_L*Delta_T*exp((V-V_T)/Delta_T) -g_e(t)(V-E_e) -g_i(t)(V-E_i) -w1 -w2 +I_e

and
tau_w1 * dw1/dt= a(V +V_woff -E_L) -w1
tau_w2 * dw2/dt= a(V +V_woff -E_L) -w2

Parameters: 
The following parameters can be set in the status dictionary.

Dynamic state variables:
  V_m        double - Membrane potential in mV
  g_ex       double - Excitatory synaptic conductance in nS.
  g_in       double - Inhibitory synaptic conductance in nS.
  w1         double - Short-term spike-adaptation current in pA.
  w2         double - Long-term spike-adaptation current in pA.

Membrane Parameters:
  C_m        double - Capacity of the membrane in pF
  t_ref      double - Duration of refractory period in ms. 
  V_peak     double - Spike detection threshold in mV.
  V_reset    double - Reset value for V_m after a spike. In mV.
  E_L        double - Leak reversal potential in mV. 
  g_L        double - Leak conductance in nS.
  I_e        double - Constant external input current in pA.

Spike adaptation parameters:
  a1         double - Short-term subthreshold adaptation in nS.
  b1         double - Short-term spike-triggered adaptation in pA.
  tau_w1     double - Short-term adaptation time constant in ms
  a2         double - Long-term subthreshold adaptation in nS.
  b2         double - Long-term spike-triggered adaptation in pA.
  tau_w2     double - Long-term adaptation time constant in ms
  Delta_T    double - Slope factor in mV
  V_t        double - Spike initiation threshold in mV (V_th can also be used for compatibility).

Synaptic parameters
  E_ex       double - Excitatory reversal potential in mV.
  tau_syn_ex double - Rise time of excitatory synaptic conductance in ms (exponential function).
  E_in       double - Inhibitory reversal potential in mV.
  tau_syn_in double - Rise time of the inhibitory synaptic conductance in ms (exponential function).

Hardware parameters that can be set
  V_woff     double - Voltage offset in the adaptation curcuit in mV
  V_min      double - Desirable minimal voltage in the simulation in mV
  s_hw       double - Hardware voltage scaling factor restricting maximum voltage range

Hardware parameters that cannot be set
  V_rang     double - Hardware voltage range, internal parameter in mV

Author: March 2010 A. Kononov, T. Pfeil

Sends: SpikeEvent

Receives: SpikeEvent, CurrentEvent, DataLoggingRequest

References: *  Brette R and Gerstner W (2005) Adaptive Exponential Integrate-and-Fire Model as 
               an Effective Description of Neuronal Activity. J
               Neurophysiol 94:3637-3642

            ** Mensi S, Naud R, Avermann M, Petersen C and Gerstner W (2010).
               Complexity and performance in simple neuron models.
               Conference Abstract: Computational and Systems Neuroscience 2010.
               doi: 10.3389/conf.fnins.2010.03.00064 

SeeAlso: a2eif_cond_exp, aeif_cond_alpha, aeif_cond_exp
*/

namespace nest
{
  /**
   * Function computing right-hand side of ODE for GSL solver.
   * @note Must be declared here so we can befriend it in class.
   * @note Must have C-linkage for passing to GSL. Internally, it is
   *       a first-class C++ function, but cannot be a member function
   *       because of the C-linkage.
   * @note No point in declaring it inline, since it is called
   *       through a function pointer.
   * @param void* Pointer to model neuron instance.
   */
  extern "C"
  int a2eif_cond_exp_HW_dynamics (double, const double*, double*, void*);

  class a2eif_cond_exp_HW:
  public Archiving_Node
  {
    
  public:        
    
    a2eif_cond_exp_HW();
    a2eif_cond_exp_HW(const a2eif_cond_exp_HW&);
    ~a2eif_cond_exp_HW();

    /**
     * Import sets of overloaded virtual functions.
     * @see Technical Issues / Virtual Functions: Overriding, Overloading, and Hiding
     */

    using Node::connect_sender;
    using Node::handle;

    port check_connection(Connection &, port);

    void handle(SpikeEvent &);
    void handle(CurrentEvent &);
    void handle(DataLoggingRequest &);

    port connect_sender(SpikeEvent &, port);
    port connect_sender(CurrentEvent &, port);
    port connect_sender(DataLoggingRequest &, port);
    
    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);
    
  private:
    
    void init_state_(const Node &proto);
    void init_buffers_();
    void calibrate();
    void update(const Time &, const long_t, const long_t);

    // END Boilerplate function declarations ----------------------------

    // Friends --------------------------------------------------------
    // make dynamics function quasi-member
    friend int a2eif_cond_exp_HW_dynamics(double, const double*, double*, void*);

    // The next two classes need to be friends to access the State_ class/member
    friend class RecordablesMap<a2eif_cond_exp_HW>;
    friend class UniversalDataLogger<a2eif_cond_exp_HW>;

  private:
    // ---------------------------------------------------------------- 

    //! Independent parameters
    struct Parameters_
    {
      double_t V_peak_;     //!< Spike detection threshold in mV
      double_t V_reset_;    //!< Reset Potential in mV
      double_t t_ref_;      //!< Refractory period in ms
      double_t g_L;         //!< Leak Conductance in nS
      double_t C_m;         //!< Membrane Capacitance in pF
      double_t E_ex;        //!< Excitatory reversal Potential in mV
      double_t E_in;        //!< Inhibitory reversal Potential in mV
      double_t E_L;         //!< Leak reversal Potential (aka resting potential) in mV
      double_t Delta_T;     //!< Slope faktor in ms
      double_t tau_w1;      //!< Short-term adaptation time-constant in ms
      double_t a1;          //!< Short-term subthreshold adaptation in nS
      double_t b1;          //!< Short-term spike-triggered adaptation in pA
      double_t tau_w2;      //!< Long-term adaptation time-constant in ms
      double_t a2;          //!< Long-term subthreshold adaptation in nS
      double_t b2;          //!< Long-term spike-triggered adaptation in pA
      double_t V_th;        //!< Spike threshold in mV
      double_t t_ref;       //!< Refractory period in ms
      double_t tau_syn_ex;  //!< Excitatory synaptic rise time in ms
      double_t tau_syn_in;  //!< Inhibitory synaptic rise time in ms
      double_t I_e;         //!< Intrinsic current in pA
      double_t V_woff;      //!< Offset in adaptation term in mV
      double_t V_min;       //!< Desirable overall minimal voltage in mV
      double_t s_hw;        //!< Hardware voltage scaling factor
      double_t V_rang;      //!< Voltage range in hardware in mV

      Parameters_();  //!< Sets default parameter values

      void get(DictionaryDatum &) const;  //!< Store current values in dictionary
      void set(const DictionaryDatum &);  //!< Set values from dicitonary
    };

  public:
    // ---------------------------------------------------------------- 

    /**
     * State variables of the model.
     * @note Copy constructor and assignment operator required because
     *       of C-style array.
     */
    struct State_
    {
      /**
       * Enumeration identifying elements in state array State_::y_.
       * The state vector must be passed to GSL as a C array. This enum
       * identifies the elements of the vector. It must be public to be
       * accessible from the iteration function.
       */  
      enum StateVecElems
      {
  V_M   = 0,
  G_EXC    ,  // 1
  G_INH    ,  // 2
  W1       ,  // 3
  W2       ,  // 4
	STATE_VEC_SIZE
      };

      double_t y_[STATE_VEC_SIZE];  //!< neuron state, must be C-array for GSL solver
      int_t    r_;           //!< number of refractory steps remaining

      State_(const Parameters_ &);  //!< Default initialization
      State_(const State_ &);
      State_& operator = (const State_ &);

      void get(DictionaryDatum &) const;
      void set(const DictionaryDatum &, const Parameters_ &);
    };

    // ---------------------------------------------------------------- 

    /**
     * Buffers of the model.
     */
    struct Buffers_
    {
      Buffers_(a2eif_cond_exp_HW &);
      Buffers_(const Buffers_ &, a2eif_cond_exp_HW &);

      UniversalDataLogger<a2eif_cond_exp_HW> logger_;    //! Logger for all analog data

      /** buffers and sums up incoming spikes/currents */
      RingBuffer spike_exc_;
      RingBuffer spike_inh_;
      RingBuffer currents_;

      /** GSL ODE stuff */
      gsl_odeiv_step*    s_;    //!< stepping function
      gsl_odeiv_control* c_;    //!< adaptive stepsize control function
      gsl_odeiv_evolve*  e_;    //!< evolution function
      gsl_odeiv_system   sys_;  //!< struct describing system

      // IntergrationStep_ should be reset with the neuron on ResetNetwork,
      // but remain unchanged during calibration. Since it is initialized with
      // step_, and the resolution cannot change after nodes have been created,
      // it is safe to place both here.
      double_t step_;           //!< step size in ms
      double   IntegrationStep_;//!< current integration time step, updated by GSL

      /** 
       * Input current injected by CurrentEvent.
       * This variable is used to transport the current applied into the
       * _dynamics function computing the derivative of the state vector.
       * It must be a part of Buffers_, since it is initialized once before
       * the first simulation, but not modified before later Simulate calls.
       */
      double_t I_stim_;
    };

    // ---------------------------------------------------------------- 

    /**
     * Internal variables of the model.
     */
    struct Variables_
    {
      int_t RefractoryCounts_;
    };

    // Access functions for UniversalDataLogger -------------------------------

    //! Read out state vector elements, used by UniversalDataLogger
    template <State_::StateVecElems elem>
    double_t get_y_elem_() const { return S_.y_[elem]; }

    // ---------------------------------------------------------------- 

    Parameters_ P_;
    State_      S_;
    Variables_  V_;
    Buffers_    B_;

    //! Mapping of recordables names to access functions
    static RecordablesMap<a2eif_cond_exp_HW> recordablesMap_;
  };


  inline  
  port a2eif_cond_exp_HW::check_connection(Connection &c, port receptor_type)
  {
    SpikeEvent e;
    e.set_sender(*this);
    c.check_event(e);
    return c.get_target()->connect_sender(e, receptor_type);
  }

  inline
  port a2eif_cond_exp_HW::connect_sender(SpikeEvent &, port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
 
  inline
  port a2eif_cond_exp_HW::connect_sender(CurrentEvent &, port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
 
  inline
  port a2eif_cond_exp_HW::connect_sender(DataLoggingRequest &dlr, 
				      port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return B_.logger_.connect_logging_device(dlr, recordablesMap_);
  }

  inline
  void a2eif_cond_exp_HW::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    S_.get(d);
    Archiving_Node::get_status(d);

    (*d)[names::recordables] = recordablesMap_.get_list();
  }

  inline
  void a2eif_cond_exp_HW::set_status(const DictionaryDatum &d)
  {
    Parameters_ ptmp = P_;  // temporary copy in case of errors
    ptmp.set(d);            // throws if BadProperty
    State_      stmp = S_;  // temporary copy in case of errors
    stmp.set(d, ptmp);      // throws if BadProperty

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

#endif //HAVE_GSL_1_11
#endif //A2EIF_COND_EXP_HW_H
