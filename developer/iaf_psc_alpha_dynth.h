/*
 *  iaf_psc_alpha_dynth.h
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2004-2011 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#ifndef IAF_PSC_ALPHA_DYNTH_H
#define IAF_PSC_ALPHA_DYNTH_H

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"
#include "recordables_map.h"

/* BeginDocumentation
Name: iaf_psc_alpha_dynth - Leaky integrate-and-fire with dynamic threshold neuron model

Description:

  iaf_psc_alpha_dynth is a copy of iaf_psc_alpha_dynth neuron model with addition
  of a dynamic threshold (see figure 4.5 in [1]). One of the earliest extensions 
  to the integrate-and-fire model, dynamic threshold introduces simple spike 
  frequency adaptation by increasing the firing threshold by some delta and then
  decaying it back exponentially to the initial threshold value with a specified
  time constant. Specifically, after each emitted spike the threshold is increased 
  (cumulatively) by ’V_th_delta’ and decayed with time constant ’tau_th’.

Parameters: 

  The following parameters can be set in the status dictionary.

  V_m        double - Membrane potential in mV 
  E_L        double - Resting membrane potential in mV. 
  C_m        double - Capacity of the membrane in pF
  tau_m      double - Membrane time constant in ms.
  t_ref      double - Duration of refractory period in ms. 
  V_th       double - Spike threshold in mV.
  V_th_delta double - Post-spike dynamic threshold change in mV. // MOD
  V_reset    double - Reset potential of the membrane in mV.
  tau_th     double - Dynamic threshold time constant in ms. // MOD
  tau_syn_ex double - Rise time of the excitatory synaptic alpha function in ms.
  tau_syn_in double - Rise time of the inhibitory synaptic alpha function in ms.
  I_e        double - Constant external input current in pA.
  V_min      double - Absolute lower value for the membrane potential.
 
References:
  [1] Gerstner W & Kistler W (2002) Spiking Neuron Models. Single Neurons, 
      Populations, Plasticity. Cambridge University Press

Sends: SpikeEvent

Receives: SpikeEvent, CurrentEvent, DataLoggingRequest

Author:  June 2011, Kazeka
SeeAlso: iaf_psc_alpha_dynth
*/

namespace nest
{
  class Network;

  /**
   * Leaky integrate-and-fire neuron with alpha-shaped PSCs and dynamic threshold.
   */
  class iaf_psc_alpha_dynth : public Archiving_Node
  {
    
  public:        
    
    typedef Node base;
    
    iaf_psc_alpha_dynth();
    iaf_psc_alpha_dynth(const iaf_psc_alpha_dynth&);

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
    void handle(CurrentEvent &);
    void handle(DataLoggingRequest &); 
    
    port connect_sender(SpikeEvent&, port);
    port connect_sender(CurrentEvent&, port);
    port connect_sender(DataLoggingRequest &, port);

    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);

  private:

    void init_state_(const Node& proto);
    void init_buffers_();
    void calibrate();

    void update(Time const &, const long_t, const long_t);

    // The next two classes need to be friends to access the State_ class/member
    friend class RecordablesMap<iaf_psc_alpha_dynth>;
    friend class UniversalDataLogger<iaf_psc_alpha_dynth>;

    // ---------------------------------------------------------------- 

    /** 
     * Independent parameters of the model. 
     */
    struct Parameters_ {
  
      /** Membrane time constant in ms. */
      double_t Tau_; 

      /** Membrane capacitance in pF. */
      double_t C_;
    
      /** Refractory period in ms. */
      double_t TauR_;

      /** Resting potential in mV. */
      double_t U0_;

      /** External current in pA */
      double_t I_e_;

      /** reset value of the membrane potential */
      double_t V_reset_;

      /** Threshold, RELATIVE TO RESTING POTENTIAL(!).
          I.e. the real threshold is (U0_+Theta_). */
      double_t Theta_;

      /** Post-spike dynamic threshold change, RELATIVE TO THRESHOLD (!).
          I.e. the immeate post-spike threshold is (ThetaDelta_+Theta_+U0_). */
      double_t ThetaDelta_; // MOD

      /** Lower bound, RELATIVE TO RESTING POTENTIAL(!).
          I.e. the real lower bound is (U0_+LowerBound_). */
      double_t LowerBound_;

      /** Time constant of excitatory synaptic current in ms. */
      double_t tau_ex_;

      /** Time constant of inhibitory synaptic current in ms. */
      double_t tau_in_;

      /** Time constant dynamic threshold change in ms. */
      double_t tau_th_; // MOD

      Parameters_();  //!< Sets default parameter values

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary

      /** Set values from dictionary.
       * @returns Change in reversal potential E_L, to be passed to State_::set()
       */
      double set(const DictionaryDatum&);
    };
    
    // ---------------------------------------------------------------- 

    /**
     * State variables of the model.
     */
    struct State_ {
      double_t y0_; //!< Constant current
      double_t y1_ex_;  
      double_t y2_ex_;
      double_t y1_in_;
      double_t y2_in_;
      double_t y3_; //!< This is the membrane potential RELATIVE TO RESTING POTENTIAL.
      double_t y4_; //!< This is the effective change to threshold. // MOD

      int_t    r_;  //!< Number of refractory steps remaining

      State_();  //!< Default initialization
      
      void get(DictionaryDatum&, const Parameters_&) const;

      /** Set values from dictionary.
       * @param dictionary to take data from
       * @param current parameters
       * @param Change in reversal potential E_L specified by this dict
       */
      void set(const DictionaryDatum&, const Parameters_&, double);
    };    

    // ---------------------------------------------------------------- 

    /**
     * Buffers of the model.
     */
    struct Buffers_ {
      Buffers_(iaf_psc_alpha_dynth&);
      Buffers_(const Buffers_&, iaf_psc_alpha_dynth&);

      /** buffers and summs up incoming spikes/currents */
      RingBuffer ex_spikes_;
      RingBuffer in_spikes_;
      RingBuffer currents_;

      //! Logger for all analog data
      UniversalDataLogger<iaf_psc_alpha_dynth> logger_;
    };
    
    // ---------------------------------------------------------------- 

    /**
     * Internal variables of the model.
     */
    struct Variables_ { 
      /** Amplitude of the synaptic current.
	        This value is chosen such that a post-synaptic potential with
	        weight one has an amplitude of 1 mV.
       */
      double_t EPSCInitialValue_;
      double_t IPSCInitialValue_;
      int_t    RefractoryCounts_;
    
      double_t P11_ex_;
      double_t P21_ex_;
      double_t P22_ex_;
      double_t P31_ex_;
      double_t P32_ex_;
      double_t P11_in_;
      double_t P21_in_;
      double_t P22_in_;
      double_t P31_in_;
      double_t P32_in_;
      double_t P30_;
      double_t P33_;
      double_t P44_; // MOD
      double_t expm1_tau_m_;

      double_t weighted_spikes_ex_;
      double_t weighted_spikes_in_;
    };

    // Access functions for UniversalDataLogger -------------------------------

    //! Read out the real membrane potential
    double_t get_V_m_() const { return S_.y3_ + P_.U0_; }

    double_t get_weighted_spikes_ex_() const { return V_.weighted_spikes_ex_; }
    double_t get_weighted_spikes_in_() const { return V_.weighted_spikes_in_; }

    //! Read out the effective threshold
    double_t get_effective_th_() const { return S_.y4_ + P_.U0_ + P_.Theta_; } // MOD

    // Data members ----------------------------------------------------------- 
    
   /**
    * @defgroup iaf_psc_alpha_dynth_data
    * Instances of private data structures for the different types
    * of data pertaining to the model.
    * @note The order of definitions is important for speed.
    * @{
    */   
   Parameters_ P_;
   State_      S_;
   Variables_  V_;
   Buffers_    B_;
   /** @} */

   //! Mapping of recordables names to access functions
   static RecordablesMap<iaf_psc_alpha_dynth> recordablesMap_;
  };

inline
port nest::iaf_psc_alpha_dynth::check_connection(Connection& c, port receptor_type)
{
  SpikeEvent e;
  e.set_sender(*this);
  c.check_event(e);
  return c.get_target()->connect_sender(e, receptor_type);
}
  
inline
port iaf_psc_alpha_dynth::connect_sender(SpikeEvent&, port receptor_type)
{
  if (receptor_type != 0)
    throw UnknownReceptorType(receptor_type, get_name());
  return 0;
}
 
inline
port iaf_psc_alpha_dynth::connect_sender(CurrentEvent&, port receptor_type)
{
  if (receptor_type != 0)
    throw UnknownReceptorType(receptor_type, get_name());
  return 0;
}
 
inline
port iaf_psc_alpha_dynth::connect_sender(DataLoggingRequest& dlr, 
                                   port receptor_type)
{
  if (receptor_type != 0)
    throw UnknownReceptorType(receptor_type, get_name());
  return B_.logger_.connect_logging_device(dlr, recordablesMap_);
}

inline
void iaf_psc_alpha_dynth::get_status(DictionaryDatum &d) const
{
  P_.get(d);
  S_.get(d, P_);
  Archiving_Node::get_status(d);

  (*d)[names::recordables] = recordablesMap_.get_list();
}

inline
void iaf_psc_alpha_dynth::set_status(const DictionaryDatum &d)
{
  Parameters_ ptmp = P_;  // temporary copy in case of errors
  const double delta_EL = ptmp.set(d);                       // throws if BadProperty
  State_      stmp = S_;  // temporary copy in case of errors
  stmp.set(d, ptmp, delta_EL);                 // throws if BadProperty

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

#endif /* #ifndef IAF_PSC_ALPHA_DYNTH_H */
