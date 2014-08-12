#ifndef IAF_PSC_EXP_CANON_H
#define IAF_PSC_EXP_CANON_H

#include "config.h"

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "../precise/slice_ring_buffer.h"
#include "connection.h"
#include "numerics.h"
#include "analog_data_logger.h"

#include <vector>

/* BeginDocumentation
Name: iaf_psc_exp_canon - Leaky integrate-and-fire neuron 
with exponentially shaped postsynaptic currents; canoncial implementation.

Description:
iaf_psc_exp_canon is the "canonical" implementatoin of the leaky 
integrate-and-fire model neuron with exponentially shaped postsynaptic 
currents in the sense of [1].

The canonical implementation handles neuronal dynamics in a locally 
event-based manner within coarse time grid defined by the minimum 
delay in the network, see [1]. Incoming spikes are applied at the 
precise moment of their arrival, while the precise time of outgoing 
spikes is determined by interpolation once a threshold crossing has 
been detected. Return from refractoriness occurs precisly at spike 
time plus refractory period.

This implementation is more complex than the plain iaf_psc_exp 
neuron, but achieves much higher precision. In particular, it does not 
suffer any binning of spike times to grid points. Depending on your 
application, the canonical application may provide superior overall 
performance given an accuracy goal; see [1] for details. Subthreshold 
dynamics are integrated using exact integration between events [2].

Remarks:
  The iaf_psc_exp_canon neuron does not accept CurrentEvent connections.
  This is because the present method for transmitting CurrentEvents in 
  NEST (sending the current to be applied) is not compatible with off-grid 
  currents, if more than one CurrentEvent-connection exists. Once CurrentEvents 
  are changed to transmit change-of-current-strength, this problem will 
  disappear and the canonical neuron will also be able to handle CurrentEvents.
  For now, the only way to inject a current is the built-in current I_e.

Parameters:
  The following parameters can be set in the status dictionary.

  E_L             double - Resting membrane potential in mV.
  c_m             double - Specific capacitance of the membrane in pF/mum^2
  tau_m           double - Membrane time constant in ms.
  t_ref           double - Duration of refractory period in ms.
  U_th            double - Spike threshold in mV.
  I_e             double - Constant input current in pA.
  U_min           double - Absolute lower value for the membrane potential.
  Interpol_Order  int    - Interpolation order for spike time:
                           0-none, 1-linear, 2-quadratic, 3-cubic
References:
  [1] Morrison A, Straube S, Plesser HE, Diesmann M. Exact subthreshold 
      integration with continuous spike times in discrete-time neural network 
      simulations. Neural Comput. 2007 Jan;19(1):47-79.
  [2] Rotter S, Diesmann M. Exact digital simulation of time-invariant linear 
      systems with applications to neuronal modeling. Biol Cybern. 1999 Nov;
      81(5-6):381-402.

Author: Hanuschkin, Kunkel

SeeAlso: iaf_psc_exp, iaf_psc_alpha_canon
*/

namespace nest
{
  /**
   * Leaky iaf neuron, exponential PSC synapses, canonical implementation,
   * exact timing of outgoing spikes.
   * @note Inherit privately from Node, so no classes can be derived
   * from this one.
   * @todo Implement current input in consistent way.
   */
  class iaf_psc_exp_canon : public Node
  {
  class Network;

  public:

    /**
     * Basic constructor.
     * This constructor should only be used by GenericModel to create
     * model prototype instances.
     */
    iaf_psc_exp_canon();

    /**
     * Copy constructor.
     * GenericModel::allocate_() uses the copy constructor to clone
     * actual model instances from the prototype instance.
     * @note The copy constructor MUST NOT be used to create nodes based
     * on nodes that have been placed in the network.
     */
    iaf_psc_exp_canon(iaf_psc_exp_canon const &);

    /**
     * Import sets of overloaded virtual functions.
     * @see Technical Issues / Virtual Functions: Overriding, Overloading, and Hiding
     */

    using Node::connect_sender;
    using Node::handle;

    port check_connection(Connection &, port);
    port connect_sender(SpikeEvent &, port);
    port connect_sender(PotentialRequest &, port);

    void get_status(DictionaryDatum &) const;
    void set_status(DictionaryDatum const &);

    // uses off_grid events
    bool is_off_grid() const
    {
      return true;
    }

    void handle(SpikeEvent &);
    void handle(PotentialRequest &);

  private:

    /**
     * @name Interface functions
     * @note These functions are private, so that they can be accessed
     * only through a Node*.
     */
    //@{

    void init_state_(Node const &proto);
    void init_buffers_();
    void calibrate();

    /**
     * Time Evolution Operator.
     *
     * update() promotes the state of the neuron from origin+from to origin+to.
     * It does so in steps of the resolution h.  Within each step, time is
     * advanced from event to event, as retrieved from the spike queue.
     *
     * Return from refractoriness is handled as a special event in the
     * queue, which is marked by a weight that is GSL_NAN. This greatly simplifies
     * the code.
     *
     * For steps, during which no events occur, the precomputed propagator matrix
     * is used.  For other steps, the propagator matrix is computed as needed.
     *
     * While the neuron is refractory, membrane potential (y2_) is
     * clamped to U_reset_.
     */
    void update(Time const &origin, long_t const from, long_t const to);

    //@}

    Time get_spiketime() const;
    void set_spiketime(Time const &);

    /**
     * Propagate neuron state.
     * Propagate the neuron's state by dt.
     * @param dt Interval over which to propagate
     */
    void propagate_(double_t const dt);

    /**
     * Emit a single spike caused by DC current in absence of spike input.
     * Emits a single spike and reset neuron given that the membrane
     * potential was below threshold at the beginning of a mini-timestep
     * and above afterwards.
     *
     * @param origin  Time stamp at beginning of slice
     * @param lag     Time step within slice
     * @param t0      Beginning of mini-timestep
     * @param dt      Duration of mini-timestep
     */
    void emit_spike_(Time const &origin, long_t const lag, 
                       double_t const t0, double_t const dt);

    /**
     * Emit a single spike at a precisely given time.
     *
     * @param origin        Time stamp at beginning of slice
     * @param lag           Time step within slice
     * @param spike_offset  Time offset for spike
     */
    void emit_instant_spike_(Time const &origin, long_t const lag, 
                               double_t const spike_offset);

    /**
     * @name Threshold-crossing interpolation
     * These functions determine the time of threshold crossing using
     * interpolation, one function per interpolation
     * order. thresh_find() is the driver function and the only one to
     * be called directly. 
     */
    //@{

    /** Interpolation orders. */
    enum interpOrder { NO_INTERPOL, LINEAR, QUADRATIC, CUBIC, END_INTERP_ORDER };

    /**
     * Localize threshold crossing.
     * Driver function to invoke the correct interpolation function
     * for the chosen interpolation order.
     * @param   double_t length of interval since previous event
     * @returns time from previous event to threshold crossing
     */
    double_t thresh_find_(double_t const) const;
    double_t thresh_find1_(double_t const) const;
    double_t thresh_find2_(double_t const) const;
    double_t thresh_find3_(double_t const) const;

    //@}

    // ---------------------------------------------------------------- 

    /** 
     * Independent parameters of the model. 
     */
    struct Parameters_
    {
      double_t tau_m_;       //!< Membrane time constant in [ms]
      double_t tau_ex_;      //!< Time constant of exc. synaptic current in [ms]
      double_t tau_in_;      //!< Time constant of inh. synaptic current in [ms]
      double_t c_m_;         //!< Membrane capacitance in [pF]
      double_t t_ref_;       //!< Refractory period in [ms]
      double_t E_L_;         //!< Resting potential in [mV]
      double_t I_e_;         //!< External DC in [pA]
      double_t U_th_;        //!< Threshold (rel. to resting potential)
      double_t U_min_;       //!< Lower bound (rel. to resting potential)
      double_t U_reset_;     //!< Reset potential (rel. to resting potential)
      interpOrder Interpol_; //!< Interpolation order

      Parameters_(); //!< Sets default parameter values

      void get(DictionaryDatum &) const; //!< Store current values in dictionary
      void set(DictionaryDatum const &); //!< Set values from dicitonary
    };

    // ---------------------------------------------------------------- 

    /**
     * State variables of the model.
     */
    struct State_
    {
      double_t y1_ex_; //!< Exc. exponetial current
      double_t y1_in_; //!< Inh. exponetial current
      double_t y2_;    //!< Membrane potential (rel. to resting potential)

      bool     is_refractory_;     //!< True while refractory
      long_t   last_spike_step_;   //!< Time stamp of most recent spike
      double_t last_spike_offset_; //!< Offset of most recent spike 

      State_(); //!< Default initialization
      
      void get(DictionaryDatum &, Parameters_ const &) const;
      void set(DictionaryDatum const &, Parameters_ const &);
    };

    // ---------------------------------------------------------------- 

    /**
     * Buffers of the model.
     */
    struct Buffers_
    {
      /**
       * Queue for incoming events.
       * @note Handles also pseudo-events marking return from refractoriness.
       */
      SliceRingBuffer events_;

      /**
       * Buffer for membrane potential.
       * @note Rows are written alternatingly.
       */
      AnalogDataLogger<PotentialRequest> potentials_; 
    };

    // ---------------------------------------------------------------- 

    /**
     * Internal variables of the model.
     */
    struct Variables_
    { 
      double_t h_ms_;             //!< Time resolution in [ms]
      long_t   refractory_steps_; //!< Refractory time in steps
      double_t expm1_tau_m_;      //!< exp(-h/tau_m) - 1 
      double_t expm1_tau_ex_;     //!< exp(-h/tau_ex) - 1
      double_t expm1_tau_in_;     //!< exp(-h/tau_in) - 1
      double_t P20_;              //!< Progagator matrix element, 2nd row
      double_t P21_ex_;           //!< Progagator matrix element, 2nd row, exc.
      double_t P21_in_;           //!< Progagator matrix element, 2nd row, inh.
      double_t y1_ex_before_;     //!< y1_ex at beginning of ministep, for interpolation
      double_t y1_in_before_;     //!< y1_in at beginning of ministep, for interpolation
      double_t y2_before_;        //!< y2_ at beginning of ministep, for interpolation
    };
    
    // ---------------------------------------------------------------- 

     /**
      * @defgroup iaf_psc_exp_data
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

  };

inline
port iaf_psc_exp_canon::check_connection(Connection &c, port receptor_type)
{
  SpikeEvent e;

  e.set_sender(*this);
  c.check_event(e);
  return c.get_target()->connect_sender(e, receptor_type);
}

inline
port iaf_psc_exp_canon::connect_sender(SpikeEvent &, port receptor_type)
{
  if ( receptor_type != 0 )
    throw UnknownReceptorType(receptor_type, get_name());
  return 0;
}

inline
port iaf_psc_exp_canon::connect_sender(PotentialRequest &pr, port receptor_type)
{
  if ( receptor_type != 0 )
    throw UnknownReceptorType(receptor_type, get_name());
    B_.potentials_.connect_logging_device(pr);
  return 0;
}

inline
void iaf_psc_exp_canon::get_status(DictionaryDatum &d) const
{
  P_.get(d);
  S_.get(d, P_);
}

inline
void iaf_psc_exp_canon::set_status(const DictionaryDatum &d)
{
  Parameters_ ptmp = P_; // temporary copy in case of errors
  ptmp.set(d);           // throws if BadProperty

  State_ stmp = S_;  // temporary copy in case of errors
  stmp.set(d, ptmp); // throws if BadProperty

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
  S_ = stmp;
}

} // namespace

#endif //IAF_PSC_EXP_CANON_H
