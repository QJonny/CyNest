/*
 *  ac_gamma_generator.h
 *
 *  Copyright (C) 2010 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  This file is proprietary material of the NEST Initiative
 *
 */

#ifndef AC_GAMMA_GENERATOR_H
#define AC_GAMMA_GENERATOR_H

#include "config.h"

#ifdef HAVE_GSL

#include "nest.h"
#include "event.h"
#include "node.h"
#include "stimulating_device.h"
#include "connection.h"
#include "universal_data_logger.h"

#include <vector>

namespace nest{

  class Network;

  /* BeginDocumentation
	Name: ac_gamma_generator - Generates sinusoidally modulated Gamma train.

	Description:
	  ac_gamma_generator generates a sinusoidally modulate Gamma spike train.
	  The instantaneous rate of the process is given by

			   f(t) = dc + ac sin ( 2 pi freq t + phi )

	Parameters:
	  The following parameters can be set in the status dictionary.

	  order      double - the gamma order, >= 1
	  dc         double - Mean firing rate in spikes/second, default: 0 s^-1
	  ac         double - Firing rate modulation amplitude in spikes/second (0 <= ac < dc), default: 0 s^-1
	  freq       double - Modulation frequency in Hz, default: 0 Hz
	  phi        double - Modulation phase in radian, default: 0

	Remarks:
	  - THIS GENERATOR IS NOT SUITABLE FOR PARALLEL SIMULATION,
		since targets on different virtual process will receive different spike
		trains.
	  - This generator is at present not tuned for speed.
	  - The generator will not give precise results for very high firing
		rates relative to the resolution of the simulation.
	  - The generator does not support precise spike timing.
	  - The generator sends the same spike train to all of its targets.
	  - The generator will create different trains if run at different
		temporal resolutions.

	Receives: DataLoggingRequest

	Sends: SpikeEvent

	References: Barbieri et al, J Neurosci Methods 105:25-37 (2001)
	FirstVersion: October 2007
	Author: Hans E Plesser, Thomas Heiberg
	SeeAlso: smp_generator
   */

  /**
   * AC Gamma Generator.
   * Generates AC-modulated inhomogeneous gamma process.
   * @todo The implementation is very quick and dirty and not tuned for performance at all. 
   * @note  The simulator works by calculating the hazard h(t) for each time step
   *  and comparing h(t) dt to a [0,1)-uniform number. The hazard is given by
   * $[
   *     h(t) = \frac{a \lambda(t) \Lambda(t)^{a-1} e^{-\Lambda(t)}}{\Gamma(a, \Lambda(t))}
   * $]
   * with
   * $[  \lambda(t) = dc + ac \sin ( 2 \pi f t + \phi ) $]
   * $[  \Lambda(t) = a \int_0^t \lambda(s) ds $]
   * and the incomplete Gamma function $\Gamma(a,z)$; $a$ is the order of the gamma function.
   *
   * @note This implementation includes an additional $a$ factor in the calculation of $\Lambda(t)$
   * and $h(t)$ in order to keep the mean rate constant with varying $a$
   */
  class ac_gamma_generator: public Node
  {
    
  public:        
    
    ac_gamma_generator();
    ac_gamma_generator(const ac_gamma_generator&);
    
    port check_connection(Connection&, port);

    void handle(DataLoggingRequest &);
    port connect_sender(DataLoggingRequest &, port);

    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &) ;

  private:
	void init_state_(const Node&);
	void init_buffers_();
	void calibrate();

	void update(Time const &, const long_t, const long_t);

	struct Parameters_ {
		/** Frequency in radian */
		double_t om_;

		/** phase in radian */
		double_t phi_;

		/** gamma order */
		double_t order_;

		/** DC amplitude */
		double_t dc_;

		/** AC amplitude */
		double_t ac_;

		Parameters_();  //!< Sets default parameter values
		Parameters_(const Parameters_& );
		Parameters_& operator=(const Parameters_& p); // Copy constructor EN

		void get(DictionaryDatum&) const;  //!< Store current values in dictionary

		/**
		* Set values from dicitonary.
		* @note State is passed so that the position can be reset if the
		*       spike_times_ vector has been filled with new data.
		*/
		void set(const DictionaryDatum&);
	};

	struct State_ {

      double_t rate_;    //!< current rate, kept for recording
      Time last_spike_;  //!< time stamp of most recent spike fired.

      State_();  //!< Sets default state value

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary
      void set(const DictionaryDatum&, const Parameters_&);  //!< Set values from dicitonary
    };
    // ------------------------------------------------------------

    // These friend declarations must be precisely here.
    friend class RecordablesMap<ac_gamma_generator>;
    friend class UniversalDataLogger<ac_gamma_generator>;

    // ----------------------------------------------------------------

    /**
     * Buffers of the model.
     */
    struct Buffers_ {
      Buffers_(ac_gamma_generator&);
      Buffers_(const Buffers_&, ac_gamma_generator&);
      UniversalDataLogger<ac_gamma_generator> logger_;
      double_t Lambda_hist_;   //!< integral over lambda(t) from S_.last_spike_ to Lamdba_t0_
      Time Lambda_t0_; //!< time to which Lamdba_hist_ is integrated
      Parameters_ P_prev_;  //!< parameter values prior to last SetStatus
    };

    // ------------------------------------------------------------

    struct Variables_ {
    };

    double_t get_rate_() const { return 1000.0 * S_.rate_; }

    // compute deltaLambda for given parameters from ta to tb
    double_t deltaLambda_(const Parameters_&, double_t, double_t);

    StimulatingDevice<SpikeEvent> device_;
    static RecordablesMap<ac_gamma_generator> recordablesMap_;

    Parameters_ P_;
    State_      S_;
    Variables_  V_;
    Buffers_    B_;
};

  inline
    port ac_gamma_generator::check_connection(Connection& c, port receptor_type)
  {
    SpikeEvent e;
    e.set_sender(*this);
    c.check_event(e);
    return c.get_target()->connect_sender(e, receptor_type);
  }
  inline
    port ac_gamma_generator::connect_sender(DataLoggingRequest& dlr,
				      port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return B_.logger_.connect_logging_device(dlr, recordablesMap_);
  }

  inline
    void ac_gamma_generator::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    S_.get(d);
    device_.get_status(d);
    (*d)[names::recordables] = recordablesMap_.get_list();
  }

  inline
    void ac_gamma_generator::set_status(const DictionaryDatum &d)
  {
    Parameters_ ptmp = P_;  // temporary copy in case of errors

    ptmp.set(d);                       // throws if BadProperty
    // We now know that ptmp is consistent. We do not write it back
    // to P_ before we are also sure that the properties to be set
    // in the parent class are internally consistent.
    device_.set_status(d);

    // if we get here, temporaries contain consistent set of properties
    P_ = ptmp;
  }

} // namespace

#endif // AC_GAMMA_GENERATOR_H

#endif //HAVE_GSL
