/*
 *  ac_poisson_generator.h
 *
 *  (C) 2006-2008 by The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  This file is proprietary material of the NEST Initiative
 *
 */

#ifndef AC_POISSON_GENERATOR_H
#define AC_POISSON_GENERATOR_H

#include "nest.h"
#include "event.h"
#include "node.h"
#include "ring_buffer.h"
#include "stimulating_device.h"
#include "poisson_randomdev.h"
#include "connection.h"
#include "analog_data_logger.h"
#include <valarray>


namespace nest{

  class Network;

  /* BeginDocumentation
     Name: ac_poisson_generator - Generates sinusoidally modulated Poisson train.

     Description:
     ac_poisson_generator generates a sinusoidally modulate Poisson spike train.
     The instantaneous rate of the process is given by

     f(t) = A_0 + Sum_k A_k sin ( 2 pi freq_k t + phi_k )

     Parameters: 
     The following parameters can be set in the status dictionary.

     DC         double - DC input, scalar ( = A_0 ) 
     AC         array  - AC input, one value per component ( = A_k )
     Freq       array  - AC frequency, same size as AC, in Hz ( = freq_k )
     Phi        array  - AC phase, same size as AC, in Radian ( = phi_k )
  
 
     Remarks:
     - THIS GENERATOR IS NOT SUITABLE FOR PARALLEL SIMULATION, since
       targets on different virtual process will receive different
       spike trains.
     - AC, Freq, and Phi must be arrays of equal size.
     - The state of the generator is reset on calibration.
     - The generator does not support precise spike timing.
     - You can use the voltmeter to sample the rate of the generator.
     - The generator sends the same spike train to all of its targets.
     - The generator will create different trains if run at different
       temporal resolutions.

     Author:  July 2006, Hans E Plesser
     SeeAlso: ac_generator, ReichSpikeGen
  */

  class ac_poisson_generator: public Node
  {
    
  public:        
    
    ac_poisson_generator();
    ac_poisson_generator(const ac_poisson_generator&);

    // Import overloaded virtual function set to local scope. 
    using Node::connect_sender;
    using Node::handle;

    bool has_proxies() const {return false;}
    bool local_receiver() const { return true; }  // receives PotentialRequests

    port check_connection(Connection&, port);

    void handle(PotentialRequest &);
    port connect_sender(PotentialRequest &, port);

    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &) ;

  private:
    void init_state_(const Node&);
    void init_buffers_();
    void calibrate();

    void update(Time const &, const long_t, const long_t);
    
    struct Parameters_ {
      /** temporal frequency in radian/ms. */
      std::valarray<double_t> om_;
            
      /** phase in radian */
      std::valarray<double_t> phi_;
      
      /** DC amplitude */
      double_t dc_;
      
      /** AC amplitude */
      std::valarray<double_t> ac_;
      
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

      /**
       * Extract array from a dictionary.
       * @param d     dictionary to read from
       * @param dname name in dictionary
       * @param data  destination for data
       * @returns true if dname was in dictionary
       */
      bool extract_array_(const DictionaryDatum &d, 
			  const std::string& dname,
			  std::valarray<double>& data) const;
      
      /**
       * Store valarray in a dictionary.
       * @param d     dictionary to store in
       * @param dname name in dictionary
       * @param data  data to be stored; 
       *              copy by value to facilitate on-the-fly conversions
       */
      void store_array_(DictionaryDatum &d, 
			const std::string& dname,
			std::valarray<double> data) const;
    };

    // ------------------------------------------------------------

    struct Variables_ {
      librandom::PoissonRandomDev poisson_dev_;  //!< random deviate generator
      
      std::valarray<double_t> sins_;  //!< sin(h om) in propagator
      std::valarray<double_t> coss_;  //!< cos(h om) in propagator

    };

    // ------------------------------------------------------------
    
    /**
     * State
     */
    struct State_ {      
      Time last_spike_;  //!< time stamp of most recent spike fired.

      State_();  //!< Sets default state value

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary
      void set(const DictionaryDatum&, const Parameters_&);  //!< Set values from dicitonary
    };

    // ---------------------------------------------------------------- 

    /**
     * Buffers of the model.
     */
    struct Buffers_ {
      /** Number of oscillators */
      unsigned int N_osc_;

      /** 
       * Oscillators representing rate.
       * These must not be variables, since they must persist across 
       * repeated simulate calls.
       */
      std::valarray<double_t> state_oscillators_; 

      //! We abuse potential logging for making rate information available
      AnalogDataLogger<PotentialRequest> rates_;
    };

    // ------------------------------------------------------------

    StimulatingDevice<SpikeEvent> device_;

    Parameters_ P_;
    State_      S_;
    Variables_  V_;
    Buffers_    B_;

  };

  inline
    port ac_poisson_generator::check_connection(Connection& c, port receptor_type)
    {
      SpikeEvent e;
      e.set_sender(*this);
      c.check_event(e);
      return c.get_target()->connect_sender(e, receptor_type);
    }

  inline
    port ac_poisson_generator::connect_sender(PotentialRequest& pr, port receptor_type)
    {
      if (receptor_type != 0)
	throw UnknownReceptorType(receptor_type, get_name());
      B_.rates_.connect_logging_device(pr);
      return 0;
    }

  inline
    void ac_poisson_generator::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    S_.get(d);
    device_.get_status(d);
  }

  inline
    void ac_poisson_generator::set_status(const DictionaryDatum &d)
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

#endif // AC_POISSON_GENERATOR_H
