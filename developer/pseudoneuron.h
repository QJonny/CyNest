/*
 *  pseudoneuron.h
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2004 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  This file is CONFIDENTIAL and PROPRIETARY.
 *
 */

#ifndef PSEUDONEURON_H
#define PSEUDONEURON_H

#include <vector>
#include "nest.h"
#include "event.h"
#include "node.h"
#include "device.h"
#include "scheduler.h"
#include "exp_randomdev.h"

/*BeginDocumentation
Name: pseudoneuron - spit out and accept but forget incoming spikes

Description:
This is a pseudoneuron for debugging purposes only.  It emits a Poisson 
spike train at a given rate, but also receives incoming spikes.   The spikes
are just counted, though.
    
Parameters:
   The following parameters appear in the element's status dictionary:

   rate     - mean firing rate [Hz]. (double, var)

   origin   - Time origin for device timer in realtime. (double, var)
   start    - begin of device application in realtime with resp. to origin. (double, var)
   stop     - end of device application in realtime with resp. to origin. (double, var)

   spikes_rcvd - number of spikes received [read only]

Author: Plesser

SeeAlso: poisson_generator
*/

namespace nest{

  /** 
   * Pseudoneuron emitting spikes like a Poisson generator and counting 
   * incoming spikes.
   * @note The pseudoneuron emits spikes through the queue.
   */
  class pseudoneuron: public Device
  {
    
  public:        
    
    typedef Node base;
    
    pseudoneuron();
    pseudoneuron(const pseudoneuron&);

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
    
    void handle(SpikeEvent &);
    
    port connect_sender(SpikeEvent &) {return 0;}

  protected:
    
    void init();
    void calibrate();
    void update(thread, Time const &, const long_t, const long_t);
    
    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);

  private:
    librandom::ExpRandomDev exp_dev_;  //!< random deviate generator
    double_t rate_;   //!< process rate in Hz
    double_t h_;      //!< resolution of simulation
    std::vector<double_t> overflow_;      //!< overflows into next time bin
    long_t spikes_rcvd_;  //!< number of spikes received by neuron
};
 
}
#endif //PSEUDONEURON_H
