/*
 *  environment_gauss.h
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2004-2008 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */


/*BeginDocumentation
Name: environment_gauss - provides DC input current

Description: The DC-Generator provides a constant DC Input
to the connected node. The unit of the current is pA.

Examples: The dc current can be altered in the following way:
   /environment_gauss Create /dc_gen Set         % Creates a environment_gauss, which is a node
   dc_gen GetStatus info                    % View properties (amplitude is 0)
   dc_gen << /amplitude 1500. >> SetStatus
   dc_gen GetStatus info                    % amplitude is now 1500.0

Note: The environment_gauss is rather inefficient, since it needs to
      send the same current information on each time step. If you
      only need a constant bias current into a neuron, you should 
      set it directly in the neuron, e.g., iaf_neuron.
      
Author: docu by Sirko Straube

SeeAlso: Device, StimulatingDevice

*/ 

#ifndef ENVIRONMENT_GAUSS_H
#define ENVIRONMENT_GAUSS_H

#include <vector>
#include <fstream>
#include "nest.h"
#include "event.h"
#include "node.h"
#include "ring_buffer.h"
#include "connection.h"
//#include "fstream.h"
//#include "exeptions.h"
//#include "stimulating_device.h"
#include "archiving_node.h"
#include "namedatum.h"


namespace nest
{

  class Network; 

  typedef lockPTR<std::ostream> ostreamPtr;
  /**
   * DC current generator.
   *
   * @ingroup Devices
   */
  class environment_gauss : public Archiving_Node//, public StimulatingDevice<CurrentEvent>
  {
    
  public:        
    
    typedef Node base;

    environment_gauss();
    environment_gauss(const environment_gauss&);

    bool has_proxies() const {return false;}
    bool local_receiver() const {return false;}

    using Node::connect_sender;
    using Node::handle;
    port check_connection(Connection&, port);
    port connect_sender(SpikeEvent &, port);
    
    void handle(SpikeEvent &);

    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);

  protected:
    void init_parameters_(Node const*); //?
    //void calibrate();
    void event_hook(DSCurrentEvent&);
    void init_dynamic_state_(Node const*); //?
      //void update(Time const &, const long_t, const long_t);
    
  private:

    void move(int_t direction);
    void go_north();
    void go_south();
    void go_east();
    void go_west();
    bool pos_valid();
    int_t pos(int_t x, int_t y, int_t columns);
    void special_moves();
    void random_position();
    void put_agent(int x, int y);
    void plot_agent();

    void init_state_(const Node& proto);
    void init_buffers_();
    void calibrate();
    void finalize();

    void update(Time const &, const long_t, const long_t);
    const std::string build_filename_() const;
    //bool to_file()   const { return P_.to_file_;   }

    //-----------------------------------------------------------
    /**
     * Independent parameters of the model.
     */

    struct Parameters_ {
      /** Max amplitude in [pA] for stimulating the state neurons. */
      double_t    amp_;

      /** width in rad of gaussian stimulation of the state neurons. */
      double_t    sigma_;

      /** Amplitude in [pA] for inhibiting the actor neurons. */
      double_t I_inhib_;
    
      /** Number of neurons belonging to one state.*/
      long_t nps_;

      /** Number of states the agent jumps in one movie.*/
      int_t agent_step_;

      /** x_dim of the grid. */
      int_t x_dim_;

      /** y_dim of the grid. */
      int_t y_dim_;

      /** dim of the border of the grid (how many states/neurons form border). */
      int_t border_dim_;

      /** Time period of inhibiting the actor neurons. */
      double_t refractor_;

      /** Vector of neuron id's representing the state. */ 

      std::vector<long_t> input_neurons_; 
      
      long_t neuron_west_;
      long_t neuron_east_;
      long_t neuron_south_;
      long_t neuron_north_;

      /** Vector of neuron id's representing the critic. */
      std::vector<long_t> critic_neurons_;

      /** x and y position of reward state */
      long_t specialx_, specialy_;

      /** x and y position to which agent is put after it found the reward. -1,-1 means that it is put to a random position */
      long_t warpx_, warpy_; 

      /** amount of reward in pA */
      double_t special_reward_;

      /** seed */
      long_t seed_;

      /** if true, close stream in init_buffers()*/
      bool close_on_reset_;

      /** the filename for recording the steps performed by the agent*/
      std::string filename_;

      /** a user-defined label for filename*/
      std::string label_;

      /** if true, finalize() shall flush the stream */
      bool flush_after_simulate_;

      /** if true, finalize() shall close the stream */
      bool close_after_simulate_;



      Parameters_();  //!< Sets default parameter values

      void get(DictionaryDatum&) const; //!< Store current values in dictionary
      void set(const DictionaryDatum&);  //!< Set values from dicitonary
      //bool to_file_;       //!< true if recorder writes its output to a file
    };

    //----------------------------------------------------------------

    /** 
     * Buffers of the model.
     */
    
    struct Buffers_ {
      /** buffer of incoming spikes */
      RingBuffer spikes_;
      std::ofstream fs_; //!< the file to write the recorded data to
      int_t state_; //current state
      int_t prev_state_; //previous state 
      double_t ref_period_;
      double_t h;
      bool wall_;
      int_t agentx_;
      int_t agenty_;
      int_t actx_;
      int_t acty_;
      double_t Gausscurr_;
      double_t sumGausscurr_;
      int_t old_agent_[2];
      double_t reward_;
      double_t effamp_;
      double_t sumGauss_;     
      //std::vector<double_t> grid_mat_;
      std::vector<std::vector<double_t> > grid_mat_; /** Input Matrix. */ 
    };

    //-----------------------------------------------------------------

    /**
     * Internal variables of the model.
     */

    struct Variables_ {
      void open_file_();
      bool new_file_;
    };

    // ---------------------------------------------------------------

    Parameters_ P_;
    Variables_ V_;
    Buffers_ B_;

    /** 
     * Node to which environment_gauss instance belongs.
     * @todo This reference must be replaced by a pointer, see #299.
     */
    const Node& node_;
  };

  

inline  
port environment_gauss::check_connection(Connection& c, port receptor_type)
    {
      DSCurrentEvent e;
      e.set_sender(*this);
      c.check_event(e);
      return c.get_target()->connect_sender(e, receptor_type);
    }

inline
  port environment_gauss::connect_sender(SpikeEvent&, port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return 0;
  }


inline
void environment_gauss::get_status(DictionaryDatum &d) const
{
  P_.get(d);
  (*d)[names::type] = LiteralDatum(names::other);
}

/*
inline
void environment_gauss::set_status(const DictionaryDatum &d)
{
  Parameters_ ptmp = P_;  // temporary copy in case of errors
  ptmp.set(d);                       // throws if BadProperty
  
  P_ = ptmp;
 
}
*/
  
  
  
} // namespace

#endif /* #ifndef ENVIRONMENT_GAUSS_H */

