/*
 *  cython_neuron.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "exceptions.h"
#include "cython_neuron.h"
#include "network.h"
#include "dict.h"
#include "dictutils.h"
#include "numerics.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "universal_data_logger_impl.h"
#include "dictstack.h"

#include <string>
#include <limits>

#include <stdio.h>
/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */
nest::RecordablesMap<nest::cython_neuron> nest::cython_neuron::recordablesMap_;

namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
  template <>
  void RecordablesMap<cython_neuron>::create()
  {
    // use standard names whereever you can for consistency!
    insert_(names::V_m, &cython_neuron::get_V_m_);
  }
}

nest::cython_neuron::Buffers_::Buffers_(cython_neuron &n)
  : logger_(n)
{}

nest::cython_neuron::Buffers_::Buffers_(const Buffers_ &, cython_neuron &n)
  : logger_(n)
{}

/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */
nest::cython_neuron::cython_neuron()
  : Archiving_Node(),
    state_(new Dictionary()),
    B_(*this)
{
  recordablesMap_.create();
  pyObj = NULL;
}

nest::cython_neuron::cython_neuron(const cython_neuron& n)
  : Archiving_Node(n),
    state_(new Dictionary(*n.state_)),
    B_(n.B_, *this)
{
  init_state_(n);
  pyObj = NULL;
}

nest::cython_neuron::~cython_neuron()
{
   delete pyObj;
}


/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::cython_neuron::init_state_(const Node& proto)
{
  const cython_neuron& pr = downcast<cython_neuron>(proto);
  state_= DictionaryDatum(new Dictionary(*pr.state_));
}

void nest::cython_neuron::init_buffers_()
{
  B_.ex_spikes_.clear();       // includes resize
  B_.in_spikes_.clear();       // includes resize
  B_.currents_.clear();        // includes resize
  B_.logger_.reset(); 	       // includes resize
  Archiving_Node::clear_history();
}



void nest::cython_neuron::calibrate()
{
  B_.logger_.init();

  if(pyObj != NULL) {	  
	// Pointers to Standard Parameters passing
	pyObj->putStdParams(&currents, &in_spikes, &ex_spikes, &t_lag, &spike, &current_value);
	*spike = 0;
	*current_value = 0.0;
	pyObj->call_method("calibrate");

	if(state_->known(Name("optimized")) && (*state_)[Name("optimized")]) {
		optimized = true;
	}
	else {
		optimized = false;
	}
  }
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 */
void nest::cython_neuron::update(Time const & origin, const long_t from, const long_t to)
{
  for ( long_t lag = from ; lag < to ; ++lag )
  {
	*currents = B_.currents_.get_value(lag);
    *in_spikes = B_.in_spikes_.get_value(lag); // in spikes arriving at right border
    *ex_spikes = B_.ex_spikes_.get_value(lag); // ex spikes arriving at right border
    *t_lag = lag;


    if(this->optimized) {
		pyObj->call_update_optimized();
    }
    else {
		pyObj->call_update();
    }

    // threshold crossing
    if (*spike)
    {
      set_spiketime(Time::step(origin.get_steps()+lag+1));
      SpikeEvent se;
      network()->send(*this, se, lag);
    }

    if(*current_value != 0.0)
    {
        CurrentEvent ce;
		ce.set_current(*current_value);
		network()->send(*this, ce, lag);
    }

    B_.logger_.record_data(origin.get_steps()+lag);
  }
}


void nest::cython_neuron::handle(SpikeEvent & e)
{
  assert(e.get_delay() > 0);

  if(e.get_weight() > 0.0)
    B_.ex_spikes_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
                            e.get_weight() * e.get_multiplicity() );
  else
    B_.in_spikes_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
                            e.get_weight() * e.get_multiplicity() );
}

void nest::cython_neuron::handle(CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  const double_t I = e.get_current();
  const double_t w = e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
		                     w * I);
}

void nest::cython_neuron::handle(DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}

void nest::cython_neuron::setStatusCython()
{
	// The first setStatus iteration will remove the pyobject from the dictionary and put it in a different object
	static Name pyObjectName = Name("pyobject");
    if(pyObj == NULL && state_->known(pyObjectName) ) {
       ((PyObjectDatum*)(&(*(*state_)[pyObjectName])))->call_status_method(SET_STATUS_METHOD, &state_);
       pyObj = ((PyObjectDatum*)(&(*(*state_)[pyObjectName])))->clone();
       state_->remove(pyObjectName);
    }
    else if(pyObj != NULL) {
       pyObj->call_status_method(SET_STATUS_METHOD, &state_);
    }
}

void nest::cython_neuron::getStatusCython() const
{
    if(pyObj != NULL) {
		pyObj->call_status_method(GET_STATUS_METHOD, &state_);
    }
}


void nest::register_cython_model(nest::Network *net, std::string model)
{
	nest::register_model<cython_neuron>(*net, model.c_str());
}
