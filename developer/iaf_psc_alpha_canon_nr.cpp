/*
 *  iaf_psc_alpha_canon_nr.cpp
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2004-2006 by
 *  The NEST Initiative 
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#include "iaf_psc_alpha_canon_nr.h"

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include "analog_data_logger_impl.h"

#include <limits>

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
nest::iaf_psc_alpha_canon_nr::Parameters_::Parameters_()
  : tau_m_  ( 10.0    ),  // ms
    tau_syn_(  2.0    ),  // ms
    c_m_    (250.0    ),  // pF
    t_ref_  (  2.0    ),  // ms
    E_L_    (-70.0    ),  // mV
    I_e_    (  0.0    ),  // pA
    U_th_   (-55.0-E_L_),  // mV, rel to E_L_
    U_min_  (-std::numeric_limits<double_t>::infinity()),  // mV
    U_reset_(-70.0-E_L_),  // mV, rel to E_L_
    Interpol_(iaf_psc_alpha_canon_nr::CUBIC)
{}

nest::iaf_psc_alpha_canon_nr::State_::State_()
  : y0_(0.0),
    y1_(0.0),
    y2_(0.0),
    y3_(0.0),  
    is_refractory_(false),
    last_spike_step_(-1),
    last_spike_offset_(0.0)
{}

/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::iaf_psc_alpha_canon_nr::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d, names::E_L, E_L_);
  def<double>(d, names::I_e, I_e_);
  def<double>(d, names::V_th, U_th_+E_L_);
  def<double>(d, names::V_min, U_min_+E_L_);
  def<double>(d, names::V_reset, U_reset_+E_L_);
  def<double>(d, names::C_m, c_m_);
  def<double>(d, names::tau_m, tau_m_);
  def<double>(d, names::tau_syn, tau_syn_);
  def<double>(d, names::t_ref, t_ref_);
  def<long>(d, names::Interpol_Order, Interpol_);
}

void nest::iaf_psc_alpha_canon_nr::Parameters_::set(const DictionaryDatum& d)
{
  updateValue<double>(d, names::tau_m, tau_m_);
  updateValue<double>(d, names::tau_syn, tau_syn_);
  updateValue<double>(d, names::C_m, c_m_);
  updateValue<double>(d, names::t_ref, t_ref_);
  updateValue<double>(d, names::E_L, E_L_);
  updateValue<double>(d, names::I_e, I_e_);

  if (updateValue<double>(d, names::V_th, U_th_)) 
    U_th_ -= E_L_;

  if (updateValue<double>(d, names::V_min, U_min_)) 
    U_min_ -= E_L_;

  if (updateValue<double>(d, names::V_reset, U_reset_)) 
    U_reset_ -= E_L_;
  
  long_t tmp;
  if ( updateValue<long_t>(d, names::Interpol_Order, tmp) )
  {
    if ( NO_INTERPOL <= tmp && tmp < END_INTERP_ORDER )
      Interpol_ = static_cast<interpOrder>(tmp);
    else
      throw BadProperty("Invalid interpolation order. "
                        "Valid orders are 0, 1, 2, 3, 4.");
  }

  if ( U_reset_ >= U_th_ )
    throw BadProperty("Reset potential must be smaller than threshold.");

  if ( U_reset_ < U_min_ )
    throw BadProperty("Reset potential must be greater equal minimum potential.");
    
  if ( c_m_ <= 0 )
    throw BadProperty("Capacitance must be strictly positive.");
    
  if ( t_ref_ < 0 )
    throw BadProperty("Refractory time must not be negative.");
    
  if ( tau_m_ <= 0 || tau_syn_ <= 0 )
    throw BadProperty("All time constants must be strictly positive.");
}

void nest::iaf_psc_alpha_canon_nr::State_::get(DictionaryDatum &d, 
                                            const Parameters_& p) const
{
  def<double>(d, names::V_m, y3_ + p.E_L_); // Membrane potential
  def<double>(d, names::t_spike, Time(Time::step(last_spike_step_)).get_ms());
  def<double>(d, names::offset, last_spike_offset_);
  def<bool>(d, names::is_refractory, is_refractory_);
}

void nest::iaf_psc_alpha_canon_nr::State_::set(const DictionaryDatum& d, const Parameters_& p)
{
  if ( updateValue<double>(d, names::V_m, y3_) )
    y3_ -= p.E_L_;
}

/* ---------------------------------------------------------------- 
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::iaf_psc_alpha_canon_nr::iaf_psc_alpha_canon_nr()
  : Node(), 
    P_(), 
    S_()
{}

nest::iaf_psc_alpha_canon_nr::iaf_psc_alpha_canon_nr(const iaf_psc_alpha_canon_nr& n)
  : Node(n), 
    P_(n.P_), 
    S_(n.S_)
{}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::iaf_psc_alpha_canon_nr::init_state_(const Node& proto)
{
  const iaf_psc_alpha_canon_nr& pr = downcast<iaf_psc_alpha_canon_nr>(proto);
  S_ = pr.S_;
}

void nest::iaf_psc_alpha_canon_nr::init_buffers_()
{
  B_.events_.resize();
  B_.events_.clear(); 
  B_.currents_.clear();  // includes resize
  B_.potentials_.clear_data(); // includes resize
}

void nest::iaf_psc_alpha_canon_nr::calibrate()
{
  V_.h_ms_ = Time::get_resolution().get_ms();

  V_.PSCInitialValue_ = 1.0 * numerics::e / P_.tau_syn_;

  V_.gamma_    = 1/P_.c_m_ /  ( 1/P_.tau_syn_ - 1/P_.tau_m_ );
  V_.gamma_sq_ = 1/P_.c_m_ /( ( 1/P_.tau_syn_ - 1/P_.tau_m_ ) 
                            * ( 1/P_.tau_syn_ - 1/P_.tau_m_ ) );

  // pre-compute matrix for full time step
  V_.expm1_tau_m_   = numerics::expm1(-V_.h_ms_/P_.tau_m_);
  V_.expm1_tau_syn_ = numerics::expm1(-V_.h_ms_/P_.tau_syn_);
  V_.P30_ = -P_.tau_m_ / P_.c_m_ * V_.expm1_tau_m_;
  V_.P31_ = V_.gamma_sq_ * V_.expm1_tau_m_ - V_.gamma_sq_ * V_.expm1_tau_syn_  
           - V_.h_ms_ * V_.gamma_ * V_.expm1_tau_syn_ - V_.h_ms_ * V_.gamma_; 
  V_.P32_ = V_.gamma_ * V_.expm1_tau_m_ - V_.gamma_ * V_.expm1_tau_syn_;

  // t_ref_ is the refractory period in ms
  // refractory_steps_ is the duration of the refractory period in whole
  // steps, rounded down
  V_.refractory_steps_ = Time(Time::ms(P_.t_ref_)).get_steps();
  assert(V_.refractory_steps_ >= 0);  // since t_ref_ >= 0, this can only fail in error
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void nest::iaf_psc_alpha_canon_nr::update(Time const & origin, 
				       const long_t from, const long_t to)
{
  assert(to >= 0);
  assert(static_cast<delay>(from) < Scheduler::get_min_delay());
  assert(from < to);
  
  // at start of slice, tell input queue to prepare for delivery
  if ( from == 0 )
    B_.events_.prepare_delivery();

  /* Neurons may have been initialized to superthreshold potentials.
     We need to check for this here and issue spikes at the beginning of
     the interval.
  */
  if ( S_.y3_ >= P_.U_th_ )
    emit_instant_spike_(origin, from, 
			V_.h_ms_*(1-std::numeric_limits<double_t>::epsilon()));

  for ( long_t lag = from ; lag < to ; ++lag )
    {
      // time at start of update step
      const long_t T = origin.get_steps() + lag;
      // if neuron returns from refractoriness during this step, place
      // pseudo-event in queue to mark end of refractory period
      if (S_.is_refractory_ && (T + 1 - S_.last_spike_step_ == V_.refractory_steps_))
	      B_.events_.add_refractory(T, S_.last_spike_offset_);

      // save state at beginning of interval for spike-time interpolation
      V_.y0_before_ = S_.y0_;
      V_.y1_before_ = S_.y1_;
      V_.y2_before_ = S_.y2_;
      V_.y3_before_ = S_.y3_;

      // get first event
      double_t ev_offset;
      double_t ev_weight;
      bool     end_of_refract;

      if ( !B_.events_.get_next_spike(T, ev_offset, ev_weight, end_of_refract) )
	{ // No incoming spikes, handle with fixed propagator matrix.
	  // Handling this case separately improves performance significantly
	  // if there are many steps without input spikes.

	  // update membrane potential
	  if ( !S_.is_refractory_ )
	    {
	      S_.y3_ = V_.P30_*(P_.I_e_+S_.y0_) + V_.P31_*S_.y1_ + V_.P32_*S_.y2_ +
		      V_.expm1_tau_m_ * S_.y3_  + S_.y3_;

	      // lower bound of membrane potential
	      S_.y3_ = ( S_.y3_ < P_.U_min_ ? P_.U_min_ : S_.y3_); 
	    }

	  // update synaptic currents
	  S_.y2_ = V_.expm1_tau_syn_ * V_.h_ms_ * S_.y1_  
	        + V_.expm1_tau_syn_ * S_.y2_ + V_.h_ms_ * S_.y1_  +  S_.y2_;
	  S_.y1_ = V_.expm1_tau_syn_ * S_.y1_ + S_.y1_;
	  
	  /* The following must not be moved before the y1_, y2_ update,
	     since the spike-time interpolation within emit_spike_ depends 
	     on all state variables having their values at the end of the
	     interval. 
	  */
	  if ( S_.y3_ >= P_.U_th_ )
	    emit_spike_(origin, lag, 0, V_.h_ms_);
	}
      else
	{
	  // We only get here if there is at least on event, 
	  // which has been read above.  We can therefore use 
	  // a do-while loop.

	  // Time within step is measured by offsets, which are h at the beginning
	  // and 0 at the end of the step.
	  double_t last_offset = V_.h_ms_;  // start of step

	  do {
	    // time is measured backward: inverse order in difference
	    const double_t ministep = last_offset - ev_offset;

	    propagate_(ministep);

	    // check for threshold crossing during ministep
	    // this must be done before adding the input, since
	    // interpolation requires continuity
	    if ( S_.y3_ >= P_.U_th_ ) 
	      emit_spike_(origin, lag, V_.h_ms_-last_offset, ministep);

	    // handle event
	    if ( end_of_refract )
	      S_.is_refractory_ = false;   // return from refractoriness
	    else
	      S_.y1_ += V_.PSCInitialValue_ * ev_weight;  // spike input

  	    // store state    //// STATE UPDATE HERE?
            V_.y0_before_ = S_.y0_;
            V_.y1_before_ = S_.y1_;
            V_.y2_before_ = S_.y2_;
            V_.y3_before_ = S_.y3_;

//old	    V_.y2_before_ = S_.y2_;
//	    V_.y3_before_ = S_.y3_;
	    last_offset = ev_offset;

	  } while ( B_.events_.get_next_spike(T, ev_offset, ev_weight, 
					   end_of_refract) );

	  // no events remaining, plain update step across remainder 
	  // of interval
	  if ( last_offset > 0 )  // not at end of step, do remainder
	    {
	      propagate_(last_offset);
	      if ( S_.y3_ >= P_.U_th_ ) 
		emit_spike_(origin, lag, V_.h_ms_-last_offset, last_offset);
	    }
	}  // else

      // Set new input current. The current change occurs at the
      // end of the interval and thus must come AFTER the threshold-
      // crossing interpolation
      S_.y0_  = B_.currents_.get_value(lag);


      // voltage logging
      B_.potentials_.record_data(origin.get_steps()+lag, S_.y3_ + P_.E_L_);
    }  // from lag = from ...
  
}


//function handles exact spike times
void nest::iaf_psc_alpha_canon_nr::handle(SpikeEvent & e)
{
  assert(e.get_delay() > 0 );

  /* We need to compute the absolute time stamp of the delivery time
     of the spike, since spikes might spend longer than min_delay_
     in the queue.  The time is computed according to Time Memo, Rule 3.
  */
  const long_t Tdeliver = e.get_stamp().get_steps() + e.get_delay() - 1;
  B_.events_.add_spike(e.get_rel_delivery_steps(network()->get_slice_origin()), 
		    Tdeliver, e.get_offset(), e.get_weight() * e.get_multiplicity());
}

void nest::iaf_psc_alpha_canon_nr::handle(CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  const double_t c=e.get_current();
  const double_t w=e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
		      w * c);
}

void nest::iaf_psc_alpha_canon_nr::handle(PotentialRequest& e)
{
  B_.potentials_.handle(*this, e);
}

// auxiliary functions ---------------------------------------------

inline 
void nest::iaf_psc_alpha_canon_nr::set_spiketime(Time const & now)
{
  S_.last_spike_step_ = now.get_steps();
}

void nest::iaf_psc_alpha_canon_nr::propagate_(const double_t dt)
{
  const double_t ps_e_TauSyn = numerics::expm1(-dt/P_.tau_syn_); // needed in any case 

  // y3_ remains unchanged at 0.0 while neuron is refractory
  if ( !S_.is_refractory_ )
    {
      const double_t ps_e_Tau = numerics::expm1(-dt/P_.tau_m_);
      const double_t ps_P30   = -P_.tau_m_/P_.c_m_*ps_e_Tau;
      const double_t ps_P31   =  V_.gamma_sq_ * ps_e_Tau - V_.gamma_sq_ * ps_e_TauSyn 
	                         - dt*V_.gamma_*ps_e_TauSyn - dt*V_.gamma_;
      const double_t ps_P32   =  V_.gamma_*ps_e_Tau - V_.gamma_* ps_e_TauSyn;
      S_.y3_ = ps_P30 * (P_.I_e_ + S_.y0_) + ps_P31 * S_.y1_ 
                 + ps_P32 * S_.y2_ + ps_e_Tau * S_.y3_  + S_.y3_;

      // lower bound of membrane potential
      S_.y3_ = ( S_.y3_< P_.U_min_ ? P_.U_min_ : S_.y3_); 
    }

  // now the synaptic components
  S_.y2_ = ps_e_TauSyn * dt * S_.y1_  + ps_e_TauSyn * S_.y2_ + dt * S_.y1_  +  S_.y2_;
  S_.y1_ = ps_e_TauSyn * S_.y1_ + S_.y1_;

  return;
}

void nest::iaf_psc_alpha_canon_nr::emit_spike_(Time const &origin, const long_t lag, 
					    const double_t t0,  const double_t dt)
{
  // we know that the potential is subthreshold at t0, super at t0+dt

  // compute spike time relative to beginning of step
  const double_t spike_offset = V_.h_ms_ - (t0 + thresh_find_(dt));
  set_spiketime(Time::step(origin.get_steps() + lag + 1));
  S_.last_spike_offset_ = spike_offset;  

  // reset neuron and make it refractory
  S_.y3_ = P_.U_reset_;
  S_.is_refractory_ = true; 

  // send spike
  SpikeEvent se;
  se.set_offset(spike_offset);
  network()->send(*this, se, lag);

  return;
}

void nest::iaf_psc_alpha_canon_nr::emit_instant_spike_(Time const &origin, const long_t lag,
						    const double_t spike_offs) 
{
  assert(S_.y3_ >= P_.U_th_);  // ensure we are superthreshold
 
  // set stamp and offset for spike
  set_spiketime(Time::step(origin.get_steps() + lag + 1));
  S_.last_spike_offset_ = spike_offs;

  // reset neuron and make it refractory
  S_.y3_ = P_.U_reset_;
  S_.is_refractory_ = true; 

  // send spike
  SpikeEvent se;
  se.set_offset(S_.last_spike_offset_);
  network()->send(*this, se, lag);

  return;
}

// finds threshpassing
inline
nest::double_t nest::iaf_psc_alpha_canon_nr::thresh_find_(double_t const dt) const
{
  switch (P_.Interpol_) {
    case NO_INTERPOL: return dt;
    case LINEAR     : return thresh_find1_(dt);
    case QUADRATIC  : return thresh_find2_(dt);
    case CUBIC      : return thresh_find3_(dt);
    case NR         : return thresh_find4_(dt);
    default: 
      throw BadProperty("Invalid interpolation order in iaf_psc_alpha_canon_nr.");
  }
  return 0;
}

// finds threshpassing via linear interpolation
nest::double_t nest::iaf_psc_alpha_canon_nr::thresh_find1_(double_t const dt) const
{
  double_t tau = ( P_.U_th_ - V_.y3_before_ ) * dt / ( S_.y3_ - V_.y3_before_ );
  return tau;
}
  
// finds threshpassing via quadratic interpolation
nest::double_t nest::iaf_psc_alpha_canon_nr::thresh_find2_(double_t const dt) const
{
  const double_t h_sq = dt * dt;
  const double_t derivative = - V_.y3_before_/P_.tau_m_ + (P_.I_e_ + V_.y0_before_ + V_.y2_before_)/P_.c_m_;
  
  const double_t a = (-V_.y3_before_/h_sq) + (S_.y3_/h_sq) - (derivative/dt);
  const double_t b = derivative;
  const double_t c = V_.y3_before_;
  
  const double_t sqr_ = std::sqrt(b*b - 4*a*c + 4*a*P_.U_th_);
  const double_t tau1 = (-b + sqr_) / (2*a);
  const double_t tau2 = (- b - sqr_) / (2*a);
      
  if (tau1 >= 0)
    return tau1;
  else if (tau2 >= 0)
    return tau2;
  else
    return thresh_find1_(dt);
}
 
nest::double_t nest::iaf_psc_alpha_canon_nr::thresh_find3_(double_t const dt) const
{
  const double_t h_ms = dt;
  const double_t h_sq = h_ms*h_ms;
  const double_t h_cb = h_sq*h_ms;
  
  const double_t deriv_t1 = - V_.y3_before_/P_.tau_m_ + (P_.I_e_ + V_.y0_before_ + V_.y2_before_)/P_.c_m_;
  const double_t deriv_t2 = - S_.y3_/P_.tau_m_ + (P_.I_e_ + S_.y0_ + S_.y2_)/P_.c_m_;
      
  const double_t w3_ = (2 * V_.y3_before_ / h_cb) - (2 * S_.y3_ / h_cb) 
                       + ( deriv_t1 / h_sq) + ( deriv_t2 / h_sq) ;
  const double_t w2_ = - (3 * V_.y3_before_ / h_sq) + (3 * S_.y3_ / h_sq) 
                       - ( 2 * deriv_t1 / h_ms) - ( deriv_t2 / h_ms) ;
  const double_t w1_ = deriv_t1;
  const double_t w0_ = V_.y3_before_;

  //normal form :    x^3 + r*x^2 + s*x + t with coefficients : r, s, t
  const double_t r = w2_ / w3_;
  const double_t s = w1_ / w3_;
  const double_t t = (w0_ - P_.U_th_) / w3_;
  const double_t r_sq= r*r;

  //substitution y = x + r/3 :  y^3 + p*y + q == 0 
  const double_t p = - r_sq / 3 + s;
  const double_t q = 2 * ( r_sq * r ) / 27 - r * s / 3 + t;

  //discriminante
  const double_t D = std::pow( (p/3), 3) + std::pow( (q/2), 2);

  double_t tau1;
  double_t tau2;
  double_t tau3;
      
  if(D<0){
    const double_t roh = std::sqrt( -(p*p*p)/ 27 );
    const double_t phi = std::acos( -q/ (2*roh) );
    const double_t a = 2 * std::pow(roh, (1.0/3.0));
    tau1 = (a * std::cos( phi/3 )) - r/3;
    tau2 = (a * std::cos( phi/3 + 2* numerics::pi/3 )) - r/3;
    tau3 = (a * std::cos( phi/3 + 4* numerics::pi/3 )) - r/3;
  }
  else{
    const double_t sgnq = (q >= 0 ? 1 : -1);
    const double_t u = -sgnq * std::pow(std::fabs(q)/2.0 + std::sqrt(D), 1.0/3.0);
    const double_t v = - p/(3*u);
    tau1= (u+v) - r/3;
    if (tau1 >= 0) {
      return tau1;
    }
    else {
      return thresh_find2_(dt);
    }
  }
      
  //set tau to the smallest root above 0
      
  double tau = (tau1 >= 0) ? tau1 : 2*h_ms;
  if ((tau2 >=0) && (tau2 < tau)) tau = tau2;
  if ((tau3 >=0) && (tau3 < tau)) tau = tau3;      
  return (tau <= V_.h_ms_) ? tau : thresh_find2_(dt);
}

inline
nest::double_t nest::iaf_psc_alpha_canon_nr::membrane_potential_th_(const double_t dt) const  // f(x)
{
  double_t gamma_    	= 1/P_.c_m_ /  ( 1/P_.tau_syn_ - 1/P_.tau_m_ );
  double_t gamma_sq_ 	= 1/P_.c_m_ / ( ( 1/P_.tau_syn_ - 1/P_.tau_m_ )
                            * ( 1/P_.tau_syn_ - 1/P_.tau_m_ ) );
  double_t expm1_tau_m_ = numerics::expm1(-1.*dt/P_.tau_m_);
  double_t expm1_tau_syn_ = numerics::expm1(-1.*dt/P_.tau_syn_);
  double_t P30_ 	= -P_.tau_m_ / P_.c_m_ * expm1_tau_m_;
  double_t P31_ 	= gamma_sq_ * expm1_tau_m_ - gamma_sq_ * expm1_tau_syn_
		           - dt * gamma_ * expm1_tau_syn_ - dt * gamma_;
  double_t P32_ 	=  gamma_ * expm1_tau_m_ - gamma_ * expm1_tau_syn_;
  return P30_*(P_.I_e_+V_.y0_before_) + P31_*V_.y1_before_ + P32_*V_.y2_before_ +
                      expm1_tau_m_ * V_.y3_before_  + V_.y3_before_-P_.U_th_;
}

inline
nest::double_t nest::iaf_psc_alpha_canon_nr::membrane_potential_d1_(const double_t dt) const  //  f'(x)
{
  double_t gamma_    	= 1/P_.c_m_ /  ( 1/P_.tau_syn_ - 1/P_.tau_m_ );
  double_t gamma_sq_ 	= 1/P_.c_m_ / ( ( 1/P_.tau_syn_ - 1/P_.tau_m_ )
                            * ( 1/P_.tau_syn_ - 1/P_.tau_m_ ) );
  double_t exp_tau_m_ 	= exp(-1.*dt/P_.tau_m_);
  double_t exp_tau_syn_	= exp(-1.*dt/P_.tau_syn_);
  double_t expm1_tau_syn_ = numerics::expm1(-1.*dt/P_.tau_syn_);

  double_t P30_		= 1./ P_.c_m_ * exp_tau_m_;
  double_t P31_ 	= -1.*gamma_sq_ * 1./P_.tau_m_ * exp_tau_m_ + gamma_sq_ * 1./P_.tau_syn_ * exp_tau_syn_ - gamma_*expm1_tau_syn_
			  + dt*gamma_* 1./P_.tau_syn_ * exp_tau_syn_ - gamma_;
  double_t P32_         = gamma_*  1./P_.tau_syn_ * exp_tau_syn_ - gamma_*  1./P_.tau_m_  * exp_tau_m_;
  return  P30_*(P_.I_e_+V_.y0_before_) + P31_*V_.y1_before_ + P32_*V_.y2_before_ -  1./P_.tau_m_  * exp_tau_m_ * V_.y3_before_; 
}

// finds threshpassing via NR method
nest::double_t nest::iaf_psc_alpha_canon_nr::thresh_find4_(double_t const dt) const
{
  double_t last_root=0.;
  double_t e;

  do
  {
    double_t root=last_root-membrane_potential_th_(last_root)/membrane_potential_d1_(last_root);  // x=x_0-f(x_0)/f'(x_0)

    assert(-1e-08<=root && root<=dt+1e-08);
    e=fabs(last_root-root);
    last_root=root;
  } while(e>=1e-08);
  return last_root;
}

