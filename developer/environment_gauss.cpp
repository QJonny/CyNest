/*
 *  environment_gauss.cpp
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

#include "environment_gauss.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
//#include "fstream.h"
#include <iostream> //Larissa 7.7.09 braucht man für cout
#include <cmath> //Larissa 7.7.09
#include <vector> //Larissa 27.7.09 ??
#include "network.h"


#include "fdstream.h"

namespace nest
{

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters
 * ---------------------------------------------------------------- */

  nest::environment_gauss::Parameters_::Parameters_()
    : amp_           (0.0), // pA
      sigma_         (1.0), // width of Gaussian 
      I_inhib_       (0.0), // pA
      nps_           (0.0),
      agent_step_    (1),
      x_dim_         (0),
      y_dim_         (0),
      border_dim_    (0),
      refractor_     (0.), // ms
      input_neurons_ (0),
      neuron_west_   (0),
      neuron_east_   (0),
      neuron_south_  (0),
      neuron_north_  (0),
      critic_neurons_(0),
      specialx_      (0),
      specialy_      (0),
      warpx_         (-1),
      warpy_         (-1),
      special_reward_(0),
      seed_(0),
      close_on_reset_(true),
      filename_(),
      label_(),
      flush_after_simulate_(true),
      close_after_simulate_(false)
  {}

 /* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

 void nest::environment_gauss::Parameters_::get(DictionaryDatum &d) const
{ 
  def<double>(d, names::amplitude, amp_);
  def<double>(d, "sigma", sigma_);
  def<double>(d, "I_inhib", I_inhib_);
  def<long_t>(d, "nps", nps_);
  def<double>(d, "ref_period", refractor_);
  def<int_t> (d, "agent_step", agent_step_);
  def<int_t> (d, "x_dim", x_dim_);
  def<int_t> (d, "y_dim", y_dim_);
  def<int_t> (d, "border_dim",  border_dim_);
  def<long_t> (d, "special_x", specialx_);
  def<long_t> (d, "special_y", specialy_);
  def<double_t> (d, "special_reward", special_reward_);
  def<long_t> (d, "warpx", warpx_);
  def<long_t> (d, "warpy", warpy_);
  def<long_t> (d, "neuron_west", neuron_west_);
  def<long_t> (d, "neuron_east", neuron_east_);
  def<long_t> (d, "neuron_north", neuron_north_);
  def<long_t> (d, "neuron_south", neuron_south_);
  def<long_t> (d, "seed", seed_);
  (*d)[names::close_on_reset] = close_on_reset_;
  if (!filename_.empty() )
    (*d)[names::filename] = filename_;
  (*d)[names::label] = label_;
  (*d)[names::close_after_simulate] = close_after_simulate_;
  (*d)[names::flush_after_simulate] = flush_after_simulate_;
}


void nest::environment_gauss::Parameters_::set(const DictionaryDatum& d)
{
  updateValue<double>(d, names::amplitude, amp_);
  updateValue<double>(d, "sigma", sigma_);
  updateValue<double>(d, "I_inhib", I_inhib_);
  updateValue<long_t>(d, "nps", nps_);
  updateValue<double>(d, "ref_period", refractor_);
  updateValue<long_t>(d, "agent_step", agent_step_);
  updateValue<long_t> (d, "x_dim", x_dim_);
  updateValue<long_t> (d, "y_dim", y_dim_);
  updateValue<long_t> (d, "border_dim",  border_dim_);
  updateValue<long_t> (d, "special_x", specialx_);
  updateValue<long_t> (d, "special_y", specialy_);
  updateValue<double_t> (d, "special_reward", special_reward_);
  updateValue<long_t> (d, "seed", seed_);
  updateValue<long_t> (d, "warpx", warpx_);
  updateValue<long_t> (d, "warpy", warpy_);
  updateValue<bool>(d, names::close_on_reset, close_on_reset_);
  updateValue<std::string>(d, names::label, label_);
  updateValue<bool>(d, names::close_after_simulate, close_after_simulate_);
  updateValue<bool>(d, names::flush_after_simulate, flush_after_simulate_);
  updateValue<long_t>(d, "neuron_west", neuron_west_);
  updateValue<long_t>(d, "neuron_east", neuron_east_);
  updateValue<long_t>(d, "neuron_south", neuron_south_);
  updateValue<long_t>(d, "neuron_north", neuron_north_);
  
  std::srand(seed_);
  
  Name input_neurons_name("input_neurons");

  if(d->known(input_neurons_name))
    {
      const Token& t = d->lookup(input_neurons_name);
      ArrayDatum *ad = dynamic_cast<ArrayDatum*>(t.datum());
      if(ad != 0)
	{
	  ad->toVector(input_neurons_);
	}
      //else
      //throw TypeMismatch();
    }

  Name critic_neurons_name("critic_neurons");

  if(d->known(critic_neurons_name))
    {
      const Token& t = d->lookup(critic_neurons_name);
      ArrayDatum *ad = dynamic_cast<ArrayDatum*>(t.datum());
      if(ad != 0)
	{
	  ad->toVector(critic_neurons_);
	}
      //else
      //throw TypeMismatch();
    }

}

/* ---------------------------------------------------------------- 
 * Default and copy constructor for environment_gauss
 * ---------------------------------------------------------------- */

  nest::environment_gauss::environment_gauss()
    : Archiving_Node(),// StimulatingDevice<CurrentEvent>(), 
      P_(),
      node_(*this)
  {}

  nest::environment_gauss::environment_gauss(const environment_gauss& n)
    : Archiving_Node(n), 
      P_(n.P_),
      node_(n)    
  {}

/* ---------------------------------------------------------------- 
 * environment_gauss initialization functions
 * ---------------------------------------------------------------- */

void nest::environment_gauss::init_state_(const Node&)
{
}

void nest::environment_gauss::init_buffers_()
{
  B_.spikes_.clear(); // includes resize
  // Archiving_Node::clear_history();
  // we only close files here, opening is left to calibrate()
   if ( P_.close_on_reset_ && B_.fs_.is_open() )
   {
     B_.fs_.close();
     P_.filename_.clear();  // filename_ only visible while file open
   }
   B_.ref_period_ = 0;
   B_.h = Time::get_resolution().get_ms(); //time step
   B_.agentx_ = P_.border_dim_; //40
   B_.agenty_= P_.border_dim_; //30
   B_.actx_ = 0; //only initialization, will be changed before first trial
   B_.acty_ = 0;
   B_.state_ = pos(B_.actx_,B_.acty_,(P_.x_dim_+P_.border_dim_*2)); //current state
   random_position();
   std::cout << "agentx/agenty: " << B_.agentx_ << "/" << B_.agenty_ << std::endl;
   B_.Gausscurr_ = 0.;
   B_.sumGausscurr_ = 0.;
   B_.effamp_ = 0.;
   B_.sumGauss_ = 0.;
   B_.grid_mat_.resize(P_.y_dim_+2*P_.border_dim_);
   for (int i = 0; i < P_.y_dim_+2*P_.border_dim_; i++) {
		// set the number of columns and initialize
		// all elements to 0;
		B_.grid_mat_[i].resize(P_.x_dim_+2*P_.border_dim_); }
}

void nest::environment_gauss::calibrate()
{
  // do we need to (re-)open the file
  bool newfile = false;
     

  if ( !B_.fs_.is_open() )
    {
      newfile = true;   // no file from before
      P_.filename_ = build_filename_();
    }
  else
    {
      std::string newname = build_filename_();
      if ( newname != P_.filename_ )
      	{
	  //Node::network()->info("RecordingDevice::calibrate()",
	  //	"Closing file " + P_.filename_ +
	  //	", opening file " + newname);
	  
	  B_.fs_.close(); // close old file
	  P_.filename_ = newname;
	  newfile = true;
	}
    }

  
  if(newfile)
    {
      assert(!B_.fs_.is_open());
      if(Node::network()->overwrite_files())
	B_.fs_.open(P_.filename_.c_str());
      else
	{// try opening for reading
         std::ifstream test(P_.filename_.c_str());
         if ( test.good() )
         {
           Node::network()->message(SLIInterpreter::M_ERROR, "RecordingDevice::calibrate()",
             "The device file " + P_.filename_ + " exists already and will not be overwritten.\n"
             "Please change data_path, data_prefix or label, or set overwrite_files." );
           throw IOError();
         }
         else 
           test.close();
           
         // file does not exist, so we can open
         B_.fs_.open(P_.filename_.c_str());
       }
     }

    if ( !B_.fs_.good() )
     {
       Node::network()->message(SLIInterpreter::M_ERROR, "RecordingDevice::calibrate()", 
                              "I/O error while opening file " + P_.filename_);
       if ( B_.fs_.is_open() )
         B_.fs_.close();                        
       P_.filename_.clear();
       throw IOError();
     }
  
    

    //V_.h=Time::get_resolution().get_ms();
  //V_.new_file_ = false;

  //if(V_.new_file_){
  //open_file_();
    //  V_.new_file_ = false;
}


void nest::environment_gauss::finalize()
 {
   if ( B_.fs_.is_open() )
   {
     if ( P_.close_after_simulate_ ) 
     {
       B_.fs_.close();
       return;
     }
     
     if ( P_.flush_after_simulate_ )
       B_.fs_.flush();
       
     if ( !B_.fs_.good() )
     {
       Node::network()->message(SLIInterpreter::M_ERROR, "RecordingDevice::finalize()", 
                              "I/O error while writing to file " + P_.filename_);
       throw IOError();
     }
   }
 }



void nest::environment_gauss::set_status(const DictionaryDatum &d)
{
  Parameters_ ptmp = P_;  // temporary copy in case of errors
  ptmp.set(d);                       // throws if BadProperty
  
  P_ = ptmp;
  std::srand(99);
}

    
  


const std::string nest::environment_gauss::build_filename_() const
{

  //could be improved!
  //at the moment a new file is created for each rank
  // however the data is the same on each rank, so we would only need one file
//   const int vpdigits = static_cast<int>(std::floor(std::log10(static_cast<float>(Communicator::n_vps))) + 1);
  std::ostringstream basename;
  const std::string& path = Node::network()->get_data_path();
  if (!path.empty())
    basename << path << '/';
  if (!P_.label_.empty())
    basename << P_.label_;
  else 
    basename << "steps";

  basename <<"-"<<Communicator::rank;
  return basename.str() + '.' + "gdf";
}
 

void nest::environment_gauss::update(Time const &, const long_t from, const long_t to)
{
//   long_t start = origin.get_steps();

  
     
  for ( long_t offs = from ; offs < to ; ++offs )
    {
  
      //reads out current action
      if(B_.ref_period_>0){B_.ref_period_-=B_.h;}
 
      //reward decays exponentially
      /*
      if(reward_!=0)
      {
        reward_=special_reward_*std::exp(-(refractor_-ref_period_)/(refractor_/10.0));
      }
      */

      //reward given for delay time
      if(B_.reward_ != 0)
	{
	  if(B_.ref_period_<800.0)
	    B_.reward_=0;
	}

      long_t current_spike=(long_t) B_.spikes_.get_value(offs);
	
      //if(current_spike != 0)
      //  B_.ref_period_=P_.refractor_-B_.h; //in ms
      // then decide direction without else

      if ((P_.neuron_east_==current_spike)||(P_.neuron_west_==current_spike)||(P_.neuron_north_==current_spike)||(P_.neuron_south_==current_spike))
	 {
            int_t direction = 0;

	    B_.ref_period_=P_.refractor_-B_.h; //in ms
	    if (P_.neuron_east_==current_spike)
	      direction=0;
	    else if (P_.neuron_west_==current_spike)
	      direction=1;
	    else if(P_.neuron_north_==current_spike)
	      direction=2;
	    else if(P_.neuron_south_==current_spike)
	      direction=3;

	    move(direction);
	    plot_agent();
	 
	 }

      DSCurrentEvent ce; 
     //stimulate state, critic and action neurons (amplitude is set in event_hook)
      network()->send(*this, ce, offs); 
    }
  
}





void nest::environment_gauss::event_hook(DSCurrentEvent& e)
{
//-------------find absolute amplitude of Gauss curve so that overall input is "amp"-------------------------------------
  if(B_.effamp_ == 0.){ 	
    for(int i = 1; i <= P_.input_neurons_.size()/P_.nps_; i++){
	B_.actx_ =  i%(P_.x_dim_+P_.border_dim_*2);
	B_.acty_ =  int(i/(P_.x_dim_+P_.border_dim_*2));
	B_.sumGauss_ += exp(-(pow((B_.actx_-round((P_.x_dim_+P_.border_dim_*2)/2.)),2.)+pow((B_.acty_-round((P_.y_dim_+P_.border_dim_*2)/2.)),2.))/(2.*pow(P_.sigma_,2.))); 
	/*if(i%P_.nps_ == 0){
		//std::cout << "x: " << B_.actx_ << " y: " << B_.acty_ << " Gauss: " << exp(-(pow((B_.actx_-round((P_.x_dim_+P_.border_dim_*2)/2.)),2.)+pow((B_.acty_-round((P_.y_dim_+P_.border_dim_*2)/2.)),2.))/(2.*pow(P_.sigma_,2.)))  << std::endl;
	}*/
    } 
  std::cout << " " << std::endl;
  B_.effamp_ = P_.amp_/B_.sumGauss_;
  std::cout << "Gauss Summe: " << B_.sumGauss_ << "	 eff. Amp: " << B_.effamp_ << std::endl;
  }  
//-------------Gaussian distribution of inputs from environment to state neurons-----------------------------------------
  //std::cout << "Number of neuron: " << static_cast<long_t>(e.get_receiver().get_gid()) << std::endl;
  //std::cout << "size of P_.input_neurons_: " << P_.input_neurons_.size() << std::endl;
  if((static_cast<long_t>(e.get_receiver().get_gid()) >= P_.input_neurons_[0])&&
     (static_cast<long_t>(e.get_receiver().get_gid()) < (P_.input_neurons_.size()+P_.input_neurons_[0])))
  {  
      B_.actx_ =  ((static_cast<long_t>(e.get_receiver().get_gid())-P_.input_neurons_[0])/P_.nps_)%(P_.x_dim_+P_.border_dim_*2);
      B_.acty_ =  int((((static_cast<long_t>(e.get_receiver().get_gid()))-P_.input_neurons_[0])/P_.nps_)/(P_.x_dim_+P_.border_dim_*2));
      B_.Gausscurr_ = B_.effamp_*exp(-(pow((B_.actx_-B_.agentx_),2.)+pow((B_.acty_-B_.agenty_),2.))/(2.*pow(P_.sigma_,2.))); 
      e.set_current(B_.Gausscurr_);
      if((network()->get_slice_origin().get_steps() == 0) && (((static_cast<long_t>(e.get_receiver().get_gid())-P_.input_neurons_[0])%P_.nps_) == 0)){
	//std::cout << "State: " << B_.actx_ << "/" << B_.acty_ << " " << pos(B_.actx_,B_.acty_,(P_.x_dim_+P_.border_dim_*2)) << "	 Neuron: " << static_cast<long_t>(e.get_receiver().get_gid())-P_.input_neurons_[0] << "  Gauss_amp: " << B_.Gausscurr_ << std::endl;
        //std::cout << B_.acty_+P_.border_dim_ << " " << B_.actx_+P_.border_dim_ << "\n";
        B_.sumGausscurr_ += B_.Gausscurr_;
        B_.grid_mat_[B_.acty_][B_.actx_] = B_.Gausscurr_;
        if((B_.actx_ == P_.x_dim_+2*P_.border_dim_-1) && (B_.acty_ == P_.y_dim_+2*P_.border_dim_-1)){
	   //std::ofstream f("/home/larissa/runPython/test_fullrunGauss/grid_input.txt",std::ios::app);
	   std::ofstream f("/home/larissa/runPython/TestInari2/grid_input.txt",std::ios::trunc);	      	      
	   std::cout << "B_sumGausscurrent " << B_.sumGausscurr_ << std::endl;
           for(int i=0; i<(P_.y_dim_+2*P_.border_dim_); i++){
	       	for(int j=0; j<(P_.x_dim_+2*P_.border_dim_); j++){ 
   		  //std::cout << std::fixed << std::setprecision(1);
           	  //std::cout << std::setw(3) << "  " << B_.grid_mat_[i][j];
		  
    		  f<< std::fixed << std::setprecision(1);
            	  f<< std::setw(3) << "  " << B_.grid_mat_[i][j];
                  }
		//std::cout << "\n";
		f<< "\n";
      	   }
	f.close();
	}
      }
  }
  else e.set_current(0.0);

  //inhibit all output neurons for ref_period
   if (B_.ref_period_ > 0)
    {
      if (static_cast<long_t>(e.get_receiver().get_gid())==P_.neuron_west_)
	e.set_current(P_.I_inhib_);
      if (static_cast<long_t>(e.get_receiver().get_gid())==P_.neuron_east_)
	e.set_current(P_.I_inhib_);
      if (static_cast<long_t>(e.get_receiver().get_gid())==P_.neuron_north_)
	e.set_current(P_.I_inhib_);
      if (static_cast<long_t>(e.get_receiver().get_gid())==P_.neuron_south_)
	e.set_current(P_.I_inhib_);
                 
    }
 

  //stimulate critic neurons if there is a reward
  for (index i=0; i<P_.critic_neurons_.size(); i++)
    {
      if(static_cast<long_t>(e.get_receiver().get_gid())==P_.critic_neurons_[i])
	e.set_current(B_.reward_);
    }

  e.get_receiver().handle(e);
}





void nest::environment_gauss::handle(SpikeEvent & e)
{
  
   if(e.get_multiplicity()!=0)
    {
     
      B_.spikes_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin())-e.get_delay()+1,
			(double_t) e.get_sender().get_gid());
     
    }
   else
    B_.spikes_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin())-e.get_delay()+1,
			(double_t) e.get_multiplicity());

  // -e.get_delay()+1 ensures that the spike delivered here can be read out without delay time
  
  
}

void nest::environment_gauss::move(int_t direction){
  B_.old_agent_[0]=B_.agentx_;
  B_.old_agent_[1]=B_.agenty_;
  if(P_.y_dim_>1){
    switch(direction){
    case 0 : go_east(); break;
    case 1 : go_west(); break;
    case 2 : go_north(); break;
    case 3 : go_south(); break;
      //default : std::cout<< "environment_gauss->Move-Error: Ungueltige Bewegung or two spikes at same time!! direction:"<< direction<<std::endl;
    }
  }
  else{ // only 1 dimension
    switch(direction){
    case 0 : go_east(); break;
    case 1 : go_west(); break;
      // default: std::cout<< "environment_gauss->Move-Error: Ungueltige Bewegung or two spikes at same time!! direction:"<< direction<<std::endl;
    }
  }
  //set_reward();
  if(P_.special_reward_!=0) special_moves();
  B_.prev_state_=B_.state_;
  B_.state_=pos(B_.agentx_,B_.agenty_,(P_.x_dim_+P_.border_dim_*2));
  
}

void nest::environment_gauss::go_south(){
  B_.agenty_ += P_.agent_step_;
  if(pos_valid()) B_.wall_=false;
  else {
    B_.agenty_-= P_.agent_step_;
    B_.wall_=true;
  }
}

void nest::environment_gauss::go_north(){
  B_.agenty_ -= P_.agent_step_;
  if(pos_valid()) B_.wall_=false;
  else {
    B_.agenty_ += P_.agent_step_;
    B_.wall_=true;
  }
}

void nest::environment_gauss::go_west(){
  B_.agentx_ -= P_.agent_step_;
  if(pos_valid()) B_.wall_=false;
  else {
    B_.agentx_ += P_.agent_step_;
    B_.wall_=true;
  }
}

void nest::environment_gauss::go_east(){
  B_.agentx_ += P_.agent_step_;
  if(pos_valid()) B_.wall_=false;
  else{
    B_.agentx_ -= P_.agent_step_;
    B_.wall_=true;
  }
}


bool nest::environment_gauss::pos_valid(){
  if (B_.agentx_ < P_.border_dim_ || B_.agentx_ >= P_.x_dim_+P_.border_dim_) return false;
  if(P_.y_dim_>1){
    if(B_.agenty_ < P_.border_dim_ || B_.agenty_ >= P_.y_dim_+P_.border_dim_) return false;
  }
  return true; 
}


int nest::environment_gauss::pos(int_t x, int_t y, int_t columns){
  return y*columns+x;
}

void nest::environment_gauss::special_moves(){
  //warp, if ran into special
  if(abs(B_.old_agent_[0]-P_.specialx_) < 5 && abs(B_.old_agent_[1]-P_.specialy_) < 5){
//  if(B_.old_agent_[1]==P_.specialy_){
//    if(B_.old_agent_[0]==P_.specialx_){
      B_.reward_=0;
      //random warp
      if(P_.warpx_<0 || P_.warpy_<0)
	{
	  random_position();
	}
      else { //warp to specific position
	put_agent(P_.warpx_,P_.warpy_);
      }   
//    }
  }

  if(sqrt(pow((B_.agentx_-P_.specialx_),2.)+pow((B_.agenty_-P_.specialy_),2.)) < 4.5){
//  if(B_.agenty_==P_.specialy_){
//    if(B_.agentx_==P_.specialx_){
      B_.reward_=P_.special_reward_;
    }
//  }
}

void nest::environment_gauss::random_position(){
  double_t xpos, ypos;
  // before with C++ rng, calling only rand()
  // CHANGEME!!! check RAND_MAX with NEST rng
  // librandom::RngPtr rng = net_->get_rng(get_thread());
  // xpos = ((double) P_.x_dim_*rng->rand())/(RAND_MAX+1.0);
  xpos = ((double) P_.x_dim_ *rand())/(RAND_MAX+1.0) + P_.border_dim_;
  ypos = ((double) P_.y_dim_ *rand())/(RAND_MAX+1.0) + P_.border_dim_; 
  std::cout<<"random pos"<<"\t"<<(int)xpos<<"\t"<<(int)ypos<<std::endl;
  put_agent((int)xpos,(int)ypos);
}

void nest::environment_gauss::put_agent(int x, int y){
  B_.agentx_=x;
  if(P_.y_dim_>1) B_.agenty_=y;

  if(pos_valid())
    return;
  else{
    //std::cout<<"error: agent's position does not exist!!"<<std::endl;
    //put_agent(P_.border_dim_,P_.border_dim_); Larissa: This leads to increased probability to start at 0,0 !!
    random_position(); //Larissa: Not sure if this works 
  }
  B_.prev_state_=B_.state_;
  B_.state_=pos(x,y,(P_.x_dim_+P_.border_dim_*2));

}


void nest::environment_gauss::plot_agent(){
   // file_<< network()->get_slice_origin().get_steps() << "\t"<< B_.state_ <<"\t"<< reward_<< std::endl;
  B_.fs_ << network()->get_slice_origin().get_steps()<< '\t'<< B_.state_ << '\t'<< B_.reward_ << '\n';
  // fs_ <<  B_.state_ << '\t'<<I_ reward_ << '\n';
}
}

/* Larissa:
void nest::environment_gauss::cout_grid_mat(int rows, int cols){
   for (int i = 0; i < rows; i++) {		
	for (int j = 0; j < cols; j++) {
		std::cout << fixed << setprecision(6);
		std::cout << setw(12) << "  " << getPoolrate(j);
	}
	cout << "\n";
}
*/
  /*
void nest::environment_gauss::open_file_()
{
  std::ostringstream basename;
  basename << "steps";
  basename << "-" << Communicator::rank << ".gdf";
  //std::string filename_;
  filename_ = basename.str();
  // fs_ = ostreamPtr(new ofdstream(filename_.c_str()));

 
}
  */
