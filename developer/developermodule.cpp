/*
 *  developermodule.cpp
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

/*
    This file is part of NEST

    developermodule.cpp -- module providing unreleased developer modules

    Author(s):
    Hans Ekkehard Plesser,
    based on work by M.O.Gewaltig, E. Mueller and M. Helias

    First Version: June 2006
*/

/*
  Symbol to be accessible from outside the library through libltdl.
  Its name has to be composed like <modulename>_LTX_<symbolname>. It is then
  accessible by dlopening modulename and seaching for the symbol <symbolname>.
*/

#include <vector>
#include "config.h"
#include "developermodule.h"
#include "network.h"
#include "model.h"
#include "genericmodel.h"
#include <string>
#include "generic_connector_model.h"
#include "generic_connector.h"
#include "common_synapse_properties.h"
#include "booldatum.h"
#include "integerdatum.h"
#include "tokenarray.h"
#include "exceptions.h"
#include "sliexceptions.h"
#include "network.h"

#include "selective_connector.h"
#include "static_connection.h"
#include "static_connection_hom_wd.h"
#include "annealing_connection.h"

//
// The following models are commented out for a reason! Please see the mail at:
// http://ken.brainworks.uni-freiburg.de/cgi-bin/mailman/private/nest_developer/2008-July/001765.html
//

// Node headers
//#include "ac_poisson_generator.h"

#include "iaf_psc_alpha_mod.h"
#include "iaf_psc_alpha_norec.h"
#include "iaf_psc_alpha_dynth.h"
#include "theta_neuron.h"
#include "theta_neuron_ps.h"
#include "iaf_chs_2007.h"

#include "iaf_psc_delta_canon_cvv.h"
#include "iaf_psc_delta_nodelay.h"
#include "pif_psc_delta_canon_cvv.h"
#include "tsodyks2_stdp_doublet_connection.h"
// #include "iaf_psc_delta_canon_stepcurr.h"
// #include "iaf_psc_alpha_canon_nr.h"
// #include "iaf_psc_exp_canon.h"

#include "spike_dilutor.h"
#include "individual_spike_generator.h"


#ifdef HAVE_GSL
#include "inh_gamma_generator.h"
#include "ac_gamma_generator.h"
#include "a2eif_cond_exp.h"
#include "a2eif_cond_exp_HW.h"
#include "random_datums.h"
#include "iaf_chxk_2008.h"
#include "ht_neuron_fs.h"
#endif

// for RandomPopulationConnectD
#include "binomial_randomdev.h"
#include "normal_randomdev.h"
#include "fdstream.h"

// connector headers
#include "maturing_connection.h"
#include "maturing_connector.h"
#include "maturing_connection_fr.h"
#include "maturing_connector_fr.h"

#include "markram_connection.h"
#include "schemmel_s1_connection.h"

#include "stdp_connection.h"
#include "stdp_connection_hom.h"
#include "stdp_pl_connection_hom.h"

#include "stdp_triplet_connection.h"

// #include "dopamine_connection.h"
// #include "policy_connection.h"
#include "environment.h"
#include "environment_cliff.h"
// #include "dopamine_th_connection.h"
// #include "policy_th_connection.h"
// #include "iaf_cond_delta.h"
#include "static_connection_mult0.h"

#include "stdp_connection_rsnn_spikepairing_hom.h"

#include "lossy_connection.h"

namespace nest
{
  Network *DeveloperModule::net_ = 0;

  /* At the time when DeveloperModule is constructed, the SLI Interpreter
     must already be initialized. DeveloperModule relies on the presence of
     the following SLI datastructures: Name, Dictionary
  */

  DeveloperModule::DeveloperModule(Network& net)
  {
    assert(net_ == 0);
    net_ = &net;
  }

  DeveloperModule::~DeveloperModule()
  {
  }

  const std::string DeveloperModule::name(void) const
  {
    return std::string("NEST Developer Models Module"); // Return name of the module
  }

  const std::string DeveloperModule::commandstring(void) const
  {
    return std::string("/developer-init /C++ ($Revision: 7766 $) provide-component "
		       "/developer-init /SLI ($Revision: 7766 $) require-component");
  }

#ifdef HAVE_GSL
/* BeginDocumentation
   Name: RandomConnectSynfireInDegreeOutdegree - Connects a network with same indegree-outdegree as a network of synfire chains.

   Synopsis:
   [gid] w l m           RandomConnectSynfireInDegreeOutdegree -> ds dt
   [gid] w l m /synmodel RandomConnectSynfireInDegreeOutdegree -> ds dt

   Parameters:
   [gid]     - array of neurons to be wired
   w         - width of one chain
   l         - length of one chain
   m         - number of chains
   /synmodel - The synapse model for the connection (see Options below)

   Options:
   If not given, the synapse model is taken from the Options dictionary
   of the Connect command.

   Description:
   RandomConnectSynfireInDegreeOutdegree connects a network (given as
   [gid]) to have the same indegree-outdegree distributions than a
   network of randomly embedded synfire chains. Synfire chains are NOT constructed,
   only their synaptic distribution. Therefore the with and length of the
   (hypothetical) chain needs to be given. An estimate of in/out degree
   not machted is put on the stack.

   Examples:
   modeldict begin
   subnet Create /s Set
   s ChangeSubnet
   iaf_neuron 2000 Create
   [0] ChangeSubnet s GetGlobalNodes /neurons Set ;
   neurons 100 20 50 RandomConnectSynfireInDegreeOutdegree
   /dt Set
   /ds Set
   (number of neurons with not enough sources for indegree: ) =only ds =
   (number of neurons with not enough targets for indegree: ) =only dt =

   Author: Schrader, obtained from paranel
   SeeAlso: Connect
*/
  void DeveloperModule::RandomConnectSynfireInDegreeOutdegree_a_i_i_i_lFunction::execute(SLIInterpreter *i) const
  {
    i->assert_stack_load(5);

    TokenArray neuron_adr = getValue<TokenArray>(i->OStack.pick(4));

    const long w = getValue<long>(i->OStack.pick(3));
    const long L = getValue<long>(i->OStack.pick(2));
    const long M = getValue<long>(i->OStack.pick(1));
    const Name synmodel_name = getValue<std::string>(i->OStack.pick(0));

    const Token synmodel = get_network().get_synapsedict().lookup(synmodel_name);
    if ( synmodel.empty() )
      throw UnknownSynapseType(synmodel_name.toString());
    const index synmodel_id = static_cast<index>(synmodel);

    index source;
    index target;
    const size_t n = neuron_adr.size();
    const double p = static_cast<double>(w*L)/n;

    librandom::RngPtr rng = get_network().get_grng();
    librandom::BinomialRandomDev bino(rng, p,M);

    vector<long> indegree(n);
    vector<long> outdegree(n);

    const int maxiter=50000;

    for(index j=0;j<n;++j)
    {
      indegree[j]  = w* bino.uldev();
      outdegree[j] =0;
    }

    int jj=0;
    for(index to=0;to<n;++to)
    {
      int j=1;
      int c=indegree[to];

      while ( c>0 && j<maxiter)
      {
        int k=rng->ulrand(n);

        if ( outdegree[k]<indegree[k] )
        {
          outdegree[k]++;
          c--;

          source = neuron_adr[k];
          target = neuron_adr[to];

          try
          {
            get_network().connect(source, target, synmodel_id);
          }
          catch(KernelException &e)
          {
            i->message(SLIInterpreter::M_ERROR,
		       "Connect", "Operation caused an error.");
            i->raiseerror(e.what());
            return;
          }
          catch(...)
          {
            i->message(SLIInterpreter::M_ERROR,
		       "Connect", "Operation caused an error.");
            return;
          }
        }
        else
          j++;
      }
      if (j==maxiter)
        jj++;
    }

    long dt=0;
    for(size_t j=0;j<n;j++)
    {
      if ( outdegree[j]-indegree[j] < 0 )
        dt++;
    }

    i->OStack.pop(5);

    // number of neurons with not enough sources for indegree
    i->OStack.push(IntegerDatum(jj));
    // number of neurons with not enough targets for indegree
    i->OStack.push(IntegerDatum(dt));

    i->EStack.pop();
  }
#endif //HAVE_GSL


/*
BeginDocumentation

   Name: RandomPopulationConnectD - Randomly connect a population of
   source nodes to a population of target nodes using normally distributed
   parameters.

   Synopsis:
   source_array target_array N param_dict synapse_model RandomPopulationConnectD

   Parameters:
   source_array       - array containing the gids of the presynaptic populations
   target_array       - array containing the gids of the postsynaptic populations
   N                  - number of connections to be established
   param_dict         - dictionary for discribing the parameter
                        distributions and the property switches of
			the connection routine
   synapse_model      - synapse model of the connection

   Options:

   Description:
   This connection routine randomly connects the pre- and
   postsynaptic populations by randomly drawing N connections
   from the given populations. The parameters of the connection
   are given in a specialized dictionary.

   Necessary entries are:

   Optional entries are:

   /!\ Currently this high-level connection routine can cope only with
   static synapses and tsodyks synapses. /!\

  Author: Tobias Potjans, Markus Diesmann
  SeeAlso: RandomDivergentConnect, Connect
*/
   void DeveloperModule::RPopulationConnect_ia_ia_i_d_lFunction::execute(SLIInterpreter *i) const
   {
     i->assert_stack_load(5);

     TokenArray source_adr = getValue<TokenArray>(i->OStack.pick(4));
     TokenArray target_adr = getValue<TokenArray>(i->OStack.pick(3));
     ulong_t N = (ulong_t) getValue<long_t>(i->OStack.pick(2));
     DictionaryDatum param_dict = getValue<DictionaryDatum>(i->OStack.pick(1));

     const Name synmodel_name = getValue<std::string>(i->OStack.pick(0));
     const Token synmodel = get_network().get_synapsedict().lookup(synmodel_name);
     if ( synmodel.empty() )
       throw UnknownSynapseType(synmodel_name.toString());
     const index synmodel_id = static_cast<index>(synmodel);

     int_t M = Communicator::get_num_virtual_processes();
     int_t proc = Communicator::get_rank();

     long_t size_sources = source_adr.size();
     long_t size_targets = target_adr.size();
     long_t last_gid_sources = source_adr[size_sources-1];
     long_t last_gid_targets = target_adr[size_targets-1];
     long_t first_gid_sources = last_gid_sources - size_sources + 1;
     long_t first_gid_targets = last_gid_targets - size_targets + 1;

     uint_t offset_targets = first_gid_targets % M;

     // get switches
     bool static_type = param_dict->known("static_type") ? getValue<bool>(param_dict, "static_type") : false;
     bool tsodyks_type = param_dict->known("tsodyks_type") ? getValue<bool>(param_dict, "tsodyks_type") : false;
     bool distribute_weights = param_dict->known("distribute_weights") ? getValue<bool>(param_dict, "distribute_weights") : false;
     bool distribute_delays = param_dict->known("distribute_delays") ? getValue<bool>(param_dict, "distribute_delays") : false;
     bool distribute_stp_params = param_dict->known("distribute_stp_params") ? getValue<bool>(param_dict, "distribute_stp_params") : false;
     bool initialize_stp_synapse = param_dict->known("initialize_stp_synapse") ? getValue<bool>(param_dict, "initialize_stp_synapse") : false;
     bool distribute_stp_synapse_initialization = param_dict->known("distribute_stp_synapse_initialization") ? getValue<bool>(param_dict, "distribute_stp_synapse_initialization") : false;
     //bool num_targets_per_source = param_dict->known("num_targets_per_source") ? getValue<bool>(param_dict, "num_targets_per_source") : false;
     bool num_sources_per_target = param_dict->known("num_sources_per_target") ? getValue<bool>(param_dict, "num_sources_per_target") : false;
     bool num_syn_per_conn = param_dict->known("num_syn_per_conn") ? getValue<bool>(param_dict, "num_syn_per_conn") : false;
     bool exact_multinomial = param_dict->known("exact_multinomial") ? getValue<bool>(param_dict, "exact_multinomial") : true;

     if ((!static_type) && (!tsodyks_type))
     {
       i->message(SLIInterpreter::M_ERROR, "RandomPopulationConnectD", "static_type or tsodyks_type have to be set to true");
       throw BadProperty();
     }

     double_t res = Time::get_resolution().get_ms();

     // extraction of parameters from param_dict
     // 1. check for setting of all required parameters
     // 2. declaration and default initialization
     // 3. extraction

     // 1. check for required dictionary elements
     try
     {
       getValue<double_t>(param_dict, "weight_m");
       if (getValue<double_t>(param_dict, "delay_m") < res)
       {
	 i->message(SLIInterpreter::M_ERROR, "RandomPopulationConnectD", "delay_m specified to be smaller than resolution.");
	 return;
       }
       if (distribute_weights)
	 getValue<double_t>(param_dict, "weight_s");
       if (distribute_delays)
	 getValue<double_t>(param_dict, "delay_s");
     }
     catch (DictError &e)
     {
       i->message(SLIInterpreter::M_ERROR, "RandomPopulationConnectD", "One or more necessary parameters are not specified in param_dict. See documentation for list of necessary parameters.");
       i->raiseerror(e.what());
       return;
     }

     // 2. declarations and default value initialization
     double_t weight;
     double_t weight_m;
     double_t weight_s;
     bool insert_weights;
     double_t delay;
     double_t delay_m;
     double_t delay_s;
     bool insert_delays;
     double_t U;
     double_t tau_rec;
     double_t tau_fac;
     double_t tau_psc;
     double_t U_m;
     double_t tau_rec_m;
     double_t tau_fac_m;
     double_t U_s;
     double_t tau_rec_s;
     double_t tau_fac_s;
     bool insert_U;
     bool insert_tau_rec;
     bool insert_tau_fac;
     bool insert_tau_psc;
     double_t u;
     double_t x;
     double_t y;
     double_t u_m;
     double_t x_m;
     double_t y_m;
     double_t u_s;
     double_t x_s;
     double_t y_s;
     bool insert_u;
     bool insert_xy;
     //     ostreamPtr os_ntps;
     //ostreamPtr os_nspt;
     //ostreamPtr os_nspc;
     //std::ofstream os_ntps;
     std::ofstream os_nspt;
     std::ofstream os_nspc;
     std::ostringstream filename;
     std::ostringstream filename_mul;
     std::vector<ulong_t> nspt;
     std::vector<ulong_t> nspc;
     ulong_t nspc_samples;

     weight_m = 0.0;
     weight_s = 0.0;
     delay_m = res;
     delay_s = 0.0;
     U_m = 0.0;
     U_s = 0.0;
     tau_rec_m = 0.0;
     tau_rec_s = 0.0;
     tau_fac_m = 0.0;
     tau_fac_s = 0.0;
     tau_psc = 0.0;
     u_m = 0.0;
     u_s = 0.0;
     x_m = 0.0;
     x_s = 0.0;
     y_m = 0.0;
     y_s = 0.0;
     insert_U = false;
     insert_tau_rec = false;
     insert_tau_fac = false;
     insert_tau_psc = false;
     insert_u = false;
     insert_xy = false;

     // 3. extraction of parameters from dictionary param_dict
     if ( param_dict->known("weight_m") )
     {
       insert_weights = true;
       if ( distribute_weights )
       {
	 weight_m = getValue<double_t>(param_dict, "weight_m");
	 weight_s = param_dict->known("weight_s") ? getValue<double_t>(param_dict, "weight_s") : 0.0;
       }
       else
	 weight = getValue<double_t>(param_dict, "weight_m");
     }
     else
       insert_weights = false;

     if ( param_dict->known("delay_m") )
     {
       insert_delays = true;
       if ( distribute_delays )
       {
	 delay_m = getValue<double_t>(param_dict, "delay_m");
	 delay_s = param_dict->known("delay_s") ? getValue<double_t>(param_dict, "delay_s") : 0.0;
       }
       else
	 delay = getValue<double_t>(param_dict, "delay_m");
     }
     else
       insert_delays = false;

     if (tsodyks_type)
     {
       insert_U = param_dict->known("U_m");
       insert_tau_rec = param_dict->known("tau_rec_m");
       insert_tau_fac = param_dict->known("tau_fac_m");
       insert_tau_psc = param_dict->known("tau_psc");
       if ( insert_U )
       {
	 U_m = getValue<double_t>(param_dict, "U_m");
	 U_s = param_dict->known("U_s") ? getValue<double_t>(param_dict, "U_s") : 0.0;
       }
       if ( insert_tau_rec )
       {
	 tau_rec_m = getValue<double_t>(param_dict, "tau_rec_m");
	 tau_rec_s = param_dict->known("tau_rec_s") ? getValue<double_t>(param_dict, "tau_rec_s") : 0.0;
       }
       if ( insert_tau_fac )
       {
	 tau_fac_m = getValue<double_t>(param_dict, "tau_fac_m");
	 tau_fac_s = param_dict->known("tau_fac_s") ? getValue<double_t>(param_dict, "tau_fac_s") : 0.0;
       }
       if ( insert_tau_psc )
	 tau_psc = getValue<double_t>(param_dict, "tau_psc");

       if (initialize_stp_synapse)
       {
	 insert_u = param_dict->known("u_m");
	 insert_xy = (param_dict->known("x_m") && param_dict->known("y_m"));

	 if ( insert_u )
	 {
	   u_m = getValue<double_t>(param_dict, "u_m");
	   u_s = param_dict->known("u_s") ? getValue<double_t>(param_dict, "u_s") : 0.0; // 1.0?
	 }
	 if ( insert_xy )
	 {
	   x_m = getValue<double_t>(param_dict, "x_m");
	   x_s = param_dict->known("x_s") ? getValue<double_t>(param_dict, "x_s") : 0.0; // 1.0?
	   y_m = getValue<double_t>(param_dict, "y_m");
	   y_s = param_dict->known("y_s") ? getValue<double_t>(param_dict, "y_s") : 0.0;
	 }
       }
     }

     // drawing connection ids

     // unnormalized distributions of nodes on processes will be used
     // as probability distributions for drawing partitions of
     // connections for virtual processes from multinomial
     // distribution
     std::vector<ulong_t> target_distribution(size_targets, size_targets/M);
     // a generalization to create target_distribution would be a
     // function called locate_neurons that might be expanded to usage
     // of arrays as populations for now, we create the
     // target_distribution here
     int_t k;
     // correction for nodes not uniformly distributed
     // due to size_population%num_processes!=0
     for( k = offset_targets; k < (int_t) ( size_targets % M ) + (int_t) offset_targets; ++k )
     {
       target_distribution[k]++;
     }

     // first_gid_targets_distribution contains at position vp the
     // first (target) gid of all (target) gids belonging to that
     // virtual process. This will be used to map later random numbers
     // to local target gids.
     std::vector<ulong_t> first_gid_targets_distribution(M, 0);
     for( k = 0; k < M; ++k )
     {
       // fill with correct values
       first_gid_targets_distribution[ ( first_gid_targets + k ) % M ] = first_gid_targets + k;
     }

     // We use the multinomial distribution to determine the number of
     // connections that will be made on one virtual process, i.e. we
     // partition the set of edges into n_vps subsets. The number of
     // edges on one virtual process is binomially distributed with
     // the boundary condition that the sum of all edges over virtual
     // processes is the total number of processes.
     // To obtain the target_partitioning we adapt the gsl
     // implementation of the multinomial distribution.

     // K from gsl is equivalent to M = n_vps
     // N is already taken from stack
     // p[] is target_distribution
     std::vector<ulong_t> target_partitioning; // corresponds to n[]

     // calculate exact multinomial distribution
     if (exact_multinomial)
     {
       // get global rng that is tested for synchronization for all threads
       librandom::RngPtr grng = get_network().get_grng();

       // begin code adapted from gsl 1.8 //
       double_t sum_dist = 0.0; // corresponds to sum_p
       //norm is equivalent to size_targets
       uint_t sum_partitions = 0; // corresponds to sum_n

       for ( k = 0; k < M; k++ )
       {
	 if (target_distribution[k] > 0)
         {
	   // substituting gsl_ran call
	   librandom::BinomialRandomDev bino(grng, ((double_t) target_distribution[k] ) / (size_targets - sum_dist), N - sum_partitions);
	   target_partitioning.push_back(bino.uldev());
	 }
	 else
	 {
	   target_partitioning.push_back(0);
	 }

	 sum_dist += (double_t) target_distribution[k];
	 sum_partitions += (uint_t) target_partitioning[k];
       }
       // end code adapted from gsl 1.8
     }
     else
     {
       for ( k = 0; k < M; k++ )
       {
	 target_partitioning.push_back(0);
       }
     }
     // rng for local vp random numbers
     librandom::RngPtr rng;

     // loop over local threads on one machine
     int_t p;
     for (p = 0; p < M / Communicator::get_num_processes(); p++)
     {
       // current vp corresponding to proc
       int_t j;

       j = p * Communicator::get_num_processes() + proc;

       if(get_network().is_local_vp(j))
       {
	 rng = get_network().get_rng(get_network().vp_to_thread(j));
	 // normal deviate for parameter grabbing
	 librandom::NormalRandomDev norm(rng);

	 // if not exact multinomial mode, then calculate here number of neuron
	 // on this process
	 if (!exact_multinomial)
	 {
	   if (target_distribution[j] > 0)
	   {
	     librandom::BinomialRandomDev bino(rng, ((double_t) target_distribution[j] ) / size_targets, N);
	     target_partitioning[j] = bino.uldev();
	   }
	 }

	 // prepare filenames and initialize vectors for graph property readouts
	 if (num_sources_per_target)
	 {
	   filename << getValue<std::string>(param_dict, "num_sources_per_target_fstem") << "-" << Communicator::get_rank() << "-" << p << ".dat";
	   //os_nspt = ostreamPtr(new ofdstream(filename.str().c_str()));
	   os_nspt.open(filename.str().c_str());
	   if ( !os_nspt.good() )
	   {
	     Node::network()->message(SLIInterpreter::M_ERROR,
				    "RandomPopulationConnect",
				    "I/O error while opening file " + filename.str());
	     if ( os_nspt.is_open() )
	       os_nspt.close();
	     filename.str().clear();
	     throw IOError();
	   }

	   nspt = vector<ulong_t>(target_distribution[j]);
	 }
	 if (num_syn_per_conn)
	 {
	   filename_mul << getValue<std::string>(param_dict, "num_syn_per_conn_fstem") << "-" << Communicator::get_rank() << "-" << p << ".dat";
	   nspc_samples = param_dict->known("nspc_samples") ? getValue<long_t>(param_dict, "nspc_samples") : 1;
	   //os_nspc = ostreamPtr(new ofdstream(filename_mul.str().c_str()));
	   os_nspc.open(filename.str().c_str());
	   if ( !os_nspc.good() )
	   {
	     Node::network()->message(SLIInterpreter::M_ERROR,
				      "RandomPopulationConnect",
				      "I/O error while opening file " + filename.str());
	     if ( os_nspc.is_open() )
	       os_nspc.close();
	     filename.str().clear();
	     throw IOError();
	   }
	   nspc = vector<ulong_t>(nspc_samples * target_distribution[j]);
	 }
	 else
	 {
	   nspc_samples = 0;
	 }



	 while( target_partitioning[j] > 0 )
	 {
	   // draw random numbers for source node from all source neurons
	   long_t s_index  = rng->ulrand(size_sources);
	   // draw random numbers for target node from
	   // target_distribution on this virtual process
	   long_t t_index  = rng->ulrand(target_distribution[j]);

	   // use s_index and t_index to fill graph readout variables
	   if (num_sources_per_target)
	   {
	     nspt[t_index]++;
	   }
	   if (num_syn_per_conn)
	   {
	     for ( size_t nspc_i = 0; nspc_i < nspc_samples; ++nspc_i )
	     {
	       if (t_index == s_index + (long_t) nspc_i)
	       {
		 nspc[t_index + nspc_i * target_distribution[j]]++;
	       }
	     }
	   }
	   // map random number of source node to gid corresponding to
	   // first_gid_sources
	   long_t source_gid = first_gid_sources + s_index;
	   // map random number of target node to gid using the
	   // first_gid_targets_distribution corresponding to the
	   // modulo distribution of the nodes on virtual processes
	   long_t target_gid = first_gid_targets_distribution[j] + t_index * (long_t) M;

	   // local part of connection settings:
	   // - draw random number for every non-constant parameter
	   // - make all parameters available for connect via substituting
	   //   SetSynapseDefaults call

	   // define parameters specified in param_dict for this connection
	   // define connection parameter dict cpdict
	   DictionaryDatum cpdd(new Dictionary());
	   if (insert_weights)
	   {
	     if (distribute_weights)
	     {
	       do {
		 weight = norm()*weight_s + weight_m;
	       } while ((( weight_m >= 0.0 ) && ( weight < 0.0 )) || (( weight_m < 0.0 ) && ( weight > 0.0 )));
	       // weight_m == 0.0 is interpreted as exc. if weights are distributed
	     }
	     def<double_t>(cpdd, "weight", weight);
	   }
	   if (insert_delays)
	   {
	     if (distribute_delays)
	     {
	       do {
		 delay = distribute_delays ? round((norm()*delay_s + delay_m)/res)*res : round(delay_m/res)*res;
	       } while (delay < 0.0);
	       delay = (delay == 0.0) ? res : delay;
	     }
	     def<double_t>(cpdd, "delay", delay);
	   }
	   if (tsodyks_type)
	   {
	     if (distribute_stp_params)
	     {
	       if (insert_U)
	       {
		 do {
		   U = U_s > 0.0 ? norm()*U_s + U_m : U_m;
		 } while (U < 0.0);
		 def<double_t>(cpdd, "U", U);
	       }
	       if (insert_tau_rec)
	       {
		 do {
		   tau_rec = tau_rec_s > 0.0 ? norm()*tau_rec_s + tau_rec_m : tau_rec_m;
		 } while (tau_rec < 0.0);
		 def<double_t>(cpdd, "tau_rec", tau_rec);
	       }
	       if (insert_tau_fac)
	       {
		 do {
		   tau_fac = tau_fac_s > 0.0 ? norm()*tau_fac_s + tau_fac_m : tau_fac_m;
		 } while (tau_fac < 0.0);
		 def<double_t>(cpdd, "tau_fac", tau_fac);
	       }
	     }
	     else // don't distribute_stp_params
	     {
	       if (insert_U)
	       {
		 U = U_m;
		 def<double_t>(cpdd, "U", U);
	       }
	       if (insert_tau_rec)
	       {
		 tau_rec = tau_rec_m;
		 def<double_t>(cpdd, "tau_rec", tau_rec);
	       }
	       if (insert_tau_fac)
	       {
		 tau_fac = tau_fac_m;
		 def<double_t>(cpdd, "tau_fac", tau_fac);
	       }
	     }
	     if (insert_tau_psc) //not distributed for consistency
	                         //with neuron model, see warning above
	     {
	       def<double_t>(cpdd, "tau_psc", tau_psc);
	     }
	     if (initialize_stp_synapse)
	     {
	       if (distribute_stp_synapse_initialization)
	       {
		 if (insert_u)
		 {
		   u = u_s > 0.0 ? norm()*u_s + u_m : u_m;
		   do
		   {
		     u = norm()*u_s + u_m;
		   } while ( u < 0.0 || u > 1.0 );
		   def<double_t>(cpdd, "u", u);
		 }
		 if (insert_xy)
		 {
		   do
		   {
		     do {
		       x = x_s > 0.0 ? norm()*x_s + x_m : x_m;
		     } while (x < 0.0);
		     do {
		       y = y_s > 0.0 ? norm()*y_s + y_m : y_m;
		     } while (y < 0.0);
		   } while ( x + y > 1.0 );
		   def<double_t>(cpdd, "x", x);
		   def<double_t>(cpdd, "y", y);
		 }
	       }
	       else
	       {
		 if (insert_u)
		 {
		   u = u_m;
		   def<double_t>(cpdd, "u", u);
		 }
		 if (insert_xy)
		 {
		   x = x_m;
		   def<double_t>(cpdd, "x", x);
		   y = y_m;
		   def<double_t>(cpdd, "y", y);
		 }
	       }
	     }
	   }

	   // connect source and target with dictionary
	   try
	   {
	     get_network().connect(source_gid, target_gid, cpdd, synmodel_id);
	   }
	   catch (UnknownSynapseType &e)
	   {
	     i->message(SLIInterpreter::M_ERROR, "RandomPopulationConnectD", "Specified synapse type does not exist.");
	     i->raiseerror(e.what());
	     return;
	   }
	   catch(BadProperty &e)
	   {
	     i->message(SLIInterpreter::M_ERROR, "RandomPopulationConnectD", "One of the properties does not exist or has wrong type.");
	     i->raiseerror(e.what());
	     return;
	   }
	   catch(BadDelay &e)
	   {
	     i->message(SLIInterpreter::M_ERROR, "RandomPopulationConnectD", "The delay has to be greater than or equal to the resolution.");
	     i->raiseerror(e.what());
	     return;
	   }
	   catch(DimensionMismatch &e)
	   {
	     i->message(SLIInterpreter::M_ERROR, "RandomPopulationConnectD", "One of the properties was given with wrong dimensions.");
	     i->raiseerror(e.what());
	     return;
	   }

	   target_partitioning[j]--;
	 }
	 // write graph readout to file
	 if (num_sources_per_target)
	 {
	   for ( size_t nspt_i = 0; nspt_i < nspt.size(); ++nspt_i )
	   {
	     os_nspt << nspt[nspt_i] << "\n";
	   }
	 }
	 if (num_syn_per_conn)
	 {
	   std::vector<ulong_t> nspc_hist = std::vector<ulong_t>(0);
	   for ( size_t nspc_i = 0; nspc_i < nspc.size(); ++nspc_i )
	   {
	     while (nspc[nspc_i] > nspc_hist.size())
	     {
	       nspc_hist.push_back(0);
	     }
	     nspc_hist[nspc[nspc_i]]++;
	   }
	   for ( size_t nspc_hist_i = 0; nspc_hist_i < nspc_hist.size(); ++nspc_hist_i )
	   {
	     os_nspc << nspc_hist_i << "\t" << nspc_hist[nspc_hist_i] << "\n";
	     //std::cout << nspc_hist_i << "\t" << nspc_hist[nspc_hist_i] << "\n";
	   }
	 }
       }
     }
     i->OStack.pop(5);
     i->EStack.pop();
   }




  //-------------------------------------------------------------------------------------

  void DeveloperModule::init(SLIInterpreter *i)
  {
    assert(net_ != 0);

    // register models
	//    register_model<ac_poisson_generator>(*net_, "ac_poisson_generator");
    register_model<environment>(*net_, "environment");
    register_model<environment_cliff>(*net_, "environment_cliff");

    register_model<iaf_psc_alpha_mod>(*net_, "iaf_psc_alpha_mod");
    register_model<iaf_psc_alpha_norec>(*net_, "iaf_psc_alpha_norec");
    register_model<iaf_psc_alpha_dynth>(*net_, "iaf_psc_alpha_dynth");
    register_model<iaf_chs_2007>(*net_, "iaf_chs_2007");
    //register_model<iaf_psc_delta_cvv>(*net_, "iaf_psc_delta_cvv");

    register_model<iaf_psc_delta_canon_cvv>(*net_, "iaf_psc_delta_canon_cvv");
    register_model<iaf_psc_delta_nodelay>(*net_, "iaf_psc_delta_nodelay");
    register_model<pif_psc_delta_canon_cvv>(*net_, "pif_psc_delta_canon_cvv");


    // register_model<iaf_psc_delta_canon_stepcurr>(*net_, "iaf_psc_delta_canon_stepcurr");

    register_model<theta_neuron>(*net_, "theta_neuron");
    register_model<theta_neuron_ps>(*net_, "theta_neuron_ps");
    // register_model<iaf_psc_alpha_canon_nr>(*net_, "iaf_psc_alpha_canon_nr");
    // register_model<iaf_psc_exp_canon>(*net_, "iaf_psc_exp_canon");

    register_model<spike_dilutor>(*net_, "spike_dilutor");
    register_model<individual_spike_generator>(*net_, "individual_spike_generator");

    i->createcommand("RandomPopulationConnect_ia_ia_i_d_l",
		     &rpopulationconnect_ia_ia_i_d_lfunction);

#ifdef HAVE_GSL
    i->createcommand("RandomConnectSynfireInDegreeOutdegree_a_i_i_i_l", &randomconnectsynfireindegreeoutdegree_a_i_i_i_lfunction);

    // register_model<iaf_cond_delta>(*net_, "iaf_cond_delta");
    register_model<iaf_chxk_2008>(*net_, "iaf_chxk_2008");
    register_model<inh_gamma_generator>(*net_, "inh_gamma_generator");
    register_model<ac_gamma_generator>(*net_, "ac_gamma_generator");
#endif

#ifdef HAVE_GSL_1_11
    register_model<a2eif_cond_exp>(*net_, "a2eif_cond_exp");
    register_model<a2eif_cond_exp_HW>(*net_, "a2eif_cond_exp_HW");
    register_model<ht_neuron_fs>(*net_, "ht_neuron_fs");
#endif

    // register additional synapse types

    register_prototype_connection_connector_commonproperties< MaturingConnection,
                                                              MaturingConnector,
                                                              MaturingCommonProperties >
      (*net_, "maturing_synapse");



    register_prototype_connection_connector_commonproperties< MaturingConnectionFr,
                                                              MaturingConnectorFr,
                                                              MaturingCommonPropertiesFr >
      (*net_, "maturing_synapse_fr");


    // Register SelectiveConnector variants of synapse types
    net_->register_synapse_prototype(
       new GenericConnectorModel< MarkramConnection, CommonSynapseProperties,
       SelectiveConnector<MarkramConnection, CommonSynapseProperties> >(*net_, "markram_synapse_S") );

    net_->register_synapse_prototype(
       new GenericConnectorModel< SchemmelS1Connection, CommonSynapseProperties,
       SelectiveConnector<SchemmelS1Connection, CommonSynapseProperties> >(*net_, "schemmel_s1_synapse_S") );

    net_->register_synapse_prototype(
       new GenericConnectorModel< StaticConnection, CommonSynapseProperties,
       SelectiveConnector<StaticConnection, CommonSynapseProperties> >(*net_, "static_synapse_S") );

    net_->register_synapse_prototype(
       new GenericConnectorModel< STDPConnection, CommonSynapseProperties,
       SelectiveConnector<STDPConnection, CommonSynapseProperties> >(*net_, "stdp_synapse_S") );

    net_->register_synapse_prototype(
       new GenericConnectorModel< STDPTripletConnection, CommonSynapseProperties,
       SelectiveConnector<STDPTripletConnection, CommonSynapseProperties> >(*net_, "stdp_triplet_synapse_S") );



    net_->register_synapse_prototype(
       new GenericConnectorModel< STDPConnectionHom, STDPHomCommonProperties,
       SelectiveConnector<STDPConnectionHom, STDPHomCommonProperties> >(*net_, "stdp_synapse_hom_S") );

    net_->register_synapse_prototype(
       new GenericConnectorModel< STDPPLConnectionHom, STDPPLHomCommonProperties,
       SelectiveConnector<STDPPLConnectionHom, STDPPLHomCommonProperties> >(*net_, "stdp_pl_synapse_hom_S") );


    // Register non-selective variants
    register_prototype_connection<MarkramConnection>(*net_, "markram_synapse");

    register_prototype_connection<SchemmelS1Connection>(*net_, "schemmel_s1_synapse");

    register_prototype_connection<STDPTripletConnection>(*net_, "stdp_triplet_synapse");

    register_prototype_connection<Tsodyks2STDPDoubletConnection>(*net_, "tsodyks2_stdp_doublet_synapse");

    // register_prototype_connection_commonproperties<PolicyConnection, PolicyCommonProperties>(*net_, "policy_synapse");

    // register_prototype_connection_commonproperties<STDPDopaConnection, STDPDopaCommonProperties>(*net_, "stdp_dopamine_synapse");

    // register_prototype_connection_commonproperties<DOPAMINE_TH_Connection, DOPAMINE_TH_CommonProperties>(*net_, "dopamine_th_synapse");

    // register_prototype_connection_commonproperties<Policy_TH_Connection, Policy_TH_CommonProperties>(*net_, "policy_th_synapse");

    // static connection with weight, delay, rport, target
    register_prototype_connection<StaticConnectionMult0>(*net_,    "static_synapse_mult0");

    register_prototype_connection_commonproperties<STDPRSNNSpikePairingConnectionHom, STDPRSNNSpikePairingHomCommonProperties>(*net_, "stdp_rsnn_spikepairing_synapse_hom");

    register_prototype_connection<LossyConnection>(*net_, "lossy_synapse");

    register_prototype_connection_commonproperties <AnnealingConnection, AnnealingCommon > (*net_, "annealing_synapse");

  }  // DeveloperModule::init()

} // namespace nest
