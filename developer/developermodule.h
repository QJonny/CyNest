/*
 *  developermodule.h
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

#ifndef DEVELOPERMODULE_H
#define DEVELOPERMODULE_H
/* 
    This file is part of NEST

    developermodule.h -- Header to the developermodule
    (see cpp file for details)


    Author: Marc-Oliver Gewaltig (marc-oliver.gewaltig@honda-ri.de)
            Hans Ekkehard Plesser (hans.ekkehard.plesser@umb.no)

    $Date$
    Last change: $Author$
    $Revision$
*/

#include "slimodule.h"
#include "slifunction.h"
#include <cassert>

namespace nest
{
  class Network;
  class Node;

  /**
   * Module supplying models not released to the public.
   */
  class DeveloperModule: public SLIModule
  {
    public:

    DeveloperModule(Network&);
    ~DeveloperModule();

    /**
     * Initialize module by registering models with the network.
     * @param SLIInterpreter* SLI interpreterm, must know modeldict
     * @param nest::Network&  Network with which to register models
     */
    void init(SLIInterpreter*);

    const std::string name(void) const;
    const std::string commandstring(void) const;

    /**
     * Return a reference to the network
     */
    static Network &get_network();
    
    class RPopulationConnect_ia_ia_i_d_lFunction: public SLIFunction
    {
    public:
      void execute(SLIInterpreter *) const;
    } rpopulationconnect_ia_ia_i_d_lfunction;

#ifdef HAVE_GSL
    class RandomConnectSynfireInDegreeOutdegree_a_i_i_i_lFunction: public SLIFunction
    { 
     public:
      void execute(SLIInterpreter *) const;
    } randomconnectsynfireindegreeoutdegree_a_i_i_i_lfunction;
#endif //HAVE_GSL

   private:
    /** @todo Should be reference, might not need to be static. */
    static Network* net_;
  };

inline
Network& DeveloperModule::get_network()
{
  assert(net_ != 0);
  return *net_;
}

} // namespace

#endif
