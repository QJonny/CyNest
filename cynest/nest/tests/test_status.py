#! /usr/bin/env python
#
# test_status.py
#
# This file is part of cynest.
#
# Copyright (C) 2004 The cynest Initiative
#
# cynest is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# cynest is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cynest.  If not, see <http://www.gnu.org/licenses/>.
"""
Test if Set/GetStatus work properly
"""

import unittest
import cynest
import sys

if "sample_neuron" not in cynest.Models():
    cynest.RegisterNeuron("sample_neuron")

class StatusTestCase(unittest.TestCase):
    """Tests of Set/GetStatus"""

    def test_GetDefaults(self):
        """GetDefaults"""

        model = "sample_neuron"
        
        cynest.ResetKernel()
        d = cynest.GetDefaults(model)        


    def test_SetDefaults(self):
        """SetDefaults"""

        m = "sample_neuron"
        
        cynest.ResetKernel()
        cynest.SetDefaults(m,{'V_m':-1.})
        self.assertEqual(cynest.GetDefaults(m)['V_m'], -1.)

        try:
            cynest.SetDefaults(m,{'DUMMY':0})
        except:
            info = sys.exc_info()[1]
            if not "DictError" in info.__str__():
                self.fail('wrong error message')


    def test_GetStatus(self):
        """GetStatus"""

        m = "sample_neuron"
       
        cynest.ResetKernel()
        n  = cynest.Create(m)
        d  = cynest.GetStatus(n)
        v1 = cynest.GetStatus(n)[0]['param']
        v2 = cynest.GetStatus(n,'param')[0]
        self.assertEqual(v1,v2)
        n  = cynest.Create(m,10)
        self.assertEqual(len(cynest.GetStatus(n,'param')), 10)


    def test_SetStatus(self):
        """SetStatus with dict"""

        m = "sample_neuron"
      
        cynest.ResetKernel()
        n = cynest.Create(m)
        cynest.SetStatus(n,{'V_m':1.})
        self.assertEqual(cynest.GetStatus(n,'V_m')[0], 1.)

        
    def test_SetStatusList(self):
        """SetStatus with list"""

        m = "sample_neuron"

        cynest.ResetKernel()
        n  = cynest.Create(m)
        cynest.SetStatus(n,[{'V_m':2.}])
        self.assertEqual(cynest.GetStatus(n,'V_m')[0], 2.)


    def test_SetStatusParam(self):
        """SetStatus with parameter"""

        m = "sample_neuron"
        
        cynest.ResetKernel()
        n  = cynest.Create(m)
        cynest.SetStatus(n,'V_m',3.)
        self.assertEqual(cynest.GetStatus(n,'V_m')[0], 3.)


    def test_SetStatusVth_E_L(self):
        """SetStatus of reversal and threshold potential """

        m = "sample_neuron"

        cynest.ResetKernel()

        neuron1 = cynest.Create(m)
        neuron2 = cynest.Create(m)

        # must not depend on the order.
        new_EL = -90.
        new_Vth= -10.
        cynest.SetStatus(neuron1,{'E_L': new_EL})
        cynest.SetStatus(neuron2,{'V_th': new_Vth})

        cynest.SetStatus(neuron1,{'V_th': new_Vth})
        cynest.SetStatus(neuron2,{'E_L': new_EL})

        # next three lines for debugging
        vth1, vth2 = cynest.GetStatus(neuron1,'V_th'), cynest.GetStatus(neuron2,'V_th')
        if vth1 != vth2:
            print (m, vth1, vth2, cynest.GetStatus(neuron1,'E_L'), cynest.GetStatus(neuron2,'E_L'))

        assert(cynest.GetStatus(neuron1,'V_th')==cynest.GetStatus(neuron2,'V_th'))
    



def suite():

    suite = unittest.makeSuite(StatusTestCase,'test')
    return suite


def run():
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
