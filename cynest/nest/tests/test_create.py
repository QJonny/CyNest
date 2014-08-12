#! /usr/bin/env python
#
# test_create.py
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
Creation tests
"""

import unittest
import cynest
import sys

if "sample_neuron" not in cynest.Models():
    cynest.RegisterNeuron("sample_neuron")

class CreateTestCase(unittest.TestCase):


    def test_ModelCreate(self):
        """Model Creation"""       

        cynest.ResetKernel()
        node = cynest.Create("sample_neuron")


    def test_ModelCreateN(self):
        """Model Creation with N"""       

        cynest.ResetKernel()
        node = cynest.Create("sample_neuron",10)


    def test_ModelCreateNdict(self):
        """Model Creation with N and dict"""       

        cynest.ResetKernel()
        node = cynest.Create("sample_neuron",10, {})


    def test_ModelDict(self):
        """sample_neuron Creation with N and dict"""       

        cynest.ResetKernel()

        n = cynest.Create("sample_neuron", 10, [{'V_m':12.0}])
        V_m = [12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0]     
        self.assertEqual(cynest.GetStatus(n,'V_m'), V_m)
        self.assertEqual([key['V_m'] for key in cynest.GetStatus(n)], V_m)

        
    def test_ModelDicts(self):
        """sample_neuron Creation with N and dicts"""       

        cynest.ResetKernel()

        V_m = [0.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.]
        n   = cynest.Create("sample_neuron", 10, [{'V_m':v} for v in V_m])
        self.assertEqual(cynest.GetStatus(n,'V_m'), V_m)
        self.assertEqual([key['V_m'] for key in cynest.GetStatus(n)], V_m)


    def test_CopyModel(self):
        """CopyModel"""

        cynest.ResetKernel()        
        cynest.CopyModel("sample_neuron",'new_neuron',{'V_m':10.0})
        vm = cynest.GetDefaults('new_neuron')['V_m']
        self.assertEqual(vm,10.0)
        
        n = cynest.Create('new_neuron',10)
        vm = cynest.GetStatus([n[0]])[0]['V_m']
        self.assertEqual(vm,10.0)
        
        cynest.CopyModel('static_synapse', 'new_synapse',{'weight':10.})
        cynest.Connect([n[0]],[n[1]],model='new_synapse')
        w = cynest.GetDefaults('new_synapse')['weight']
        self.assertEqual(w,10.0)
        
        try:
            cynest.CopyModel("sample_neuron",'new_neuron') # shouldn't be possible a second time
            self.fail('an error should have risen!')  # should not be reached
        except: 
            info = sys.exc_info()[1]
            if not "NewModelNameExists" in info.__str__():
                self.fail('could not pass error message to cynest!')    

        
def suite():
    suite = unittest.makeSuite(CreateTestCase,'test')
    return suite


def run():
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
