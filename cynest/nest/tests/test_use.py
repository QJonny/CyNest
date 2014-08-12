#! /usr/bin/env python
#
# test_use.py
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

def testFct():
    return -1

class testClass:
    def __init__(self):
        pass
    def m1(self):
        return -2

class CreateTestCase(unittest.TestCase):


    def test_ModelCreateSimulate(self):
        """Model Creation and Simulation"""       

        cynest.ResetKernel()

        cynest.SetDefaults("sample_neuron", {"param":20})
        node = cynest.Create("sample_neuron")
        cynest.Simulate(1)
        
        self.assertEqual(cynest.GetStatus(node)[0]["param"], 30)   

    def test_Instance(self):
        """Test of instance object"""

        cynest.ResetKernel()

        t = testClass()

        node = cynest.Create("sample_neuron")

        try:
            cynest.SetStatus(node, { "param2":t})

        except:
            info = sys.exc_info()[1]


    def test_Function(self):
        """Test of function object"""

        cynest.ResetKernel()

        node = cynest.Create("sample_neuron")

        try:
            cynest.SetStatus(node, { "param2":testFct})

        except:
            info = sys.exc_info()[1]



    def test_ParamTypes(self):
        """Test of all parameter types"""    

        cynest.ResetKernel()

        node = cynest.Create("sample_neuron")
        cynest.SetStatus(node, {"param1":2, "param2":6.5, "param3": True, "param4":"sd", "param5":'r',  \
        "param6":[1,2,3,4], "param7":{"c1":1, "c2":2}, "param8":{"e1":{"s1":1}}})

        cynest.Simulate(1)

        s = cynest.GetStatus(node)[0]

        self.assertEqual(s["param1"], 2)
        self.assertEqual(s["param2"], 6.5)
        self.assertEqual(s["param3"], True)
        self.assertEqual(s["param4"], "sd")
        self.assertEqual(s["param5"], 'r')
        self.assertEqual(s["param6"][0], 1)
        self.assertEqual(s["param6"][1], 2)
        self.assertEqual(s["param6"][2], 3)
        self.assertEqual(s["param6"][3], 4)
        self.assertEqual(s["param7"]["c1"], 1)
        self.assertEqual(s["param7"]["c2"], 2)
        self.assertEqual(s["param8"]["e1"]["s1"], 1)


def suite():
    suite = unittest.makeSuite(CreateTestCase,'test')
    return suite


def run():
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
