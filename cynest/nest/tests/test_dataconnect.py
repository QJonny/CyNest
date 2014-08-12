#! /usr/bin/env python
#
# test_create.py
#
# This file is part of cycynest.
#
# Copyright (C) 2004 The cycynest Initiative
#
# cycynest is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# cycynest is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cycynest.  If not, see <http://www.gnu.org/licenses/>.
"""
Creation tests
"""

import unittest
import cynest
import sys
import time

class DataConnectTestCase(unittest.TestCase):

    def test_DataConnect1(self):
        """DataConnect"""

        cynest.ResetKernel()
        
        a=cynest.Create("iaf_neuron", 100000)
        sources=[1]
        target=[1.0 * x for x in list(range(2,100000))]
        weight=[2.0 * x for x in target]
        delay=[1.0 * x for x in target]
        connections=[{'target':target, 'weight':weight, 'delay':delay} for t in target ]
        cynest.DataConnect(sources,connections, "static_synapse", 1)
        conn1=cynest.GetConnections(sources)
        stat1=cynest.GetStatus(conn1)
        target1=[ d['target'] for d in stat1]
        self.assertEqual(target, target1)

    def test_DataConnect2(self):
        """DataConnect"""

        cynest.ResetKernel()
        
        a=cynest.Create("iaf_neuron", 100000)
        sources=[1]
        target=[x for x in list(range(2,100000))]
        weight=[2.0 * x for x in target]
        delay=[1.0 * x for x in target]
        connections=[{'target':target, 'weight':weight, 'delay':delay} for t in target ]
        cynest.DataConnect(sources,connections, "static_synapse", 2)
        conn1=cynest.GetConnections(sources)
        stat1=cynest.GetStatus(conn1)
        target1=[ d['target'] for d in stat1]
        self.assertEqual(target, target1)

    def test_ConvergentConnect(self):
        """ConvergentConnect"""

        cynest.ResetKernel()
        
        a=cynest.Create("iaf_neuron", 10)
        target=[1]
        sources=[1 * x for x in list(range(2,11))]

        cynest.ConvergentConnect(sources, target, [1.0], [1.0])
        conn1=cynest.GetConnections(sources)
        stat1=cynest.GetStatus(conn1)
        target1=[ d['source'] for d in stat1]
        self.assertEqual(sources, target1)

    def test_DivergentConnect(self):
        """DivergentConnect"""

        cynest.ResetKernel()
        
        a=cynest.Create("iaf_neuron", 10)
        source=[1]
        targets=[1 * x for x in list(range(2,11))]

        cynest.DivergentConnect(source, targets, [1.0], [1.0])
        conn1=cynest.GetConnections(source)
        stat1=cynest.GetStatus(conn1)
        target1=[ d['target'] for d in stat1]
        self.assertEqual(targets, target1)

    def test_RandomDivergentConnect(self):
        """RandomDivergentConnect"""

        cynest.ResetKernel()
        
        a=cynest.Create("iaf_neuron", 1000)
        source=[1]
        targets=range(2,1000)

        cynest.RandomDivergentConnect(source,targets, 500, [1.0], [1.0])
        conn1=cynest.GetConnections(source)
        self.assertEqual(len(conn1), 500)

    def test_RandomConvergentConnect(self):
        """RandomConvergentConnect"""

        cynest.ResetKernel()
        
        a=cynest.Create("iaf_neuron", 1000)
        sources=list(range(2,1000))
        target=[1]

        cynest.RandomConvergentConnect(sources,target, 500, [1.0], [1.0])
        conn1=cynest.GetConnections(sources)
        self.assertEqual(len(conn1), 500)

def suite():
    suite = unittest.makeSuite(DataConnectTestCase,'test')
    return suite


def run():
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())


