# -*- coding: utf-8 -*-
#
# brunel2000_classes.py
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

import cynest
import cynest.raster_plot
import pylab
import os
import sys
from subprocess import call

cynest.ResetKernel()

if "cython_iaf_psc_delta_c_members" not in cynest.Models():
    cynest.RegisterNeuron("cython_iaf_psc_delta_c_members")
    
if "testmodel" not in cynest.Models():
    cynest.RegisterNeuron("testmodel")

class Brunel2000:
    """
    Implementation of the sparsely connected random network, described
    by Brunel (2000) J. Comp. Neurosci.  Parameters are chosen for the
    asynchronous irregular state (AI).
    """
    g     = 5.0    # Ratio of IPSP to EPSP amplitude: J_I/J_E
    eta   = 2.0    # rate of external population in multiples of threshold rate
    delay = 1.5    # synaptic delay in ms
    tau_m = 20.0   # Membrane time constant in mV
    V_th  = 20.0   # Spike threshold in mV
    N_E   = 800 # to change
    N_I   = 200 # to change
    J_E   = 0.1
    threads=1
    N_rec   = 50    # Number of neurons to record from
    built=False
    connected=False

    def __init__(self):
        """
        Initialize an object of this class.
        """
        self.name=self.__class__.__name__
        self.data_path=self.name+"/"
        if not os.path.exists(self.data_path):
            os.makedirs(self.data_path)
        print ("Writing data to: "+self.data_path)
        cynest.ResetKernel()
        cynest.SetKernelStatus({"data_path": self.data_path})

    def calibrate(self):
        """
        Compute all parameter dependent variables of the
        model.
        """
        self.N_neurons = self.N_E+self.N_I
        self.C_E    = self.N_E//10
        self.C_I    = self.N_I//10
        self.J_I    = -self.g*self.J_E
        self.nu_ex  = self.eta* self.V_th/(self.J_E*self.C_E*self.tau_m)
        self.p_rate = 1000.0*self.nu_ex*self.C_E
        cynest.SetKernelStatus({"print_time": True,
                              "local_num_threads":self.threads})


    def build(self, neuronName, custom):
        """
        Create all nodes, used in the model.
        """
        if self.built==True: return
        self.calibrate()


        self.nodes   = cynest.Create(neuronName, self.N_neurons)

        self.noise=cynest.Create("poisson_generator",1,{"rate": self.p_rate})
        self.spikes=cynest.Create("spike_detector",2, 
                                [{"label": "brunel-py-ex"},
                                 {"label": "brunel-py-in"}])
        if custom == True:
            cynest.SetStatus(self.nodes, [{"optimized":True}])

        self.nodes_E= self.nodes[:self.N_E]
        self.nodes_I= self.nodes[self.N_E:]
        self.spikes_E=self.spikes[:1]
        self.spikes_I=self.spikes[1:]
        self.built=True

    def connect(self, neuronName, custom):
        """
        Connect all nodes in the model.
        """
        if self.connected: return
        if not self.built:
            self.build(neuronName, custom)

        cynest.CopyModel("static_synapse_hom_wd",
                       "excitatory",
                       {"weight":self.J_E, 
                        "delay":self.delay})
        cynest.CopyModel("static_synapse_hom_wd",
                       "inhibitory",
                       {"weight":self.J_I, 
                        "delay":self.delay})

        cynest.RandomConvergentConnect(self.nodes_E, self.nodes, 
                                     self.C_E, model="excitatory")
        cynest.RandomConvergentConnect(self.nodes_I, self.nodes, 
                                     self.C_I, model="inhibitory")
        cynest.DivergentConnect(self.noise,self.nodes,model="excitatory")
        cynest.ConvergentConnect(self.nodes_E[:self.N_rec],self.spikes_E)
        cynest.ConvergentConnect(self.nodes_I[:self.N_rec],self.spikes_I)
        self.connected=True

    def run(self, neuronName, custom, simtime=40.):
        """
        Simulate the model for simtime milliseconds and print the
        firing rates of the network during htis period.  
        """
        if not self.connected:
            self.connect(neuronName, custom)
        cynest.Simulate(simtime)
        events = cynest.GetStatus(self.spikes,"n_events")

        if sys.version_info >= (3, 0):
            self.rate_ex= events[0]//simtime*1000.0/self.N_rec
            self.rate_in= events[1]//simtime*1000.0/self.N_rec
        else:
            self.rate_ex= events[0]/simtime*1000.0/self.N_rec
            self.rate_in= events[1]/simtime*1000.0/self.N_rec

        print ("Excitatory rate   : %.2f Hz" % self.rate_ex)
        print ("Inhibitory rate   : %.2f Hz" % self.rate_in)
        cynest.raster_plot.from_device(self.spikes_E, hist=True)
        pylab.show()



def runNeurons(ms, opt = True):
    print ("Running native, SLI and cython neurons for " + str(ms) + " ms")

    print ("\n\nRunning native neurons")
    # native neuron
    b = Brunel2000()
    b.run("iaf_psc_delta", False, ms)
    NativRTF = cynest.GetKernelStatus()["realtime factor"]

    cynest.ResetKernel()

    print ("\n\nRunning cython neurons")
    # cython neuron
    b = Brunel2000()
    
    b.run("cython_iaf_psc_delta_c_members", opt, ms)
    #b.run("testmodel", opt, ms)

        
    CythonRTF = cynest.GetKernelStatus()["realtime factor"]

    print ("Faster factor (native / cython) : " + str(NativRTF / CythonRTF))



def start(v, t=40):
    runNeurons(t, v)

def run():
    start(True)
