#!/usr/bin/env python

# This script compares the two variants of the Tsodyks/Markram synapse in NEST.

import cynest as nest
import cynest.voltage_trace

nest.ResetKernel()

#Parameter set for depression
dep_params={"U":0.67, "weight":250.}

# parameter set for facilitation
fac_params={"U":0.1, "tau_fac":1000.,"tau_rec":100.,"weight":250.}

# Here we assign the parameter set to the synapse models 
t1_params=fac_params       # for tsodyks_synapse
t2_params=t1_params.copy() # for tsodyks2_synapse

nest.SetDefaults("tsodyks_synapse",t1_params)
nest.SetDefaults("tsodyks2_synapse",t2_params)
nest.SetDefaults("iaf_psc_exp",{"tau_syn_ex": 3.})

neuron = nest.Create("iaf_psc_exp",3)

nest.Connect([neuron[0]],[neuron[1]],model="tsodyks_synapse")
nest.Connect([neuron[0]],[neuron[2]],model="tsodyks2_synapse")
voltmeter = nest.Create("voltmeter",2)
nest.SetStatus(voltmeter, {"withgid": True, "withtime": True})
nest.Connect([voltmeter[0]], [neuron[1]])
nest.Connect([voltmeter[1]], [neuron[2]])

nest.SetStatus([neuron[0]], "I_e", 376.0)
nest.Simulate(500.0)
nest.SetStatus([neuron[0]], "I_e", 0.0)
nest.Simulate(800.0)
nest.SetStatus([neuron[0]], "I_e", 376.0)
nest.Simulate(500.0)
nest.SetStatus([neuron[0]], "I_e", 0.0)
nest.Simulate(100.0)

nest.voltage_trace.from_device([voltmeter[0]])
nest.voltage_trace.from_device([voltmeter[1]])
nest.voltage_trace.show()
