#! /usr/bin/env python

import matplotlib
# matplotlib.use("macosx")
import pylab

import cynest as nest
import cynest.voltage_trace

weight=20.0
delay=1.0
stim=1000.0

neuron1 = nest.Create("iaf_neuron")
neuron2 = nest.Create("iaf_neuron")
voltmeter = nest.Create("voltmeter")

nest.SetStatus(neuron1, {"I_e": stim})
nest.Connect(neuron1,neuron2,weight,delay)
nest.Connect(voltmeter, neuron2)

nest.Simulate(100.0)

nest.voltage_trace.from_device(voltmeter)
nest.voltage_trace.show()
