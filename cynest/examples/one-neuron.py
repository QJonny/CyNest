#!/usr/bin/env python

import cynest as nest
import cynest.voltage_trace

nest.ResetKernel()

neuron = nest.Create("iaf_neuron")
nest.SetStatus(neuron, "I_e", 376.0)

voltmeter = nest.Create("voltmeter")
nest.SetStatus(voltmeter, {"withgid": True, "withtime": True})

nest.Connect(voltmeter, neuron)

nest.Simulate(1000.0)

nest.voltage_trace.from_device(voltmeter)
nest.voltage_trace.show()
