# This is a simple test script which gererates independant spike train
# input to two spike detectors with the the inh_gamma_generator,
# Once via two generators and two connections
# Second via 1 generator and two _S (selective connections)

# The spike trains for the two strategies are then checked to be identical.

# Author: Eilif Muller

# See also: bench_ng.py, and bench.py for a comparison of the performence of
# the two approaches

# On my Operton system, bench_ng.py is more than 10 times faster than bench.py
# (These are equivalent simulations in terms of output).

import nest
import numpy
import time

bg_e = 6.0

ige_Params = {}
ige_Params['tbins'] = numpy.array([0.0])
ige_Params['a'] = numpy.array([1.0])
ige_Params['b'] = 1.0/numpy.array([bg_e])

nest.ResetKernel()

g = nest.Create('inh_gamma_generator',2)
nest.SetStatus(g,ige_Params)

sd = nest.Create('spike_detector',2)
nest.Connect(g,sd)

print "Simulating inh_gamma_generator control." 
t1 = time.time()
nest.Simulate(10*1000)
t2 = time.time()
print "Elapsed: ", t2-t1, " seconds."

spikes1 = nest.GetStatus(sd,'events')[0]['times']

nest.ResetKernel()

# Now we create only one inh_gamma_generator
# Connecting it with _S (selective connectors)
# Allows it to send different realizations to its targets
g = nest.Create('inh_gamma_generator',1)
nest.SetStatus(g,ige_Params)

sd = nest.Create('spike_detector',2)
nest.Connect(g*2,sd,model='static_synapse_S')

print "Simulating inh_gamma_generator_ng test." 
t1 = time.time()
nest.Simulate(10*1000)
t2 = time.time()
print "Elapsed: ", t2-t1, " seconds."

spikes2 = nest.GetStatus(sd,'events')[0]['times']


if (spikes1[0] == spikes2[0]).all() and (spikes1[1] == spikes2[1]).all():
    print "Spikes for multiple targets agree (inh_gamma_generator_ng == inh_gamma_generator)"
else:
    print "ERROR: spikes do not agree (inh_gamma_generator_ng != inh_gamma_generator)"
    
