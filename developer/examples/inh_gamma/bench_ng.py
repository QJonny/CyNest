
# a benchmark of the inh_gamma_generator independent multiple targets
# feature.

# Compare simulation times to the equivalent multiple generator
# objects for independent realizations (bench.py).

# Test system: Intel core2quad Q6600 @ 2.4 Ghz
#
# On my test system, the timed simulation step is almost
# 20 times faster (than the bench.py implementation)
# and 7 times faster than bench_poisson.py.
# However, for STDP and dynamic synapses, one must use a
# configuration as in bench_poisson_parrot.py, due to the callback
# problem with poisson_generators.  We are 13 times faster
# than bench_poisson_parrot.py
#
# When pooling to one Poisson process is possible, as in bench_single_poisson.py,
# one is faster by a margin of about 3 times over the approach here.
# However, pooling is NOT possible for STDP, or dynamic synapses.

# So for STDP, and dynamics synapses, the inh_gamma_generator is a big
# improvement over the presently best approach, bench_poisson_parrot.py,
# by a margin of 13 times.


# Author: Eilif Muller

# See also: test_ng.py for a verification that this feature of
# inh_gamma_generator approachs yields identical spike trains to
# the equivalent multiple generator objects for independent realizations.


import nest
import numpy

import time

nest.ResetKernel()

dt = 0.1 # milliseconds
tsim = 100 # seconds
bg_e = 6.0 # Hz
bg_i = 10.2 # Hz

nest.SetStatus([0],{'resolution':dt})

# neuron to receive input 

iaf = nest.Create('iaf_cond_exp_sfa_rr',1)

iafParams = {'V_th':-57.0, 'V_reset': -70.0, 't_ref': 0.5, 'g_L':28.95,
'C_m':289.5, 'E_L' : -70.0, 'E_ex': 0.0, 'E_in': -75.0, 'tau_syn_ex':1.5,
'tau_syn_in': 10.0, 'E_sfa': -70.0, 'q_sfa':14.48,
'tau_sfa':110.0, 'E_rr': -70.0, 'q_rr': 3214.0,
'tau_rr': 1.97}

nest.SetStatus(iaf,iafParams)

# excitatory input

ige = nest.Create('inh_gamma_generator',1)

# inhihibtory input

igi = nest.Create('inh_gamma_generator',1)

# set inh_gamma_generator parameters

# excitatory

ige_Params = {}
ige_Params['tbins'] = numpy.array([0.0])
ige_Params['a'] = numpy.array([1.0])
ige_Params['b'] = 1.0/numpy.array([bg_e])
nest.SetStatus(ige,ige_Params)

# inhibitory

igi_Params = {}
igi_Params['tbins'] = numpy.array([0.0])
igi_Params['a'] = numpy.array([1.0])
igi_Params['b'] = 1.0/numpy.array([bg_i])
nest.SetStatus(igi,igi_Params)
    
# connect

# note list multiplication : [0]*4 = [0,0,0,0]
# thus we have 1000 ex, 250 inh independent connections to iaf  

nest.ConvergentConnect(ige*1000,iaf,[2.0],[dt],model='static_synapse_S')
nest.ConvergentConnect(igi*250,iaf,[-2.0],[dt],model='static_synapse_S')

# record spikes

spike_detector = nest.Create('spike_detector',1)

nest.Connect(iaf,spike_detector,model='static_synapse_S')

# simulate

t1 = time.time()
nest.Simulate(tsim*1000)
t2 = time.time()
print "Elapsed: ", t2-t1, " seconds."

# get iaf spike output

spikes = nest.GetStatus(spike_detector,'events')[0]['times']
espikes = spikes.astype(float)

espikes = espikes.compress(espikes>500.0)

erate = float(len(espikes))/(tsim-0.5)

print "IAF firing rate = ",erate, ' Hz.'






