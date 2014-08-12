
# See bench_ng.py
# This is an equivalent implementation of bench_ng.py without using
# the inh_gamma_generator_ng independent multiple targets feature.

# Compare simulation times to the equivalent inh_gamma_generator_ng
# (single generator object delivers independent trains to multiple
# targets) implementation (bench_ng.py)
#
# On my Opteron system, the timed simulation step here is more than 10
# times slower (than bench_ng.py).
#

# Author: Eilif Muller

# See also: test_ng.py for a verification that the
# inh_gamma_generator_ng approachs yields identical spike trains to
# the equivalent inh_gamma_generator (multiple generator objects for
# independence) implementation.


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

ige = nest.Create('inh_gamma_generator',1000)

# inhihibtory input

igi = nest.Create('inh_gamma_generator',250)

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

nest.ConvergentConnect(ige,iaf,[2.0],[dt])
nest.ConvergentConnect(igi,iaf,[-2.0],[dt])

# record spikes

spike_detector = nest.Create('spike_detector',1)

nest.Connect(iaf,spike_detector)



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






