

# Simulate the PSTH of a stepped inh_gamma_generator_ng

# The method given here is however applicable to any arbitrary modulation
# (specified by a time histogram).

# Author: Eilif Muller 

# NOTE:
# plotting - requires "ipython -pylab" interactive plotting.


import nest
import numpy
import time

bg_e = 6.0 # in Hz
dt = 0.1 # in ms
tsim = 1.5 # in s
trials = 100000 # ummm trials.




# the equilibrium firing rate function ...

tbins = numpy.array([0.0,0.5,1.0,1.5])

f = numpy.array([6.0,10.0,6.0,6.0])

# gamma order of 3.0
gamma_a = 3.0
# gamma order of 1.0
#gamma_a = 1.0

a = numpy.ones_like(tbins)*gamma_a
# need to skip first bin in f
# inh_gamma_generator time histogram for a,b
# defines a(t) = a_i on [t_i,t_{i+1})
b = 1.0/(a*f)

# tbins in ms
ig_Params = {}
ig_Params['tbins'] = tbins*1000.0
ig_Params['a'] = a
ig_Params['b'] = b
#ig_Params['stop'] = 800.0


nest.ResetKernel()

g = nest.Create('inh_gamma_generator',1)
nest.SetStatus(g,ig_Params)

sd = nest.Create('spike_detector',1)
nest.ConvergentConnect(g*trials,sd,model='static_synapse_S')

print "Simulating inh_gamma_generator test." 
t1 = time.time()
nest.Simulate(tsim*1000.0)
t2 = time.time()
print "Elapsed: ", t2-t1, " seconds."

spikes = nest.GetStatus(sd,'events')[0]['times']

# "nice" the spikes, convert to seconds.
spikes = spikes.astype(float)/1000.0

# plot PSTH

h_dt = 0.001 # in s

h_t = numpy.arange(0,tsim,h_dt) 

psth = numpy.histogram(spikes,h_t)[0]

# normalize

psth = psth.astype(float)/(h_dt*trials)

# plotting - requires "ipython -pylab" interactive plotting

figure()
plot(h_t[:-1],psth,ls='steps',lw=2)
plot(tbins,[0.0]+list(f[:-1]),'r-',ls='steps',lw=3)

legend(('gamma psth, a=%.2f' % gamma_a,'poisson (a=1.0)'),loc='lower right')

xticks(size=18)
yticks(size=18)

xlabel('time [s]',size=20)
ylabel('"firing rate" [Hz]',size=20)
