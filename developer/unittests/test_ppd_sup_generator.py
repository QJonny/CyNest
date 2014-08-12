#! /usr/bin/env python

# A script to test the sanity of the output of the ppd_sup_generator.
#
# Author: Moritz Deger
#
# try:
# $ ipython -pylab
# >>> run test_ppd_sup_generator
# produces nice plot, hopefully!

import nest
import numpy
import pylab



### parameters ###
params = {
  'amplitude': 0.1,
  'dead_time': 30.0,
  'frequency': 1.0,
  'n_proc': 1000,
  'rate': 10.0
  }
T = 10000. # [ms]
dt = 0.1 # [ms]

### info
prob = params['rate'] * dt / 1000.0
print 'transition probability is p=', prob
print 'number of processes is n=', params['n_proc']
occ0 = params['rate'] * params['dead_time']/1000.0 * params['n_proc']
print 'zero occupation in eq. is', occ0
print 'act. occupation in eq. is', params['n_proc']-occ0
print 'product with transition prob is',  prob * (params['n_proc']-occ0)


### simulation ###
nest.ResetKernel()
nest.SetStatus([0], {'resolution':dt})
gg = nest.Create('ppd_sup_generator')
sd = nest.Create('spike_detector',2)
nest.DivergentConnect(gg, sd)
nest.SetStatus(gg, params)
nest.SetStatus(sd, {'withgid': False})
nest.Simulate(T)


### analysis ###
def power_spectrum(st, om, T):
    j = numpy.complex(0,1)
    T = st[-1] - st[0]
    st_T = st - T/2.0
    contribs = numpy.exp(-j*om*st)
    return numpy.abs(numpy.sum(contribs))**2/T


spikes1 = nest.GetStatus(sd)[0]['events']['times']/1000.0
spikes2 = nest.GetStatus(sd)[1]['events']['times']/1000.0
dbin = T/1000.0/100.0
bins = numpy.arange(0., T/1000.0, dbin)
hist1 = numpy.histogram(spikes1, bins=bins)[0] / dbin / params['n_proc']
hist2 = numpy.histogram(spikes2, bins=bins)[0] / dbin / params['n_proc']
rate1 = len(spikes1) / (T/1000.0) / params['n_proc']
rate2 = len(spikes2) / (T/1000.0) / params['n_proc']
freq = numpy.arange(0.0, 100.0, 0.5)
ps1 = numpy.array([power_spectrum(spikes1, 2.*numpy.pi*f, T/1000.0) for f in freq])
ps2 = numpy.array([power_spectrum(spikes2, 2.*numpy.pi*f, T/1000.0) for f in freq])



### show some results
print 'spike rates [Hz]:', rate1, rate2
pylab.figure(1)
pylab.clf()
pylab.hist(numpy.diff(spikes1), bins=50, histtype='step', normed=True)
pylab.hist(numpy.diff(spikes2), bins=50, histtype='step', normed=True)
pylab.title('ISI histogram')
pylab.xlabel('time [s]')
pylab.ylabel('relative frequency')

pylab.figure(2)
pylab.clf()
pylab.plot(freq, ps1)
pylab.plot(freq, ps2)
pylab.title('power spectrum estimate')
pylab.xlabel('frequency [Hz]')
pylab.ylabel('power')

pylab.figure(3)
pylab.clf()
pylab.plot(bins[:-1], hist1)
pylab.plot(bins[:-1], hist2)
pylab.title('rate over time')
pylab.xlabel('time [s]')
pylab.ylabel('average rate [Hz]')



pylab.show()


