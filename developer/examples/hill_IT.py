"""
Current stimulation protocol to test intrinsic currents 
with current stimulation protocol in Hill model.
"""

import nest
import numpy as np
import pylab as pl

nest.ResetKernel()

nest.SetKernelStatus({'overwrite_files': True})

# turn off other intrinsic currents to avoid confusion
n = nest.Create('hill_neuron_hT', params={'KNa_g_peak': 0.0, 'NaP_g_peak': 0.0, 'T_g_peak': 0.0,
                                          'h_g_peak': 20.0})
ih_defs = nest.GetStatus(n)[0]
r = ih_defs['recordables']
s = ih_defs['receptor_types']

# cannot use dict as params = in Create()
m = nest.Create('multimeter')
nest.SetStatus(m, {'to_file': True, 'interval': 0.1,
                   'record_from': ['V_m', 'I_T', 'I_NaP', 'I_KNa', 'I_h', 'I_stim']})

g = nest.Create('dc_generator')

nest.Connect(g, n, params = {'weight': 1.0, 'delay': 1.0})
nest.Connect(m, n)

for t in [50, 100, 200, 400, 800]:
    nest.SetStatus(g, {'amplitude': 20.})
    nest.Simulate(200)
    nest.SetStatus(g, {'amplitude': 10.})
    nest.Simulate(t)

events = nest.GetStatus(m)[0]['events']
t = events['times'];

pl.clf()
Vlines=pl.plot(t, events['V_m'], 'b')
vax = pl.gca()
nax = vax.twinx()
Ilines=nax.plot(t, events['I_h'], 'purple', t, events['I_stim'], 'aqua') #, t, events['I_KNa'], 'maroon', t, events['I_T'], 'orange')
vax.set_xlabel('Time [ms]')
vax.set_ylabel('Membrane potential [mV]')
nax.set_ylabel('Current [pA]')
pl.legend(Vlines+Ilines, ['V_m', 'I_h', 'I_stim'])
pl.title('Hill-Tononi hill_neuron_hT, Current stimulation protocol')

#pl.legend(('V_m',  'I_NaP'))
#pl.show()

