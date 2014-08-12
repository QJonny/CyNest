"""
Bombard Hill model with poisson input
"""

import nest
import numpy as np
import pylab as pl

nest.ResetKernel()

nest.SetKernelStatus({'overwrite_files': True})

p = nest.Create('parrot_neuron', 4)
nhT = nest.Create('hill_neuron_hT')
n = nest.Create('hill_neuron')

ih_defs = nest.GetStatus(n)[0]
r = ih_defs['recordables']
syns = ih_defs['receptor_types']

# cannot use dict as params = in Create()
mhT = nest.Create('multimeter')
nest.SetStatus(mhT, {'to_file': False, 'interval': 0.1,
                     'record_from': ['V_m', 'I_NaP', 'I_KNa', 'I_T', 'I_h']})
m = nest.Create('multimeter')
nest.SetStatus(m, {'to_file': False, 'interval': 0.1,
                     'record_from': ['V_m', 'I_NaP', 'I_KNa']})

g = nest.Create('poisson_generator', 4, params={'rate': 100.})

for s in zip(g, p, [('AMPA', 500.0), ('NMDA', 50.), ('GABA_A', 250.), ('GABA_B', 100.0)]):
    nest.Connect([s[0]], [s[1]])
    nest.Connect([s[1]], n, params = {'weight': s[2][1], 'receptor_type': syns[s[2][0]]},
                 model= 'hill_synapse')
    nest.Connect([s[1]], nhT, params = {'weight': s[2][1], 'receptor_type': syns[s[2][0]]},
                 model= 'hill_synapse')

nest.Connect(m, n)
nest.Connect(mhT, nhT)

nest.Simulate(1000)

events = nest.GetStatus(m)[0]['events']
eventshT = nest.GetStatus(mhT)[0]['events']
t = events['times'];

pl.clf()
vax = pl.subplot(211)
Vlines = pl.plot(t, events['V_m'], 'b')
nax = pl.gca().twinx()
Ilines = nax.plot(t, events['I_NaP'], 'r', t, events['I_KNa'], 'yellow')
pl.title('hill_neuron: Poisson stimulation via all synpases')
pl.legend(Vlines+Ilines, ['V_m', 'I_NaP', 'I_KNa'])
vax.set_ylabel('Membrane potential [mV]')
nax.set_ylabel('Current [pA]')


vax = pl.subplot(212)
Vlines = pl.plot(t, eventshT['V_m'], 'b')
nax = pl.gca().twinx()
Ilines = nax.plot(t, eventshT['I_NaP'], 'r', t, eventshT['I_KNa'], 'yellow', t, eventshT['I_T'], 'pink',
         t, eventshT['I_h'], 'black')
pl.title('hill_neuron_hT: Poisson stimulation via all synpases')
pl.legend(Vlines+Ilines, ['V_m', 'I_NaP', 'I_KNa', 'I_T', 'I_h'])
vax.set_xlabel('Time [ms]')
vax.set_ylabel('Membrane potential [mV]')
nax.set_ylabel('Current [pA]')
