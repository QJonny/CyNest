import nest
import numpy as np
import pylab as pl

nest.ResetKernel()

nest.SetKernelStatus({'overwrite_files': True})

n = nest.Create('hill_neuron_hT')
ih_defs = nest.GetStatus(n)[0]
r = ih_defs['recordables']
s = ih_defs['receptor_types']

# cannot use dict as params = in Create()
m = nest.Create('multimeter')
nest.SetStatus(m, {'to_file': True, 'interval': 0.1,
                   'record_from': ['V_m', 'I_NaP', 'Theta', 'I_KNa', 'I_T', 'I_h']})

g = nest.Create('spike_generator',
                params = {'spike_times': 
                          sorted(np.ravel([np.linspace(10,30,num=15), np.linspace(50,70,num=15),
                                    np.linspace(15,25,num=15),
                                    np.linspace(120,140,num=15)]))})

nest.Connect(g, n, params = {'weight': 10.0, 'delay': 1.0, 'receptor_type': s['AMPA']},
                             model= 'hill_synapse')
nest.Connect(g, n, params = {'weight': 10.0, 'delay': 1.5, 'receptor_type': s['GABA_A']})
nest.Connect(m, n)

nest.Simulate(1000)

events = nest.GetStatus(m)[0]['events']
t = events['times'];

pl.clf()
Vlines = pl.plot(t, events['V_m'], 'b', t, events['Theta'], 'g')
vax = pl.gca()
nax = pl.gca().twinx()
Ilines = nax.plot(t, events['I_NaP'], 'r', t, events['I_KNa'], 'maroon', t, events['I_T'], 'orange',
                  t, events['I_h'], 'purple')
pl.legend(Vlines+Ilines, ['V_m', 'Theta', 'I_NaP', 'I_KNa', 'I_T', 'I_h'])
vax.set_xlabel('Time [ms]')
vax.set_ylabel('Membrane potential [mV]')
nax.set_ylabel('Current [pA]')
pl.title('hill_neuron_hT: Crafted AMPA + delayed GABA_A input')
#pl.legend(('V_m',  'I_NaP'))
#pl.show()

