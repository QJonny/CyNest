# Simulate iaf_cond_alpha using multimeter, voltmeter and conductancemeter 
# for comparative measurement

import nest
import numpy as np
import pylab as pl

nest.ResetKernel()

nest.SetKernelStatus({'overwrite_files': True})

n = nest.Create('iaf_cond_alpha')

# cannot use dict as params = in Create()
m = nest.Create('multimeter')
nest.SetStatus(m, {'interval': 0.5, 'start': 0.0, 'stop': 45.0,
                   'record_from': ['V_m', 'g_ex', 'g_in']})

v = nest.Create('voltmeter', params = {'interval': 0.5, 'withtime': True, 'start': 0.0, 'stop': 45.0})
c = nest.Create('conductancemeter', params = {'interval': 0.5,'withtime': True, 'start': 0.0, 'stop': 45.0})

gex = nest.Create('spike_generator',
                  params = {'spike_times': np.array([10.0, 20.0, 50.0])})
gin = nest.Create('spike_generator',
                  params = {'spike_times': np.array([15.0, 25.0, 55.0])})

nest.Connect(gex, n, params={'weight': 200.0})
nest.Connect(gin, n, params={'weight': -10.0})
nest.Connect(m, n)
nest.Connect(v, n)
nest.Connect(c, n)

nest.Simulate(100)

me = nest.GetStatus(m)[0]['events']
mt = me['times'];

vp=nest.GetStatus(v)[0]['events']['potentials']
vt=nest.GetStatus(v)[0]['events']['times']
ce=nest.GetStatus(c)[0]['events']['exc_conductance']
ci=nest.GetStatus(c)[0]['events']['inh_conductance']
#pl.clf()
#pl.plot(t, events['V_m'], t, 0.1 * events['g_ex'] - 70, t, 0.1 * events['g_in'] - 70)
#pl.legend(('V_m', 'g_ex', 'g_in'))
#pl.show()

