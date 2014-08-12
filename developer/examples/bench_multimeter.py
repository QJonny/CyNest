# benchmark of multimeter vs voltmeter

import nest
import time

nrns = 10000
sim = 1000

nest.ResetKernel()
n = nest.Create('iaf_hill', nrns)
r = nest.GetStatus(n)[0]['recordables']
m = nest.Create('multimeter')
nest.SetStatus(m, {'record_from': [r['V_m']], 
                   'to_file': False, 'withtime': False, 'withgid': False})
nest.ConvergentConnect(m, n)

mm_start = time.time()
nest.Simulate(sim)
mm_stop = time.time()


nest.ResetKernel()
n = nest.Create('iaf_hill', nrns)
r = nest.GetStatus(n)[0]['recordables']
v = nest.Create('voltmeter')
nest.SetStatus(v, {'to_file': False, 'withtime': False, 'withgid': False})
nest.ConvergentConnect(v, n)

vm_start = time.time()
nest.Simulate(sim)
vm_stop = time.time()

print 'Multimeter: %8.3f' % (mm_stop-mm_start)
print 'Voltmeter: %8.3f' % (vm_stop-vm_start)
