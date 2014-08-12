#! /usr/bin/env python

# A script to test triplet STDP (stdp_triplet_synapse)

# Presently, both potentiation and depression are additive

# Author: Eilif Muller

# try:
# $ ipython -pylab
# >>> execfile('test_triplet_stdp.py')
# produces nice plot, hopefully!

import numpy
import numpy.random as nr

import time

dt = 0.05
tsim = 10.0
results = []
stdp_params = {}

import nest

import numpy

import math

# avoid dependence on pylab
def load(filename):
    
    fo = file(filename,'r')
    
    # predicted size of file
    
    len_inc = 100000
    length = len_inc
    
    line = fo.readline()
    
    dim = len(line.split())
    
    # preallocate memory
    d = numpy.zeros((dim,length),float)
    
    # go to start of file
    
    fo.seek(0)
    i = 0
    
    for line in fo:
        if i>=length:
            # we need to make the arrays larger
            # add size increment
            new_d = numpy.zeros((dim,length+len_inc),float)
            new_d[:,:length] = d[:,:]
            length = length+len_inc
            d = new_d
            
        d[:,i] = numpy.array(map(float,line.split()))
        i+=1
    
    # set array size to number of elements read
    
    fo.close()
    
    return d[:,:i]


def setup_stdp():

    global eN, eP, eP1, cD, sD, stdp_params
    
    global dt

    global espikes

    nest.ResetKernel()
    
    nest.SetKernelStatus({'resolution':dt, 'overwrite_files': True})
    
    [eN]= nest.Create('iaf_cond_exp_sfa_rr',1)
        
    # set neuron model params
    
    eParams = {'V_th':-57.0, 'V_reset': -70.0, 't_ref': 0.5, 'g_L':28.95,
    'C_m':289.5, 'E_L' : -70.0, 'E_ex': 0.0, 'E_in': -75.0, 'tau_syn_ex':1.5,
    'tau_syn_in': 10.0, 'E_sfa': -70.0, 'q_sfa':14.48,
    'tau_sfa':110.0, 'E_rr': -70.0, 'q_rr': 3214.0,
    'tau_rr': 1.97,'tau_minus':33.7,'tau_minus_triplet':110.0}
    
    nest.SetStatus([eN],eParams)
    

    # spike generator input


    [eP,eP1] = nest.Create('spike_generator',2)
    
    nest.SetStatus([eP],
                   {'spike_times':numpy.array([50.0])})

    # static synapse to stimulate neuron to spike

    nest.Connect([eP],[eN],100.0,dt, model='static_synapse')
                         
    # hippocampus

    ee_synParams = {'tau_plus':16.8,'tau_x':946.0,
                    'A_2p':6.1e-3,'A_3p':6.7e-3,
                    'A_2m':1.6e-3/2.0,'A_3m':1.4e-3/2.0}

    neuron_stdpParams = {'tau_minus':33.7,'tau_minus_triplet':125.0}


    # cortex

    #ee_synParams = {'tau_plus':16.8,'tau_x':100.0,
    #                'A_2p':0.0,'A_3p':6.2e-3,
    #                'A_2m':7.0e-3,'A_3m':2.3e-4}

    #nest.SetStatus([eN],{'tau_minus':33.7,'tau_minus_triplet':27.0})

    ee_synParams['weight'] = 2.0
    ee_synParams['delay'] = dt
    
    nest.SetDefaults('stdp_triplet_synapse_S',ee_synParams)

    
    nest.SetStatus([eN],neuron_stdpParams)

    # save STDP params in a global dictionary for the analytical calculations

    stdp_params.update(ee_synParams)
    stdp_params.update(neuron_stdpParams)

    # finally, make the stdp_triplet connection

    nest.Connect([eP1],[eN], model='stdp_triplet_synapse_S')

    # measure voltage in the neuron

    [cD] = nest.Create('voltmeter',1)
    nest.SetStatus([cD],[{"to_file": True, "withtime": True}])
    nest.Connect([cD],[eN], model='stdp_triplet_synapse_S')

    # spike detector

    [sD] = nest.Create('spike_detector',1)
#    nest.SetStatus([sD],{'record_events':True})

    nest.Connect([eN],[sD], model='static_synapse')


def dw_theo_quad(dt,tau):
    """ analytical calculation of the weight change induced by
    additive triplet STDP by a quad of spikes: two pre-post pairs,
    each seperated by dt, pre-to-pre seperation of tau

    globals:
    stdp_params - parameters of the triplet STDP rule

    returns the weight change
    
    """

    global stdp_params

    A_2m = stdp_params['A_2m']
    A_3m = stdp_params['A_3m']
    A_2p = stdp_params['A_2p']
    A_3p = stdp_params['A_3p']

    tau_p = stdp_params['tau_plus']
    tau_x = stdp_params['tau_x']

    tau_m = stdp_params['tau_minus']
    tau_y = stdp_params['tau_minus_triplet']


    if dt>0:

        dw_p1 = numpy.exp(-dt/tau_p)*A_2p
        dw_p2 = (numpy.exp(-tau/tau_p)+1.0)*numpy.exp(-dt/tau_p)
        dw_p2 *= A_2p + A_3p*numpy.exp(-tau/tau_y)

        dw_m1 = numpy.exp((-tau+dt)/tau_m)*(A_2m + A_3m*numpy.exp(-tau/tau_x))
        return dw_p1+dw_p2-dw_m1
    elif dt<0:
        dw_m1 = numpy.exp(+dt/tau_m)*A_2m
        dw_p1 = numpy.exp((-tau-dt)/tau_p)
        dw_p1 *= A_2p + A_3p*numpy.exp(-tau/tau_y)

        dw_m2 = (numpy.exp(-tau/tau_m)+1.0)*numpy.exp(+dt/tau_m)
        dw_m2 *= A_2m + A_3m*numpy.exp(-tau/tau_x)

        return dw_p1-dw_m1-dw_m2

    else:
        return 0.0
     


def test_quad():
    ''' a simple simulation to test triplet stdp
    Two pre-post pairs (i.e. a quad) are simulated with NEST and computed analystically using dw_theo_quad()
    to compare the NEST results to the analytical results.

    local variables:
    dt_prepost - seperation between pre and post
    tau - seperation between quads (pre to pre)

    results are in:
    dw_theo, dw_prepost, and binned by tau.

    '''

    global eP1
    global dw_prepost,dw_theo,tau

    rho = numpy.array([2.0,10.0,20.0,30.0,40.0,50.0,60.0,80.0])
    tau = 1.0/rho*1000.0

    dw_prepost = numpy.zeros_like(tau)
    dw_theo = numpy.zeros_like(tau)
    
    dt_prepost = 10.0

    # pre post protocol

    print 'pre->post'

    for i in xrange(len(tau)):

        print i

        setup_stdp()

        dw_theo[i] = dw_theo_quad(dt_prepost,tau[i])

        # pre synaptic spike time is caluclated (50.0 + 0.82 (onset delay of postsynaptic spike) - dt_prepost)
        # last spike at approx 3000.0 is there to trigger STDP calculation
        # as NEST computes STDP only at pre-synaptic spike arrivals.
        nest.SetStatus([eP1],
                       {'spike_times':numpy.array([50.82-dt_prepost,50.82+tau[i]-dt_prepost,3000.0+tau[i]])})

        # postsynaptic spikes (induced by stong input weight at eP)
        nest.SetStatus([eP],
                       {'spike_times':numpy.array([50.0,50.0+tau[i]])})

        nest.Simulate(3020.0+tau[i])

        # get voltage

        #fname = nest.GetStatus([cD],"filenames")[0][0]
        #data = load(fname)
    
        #v = data[1,:]
        #t = data[0,:]

        # get spikes

        [de] = [nest.GetStatus([sD],'events')[0]['times']]
        # convert from timesteps to ms
        espikes = de.astype(float)*dt
        del de

        # determine exact dt_prepost
        #dt_prepost[i] = espikes[0]-50.82+prepost_dts[i]-dt
        conn = nest.FindConnections([eP1], synapse_type='stdp_triplet_synapse_S')
        dw_prepost[i] = nest.GetStatus(conn,'weight')[0] -2.0

    




def test_stdp():
    ''' a simple simulation to test stdp
     - produces a plot of weight change as a function of t_post-t_pre
     i.e. the standard STDP plot.'''

    global eP1
    global dw_prepost,dt_prepost, dw_postpre, dt_postpre

    prepost_dts = arange(1.0,30.0)

    dw_prepost = numpy.zeros_like(prepost_dts)
    dt_prepost = numpy.zeros_like(prepost_dts)

    # pre post protocol

    print 'pre->post'

    for i in xrange(len(prepost_dts)):

        print i

        setup_stdp()


        # second spike needs to be there that effect due to
        # post synaptic spike is caluclated
        nest.SetStatus([eP1],
                       {'spike_times':numpy.array([50.82-prepost_dts[i],300.0])})
                         


        nest.Simulate(320.0)

        # get voltage

        fname = nest.GetStatus([cD],"filenames")[0][0]
        data = load(fname)
    
        v = data[1,:]
        t = data[0,:]

        # get spikes
    
        [de] = [nest.GetStatus([sD],'events')[0]['times']]
        # convert from timesteps to ms
        espikes = de.astype(float)*dt
        del de

        # determine exact dt_prepost
        dt_prepost[i] = espikes[0]-50.82+prepost_dts[i]-dt

        conn = nest.FindConnections([eP1], synapse_type='stdp_triplet_synapse_S')
        dw_prepost[i] = nest.GetStatus(conn,'weight')[0] -2.0


    
    postpre_dts = arange(1.0,30.0)

    dw_postpre = numpy.zeros_like(postpre_dts)
    dt_postpre = numpy.zeros_like(postpre_dts)

    # pre post protocol

    print 'post->pre'

    for i in xrange(len(postpre_dts)):

        print i

        setup_stdp()

        # no second spike needed as for pre->post

        nest.SetStatus([eP1],
                       {'spike_times':numpy.array([50.82+postpre_dts[i]])})
                         

        nest.Simulate(100.0)

        # get voltage

        fname = nest.GetStatus([cD],"filenames")[0][0]
        data = load(fname)
    
        v = data[1,:]
        t = data[0,:]

        # get spikes
        [de] = [nest.GetStatus([sD],'events')[0]['times']]
        # convert from timesteps to ms
        espikes = de.astype(float)*dt
        del de

        # determine exact dt_prepost
        dt_postpre[i] = -(postpre_dts[i]+dt - (espikes[0]-50.82))

        conn = nest.FindConnections([eP1], synapse_type='stdp_triplet_synapse_S')
        dw_prepost[i] = nest.GetStatus(conn,'weight')[0] -2.0



    figure()

    plot(dt_postpre,dw_postpre,'b',lw=2)
    plot(dt_prepost,dw_prepost,'r',lw=2)
    plot([-30.0,30.0],[0.0,0.0],'g--',lw=2)
    xticks(size=16)
    yticks(size=16)
    xlabel(r'$t_{\rm{post}}-t_{\rm{pre}}\ \rm{\left[ms\right]}$',size=20)
    ylabel(r'$\Delta w\  \left[nS\right]$',size=20)
    
    

if __name__=='__main__':

    import pylab

    test_quad()
    pylab.figure()
    pylab.plot(tau, dw_prepost,'b',lw=2)
    pylab.plot(tau, dw_theo,'r--',lw=2)
    pylab.legend(('NEST2','Analytical'))
    pylab.title('Comparison triplet STDP of quads of spikes \n (pairs pre-post dt=10.0ms, seperated by $\\tau$)')
    pylab.xlabel(r'$\tau$')
    pylab.ylabel('dw [nS]')
    pylab.show()



