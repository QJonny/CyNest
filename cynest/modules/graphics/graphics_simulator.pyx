# 
# First python script to interactively simulate nest_engine.cynest with opengl in-and-output.
# 
# @authors: csaba.ero@epfl.ch, jonny.quarta@epfl.ch
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


import os
import sys
import socket
import subprocess
import random
import time

# Note that every method begins with graphics_simulator_. This
# is done because of possible name redundancy
# in other parts of the kernel (everything is global)

# sockets
graphics_simulator_sender = graphics_simulator_listener = None
# list of added spike detectors (for event recording)
graphics_simulator_spike_detector = []
# available ids (without added spike detectors)
graphics_simulator_ids = []
# simulation_step
graphics_simulator_sim_time_delta = 0.0
graphics_simulator_verbosity_level = 0
graphics_simulator_overwrite_files = False
graphics_simulator_sim_time_total = 0.0


def graphics_simulator_socket_listen(port, host = 'localhost'):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.setblocking(0);
    s.bind((host, port))
    s.listen(1)
    
    while True:
        try:
            conn, addr = s.accept()
            conn.setblocking(0);
            break
        except:
            pass

    return conn


def graphics_simulator_socket_speak(port, host = 'localhost'):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((host, port))
    return s




def graphics_simulator_receiveKeyword(keys, length = -1):
    if length == -1:
        length = len(keys[0].encode())
        
    start_time = time.time()
    while True:
        if(time.time() - start_time >= 5):
            raise NESTError("Graphics Simulator not responding. Stop.")
            
        try:
            res = graphics_simulator_listener.recv(length)
            for i, k in enumerate(keys):
                if res == k.encode():
                    return i
        except:
            pass
            
        
def graphics_simulator_receiveKeywordNonBlocking(keys, length = -1):
    if length == -1:
        length = len(keys[0].encode())
        
    try:
        res = graphics_simulator_listener.recv(length)
        for i, k in enumerate(keys):
            if res == k.encode():
                return i
    except:
        pass
    
    return -1


def graphics_simulator_send(msg):
    graphics_simulator_sender.send(msg.encode())




# Creates a new spike detector (if not already created) and 
# attaches it to the network, then connects everything to it
def graphics_simulator_initOperations():
    global graphics_simulator_spike_detector, graphics_simulator_ids, graphics_simulator_sim_time_delta, graphics_simulator_verbosity_level, graphics_simulator_overwrite_files
    
    graphics_simulator_overwrite_files = nest_engine.cynest.GetKernelStatus()["overwrite_files"]
    graphics_simulator_verbosity_level = nest_engine.cynest.get_verbosity()
    nest_engine.cynest.set_verbosity("M_WARNING")
    nest_engine.cynest.SetKernelStatus({"overwrite_files":True})
    
    size = int(nest_engine.cynest.GetKernelStatus()['network_size'])
    
    if(len(graphics_simulator_spike_detector) == 0):
        graphics_simulator_spike_detector = nest_engine.cynest.Create("spike_detector")
    
    graphics_simulator_ids = [x['global_id'] for x in nest_engine.cynest.GetStatus(list(range(1,size))) if x['model'] != "spike_detector"]

    not_connected_ids = [x for x in graphics_simulator_ids if not graphics_simulator_spike_detector in [y[1] for y in nest_engine.cynest.GetConnections([x])] ]

    if len(not_connected_ids) > 0:
        nest_engine.cynest.ConvergentConnect(not_connected_ids, graphics_simulator_spike_detector) # why does this connection have an impact on the overal simulation???
        
    if classes.get_max_delay() > 0.0:
        graphics_simulator_sim_time_delta = float(classes.get_max_delay() * 2.0)
    else:
        graphics_simulator_sim_time_delta = float(2.0)
    



def graphics_simulator_sendPositions():
    for i in graphics_simulator_ids:
        status = nest_engine.cynest.GetStatus([i])[0]

        if 'pos-x' in status and 'pos-y' in status and 'pos-z' in status:
            msg = "[" + str(i) + "," + str(status["pos-x"]) + "," + str(status["pos-y"]) + "," + str(status["pos-z"]) + "]" # [id, x,y,z]
        else:
            msg = "[" + str(i) + ",#]"
            
        graphics_simulator_send(msg)
        graphics_simulator_receiveKeyword(["ok"])
        
    graphics_simulator_send("end")
    graphics_simulator_receiveKeyword(["ok"])
    
    
def graphics_simulator_sendConnections():    
    for i in graphics_simulator_ids:
        listMsg = list(  set([x[1] for x in nest_engine.cynest.GetConnections([i])]) - set(graphics_simulator_spike_detector)  ) # [id1,id2,id3...], no multiple values
        msg = str(listMsg)

        graphics_simulator_send("[" + str(i) + "," + str(len(listMsg)) + "," + str(len(msg)) + "]") # [id, nb connections, length next msg]
        graphics_simulator_receiveKeyword(["param_ok"])
        graphics_simulator_send(msg)
        graphics_simulator_receiveKeyword(["msg_ok"])
        
    graphics_simulator_send("end")
    graphics_simulator_receiveKeyword(["ok"])








def graphics_simulator_init():
    global graphics_simulator_sender, graphics_simulator_listener
    
    print ("Initializing graphics simulator...")
    
    port = random.randint(10000, 65500)
    
    subprocess.Popen([nest_engine.exec_dir + '/g_simulator', str(port)])
    
    graphics_simulator_listener = graphics_simulator_socket_listen(port)
    graphics_simulator_sender = graphics_simulator_socket_speak(port + 1)
    
    graphics_simulator_receiveKeyword(["ready"])
    graphics_simulator_send("[" + str(graphics_simulator_sim_time_total) + "]")
    graphics_simulator_receiveKeyword(["ok"])
    
    graphics_simulator_initOperations()
    
    graphics_simulator_sendPositions()
    graphics_simulator_sendConnections()
    
    
    
def graphics_simulator_simulate(time):
    t = 0.0
    keepGoing = True
    pt = nest_engine.cynest.GetKernelStatus()["print_time"]
    nest_engine.cynest.SetKernelStatus({"print_time":False})
    
    graphics_simulator_receiveKeyword(["simulate"])
    
    print ("Simulation started...")
       
    while t < time:
        result = graphics_simulator_receiveKeywordNonBlocking(["quit", "stop", "resume"], 6)
        
        if result == 0:
            return
        elif result == 1:
            keepGoing = False
        elif result == 2:
            keepGoing = True
        
        if keepGoing:
            nest_engine.cynest.sps(graphics_simulator_sim_time_delta)
            nest_engine.cynest.sr('ms Simulate')
        
            ev = nest_engine.cynest.GetStatus(graphics_simulator_spike_detector, "events") # {[times], [senders]}
            l = list(ev[0]["senders"]) + list(ev[0]["times"])
            msg = str(l)
            graphics_simulator_send("[" + str(len(l)) + "," + str(len(msg)) + "]" ) # [nb_spikes + nb_times, length msg]
            
            r = graphics_simulator_receiveKeyword(["ok", "quit", "stop"], 4)
        
            if r == 1:
                return
            elif r == 2:
                keepGoing = False
            
            graphics_simulator_send(msg) # [ id1, id2, ..., idn, t1, t2, ... tn]
            nest_engine.cynest.SetStatus(graphics_simulator_spike_detector, [{"n_events": 0}])
        
            t += graphics_simulator_sim_time_delta
            
                
    nest_engine.cynest.SetKernelStatus({"print_time": pt})
    graphics_simulator_send("finish")
    graphics_simulator_receiveKeyword(["ok"])
    print ("Simulation terminated")
        

def graphics_simulator_close():
    nest_engine.cynest.set_verbosity(graphics_simulator_verbosity_level)
    nest_engine.cynest.SetKernelStatus({"overwrite_files":graphics_simulator_overwrite_files})
    graphics_simulator_sender.close()
    graphics_simulator_listener.close()
    
    






def graphics_simulate(t):
    global graphics_simulator_sim_time_total
    graphics_simulator_sim_time_total = t
    
    
    graphics_simulator_init()
    graphics_simulator_simulate(t)
    graphics_simulator_close()


    
    
    
