from multiprocessing import Process, Queue, freeze_support, Pipe
import numpy as np
import time
import matplotlib.pyplot as plt
from neuron import h
#from time import sleep
from GA_RF import GA_Pre_activation , GA_StimTrajectory
from GA_h5 import GA_SaveH5
from neuron import h
import copy
import pickle
#import sys
import argparse


h.load_file('stdrun.hoc')

from GA_NEURON import NEURON_SetCells, NEURON_UpdateSynapses, NEURON_mutation, GA_RecordingVectors

h.load_file('stdrun.hoc')
# Receptive fields of Excitatory and Inhibitory presynaptic populations
pre_cell_types = ["excitation", "inhibition"]
pre_RF_components = ["center", "surround"]

global_params = {
    'randomStart':      True,   # random params or load from file
    'numPop':           10,     # population size
    'numGen':           300,      # number of generations
    'activeChannels':   '',#'Ka,CaN',   #CaN,CaL,caLinear, Ka,Km,IH,Kdr
    'cellType':         'RGC',   #['RGC', 'L23', 'L5', 'SAC', 'SAC network'],      # Cell type (RGC, L23, L5, SAC, SAC-simple morphology)
    'numSpeed':         3,      # number of probed speeds
    'numContrast':      1,      # number of probed contrasts
    'numDir':           2,      # number of probed directions, keep 2 or more

    'numPreClus':       4,      # number of different presynaptic clusters with different response profiles 
    'numPostComp':      4,      # number of different postsynaptic compartments with possible different distributions of  voltage gated channels
    'numVoltagePoints': 20,      # number of different postsynaptic voltage points to modify the voltage gated channels, starting at -80 and increse by 5mV
    #compartments with possible different distributions of  voltage gated channels
    
    'numSyn':           100,    # number of synaptic inputs  
    'distSyn':          0,     # distance in microns between synapses, set to zero or negative to use random numSyn only
    'dt':               10,     # Time step
    'mutationRate':     0.1,    # Change in param values between generations
    # RF structue
    "RF_constrains": {
        cell_type: {
            cs: {
                'sameKinetics': False,  # Presynaptic inputs that vary in their kinetics
                'sameSize': False,       # Presynaptic inputs that vary in their RF size  
                'sameAmplitude': False,       # Presynaptic inputs that vary in their strength  
                'sameOrientation': True,# Presynaptic inputs that vary in their RF orientation
                'doSurround': False,    # mutate surround (does nothing for center components)
            }
            for cs in pre_RF_components
        }
        for cell_type in pre_cell_types
    },

    'inhibition':       True,          # include inhibitory presynaptic inputs
    
    # SAC NETWORK ONLY
    'numSAClayersX':    1,		# how many SACs are arranged around the target cell; 0 will create one cell
    'numSAClayersY': 	1,	    # how many SACs are arranged around the target cell in Y; 1 is minimal number
    'distSAC': 	        25,		# Distance between cells
    'maxSAC_syndist': 	15,	    # highest possible distance between pre and post synaptic locations
    # SIMULATION
    'debugger':  {
        'run_neuron': 	    True,	    # Actually run the simulation
        'multithread': 	    False,	    # Use multuple threads or not
        'stop_mutations':   False,       # Do not mutate the models
        'set_speed': 	    1,	        # Use this speed only (set to zero/negative to disable)
        'set_dir': 	        -1,	        # Use this direction only (set to negative to disable)
        'plot_outcome': 	True,	    # Plot simulation results
        'plot_inputs':      False,      # Plot the inputs
        'plot_positions':   False,       # Plot the cell positions
        'print_dsi':        True,       # print DSI values for all cells
    },
    #'run_debugger':     True,       # RUN DEBUGGER
    'run_debugger':     False,       # TRUE  to run normal sim, otherwise apply the debugger
    'job_id':           0,
    'save_every_gen':   10,
    'numCells':	        1,		# number of cells
    'ca_present':	    False,		# ca conductance exists
    'k_present':	    False,		# k conductance exists
    'save_all':         False,      # Save all parameters (takes a lot of space)
    
}
if not global_params['run_debugger']:
    global_params['debugger']['run_neuron']= True 
    global_params['debugger']['multithread']= True 
    global_params['debugger']['stop_mutations']= False 
    global_params['debugger']['set_speed']= -1 
    global_params['debugger']['set_dir']= -1 
    global_params['debugger']['plot_outcome']= False
    global_params['debugger']['plot_inputs']= False 
    global_params['debugger']['plot_positions']= False 
    global_params['debugger']['print_dsi']= False 
    

#----------------------------- arguments!
parser = argparse.ArgumentParser(description="Simulation Configuration")
parser.add_argument('--cellType', type=str)
parser.add_argument('--activeChannels', type=str)
parser.add_argument('--numSAClayersX', type=int)
parser.add_argument('--distSAC', type=int)
parser.add_argument('--job_id', type=int)


global_params.update({key: val for key, val in vars(parser.parse_args()).items() if val is not None and key in global_params})
if len(global_params['activeChannels']) > 0 :
    print(f"active channels: {global_params['activeChannels']}")

stim_params = {
    'dt':               global_params['dt'], 
    'dx':               10, 
    'arena':            1000, 
    'tStop':            3000, 
    'speed':            1, 
    'contrast':         1, 
    'angle':            0, 
    'delay':            500, 
    'duration':         200
}

class Model:
    pass
    def __init__(self):
        self.InputE_synapses = []   # Synaptic inputs from presynaptic populations that can have difefrent RF properties
        self.InputI_synapses = []
        self.SAC_SAC_synapses = []  # Feedback inhibition from pther SACs
        self.cell = []              # Morhphology etc
        self.input_params = {}
        self.output_params = {}


def GA_Run(model, stim_params, global_params):
    h(f"forall g_pas={float(model.input_params['passive_params']['pas'])}")
    h(f"forall Ra={float(model.input_params['passive_params']['Ra'])}") 
    # Flexible active conductances
    v_vec= np.linspace(-80, -80+5*global_params['numVoltagePoints'], global_params['numVoltagePoints']) # Voltage stops
    V_Vec = h.Vector(v_vec)
    #mtau_vec= h.Vector(model.input_params['active_params']['mtau_caGA'])
    #print('presetn',V_Vec.size(), mtau_vec.size())  # Should return True

    if "caGA" in global_params['activeChannels']:     # Flexible conductance present    
        h.table_minf_caGA(h.Vector(model.input_params['active_params']['minf_caGA']), V_Vec)           
        h.table_mtau_caGA(h.Vector(model.input_params['active_params']['mtau_caGA']), V_Vec)           
        h.table_hinf_caGA(h.Vector(model.input_params['active_params']['hinf_caGA']), V_Vec)           
        h.table_htau_caGA(h.Vector(model.input_params['active_params']['htau_caGA']), V_Vec)           
    if "kGA" in global_params['activeChannels']:     # Flexible conductance present    
        h.table_ninf_kGA(h.Vector(model.input_params['active_params']['minf_kGA']), V_Vec)           
        h.table_ntau_kGA(h.Vector(model.input_params['active_params']['mtau_kGA']), V_Vec)           
        h.table_hinf_kGA(h.Vector(model.input_params['active_params']['hinf_kGA']), V_Vec)           
        h.table_htau_kGA(h.Vector(model.input_params['active_params']['htau_kGA']), V_Vec)           


    for cell in range(global_params['numCells']):
        for d, sec in enumerate(model.cell[cell].all):
            #sec.push()    
            type= model.output_params['cell'][cell]['activeCompartment']['type'][d]
            if h.ismembrane("canrgc", sec=sec):     # N-type conductance present
                for seg in sec:
                    seg.gbar_canrgc = model.input_params['active_params']['gbar_canrgc'][type].item()
                    seg.shift_canrgc = model.input_params['active_params']['shift_canrgc']
            if h.ismembrane("calrgcfix", sec=sec):     # L-type conductance present
                for seg in sec:
                    seg.gbar_calrgc = model.input_params['active_params']['gbar_calrgc'][type].item()
                    seg.shift_calrgc = model.input_params['active_params']['shift_calrgc']
            if h.ismembrane("caGA", sec=sec):     # Flexible conductance present
                for seg in sec:
                    seg.gMax_caGA = model.input_params['active_params']['gMax_caGA'][type].item()
                    seg.mN_caGA = model.input_params['active_params']['mN_caGA']

            
            if h.ismembrane("kap", sec=sec):        # A-type conductance present
                for seg in sec:
                    seg.gkabar_kap = model.input_params['active_params']['gbar_ka'][type].item()
                    seg.shift_kap = model.input_params['active_params']['shift_ka']
            if h.ismembrane("km", sec=sec):        # M-type conductance present
                for seg in sec:
                    seg.gbar_km = model.input_params['active_params']['gbar_km'][type].item()
                    seg.shift_km = model.input_params['active_params']['shift_km']
            if h.ismembrane("kSlow", sec=sec):        # DR-type conductance present
                for seg in sec:
                    seg.gkbar_kSlow = model.input_params['active_params']['gbar_kd'][type].item()
                    seg.v_shift_kSlow = model.input_params['active_params']['shift_kd']
            if h.ismembrane("ih", sec=sec):        # M-type conductance present
                for seg in sec:
                    seg.ghdbar_ih = model.input_params['active_params']['gbar_ih'][type].item()
                    seg.shift_ih = model.input_params['active_params']['shift_ih']
            if h.ismembrane("kGA", sec=sec):     # Flexible conductance present
                for seg in sec:
                    seg.gMax_kGA = model.input_params['active_params']['gMax_kGA'][type].item()
                    seg.nN_kGA = model.input_params['active_params']['nN_kGA']


    GA_RecordingVectors(global_params, stim_params, model, prep= True)


    # Run multiple simulations for a range of contrasts and velocities
    for speed in range(global_params['numSpeed']):
        # Compute speed
        stim_params['speed']= 1
        if(global_params['numSpeed'] == 5):
            stim_params['speed']=  2**(speed - 2)
        if(global_params['numSpeed'] == 3):
            stim_params['speed']=  4**(speed - 1)
        if(global_params['numSpeed'] == 2):
            stim_params['speed']=  0.25 + 0.75 * speed
        if(global_params['debugger']['set_speed'] > 0):
            stim_params['speed']= global_params['debugger']['set_speed']
        
        # Computation of activation params
        model.input_params['trajectory'] = GA_StimTrajectory(stim_params)
        for cls in range(global_params['numPreClus']):
            model.input_params['InputE_syn_time'][cls, :] = GA_Pre_activation('excitation', model.input_params['RF_params'], model.input_params['trajectory'], stim_params, pop=0, cls=cls)
            if global_params['inhibition']:
                model.input_params['InputI_syn_time'][cls, :] = GA_Pre_activation('inhibition', model.input_params['RF_params'], model.input_params['trajectory'], stim_params, pop=0, cls=cls)
        
        for contrast in range(global_params['numContrast']):
            stim_params['contrast']= (3**(contrast + 1)) / 3**global_params['numContrast']
            # Run for different directions, compute DSI and store it
            score_right= []
            score_left= []
            model.output_params['cell'][cell]['max_soma_single_run']= []
            for cell in range(global_params['numCells']):
                model.output_params['cell'][cell]['dendrite_vectors']['max_right_single_run']= []
                model.output_params['cell'][cell]['dendrite_vectors']['max_left_single_run']= []
            for dr in range (global_params['numDir']):
                
                stim_params['angle'] = (dr / (global_params['numDir']-1)) * np.pi                       
                if(global_params['debugger']['set_dir'] >= 0):
                    stim_params['angle']= global_params['debugger']['set_dir']
                #for cell in range(global_params['numCells']):
                score_right.append(np.cos(stim_params['angle']))
                score_left.append(-np.cos(stim_params['angle']))
                h.tstop= NEURON_UpdateSynapses(global_params, stim_params, model)

                if(not global_params['debugger']['run_neuron']):
                    h.tstop= speed + 10
                h.run()
                GA_RecordingVectors(global_params, stim_params, model, populate= True)


            # Compute DSI  after all dirs are done - from somatic responses for RGCs and dendritic calcium levels for SACs
            
            for cell in range(global_params['numCells']):
                if (global_params['cellType'] in ['RGC', 'L23', 'L5']): # Somatic 
                    dsi_list=np.dot(np.array(score_right), np.array(model.output_params['cell'][cell]['max_soma_single_run'])) / sum( model.output_params['cell'][cell]['max_soma_single_run']) 
                else:   # SAC
                    dsi_list= np.dot(np.array(score_right), np.array(model.output_params['cell'][cell]['dendrite_vectors']['max_right_single_run'])) / sum( model.output_params['cell'][cell]['dendrite_vectors']['max_right_single_run']) / 2
                    dsi_list+= np.dot(np.array(score_left), np.array(model.output_params['cell'][cell]['dendrite_vectors']['max_left_single_run'])) / sum( model.output_params['cell'][cell]['dendrite_vectors']['max_left_single_run']) / 2
                model.output_params['cell'][cell]['dsi_list'].append(dsi_list)

                if global_params['debugger']['print_dsi']:  # Report DSI values from all cells
                    print(np.array(model.output_params['cell'][cell]['max_soma_single_run']), dsi_list)
                
    for cell in range(global_params['numCells']):
        model.output_params['cell'][cell]['dsi']=np.mean(model.output_params['cell'][cell]['dsi_list'])
    

def GA_Start(conn, model, stim_params, global_params):
    h('v_init= -60')
    h('stdinit()')
    NEURON_SetCells(global_params, stim_params, model)   
    

    while True:
        if conn.poll():  # Check for incoming command
            msg = conn.recv()
            response = "---"
            if isinstance(msg, str):
                if msg == "input_params":
                    response = model.input_params
                if msg == 'global_params':
                    response= global_params
                if msg == 'I_syn':
                    response= len(model.SAC_SAC_synapses)
                elif msg == "quit":
                    #h.quit()
                    conn.send("end of simulation")
                    break
                elif msg == "show_inputs":
                    maxx= -10000
                    minx= 10000
                    for syn in model.InputE_synapses:
                        maxx= max(maxx, syn.x)
                        minx= min(minx, syn.x)
                    for syn in model.InputE_synapses:
                        plt.plot(syn.shifted_drive_vector.to_python(), color= ((syn.x - minx) / (maxx - minx), 0, 0))
                    plt.show()   
                    response= "done plotting"
                    #break            
            elif isinstance(msg, dict):
                # Apply the passive and active params, compute the synaptic activation and run the simulation               
                model.input_params = msg
                GA_Run(model, stim_params, global_params)
                response = model.output_params
            conn.send(response) 

if __name__ == "__main__":
    
    freeze_support()
    start_time = time.time()  # ⏱️ Start timer
    # Initialize the models
    models=[]       # Container of the NEURON models
    for pop in range(global_params['numPop']):
        model = Model()
        models.append(model)
    
    sim_procs = []
    sim_conns = []
    # Initialize the NEURON threads
    for pop in range(global_params['numPop']):      # Adjust number of NEURON simulations here
        if(global_params['debugger']['multithread']):
            parent_conn, child_conn = Pipe()
            p = Process(target=GA_Start, args=(child_conn, models[pop], stim_params, global_params))
            p.start()
            sim_procs.append(p)
            sim_conns.append(parent_conn)
        else:                   # Single thread
            h('load_file("nrngui.hoc")')
            NEURON_SetCells(global_params, stim_params, model)
            h('load_file("neuron.ses")')
            h('v_init= -60')


    # Save a copy of the created params
    if(global_params['debugger']['multithread']):
        for pop,conn in enumerate(sim_conns):
            conn.send("input_params")
            models[pop].input_params= conn.recv()
        sim_conns[0].send('global_params')
        global_params= sim_conns[0].recv()
        #sim_conns[0].send('I_syn')
        #nSynI= sim_conns[0].recv()  

    if(global_params['numCells'] > 1):
        print(f"Number of SAC-SAC synapses - {global_params['numIsyn']}  , Number of cells - {global_params['numCells']}")
    
    models[0].output_params['score']= []
    # Generation loop
    for gen in range(global_params['numGen']):    # number of generations
        gen_time = time.time()  

        # Send new input parameters and run the simulations
        if(global_params['debugger']['multithread']):
            for pop, conn in enumerate(sim_conns):    
                #print('send', pop, models[pop].input_params['RF_params']['excitation']['center']['peak'])     
                conn.send(models[pop].input_params)

            # Gather values
            for pop,conn in enumerate(sim_conns):
                models[pop].output_params= conn.recv()
                conn.send("input_params")
                models[pop].input_params= conn.recv()
                #print('back',pop, models[pop].input_params['RF_params']['excitation']['center']['peak'])     
        else:
            for pop in range(global_params['numPop']):
                GA_Run(models[pop], stim_params, global_params)
       
        # Find the model with the largest DSI
        best_dsi= -1
        best_pos= 0
        for pop, model in enumerate(models):
            if(model.output_params['cell'][0]['dsi'] > best_dsi):
                best_dsi= model.output_params['cell'][0]['dsi']
                best_pos= pop 
            if global_params['debugger']['print_dsi']:
                print(f"Model={pop}, DSI={model.output_params['cell'][0]['dsi']}. Best pos={best_pos} ({best_dsi})")            
       
        print(f"Done generation {gen}, best score= {best_dsi}, time - {(time.time()-gen_time):.4f}")
        models[0].output_params['score'].append(best_dsi) # Save progress

        #print(f"Gen - {gen}, time - {(time.time()-gen_time):.4f}")
        if((gen%(global_params['save_every_gen'])) == (global_params['save_every_gen']-1)):
            GA_SaveH5(global_params, models, final= False)
            with open(f"params/input_params_{global_params['job_id']}.pkl", "wb") as f:
                pickle.dump(models[0].input_params, f)
        
        # Duplicate best model
        if gen < global_params['numGen'] - 1:            
            for pop, model in enumerate(models): 
                if pop != best_pos:
                    models[pop] = copy.deepcopy(models[best_pos])
                    if global_params['debugger']['print_dsi']:
                        print(f'Copy {best_pos} into {pop}')
            
            # Mutations
            if (global_params['debugger']['stop_mutations'] == False):
                for pop, model in enumerate(models):
                    if(pop > 0):    # keep the first model intact
                        NEURON_mutation(global_params, model.input_params)
                        if global_params['debugger']['print_dsi']:
                            print(f'Mutate {pop}')                    

    # End generation loop
    #models[0].output_params['score']= score # Save progress
    if (global_params['debugger']['plot_inputs']):
        sim_conns[0].send("show_inputs")
        dummy= conn.recv()
    
    # Cleanup
    if(global_params['debugger']['multithread']):
        for conn in sim_conns:
            conn.send("quit")
            dummy= conn.recv()
        for p in sim_procs:
            p.join()

    print("All NEURON simulations completed.")
    end_time = time.time()  # ⏱️ End timer
    elapsed_time = end_time - start_time
    
    
    # Save to HDF5
    GA_SaveH5(global_params, models, final= True)
    print("Saved H5 File")
    with open(f"params/input_params_{global_params['job_id']}.pkl", "wb") as f:
        pickle.dump(models[0].input_params, f)
    

    
    print(f"Execution Time: {(elapsed_time / 60):.2f} minutes")  

    if (global_params['debugger']['plot_outcome']):
        #print("op ",models[0].output_params)
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 6), height_ratios=[1, 1, 1])
        ax3.plot(models[0].output_params['score'], color='black')
        for cell in range(global_params['numCells']): 
            for v in model.output_params['cell'][cell]['somaV_all_angles']:
                ax1.plot(v)
        if (global_params['cellType'] in ['RGC', 'L23', 'L5']): # Somatic 
            pass

        else:
            for ca in model.output_params['cell'][0]['dendrite_vectors']['ca_right_all_angles']:
                ax2.plot(ca, color= 'orange')

            for ca in model.output_params['cell'][0]['dendrite_vectors']['ca_left_all_angles']:
                ax2.plot(ca, color= 'black')

        plt.show() 

        #for cell in range(global_params['numCells']): 
    if (global_params['debugger']['plot_positions']):
        plt.plot(models[0].input_params['cellPos'][:][0],models[0].input_params['cellPos'][:][1])
        x_coords, y_coords = zip(*models[0].input_params['cellPos'])

        plt.scatter(x_coords, y_coords)
        plt.show() 
     
