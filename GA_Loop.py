import numpy as np
from neuron import h
import time
from multiprocessing import Process, freeze_support, Pipe
import copy
import pickle
import math 

#from GA_Loop import GA_Run, GA_Start, NEURON_mutation
from GA_NEURON import Model, NEURON_SetCells 
from GA_h5 import GA_SaveH5

pre_cell_types = ["excitation", "inhibition"]
pre_RF_components = ["center", "surround"]

from GA_NEURON import NEURON_SetCells, GA_RecordingVectors
from GA_BC import GA_UpdateInputSynapses
from GA_RF import GA_StimTrajectory, GA_Pre_activation

# Introduce mutations and check boundaries
def NEURON_mutation(global_params, input_params, output_params):
    np.random.seed()
    # RF params
    for cell_type in pre_cell_types:    # Excitation/Inhibition     
        for cs in pre_RF_components:    # Center/Surround        
            for key, value in input_params['RF_params'][cell_type][cs].items():
                # Introduce mutations
                if(key in ['widthC']):  # rotation has no special meaning for zero
                    input_params['RF_params'][cell_type][cs][key]= value + np.random.normal(0, global_params['mutationRate'] / 10, (global_params['num_input_types'], 1)) # [[0],[.7],[1.1],[1.7]]#
                else:
                    input_params['RF_params'][cell_type][cs][key]= value * np.random.normal(1, global_params['mutationRate'], (global_params['num_input_types'], 1)) + np.random.uniform(0, global_params['mutationRate'] / 10000, (global_params['num_input_types'], 1))
                
                # Impose boundaries on the RF parameters
                if(key in ['widthX', 'widthY']):
                    input_params['RF_params'][cell_type][cs][key]= np.clip(input_params['RF_params'][cell_type][cs][key], 10, 200 + 2000 * (cs == 'surround'))

                if(key in ['tau', 'tauRRP']):
                    input_params['RF_params'][cell_type][cs][key]= np.clip(input_params['RF_params'][cell_type][cs][key], global_params['stim_params']['dt'], 10000) 
                
                if((key == 'peak') and (cs == 'surround')):
                    input_params['RF_params'][cell_type][cs][key]= np.clip(input_params['RF_params'][cell_type][cs][key], 0, 1) 
          
            # Modify the RF params based on constrains (when params are shared, they get the value of the first input)
            if not global_params['RF_constrains'][cell_type][cs]['varyKinetics']:        # similar kinetics
                for key in ['tau', 'tauRRP']:
                    input_params['RF_params'][cell_type][cs][key].fill(input_params['RF_params'][cell_type][cs][key][0]) 
            
            if not global_params['RF_constrains'][cell_type][cs]['varySize']:             # similar RF size
                for key in ['widthX', 'widthY']:
                    input_params['RF_params'][cell_type][cs][key].fill(input_params['RF_params'][cell_type][cs][key][0]) 
            
            if not global_params['RF_constrains'][cell_type][cs]['varyOrientation']:             # similar RF orientation
                input_params['RF_params'][cell_type][cs]['widthC'].fill(input_params['RF_params'][cell_type][cs]['widthC'][0]) 
                for cls in range(global_params['num_input_types']):
                    input_params['RF_params'][cell_type][cs]['widthY'][cls]=input_params['RF_params'][cell_type][cs]['widthX'][cls]  # Circle and not oval
                
            if not global_params['RF_constrains'][cell_type][cs]['varyAmplitude'] :                # similar RF strength
                input_params['RF_params'][cell_type][cs]['peak'].fill(input_params['RF_params'][cell_type][cs]['peak'][0]) 
                   
        # Disable surround
        if not global_params['RF_constrains'][cell_type]['doSurround']:        # No surround
            input_params['RF_params'][cell_type]['surround']['peak']*= 0 
        
        # Short term plasticity
        for key in ['depression', 'facilitation', 'depression_tau', 'facilitation_tau']:
            if global_params['RF_constrains'][cell_type][key]:
                input_params['RF_params'][cell_type][key]*= np.random.normal(1, global_params['mutationRate'], (global_params['num_input_types'], 1)) 
                input_params['RF_params'][cell_type][key]= np.clip(input_params['RF_params'][cell_type][key], 0, 1)
            else:
                input_params['RF_params'][cell_type][key].fill(input_params['RF_params'][cell_type][key][0])   
    # Passive params
    for key, value in input_params['passive_params'].items():
        input_params['passive_params'][key]= value * np.random.normal(1, global_params['mutationRate']) + np.random.uniform(0, global_params['mutationRate'] / 10000)
        if(key == 'Ra'):
           input_params['passive_params'][key]= np.clip(input_params['passive_params'][key], 50, 300) 
    
    # Inputs
    input_params["input_exc_syn_gain"]= input_params["input_exc_syn_gain"] * np.random.normal(1, global_params['mutationRate'] ) + np.random.uniform(0, global_params['mutationRate'] / 1000)
    input_params["input_inh_syn_gain"]= input_params["input_inh_syn_gain"] * np.random.normal(1, global_params['mutationRate'] ) + np.random.uniform(0, global_params['mutationRate'] / 1000)
    
    # Gain from photoreceptors
    input_params["fromPhsynG"]= input_params["fromPhsynG"] * np.random.normal(1, global_params['mutationRate'] ) 

    
# If the cell did not spike, increase stimulation intensity
def NEURON_spike_reached(model):
    # tweaks
    if model.output_params['cell']['dsi'] == 0:
        model.input_params["input_exc_syn_gain"]*= 2
        print(f'Increasing BC input to {model.input_params["input_exc_syn_gain"]}')
        
        
# Apply parameters
def GA_Run(model, global_params):
    
    # Add active params (when spikes are implemented)
    if model.input_params['gen'] == global_params['switch_to_spiking_model']:            # APs at the soma
        print("Switching to spiking DS")
        model.cell[0].soma.push()   # Place spiking mechanism at the soma
        h('celsius = 22')  
        h('insert na16')  
        h('insert kSlow')  
        h('insert kv')  
        h('vshift_na16= -15')  

        if(global_params['cellType'] == 'single compartment'):
            h('gbar_na16= 2000')  
            h('gbar_kv= 400')  
            h('gkbar_kSlow= 200')  
        else:                                            # RGC
            h('gbar_na16= 10000')  
            h('gbar_kv= 2000')  
            h('gkbar_kSlow= 000')   
        h.pop_section()    
        

    # Passive parameters
    for d, sec in enumerate(model.cell.all):
        sec.push()    
        h(f"g_pas={float(model.input_params['passive_params']['pas'])}")
        h(f"Ra={float(model.input_params['passive_params']['Ra'])}") 
        h.pop_section()

    GA_RecordingVectors(global_params, model, prep= True)

    # Run multiple simulations for a range of contrasts and velocities
    for r in range(global_params['numRep']):                        # Multiple repeats
        for s in range(global_params['numSpeed']):                  # Multiple speeds
            # Compute speed
            speed= 1        # Single speed
            cycle= 400      # Single drifting grating level
            duration= 200   # Bar duration
            if(global_params['numSpeed'] == 5):
                speed=  2**(s - 2)
            if(global_params['numSpeed'] == 3):
                speed=  4**(s - 1)
            if(global_params['numSpeed'] == 2):
                speed=  0.25 + 0.75 * s
            if(global_params['debugger']['set_speed'] > 0):
                speed= global_params['debugger']['set_speed']
            
            # Compute cycle (when activated with drifting grating)
            if global_params['stim_params']['type'] == 'dg':
                speed=  1
                if(global_params['numSpeed'] == 5):
                    #cycles= [ 100, 200, 400, 800, 1600]
                    cycle= 100 * ( 2**s ) 
                if(global_params['numSpeed'] == 2):
                    cycle= 100 + 400 * s

            # vary bar duration
            if global_params['stim_params']['extra'] == 'vary bar duration':
                speed=  1
                if(global_params['numSpeed'] == 5):
                    duration= 50 * (2**s)
                if(global_params['numSpeed'] == 2):
                    duration= 50 + (950*s)

            global_params['stim_params']['speed'].append( speed )
            global_params['stim_params']['cycle'].append( cycle )
            global_params['stim_params']['duration'].append( duration )

            for contrast in range(global_params['numContrast']):    # Multiple contrasts
                global_params['stim_params']['contrast'].append((3**(contrast + 1)) / 3**global_params['numContrast'])
                # Run for different directions, compute DSI and store it
                score_right= []
                score_left= []
                dsi_list= []
                adjusted_dsi_list= []            
                model.output_params['cell']['max_soma_single_run']= []
                
                # Directions
                for dr in range (global_params['numDir']):
                    if global_params['numDir'] > 1:
                        angle = (dr / (global_params['numDir']-1)) * np.pi                       
                    else:
                        angle= 0
                    if(global_params['debugger']['set_dir'] >= 0):
                        angle= global_params['debugger']['set_dir']
                    global_params['stim_params']['angle'].append(angle)
                    
                    # compute the stimulus
                    GA_StimTrajectory(global_params, model.input_params, repeat= r) # Stimulus trajectory
                    
                    if global_params['fullRFcomputation']:
                        list_range= len(model.input_exc_list)           # each synapse is computed independently
                    else:
                        list_range= global_params['num_input_types']    # input is divided into clusters, vectors are time shifted for moving bar presentation

                    input_population= [model.input_exc_list]
                    input_type= ['excitation']
                    input_time=['input_exc_syn_time']
                    if global_params['pre_inhibition']:
                        input_population.append(model.input_inh_list)
                        input_type.append('inhibition')
                        input_time.append('input_inh_syn_time')
                    
                    for ei in range(len(input_type)):
                        for counter in range(list_range):
                            cluster= counter
                            xcenter= 0
                            if global_params['fullRFcomputation']:
                                cluster= input_population[ei][counter].cluster
                                xcenter= input_population[ei][counter].x
                            
                            pre_drive, spatial_sum_center, activation_center, inactivation_center, spatial_sum_surround, activation_surround, inactivation_surround= GA_Pre_activation(input_type[ei], model.input_params, global_params, cls=cluster, xcenter= xcenter)
                            #print(counter, len(model.input_exc_list), len(model.input_params[input_time[ei]]))
                            if counter >= len(model.input_params[input_time[ei]]):
                                template_item = {
                                    'drive': [],
                                    'original': [],
                                    'activation': [],
                                    'inactivation': [],
                                    'original_s': [],
                                    'activation_s': [],
                                    'inactivation_s': [],
                                }
                                model.input_params[input_time[ei]].append({k: [] for k in template_item})
                            model.input_params[input_time[ei]][counter]['drive'].append(pre_drive)                           # Save the RF temporal drive waveforms
                            model.input_params[input_time[ei]][counter]['original'].append(spatial_sum_center)               # Save the convolved(stimulus*RF) waveforms
                            model.input_params[input_time[ei]][counter]['activation'].append(activation_center)              # Save the RF activation waveforms
                            model.input_params[input_time[ei]][counter]['inactivation'].append(inactivation_center)          # Save the RF inactivation waveforms                            
                            model.input_params[input_time[ei]][counter]['original_s'].append(spatial_sum_surround)           # Save the surround convolved (stimulus, RF) waveforms
                            model.input_params[input_time[ei]][counter]['activation_s'].append(activation_surround)          # Save the RF surround activation waveforms
                            model.input_params[input_time[ei]][counter]['inactivation_s'].append(inactivation_surround)      # Save the RF surround inactivation waveforms
                            
                            if global_params['fullRFcomputation']:
                                input_population[ei][counter].drive_vector= pre_drive
                        # apply the waveforms if in cluster mode
                        for bc in input_population[ei]:
                            bc.drive_vector= model.input_params[input_time[ei]][bc.cluster]['drive'][-1]
                            bc.input_vector= model.input_params[input_time[ei]][bc.cluster]['original'][-1]                        

                    score_right.append(math.cos(global_params['stim_params']['angle'][-1]))
                    score_left.append(-math.cos(global_params['stim_params']['angle'][-1]))
                    GA_UpdateInputSynapses(global_params, model)
                    h.tstop= global_params['stim_params']['tStop'][-1]
                    if(not global_params['debugger']['run_neuron']):
                        h.tstop= 10
                    h.finitialize(-60)
                    h.run()
                    GA_RecordingVectors(global_params, model, populate= True)

                # Compute DSI  after all dirs are done                
                sum_run= sum( model.output_params['cell']['max_soma_single_run'] )
                if sum_run == 0:
                    dsi_list.append(0)
                    adjusted_dsi_list.append(0)
                    if global_params['debugger']['print_dsi']:  # Report DSI values from all cells
                        print('no signal')      # for example, where there are no spikes
                else:
                    dsi_top= 0
                    for sim in range(len(score_right)):
                        dsi_top+= score_right[sim] * model.output_params['cell']['max_soma_single_run'][sim] #// sum_run
                    dsi_list.append( dsi_top / sum_run)
                    # GA is trained based on adjusted DSI that takes into account response amplitude
                    adjusted_dsi_list.append((dsi_top / sum_run) * (math.tanh(0.2 * max( model.output_params['cell']['max_soma_single_run']) )))

                    model.output_params['cell']['dsi_list'].append(dsi_list)
                    model.output_params['cell']['adjusted_dsi_list'].append(adjusted_dsi_list)
                    if global_params['debugger']['print_dsi']:  # Report DSI values from all cells
                        print(np.array(model.output_params['cell']['max_soma_single_run']), dsi_list, 'adjusted', adjusted_dsi_list)
    
    # Find the final DSI from all speeds/contrasts/etc
    model.output_params['cell']['dsi']= np.mean(model.output_params['cell']['dsi_list'])
    model.output_params['cell']['adjusted_dsi']= np.mean(model.output_params['cell']['adjusted_dsi_list'])
    
    if np.isnan(model.output_params['cell']['dsi']):    # No data
        model.output_params['cell']['dsi']= -1
        model.output_params['cell']['adjusted_dsi']= -1

# Start configurations
def GA_Start(conn, model, global_params):
    h('v_init= -60')
    h('stdinit()')
    NEURON_SetCells(global_params, model)   
    if global_params['randomStart']:    # Apply boundaries when starting from random initial values
        NEURON_mutation(global_params, model.input_params, model.output_params)  

    while True:
        if conn.poll():  # Check for incoming command
            msg = conn.recv()
            response = "---"
            if isinstance(msg, str):
                if msg == "input_params":
                    response = model.input_params
                if msg == "global_params":
                    response = global_params
                if msg == 'I_syn':
                    response= len(model.SAC_SAC_I_synapses)
                elif msg == "quit":
                    #h.quit()
                    conn.send("end of simulation")
                    break
        
            elif isinstance(msg, dict):
                # Apply the passive and active params, compute the synaptic activation and run the simulation               
                model.input_params = msg
                GA_Run(model,  global_params)
                response = model.output_params
            conn.send(response) 

# Execute the main loop
def GA_Execute_loop(models,  global_params):
    freeze_support()
    start_time = time.time()  # ⏱️ Start timer
    # Initialize the models
    for pop in range(global_params['numPop']):
        model = Model()
        models.append(model)
    
    sim_procs = []
    sim_conns = []
    # Initialize the NEURON threads
    for pop in range(global_params['numPop']):      # Adjust number of NEURON simulations here
        if(global_params['debugger']['multithread']):
            parent_conn, child_conn = Pipe()
            p = Process(target=GA_Start, args=(child_conn, models[pop],  global_params))
            p.start()
            sim_procs.append(p)
            sim_conns.append(parent_conn)
        else:                   # Single thread
            h('load_file("nrngui.hoc")')
            NEURON_SetCells(global_params,  models[pop])
            h('load_file("neuron.ses")')
            h('v_init= -60')


    # Save a copy of the created params
    if(global_params['debugger']['multithread']):
        for pop,conn in enumerate(sim_conns):
            conn.send("input_params")
            models[pop].input_params= conn.recv()
        sim_conns[0].send('global_params')
        global_params= sim_conns[0].recv()

    
    score= []
    # Generation loop
    for global_params['gen'] in range(global_params['numGen']):    # number of generations
        gen_time = time.time()  
        
        for pop in range(global_params['numPop']):
            models[pop].input_params['gen']= global_params['gen']

        # Send new input parameters and run the simulations
        if(global_params['debugger']['multithread']):
            for pop, conn in enumerate(sim_conns):    
                conn.send(models[pop].input_params)

            # Gather values
            for pop,conn in enumerate(sim_conns):
                models[pop].output_params= conn.recv()
                conn.send("input_params")
                models[pop].input_params= conn.recv()
        else:
            for pop in range(global_params['numPop']):
                GA_Run(models[pop], global_params)
       
        # Find the model with the largest DSI
        best_adjusted_dsi= -1
        best_dsi= -1
        best_pos= 0
        for pop, model in enumerate(models):
            if(model.output_params['cell']['adjusted_dsi'] > best_adjusted_dsi):
                best_adjusted_dsi= model.output_params['cell']['adjusted_dsi']
                best_dsi= model.output_params['cell']['dsi']
                best_pos= pop 
            if global_params['debugger']['print_dsi']:
                print(f"Model={pop}, DSI={model.output_params['cell']['dsi']}. Best pos={best_pos} ({best_dsi})")            
       
        print(f"Done generation {global_params['gen']} out of {global_params['numGen']}, best score= {best_dsi} [{best_adjusted_dsi}], time - {(time.time()-gen_time):.4f}")
        score.append(best_dsi) # Save progress


        if((global_params['gen']%(global_params['save_every_gen'])) == (global_params['save_every_gen']-1)):
            models[0].output_params['score']= score
            GA_SaveH5(global_params, models, final= False)
        if False:   # optinal-save input params
            with open(f"params/input_params_{global_params['job_id']}.pkl", "wb") as f:
                pickle.dump(models[0].input_params, f)
        
        # Duplicate best model
        if global_params['gen'] < global_params['numGen'] - 1:            
            for pop, model in enumerate(models): 
                if pop != best_pos:
                    models[pop] = copy.deepcopy(models[best_pos])

            
            # Mutations
            if (global_params['debugger']['stop_mutations'] == False):
                for pop, model in enumerate(models):
                    if(pop > 0):    # keep the first model intact
                        NEURON_mutation(global_params, model.input_params, model.output_params)
                    NEURON_spike_reached(model) # Check that the cell spikes (if in spiking regime)                   

    # End generation loop
    models[0].output_params['score']= score # Save progress   
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