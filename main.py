from neuron import h
import argparse

# Local imports
from GA_Loop import GA_Execute_loop 


h.load_file('stdrun.hoc')

#-------Experiment description---------#
# First must be cell type ['single compartment', 'RGC'], 
# Other flags are   'pre_inhibition', 'drifting grating-steady state'

#experiment_type= ['RGC', 'pre_inhibition']
experiment_type= ['RGC']
#--------------------------------------#

# Receptive fields of Excitatory and Inhibitory presynaptic populations
pre_cell_types = ["excitation", "inhibition"]
pre_RF_components = ["center", "surround"]
#----------------------------- arguments (when exist override defaults)
parser = argparse.ArgumentParser(description="Simulation Configuration")
parser.add_argument('--cellType', type=str, default=experiment_type[0], nargs='?', const=0)
parser.add_argument('--job_id', type=int)

experiment_type[0]= parser.parse_args().cellType

global_params = {
    'randomStart':      True,   # random params or load from file
    'numPop':           10,      # population size [10]
    'numGen':           100,      # number of generations

    'cellType':         experiment_type[0],   
    'numSpeed':         1,           # number of probed speeds [5]
    'numContrast':      1,           # number of probed contrasts [2]
    'numRep':           1,      # number of repeats (for example when the stimulus is noisy)
    'numDir':           2,      # number of probed directions, keep 2 or more

    'num_input_types':       4 if ('RGC' in experiment_type) else 2,      # number of different presynaptic clusters with different response profiles 
    
    'switch_to_spiking_model':          1000,   # Set to a number < numGen to produce somatic spikes
    
    'num_syn_replace_input_cell_dist':     100,               # Presynaptic inputs are determined by number of synapses, to cancel - set to zero
    # If    num_syn_replace_input_cell_dist==0, provide the following information:
    'dist_between_input_cells':            28,           # Distance between BCs
    'dist_between_input_cells_SD':         10,            # Some variability so not all BCs are on the same grid

    'mutationRate':      0.1,     # Change in param values between generations
    'analysis_time':     0.7 if ('drifting grating-steady state' in experiment_type) else 0,    # When to start DSI analysis (to avoid initial bumps). Should be zero for moving bars. 
    # RF structue
    "RF_constrains": {
        cell_type: {
            **{
                cs: {
                    'varyAmplitude': False,          # Presynaptic inputs that vary in their strength  
                    'varyKinetics': False,          # Presynaptic inputs that vary in their kinetics
                    'varySize': False,               # Presynaptic inputs that vary in their RF size  
                    'varyOrientation': False,# Presynaptic inputs that vary in their RF orientation
                }
                for cs in pre_RF_components
            },
            'doSurround': False,    # mutate surround (does nothing for center components)
            'depression': False,
            'facilitation': False,
            'depression_tau': False,
            'facilitation_tau': False,
        }
        for cell_type in pre_cell_types
    },

    'pre_inhibition':       'pre_inhibition' in experiment_type,          # include inhibitory presynaptic inputs
    
    # SIMULATION
    'debugger':  {
        'run_neuron': 	    True,	    # Actually run the simulation
        'multithread': 	    False,	    # Use multuple threads or not
        'stop_mutations':   False,       # Do not mutate the models
        'set_speed': 	    -1,	        # Use this speed only (set to zero/negative to disable)
        'set_dir': 	        -1,	        # Use this direction only (set to negative to disable)
        'print_dsi':        True,       # print DSI values 
    },
    'run_debugger':     False,           # RUN DEBUGGER
    
    'job_id':           0,
    'save_every_gen':   10,             # Save RF parameters
    'gen':              0,              # Generation number (set by the GA)
    'fullRFcomputation':False,		    # compute RF activation from exact location info vs using an offset
    'save_all':         False,          # Save all parameters (takes a lot of space)


    'stim_params': {
        'type':             'dg' if ('drifting grating-steady state' in experiment_type) else 'bar', # [bar, dg] # Moving bar, drifting grating,
        'extra':            '',     # Currently of 'vary bar duration' is implemented
        'dt':               2, # if ('drifting grating-steady state' in experiment_type) else 10,       # Time step
        'dx':               10,                     # Spatial step
        'arena':            1000,                   # Visual presentation arena size in microns
        'tStop':            [],                     # Simulation duration
        'speed':            [],      
        'contrast':         [], 
        'angle':            [], 
        'delay':            [200], 
        'duration':         [500] if ('drifting grating-steady state' in experiment_type) else [200],
        'cycle':            [],
    }
}

# Hassenstein-Reichardt model
global_params["RF_constrains"]["excitation"]["center"]["varyKinetics"]= True


if not global_params['run_debugger']:
    global_params['debugger']['run_neuron']= True 
    global_params['debugger']['multithread']= True 
    global_params['debugger']['stop_mutations']= False 
    global_params['debugger']['set_speed']= -1 
    global_params['debugger']['set_dir']= -1 
    global_params['debugger']['print_dsi']= False 


# Reapply arguments
global_params.update({key: val for key, val in vars(parser.parse_args()).items() if val is not None and key in global_params})


if global_params['stim_params']['type'] != 'bar':   # Moving bar can be computed with exact delays, otherwise run full RF computation
    global_params['fullRFcomputation']= True

if not global_params['randomStart']:                # Alert the user that params are pulled from previous run
    print('Not randomized initial parameters!')

if ('single compartment' in experiment_type):   
    global_params['num_syn_replace_input_cell_dist']= 0

if __name__ == "__main__":
    models=[]       # Container of the NEURON models
    GA_Execute_loop(models,  global_params)

