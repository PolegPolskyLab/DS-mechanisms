import numpy as np
from neuron import h
import pickle

#from GA_RF import GA_StimTrajectory

pre_cell_types = ["excitation", "inhibition"]
pre_RF_components = ["center", "surround"]

# Set model params
def GA_CreateParams(global_params, cellPos, model):
    model.input_params = {                  # Parameters that set the inputs the simulation will act upon
        "RF_params": {
            cell_type: {    # Excitation/Inhibition
                **{
                    cs: {       # Center/Surround
                        'peak': np.random.uniform(0.01, 1, (global_params['num_input_types'], 1)),      # Strength of the component
                        # Spatial RF
                        'widthX': np.random.uniform(10, 200, (global_params['num_input_types'], 1)),    # Width (x)
                        'widthY': np.random.uniform(10, 200, (global_params['num_input_types'], 1)),    # Width (y)
                        'widthC': np.random.uniform(0, 1.7, (global_params['num_input_types'], 1)),                      # RF rotation
                        # Temporal RF
                        'tau': np.random.uniform(10, 100, (global_params['num_input_types'], 1)),      # Activation tau
                        'tauRRP': np.random.uniform(10, 1000, (global_params['num_input_types'], 1)),   # Inactivation tau
                    }
                    for cs in pre_RF_components
                },
                'depression': np.random.uniform(0.0, 1, (global_params['num_input_types'], 1)),
                'depression_tau': np.random.uniform(10, 10, (global_params['num_input_types'], 1)),
                'facilitation': np.random.uniform(0.0, 1, (global_params['num_input_types'], 1)),
                'facilitation_tau': np.random.uniform(10, 10, (global_params['num_input_types'], 1)),
            }
            for cell_type in pre_cell_types
        },

        "passive_params" : {                # Passive parameters
            'pas': np.random.uniform(1e-5, 1e-4),
            'Ra': np.random.uniform(50, 300),
        },

        "trajectory" :      [],     # Stimulus to use
        "zero_trajectory" : [],     # Stimulus before rotation

        "input_exc_syn_time" :     [],
        "input_inh_syn_time" :     [],

        "input_exc_syn_gain": np.random.uniform(0, 1),       # Strength of input signals (excitatory)
        "input_inh_syn_gain": np.random.uniform(1, 10),       # Strength of input signals (inhibitory)
        "fromPhsynG": np.random.uniform(0, 1),       # Strength of BC signals
        "gen": 0,
    }

    # Load parametrs
    if not global_params['randomStart'] :
        with open(f"params/input_params_0.pkl", "rb") as f:
            model.input_params = pickle.load(f)

    # Output of the detector
    model.output_params = {                     # Return info - voltage and calcium
        "cell" : {
            'somaV':                        [],         # Somatic voltage
            'somaV_all_angles':             [],         # Somatic voltage for all probed angles
            'max_soma_single_run':          [],

            'dsi':                      -1,           #computed DSI index
            'adjusted_dsi':             -1,           #DSI index adjusted by the peak of the pref response
            'dsi_list':                 [],
            'adjusted_dsi_list':        [],
            
            'morph_x':                  [],
            'morph_y':                  [],

        } ,           
        "input_exc_cell" : [],
        "input_inh_cell" : [],     
        'score' : [],
    }