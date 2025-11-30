import numpy as np
from neuron import h
import math 

from GA_Params import GA_CreateParams
from GA_BC import GA_Input_SetSynapses
pre_cell_types = ["excitation", "inhibition"]
pre_RF_components = ["center", "surround"]

# Main module describing the cell of interestes (typically SAC / DSGCs)
class Model:
    pass
    def __init__(self):
        self.cell = []              # Morhphology etc
        self.input_params = {}
        self.output_params = {}
        self.name = []
        # input
        self.input_exc_list= []
        self.input_inh_list= []
 
        # Internal
        self.max_input_extent= []    # Size of the InputDends (from soma or along the x axis)
        self.min_input_extent= []    # Size of the InputDends (from soma or along the x axis)
        self.APcounter= []

# Prepare and populate recording vectors
def GA_RecordingVectors(global_params, model, prep= False, populate= False):
    if prep:
        # Prepare the recording waves
        input_time= ['input_exc_syn_time', 'input_inh_syn_time']
        for ei in input_time:
            for cls in model.input_params[ei]:
                cls['drive']= []
                cls['original']= []
                cls['activation']= []
                cls['inactivation']= []
                cls['original_s']= []
                cls['activation_s']= []
                cls['inactivation_s']= []        

        for bc in range(len(model.input_exc_list)):
            model.output_params['input_exc_cell'][bc]['somaV_all_angles']= []
            model.output_params['input_exc_cell'][bc]['cell_input_all_angles']= []
            model.output_params['input_exc_cell'][bc]['cell_output_all_angles']= []
        for ac in range(len(model.input_inh_list)):
            model.output_params['input_inh_cell'][ac]['somaV_all_angles']= []
            model.output_params['input_inh_cell'][ac]['cell_input_all_angles']= []
            model.output_params['input_inh_cell'][ac]['cell_output_all_angles']= []

        model.output_params['cell']['somaV_all_angles']= []
        model.output_params['trajectory']= []
        model.output_params['cell']['dsi_list']= []
        model.output_params['cell']['adjusted_dsi_list']= []


    if populate:
        if global_params['save_all']:
            model.output_params['trajectory'].append(model.input_params['trajectory'])
        # Save the input responses
        exc_drive= []
        inh_drive= []
        for i, bc in enumerate(model.input_exc_list):
            model.output_params['input_exc_cell'][i]['somaV_all_angles'].append(model.output_params['input_exc_cell'][i]['somaV'].to_python())
            model.output_params['input_exc_cell'][i]['cell_output_all_angles'].append(bc.shifted_drive_vector)
            model.output_params['input_exc_cell'][i]['cell_input_all_angles'].append(bc.input_vector)
            exc_drive.append(bc.shifted_drive_vector)

        for i, ac in enumerate(model.input_inh_list):
            model.output_params['input_inh_cell'][i]['somaV_all_angles'].append(model.output_params['input_inh_cell'][i]['somaV'].to_python())
            model.output_params['input_inh_cell'][i]['cell_output_all_angles'].append(ac.shifted_drive_vector)        
            model.output_params['input_inh_cell'][i]['cell_input_all_angles'].append(ac.input_vector)      
            inh_drive.append(ac.shifted_drive_vector)  
        
        # Save the somatic responses

            
        model.output_params['cell']['somaV_all_angles'].append(model.output_params['cell']['somaV'].to_python())
        analysis_time= int(model.output_params['cell']['somaV'].size() * global_params['analysis_time'])

        sum_vec = np.sum(exc_drive, axis=0)
        model.output_params['cell']['max_conductance_single_run'].append(np.max(sum_vec[analysis_time : ]) - np.min(sum_vec[analysis_time : ]))
            
        if(model.input_params['gen'] >= global_params['switch_to_spiking_model']):
            model.output_params['cell']['max_soma_single_run'].append(model.cell.APcounter.size())
        else:
            # dV response, penalty for depolarized potentials
            max_response= np.max(model.output_params['cell']['somaV'].to_python()[analysis_time : ] )-np.min(model.output_params['cell']['somaV'].to_python()[analysis_time : ])
            max_response*= 1-math.tanh((60 + np.min(model.output_params['cell']['somaV'].to_python()[analysis_time : ]))/60) #np.clip(1 - ((60 + np.min(model.output_params['cell'][cell]['somaV'].to_python()[analysis_time : ]))/60), 0, 1)
            model.output_params['cell']['max_soma_single_run'].append(max_response)


# Create cells and synaptic connections
def NEURON_SetCells(global_params, model):
    cellPos= [(0,0)]

    GA_CreateParams(global_params, cellPos, model)


    # Create the detector in NEURON       
    h(f"objref Cell") 
    if(global_params['cellType'] == 'RGC'):					# RGC	       
        h('load_file("RGCmodel.hoc")')
        h("Cell= new Ganglion()")       
    else:                  
        h('load_file("single_compartment.hoc")')
        h("Cell= new s_comp()")	

    model.cell=h.Cell
    # Save morphology
    max_dist= 0
    for sec in model.cell.all:
        sec.push()   
        for i in range(int(h.n3d())):      
            model.output_params['cell']['morph_x'].append(h.x3d(i))
            model.output_params['cell']['morph_y'].append(h.y3d(i))
        model.output_params['cell']['morph_x'].append(np.nan)
        model.output_params['cell']['morph_y'].append(np.nan)  
        max_dist= max(max_dist, h.distance(model.cell.soma(0.5), sec(1)))
        h.pop_section()  

    # Add passive params
    h('forall insert pas')  
    h('forall e_pas= -60')  
    h('forall nseg= 1')  

    # Connect presynaptic inputs
    GA_Input_SetSynapses(global_params, model)
    
    input_population= [model.input_exc_list]
    if global_params['pre_inhibition']:
        input_population.append(model.input_inh_list) 

    for ei in range(len(input_population)):
        input_cell= 'input_exc_cell'
        if ei == 1: # inhibition
            input_cell= 'input_inh_cell'
        for bc in input_population[ei]:
            post_x= []
            post_y= []
            for syn in bc.output:
                post_x.append(syn.post_x)
                post_y.append(syn.post_y)
            model.output_params[input_cell].append({
                'somaV':                    [],    # Voltage 
                'somaV_all_angles':         [],
                'cell_input_all_angles':    [],
                'cell_output_all_angles':    [],
                'pre_x':                    bc.x,
                'pre_y':                    bc.y,
                'cluster':                  bc.cluster,
                'post_x':                   post_x,
                'post_y':                   post_y,            
            })


    # Record signals
    def NEURON_RecordSignals(global_params, model):
        # Bipolar cells
        for bc in range(len(model.input_exc_list)):
            model.output_params['input_exc_cell'][bc]['somaV'] = h.Vector().record(model.input_exc_list[bc].soma(0.5)._ref_v, global_params['stim_params']['dt'])
        for ac in range(len(model.input_inh_list)):
            model.output_params['input_inh_cell'][ac]['somaV'] = h.Vector().record(model.input_inh_list[ac].soma(0.5)._ref_v, global_params['stim_params']['dt'])

        # Record somatic signals
        if model.input_params['gen'] >= global_params['switch_to_spiking_model']:
            model.output_params['cell']['somaV'] = h.Vector().record(model.cell.soma(0.5)._ref_v, .1)
        else:
            model.output_params['cell']['somaV'] = h.Vector().record(model.cell.soma(0.5)._ref_v, global_params['stim_params']['dt'])

    
    # Record signals
    NEURON_RecordSignals(global_params, model)
