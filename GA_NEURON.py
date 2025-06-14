import numpy as np
from neuron import h
from GA_RF import GA_Pre_activation, GA_StimTrajectory, quick_plot
#from h5 import GA_LoadInputsFromH5
import pickle
import copy
pre_cell_types = ["excitation", "inhibition"]
pre_RF_components = ["center", "surround"]

class Input_Synapse:
    def __init__(self, x, y, cls, drive_vector, sec, cell_num, e=0, dt=10, counter=0):
        self.x = x
        self.y = y
        self.cls= cls
        self.drive_vector= drive_vector    # assumed to be a numpy array
        self.syn = h.SynVec(sec(0.5))  # a NEURON point process instance
        self.cell = cell_num
        self.syn.e = e
        self.counter= counter
        self.shifted_drive_vector= h.Vector().from_python(drive_vector)        
        self.shifted_drive_vector.play(self.syn._ref_g, dt)

class Feedback_Synapse:
    def __init__(self, preX, preY, postX, postY, post_sec, pre_sec, postcell, precell, e=-60):
        self.preX = preX
        self.preY = preY
        self.postX = postX
        self.postY = postY

        h(f'objref TempGi,TempSyn')
        h(f'TempGi= new Vector()')
        post_sec.push()                       
        h(f'TempSyn= new SynPointer(.5)')
        h.pop_section()  
        pre_sec.push()
        h(f"setpointer TempSyn.g, cai(.5)")
        h.pop_section()  
        h(f'TempGi.record(&TempSyn.g, 1)')

        self.syn = h.TempSyn    #SynPointer(post_sec(0.5))  # a NEURON point process instance
        self.precell = precell
        self.postcell = postcell
        self.syn.e = e
        self.presynapticCa= h.Vector().record(pre_sec(.5)._ref_cai, 1) 
        self.conductance= h.TempGi       
        h(f"setpointer {self.syn}.g, {pre_sec}.cai(.5)")
    
# Prepare and populate recording vectors
def GA_RecordingVectors(global_params, stim_params, model, prep= False, populate= False):
    if prep:
        # Prepare the recording waves
        for cell in range(global_params['numCells']):
            model.output_params['cell'][cell]['somaV_all_angles']= []
            model.output_params['cell'][cell]['conductance']['ca_total']= []
            model.output_params['cell'][cell]['conductance']['k_total']= []

            model.output_params['cell'][cell]['dendrite_vectors']['ca_right_all_angles']= []
            model.output_params['cell'][cell]['dendrite_vectors']['ca_left_all_angles']= []
            model.output_params['cell'][cell]['dendrite_vectors']['v_right_all_angles']= []
            model.output_params['cell'][cell]['dendrite_vectors']['v_left_all_angles']= []
            model.output_params['cell'][cell]['dendrite_vectors']['max_right']= []
            model.output_params['cell'][cell]['dendrite_vectors']['max_left']= []
            model.output_params['cell'][cell]['dendrite_vectors']['speed']= []
            model.output_params['cell'][cell]['dendrite_vectors']['contrast']= []
            model.output_params['cell'][cell]['dendrite_vectors']['angle']= []

            model.output_params['cell'][cell]['peak_ca']['peak']= []
            model.output_params['cell'][cell]['peak_ca']['x']= []
            model.output_params['cell'][cell]['peak_ca']['y']= []
            model.output_params['cell'][cell]['peak_ca']['speed']= []
            model.output_params['cell'][cell]['peak_ca']['contrast']= []
            model.output_params['cell'][cell]['peak_ca']['angle']= []


            model.output_params['cell'][cell]['dsi_list']= []

        model.output_params['InputE_syn_time_copy']= []
        model.output_params['InputI_syn_time_copy']= []
        model.output_params['I_syn_g']= []
        model.output_params['I_syn_data']['postX']= []
        model.output_params['I_syn_data']['postY']= []
        model.output_params['I_syn_data']['precell']= []
        model.output_params['I_syn_data']['postcell']= []
        model.output_params['I_syn_data']['speed']= []
        model.output_params['I_syn_data']['contrast']= []
        model.output_params['I_syn_data']['angle']= []


    if populate:

        # Save the somatic responses
        for cell in range(global_params['numCells']):
            model.output_params['cell'][cell]['dendrite_vectors']['speed'].append(stim_params['speed'])
            model.output_params['cell'][cell]['dendrite_vectors']['contrast'].append(stim_params['contrast'])
            model.output_params['cell'][cell]['dendrite_vectors']['angle'].append(stim_params['angle'])
            
            model.output_params['cell'][cell]['somaV_all_angles'].append(model.output_params['cell'][cell]['somaV'].to_python())
            model.output_params['cell'][cell]['max_soma_single_run'].append(np.max(model.output_params['cell'][cell]['somaV'].to_python())-np.min(model.output_params['cell'][cell]['somaV'].to_python()))
            # Compute the mean calcium/voltage responses
            if(len(model.output_params['cell'][cell]['dendrite_vectors']['ca_right']) > 0):
                # Calcium - all dendrites
                for i, ca in enumerate( model.output_params['cell'][cell]['peak_ca']['ca_all'] ):                        
                    model.output_params['cell'][cell]['peak_ca']['peak'].append(np.max(ca.to_python()))
                    model.output_params['cell'][cell]['peak_ca']['x'].append(model.output_params['cell'][cell]['peak_ca']['pos_rec_all_x'][i])
                    model.output_params['cell'][cell]['peak_ca']['y'].append(model.output_params['cell'][cell]['peak_ca']['pos_rec_all_y'][i])
                    model.output_params['cell'][cell]['peak_ca']['speed'].append(stim_params['speed'])
                    model.output_params['cell'][cell]['peak_ca']['contrast'].append(stim_params['contrast'])
                    model.output_params['cell'][cell]['peak_ca']['angle'].append(stim_params['angle'])
                # Calcium - important dendrites
                response_right= np.mean(np.stack(model.output_params['cell'][cell]['dendrite_vectors']['ca_right']), axis=0)
                response_left= np.mean(np.stack(model.output_params['cell'][cell]['dendrite_vectors']['ca_left']), axis=0)
                # Conductances
                if(len(model.output_params['cell'][cell]['conductance']['ca_all'])>0):
                    model.output_params['cell'][cell]['conductance']['ca_total'].append(np.sum( model.output_params['cell'][cell]['conductance']['ca_all'], axis=0))
                if(len(model.output_params['cell'][cell]['conductance']['k_all'])>0):
                    model.output_params['cell'][cell]['conductance']['k_total'].append(np.sum( model.output_params['cell'][cell]['conductance']['k_all'], axis=0))

                # Save params
                #response_right[0]= stim_params['speed']
                #response_right[1]= stim_params['contrast']
                #response_right[2]= stim_params['angle']
                
                model.output_params['cell'][cell]['dendrite_vectors']['ca_right_all_angles'].append(response_right)
                model.output_params['cell'][cell]['dendrite_vectors']['ca_left_all_angles'].append(response_left)
                model.output_params['cell'][cell]['dendrite_vectors']['max_right'].append(np.max(response_right))
                model.output_params['cell'][cell]['dendrite_vectors']['max_left'].append(np.max(response_left))
                model.output_params['cell'][cell]['dendrite_vectors']['max_right_single_run'].append(np.max(response_right))
                model.output_params['cell'][cell]['dendrite_vectors']['max_left_single_run'].append(np.max(response_left))
                # Voltage
                model.output_params['cell'][cell]['dendrite_vectors']['v_right_all_angles'].append(np.mean(np.stack(model.output_params['cell'][cell]['dendrite_vectors']['v_right']), axis=0))
                model.output_params['cell'][cell]['dendrite_vectors']['v_left_all_angles'].append(np.mean(np.stack(model.output_params['cell'][cell]['dendrite_vectors']['v_left']), axis=0))
            
        # Synapses
        for syn in model.SAC_SAC_synapses:
            #model.output_params['I_syn_time_copy'].append(syn.presynapticCa.c())
            model.output_params['I_syn_g'].append(syn.conductance.c())
            model.output_params['I_syn_data']['postX'].append(syn.postX)
            model.output_params['I_syn_data']['postY'].append(syn.postY)
            model.output_params['I_syn_data']['precell'].append(syn.precell)
            model.output_params['I_syn_data']['postcell'].append(syn.postcell)
            model.output_params['I_syn_data']['speed'].append(stim_params['speed'])
            model.output_params['I_syn_data']['contrast'].append(stim_params['contrast'])
            model.output_params['I_syn_data']['angle'].append(stim_params['angle'])
        for syn in model.InputE_synapses:             
            model.output_params['InputE_syn_time_copy'].append(syn.shifted_drive_vector.to_python()) 
        for syn in model.InputI_synapses:             
            model.output_params['InputI_syn_time_copy'].append(syn.shifted_drive_vector.to_python()) 

# Introduce mutations and check boundaries
def NEURON_mutation(global_params, input_params):
    np.random.seed()
    # RF params
    for cell_type in pre_cell_types:    # Excitation/Inhibition
        for cs in pre_RF_components:    # Center/Surround
            for key, value in input_params['RF_params'][cell_type][cs].items():
                input_params['RF_params'][cell_type][cs][key]= value * np.random.normal(1, global_params['mutationRate'], (global_params['numPreClus'],1)) + np.random.uniform(0, global_params['mutationRate'] / 10000, (global_params['numPreClus'], 1))
                if(key in ['widthX', 'widthY']):
                    input_params['RF_params'][cell_type][cs][key]= np.clip(input_params['RF_params'][cell_type][cs][key], 10, 200 + 2000 * (cs == 'surround'))

                if(key in ['tau', 'tauRRP']):
                    input_params['RF_params'][cell_type][cs][key]= np.clip(input_params['RF_params'][cell_type][cs][key], global_params['dt'], 10000) 
                if((key == 'peak') and (cs == 'surround')):
                    input_params['RF_params'][cell_type][cs][key]= np.clip(input_params['RF_params'][cell_type][cs][key], 0, 1) 
          
            # modify the params based on constrains
            if(global_params['RF_constrains'][cell_type][cs]['sameKinetics'] ):        # similar kinetics
                for key in ['tau', 'tauRRP']:
                    input_params['RF_params'][cell_type][cs][key].fill(input_params['RF_params'][cell_type][cs][key][0]) 
            if(global_params['RF_constrains'][cell_type][cs]['sameSize'] ):             # similar RF size
                for key in ['widthX', 'widthY']:
                    input_params['RF_params'][cell_type][cs][key].fill(input_params['RF_params'][cell_type][cs][key][0]) 
            if(global_params['RF_constrains'][cell_type][cs]['sameOrientation'] ):             # similar RF orientation
                input_params['RF_params'][cell_type][cs]['widthC'].fill(input_params['RF_params'][cell_type][cs]['widthC'][0]) 
                for cls in range(global_params['numPreClus']):
                    input_params['RF_params'][cell_type][cs]['widthX'][cls]=input_params['RF_params'][cell_type][cs]['widthY'][cls]  # Circle and not oval
                
            if(global_params['RF_constrains'][cell_type][cs]['sameAmplitude'] ):                # similar RF strength
                input_params['RF_params'][cell_type][cs]['peak'].fill(input_params['RF_params'][cell_type][cs]['peak'][0]) 
                   

        if not global_params['RF_constrains'][cell_type]['surround']['doSurround']:        # No surround
            input_params['RF_params'][cell_type]['surround']['peak']*= 0 
            
    # Passive params
    for key, value in input_params['passive_params'].items():
        input_params['passive_params'][key]= value * np.random.normal(1, global_params['mutationRate']) + np.random.uniform(0, global_params['mutationRate'] / 10000)
        if(key == 'Ra'):
           input_params['passive_params'][key]= np.clip(input_params['passive_params'][key], 50, 300) 
    
    # Active params
    for key, value in input_params['active_params'].items():
        if key in ['shift', 'mN_caGA', 'nN_kGa' ]:  # integer
            input_params['active_params'][key]= value + np.random.normal(0, global_params['mutationRate'] / 10)
            if key in ['mN_caGA', 'nN_kGa']:
                input_params['active_params'][key]= np.clip(input_params['active_params'][key], 1, 8) # Nonlinear gating (multiple subunits)
        else:
            if key in ['minf_caGA', 'ninf_kGA', 'hinf_caGA', 'hinf_kGA']:
                input_params['active_params'][key]= np.clip(input_params['active_params'][key], 0, 1) # Gating states
                input_params['active_params'][key]= value * np.random.normal(1, global_params['mutationRate'], (input_params['active_params'][key].shape[0], 1)) #+ np.random.uniform(0, global_params['mutationRate'] / 10000, (global_params['numPreClus'], 1))
            elif key in ['mtau_caGA', 'ntau_kGA', 'htau_caGA', 'htau_kGA']:
                input_params['active_params'][key]= np.clip(input_params['active_params'][key], 0.1, 100000) # Activation/Inactivation 
                input_params['active_params'][key]= value * np.random.normal(1, global_params['mutationRate'], (input_params['active_params'][key].shape[0], 1)) #+ np.random.uniform(0, global_params['mutationRate'] / 10000, (global_params['numPreClus'], 1))
            else:
                input_params['active_params'][key]= value * np.random.normal(1, global_params['mutationRate'], (global_params['numPostComp'], 1)) #+ np.random.uniform(0, global_params['mutationRate'] / 10000, (global_params['numPreClus'], 1))

    # SAC - SAC syn          
    input_params["SAC_SACsynG"]= input_params["SAC_SACsynG"] * np.random.normal(1, global_params['mutationRate'] ) + np.random.uniform(0, global_params['mutationRate'] / 10000)

# Set model params
def NEURON_CreateParams(global_params, stim_params, cellPos, model):
    trajectory = GA_StimTrajectory(stim_params)

    model.input_params = {                  # Parameters that set the inputs the simulation will act upon
        "RF_params": {
            cell_type: {    # Excitation/Inhibition
                cs: {       # Center/Surround
                    'peak': np.random.uniform(0., 0.1, (global_params['numPreClus'], 1)),      # Strength of the component
                    # Spatial RF
                    'widthX': np.random.uniform(10, 200, (global_params['numPreClus'], 1)),    # Width (x)
                    'widthY': np.random.uniform(10, 200, (global_params['numPreClus'], 1)),    # Width (y)
                    'widthC': np.zeros((global_params['numPreClus'], 1)),                      # RF rotation
                    # Temporal RF
                    'tau': np.random.uniform(10, 1000, (global_params['numPreClus'], 1)),      # Activation tau
                    'tauRRP': np.random.uniform(10, 1000, (global_params['numPreClus'], 1)),   # Inactivation tau
                    'del': np.zeros((global_params['numPreClus'], 1)),                         # Delay from visual stimulus
                }
                for cs in pre_RF_components
            }
            for cell_type in pre_cell_types
        },

        "passive_params" : {                # Passive parameters
            'pas': np.random.uniform(1e-5, 1e-4),
            'Ra': np.random.uniform(50, 300),
        },

        "active_params" : {
            # N type calcium
            'gbar_canrgc' : np.random.uniform(0.0001, .0005, (global_params['numPostComp'], 1)),
            'shift_canrgc' : np.random.uniform(-1, 1),
            # L type calcium
            'gbar_calrgc' : np.random.uniform(0.0001, .0005, (global_params['numPostComp'], 1)),
            'shift_calrgc' : np.random.uniform(-1, 1),
            # genetic algorithm trainable 
            'gMax_caGA' : np.random.uniform(0.00001, .0001, (global_params['numPostComp'], 1)),
            'mN_caGA' : np.random.uniform(1, 3),
            'minf_caGA' : np.random.uniform(0.0001, 1, (global_params['numVoltagePoints'], 1)),
            'hinf_caGA' : np.random.uniform(0.0001, 1, (global_params['numVoltagePoints'], 1)),
            'mtau_caGA' : np.random.uniform(0.1, 10, (global_params['numVoltagePoints'], 1)),
            'htau_caGA' : np.random.uniform(0.1, 100, (global_params['numVoltagePoints'], 1)),


            # "Linear" calcium"
            'gbar_caLinear' : np.random.uniform(0.0001, .001, (global_params['numPostComp'], 1)),

            # A type potassium
            'gbar_ka' : np.random.uniform(0.001, .01, (global_params['numPostComp'], 1)),
            'shift_ka' : np.random.uniform(-1, 1),                      
            # M type potassium
            'gbar_km' : np.random.uniform(0.1, 1, (global_params['numPostComp'], 1)),
            'shift_km' : np.random.uniform(-1, 1),                      
            # DR slow type potassium
            'gbar_kd' : np.random.uniform(0.10, 0.1, (global_params['numPostComp'], 1)),
            'shift_kd' : np.random.uniform(-1, 1),                      
            # IH
            'gbar_ih' : np.random.uniform(0.001, 0.01, (global_params['numPostComp'], 1)),
            'shift_ih' : np.random.uniform(-1, 1),                      

            # genetic algorithm trainable 
            'gMax_kGA' : np.random.uniform(0.0001, .001, (global_params['numPostComp'], 1)),
            'nN_kGA' : np.random.uniform(1, 3),
            'ninf_kGA' : np.random.uniform(0.0001, 1, (global_params['numVoltagePoints'], 1)),
            'hinf_kGA' : np.random.uniform(0.0001, 1, (global_params['numVoltagePoints'], 1)),
            'ntau_kGA' : np.random.uniform(0.0001, 1, (global_params['numVoltagePoints'], 1)),
            'htau_kGA' : np.random.uniform(0.0001, 1, (global_params['numVoltagePoints'], 1)),

        },

        "trajectory" : trajectory,

        "InputE_syn_time" : np.zeros((global_params['numPreClus'], int(stim_params['tStop'] / stim_params['dt']))),
        "InputI_syn_time" : np.zeros((global_params['numPreClus'], int(stim_params['tStop'] / stim_params['dt']))),

        "SAC_SACsynG": np.random.uniform(00, 100),
        
    }
    # Load parametrs
    if(global_params['randomStart'] == False):
        with open(f"params/input_params_{global_params['job_id']}.pkl", "rb") as f:
            model.input_params = pickle.load(f)
    
    model.input_params['cellPos']= cellPos
    model.input_params['cell_syn_loc'] = [{
        'x': [],
        'y': [],
        'section': [],
        'offset': [],
    } for _ in range(global_params['numCells'])]

    model.output_params = {                     # Return info - voltage and calcium
        "cell" : [{
            'somaV':                        h.Vector(), # Somatic voltage
            'somaV_all_angles':             [],          # Somatic voltage for all probed angles
            'max_soma_single_run':          [],
            'conductance':                  {
                'ca_total':                    [],  
                'k_total':                     [],
                'ca_all':                      [],
                'k_all':                       [],
            },

            'dendrite_vectors':{
                'v_right':                 [],
                'v_right_all_angles':      [],
                'v_left':                  [],
                'v_left_all_angles':       [],

                'ca_right':                [],          # SAC - calcium signals on the right side of the cell
                'ca_right_all_angles':     [],          # SAC - MEAN calcium signals on the right side of the cell for different direction
                'ca_left':                 [],
                'ca_left_all_angles':      [], 

                'pos_rec_right_x' :        [],
                'pos_rec_right_y' :        [],
                'pos_rec_left_x' :         [],
                'pos_rec_left_y' :         [],

                'max_right':               [],
                'max_left':                [],

                #'score_right':               [],
                #'score_left':                [],

                'speed':                   [],
                'contrast':                [],
                'angle':                   [],  

                'max_right_single_run':     [],
                'max_left_single_run':     [],
            },

            'dsi':                     0,           #computed DSI index
            'dsi_list':                 [],
            #'speed':                   [],
            #'contrast':                [],
            #'angle':                   [],              
            'morph_x':                  [],
            'morph_y':                  [],
            'activeCompartment':        {
                'x':                        [],
                'y':                        [],
                'type':                     [],
            },

            'peak_ca':{
                'ca_all' :                  [],
                'pos_rec_all_x' :           [],
                'pos_rec_all_y' :           [],
                'peak':                     [],
                'x':                        [],
                'y':                        [],
                'speed':                    [],
                'contrast':                 [],
                'angle':                    [],                  
            }

        }for _ in range(global_params['numCells'])]  ,      
        # Synaptic inputs: Excitatory
        'InputE_syn_time_copy' :        [],  #Temporal responses #
        'InputI_syn_time_copy' :        [],  #Temporal responses 
        # Synaptic inputs: Inhibitory
        #'I_syn_time_copy' :        [],  # Temporal responses
        'I_syn_g' :        [],  # Temporal responses
        'I_syn_data'    :          {
            "postX" :       [],
            "postY" :       [],
            "postcell" :    [],
            "precell" :     [],
            'speed':        [],
            'contrast':     [],
            'angle':        [],          
        },
        'score' : [],
    }

# Create cells and synaptic connections
def NEURON_SetCells(global_params, stim_params, model):
    cellPos= [(0,0)]
    #sacPosY= [0]
    global_params['numCells']= 1
    if(global_params['numSAClayersY'] < 1):
        global_params['numSAClayersY']= 1

    if(global_params['cellType'] == 'SAC network'):	# SAC network
        for layerX in range(global_params['numSAClayersX']):
            for layerY in range(global_params['numSAClayersY']):
                pos= ((layerX - np.floor(global_params['numSAClayersX'] / 2) + 0.5 * abs((layerY + (global_params['numSAClayersY'] - 1) / 2) % 2)) * global_params['distSAC'] , (layerY - np.floor(global_params['numSAClayersY'] / 2 )) * global_params['distSAC'] * np.cos(60/180*np.pi))
                if pos not in cellPos:
                    cellPos.append(pos)

 
        global_params['numCells'] = len(cellPos)		#---number of cells
    #print(len(cellPos), cellPos)
    # Create the 'models' container
    NEURON_CreateParams(global_params, stim_params, cellPos, model)

    """
    # Manually load saved model params
    if(False):
        # Use specific params 
        model.input_params['passive_params']['Ra']= 50
        model.input_params['passive_params']['pas']= 8e-6
        model.input_params['active_params']['gbar_calrgc'][0]=0.0001
        model.input_params['active_params']['gbar_calrgc'][1]=0.0001882045831671943
        model.input_params['active_params']['gbar_calrgc'][2]=0.0003609851401634234
        model.input_params['active_params']['gbar_calrgc'][3]=0.002841117527129199
        model.input_params['active_params']['gbar_ka'][0]=0.005404519983427247
        model.input_params['active_params']['gbar_ka'][1]=0.001692975083387725
        model.input_params['active_params']['gbar_ka'][2]=0.0006458936041504057
        model.input_params['active_params']['gbar_ka'][3]=0.00251459867737246
        model.input_params['SAC_SACsynG']=16

    if(global_params['randomStart'] == False):
        with open("input_params.pkl", "rb") as f:
            model.input_params = pickle.load(f)
        model.input_params['cellPos']= cellPos
    """    
    # Create the cells in NEURON       
    h(f"objref Cell[{global_params['numCells']}]") 

    if(global_params['cellType'] == 'RGC'):					# RGC	
        
        h('load_file("RGCmodel.hoc")')
        h("Cell= new Ganglion()")       
    else:
        if(global_params['cellType'] == 'L23'):				# L23	
           
            h('load_file("layer_2_3.hoc")')
            h("Cell= new CorticalL23()")
        else:
            if(global_params['cellType'] == 'L5'):           # L5	
                
                h('load_file("layerV1.hoc")')
                h("Cell= new CorticalL5(0,0)")	
            else:
                if(global_params['cellType'] == 'SAC'):       # SAC
                    
                    h('load_file("SAC1.hoc")')
                    h("Cell= new star(0,0)")	    
                else:
                    if(global_params['cellType'] == 'SAC network'):   # SAC network
                        
                        h('load_file("SACnet.hoc")')
                        for cell in range(global_params['numCells']):
                            h(f"Cell[{cell}] = new StarNet({cellPos[cell][0]},{cellPos[cell][1]})")                       
                    else:               # single compartment
                        
                        h('load_file("single_compartment.hoc")')
                        h("Cell= new s_comp()")	
  
    for cell in range(global_params['numCells']): 
        model.cell.append(h.Cell[cell]) 
        # Save morphology
        #model.cell[cell].soma.push()
        #h('distance()')
        #h.pop_section() 
        max_dist= 0
        for sec in model.cell[cell].all:
            sec.push()   
            for i in range(int(h.n3d())):      
                model.output_params['cell'][cell]['morph_x'].append(h.x3d(i))
                model.output_params['cell'][cell]['morph_y'].append(h.y3d(i))
            model.output_params['cell'][cell]['morph_x'].append(np.nan)
            model.output_params['cell'][cell]['morph_y'].append(np.nan)  
            max_dist= max(max_dist, h.distance(model.cell[cell].soma(0.5), sec(1)))
            h.pop_section()  
        # Find the positions of the active compartments (based on distance from soma or along the null-preferred axis)
        for d, sec in enumerate(model.cell[cell].all):
            sec.push()    
            model.output_params['cell'][cell]['activeCompartment']['x'].append(h.x3d(int(h.n3d()/2)))
            model.output_params['cell'][cell]['activeCompartment']['y'].append(h.y3d(int(h.n3d()/2)))
            type= min(global_params['numPostComp'] - 1, int(h.distance(model.cell[cell].soma(0.5), sec(1)) / max_dist * global_params['numPostComp']))
            model.output_params['cell'][cell]['activeCompartment']['type'].append(type)
            #print(int(h.distance(model.cell[cell].soma(0.5), sec(1)) / max_dist * global_params['numPostComp']), h.distance(model.cell[cell].soma(0.5), sec(1)))
            h.pop_section()             
    
    # Add passive params
    h('forall insert pas')  
    h('forall e_pas= -60')  
   
    # Add active params
    if 'CaN' in global_params['activeChannels']:            # N-type
       h('forall insert canrgc')  
       global_params['ca_present']= True
    else:
       model.input_params["active_params"]["gbar_canrgc"]= 0 
    
    if 'CaL' in global_params['activeChannels']:            # L-type
       h('forall insert calrgcfix')  
       global_params['ca_present']= True
    else:
       model.input_params["active_params"]["gbar_calrgc"]= 0 

    if 'caGA' in global_params['activeChannels']:            #Flexible ca++
       h('forall insert caGA')  
       global_params['ca_present']= True
    else:
       model.input_params["active_params"]["gMax_caGA"]= 0 

    if 'caLinear' in global_params['activeChannels']:            # Linear-type
       h('forall insert caLinear')  
       global_params['ca_present']= True
    else:
       model.input_params["active_params"]["gbar_caLinear"]= 0 

    if global_params['ca_present']:                         # Diffusion
       h('forall insert cadiff')

       

    if 'Ka' in global_params['activeChannels']:            # A-type
       h('forall insert kap')  
       global_params['k_present']= True
    else:
       model.input_params["active_params"]["gbar_ka"]= 0 

    if 'Km' in global_params['activeChannels']:            # M-type
       h('forall insert km')  
       global_params['k_present']= True
    else:
       model.input_params["active_params"]["gbar_km"]= 0 

    if 'Kdr' in global_params['activeChannels']:            # DR-type
       h('forall insert kSlow')  
       global_params['k_present']= True
    else:
       model.input_params["active_params"]["gbar_kd"]= 0 

    if 'IH' in global_params['activeChannels']:            # DR-type
       h('forall insert ih')  
       global_params['ih_present']= True
    else:
       model.input_params["active_params"]["gbar_ih"]= 0 

    if 'kGA' in global_params['activeChannels']:            #Flexible k+
       h('forall insert kGA')  
       global_params['k_present']= True
    else:
       model.input_params["active_params"]["gMax_kGA"]= 0 

    # Record signals (mostly from central cell)
    def NEURON_RecordSignals(global_params, model):
        for cell in range(global_params['numCells']):
            # Record somatic signals
            model.output_params['cell'][cell]['somaV'] = h.Vector().record(model.cell[cell].soma(0.5)._ref_v, 1)
            # Record conductances
            for sec in model.cell[cell].all:
                # Potassium
                if 'Ka' in global_params['activeChannels']:
                    model.output_params['cell'][cell]['conductance']['k_all'].append(h.Vector().record(sec(0.5)._ref_gka_kap, 1))            
                if 'Km' in global_params['activeChannels']:
                    model.output_params['cell'][cell]['conductance']['k_all'].append(h.Vector().record(sec(0.5)._ref_gk_km, 1))            
                if 'Kdr' in global_params['activeChannels']:
                    model.output_params['cell'][cell]['conductance']['k_all'].append(h.Vector().record(sec(0.5)._ref_gk_kSlow, 1))            
                if 'IH' in global_params['activeChannels']:
                    model.output_params['cell'][cell]['conductance']['k_all'].append(h.Vector().record(sec(0.5)._ref_gq_ih, 1))            
                # Calcium
                if 'CaN' in global_params['activeChannels']:
                    model.output_params['cell'][cell]['conductance']['ca_all'].append(h.Vector().record(sec(0.5)._ref_gca_canrgc, 1))            
                if 'CaL' in global_params['activeChannels']:
                    model.output_params['cell'][cell]['conductance']['ca_all'].append(h.Vector().record(sec(0.5)._ref_gca_calrgc, 1))    
                #if 'CaLinear' in global_params['activeChannels']:
                #    model.output_params['cell'][cell]['conductance']['ca_all'].append(h.Vector().record(sec(0.5)._ref_ica_caLinear, 1))  
            if(global_params['cellType'] in ['SAC', 'SAC network']):  # Dendritic calcium
                
                for sec in model.cell[cell].OutputDends:
                    sec.push() 
                    model.output_params['cell'][cell]['peak_ca']['ca_all'].append(h.Vector().record(sec(0.5)._ref_cai, 1)) # All output dendrites
                    model.output_params['cell'][cell]['peak_ca']['pos_rec_all_x'].append((h.x3d(0)+h.x3d(h.n3d()-1))/2)
                    model.output_params['cell'][cell]['peak_ca']['pos_rec_all_y'].append((h.y3d(0)+h.y3d(h.n3d()-1))/2)
                    if abs(h.y3d(0) - model.input_params['cellPos'][cell][1]) < 20:  # Close to horizontal axis
                        
                        if h.x3d(0) - model.input_params['cellPos'][cell][0] > 90:
                            model.output_params['cell'][cell]['dendrite_vectors']['ca_right'].append(h.Vector().record(sec(0.5)._ref_cai, 1))
                            model.output_params['cell'][cell]['dendrite_vectors']['v_right'].append(h.Vector().record(sec(0.5)._ref_v, 1))
                            model.output_params['cell'][cell]['dendrite_vectors']['pos_rec_right_x'].append((h.x3d(0)+h.x3d(h.n3d()-1))/2)
                            model.output_params['cell'][cell]['dendrite_vectors']['pos_rec_right_y'].append((h.y3d(0)+h.y3d(h.n3d()-1))/2)
                            
                        if h.x3d(0) - model.input_params['cellPos'][cell][0] < -70 - 20*(global_params['cellType'] == 'SAC network'):
                            model.output_params['cell'][cell]['dendrite_vectors']['ca_left'].append(h.Vector().record(sec(0.5)._ref_cai, 1))
                            model.output_params['cell'][cell]['dendrite_vectors']['v_left'].append(h.Vector().record(sec(0.5)._ref_v, 1))
                            model.output_params['cell'][cell]['dendrite_vectors']['pos_rec_left_x'].append((h.x3d(0)+h.x3d(h.n3d()-1))/2)
                            model.output_params['cell'][cell]['dendrite_vectors']['pos_rec_left_y'].append((h.y3d(0)+h.y3d(h.n3d()-1))/2)
                    h.pop_section()           
        #print(f"#dends on left: {len(model.output_params['cell'][cell]['dendrite_vectors']['ca_left'])}, #dends of right: {len(model.output_params['cell'][cell]['dendrite_vectors']['ca_right'])}")
    
    # Distribute synaptic inputs
    def NEURON_SetSynapses(global_params, model):
        np.random.seed(global_params['job_id'])
        # Excitatory inputs (coming from bipolar cells)
        # Potentially recompute number of synapses
        x_dist= []
        y_dist= []
        sec_dist= []        

        if(global_params['distSyn'] > 0):   # Use distances on dendrite and not random positions
            for sec_num, sec in enumerate(model.cell[0].InputDends):
                sec.push()                
                numSyn= max(1, int(sec.L/global_params['distSyn']) ) # Calc # inputs if separate input by global_params['distSyn'] micron
                for _ in range(numSyn):                     # Over all Synapses
                    x_dist.append((h.x3d(0)+h.x3d(int(h.n3d()-1)))/2 )  # mid point
                    y_dist.append((h.y3d(0)+h.y3d(int(h.n3d()-1)))/2 )
                    sec_dist.append(sec_num)                            # section stored        
                    
                h.pop_section()
            global_params['numSyn']= len(sec_dist)

        # Go over all cells
        for cell in range(global_params['numCells']): 
            # Initialize the synaptic input lists
            model.input_params["cell_syn_loc"][cell]["x"]= []
            model.input_params["cell_syn_loc"][cell]["y"]= []
            model.input_params["cell_syn_loc"][cell]["section"]= []
            model.input_params["cell_syn_loc"][cell]["offset"]= []
    
            # Compute the number of positions (length of input dendrites)
            x_coords = []
            y_coords = []
            sec_coords = []
            for sec_num, sec in enumerate(model.cell[cell].InputDends):
                sec.push()
                for i in range(int(h.n3d())):   # Range of x and y positions
                    x_coords.append(h.x3d(i))
                    y_coords.append(h.y3d(i))
                    sec_coords.append(sec)
                h.pop_section()
                
            xmin, xmax = min(x_coords), max(x_coords)   # Extend of dendrites (for random positions and presynaptic poulations)
            ymin, ymax = min(y_coords), max(y_coords)
            
            for syn_i in range(global_params['numSyn']):
                if(global_params['distSyn'] > 0):       # Load the computed distribution
                    x= x_dist[syn_i] + cellPos[cell][0]
                    y= y_dist[syn_i] + cellPos[cell][1]
                    close_sec= sec_dist[syn_i]
                else:                                   # Create new random distribution
                    closest= 10000                
                    close_sec= 0
                    while (closest > 10):
                        x= np.random.uniform(xmin, xmax)
                        y= np.random.uniform(ymin, ymax)
                        dist= 10000
                        # Find closest dend
                        for post in range(len(sec_coords)):
                            dist= abs(x - x_coords[post]) + abs(y - y_coords[post])           
                            if(dist < closest):
                                closest= dist
                                close_sec= post

                model.input_params["cell_syn_loc"][cell]["x"].append(x) 
                model.input_params["cell_syn_loc"][cell]["y"].append(y) 
                model.input_params["cell_syn_loc"][cell]["section"].append(close_sec) 
                model.input_params["cell_syn_loc"][cell]["offset"].append(0)    

                # Found presynaptic type and append a synapse
                cls= int((x - xmin) / (xmax - xmin) * global_params['numPreClus'] )# Presynaptic cluster type is based on x position
                model.InputE_synapses.append(Input_Synapse(x, y , cls, model.input_params['InputE_syn_time'][cls], sec_coords[close_sec], cell, e=0, dt=global_params['dt'], counter=syn_i))   
                if(global_params['inhibition']):
                    model.InputI_synapses.append(Input_Synapse(x, y , cls, model.input_params['InputI_syn_time'][cls], sec_coords[close_sec], cell, e=-60, dt=global_params['dt'], counter=syn_i))   


            # Inhibitory SAC-SAC recurrent (feedback) connections
            if(global_params['cellType'] == 'SAC network' and global_params['numCells'] > 1):         # Connect inhibitory SAC-SAC synapses
                h('objref post_secref, pre_secref')
                for post_sec in model.cell[cell].GABAInputDends:
                    post_sec.push()

                    # Find all close postsynaptic locations
                    postX= (h.x3d(0)+h.x3d(1))/2
                    postY= (h.y3d(0)+h.y3d(1))/2          
                    #h('post_secref= new SectionRef()')
                    h.pop_section()  
                    for precell in range(global_params['numCells']):
                        if(cell != precell):
                            mindist= 10000

                            for pre_sec in model.cell[precell].OutputDends: # Find all close presynaptic locations
                                pre_sec.push()
                                preX_= (h.x3d(0)+h.x3d(1))/2
                                preY_= (h.y3d(0)+h.y3d(1))/2
                                dist= np.sqrt((preX_ - postX)**2 + (preY_- postY)**2)
                                if(dist < mindist):
                                    mindist= dist
                                    pre_sec_save= pre_sec
                                    #h('pre_secref= new SectionRef()')
                                    preX= preX_
                                    preY= preY_
                                    
                                h.pop_section()
                            
                            if(mindist < global_params['maxSAC_syndist']):
                                #h('access post_secref.sec') 		# Postsynaptic
							    # Create an inhibitory synapse
                                
                                model.SAC_SAC_synapses.append(Feedback_Synapse(preX, preY, postX, postY, post_sec, pre_sec_save, cell, precell))   
                                """
                                h('TempSyn= new SynPointer(.5)')
                                h('TempSyn.e= -60')	
                                h(f'TempSyn.postX= {postX}')
                                h(f'TempSyn.postY= {postY}')

                                h('access pre_secref.sec') 	        # Presynaptic
                                h(f'TempSyn.preX= {preX}')
                                h(f'TempSyn.preY= {preY}')
                                h(f'TempSyn.cellNum= {precell}')
                                h('setpointer TempSyn.g, cai(.5)')
                                h(f'Cell[{cell}].GABA_Syn_List.append(TempSyn)')

                                h(f'TempGi= new Vector(GA_timePnts)')
                                h(f'TempGi.record(&TempSyn.g, 1)')
                                h(f'Cell[{cell}].GABA_postG_List.append(TempGi)')
                                """
                global_params['numIsyn']= len(model.SAC_SAC_synapses)
    # Record signals, central cell
    NEURON_RecordSignals(global_params, model)
    # Distribute synapses across all cells
    NEURON_SetSynapses(global_params, model)    



# Shifts the activation of the synaptic inputs based on the speed and direction of activation
def NEURON_UpdateSynapses(global_params, stim_params, model):  
    
    tstop= 2000
    # Make sure that the responses are within the simulation period by examining the 10% indices
    offset= []
    for syn in model.InputE_synapses:
        offset.append( int(syn.x / stim_params['speed'] * np.cos(stim_params['angle']) / global_params['dt'] + syn.y / stim_params['speed'] * np.sin(stim_params['angle']) / global_params['dt']) )
    
    min_10= stim_params['delay']
    max_10= tstop
    for cls in range(global_params['numPreClus']):
        syn_time= model.input_params['InputE_syn_time'][cls, :] 
        peak_val = np.max(syn_time) * global_params['dt']
        at10_thresh = np.where(syn_time >= peak_val / 10)[0]
        if len(at10_thresh) > 0:
            min_10= min(min_10, at10_thresh[0])
        if len(at10_thresh) >= 2:
            max_10= min(max_10, at10_thresh[-1])

    min_offset= min_10 / global_params['dt'] + min(offset)
    if(min_offset < stim_params['delay'] / global_params['dt'] ):
        offset = [int(x + (stim_params['delay'] / global_params['dt'] - min_offset))  for x in offset]

    tstop= max(tstop, max_10 + 500, max(offset) * global_params['dt'] + 500 + stim_params['delay'])
    
    # Go over all input synapses (excitatory and inhibitory)
    pre_syn_pop= (model.InputE_synapses,)
    if(global_params['inhibition']):
        pre_syn_pop= (model.InputE_synapses, model.InputI_synapses)
    
    for type, synpop in enumerate(pre_syn_pop):
        for i, syn in enumerate(synpop):
            if(type==0):
                syn.drive_vector= model.input_params['InputE_syn_time'][syn.cls, :]
            else:
                syn.drive_vector= model.input_params['InputI_syn_time'][syn.cls, :]

            vec_size= int(tstop / global_params['dt'])
            syn.shifted_drive_vector.resize(vec_size)
            syn.shifted_drive_vector.fill(0)
            syn.offset= offset[i]
            model.input_params["cell_syn_loc"][syn.cell]["offset"][syn.counter]= offset[i]
            for tt in range(len(syn.drive_vector)):
                if (offset[i] + tt < vec_size) and (offset[i] + tt >= 0):
                    syn.shifted_drive_vector.x[offset[i] + tt]= syn.drive_vector[tt] * stim_params['contrast']

    h(f"global_gain_SynPointer= {model.input_params['SAC_SACsynG']}")
    return tstop
