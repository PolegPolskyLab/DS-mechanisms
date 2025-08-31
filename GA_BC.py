import numpy as np
from neuron import h
import math 
""" 
    Bipolar cell models
    Get input from photoreceptors
"""
# Input synapses
class Input_Synapse:
    def __init__(self, bc, post_cell, post_sec, post_x, post_y, dend_pos, global_params, e=0):
        self.post_x = post_x
        self.post_y = post_y
        self.post_cell = post_cell
        post_sec.push() 
        h(f'objref bcSyn')
        h(f'bcSyn= new SynVec({dend_pos})')         # simple model, no internal dynamics
        bc.shifted_drive_vector.play(h.bcSyn._ref_g, global_params['stim_params']['dt'])      # Ph input (represents the Ph drive to BC)
        h.bcSyn.e= e  
        h.pop_section()
        self.syn = h.bcSyn

# Input (E or I cell)
class Input_Cell:
    def __init__(self, cluster, x, y, counter, global_params, ei):
        self.cluster= cluster       # Presynaptic cluster type
        self.x= x                   # Position
        self.y= y
        self.offset= 0
        # Output synapses
        self.output= []
        if ei == 0: # excitatory
            self.soma= h.Section(name=f'input_exc_soma_{counter}')
        else:
            self.soma= h.Section(name=f'input_inh_soma_{counter}')

        h.pt3dclear(sec=self.soma)
        h.pt3dadd(x+.1, y-.1, -10-ei*20, 10,sec=self.soma)
        h.pt3dadd(x-.1, y+.1, 0-ei*20, 10,sec=self.soma)

        # Connect input synapse from photoreceptor
        self.drive_vector= []                             # Processed BC signal coming from the convolution of trajectory and RF       
        self.shifted_drive_vector= h.Vector()
        self.input_vector= []

        self.soma.nseg= 1
        self.soma.insert('pas')  
        self.soma.e_pas= -60
        self.soma.g_pas= 0.0001

        
# Place BC/AC cells and synapses
def GA_Input_SetSynapses(global_params, model):
    np.random.seed(0)#global_params['job_id'])
    # Excitatory inputs (coming from bipolar cells), these will tile the retina 'num_input_types' fold, every 'dist_between_input_cells'
    # Potentially recompute number of synapses
    x_dist= []
    y_dist= []

    max_input_extent= 0
    for  sec in model.cell.InputDends: # Dendrites that get the inputs
        sec.push()                
        x_dist.append(h.x3d(0))
        x_dist.append(h.x3d(int(h.n3d()-1)))
        y_dist.append(h.y3d(0))
        y_dist.append(h.y3d(int(h.n3d()-1)))
        # Find max extent of the dendrites 
        max_input_extent= max(max_input_extent, h.distance(model.cell.soma(0.5), sec(1)))
        h.pop_section()
    model.max_input_extent=(max(x_dist)) # Linear distance along the preferred axis
    model.min_input_extent=(min(x_dist)) # Linear distance along the preferred axis


    # These are the input populations to be added
    input_population= [model.input_exc_list]
    if global_params['pre_inhibition']:
        input_population.append(model.input_inh_list)


    if (global_params['num_syn_replace_input_cell_dist'] > 0): # place inputs based on synapses
        total_dist=0
        for sec in model.cell.InputDends: # Dendrites that get the inputs
            total_dist+= sec.L
        x_list= []

        for syn in range(global_params['num_syn_replace_input_cell_dist']):                
            new_pos= np.random.uniform(0,int(total_dist))    
            find_pos=0
            for sec in model.cell.InputDends: # Dendrites that get the inputs
                if(new_pos>=find_pos)and(new_pos<=find_pos+sec.L ): # Found synaptic placement
                    sec.push()   
                    post_x= (h.x3d(0) + h.x3d(h.n3d() - 1))/2
                    post_y= (h.y3d(0) + h.y3d(h.n3d() - 1))/2
                    x_list.append(post_x)
                    cluster= -1 # Just  a placeholder
                    
                    for ei in range(len(input_population)): # go over excitatory and potentially inhibitory inputs
                        bc= Input_Cell(cluster, post_x, post_y, len(input_population[ei])+1, global_params, ei) 
                        input_population[ei].append(bc)                            
                        bc.output.append(Input_Synapse(bc, 0, sec, post_x, post_y, .5, global_params, e= -60*ei)) # inhibitory synapses reverse at -60 mV

                    h.pop_section()
                find_pos+= sec.L     
            # get clusters from x pos distribution
            for c in range(global_params['num_input_types']):
                threshold_min= np.percentile(x_list, 100 * (c)/ global_params['num_input_types'] )
                threshold_max= np.percentile(x_list, 100 * (c + 1)/ global_params['num_input_types'] )
                for ei in range(len(input_population)): # go over excitatory and potentially inhibitory inputs
                    for input_cell in input_population[ei]:
                        if(input_cell.x >= threshold_min) and (input_cell.x <= threshold_max):
                            input_cell.cluster= c
    else:   # Place inputs based on bipolar cells     
        for sec in model.cell.InputDends: # Dendrites that get the inputs
            sec.push()                
            numSyn= max(1, int(sec.L/global_params['dist_between_input_cells']) ) # Calc # inputs if separate input by global_params['distSyn'] micron
            if(model.name == 'single compartment'):
                numSyn=global_params['num_input_types']
            
            for pos in range(numSyn):                       # Over all Synapses
                dend_pos= (pos + 1) / ( numSyn + 1)         # Synaptic location for 1 input: 0.5, for 2 inputs: .33, .67 etc
                cum_pos= dend_pos * sec.L                   # Cumulative position
                dend_dist_so_far= 0
                dend_i= 0
                while (dend_i + 2 < h.n3d()):               # Find the location on the dendrite for the synapse
                    dend_dist_so_far+= ((h.x3d(dend_i) - h.x3d(dend_i + 1))**2 + (h.y3d(dend_i) - h.y3d(dend_i + 1))**2)**0.5
                    if(dend_dist_so_far) >= cum_pos:
                        break
                    dend_i+= 1
                # Postsynaptic/Presynaptic positions
                post_x= (h.x3d(dend_i) + h.x3d(dend_i + 1))/2
                post_y= (h.y3d(dend_i) + h.y3d(dend_i + 1))/2
                if(model.name == 'single compartment'):
                    post_x= -100 + 200 / max(1, global_params['num_input_types'] - 1) * pos
                    post_y= 0
                # Find the nearest presynaptic cell
                pre_x= np.random.normal((post_x / global_params['dist_between_input_cells']) * global_params['dist_between_input_cells'] , global_params['dist_between_input_cells_SD'] )
                pre_y= np.random.normal((post_y / global_params['dist_between_input_cells']) * global_params['dist_between_input_cells'] , global_params['dist_between_input_cells_SD'] )

                # Type of presynaptic cluster
                cluster= -1#math.floor(global_params['num_input_types'] * h.distance(model.cell[cell].soma(0.5), sec(dend_pos)) / model.max_input_extent[cell])
                if(model.name == 'RGC'):
                    cluster= max(0,min(global_params['num_input_types'] - 1, math.floor(global_params['num_input_types'] * (post_x - model.min_input_extent) / (model.max_input_extent - model.min_input_extent))))
                if(model.name == 'single compartment'):
                    cluster= pos
                                
                # Find input cells with same location and cluster type

                
                found_in= False
                for ei in range(len(input_population)): # go over excitatory and potentially inhibitory inputs
                    for bc in input_population[ei]:
                        if( ( (bc.x - pre_x)**2 + (bc.y - pre_y)**2 < (global_params['dist_between_input_cells'] / 2)**2) and (bc.cluster == cluster) ):
                            found_in= True
                            break
                    if not(found_in):   # create a new presynaptic cell
                        bc= Input_Cell(cluster, pre_x, pre_y, len(input_population[ei])+1, global_params, ei) #self, cluster, x, y, drive_vector, counter, dt
                        input_population[ei].append(bc)
                        
                    bc.output.append(Input_Synapse(bc, 0, sec, post_x, post_y, dend_pos, global_params,e= -60*ei)) # type: ignore # inhibitory synapses reverse at -60 mV

            h.pop_section()
    #print(len(model.input_exc_list), len(model.input_inh_list))

# Shifts the activation of the synaptic inputs based on the speed and direction of activation
def GA_UpdateInputSynapses(global_params, model):     
    # Make sure that the responses are within the simulation period by examining the 10% indices
    offset= []
    for bc in model.input_exc_list:
        offset.append( int(bc.x / global_params['stim_params']['speed'][-1] * np.cos(global_params['stim_params']['angle'][-1]) / global_params['stim_params']['dt'] + bc.y / global_params['stim_params']['speed'][-1] * np.sin(global_params['stim_params']['angle'][-1]) / global_params['stim_params']['dt']) )
        for syn in bc.output:
            syn.syn.gain= model.input_params["input_exc_syn_gain"]  
    
    if global_params['pre_inhibition']:
        for ac in model.input_inh_list:
            for syn in ac.output:
                syn.syn.gain= model.input_params["input_inh_syn_gain"] 


    min_10= global_params['stim_params']['delay'][-1]
    max_10= global_params['stim_params']['tStop'][-1]
    for cls in range(global_params['num_input_types']):
        syn_time= model.input_params['input_exc_syn_time'][cls]['drive'][-1] # make sure that the excitatory signal s within the stimulation duration
        peak_val = np.max(syn_time) * global_params['stim_params']['dt']
        at10_thresh = np.where(syn_time >= peak_val / 10)[0]
        if len(at10_thresh) > 0:
            min_10= min(min_10, at10_thresh[0])
        if len(at10_thresh) >= 2:
            max_10= min(max_10, at10_thresh[-1])

    min_offset= min_10 / global_params['stim_params']['dt'] + min(offset)
    if(min_offset < global_params['stim_params']['delay'][-1] / global_params['stim_params']['dt'] ):
        offset = [int(x + (global_params['stim_params']['delay'][-1] / global_params['stim_params']['dt'] - min_offset))  for x in offset]

    # Go over all input synapses (excitatory and inhibitory)
    pre_syn_pop= [model.input_exc_list]
    if(global_params['pre_inhibition']):
        pre_syn_pop.append( model.input_inh_list )
    
    for ei in pre_syn_pop:
        for i, bc in enumerate(ei):
            vec_size= int(global_params['stim_params']['tStop'][-1] / global_params['stim_params']['dt'])
            bc.shifted_drive_vector.resize(vec_size)
            bc.shifted_drive_vector.fill(0)
            if global_params['fullRFcomputation']:
                offset[i]= 0    # No shortcuts for drifting gratings  

            bc.offset= offset[i] 
            for tt in range(len(bc.drive_vector)):
                if (offset[i] + tt < vec_size) and (offset[i] + tt >= 0):
                    bc.shifted_drive_vector.x[offset[i] + tt]= bc.drive_vector[tt] * global_params['stim_params']['contrast'][-1]


