import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
import time

counter = 0
def show2d(matrix):
    plt.imshow(matrix, aspect='auto', origin='lower')
    plt.colorbar()
    plt.show()

def show3d(image):
    fig, ax = plt.subplots()
    im = ax.imshow(image[:, :, 0], cmap='gray', vmin=0, vmax=1)  # create once

    for i in range(image.shape[2]):
        im.set_data(image[:, :, i])  # update data, don't redraw
        ax.set_title(f"Slice {i}")
        plt.pause(0.01)

    plt.show()    

def GA_StimTrajectory(global_params, input_params, repeat= 0):
    start_time = time.time()
    # Compute sim duration
    tstop= max(1000, 300 + global_params['stim_params']['arena'] / global_params['stim_params']['speed'][-1] + global_params['stim_params']['delay'][-1] + global_params['stim_params']['duration'][-1])
    if global_params['stim_params']['type'] == 'dg':
        tstop = 1500
    global_params['stim_params']['tStop'].append(tstop)

    adj_tStop= int(global_params['stim_params']['tStop'][-1] / global_params['stim_params']['dt'])
    adj_size= int(global_params['stim_params']['arena'] / global_params['stim_params']['dx'])
    adj_del= int(global_params['stim_params']['delay'][-1] / global_params['stim_params']['dt'])    
    adj_dur= int(global_params['stim_params']['duration'][-1] / global_params['stim_params']['dt'])    
    adj_speed= global_params['stim_params']['speed'][-1] * global_params['stim_params']['dt'] / global_params['stim_params']['dx']

    if(global_params['stim_params']['angle'][-1] == 0): # Compute the stimulus
        trajectory = np.zeros((adj_size, adj_size, adj_tStop)) # Space, Space, Time
        
        # DS stimulus (moving bar)
        if (global_params['stim_params']['type'] == 'bar'):
            for x in range(adj_size):
                t_min= min(adj_tStop, max(0, int(x / adj_speed )))
                t_max= min(adj_tStop, max(0, int(adj_dur + x / adj_speed)))
                trajectory[x, :, t_min :  t_max]= 1


        # Drifting sine wave pattern        
        if (global_params['stim_params']['type'] == 'dg'): 
            trajectory+= .5 # Start from grayscale
            for t in range(adj_del, adj_tStop):
                for pos in range(0, max(0, min(adj_size , int((t - adj_del) * global_params['stim_params']['speed'][-1] * global_params['stim_params']['dt'])))):
                    trajectory[ pos , :, t]= np.clip(.5+.5*np.sin(2 * np.pi * (pos - (t - adj_del) * global_params['stim_params']['speed'][-1] * global_params['stim_params']['dt'])/ (global_params['stim_params']['cycle'][-1] / global_params['stim_params']['dx'])) , 0, 1)

        
        input_params['zero_trajectory']= trajectory
    else:   # Rotate saved stimulus
        trajectory = rotate(input_params['zero_trajectory'], angle=global_params['stim_params']['angle'][-1] / np.pi * 180, reshape=False)

    input_params['trajectory']= trajectory
   # show3d(trajectory)

def GA_Pre_activation(cell_type, input_params, global_params, cls, xcenter=0, ycenter=0):   
    pre_RF_components = ["center"]                                          # Compute RF center (always)
    if global_params['RF_constrains'][cell_type]['doSurround']: # Also do RF surround
        pre_RF_components = ["center", "surround"]
       
    dt= global_params['stim_params']['dt']
    dx= global_params['stim_params']['dx']
    arena= global_params['stim_params']['arena']

    conv = 2 * np.sqrt(np.log(2))
    RF_params= input_params['RF_params']
    time_steps= input_params['trajectory'].shape[2]
    spatial_points = int(arena / dx)

    def gauss2d(xc, yc, wx, wy, theta, arena, spatial_points):
        """
        Returns a 2D Gaussian centered at (xc, yc) with widths wx and wy.

        Parameters:
        - xc, yc: Center coordinates of the Gaussian
        - wx, wy: Standard deviations along x and y axes
        - arena: Grid size (creates a spatial_points x spatial_points array)

        Returns:
        - 2D numpy array of shape (spatial_points, spatial_points)
        """

        x = np.linspace(-arena / 2, arena / 2, spatial_points)
        y = np.linspace(-arena / 2, arena / 2, spatial_points)
        X, Y = np.meshgrid(x, y, indexing= 'ij')
        # Shift coordinates to center the Gaussian
        X_shift = X - xc
        Y_shift = Y - yc

        # Apply rotation
        X_rot = X_shift * np.cos(theta) + Y_shift * np.sin(theta)
        Y_rot = -X_shift * np.sin(theta) + Y_shift * np.cos(theta)

        # Compute Gaussian
        Z = np.exp(-((X_rot ** 2) / (2 * wx ** 2) + (Y_rot ** 2) / (2 * wy ** 2)))
        return Z

    spatial_sum_surround= []
    activation_surround= []
    inactivation_surround= []
    drive= []
    for cs in pre_RF_components:
        # 2D RF structure  
        RF_space = gauss2d(xcenter, ycenter, RF_params[cell_type][cs]['widthX'][cls] / conv, RF_params[cell_type][cs]['widthY'][cls] / conv, RF_params[cell_type][cs]['widthC'][cls], arena, spatial_points)

        # Convolve the spatial RF with the stimulus to get the fraction of the RF exposed to the stimulus
        RF_space_time = input_params['trajectory'] * RF_space[:, :, None]   # Convolve RF x stimulus
        spatial_sum= np.sum(RF_space_time, axis= (0,1))         # Extract spatial activation intensity at each frame
        if global_params['normalize_to_max']:
            spatial_sum /= np.max(spatial_sum)                  # Normalize to peak
        else:
            spatial_sum /= np.sum(RF_space)                     # normalize by the shape of the RF
        #spatial_sum*= input_params["fromPhsynG"]               # gain from photoreceptors 
        
        # Compute the temporal activation
        activation = np.zeros(time_steps)           # Not active
        inactivation = np.zeros(time_steps)          # Full RRP

        adjtau = (1 - 1 / RF_params[cell_type][cs]['tau'][cls]) ** dt
        adjtauRRP = (1 - 1 / float(RF_params[cell_type][cs]['tauRRP'][cls])) ** dt

        for tt in range(1, time_steps):   # Dynamics of the activation and inactivation
            activation[tt] = (( spatial_sum[tt - 1]- activation[tt - 1]) * (1 - adjtau) + activation[tt - 1]) * (1 - inactivation[tt - 1] )
            if global_params['RF_constrains'][cell_type]['inactivation']: # model RRP
                inactivation[tt]+= inactivation[tt - 1] + activation[tt] 
                inactivation[tt]= np.clip(inactivation[tt] * (adjtauRRP), 0, 1)
    
        if cs == 'surround': 
            activation*= RF_params[cell_type][cs]['peak'][cls]
        #show2d(RF_space)

        if(cs == 'center'):                         # Compute center response first
            drive= activation.copy()
            spatial_sum_center= spatial_sum
            activation_center= activation
            inactivation_center= inactivation
        else:               # Remove the surround from the center
            if global_params['relu']:
                drive = np.clip(drive - activation, 0, None ) # type: ignore
            else:
                drive = drive - activation
            spatial_sum_surround= spatial_sum
            activation_surround= activation 
            inactivation_surround= inactivation

    # combined amplitude 
    drive*= RF_params[cell_type]['center']['peak'][cls]
   # Synaptic plasticity
    facilitation= np.zeros(time_steps)
    depression= np.ones(time_steps)		# Start with full RRP, no depression
    STP= drive.copy()
    if(global_params['RF_constrains'][cell_type]['depression'] or global_params['RF_constrains'][cell_type]['facilitation']):	
        for tt in range(1, time_steps):
            if(global_params['RF_constrains'][cell_type]['depression'])	:											# depression:
                depression[tt]= (1 - RF_params[cell_type]['depression'][cls] *  drive[tt]) * depression[tt-1]
                depression[tt]= np.clip(depression[tt] + (1 - depression[tt]) / (RF_params[cell_type]['depression_tau'][cls] / dt), 0, 1)

            if(global_params['RF_constrains'][cell_type]['facilitation']):										 # facilitation
                facilitation[tt]= RF_params[cell_type]['facilitation'][cls] * drive[tt] + facilitation[tt-1]
                facilitation[tt]= np.clip(facilitation[tt] - (facilitation[tt]) / (RF_params[cell_type]['facilitation_tau'][cls] / dt), 0 ,1)

            STP[tt]=  drive[tt] * depression[tt] * (10**facilitation[tt])           # change in synaptic response
        #plt.plot(drive)
        #plt.plot(STP)
        #plt.plot(depression)        
        #plt.show()    
    return STP, spatial_sum_center, activation_center, inactivation_center, spatial_sum_surround, activation_surround, inactivation_surround # type: ignore
    
    
