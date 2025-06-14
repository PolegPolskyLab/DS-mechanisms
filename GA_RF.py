import numpy as np
import time
import matplotlib.pyplot as plt

#from scipy.ndimage import rotate, uniform_filter1d
counter = 0
def quick_plot(matrix, title="Matrix"):
    import matplotlib.pyplot as plt
    plt.imshow(matrix, aspect='auto', origin='lower')
    plt.title(title)
    plt.colorbar()
    plt.show()

def GA_StimTrajectory(stim_params, flash=0, noise1D=0, DS_noise=0, fuzzy=0, zigzag=0,
                      stim_duration=1, stim_width=1, stim_contrast=1.0):

    if stim_params is None:
        raise ValueError("stim_params must be provided.")

    dt = stim_params['dt']
    dx = stim_params['dx']
    arena = stim_params['arena']
    tStop = stim_params['tStop']
    speed = stim_params['speed']
    delay = stim_params['delay']
    duration = stim_params['duration']

    # Create Trajectory grid
    x = np.linspace(-arena / 2, arena / 2, int(arena / dx))
    y = np.linspace(0, tStop, int(tStop / dt))
    X, Y = np.meshgrid(x, y, indexing='ij')  

    Trajectory = np.zeros_like(X)

    # DS stimulus (moving bar)
    Trajectory = ((X / speed + delay >= Y - duration) & (X / speed + delay < Y)).astype(float)
    
    # Flash stimulus (overrides DS if flash is on)
    if flash:
        Trajectory = ((Y >= delay) & (Y < delay + duration)).astype(float)

    return Trajectory

def GA_Pre_activation(cell_type, RF_params, Trajectory, stim_params, pop, cls, xcenter=0, ycenter=0):   
    pre_RF_components = ["center", "surround"]
    dt = stim_params['dt']
    dx = stim_params['dx']
    arena = stim_params['arena']
    #tStop = stim_params['tStop']

    conv = 2 * np.sqrt(np.log(2))

    time_steps = Trajectory.shape[1]
    spatial_points = int(arena / dx)

    Pre_spatialSum = np.zeros(time_steps)


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

    for cs in pre_RF_components:
        # 2D RF structure
        RFspace = gauss2d(xcenter, ycenter, RF_params[cell_type][cs]['widthX'][cls][pop] / conv, RF_params[cell_type][cs]['widthY'][cls][pop] / conv, RF_params[cell_type][cs]['widthC'][cls][pop], arena, spatial_points)
        # Convolve the spatial RF with the stimulus to get the fraction of the RF exposed to the stimulus
        for tt in range(time_steps):
            if np.sum(Trajectory[:, tt]) > 0:
                stim = Trajectory[:, tt][:, None] * np.ones_like(RFspace)
                Pre_spatialSum[tt] = np.sum(stim * RFspace)

        Pre_spatialSum /= np.max(Pre_spatialSum)    # Normalize to peak
        
        # Compute the temporal activation
        PreTime = np.zeros(time_steps)      # Not active
        PreRRP = np.ones(time_steps)        # Full RRP

        adjtau = (1 - 1 / RF_params[cell_type][cs]['tau'][cls][pop]) ** dt
        adjtauRRP = (1 - 1 / RF_params[cell_type][cs]['tauRRP'][cls][pop]) ** dt

        for tt in range(1, time_steps):   # Dynamics of the activation and inactivation
            PreTime[tt] = (( Pre_spatialSum[tt - 1]- PreTime[tt - 1]) * (1 - adjtau) + PreTime[tt - 1]) * PreRRP[tt - 1]
            PreRRP[tt] = max(0, PreRRP[tt - 1] - PreTime[tt] * (1 - adjtauRRP))

        if(cs == 'center'): # Compute center response first
            PreDrive = np.zeros(time_steps)
            PreDrive = PreTime
        else:               # Remove the surround from the center
            PreDrive = PreDrive - PreTime * RF_params[cell_type][cs]['peak'][cls][pop]
    return PreDrive
    
    
