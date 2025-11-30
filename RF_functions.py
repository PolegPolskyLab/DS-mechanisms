import numpy as np
class Input_cell():
    def set_RF(self, width):
        # create a 2d gaussian for RF component (center or surround) spatial extent
        arena_size = int(self.stimulus.shape[0])
        x=np.arange(-arena_size//2,(arena_size//2),1)
        y=np.arange(-arena_size//2,(arena_size//2),1)
        x,y=np.meshgrid(x,y)
        gauss=np.exp(-((x-self.x)**2 / (width**2) + (y-self.y)**2 / (width)**2))
        gauss=gauss/np.sum(gauss)
        return gauss

    def lowpass(self, x, tau):
        # return a convolved response
        result = np.zeros_like(x)
        if tau <= 1:
            result=x
        else:
            result[0]=x[0]
            for i in range(0,x.shape[0]-1):
                result[i+1]=1.0/tau*(x[i]-result[i])+result[i]     
        return result
     
    def update_RF(self):
        # compute the spatial RF
        self.RF_center = self.set_RF(self.RF_center_width)
        self.RF_surround = self.set_RF(self.RF_surround_width)
        
        # extract the mean activation of the RF by the stimulus at each frame
        RF_center_time = self.stimulus *  self.RF_center[:, :, None]   # Convolve RF x stimulus
        self.RF_center_spatial_sum = np.sum(RF_center_time, axis= (0,1))     # Extract spatial activation intensity at each frame
        RF_surround_time = self.stimulus *  self.RF_surround[:, :, None]   # Convolve RF x stimulus
        self.RF_surround_spatial_sum = np.sum(RF_surround_time, axis= (0,1))     # Extract spatial activation intensity at each frame
        
        # compute the temporal responses
        if self.RF_center_hpass:
            self.RF_center_response = self.highpass(self.RF_center_spatial_sum, self.RF_center_tau)
        else:
            self.RF_center_response = self.lowpass(self.RF_center_spatial_sum, self.RF_center_tau)
        if self.RF_surround_hpass:
            self.RF_surround_response = self.highpass(self.RF_surround_spatial_sum, self.RF_surround_tau)
        else:
            self.RF_surround_response = self.lowpass(self.RF_surround_spatial_sum, self.RF_surround_tau)

        # the combined RF response is given by the difference between the center and the surround
        self.RF_response = self.RF_amp * (self.RF_center_response - self.RF_surround_response * self.RF_surround_amp)
        if self.norm:
            self.RF_response = (self.RF_response - np.nanmin(self.RF_response)) /  (np.nanmax(self.RF_response) - np.nanmin(self.RF_response))
        if self.rect:
            self.RF_response = np.clip(self.RF_response, 0, None)

    def __init__(self, stimulus):
        self.x = 0
        self.y = 0
        self.norm = False
        self.rect = False
        self.RF_amp = 1
        self.RF_center_width = 10
        self.RF_center_tau = 1
        self.RF_center_hpass = False
        self.RF_surround_width = 30
        self.RF_surround_tau = 5        
        self.RF_surround_hpass = False
        self.RF_surround_amp = 0    
        self.stimulus = stimulus    
        self.update_RF()

def calc_sinegrating(tf, arena_size, spat_freq):
    # drifting gratings stimulus
    def calc_singlewave(img_size, spat_freq, x, y):
        sinewave=np.sin((np.linspace(0,img_size-1,img_size) - x + y) / img_size*2.0*np.pi * spat_freq)
        return sinewave

    n=tf.shape
    maxtime=n[0]
    movie=np.zeros((arena_size,arena_size,maxtime))
    spat_wlength = arena_size/spat_freq
    velo=tf*spat_wlength/100.0   
    image=np.zeros((arena_size,arena_size))
    interim=calc_singlewave(arena_size,spat_freq,0,0)
    for i in range(arena_size):
        image[i,0::]=interim 
    for i in range(maxtime):
        movie[:,:,i]=np.roll(image,int(sum(velo[0:i])),axis=1)      

    movie=movie*0.5+0.5
    return movie