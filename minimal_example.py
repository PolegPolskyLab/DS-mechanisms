import numpy as np
import RF_functions as rf
import tifffile
import matplotlib.pyplot as plt

# stimulation parameters
maxtime=1000    # simulation duration
arena_size=200  # the size (x,y) of the arena
deltat=1       # time step

# temporal frequency
myt=np.linspace(0,maxtime*deltat*0.001,maxtime)
mytf=1          # temporal frequency 
tf=np.zeros(maxtime)

tf[50:450] = mytf       # one direction 
tf[550:950] = -mytf     # the other direction

# create the stimulus, which has periods of stability and then motion (temporal frequency set by mytf); first in one direction then in the another
stimulus=rf.calc_sinegrating(tf, arena_size, spat_freq=2)

# 2 inputs driving the directional detector
# the inputs are positioned near each other (at x=-5 and x=5 coordinates)
input1 = rf.Input_cell(stimulus)
input1.x = -5
input2 = rf.Input_cell(stimulus)
input2.x = 5 

# different experiments
experiment = 'center_width'  #'center_width' 'amplitude' etc
if 'amplitude' in experiment:           # assymetry in weight to the postsynaptic cell, should not produce direction selective responses
    input2.RF_amp = 5

if 'center_kinetics' in experiment:     # different kinetics, known solution
    input2.RF_center_tau = 25

if 'center_width' in experiment:        # different extent of RFs, works when surround is incluyded (the surround has the same properties for both inputs)
    input1.RF_center_width = 50
    input1.RF_surround_amp = 1
    input2.RF_surround_amp = 1

if 'surround_amp' in experiment:        # difference in surround contribution to the RF
    input1.RF_surround_amp = 3
    input2.RF_surround_amp = 0

if 'surround_kinetics' in experiment:   # difference in surround kinetics
    input1.RF_surround_tau = 1
    input1.RF_surround_amp = 3
    input2.RF_surround_amp = 3

if 'surround_width' in experiment:      # difference in surround width
    input2.RF_surround_width = 50
    input1.RF_surround_amp = 3
    input2.RF_surround_amp = 3

# update the inputs
input1.update_RF()
input2.update_RF()

# linear summation at the detector
detector = input2.RF_response + input1.RF_response


# plot params
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)

mylegendsize=10

xpos=0.1
ypos=0.75
xsize=.8
ysize=0.2
ystep=0.23

def setmyaxes(myxpos,myypos,myxsize,myysize):
    ax=plt.axes([myxpos,myypos,myxsize,myysize])
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')  
    
plt.figure(figsize=(7.5,10.0))

#plt.imshow(stimulus[:,:,0])

# first cell, inputs and RF components
setmyaxes(xpos,ypos,xsize,ysize)
plt.plot(myt, input1.RF_center_spatial_sum,color='black', label='center drive')
if input1.RF_surround_amp > 0:
    plt.plot(myt, input1.RF_center_response,color='green', label='center response')
    plt.plot(myt, input1.RF_surround_spatial_sum,color='gray', label='surround drive')
    plt.plot(myt, input1.RF_surround_response,color=(0,1,0), label=' surround response')
plt.legend(loc=4, frameon=True, fontsize=mylegendsize)
plt.ylabel('input 1')
plt.title("Null direction               Preferred direction")

# second cell, inputs and RF components
setmyaxes(xpos,ypos-1*ystep,xsize,ysize)
plt.plot(myt, input2.RF_center_spatial_sum,color='black', label='center drive')
if   input2.RF_surround_amp > 0:
    plt.plot(myt, input2.RF_center_response,color='green', label='center response')
    plt.plot(myt, input2.RF_surround_spatial_sum,color='gray', label='surround  drive')
    plt.plot(myt, input2.RF_surround_response,color=(0,1,0), label='surround response')
plt.legend(loc=4, frameon=True, fontsize=mylegendsize)
plt.ylabel('input 2')

# first cell, total response
setmyaxes(xpos,ypos-2*ystep,xsize,ysize)
plt.plot(myt, input1.RF_response,color='blue', label='input 1')
plt.plot(myt, input2.RF_response,color='red', label='input 2')
plt.legend(loc=4, frameon=True, fontsize=mylegendsize)
plt.ylabel('full RF response')

# response of the detector
setmyaxes(xpos,ypos-3*ystep,xsize,ysize)
plt.plot(myt, detector,color='black', label='summed signal')
plt.legend(loc=4, frameon=True, fontsize=mylegendsize)
plt.ylabel('response')
plt.xlabel('time [s]')

plt.show()

# compute and print the direction selectivity index
dir1 = np.max(detector[100:400]) - np.min(detector[100:400])
dir2 = np.max(detector[600:900]) - np.min(detector[600:900])
dsi = np.abs(dir1 - dir2) / (dir1 + dir2)
print ('DSI=', dsi)


# save the stimulus
# data = np.transpose(stimulus, (2, 0, 1))
# tifffile.imwrite('image_stack.tif', data) 



        
    
