README: Genetic Algorithm Retinal Simulation

This repository contains code to run genetic algorithm (GA)–driven optimization of receptive field (RF) parameters in a multicompartmental NEURON model. 
The framework allows exploring excitatory and inhibitory presynaptic populations, different stimulus conditions (moving bars, drifting gratings), and synaptic mechanisms such as short-term depression and facilitation.

Requirements:

Python 3.8+

NEURON simulator
 with Python bindings (pip install neuron)

numpy
matplotlib
argparse (included in Python stdlib)
Multiprocessing support (multiprocessing, included in stdlib)

Running the simulation:

Make sure to compile the MOD files (located in the 'MOD' subfolder
Create 'results' or 'output' directories in the current folder if those do not exist

run:
python main.py

This runs the GA with default parameters as defined in global_params.


Command-line arguments:
You can override defaults using flags:

--cellType
Selects the cell type. Options include:

"RGC" (default)
"single compartment"

--job_id
Assigns a numeric job ID to the simulation (useful for cluster runs).

Example:

python main.py --cellType RGC --job_id 42


Key Concepts:

At the top of the file you can set:

experiment_type = ['RGC']  
# or ['RGC', 'pre_inhibition']  
# or ['single compartment']


RGC: full retinal ganglion cell model
single compartment: simplified soma-only model
pre_inhibition: include inhibitory presynaptic populations


Global parameters:

All tunable parameters are in global_params. Some important ones:
numPop: population size (default 10)
numGen: number of generations (default 100)
mutationRate: mutation rate (default 0.1)

Circuit/Cell parameters:
num_input_types: number of presynaptic clusters (default is 4)
RF_constrains: defines which RF parameters can vary (amplitude, kinetics, size, orientation, etc.)
depression, facilitation: whether short-term plasticity is active

Simulation control:
switch_to_spiking_model: generation at which model produces spikes (set to high value to avoid spikes)
stim_params: visual stimulus parameters (bar or drifting grating)
debugger: flags for running NEURON, multithreading, etc.

Workflow:
Define experiment type (experiment_type list).
Adjust global parameters as needed in global_params.
Run simulation using python main.py.
The genetic algorithm (GA_loop) evolves RF parameters to optimize model performance.
Results are saved to H5 file in the 'results' folder.