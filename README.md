# working-manuscript-TopologicalDataAnalysisOnABrainNetwork
working-manuscript-TopologicalDataAnalysisOnABrainNetwork - for Kieran Barber

# General_biol
Requirements to run files:
pip install h5py, pip3 install pyflagser, cons_locs_pathways_mc6_Column 
which can be obtained from BBP website downloads: https://bbp.epfl.ch/nmc-portal/downloads

Sample run:
python general_biol_swap.py
OUTPUT: csv files stored in directory ../reconstruction

python3 general_biol_reconstruction.py
OUTPUT: npy file saved for future use and simplicial counts
[31346, 7822274, 35959272, 5578631, 47415, 31]

# Original (Bio-M)
Requirements to run as above (General_biol)

sample run: as above

# MMD_neuron
Requirements to run:
h5py, pyflagsercontain

sample run:
python MMD_neuron.py
OUTPUT: 3.3392330383480826
4.5146315509444
4.862046771079168
4.9274484113986246
5.017607294450558

# threeDsimplices
Requirements to run as above

sample run: as above (MMD_neuron)
OUTPUT: csv files to be used in 'Figure3.py'
