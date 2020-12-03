'''
Script takes the h5 files which contains information about the neurons and their locations, 
connectivity etc. We specifically take the file that is the average of all the instances.

TASK
What this script ultimately does is to rearrange connection instances in the whole matrix
according to distance dependence which is done through each sub matrix.

HOW
To achieve this, we first read in the relevant file (mc_file) and move into the locations 
directory. We make notes of all the morphological types. We want to create two matrices 
(L_a, L_b) consisting of neuron locations for both the pre and post synaptic neuron types. 

Next, we want to create a matrix which holds all distance values (D_) for the relevant neurons.
These values then get binned (C_) according to their distances from each other.

With this we then want to have ready the original connection matrix (a) and for each sub matrix
we shuffle the connections according to which bin the entries are in. Returned is a list of
CSV files which contain the shuffled sub matrices which then need to be put back together
to form the shuffled connectome.

NOTE: run using 'python general_biol_swap.py'
'''
import numpy as np
import h5py
import pandas as pd
import collections
import itertools
from scipy.spatial import distance
import scipy
import random
###########################################################################################
mc_file = h5py.File('../data/cons_locs_pathways_mc6_Column.h5', 'r')
populations = mc_file.get('populations')
###########################################################################################
m_type = [
    'L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 'L23_BP', 'L23_BTC', 
    'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC', 'L23_PC', 'L23_SBC', 
    'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC', u'L4_MC', 'L4_NBC', 'L4_NGC', 'L4_PC', 
    'L4_SBC', 'L4_SP', 'L4_SS', 'L5_BP', 'L5_BTC', 'L5_ChC', 'L5_DBC', 'L5_LBC', 'L5_MC', 
    'L5_NBC', 'L5_NGC', 'L5_SBC', 'L5_STPC', 'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC', 'L6_BP', 
    'L6_BPC', 'L6_BTC', 'L6_ChC', 'L6_DBC', 'L6_IPC', 'L6_LBC', 'L6_MC', 'L6_NBC', 'L6_NGC', 
    'L6_SBC', 'L6_TPC_L1', 'L6_TPC_L4', 'L6_UTPC'
    ]
###########################################################################################
for M_a in m_type:
    for M_b in m_type:
        # spacial coordinates of the neurons in each neuronal m-type
        L_a = pd.DataFrame(np.matrix(populations[M_a]['locations']), columns = ['x', 'y', 'z'])
        L_b = pd.DataFrame(np.matrix(populations[M_b]['locations']), columns = ['x', 'y', 'z'])
###########################################################################################        
        D_ = scipy.spatial.distance.cdist(L_a, L_b, 'euclidean')
###########################################################################################
        bins = np.arange(1, D_.max(), 100) - np.concatenate([[0], np.array(np.ones(len(np.arange(1, D_.max(), 100)) - 1))])
###########################################################################################
        C_ = np.array(np.digitize(D_, bins))
###########################################################################################
        groups = np.array(range(len(bins))) + 1
###########################################################################################
        connections = mc0_file.get('connectivity')
        a = np.array(connections[M_a][M_b]['cMat'])
###########################################################################################
        for aw in groups:
            b = a[C_ == aw]
            np.random.shuffle(b)
            iz, jz = np.nonzero(C_ == aw)
            for i, j, v in zip(iz, jz, b):
                a[i, j] = v
        ab = pd.DataFrame(a)
        ab.to_csv("../output/reconstruction/" + str(M_a) + str(M_b) + ".csv", header = False, index = False)
###########################################################################################
