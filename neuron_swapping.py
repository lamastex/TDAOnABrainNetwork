import numpy as np
import h5py
import pandas as pd
import collections
from mpl_toolkits.mplot3d import Axes3D
import itertools
from scipy.spatial import distance
import scipy
import random

# Read in H5 file
mc0_file = h5py.File('../../pathway_average_files/cons_locs_pathways_mc0_Column.h5', 'r')

###########################################################################################
#                                 Neuron Type Locations                                   #
###########################################################################################
populations = mc0_file.get('populations')
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
        # distances between each neuron pathway group
        D_ = scipy.spatial.distance.cdist(L_a, L_b, 'euclidean')
###########################################################################################
        # Bins
        bins = np.arange(1, D_.max(), 100) - np.concatenate([[0], np.array(np.ones(len(np.arange(1, D_.max(), 100)) - 1))])
###########################################################################################
        # Matrix of distance bins
        B_ = np.digitize(D_, bins)
        C_ = np.array(pd.DataFrame((B_)))
###########################################################################################
        # Bin groups in matrix
        groups = np.array(range(len(bins))) + 1
###########################################################################################
        # Actual connections matrix
        connections = mc0_file.get('connectivity')
        a = np.array(connections[M_a][M_b]['cMat'])
###########################################################################################
        # Shuffle of each matrix
        for aw in groups:
            b = a[C_ == aw]
            np.random.shuffle(b)
            iz, jz = np.nonzero(C_ == aw)
            for i, j, v in zip(iz, jz, b):
                a[i, j] = v
        ab = pd.DataFrame(a)
        ab.to_csv("../reconstruction/" + str(M_a) + str(M_b) + ".csv", header = False, index = False)
###########################################################################################
# script to convert all h5 files into csv files after distance dependent shuffling. 
# This creates csv files according to the general-bio model
