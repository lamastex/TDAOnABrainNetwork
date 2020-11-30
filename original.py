import numpy as np
import h5py
import pandas as pd
import collections
from mpl_toolkits.mplot3d import Axes3D
import itertools
from scipy.spatial import distance
import scipy
import random
import os

# Read in H5 file
mc0_file = h5py.File('../../pathway_average_files/cons_locs_pathways_mc0_Column.h5', 'r')
#############################################################################################
#                                 Neuron Type Locations                                     #
#############################################################################################
populations = mc0_file.get('populations')
connection = mc0_file.get('connectivity')

# # M_types
m_type = ['L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 'L23_BP', 'L23_BTC', 'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC', 'L23_PC', 'L23_SBC', 'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC', u'L4_MC', 'L4_NBC', 'L4_NGC', 'L4_PC', 'L4_SBC', 'L4_SP', 'L4_SS', 'L5_BP', 'L5_BTC', 'L5_ChC', 'L5_DBC', 'L5_LBC', 'L5_MC', 'L5_NBC', 'L5_NGC', 'L5_SBC', 'L5_STPC', 'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC', 'L6_BP', 'L6_BPC', 'L6_BTC', 'L6_ChC', 'L6_DBC', 'L6_IPC', 'L6_LBC', 'L6_MC', 'L6_NBC', 'L6_NGC', 'L6_SBC', 'L6_TPC_L1', 'L6_TPC_L4', 'L6_UTPC']

# ###########################################################################################
# #                                        All cases                                        #
# ###########################################################################################
for i in m_type:
    for j in m_type:
        connect4 = pd.DataFrame(np.matrix(connection[i][j]['cMat']))
        connect4.to_csv("../reconstruction/" + str(i) + str(j) + ".csv", header = False, index = False)
