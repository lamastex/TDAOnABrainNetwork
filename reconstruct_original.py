import numpy as np
import h5py
import pandas as pd
import collections
from mpl_toolkits.mplot3d import Axes3D
import itertools
from scipy.spatial import distance
import scipy
import random
import glob
import os
import pyflagser
#######################################################################
directory = "../reconstruction/"
#######################################################################
m_type = [
    'L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 
    'L23_BP', 'L23_BTC', 'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 
    'L23_NBC', 'L23_NGC', 'L23_PC', 'L23_SBC', 'L4_BP', 'L4_BTC', 
    'L4_ChC', 'L4_DBC', 'L4_LBC', u'L4_MC', 'L4_NBC', 'L4_NGC', 'L4_PC', 
    'L4_SBC', 'L4_SP', 'L4_SS', 'L5_BP', 'L5_BTC', 'L5_ChC', 'L5_DBC', 
    'L5_LBC', 'L5_MC', 'L5_NBC', 'L5_NGC', 'L5_SBC', 'L5_STPC', 
    'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC', 'L6_BP', 'L6_BPC', 'L6_BTC', 
    'L6_ChC', 'L6_DBC', 'L6_IPC', 'L6_LBC', 'L6_MC', 'L6_NBC', 'L6_NGC', 
    'L6_SBC', 'L6_TPC_L1', 'L6_TPC_L4', 'L6_UTPC'
    ]
#######################################################################
m_type_size = np.array([
    59, 24, 91, 71, 52, 43, 28, 105, 61, 175, 453, 332, 266, 56, 5876, 
    167, 8, 20, 8, 37, 121, 117, 96, 6, 2681, 62, 1094, 410, 34, 75, 19, 
    96, 208, 398, 202, 8, 25, 301, 2402, 1997, 343, 70, 3179, 55, 16, 31, 
    3470, 466, 333, 198, 17, 66, 1640, 1441, 1737
    ], dtype = np.int64)
#######################################################################
N = m_type_size.sum()
M = np.zeros((N, N), dtype = np.int8)
#######################################################################
for i in range(len(m_type)):
    for j in range(len(m_type)):
        filename = m_type[i] + m_type[j] + ".csv"
        M_ij = np.array(pd.read_csv(os.path.join(directory, filename), header = None, index_col = False), dtype = np.int8)
        k_start = m_type_size[:i].sum()
        l_start = m_type_size[:j].sum()
        for k in range(m_type_size[i]):
            for l in range(m_type_size[j]):
                M[k_start + k,l_start + l] = M_ij[k,l]
#######################################################################
# np.save("Original_MC0", M)
#######################################################################
M = np.load("Original_MC0.npy")
print(pyflagser.flagser_count_unweighted(M, min_dimension=0, max_dimension = 7, directed=True))
#######################################################################
