import numpy as np
import pandas as pd
import os
import sys

try:
  ins = sys.argv[1]
  mtype = sys.argv[2]
except IndexError:
    raise SystemExit(f"usage: {sys.argv[0]} <mc 0/1/2/3/4/5/6> <m_types: f/e/i>")
###############################################################################
directory = "output/reconstruction/"
###############################################################################
full = [
    'L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 'L23_BP',
     'L23_BTC', 'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC',
      'L23_PC', 'L23_SBC', 'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC',
      'L4_MC', 'L4_NBC', 'L4_NGC', 'L4_PC', 'L4_SBC', 'L4_SP', 'L4_SS', 'L5_BP',
      'L5_BTC', 'L5_ChC', 'L5_DBC', 'L5_LBC', 'L5_MC', 'L5_NBC', 'L5_NGC', 'L5_SBC',
      'L5_STPC', 'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC', 'L6_BP', 'L6_BPC', 'L6_BTC',
      'L6_ChC', 'L6_DBC', 'L6_IPC', 'L6_LBC', 'L6_MC', 'L6_NBC', 'L6_NGC', 'L6_SBC',
      'L6_TPC_L1', 'L6_TPC_L4', 'L6_UTPC'
      ]

excitatory = [
        'L23_PC','L4_PC','L4_SP','L4_SS','L5_STPC','L5_TTPC1','L5_TTPC2','L5_UTPC','L6_BPC',
        'L6_IPC','L6_TPC_L1','L6_TPC_L4','L6_UTPC'
]

inhibitory = [
    'L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 'L23_BP', 'L23_BTC',
    'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC', 'L23_SBC',
    'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC', 'L4_MC', 'L4_NBC', 'L4_NGC',
    'L4_SBC',   'L5_BP', 'L5_BTC', 'L5_ChC', 'L5_DBC', 'L5_LBC', 'L5_MC', 'L5_NBC',
    'L5_NGC', 'L5_SBC',  'L6_BP', 'L6_BTC', 'L6_ChC', 'L6_DBC',  'L6_LBC', 'L6_MC',
    'L6_NBC', 'L6_NGC', 'L6_SBC'
]
dictionary = {"f": full, "e": excitatory, "i": inhibitory}
###############################################################################
# Code for computing m_TYPE_SIZE
m_type_size = []
for f in dictionary[mtype]:
    filename = f + f + ".csv"
    a = pd.read_csv(os.path.join(directory, filename), header = None, index_col = False)
    b = pd.DataFrame(a).shape[0]
    m_type_size.append(b)
###############################################################################
N = np.array(m_type_size).sum()
M = np.zeros((N, N), dtype = np.int8)
###############################################################################
for i in range(len(dictionary[mtype])):
    for j in range(len(dictionary[mtype])):
        filename = dictionary[mtype][i] + dictionary[mtype][j] + ".csv"
        M_ij = np.array(pd.read_csv(os.path.join(directory, filename),
        header = None, index_col = False), dtype = np.int8)

        k_start = np.array(m_type_size)[:i].sum()
        l_start = np.array(m_type_size)[:j].sum()
        for k in range(m_type_size[i]):
            for l in range(m_type_size[j]):
                M[k_start + k,l_start + l] = M_ij[k,l]
##############################################################################
np.save("output/array_mc' + str(ins) + '_' + str(mtype) + '.npy', M)
##############################################################################
