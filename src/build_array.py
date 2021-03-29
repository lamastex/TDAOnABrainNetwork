import numpy as np
import h5py
import sys

try:
    mc = sys.argv[1]
except IndexError:
    raise SystemExit(f"usage: {sys.argv[0]} <mc: 0/1/2/3/4/5/6>")

mc_file = h5py.File('../data/average/cons_locs_pathways_mc' + str(mc) + '_Column.h5', 'r')
populations = mc_file.get('populations')
connections = mc_file.get('connectivity')
###########################################################################################
m_type = [
    'L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 'L23_BP',
     'L23_BTC', 'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC',
      'L23_PC', 'L23_SBC', 'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC',
      'L4_MC', 'L4_NBC', 'L4_NGC', 'L4_PC', 'L4_SBC', 'L4_SP', 'L4_SS', 'L5_BP',
      'L5_BTC', 'L5_ChC', 'L5_DBC', 'L5_LBC', 'L5_MC', 'L5_NBC', 'L5_NGC', 'L5_SBC',
      'L5_STPC', 'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC', 'L6_BP', 'L6_BPC', 'L6_BTC',
      'L6_ChC', 'L6_DBC', 'L6_IPC', 'L6_LBC', 'L6_MC', 'L6_NBC', 'L6_NGC', 'L6_SBC',
      'L6_TPC_L1', 'L6_TPC_L4', 'L6_UTPC'
      ]
###########################################################################################
m_type_size = []
for f in m_type:
        M_ij = np.array(connections[f][f]['cMat'], dtype = np.int8)
        m_type_size.append(M_ij.shape[0])
###########################################################################################
M = np.zeros((sum(m_type_size), sum(m_type_size)), dtype = np.int8)
###########################################################################################
for M_a, i in zip(m_type, range(len(m_type))):
    for M_b, j in zip(m_type, range(len(m_type))):
        M_ij = np.array(connections[M_a][M_b]['cMat'], dtype = np.int8)
        k_start = np.array(m_type_size)[:i].sum()
        l_start = np.array(m_type_size)[:j].sum()
        for k in range(m_type_size[i]):
            for l in range(m_type_size[j]):
                M[k_start + k,l_start + l] = M_ij[k,l]
###########################################################################################
np.save('../output/mc' + str(mc) + '_array.npy', M)
