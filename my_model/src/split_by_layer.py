import h5py
import numpy as np
from tqdm import tqdm

A = np.load('../output/Bio_M/model/mc6_array.npy')

mc_file = h5py.File('../data/average/cons_locs_pathways_mc6_Column.h5', 'r')
populations = mc_file.get('populations')
connections = mc_file.get('connectivity')

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

array_size_by_morph = []
for i in m_type:
  array_size_by_morph.append(np.array(connections[i][i]['cMat']).shape[0])

layer1 = np.sum(array_size_by_morph[0:6])
layer2 = np.sum(array_size_by_morph[6:16])
layer4 = np.sum(array_size_by_morph[16:28])
layer5 = np.sum(array_size_by_morph[28:41])
layer6 = np.sum(array_size_by_morph[41:])

layers = [layer1, layer2, layer4, layer5, layer6]
boundaries = np.cumsum(layers)
boundaries = np.insert(boundaries, 0, 0)

for i in range(len(boundaries)-1):
  for j in range(len(boundaries)-1):
    np.save('../output/BC/blocks/mc6_' + str(i) + str(j) + '_array.npy', A[boundaries[i]:boundaries[i+1],boundaries[j]:boundaries[j+1]])
