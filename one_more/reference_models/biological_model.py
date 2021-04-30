'''
# script to convert all h5 files into csv files after distance dependent shuffling.
# This creates csv files according to the general-bio model

'''
import numpy as np
import h5py
import pandas as pd
import collections
import itertools
from scipy.spatial import distance
import scipy
import random
import sys
from tqdm import tqdm

try:
  mc = sys.argv[1]
#   m_type = sys.argv[2]
except IndexError:
    raise SystemExit(f"usage: {sys.argv[0]} <mc 0/1/2/3/4/5/6> ")



mc_file = h5py.File('../my_model/data/average/cons_locs_pathways_mc' + str(mc) + '_Column.h5', 'r')
# mc_file = h5py.File('../my_model/data/pathway_files/pathway_individual_files/individuals/pathways_P14-17/pathways_mc' + str(mc) + '_Column.h5', 'r')
populations = mc_file.get('populations')
connections = mc_file.get('connectivity')
###########################################################################################
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
###########################################################################################
pro=0
for M_a in tqdm(full):
    for M_b in full:
        # print(pro)
        # tqdm(pro)
        pro=pro+1
        # spacial coordinates of the neurons in each neuronal m-type
        L_a = pd.DataFrame(np.matrix(populations[M_a]['locations']), columns = ['x', 'y', 'z'])
        L_b = pd.DataFrame(np.matrix(populations[M_b]['locations']), columns = ['x', 'y', 'z'])
###########################################################################################
        # distances between each neuron pathway group
        D_ = scipy.spatial.distance.cdist(L_a, L_b, 'euclidean')
###########################################################################################
        # Bins
        bins = np.arange(1, D_.max(), 100) - np.concatenate([[0],
        np.array(np.ones(len(np.arange(1, D_.max(), 100)) - 1))])
###########################################################################################
        # Matrix of distance bins
        C_ = np.array(np.digitize(D_, bins))
###########################################################################################
        # Bin groups in matrix
        groups = np.array(range(len(bins))) + 1
###########################################################################################
        # Actual connections matrix
        a = np.array(connections[M_a][M_b]['cMat'])
###########################################################################################
        ab = pd.DataFrame(a)
        ab.to_csv("../my_model/output/Bio_M/reconstruction/" + str(M_a) + str(M_b) + ".csv", header = False, index = False)
