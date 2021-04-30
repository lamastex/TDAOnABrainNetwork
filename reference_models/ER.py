import numpy as np
import pyflagser
import h5py
import matplotlib.pyplot as plt
import warnings
import sys
import time
import datetime
# import sys
# sys.path.append("../my_model/src/")
# from functions import *


print(datetime.datetime.now())
warnings.filterwarnings('ignore')

try:
    mc = sys.argv[1]
    # ins = int(sys.argv[2])
except IndexError:
    raise SystemExit(f"Usage: {sys.argv[0]} <mc> <No. of instances>")

t = time.process_time()

# for z in range(1, ins):
#   print(z)

mc_file = h5py.File('../my_model/data/average/cons_locs_pathways_mc' + str(mc) + '_Column.h5', 'r')
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
##########################################################################################
# compute full list of neuron locations. Collect array, locate connections (pre/post list)
locations = np.vstack(np.array(populations[i]['locations']) for i in m_type)
array = np.load('../my_model/output/Bio_M/model/mc' + str(mc) + '_array.npy')

N = 31346
A = np.zeros((N, N), dtype = np.int8)

for i in range(N):
  A[i,:] = np.random.rand(N) < 0.008
B = np.array(A)
np.save('../my_model/output/Erdos/model/erdos_renyi.npy', B)