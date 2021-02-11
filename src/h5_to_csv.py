import numpy as np
import h5py
import pandas as pd
import sys

try:
    mc = sys.argv[1]
    m_type = sys.argv[2]
except IndexError:
    raise SystemExit(f"usage: {sys.argv[0]} <mc: 0/1/2/3/4/5/6> <m_types: f/e/i>")

# t = time.process_time()

mc_file = h5py.File('../data/average/cons_locs_pathways_mc' + str(mc) + '_Column.h5', 'r')
populations = mc_file.get('populations')
###########################################################################################

full = [
    'L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 'L23_BP', 'L23_BTC',
    'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC', 'L23_PC', 'L23_SBC',
    'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC', 'L4_MC', 'L4_NBC', 'L4_NGC', 'L4_PC',
    'L4_SBC', 'L4_SP', 'L4_SS', 'L5_BP', 'L5_BTC', 'L5_ChC', 'L5_DBC', 'L5_LBC', 'L5_MC',
    'L5_NBC', 'L5_NGC', 'L5_SBC', 'L5_STPC', 'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC', 'L6_BP',
    'L6_BPC', 'L6_BTC', 'L6_ChC', 'L6_DBC', 'L6_IPC', 'L6_LBC', 'L6_MC', 'L6_NBC', 'L6_NGC',
    'L6_SBC', 'L6_TPC_L1', 'L6_TPC_L4', 'L6_UTPC'
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
###########################################################################################
for M_a in dictionary[m_type]:
    for M_b in dictionary[m_type]:
        connections = mc_file.get('connectivity')
        a = np.array(connections[M_a][M_b]['cMat'])
        ab = pd.DataFrame(a)
        ab.to_csv("../output/reconstruction/" + str(M_a) + str(M_b) + ".csv", header = False, index = False)
