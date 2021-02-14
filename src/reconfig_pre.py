import numpy as np
import h5py
import time
import random
import datetime
import os
import warnings
import sys

warnings.filterwarnings('ignore')

try:
    v = sys.argv[1]
    instance = sys.argv[2]
    A = int(sys.argv[3])
    B = int(sys.argv[4])
    mtype = sys.argv[5]
except IndexError:
    raise SystemExit(f"Usage: {sys.argv[0]} <mc> <instance> <bin_accuracy> <prob_accuracy (place same value as bin)> <mtype: f/e/i>")

mc_file = h5py.File('../data/average/cons_locs_pathways_mc' + str(v) + '_Column.h5', 'r')
populations = mc_file.get('populations')
connections = mc_file.get('connectivity')

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
##########################################################################################
locations = np.vstack(np.array(populations[i]['locations']) for i in dictionary[mtype])
array = np.load('../output/array_mc' + str(v) + '_' + str(mtype) + '.npy')
X = np.where(array == 1)
distance_list = np.array(np.sqrt(np.sum((locations[X[1]] - locations[X[0]])**2, axis = 1)))
pre = X[0]
post = X[1]
probability, bins = np.histogram(distance_list, bins = np.arange(0, distance_list.max(), 100))
np.set_printoptions(suppress=True)
bins = bins[:A-1]**2
probability = np.array(probability/max(probability)).astype(float)[:B]
##########################################################################################
def main_fast(pre, post, locations, probability):
    new_pre = np.empty_like(pre)
    removal_service_list = np.array(range(len(pre)))
    mask = np.zeros(len(removal_service_list), dtype = bool)
    mask_counter = 0
    mask_limit = len(mask)//2
    randint = np.random.randint
    random = np.random.random
    search = np.searchsorted
    N = len(removal_service_list)
    tries = 0

    for n in range(len(post))[:]:
        if mask_counter == mask_limit:
            removal_service_list = np.delete(removal_service_list, np.where(mask))
            mask = np.zeros(len(removal_service_list), dtype = bool)
            mask_counter = 0
            mask_limit = len(mask)//2
            N = len(removal_service_list)

        loc_post = locations[post[n],]
        done = False

        while not done:
            tries += 1
            while True:
                index_index = randint(0, N)
                if not mask[index_index]:
                    break

            index = removal_service_list[index_index]
            i = pre[index]
            loc_pre = locations[i,]

            if i != post[n] and random() < probability[search(bins, (post[2]-pre[2])**2 + (post[1] - pre[1])**2 + (post[0] - pre[0])**2) - 1]:
                done = True
                new_pre[n] = i
                mask[index_index] = True
                mask_counter += 1

    return new_pre

new_pre = main_fast(pre, post, locations, probability)
np.save('../output/instance_' + str(instance) +'_mc' + str(v) + '_' + str(mtype) + '.npy', new_pre)
print("Script 3 Complete!")
