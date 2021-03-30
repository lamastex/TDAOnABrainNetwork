import numpy as np
import h5py, time, random, os, sys, pyflagser
import matplotlib.pyplot as plt
from m_types import *

try:
    mc = sys.argv[1]
    A = int(sys.argv[2])
    B = int(sys.argv[3])
    start = int(sys.argv[4])
    stop = int(sys.argv[5])
    model = int(sys.argv[6])
except IndexError:
    raise SystemExit(f"Usage: {sys.argv[0]} <mc> <bin_accuracy> <prob_accuracy> <start> <stop> <model>")

if model == 1:

    for z in range(start, stop):
        print(z)
        mc_file = h5py.File('../data/average/cons_locs_pathways_mc' + str(mc) + '_Column.h5', 'r')
        populations = mc_file.get('populations')
        connections = mc_file.get('connectivity')
        ##########################################################################################
        # compute full list of neuron locations. Collect array, locate connections (pre/post list)
        locations = np.vstack(np.array(populations[i]['locations']) for i in m_type)
        array = np.load('../output/mc' + str(mc) + '_array.npy')
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

                    if i != post[n]:
                        done = True
                        new_pre[n] = i
                        mask[index_index] = True
                        mask_counter += 1

            return new_pre
        ##########################################################################################
        new_pre = main_fast(pre, post, locations, probability)
        ##########################################################################################
        # Rebuild configured array
        arr = np.zeros((len(array), len(array)), dtype = np.int8)

        for i, j in zip(new_pre, post):
            arr[i,j] = 1

        # Compute statistics
        simplices = pyflagser.flagser_count_unweighted(arr, directed=True)
        betti_numbers = pyflagser.flagser_unweighted(arr, min_dimension = 3, max_dimension = 5, directed = True)['betti']
        euler_char_simplices = np.sum(pyflagser.flagser_count_unweighted(arr, directed=True)[0::2]) - np.sum(pyflagser.flagser_count_unweighted(arr, directed=True)[1::2])
        
elif model == 2:
    for z in range(start, stop):
        mc_file = h5py.File('../data/average/cons_locs_pathways_mc' + str(mc) + '_Column.h5', 'r')
        populations = mc_file.get('populations')
        connections = mc_file.get('connectivity')
        ##########################################################################################
        # compute full list of neuron locations. Collect array, locate connections (pre/post list)
        locations = np.vstack(np.array(populations[i]['locations']) for i in m_type)
        array = np.load('../output/mc' + str(mc) + '_array.npy')
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

                    if i != post[n] and random() < probability[search(bins, (loc_post[2]-loc_pre[2])**2 + (loc_post[1] - loc_pre[1])**2 + (loc_post[0] - loc_pre[0])**2) - 1]:
                        done = True
                        new_pre[n] = i
                        mask[index_index] = True
                        mask_counter += 1

            return new_pre
        ##########################################################################################
        new_pre = main_fast(pre, post, locations, probability)
        ##########################################################################################
        # Rebuild configured array
        arr = np.zeros((len(array), len(array)), dtype = np.int8)

        for i, j in zip(new_pre, post):
            arr[i,j] = 1

        # Compute statistics
        simplices = pyflagser.flagser_count_unweighted(arr, directed=True)
        betti_numbers = pyflagser.flagser_unweighted(arr, min_dimension = 3, max_dimension = 5, directed = True)['betti']
        euler_char_simplices = np.sum(pyflagser.flagser_count_unweighted(arr, directed=True)[0::2]) - np.sum(pyflagser.flagser_count_unweighted(arr, directed=True)[1::2])

