import numpy as np
import h5py, time, random, os, sys, pyflagser, warnings, datetime
import matplotlib.pyplot as plt
from functions import *

warnings.filterwarnings('ignore')

try:
    mc = sys.argv[1]
    model = int(sys.argv[2])
except IndexError:
    raise SystemExit(f"Usage: {sys.argv[0]} <mc> <model>")

if model == 1:
        print('ER Model')
        ER(6)

elif model == 2:
        print('General Biological Model')
        neuron_swap(6)


# configuration
elif model == 3:
        print('Configuration Model')
        model = 'configuration'
        folder = 'configuration'
        #########################################################################################
        t = time.process_time()
        main_fast(required_inputs(6)[0], required_inputs(6)[1], required_inputs(6)[2], required_inputs(6)[3], model)
        print('runtime: ', time.process_time() - t)
        ##########################################################################################

# geometric configuration
elif model == 4:
        print('Geometric Configuration Model')
        model = 'GC'
        folder = 'geometric_configuration'
        ##########################################################################################
        t = time.process_time()
        main_fast_geometric(required_inputs(6)[0], required_inputs(6)[1], required_inputs(6)[2], required_inputs(6)[3], np.array(required_inputs(6)[4]), model, folder)
        print('runtime: ', time.process_time() - t)


# block configuration
elif model == 5:
        print('Block Configuration Model')
        complete_blocks(6)
        model = 'BC'
        folder = 'block_configuration'
        mc_file = h5py.File('../data/average/cons_locs_pathways_mc6_Column.h5', 'r')
        populations = mc_file.get('populations')
        connections = mc_file.get('connectivity')
        locations = required_inputs(6)[2]

        block_configuration(6)

# geometric block configuration
elif model == 6:
        print('Block Geometric Configuration Model')
        complete_blocks(6)
        model = 'BGC'
        folder = 'block_geometric_configuration'
        mc_file = h5py.File('../data/average/cons_locs_pathways_mc6_Column.h5', 'r')
        populations = mc_file.get('populations')
        connections = mc_file.get('connectivity')
        block_geometric_configuration(6)

elif model == 7:
        print('Bio_M Model')
        Bio_M(6)
