import numpy as np
import h5py, time, random, os, sys, pyflagser, warnings, datetime
import matplotlib.pyplot as plt
from functions import *

warnings.filterwarnings('ignore')
print(datetime.datetime.now())

try:
    mc = sys.argv[1]
    start = int(sys.argv[2])
    stop = int(sys.argv[3])
    model = int(sys.argv[4])
except IndexError:
    raise SystemExit(f"Usage: {sys.argv[0]} <mc> <start> <stop> <model>")


# configuration
if model == 3:
        model = 'configuration'
        folder = 'configuration'
        #########################################################################################
        t = time.process_time()
        print('\nConfiguration Model: ')
        main_fast(required_inputs(6)[0], required_inputs(6)[1], required_inputs(6)[2], required_inputs(6)[3], model)
        print('runtime: ', time.process_time() - t)
        ##########################################################################################

# geometric configuration
elif model == 4:
        model = 'GC'
        folder = 'geometric_configuration'
        ##########################################################################################
        t = time.process_time()
        print('\nGeometric Configuration Model: ')
        main_fast_geometric(required_inputs(6)[0], required_inputs(6)[1], required_inputs(6)[2], required_inputs(6)[3], np.array(required_inputs(6)[4]), model, folder)
        print('runtime: ', time.process_time() - t)


# block configuration
elif model == 5:
        model = 'BC'
        folder = 'block_configuration'
        mc_file = h5py.File('../data/average/cons_locs_pathways_mc6_Column.h5', 'r')
        populations = mc_file.get('populations')
        connections = mc_file.get('connectivity')
        locations = required_inputs(6)[2]
        print('\nBlock Configuration Model: ')
        block_configuration(6)

# geometric block configuration
elif model == 6:
        model = 'GBC'
        folder = 'geometric_block_configuration'
        mc_file = h5py.File('../data/average/cons_locs_pathways_mc6_Column.h5', 'r')
        populations = mc_file.get('populations')
        connections = mc_file.get('connectivity')
        print('\nBlock Geometric Configuration Model: ')
        geometric_block_configuration(6)
