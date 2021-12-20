import numpy as np
import h5py, time, random, os, sys, pyflagser, warnings, datetime
import matplotlib.pyplot as plt
from functions import *
from morphological_types import *
#------------------------------------------------------------------------------
screen_height, screen_length = os.popen('stty size', 'r').read().split()
warnings.filterwarnings('ignore')
#------------------------------------------------------------------------------
try:
    mc = sys.argv[1]
    Model = int(sys.argv[2])
    graphic__ = int(sys.argv[3])
except IndexError:
    raise SystemExit(f"Usage: {sys.argv[0]} <mc> <model>")
#------------------------------------------------------------------------------
screen_height, screen_length = os.popen('stty size', 'r').read().split()
warnings.filterwarnings('ignore')
populations, connections = king_file(mc)[0], king_file(mc)[1]
#------------------------------------------------------------------------------
def screen_refit(this_model):
        half = int(np.floor((int(screen_length) - len(this_model)))/2)
        print("-" * (half-1), this_model, '-' * (half-1))
#------------------------------------------------------------------------------
if Model == 1:
        screen_refit('Erdős–Rényi Model')
        ER(mc, m_type, populations)

elif Model == 2:
        screen_refit('General Biological Model')
        neuron_swap(m_type, populations, connections)

elif Model == 3:
        screen_refit('Configuration Model')
        model = 'configuration'
        pre, post, locations, _, _ = stat_inputs(mc, m_type, populations, connections)
        new_pre = main_fast(pre, post, locations)
        save_model(new_pre, post, model, graphic__)

elif Model == 4:
        screen_refit('Geometric Configuration Model')
        model = 'GC'
        pre, post, loc, prob, bins = stat_inputs(mc, m_type, populations, connections)
        new_pre = main_fast_geometric(pre, post, loc, prob, bins, model)
        save_model(new_pre, post, model, graphic__)

elif Model == 5:
        screen_refit('Block Configuration Model')
        subdivide_connectome(mc, m_type, populations, connections)
        block_configuration(mc, populations, connections, graphic__)

elif Model == 6:
        screen_refit('Block Geometric Configuration Model')
        subdivide_connectome(mc, m_type, populations, connections)
        block_geometric_configuration(mc, populations, connections, graphic__)

elif Model == 7:
        screen_refit('Bio-M MC')
        Bio_M(m_type, populations, connections)
