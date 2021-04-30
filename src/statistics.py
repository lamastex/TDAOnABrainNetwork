import numpy as np
import matplotlib.pyplot as plt
import sys
from functions import *

warnings.filterwarnings('ignore')
print(datetime.datetime.now())

try:
    model = int(sys.argv[1])
except IndexError:
    raise SystemExit(f"Usage: {sys.argv[0]} <model>")

mc6 = np.load('../output/Bio_M/model/mc6_array.npy')

# Erdos-Renyi
if model == 1:
  print('Erdos-Renyi: ')
  A = np.load('../output/Erdos/model/erdos_renyi.npy')
  model = 'Erdos-Renyi'
  folder = 'ER'

  # graphical_rep(A, model, folder)
  # layered_graph(A, model, folder)
  # heat_map(np.array(layer_counts(A)).reshape((5,5)), model, folder)
  # distance_distribution(required_inputs(6)[2], mc6, A, model, folder)
  # print(topological_stats(A, 3, 5))
  # degree(A, model, folder)
  # overall_degree(A, model, folder)

  model = 'Erdos-Renyi probability'
  # heat_map(np.round(np.array(layer_densities(A)).reshape((5,5)), 3), model, folder)

# general_biological
if model == 2:
  print('General Biological: ')
  A = np.load('../output/GB/gen_biol.npy')
  model = 'general_biol'
  folder = 'GB'

  # graphical_rep(A, model, folder)
  # layered_graph(A, model, folder)
  # heat_map(np.array(layer_counts(A)).reshape((5,5)), model, folder)
  # distance_distribution(required_inputs(6)[2], mc6, A, model, folder)
  # print(topological_stats(A, 3, 5))
  # degree(A, model, folder)
  # overall_degree(A, model, folder)

  # model = 'general_biol_densities'
  # heat_map(np.round(np.array(layer_densities(A)).reshape((5,5)), 3), model, folder)


# configuration
if model == 3:
  print('Configuration model: ')
  A = np.load('../output/configuration/model/configuration.npy')
  model = 'configuration'
  folder = 'configuration'

  # graphical_rep(A, model, folder)
  # layered_graph(A, model, folder)
  # heat_map(np.array(layer_counts(A)).reshape((5,5)), model, folder)
  # distance_distribution(required_inputs(6)[2], mc6, A, model, folder)
  # print(topological_stats(A, 3, 5))
  # degree(A, model, folder)
  # overall_degree(A, model, folder)


  # model = 'configuration_densities'
  # heat_map(np.round(np.array(layer_densities(A)).reshape((5,5)), 3), model, folder)

# geometric_configuration
if model == 4:
  print('Geometric Configuration model: ')
  A = np.load('../output/GC/model/GC.npy')
  model = 'GC'
  folder = 'GC'

  # graphical_rep(A, model, folder)
  # layered_graph(A, model, folder)
  # heat_map(np.array(layer_counts(A)).reshape((5,5)), model, folder)
  # distance_distribution(required_inputs(6)[2], mc6, A, model, folder)
  # print(topological_stats(A, 3, 5))
  # degree(A, model, folder)
  # overall_degree(A, model, folder)

  # model = 'GC_densities'
  # heat_map(np.round(np.array(layer_densities(A)).reshape((5,5)), 3), model, folder)

# block_configuration
if model == 5:
  print('Block Configuration model: ')
  A = np.load('../output/BC/new_blocks/block_configuration.npy')
  model = 'Block'
  folder = 'BC'

  # graphical_rep(A, model, folder)
  # layered_graph(A, model, folder)
  # heat_map(np.array(layer_counts(A)).reshape((5,5)), model, folder)
  # distance_distribution(required_inputs(6)[2], mc6, A, model, folder)
  # print(topological_stats(A, 3, 5))
  # degree(A, model, folder)
  # overall_degree(A, model, folder)

  # model = 'BC_densities'
  # heat_map(np.round(np.array(layer_densities(A)).reshape((5,5)), 3), model, folder)


# geometric_block_configuration
if model == 6:
  print('Geometric Block Configuration model: ')
  A = np.load('../output/GBC/new_blocks/geometric_block_configuration.npy')
  model = 'GBC'
  folder = 'GBC'

  # graphical_rep(A, model, folder)
  # layered_graph(A, model, folder)
  # heat_map(np.array(layer_counts(A)).reshape((5,5)), model, folder)
  # distance_distribution(required_inputs(6)[2], mc6, A, model, folder)
  # print(topological_stats(A, 3, 5))
  # degree(A, model, folder)
  # overall_degree(A, model, folder)

  # model = 'GBC_densities'
  # heat_map(np.round(np.array(layer_densities(A)).reshape((5,5)), 3), model, folder)

if model == 7:
  print('Bio-M Connectome: ')
  A = np.load('../output/Bio_M/model/mc6_array.npy')
  model = 'BioM'
  folder = 'Bio_M'

  # graphical_rep(A, model, folder)
  # layered_graph(A, model, folder)
  # heat_map(np.array(layer_counts(A)).reshape((5,5)), model, folder)
  # original_distance_distribution(required_inputs(6)[2], mc6, A, model, folder)
  # print(topological_stats(A, 5, 5))
  # degree(A, model, folder)
  # overall_degree(A, model, folder)

  # model = 'BioM_densities'
  # heat_map(np.round(np.array(layer_densities(A)).reshape((5,5)), 3), model, folder)
