import numpy as np
import matplotlib.pyplot as plt
import sys
from functions import *

warnings.filterwarnings('ignore')

try:
    model = int(sys.argv[1])
except IndexError:
    raise SystemExit(f"Usage: {sys.argv[0]} <model>")

mc6 = np.load('../output/Bio_M/model/mc6_array.npy')
mc_file = h5py.File('../data/average/cons_locs_pathways_mc6_Column.h5', 'r')
populations = mc_file.get('populations')
connections = mc_file.get('connectivity')


# Erdos-Renyi
if model == 1:
  print('Erdos-Renyi: ')
  A = np.load('../output/Erdos/model/erdos_renyi.npy')
  model = 'Erdos-Renyi'
  folder = 'ER'

  # print(distance_distribution())
  # heat_map(np.array(layer_counts(A)).reshape((5,5)), model, folder)
  # distance_distribution(required_inputs(6)[2], mc6, A, model, folder)
  print(topological_stats(A, 3, 5))
  # degree(A, model, folder)
  # overall_degree(A, model, folder)
  # cum_degree(A, model, folder)
  # graphical_rep(A, model, folder)
  # print('mean: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[1])
  # print('minimum: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[2])
  # print('maximum: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[3])
  # print('lower half: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[4])

  # model = 'Erdos-Renyi probability'
  # heat_map(np.round(np.array(layer_densities(A)).reshape((5,5)), 3), model, folder)

# general_biological
if model == 2:
  print('General Biological: ')
  A = np.load('../output/GB/gen_biol.npy')
  model = 'general_biol'
  folder = 'GB'

  # print(required_inputs(6)[5])

  # heat_map(np.array(layer_counts(A)).reshape((5,5)), model, folder)
  # distance_distribution(required_inputs(6)[2], mc6, A, model, folder)
  print(topological_stats(A, 3, 5))
  # degree(A, model, folder)
  # overall_degree(A, model, folder)
  # cum_degree(A, model, folder)
  # # print(min(morph_counts(A)))
  # heat_map_morph(np.array(morph_counts(A)).reshape((55,55)), model, folder)
  # graphical_rep(A, model, folder)
  # print('mean: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[1])
  # print('minimum: ',distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[2])
  # print('maximum: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[3])
  # print('lower half: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[4])


  # model = 'general_biol_densities'
  # heat_map(np.round(np.array(layer_densities(A)).reshape((5,5)), 3), model, folder)


# configuration
if model == 3:
  print('Configuration model: ')
  A = np.load('../output/configuration/model/configuration.npy')
  model = 'configuration'
  folder = 'configuration'

  # heat_map(np.array(layer_counts(A)).reshape((5,5)), model, folder)
  # distance_distribution(required_inputs(6)[2], mc6, A, model, folder)
  print(topological_stats(A, 3, 5))
  # degree(A, model, folder)
  # overall_degree(A, model, folder)
  # cum_degree(A, model, folder)
  # graphical_rep(A, model, folder)
  # print('mean: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[1])
  # print('minimum: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[2])
  # print('maximum: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[3])
  # print('lower half: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[4])


  # model = 'configuration_densities'
  # heat_map(np.round(np.array(layer_densities(A)).reshape((5,5)), 3), model, folder)

# geometric_configuration
if model == 4:
  print('Geometric Configuration model: ')
  A = np.load('../output/GC/model/GC.npy')
  model = 'GC'
  folder = 'GC'

  # heat_map(np.array(layer_counts(A)).reshape((5,5)), model, folder)
  distance_distribution(required_inputs(6)[2], mc6, A, model, folder)
  # print(topological_stats(A, 3, 6))
  # degree(A, model, folder)
  # overall_degree(A, model, folder)
  # cum_degree(A, model, folder)
  # graphical_rep(A, model, folder)

  # print('mean: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[1])
  # print('minimum: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[2])
  # print('maximum: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[3])
  # print('lower half: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[4])


  # model = 'GC_densities'
  # heat_map(np.round(np.array(layer_densities(A)).reshape((5,5)), 3), model, folder)

# block_configuration
if model == 5:
  print('Block Configuration model: ')
  A = np.load('../output/BC/new_blocks/block_configuration.npy')
  model = 'Block Configuration'
  folder = 'BC'

  # heat_map(np.array(layer_counts(A)).reshape((5,5)), model, folder)
  # distance_distribution(required_inputs(6)[2], mc6, A, model, folder)
  print(topological_stats(A, 3, 5))
  # degree(A, model, folder)
  # overall_degree(A, model, folder)
  # cum_degree(A, model, folder)
  # graphical_rep(A, model, folder)

  # print('mean: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[1])
  # print('minimum: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[2])
  # print('maximum: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[3])
  # print('lower half: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[4])


  # model = 'BC_densities'
  # heat_map(np.round(np.array(layer_densities(A)).reshape((5,5)), 3), model, folder)


# geometric_block_configuration
if model == 6:
  print('Block Geometric Configuration model: ')
  A = np.load('../output/BGC/new_blocks/block_geometric_configuration.npy')
  model = 'BGC'
  folder = 'GBC'

  # heat_map(np.array(layer_counts(A)).reshape((5,5)), model, folder)
  # distance_distribution(required_inputs(6)[2], mc6, A, model, folder)
  print(topological_stats(A, 3, 6))
  # degree(A, model, folder)
  # overall_degree(A, model, folder)
  # cum_degree(A, model, folder)

  # graphical_rep(A, model, folder)
  # print('mean: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[1])
  # print('minimum: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[2])
  # print('maximum: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[3])
  # print('lower half: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[4])


  # model = 'GBC_densities'
  # heat_map(np.round(np.array(layer_densities(A)).reshape((5,5)), 3), model, folder)

if model == 7:
  print('Bio-M MC: ')
  A = np.load('../output/Bio_M/model/mc6_array.npy')
  model = 'BioM'
  folder = 'Bio-M'
  graphical_rep(A, model, folder)

  # print(required_inputs(6)[5])
  # original_distance_distribution(required_inputs(6)[2], mc6, A, model, folder)
  # heat_map(np.array(layer_counts(A)).reshape((5,5)), model, folder)
  print(topological_stats(A, 5, 5))
  # degree(A, model, folder)
  # overall_degree(A, model, folder)
  # cum_degree(A, model, folder)
  # graphical_rep(A, model, folder)

  # print('mean: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[1])
  # print('minimum: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[2])
  # print('maximum: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[3])
  # print('lower half: ', distance_distribution(required_inputs(6)[2], mc6, A, model, folder)[4])


  # model = 'BioM_densities'
  # heat_map(np.round(np.array(layer_densities(A)).reshape((5,5)), 3), model, folder)
