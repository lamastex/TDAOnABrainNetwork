import numpy as np
import matplotlib.pyplot as plt
import sys
from functions import *
from morphological_types import *

warnings.filterwarnings('ignore')

try:
    mc = sys.argv[1]
    Model = int(sys.argv[2])
    graphic__ = int(sys.argv[3])
except IndexError:
    raise SystemExit(f"Usage: {sys.argv[0]} <model>")

#-----------------------------------------------------------------------------
screen_height, screen_length = os.popen('stty size', 'r').read().split()
mc6 = np.load('../output/Bio_M/model/mc6_array.npy')
def screen_refit(this_model):
  half = int(np.floor((int(screen_length) - len(this_model)))/2)
  print("-" * (half-1), this_model, '-' * (half-1))
#-----------------------------------------------------------------------------
if Model == 1:
  screen_refit('Erdős–Rényi Model Statistics')
  A = np.load('../output/Erdos/model/erdos_renyi.npy')
  block_dens = block_connection_densities(A)
  populations, connections = king_file(mc)
  _, _, locations, _, _ = stat_inputs(mc, m_type, populations, connections)
  block_counts = block_connection_counts(A)
  connect_dens, block_var = block_edge_densities(A), compute_block_var(mc6, A)

  heat_map(connect_dens, erTest, erModelName, erFolder, erModelType, layer_name, Model)
  heat_map(block_dens, erModel, erModelName, erFolder, erModelType, layer_name, Model)
  a = distance_distribution(locations, mc6, A, erModel2, erFolder, populations, erModelName)
  b = topological_stats(A, 3, 5)
  heat_map(block_counts, erModel2, erModelName, erFolder, erModelType2, layer_name, Model)
  degree(A, erModel2, erFolder, erModelName, Model)
  graphical_rep(A, erModel2, erFolder)

  print(f'{a}\n{b}\n')

#-----------------------------------------------------------------------------
if Model == 2:
  screen_refit('General Biological Model Statistics')
  A = np.load('../output/GB/gen_biol.npy')

  mtype_dens = morph_counts(A)
  populations, connections = king_file(mc)
  _, _, locations, _, _ = stat_inputs(mc, m_type, populations, connections)
  connect_dens, block_var = block_edge_densities(A), compute_block_var(mc6, A)

  heat_map(connect_dens, gbTest, gbModelName, gbFolder, gbModelType, layer_name, Model)
  heat_map(mtype_dens, gbModel, gbModelName, gbFolder, gbModelType, m_type, Model)
  a = distance_distribution(locations, mc6, A, gbModel2, gbFolder, populations, gbModelName)
  b = topological_stats(A, 3, 5)

  block_counts = block_connection_counts(A)
  heat_map(block_counts, gbModel2, gbModelName, gbFolder, gbModelType, layer_name, Model)
  block_counts = morph_counts(A)
  heat_map(block_counts, gbModel2, gbModelName, gbFolder, gbModelType, m_type, Model)

  degree(A, gbModel2, gbFolder, gbModelName, Model)
  graphical_rep(A, gbModel2, gbFolder)

  print(f'{a}\n{b}\n')

#-----------------------------------------------------------------------------
if Model == 3:
  screen_refit('Configuration Model Statistics')
  A = np.load('../output/configuration/model/configuration.npy')

  if graphic__ == 1:
    block_dens = block_connection_densities(A)
    connect_dens = block_edge_densities(A)
    populations, connections = king_file(mc)
    _, _, locations, _, _ = stat_inputs(mc, m_type, populations, connections)
    block_counts, block_var = block_connection_counts(A), compute_block_var(mc6, A)

    heat_map(connect_dens, cTest, cModelName, cFolder, cModelType, layer_name, Model)
    heat_map(block_dens, cModel, cModelName, cFolder, cModelType, layer_name, Model)
    a = distance_distribution(locations, mc6, A, cModel2, cFolder, populations, cModelName)
    b = topological_stats(A, 3, 5)
    heat_map(block_counts, cModel2, cModelName, cFolder, cModelType2, layer_name, Model)
    degree(A, cModel2, cFolder, cModelName, Model)
    print(f'{a}\n{b}\n')

  elif graphic__ == 0:
    graphical_rep(A, cModel2, cFolder)

#-----------------------------------------------------------------------------
if Model == 4:
  screen_refit('Geometric Configuration Model Statistics')
  A = np.load('../output/GC/model/GC.npy')

  if graphic__ == 1:
    block_densities = block_connection_densities(A)
    populations, connections = king_file(mc)
    _, _, locations, _, _ = stat_inputs(mc, m_type, populations, connections)
    block_counts = block_connection_counts(A)
    connect_dens, block_var = block_edge_densities(A), compute_block_var(mc6, A)

    heat_map(connect_dens, gcTest, gcModelName, gcFolder, gcModelType, layer_name, Model)
    heat_map(block_densities, gcModel, gcModelName, gcFolder, gcModelType, layer_name, Model)
    a = distance_distribution(locations, mc6, A, gcModel2, gcFolder, populations, gcModelName)
    distance_distribution(locations, mc6, A, gcModel2, gcFolder, populations, gcModelName)
    b = topological_stats(A, 3, 6)
    heat_map(block_counts, gcModel2, gcModelName, gcFolder, gcModelType2, layer_name, Model)
    degree(A, gcModel2, gcFolder, gcModelName, Model)
    print(f'{a}\n{b}\n')

  elif graphic__ == 0:
    graphical_rep(A, gcModel2, gcFolder)


#-----------------------------------------------------------------------------
if Model == 5:
  screen_refit('Block Configuration Model Statistics')
  A = np.load('../output/BC/new_blocks/block_configuration.npy')

  if graphic__ == 1:
    block_densities = block_connection_densities(A)
    connect_dens = block_edge_densities(A)
    populations, connections = king_file(mc)
    _, _, locations, _, _ = stat_inputs(mc, m_type, populations, connections)
    block_counts, block_var = block_connection_counts(A), compute_block_var(mc6, A)

    heat_map(connect_dens, bcTest, bcModelName, bcFolder, bcModelType, layer_name, Model)
    heat_map(block_densities, bcModel, bcModelName, bcFolder, bcModelType, layer_name, Model)
    a = distance_distribution(locations, mc6, A, bcModel2, bcFolder, populations, bcModelName)
    b = topological_stats(A, 3, 6)
    heat_map(block_counts, bcModel2, bcModelName, bcFolder, bcModelType2, layer_name, Model)
    degree(A, bcModel2, bcFolder, bcModelName, Model)
    print(f'{a}\n{b}\n')

  elif graphic__ == 0:
    graphical_rep(A, bcModel2, bcFolder)

#-----------------------------------------------------------------------------
if Model == 6:
  screen_refit('Block Geometric Configuration Model Statistics')
  A = np.load('../output/BGC/new_blocks/block_geometric_configuration.npy')

  if graphic__ == 1:
    block_densities = block_connection_densities(A)
    connect_dens = block_edge_densities(A)
    populations, connections = king_file(mc)
    _, _, locations, _, _ = stat_inputs(mc, m_type, populations, connections)
    block_counts, block_var = block_connection_counts(A), compute_block_var(mc6, A)

    heat_map(connect_dens, bgcTest, bgcModelName, bgcFolder, bgcModelType, layer_name, Model)
    heat_map(block_densities, bgcModel, bgcModelName, bgcFolder, bgcModelType, layer_name, Model)
    distance_distribution(locations, mc6, A, bgcModel2, bgcFolder, populations, bgcModelName)
    b = topological_stats(A, 3, 6)
    heat_map(block_counts, bgcModel2, bgcModelName, bgcFolder, bgcModelType, layer_name, Model)
    degree(A, bgcModel2, bgcFolder, bgcModelName, Model)
    print(f'{a}\n{b}\n')

  elif graphic__ == 0:
    graphical_rep(A, bgcModel2, bgcFolder)

#-----------------------------------------------------------------------------
if Model == 7:
  screen_refit('Bio-M MC Statistics')
  A = np.load('../output/Bio_M/model/mc6_array.npy')

  block_densities = block_connection_densities(A)
  connect_dens = block_edge_densities(A)
  populations, connections = king_file(mc)
  _, _, locations, _, _ = stat_inputs(mc, m_type, populations, connections)

  heat_map(connect_dens, bioTest, bioModelName, bioFolder, bioModelType3, layer_name, Model)
  heat_map(block_densities, bioModel, bioModelName, bioFolder, bioModelType, layer_name, Model)
  a = distance_distribution(locations, mc6, A, bioModel2, bioFolder, populations, bioModelName)
  b = topological_stats(A, 5, 5)
  block_counts = block_connection_counts(A)

  heat_map(block_counts, bioModel2, bioModelName, bioFolder, bioModelType2, layer_name, Model)
  degree(A, bioModel2, bioFolder, bioModelName, Model)
  graphical_rep(A, bioModel2, bioFolder)
  print(f'{a}\n{b}\n')


