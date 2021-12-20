#------------------------------------------------------------------------------------------#
import numpy as np
import pandas as pd
from scipy.spatial import distance
import scipy
import h5py, time, random, os, sys, pyflagser, warnings, datetime
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from itertools import repeat
from morphological_types import *
from pyflagser import flagser_count_unweighted as fcu
from pyflagser import flagser_unweighted as fu
import matplotlib
#------------------------------------------------------------------------------------------#
t = time.process_time()
#------------------------------------------------------------------------------------------#
''' data '''
def king_file(mc):
  mc_file = h5py.File(f'../data/average/cons_locs_pathways_mc{mc}_Column.h5', 'r')
  populations, connections = mc_file.get('populations'), mc_file.get('connectivity')
  return populations, connections

#------------------------------------------------------------------------------------------#
''' Model Constructions '''
def Bio_M(m_type, populations, connections):

        for M_a in tqdm(m_type):
          for M_b in m_type:
                # spacial coordinates of the neurons in each neuronal m-type
                L_a = pd.DataFrame(np.matrix(populations[M_a]['locations']), columns = ['x', 'y', 'z'])
                L_b = pd.DataFrame(np.matrix(populations[M_b]['locations']), columns = ['x', 'y', 'z'])

                # distances between each neuron pathway group
                D_ = scipy.spatial.distance.cdist(L_a, L_b, 'euclidean')

                # Bins
                bins = np.arange(1, D_.max(), 75) - np.concatenate([[0],
                        np.array(np.ones(len(np.arange(1, D_.max(), 75)) - 1))])

                # Matrix of distance bins
                C_ = np.array(np.digitize(D_, bins))

                # Bin groups in matrix
                groups = np.array(range(len(bins))) + 1

                # Actual connections matrix
                a = np.array(connections[M_a][M_b]['cMat'])

                ab = pd.DataFrame(a)
                ab.to_csv("../output/Bio_M/reconstruction/" + str(M_a) + str(M_b) + ".csv", header = False, index = False)
#------------------------------------------------------------------------------------------#
def ER(mc, m_type, populations):
  # compute full list of neuron locations. Collect array, locate connections (pre/post list)
  locations = np.vstack(np.array(populations[i]['locations']) for i in m_type)
  array = np.load('../output/Bio_M/model/mc' + str(mc) + '_array.npy')

  N = 31346
  A = np.zeros((N, N), dtype = np.int8)

  for i in tqdm(range(N)):
    A[i,:] = np.random.rand(N) < 0.008
  B = np.array(A)
  np.save('../output/Erdos/model/erdos_renyi.npy', B)
#------------------------------------------------------------------------------------------#
def neuron_swap(m_type, populations, connections):

        for M_a in tqdm(m_type):
          for M_b in m_type:
                  # spacial coordinates of the neurons in each neuronal m-type
                  L_a = pd.DataFrame(np.matrix(populations[M_a]['locations']), columns = ['x', 'y', 'z'])
                  L_b = pd.DataFrame(np.matrix(populations[M_b]['locations']), columns = ['x', 'y', 'z'])

                  # distances between each neuron pathway group
                  D_ = scipy.spatial.distance.cdist(L_a, L_b, 'euclidean')

                  # Bins
                  bins = np.arange(1, D_.max(), 75) - np.concatenate([[0],
                  np.array(np.ones(len(np.arange(1, D_.max(), 75)) - 1))])

                  # Matrix of distance bins
                  C_ = np.array(np.digitize(D_, bins))

                  # Bin groups in matrix
                  groups = np.array(range(len(bins))) + 1

                  # Actual connections matrix
                  a = np.array(connections[M_a][M_b]['cMat'])

                  # Shuffle of each matrix
                  for aw in groups:
                    b = a[C_ == aw]
                    np.random.shuffle(b)
                    iz, jz = np.nonzero(C_ == aw)
                    for i, j, v in zip(iz, jz, b):
                            a[i, j] = v
                  ab = pd.DataFrame(a)
                  ab.to_csv("../output/GB/general_reconstruction/" + str(M_a) + str(M_b) + ".csv", header = False, index = False)
#------------------------------------------------------------------------------------------#
def subdivide_connectome(mc, m_type, populations, connections):
  A = np.load('../output/Bio_M/model/mc' + str(mc) + '_array.npy')


  array_size_by_morph = [np.array(connections[i][i]['cMat']).shape[0] for i in m_type]

  layer1 = np.sum(array_size_by_morph[0:6])
  layer2 = np.sum(array_size_by_morph[6:16])
  layer4 = np.sum(array_size_by_morph[16:28])
  layer5 = np.sum(array_size_by_morph[28:41])
  layer6 = np.sum(array_size_by_morph[41:])

  layers = [layer1, layer2, layer4, layer5, layer6]
  boundaries = np.cumsum(layers)
  boundaries = np.insert(boundaries, 0, 0)

  for i in range(len(boundaries)-1):
    for j in range(len(boundaries)-1):
      np.save('../output/BC/blocks/mc6_' + str(i) + str(j) + '_array.npy', \
            A[boundaries[i]:boundaries[i+1],boundaries[j]:boundaries[j+1]])
#------------------------------------------------------------------------------------------#
def stat_inputs(mc, m_type, populations, connections):

    # compute full list of neuron locations. Collect array, locate connections (pre/post list)
    locations = np.vstack(np.array(populations[i]['locations']) for i in m_type)
    array = np.load('../output/Bio_M/model/mc' + str(mc) + '_array.npy')
    X = np.where(array == 1)
    distance_list = np.array(np.sqrt(np.sum((locations[X[1]] - locations[X[0]])**2, axis = 1)))
    pre, post = X[0], X[1]

    counts, bins = np.histogram(distance_list, bins = np.arange(0, distance_list.max(), 75))
    np.set_printoptions(suppress=True)
    bins1 = bins[:-1]**2

    probability = np.array(counts/max(counts)).astype(float)

    return pre, post, locations, probability, bins1
#------------------------------------------------------------------------------------------#
def main_fast(pre, post, locations):
    new_pre = np.empty_like(pre)
    vertex_removal = np.array(range(len(pre)))
    mask = np.zeros(len(vertex_removal), dtype = bool)
    mask_counter = 0
    mask_limit = len(mask)//2
    randint = np.random.randint
    random = np.random.random
    search = np.searchsorted
    N = len(vertex_removal)
    tries = 0

    for n in tqdm(range(len(post))):
        if mask_counter == mask_limit:

            vertex_removal = np.delete(vertex_removal, np.where(mask))
            mask = np.zeros(len(vertex_removal), dtype = bool)
            mask_counter = 0
            mask_limit = len(mask)//2
            N = len(vertex_removal)

        loc_post = locations[post[n],]
        done = False

        while not done:
            tries += 1
            while True:
                index_index = randint(0, N)
                if not mask[index_index]:
                    break

            index = vertex_removal[index_index]
            i = pre[index]
            loc_pre = locations[i,]

            if i != post[n]:
                done = True
                new_pre[n] = i
                mask[index_index] = True
                mask_counter += 1

    return new_pre
#------------------------------------------------------------------------------------------#
def main_fast_geometric(pre, post, locations, probability, bins, model):
    new_pre = np.empty_like(pre)
    vertex_removal = np.array(range(len(pre)))
    mask = np.zeros(len(vertex_removal), dtype = bool)
    mask_counter = 0
    mask_limit = len(mask)//2
    randint = np.random.randint
    random = np.random.random
    search = np.searchsorted
    N = len(vertex_removal)
    tries = 0

    for n in tqdm(range(len(post))):
        if mask_counter == mask_limit:
            vertex_removal = np.delete(vertex_removal, np.where(mask))
            mask = np.zeros(len(vertex_removal), dtype = bool)
            mask_counter = 0
            mask_limit = len(mask)//2
            N = len(vertex_removal)

        loc_post = locations[post[n],]
        done = False

        while not done:
            tries += 1
            while True:
                index_index = randint(0, N)
                if not mask[index_index]:
                    break

            index = vertex_removal[index_index]
            i = pre[index]
            loc_pre = locations[i,]

            if i != post[n] and random() < probability[search(bins, \
                                        (loc_post[2]-loc_pre[2])**2 + \
                                        (loc_post[1] - loc_pre[1])**2 + \
                                        (loc_post[0] - loc_pre[0])**2) - 1]:
                done = True
                new_pre[n] = i
                mask[index_index] = True
                mask_counter += 1

    return new_pre
#------------------------------------------------------------------------------------------#
def save_model(new_pre, post, model, graphic__):
  arr = np.zeros((31346, 31346), dtype=np.int8)
  if graphic__ == 1:
    for i, j in zip(new_pre, post):
        arr[i,j] += 1
  elif graphic__ == 0:
    for i, j in zip(new_pre, post):
        arr[i,j] = 1

  np.save('../output/' + str(model) + '/model/' + str(model) + '.npy', arr)
#------------------------------------------------------------------------------------------#
def block_configuration(mc, populations, connections, graphic__):
  model = 'BC'
  M = np.array([])
  subdivides = np.arange(0, 5, 1)
  for i in tqdm(subdivides):
    for j in subdivides:
      M = np.array(np.load('../output/BC/blocks/mc' + str(mc) + '_' + str(i) + str(j) + '_array.npy'))

      locations = np.vstack(populations[i]['locations'] for i in populations.keys())
      pre, post = np.where(M == 1)

      distance_list = np.array(np.sqrt(np.sum((locations[post] - locations[pre])**2, axis = 1)))
      counts, bins = np.histogram(distance_list, bins = np.arange(0, distance_list.max(), 75))

      np.set_printoptions(suppress=True)
      bins, probability = bins[:-1]**2, np.array(counts/max(counts)).astype(float)
      ##########################################################################################
      new_pre = main_fast(pre, post, locations)
      ##########################################################################################
      arr = np.zeros((M.shape), dtype = np.int8)

      if graphic__ == 1:
        for k, l in zip(new_pre, post):
                arr[k,l] += 1
      elif graphic__ == 0:
        for k, l in zip(new_pre, post):
              arr[k,l] = 1

      np.save('../output/BC/new_blocks/mc6_' + str(i) + str(j) + '.npy', arr)

  new_array = np.vstack(np.hstack(np.array(np.load('../output/BC/new_blocks/mc6_'\
                 +str(i) + str(j) + '.npy')) for j in subdivides) for i in subdivides)

  np.save('../output/BC/new_blocks/block_configuration.npy', new_array)
#------------------------------------------------------------------------------------------#
def block_geometric_configuration(mc, populations, connections, graphic__):
  model = 'BGC'
  M = np.array([])
  subdivides = np.arange(0, 5, 1)
  for i in tqdm(subdivides):
    for j in subdivides:
      M = np.array(np.load(f'../output/BC/blocks/mc{str(mc)}_{str(i)}{str(j)}_array.npy'))

      locations = np.vstack(populations[i]['locations'] for i in populations.keys())
      pre, post = np.where(M == 1)

      distance_list = np.array(np.sqrt(np.sum((locations[post] - locations[pre])**2, axis = 1)))
      counts, bins = np.histogram(distance_list, bins = np.arange(0, distance_list.max(), 75))

      np.set_printoptions(suppress=True)
      bins, probability = bins[:-1]**2, np.array(counts/max(counts)).astype(float)
      ##########################################################################################
      new_pre = main_fast_geometric(pre, post, locations, probability, bins, model)
      ##########################################################################################
      arr = np.zeros((M.shape), dtype = np.int8)
      if graphic__ == 1:
        for k, l in zip(new_pre, post):
                arr[k,l] += 1
      elif graphic__ == 0:
        for k, l in zip(new_pre, post):
                arr[k,l] = 1
      np.save('../output/BGC/new_blocks/mc6_' +str(i)+str(j)+'.npy', arr)

  new_array = np.vstack(np.hstack(np.array(np.load('../output/' + str(model)+ '/new_blocks/mc6_' + str(i) +\
                                   str(j) + '.npy')) for j in subdivides) for i in subdivides)

  np.save('../output/BGC/new_blocks/block_geometric_configuration.npy', new_array)

#------------------------------------------------------------------------------------------#
def block_connection_counts(matrix):
  al = [0, 339, 7857, 12518, 18624, 31346]
  return np.array([[np.sum(matrix[i:j, k:l]) for k, l in zip(al, al[1:])] for i, j in zip(al, al[1:])])

def block_connection_densities(matrix):
  al = [0, 339, 7857, 12518, 18624, 31346]
  return np.array([[np.sum(matrix[i:j, k:l])/((j-i)*(l-k)) for k, l in zip(al, al[1:])] for i, j in zip(al, al[1:])])

def block_edge_densities(matrix):
  al = [0, 339, 7857, 12518, 18624, 31346]
  return np.array([[np.sum(matrix[i:j, k:l])/7822274 for k, l in zip(al, al[1:])] for i, j in zip(al, al[1:])])

def compute_block_var(old, new):
  old, new = block_edge_densities(old), block_edge_densities(new)
  return np.array([np.abs(old[i] - new[i]) for i in range(len(new))])


def morph_counts(M):
  morphs = [
    0, 58, 23, 91, 72, 52, 43, 28, 105, 61, 176, 456, 332, 267, 55, 5870, 168, 8, 20, 8, 38,
    121, 118, 94, 6, 2682, 62, 1097, 407, 34, 76, 19, 96, 208, 398, 203, 8, 25, 301, 2398,
    1996, 344, 70, 3182, 55, 16, 31, 3461, 460, 335, 200, 17, 67, 1641, 1450, 1737
    ]

  msum = np.cumsum(morphs)
  counts = [[np.sum(M[i:j,k:l]) for k, l in zip(msum, msum[1:])] for i, j in zip(msum, msum[1:])]
  counts = np.nan_to_num(np.log(counts))
  counts[counts < 0] = 0
  return counts

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
''' Graphs and Statistics '''
def plot_basics(x, y, title, height):
  plt.xlabel(x), plt.ylabel(y), plt.title(title, y=height)

def heat_map(layers, model, model_name, folder, model_type, names, Model):
  indices = np.arange(0, len(names), 1)
  values = layers.flatten()
  threshold = max(values) - (max(values) - min(values))/2
  plt.imshow(layers, cmap='Blues', interpolation='nearest')

  if Model != 2:

    for i in indices:
      for j in indices:
        plt.xticks(indices, names), plt.yticks(indices, names)
        if layers[i,j] > threshold:
          text = plt.text(j,i, f'{np.round(layers[i,j],4)}', ha='center', va='center', color='white')
          # text = plt.text(j,i, f'{layers[i,j]:.2e}', ha='center', va='center', color='white')
        else:
          text = plt.text(j,i, f'{np.round(layers[i,j], 4)}', ha='center', va='center', color='black')
          # text = plt.text(j,i, f'{layers[i,j]:.2e}', ha='center', va='center', color='black')
  else:
    plt.xticks(indices, names, rotation=90), plt.yticks(indices, names)
    plt.tick_params(labelsize=5)

  if Model == 7:
      plot_basics('Post-Synaptic Layers', 'Pre-Synaptic Layers', '', 1.06)
  else:
      # plot_basics('Post-Synaptic Layers', 'Pre-Synaptic Layers', f'{str(model_name)} Model: {model_type}', 1.06)
      plot_basics('Post-Synaptic Layers', 'Pre-Synaptic Layers', '', 1.06)

  plt.savefig('../output/plots/' + str(folder) + '/heat_map_layer_' + str(model) + '.png')
  plt.clf()

def distance_distribution(locations, old_array, new_array, model, folder, populations, model_name):
  pre, post = np.where(old_array == 1)
  new_pre, new_post = np.where(new_array == 1)

  distance_lists = np.array(np.sqrt(np.sum((locations[post] - locations[pre])**2, axis = 1)))
  new_distance_lists = np.array(np.sqrt(np.sum((locations[new_post] - locations[new_pre])**2, axis = 1)))

  # Arrange data to be in form of histogram
  counts, bins = np.histogram(distance_lists, np.arange(0, distance_lists.max(), 75))
  new_counts, new_bin = np.histogram(new_distance_lists, np.arange(0, new_distance_lists.max(), 75))
  norm_counts = [i/np.sum(counts) for i in counts]
  new_norm_counts = [i/np.sum(new_counts) for i in new_counts]

  the_length = len(norm_counts) - len(new_norm_counts)
  abs_the_length = np.abs(len(norm_counts) - len(new_norm_counts))

  if the_length > 0:
    new_norm_counts = np.append(new_norm_counts, np.zeros(abs_the_length))
  elif the_length < 0:
    norm_counts = np.append(norm_counts, np.zeros(abs_the_length))
  else:
    pass

  norm_diff = [np.abs(norm_counts[i]-new_norm_counts[i]) for i in range(len(new_norm_counts))]
  ##########################################################################################
  # Plot data
  if model == 'BioM':
    plt.hist(bins[:-1], bins, weights = counts, label = 'Bio-M MC', alpha = 0.5)
    # plt.title('N2N distance for Bio-M MC', y=1.06)
  else:
    plt.hist(bins[:-1], bins, weights = counts, label = 'Bio-M MC', alpha = 0.5)
    plt.hist(new_bin[:-1], new_bin, weights = new_counts, label = str(model_name) + ' Model', alpha = 0.5, color = 'orange')
    # plt.title('N2N distance for Bio-M MC Vs ' + str(model_name) + ' Model', y=1.06)
  plt.xlabel('distance between connected neurons ($\mu m$)'),   plt.ylabel('count')
  plt.legend()
  plt.savefig('../output/plots/' + str(folder) + '/' + str(model) + '_dist_distr.png')
  plt.clf()
  return f"D-D TV distance: \t\t\t{sum(norm_diff)}"

def graphical_rep(matrix, model, folder):
  plt.rcParams.update({
    "lines.color": "white",
    "patch.edgecolor": "white",
    "text.color": "black",
    "axes.facecolor": "black",
    "axes.edgecolor": "lightgray",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "grid.color": "lightgray",
    "figure.facecolor": "white",
    "figure.edgecolor": "white",
    "savefig.facecolor": "white",
    "savefig.edgecolor": "white"})

  plt.imshow(matrix, cmap = 'gray', interpolation = 'nearest')
  plot_basics('Post-Synaptic Neurons', 'Pre-Synaptic Neurons', '', 1)

  plt.xticks([0])
  plt.yticks([0, 31346], rotation = 45)
  plt.subplot().yaxis.set_label_position("right")
  plt.subplot().yaxis.tick_right()
  plt.axes([0, 0, 1, 1])
  plt.axis('off')
  plt.savefig('../output/plots/' + str(folder) + '/matrix_' + str(model) + '.png')
  plt.clf()

def degree(matrix, model, folder, model_name, Model):
  sd = np.sum(matrix, axis = 1) - np.sum(matrix, axis = 0)
  cum_sd, size = np.cumsum(sd), np.arange(31346)

  if Model != 7:
    plt.bar(size, height = sd)
    plot_basics('Neuron', 'Signed-Degree', str(model_name) + ' Model: Signed-Degree for each Neuron', 1.06)
    plt.savefig('../output/plots/' + str(folder) + '/' + str(model) + '_sd.png'), plt.clf()

    plt.bar(size, cum_sd)
    plot_basics('Neuron', 'Cumulative Signed-Degree', str(model_name) + ' Model: Cumulative degree for all Neurons', 1.06)
    plt.savefig('../output/plots/' + str(folder) + '/cumsum_degree_' +str(model) + '.png'), plt.clf()
  else:
    plt.bar(size, height = sd)
    plot_basics('Neuron', 'Signed-Degree', str(model_name) + ' MC: Signed-Degree for each Neuron', 1.06)
    plt.savefig('../output/plots/' + str(folder) + '/' + str(model) + '_sd.png'), plt.clf()

    plt.bar(size, cum_sd)
    plot_basics('Neuron', 'Cumulative Signed-Degree', str(model_name) + ' MC: Cumulative degree for all Neurons', 1.06)
    plt.savefig('../output/plots/' + str(folder) + '/cumsum_degree_' +str(model) + '.png'), plt.clf()

def topological_stats(array, a, b):
  simplices = fcu(array,directed=True)
  Betti = fu(array, min_dimension=a, max_dimension=b, directed=True)['betti']
  euler = np.sum(fcu(array, directed=True)[0::2]) - np.sum(fcu(array, directed=True)[1::2])
  return f"Simplices: \t\t\t\t{simplices} \nBetti: Dim3+ (Dim5+ Bio-M): \t\t{Betti} \nEuler Characteristic: \t\t\t{euler}"
