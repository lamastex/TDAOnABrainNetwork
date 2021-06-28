import numpy as np
import pandas as pd
from scipy.spatial import distance
import scipy
import h5py, time, random, os, sys, pyflagser, warnings, datetime
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
#------------------------------------------------------------------------------------------#
t = time.process_time()
#------------------------------------------------------------------------------------------#

m_type = [
    'L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 'L23_BP',
    'L23_BTC', 'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC',
    'L23_PC', 'L23_SBC', 'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC',
    'L4_MC', 'L4_NBC', 'L4_NGC', 'L4_PC', 'L4_SBC', 'L4_SP', 'L4_SS', 'L5_BP',
    'L5_BTC', 'L5_ChC', 'L5_DBC', 'L5_LBC', 'L5_MC', 'L5_NBC', 'L5_NGC', 'L5_SBC',
    'L5_STPC', 'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC', 'L6_BP', 'L6_BPC', 'L6_BTC',
    'L6_ChC', 'L6_DBC', 'L6_IPC', 'L6_LBC', 'L6_MC', 'L6_NBC', 'L6_NGC', 'L6_SBC',
    'L6_TPC_L1', 'L6_TPC_L4', 'L6_UTPC'
    ]

layer1 = ['L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC']

layer2 = [
  'L23_BP','L23_BTC', 'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC',
  'L23_PC', 'L23_SBC'
]

layer4 = [
  'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC',
  'L4_MC', 'L4_NBC', 'L4_NGC', 'L4_PC', 'L4_SBC', 'L4_SP', 'L4_SS'
]

layer5 = [
  'L5_BP','L5_BTC', 'L5_ChC', 'L5_DBC', 'L5_LBC', 'L5_MC', 'L5_NBC', 'L5_NGC', 'L5_SBC',
  'L5_STPC', 'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC'
]

layer6 = [
  'L6_BP', 'L6_BPC', 'L6_BTC','L6_ChC', 'L6_DBC', 'L6_IPC', 'L6_LBC', 'L6_MC',
  'L6_NBC', 'L6_NGC', 'L6_SBC','L6_TPC_L1', 'L6_TPC_L4', 'L6_UTPC'
]

mc_file = h5py.File('../data/average/cons_locs_pathways_mc6_Column.h5', 'r')
populations = mc_file.get('populations')
connections = mc_file.get('connectivity')
#------------------------------------------------------------------------------------------#

def Bio_M(mc):
        mc_file = h5py.File('../data/average/cons_locs_pathways_mc' + str(mc) + '_Column.h5', 'r')
        populations = mc_file.get('populations')
        connections = mc_file.get('connectivity')
        ###########################################################################################
        full = [
        'L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 'L23_BP',
        'L23_BTC', 'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC',
        'L23_PC', 'L23_SBC', 'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC',
        'L4_MC', 'L4_NBC', 'L4_NGC', 'L4_PC', 'L4_SBC', 'L4_SP', 'L4_SS', 'L5_BP',
        'L5_BTC', 'L5_ChC', 'L5_DBC', 'L5_LBC', 'L5_MC', 'L5_NBC', 'L5_NGC', 'L5_SBC',
        'L5_STPC', 'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC', 'L6_BP', 'L6_BPC', 'L6_BTC',
        'L6_ChC', 'L6_DBC', 'L6_IPC', 'L6_LBC', 'L6_MC', 'L6_NBC', 'L6_NGC', 'L6_SBC',
        'L6_TPC_L1', 'L6_TPC_L4', 'L6_UTPC'
        ]
        ###########################################################################################
        for M_a in tqdm(full):
          for M_b in full:
                # spacial coordinates of the neurons in each neuronal m-type
                L_a = pd.DataFrame(np.matrix(populations[M_a]['locations']), columns = ['x', 'y', 'z'])
                L_b = pd.DataFrame(np.matrix(populations[M_b]['locations']), columns = ['x', 'y', 'z'])
        ###########################################################################################
                # distances between each neuron pathway group
                D_ = scipy.spatial.distance.cdist(L_a, L_b, 'euclidean')
        ###########################################################################################
                # Bins
                bins = np.arange(1, D_.max(), 75) - np.concatenate([[0],
                np.array(np.ones(len(np.arange(1, D_.max(), 75)) - 1))])
        ###########################################################################################
                # Matrix of distance bins
                C_ = np.array(np.digitize(D_, bins))
        ###########################################################################################
                # Bin groups in matrix
                groups = np.array(range(len(bins))) + 1
        ###########################################################################################
                # Actual connections matrix
                a = np.array(connections[M_a][M_b]['cMat'])
        ###########################################################################################
                ab = pd.DataFrame(a)
                ab.to_csv("../output/Bio_M/reconstruction/" + str(M_a) + str(M_b) + ".csv", header = False, index = False)


def ER(mc):
  mc_file = h5py.File('../data/average/cons_locs_pathways_mc' + str(mc) + '_Column.h5', 'r')
  populations = mc_file.get('populations')
  connections = mc_file.get('connectivity')

  m_type = [
      'L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 'L23_BP',
      'L23_BTC', 'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC',
        'L23_PC', 'L23_SBC', 'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC',
        'L4_MC', 'L4_NBC', 'L4_NGC', 'L4_PC', 'L4_SBC', 'L4_SP', 'L4_SS', 'L5_BP',
        'L5_BTC', 'L5_ChC', 'L5_DBC', 'L5_LBC', 'L5_MC', 'L5_NBC', 'L5_NGC', 'L5_SBC',
        'L5_STPC', 'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC', 'L6_BP', 'L6_BPC', 'L6_BTC',
        'L6_ChC', 'L6_DBC', 'L6_IPC', 'L6_LBC', 'L6_MC', 'L6_NBC', 'L6_NGC', 'L6_SBC',
        'L6_TPC_L1', 'L6_TPC_L4', 'L6_UTPC'
        ]
  ##########################################################################################
  # compute full list of neuron locations. Collect array, locate connections (pre/post list)
  locations = np.vstack(np.array(populations[i]['locations']) for i in m_type)
  array = np.load('../output/Bio_M/model/mc' + str(mc) + '_array.npy')

  N = 31346
  A = np.zeros((N, N), dtype = np.int8)

  for i in tqdm(range(N)):
    A[i,:] = np.random.rand(N) < 0.008
  B = np.array(A)
  np.save('../output/Erdos/model/erdos_renyi.npy', B)


def neuron_swap(mc):
        mc_file = h5py.File('../data/average/cons_locs_pathways_mc' + str(mc) + '_Column.h5', 'r')
        populations = mc_file.get('populations')
        connections = mc_file.get('connectivity')
        ###########################################################################################
        full = [
        'L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 'L23_BP',
        'L23_BTC', 'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC',
        'L23_PC', 'L23_SBC', 'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC',
        'L4_MC', 'L4_NBC', 'L4_NGC', 'L4_PC', 'L4_SBC', 'L4_SP', 'L4_SS', 'L5_BP',
        'L5_BTC', 'L5_ChC', 'L5_DBC', 'L5_LBC', 'L5_MC', 'L5_NBC', 'L5_NGC', 'L5_SBC',
        'L5_STPC', 'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC', 'L6_BP', 'L6_BPC', 'L6_BTC',
        'L6_ChC', 'L6_DBC', 'L6_IPC', 'L6_LBC', 'L6_MC', 'L6_NBC', 'L6_NGC', 'L6_SBC',
        'L6_TPC_L1', 'L6_TPC_L4', 'L6_UTPC'
        ]
        ###########################################################################################
        for M_a in tqdm(full):
          for M_b in full:
                  # spacial coordinates of the neurons in each neuronal m-type
                  L_a = pd.DataFrame(np.matrix(populations[M_a]['locations']), columns = ['x', 'y', 'z'])
                  L_b = pd.DataFrame(np.matrix(populations[M_b]['locations']), columns = ['x', 'y', 'z'])
          ###########################################################################################
                  # distances between each neuron pathway group
                  D_ = scipy.spatial.distance.cdist(L_a, L_b, 'euclidean')
          ###########################################################################################
                  # Bins
                  bins = np.arange(1, D_.max(), 75) - np.concatenate([[0],
                  np.array(np.ones(len(np.arange(1, D_.max(), 75)) - 1))])
          ###########################################################################################
                  # Matrix of distance bins
                  C_ = np.array(np.digitize(D_, bins))
          ###########################################################################################
                  # Bin groups in matrix
                  groups = np.array(range(len(bins))) + 1
          ###########################################################################################
                  # Actual connections matrix
                  a = np.array(connections[M_a][M_b]['cMat'])
          ###########################################################################################
                  # Shuffle of each matrix
                  for aw in groups:
                    b = a[C_ == aw]
                    np.random.shuffle(b)
                    iz, jz = np.nonzero(C_ == aw)
                    for i, j, v in zip(iz, jz, b):
                            a[i, j] = v
                  ab = pd.DataFrame(a)
                  ab.to_csv("../output/GB/general_reconstruction/" + str(M_a) + str(M_b) + ".csv", header = False, index = False)


def complete_blocks(m):
  A = np.load('../output/Bio_M/model/mc' + str(m) + '_array.npy')

  mc_file = h5py.File('../data/average/cons_locs_pathways_mc6_Column.h5', 'r')
  populations = mc_file.get('populations')
  connections = mc_file.get('connectivity')

  m_type = [
      'L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 'L23_BP',
      'L23_BTC', 'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC',
      'L23_PC', 'L23_SBC', 'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC',
      'L4_MC', 'L4_NBC', 'L4_NGC', 'L4_PC', 'L4_SBC', 'L4_SP', 'L4_SS', 'L5_BP',
      'L5_BTC', 'L5_ChC', 'L5_DBC', 'L5_LBC', 'L5_MC', 'L5_NBC', 'L5_NGC', 'L5_SBC',
      'L5_STPC', 'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC', 'L6_BP', 'L6_BPC', 'L6_BTC',
      'L6_ChC', 'L6_DBC', 'L6_IPC', 'L6_LBC', 'L6_MC', 'L6_NBC', 'L6_NGC', 'L6_SBC',
      'L6_TPC_L1', 'L6_TPC_L4', 'L6_UTPC'
      ]

  array_size_by_morph = []
  for i in m_type:
    array_size_by_morph.append(np.array(connections[i][i]['cMat']).shape[0])

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
      np.save('../output/BC/blocks/mc6_' + str(i) + str(j) + '_array.npy', A[boundaries[i]:boundaries[i+1],boundaries[j]:boundaries[j+1]])



def required_inputs(mc):
    mc_file = h5py.File('../data/average/cons_locs_pathways_mc' + str(mc) + '_Column.h5', 'r')
    populations = mc_file.get('populations')
    connections = mc_file.get('connectivity')
    ##########################################################################################
    # compute full list of neuron locations. Collect array, locate connections (pre/post list)
    locations = np.vstack(np.array(populations[i]['locations']) for i in m_type)
    array = np.load('../output/Bio_M/model/mc' + str(mc) + '_array.npy')
    X = np.where(array == 1)
    distance_list = np.array(np.sqrt(np.sum((locations[X[1]] - locations[X[0]])**2, axis = 1)))
    pre = X[0]
    post = X[1]
    probability1, bins = np.histogram(distance_list, bins = np.arange(0, distance_list.max(), 75))
    np.set_printoptions(suppress=True)
    bins = bins[:-1]**2
    probability = np.array(probability1/max(probability1)).astype(float)[:]
    mean = np.sum(distance_list)/len(distance_list)
    closer_values = np.array(np.where(distance_list < mean)).shape[1]
    minimum = min(distance_list)
    maximum = max(distance_list)
    return pre, post, locations, probability, bins, probability1, mean, closer_values, minimum, maximum


def main_fast(pre, post, locations, probability, model):
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

    print("tries = ", tries)

    arr = np.zeros((31346, 31346), dtype = np.int8)
    for i, j in zip(new_pre, post):
            arr[i,j] = 1
    np.save('../output/' + str(model) + '/model/' + str(model) + '.npy', arr)


def main_fast_geometric(pre, post, locations, probability, bins, model, folder):
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

    for n in tqdm(range(len(post))[:]):
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

            if i != post[n] and random() < probability[search(bins, (loc_post[2]-loc_pre[2])**2 + (loc_post[1] - loc_pre[1])**2 + (loc_post[0] - loc_pre[0])**2) - 1]:
                done = True
                new_pre[n] = i
                mask[index_index] = True
                mask_counter += 1

    print("tries = ", tries)

    arr = np.zeros((31346, 31346), dtype = np.int8)
    for i, j in zip(new_pre, post):
            arr[i,j] = 1
    np.save('../output/GC/model/' + str(model) + '.npy', arr)

def block_configuration(mc):
  M = np.array([])
  for i in tqdm([0, 1, 2, 3, 4]):
    for j in [0, 1, 2, 3, 4]:
      M = np.array(np.load('../output/BC/blocks/mc' + str(mc) + '_' + str(i) + str(j) + '_array.npy'))

      locations = np.vstack(populations[i]['locations'] for i in populations.keys())
      indices = np.where(M == 1)
      pre = indices[0]
      post = indices[1]

      distance_list = np.array(np.sqrt(np.sum((locations[indices[1]] - locations[indices[0]])**2, axis = 1)))
      probability, bins = np.histogram(distance_list, bins = np.arange(0, distance_list.max(), 75))

      np.set_printoptions(suppress=True)
      bins = bins[:-1]**2
      probability = np.array(probability/max(probability)).astype(float)[:]

      def main_fast(pre, post, locations, probability):
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

          for n in range(len(post))[:]:
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
      new_pre = main_fast(pre, post, locations, probability)
      ##########################################################################################
      # Rebuild configured array
      arr = np.zeros((M.shape[0], M.shape[1]), dtype = np.int8)
      for k, l in zip(new_pre, post):
        arr[k,l] += 1
      np.save('../output/BC/new_blocks/mc6_' + str(i) + str(j) + '.npy', arr)

  new_array = np.vstack(np.hstack(np.array(np.load('../output/BC/new_blocks/mc6_' +str(i) + str(j) + '.npy')) for j in [0, 1, 2, 3, 4]) for i in [0, 1, 2, 3, 4])
  np.save('../output/BC/new_blocks/block_configuration.npy', new_array)


def block_geometric_configuration(mc):
  M = np.array([])
  for i in tqdm([0, 1, 2, 3, 4]):
    for j in [0, 1, 2, 3, 4]:
      M = np.array(np.load('../output/BC/blocks/mc' + str(mc) + '_' + str(i) + str(j) + '_array.npy'))

      locations = np.vstack(populations[i]['locations'] for i in populations.keys())
      indices = np.where(M == 1)
      pre = indices[0]
      post = indices[1]

      distance_list = np.array(np.sqrt(np.sum((locations[indices[1]] - locations[indices[0]])**2, axis = 1)))
      probability, bins = np.histogram(distance_list, bins = np.arange(0, distance_list.max(), 75))

      np.set_printoptions(suppress=True)
      bins = bins[:-1]**2
      probability = np.array(probability/max(probability)).astype(float)[:]



      def main_fast(pre, post, locations, probability):
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

          for n in range(len(post)):
              if mask_counter == mask_limit:
                #   print('tries: ', tries, 'time: {:.2f}'.format(time.process_time() - t), "   r_s_list:", mask_limit)
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

                  if i != post[n] and random() < probability[search(bins, (loc_post[2]-loc_pre[2])**2 + (loc_post[1] - loc_pre[1])**2 + (loc_post[0] - loc_pre[0])**2) - 1]:
                      done = True
                      new_pre[n] = i
                      mask[index_index] = True
                      mask_counter += 1

          # print("tries = ", tries)
          return new_pre
      ##########################################################################################
      # t = time.process_time()
      new_pre = main_fast(pre, post, locations, probability)
      # print('runtime: ', time.process_time() - t)
      ##########################################################################################
      # Rebuild configured array
      arr = np.zeros((M.shape[0], M.shape[1]), dtype = np.int8)
      for k, l in zip(new_pre, post):
        arr[k,l] += 1
      np.save('../output/BGC/new_blocks/mc6_' + str(i) + str(j) + '.npy', arr)

  new_array = np.vstack(np.hstack(np.array(np.load('../output/BGC/new_blocks/mc6_' +str(i) + str(j) + '.npy')) for j in [0, 1, 2, 3, 4]) for i in [0, 1, 2, 3, 4])
  np.save('../output/BGC/new_blocks/block_geometric_configuration.npy', new_array)


#------------------------------------------------------------------------------------------#
springs = cm.get_cmap('cool', 256)
pinks = springs(np.linspace(0.4, 1, 256))
pnkcolor = ListedColormap(pinks)

def plot_examples(cms):
    """
    helper function to plot two colormaps
    """
    np.random.seed(19680801)
    data = np.random.randn(30, 30)

    fig, axs = plt.subplots(1, 2, figsize=(6, 3), constrained_layout=True)
    for [ax, cmap] in zip(axs, cms):
        psm = ax.pcolormesh(data, cmap=cmap, rasterized=True, vmin=-4, vmax=4)
        fig.colorbar(psm, ax=ax)
    plt.show()

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
  plt.title('Adjacency matrix for the ' + str(model) +' MC')
  plt.xlabel('Post Synaptic Neurons')
  plt.ylabel('Pre Synaptic Neurons')
  plt.xticks([0, 15673, 31346])
  plt.yticks([0, 15673, 31346], rotation = 45)
  plt.subplot().yaxis.set_label_position("right")
  plt.subplot().yaxis.tick_right()
  plt.axes([0, 0, 1, 1])
  plt.axis('off')
  plt.savefig('../output/plots/' + str(folder) + '/matrix_' + str(model) + '.png')
  plt.clf()

# shows same matrix with layers
def layered_graph(matrix, model, folder):
  plt.imshow(matrix, cmap = 'Greys', interpolation = 'nearest')
  plt.title('Adjacency matric Representing Connections')
  plt.axhline(339, color = 'black')
  plt.axhline(7857, color = 'black')
  plt.axhline(12518, color = 'black')
  plt.axhline(18624, color = 'black')
  plt.axvline(339, color = 'black')
  plt.axvline(7857, color = 'black')
  plt.axvline(12518, color = 'black')
  plt.axvline(18624, color = 'black')
  plt.xlabel('Post-Synaptic Neurons')
  plt.ylabel('Pre-Synaptic Neurons')
  plt.xticks([170, 4098, 10052, 15571, 24985], ['L1', 'L23', 'L4', 'L5', 'L6'])
  plt.yticks([170, 4098, 10052, 15571, 24985], ['L1', 'L23', 'L4', 'L5', 'L6'])
  plt.axes([0, 0, 1, 1])
  plt.axis('off')
  plt.savefig('../output/plots/' + str(folder) + '/matrix_' +  str(model) +'_layered.png')
  plt.clf()

def layer_counts(matrix):
    return np.sum(matrix[:339,:339]), np.sum(matrix[:339,339:7587]) ,np.sum(matrix[:339, 7587:12518]) ,np.sum(matrix[:339,12518:18624]) ,np.sum(matrix[:339, 18624:31346]) ,np.sum(matrix[339:7587, :339]), np.sum(matrix[339:7587, 339:7587]), np.sum(matrix[339:7587, 7587:12518]), np.sum(matrix[339:7587, 12518:18624]), np.sum(matrix[339:7587, 18624:31346]), np.sum(matrix[7587:12518, :339]), np.sum(matrix[7587:12518, 339:7587]), np.sum(matrix[7587:12518, 7587:12518]), np.sum(matrix[7587:12518, 12518:18624]), np.sum(matrix[7587:12518, 18624:31346]), np.sum(matrix[12518:18624, :339]), np.sum(matrix[12518:18624, 339:7587]), np.sum(matrix[12518:18624, 7587:12518]), np.sum(matrix[12518:18624, 12518:18624]), np.sum(matrix[12518:18624, 18624:31346]), np.sum(matrix[18624:31346, :339]), np.sum(matrix[18624:31346, 339:7587]), np.sum(matrix[18624:31346, 7587:12518]), np.sum(matrix[18624:31346, 12518:18624]), np.sum(matrix[18624:31346, 18624:31346])



def layer_densities(matrix):
  return np.sum(matrix[:339,:339])/(339**2), np.sum(matrix[:339,339:7587])/(339*7248) ,np.sum(matrix[:339, 7587:12518])/(339*4932), np.sum(matrix[:339,12518:18624])/(339*6106), np.sum(matrix[:339, 18624:31346])/(339*12722), np.sum(matrix[339:7587, :339])/(7248*339), np.sum(matrix[339:7587, 339:7587])/(7248**2), np.sum(matrix[339:7587, 7587:12518])/(7248*4932), np.sum(matrix[339:7587, 12518:18624])/(7248*6106), np.sum(matrix[339:7587, 18624:31346])/(7248*12722), np.sum(matrix[7587:12518, :339])/(4932*339), np.sum(matrix[7587:12518, 339:7587])/(4932*7248), np.sum(matrix[7587:12518, 7587:12518])/(4932**2), np.sum(matrix[7587:12518, 12518:18624])/(4932*6106), np.sum(matrix[7587:12518, 18624:31346])/(4932*12722), np.sum(matrix[12518:18624, :339])/(6106*339), np.sum(matrix[12518:18624, 339:7587])/(6106*7248), np.sum(matrix[12518:18624, 7587:12518])/(6106*4932), np.sum(matrix[12518:18624, 12518:18624])/(6106**2), np.sum(matrix[12518:18624, 18624:31346])/(6106*12722), np.sum(matrix[18624:31346, :339])/(12722*339), np.sum(matrix[18624:31346, 339:7587])/(12722*7248), np.sum(matrix[18624:31346, 7587:12518])/(12722*4932), np.sum(matrix[18624:31346, 12518:18624])/(12722*6106), np.sum(matrix[18624:31346, 18624:31346])/(12722**2)

def heat_map(layers, model, folder):
  size = 5
  data = np.arange(size * size).reshape((size, size))
  values = np.random.rand(size, size)

  # Limits for the extent
  x_start = -0.25
  x_end = 4.8
  y_start = -0.25
  y_end = 4.8

  extent = [x_start, x_end, y_start, y_end]

  # The normal figure
  fig = plt.figure(figsize=(16, 12))
  ax = fig.add_subplot(111)
  im = ax.imshow(layers, cmap='Oranges')

  # Add the text
  plt.rcParams.update({'font.size': 22})
  jump_x = (x_end - x_start) / (5 * size)
  jump_y = (y_end - y_start) / (5 * size)
  x_positions = np.linspace(start=x_start, stop=x_end, num=size, endpoint=False)
  y_positions = np.linspace(start=y_start, stop=y_end, num=size, endpoint=False)

  for y_index, y in enumerate(y_positions):
      for x_index, x in enumerate(x_positions):
          label = layers[y_index, x_index]
          text_x = x + jump_x
          text_y = y + jump_y
          float_formatter = "{:.4f}".format
          np.set_printoptions(formatter={'float_kind':float_formatter})


          ax.text(text_x, text_y, label, color='green', ha='center', va='center')

  plt.title(str(model) + ' of connections')
  plt.xticks([0, 1, 2, 3, 4], ['L1', 'L23', 'L4', 'L5', 'L6'])
  plt.yticks([0, 1, 2, 3, 4], ['L1', 'L23', 'L4', 'L5', 'L6'])
  plt.tick_params(labelsize = 16)
  plt.xlabel('Post Synaptic Layers', fontsize = 20)
  plt.ylabel('Pre Synaptic Layers', fontsize = 20)
  fig.colorbar(im)
  plt.savefig('../output/plots/' + str(folder) + '/heat_map_layer_' + str(model) + '.png')
  plt.clf()


def degree(matrix, model, folder):
  sd = np.sum(matrix, axis = 1) - np.sum(matrix, axis = 0)
  plt.bar(np.arange(31346), height = sd, color = pnkcolor(1.0))
  plt.title('Signed Degree for each Neuron for ' + str(model) + ' Model')
  plt.xlabel('Neuron')
  plt.ylabel('Signed Degree')
  plt.savefig('../output/plots/' + str(folder) + '/' + str(model) + '_sd.png')
  plt.clf()

def cum_degree(matrix, model, folder):
  sd = np.sum(matrix, axis = 1) - np.sum(matrix, axis = 0)
  cum_sd = np.cumsum(sd)
  plt.bar(np.arange(31346), cum_sd, color = pnkcolor(1.0))
  plt.title('Cumulative degree for all Neurons for ' + str(model) + ' Model')
  plt.xlabel('Neuron')
  plt.ylabel('Cumulative Signed Degree')
  plt.savefig('../output/plots/' + str(folder) + '/cumsum_degree_' +str(model) + '.png')
  plt.clf()



def overall_degree(matrix, model, folder):
  sd = np.sum(matrix, axis = 1) - np.sum(matrix, axis = 0)
  sd_2 = np.cumsum(sd**2)
  plt.bar(np.arange(31346), sd_2, color = pnkcolor(1.0))
  plt.title('Directionality for ' + str(model) + ' Model')
  plt.xlabel('Neuron')
  plt.ylabel('Cumulative Signed Degree')
  plt.savefig('../output/plots/' + str(folder) + '/overall_degree_' +str(model) + '_sd_2.png')
  plt.clf()



def topological_stats(array, a, b):
  return pyflagser.flagser_count_unweighted(array,directed=True), pyflagser.flagser_unweighted(array, min_dimension=a, max_dimension=b, directed=True)['betti'], np.sum(pyflagser.flagser_count_unweighted(array, directed=True)[0::2]) - np.sum(pyflagser.flagser_count_unweighted(array, directed=True)[1::2])

def distance_distribution(locations, old_array, new_array, model, folder):
  locations = np.vstack(np.array(populations[i]['locations']) for i in m_type)

  pre = np.where(old_array == 1)[0]
  post = np.where(old_array == 1)[1]
  new_pre = np.where(new_array == 1)[0]
  new_post = np.where(new_array == 1)[1]

  distance_lists = np.array(np.sqrt(np.sum((locations[post] - locations[pre])**2, axis = 1)))
  new_distance_lists = np.array(np.sqrt(np.sum((locations[new_post] - locations[new_pre])**2, axis = 1)))

  mean = np.sum(new_distance_lists)/len(new_distance_lists)
  minimum = min(new_distance_lists)
  maximum = max(new_distance_lists)
  lower_half = np.array(np.where(new_distance_lists < mean)).shape[1]

  # Arrange data to be in form of histogram
  probability, bins = np.histogram(distance_lists, np.arange(0, distance_lists.max(), 75))
  new_prob, new_bin = np.histogram(new_distance_lists, np.arange(0, new_distance_lists.max(), 75))
  ##########################################################################################
  # Plot data
  # plt.facecolor('white')
  plt.hist(bins[:-1], bins, weights = probability, label = 'Bio-M', alpha = 0.5)
  plt.hist(new_bin[:-1], new_bin, weights = new_prob, label = str(model) + ' model', alpha = 0.5, color = 'orange')
  plt.xlabel('distance between connected neurons ($\mu m$)')
  plt.ylabel('count')
  plt.title('N2N distance for Bio M model Vs ' + str(model) + ' model')
  plt.legend()
  plt.savefig('../output/plots/' + str(folder) + '/' + str(model) + '_dist_distr.png')
  plt.clf()
  return new_prob, mean, minimum, maximum, lower_half


def original_distance_distribution(locations, old_array, new_array, model, folder):
  locations = np.vstack(np.array(populations[i]['locations']) for i in m_type)

  pre = np.where(old_array == 1)[0]
  post = np.where(old_array == 1)[1]

  distance_lists = np.array(np.sqrt(np.sum((locations[post] - locations[pre])**2, axis = 1)))

  # Arrange data to be in form of histogram
  probability, bins = np.histogram(distance_lists, np.arange(0, distance_lists.max(), 75))
  ##########################################################################################
  # Plot data
  plt.hist(bins[:-1], bins, weights = probability, label = 'Bio-M', alpha = 0.5, color = 'deeppink')
  plt.xlabel('distance between connected neurons ($\mu m$)')
  plt.ylabel('count')
  plt.title('N2N distance for Bio M Connectome')
  plt.legend()
  plt.savefig('../output/plots/' + str(folder) + '/' + str(model) + '_dist_distr.png')
  plt.clf()


def morph_counts(matrix):
  morphs = [
    0, 58, 23, 91, 72, 52, 43, 28, 105, 61, 176, 456, 332, 267, 55, 5870, 168, 8, 20, 8, 38,
    121, 118, 94, 6, 2682, 62, 1097, 407, 34, 76, 19, 96, 208, 398, 203, 8, 25, 301, 2398,
    1996, 344, 70, 3182, 55, 16, 31, 3461, 460, 335, 200, 17, 67, 1641, 1450, 1737
    ]

  morphs_sum = np.cumsum(morphs)
  counts = []

  for i in range(len(morphs_sum)-1)[:]:
    for j in range(len(morphs_sum)-1)[:]:
      counts.append(len(np.where(matrix[morphs_sum[i]:morphs_sum[i+1], morphs_sum[j]:morphs_sum[j+1]] == 1)[0]))
  counts = np.nan_to_num(np.log(counts))
  counts[counts < 0] = 0
  return counts




def heat_map_morph(morphs, model, folder):
  size = 55
  data = np.arange(size * size).reshape((size, size))
  values = np.random.rand(size, size)

  # Limits for the extent
  x_start = -0.25
  x_end = 55
  y_start = -0.25
  y_end = 55

  extent = [x_start, x_end, y_start, y_end]

  # The normal figure
  fig = plt.figure(figsize=(16, 12))
  ax = fig.add_subplot(111)
  im = ax.imshow(morphs, cmap='Wistia')

  plt.title('Connections counts per morphological type (' +str(model) + ')', fontsize=24)
  plt.xticks(range(len(m_type)), m_type, rotation=90)
  plt.yticks(range(len(m_type)), m_type)
  plt.tick_params(labelsize = 10)
  plt.xlabel('Post Synaptic morphological types', fontsize = 20)
  plt.ylabel('Pre Synaptic morphological types', fontsize = 20)
  fig.colorbar(im)
  plt.savefig('../output/plots/' + str(folder) + '/heat_map_morph_layer_' + str(model) + '.png')
  plt.clf()
