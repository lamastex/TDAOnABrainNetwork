import numpy as np
import h5py, time, random, os, sys, pyflagser, warnings, datetime
import matplotlib.pyplot as plt
from tqdm import tqdm

t = time.process_time()



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
    probability, bins = np.histogram(distance_list, bins = np.arange(0, distance_list.max(), 75))
    np.set_printoptions(suppress=True)
    bins = bins[:-1]**2
    probability = np.array(probability/max(probability)).astype(float)[:]
    return pre, post, locations, probability, bins










def main_fast(pre, post, locations, probability, model):
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

    for n in tqdm(range(len(post))):
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

    print("tries = ", tries)

    arr = np.zeros((31346, 31346), dtype = np.int8)
    for i, j in zip(new_pre, post):
            arr[i,j] = 1
    np.save('../output/' + str(model) + '/model/' + str(model) + '.npy', arr)


def main_fast_geometric(pre, post, locations, probability, bins, model, folder):
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

    for n in tqdm(range(len(post))[:]):
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

    print("tries = ", tries)

    arr = np.zeros((31346, 31346), dtype = np.int8)
    for i, j in zip(new_pre, post):
            arr[i,j] = 1
    np.save('../output/GC/model/' + str(model) + '.npy', arr)



def block_configuration(mc):
  M = np.array([])
  for i in [1, 2, 4, 5, 6]:
    for j in [1, 2, 4, 5, 6]:
      print(i, j)
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
          removal_service_list = np.array(range(len(pre)))
          mask = np.zeros(len(removal_service_list), dtype = bool)
          mask_counter = 0
          mask_limit = len(mask)//2
          randint = np.random.randint
          random = np.random.random
          search = np.searchsorted
          N = len(removal_service_list)
          tries = 0

          for n in tqdm(range(len(post))[:]):
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

          print("tries = ", tries)
          return new_pre
      ##########################################################################################
      t = time.process_time()
      new_pre = main_fast(pre, post, locations, probability)
      print('runtime: ', time.process_time() - t)
      ##########################################################################################
      # Rebuild configured array
      arr = np.zeros((M.shape[0], M.shape[1]), dtype = np.int8)
      for k, l in zip(new_pre, post):
        arr[k,l] = 1
      np.save('../output/BC/new_blocks/mc6_' + str(i) + str(j) + '.npy', arr)

  new_array = np.vstack(np.hstack(np.array(np.load('../output/BC/new_blocks/mc6_' +str(i) + str(j) + '.npy')) for j in [1, 2, 4, 5, 6]) for i in [1, 2, 4, 5, 6])
  np.save('../output/BC/new_blocks/block_configuration.npy', new_array)


def geometric_block_configuration(mc):
  M = np.array([])
  for i in [1, 2, 4, 5, 6]:
    for j in [1, 2, 4, 5, 6]:
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
          removal_service_list = np.array(range(len(pre)))
          mask = np.zeros(len(removal_service_list), dtype = bool)
          mask_counter = 0
          mask_limit = len(mask)//2
          randint = np.random.randint
          random = np.random.random
          search = np.searchsorted
          N = len(removal_service_list)
          tries = 0

          for n in tqdm(range(len(post))):
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

          print("tries = ", tries)
          return new_pre
      ##########################################################################################
      t = time.process_time()
      new_pre = main_fast(pre, post, locations, probability)
      print('runtime: ', time.process_time() - t)
      ##########################################################################################
      # Rebuild configured array
      arr = np.zeros((M.shape[0], M.shape[1]), dtype = np.int8)
      for k, l in zip(new_pre, post):
        arr[k,l] = 1
      np.save('../output/GBC/new_blocks/mc6_' + str(i) + str(j) + '.npy', arr)

  new_array = np.vstack(np.hstack(np.array(np.load('../output/GBC/new_blocks/mc6_' +str(i) + str(j) + '.npy')) for j in [1, 2, 4, 5, 6]) for i in [1, 2, 4, 5, 6])
  np.save('../output/GBC/new_blocks/geometric_block_configuration.npy', new_array)


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
  plt.title('Adjacency matrix for the ' + str(model) +' model')
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
    return len(np.where(matrix[:339,:339] == 1)[0]), len(np.where(matrix[:339,339:7587] == 1)[0]) ,len(np.where(matrix[:339, 7587:12518] == 1)[0]) ,len(np.where(matrix[:339,12518:18624] == 1)[0]) ,len(np.where(matrix[:339, 18624:31346] == 1)[0]) ,len(np.where(matrix[339:7587, :339] == 1)[0]) ,len(np.where(matrix[339:7587, 339:7587] == 1)[0]) ,len(np.where(matrix[339:7587, 7587:12518] == 1)[0]) ,len(np.where(matrix[339:7587, 12518:18624] == 1)[0]) ,len(np.where(matrix[339:7587, 18624:31346] == 1)[0]) ,len(np.where(matrix[7587:12518, :339] == 1)[0]) ,len(np.where(matrix[7587:12518, 339:7587] == 1)[0]) ,len(np.where(matrix[7587:12518, 7587:12518] == 1)[0]) ,len(np.where(matrix[7587:12518, 12518:18624] == 1)[0]) ,len(np.where(matrix[7587:12518, 18624:31346] == 1)[0]) ,len(np.where(matrix[12518:18624, :339] == 1)[0]) ,len(np.where(matrix[12518:18624, 339:7587] == 1)[0]) ,len(np.where(matrix[12518:18624, 7587:12518] == 1)[0]) ,len(np.where(matrix[12518:18624, 12518:18624] == 1)[0]) ,len(np.where(matrix[12518:18624, 18624:31346] == 1)[0]) ,len(np.where(matrix[18624:31346, :339] == 1)[0]) ,len(np.where(matrix[18624:31346, 339:7587] == 1)[0]) ,len(np.where(matrix[18624:31346, 7587:12518] == 1)[0]) ,len(np.where(matrix[18624:31346, 12518:18624] == 1)[0]) ,len(np.where(matrix[18624:31346, 18624:31346] == 1)[0])


def layer_densities(matrix):
  return len(np.where(matrix[:339,:339] == 1)[0])/(339**2), len(np.where(matrix[:339,339:7587] == 1)[0])/(339*7248) ,len(np.where(matrix[:339, 7587:12518] == 1)[0])/(339*4932), len(np.where(matrix[:339,12518:18624] == 1)[0])/(339*6106), len(np.where(matrix[:339, 18624:31346] == 1)[0])/(339*12722), len(np.where(matrix[339:7587, :339] == 1)[0])/(7248*339), len(np.where(matrix[339:7587, 339:7587] == 1)[0])/(7248**2), len(np.where(matrix[339:7587, 7587:12518] == 1)[0])/(7248*4932), len(np.where(matrix[339:7587, 12518:18624] == 1)[0])/(7248*6106), len(np.where(matrix[339:7587, 18624:31346] == 1)[0])/(7248*12722), len(np.where(matrix[7587:12518, :339] == 1)[0])/(4932*339), len(np.where(matrix[7587:12518, 339:7587] == 1)[0])/(4932*7248), len(np.where(matrix[7587:12518, 7587:12518] == 1)[0])/(4932**2), len(np.where(matrix[7587:12518, 12518:18624] == 1)[0])/(4932*6106), len(np.where(matrix[7587:12518, 18624:31346] == 1)[0])/(4932*12722), len(np.where(matrix[12518:18624, :339] == 1)[0])/(6106*339), len(np.where(matrix[12518:18624, 339:7587] == 1)[0])/(6106*7248), len(np.where(matrix[12518:18624, 7587:12518] == 1)[0])/(6106*4932), len(np.where(matrix[12518:18624, 12518:18624] == 1)[0])/(6106**2), len(np.where(matrix[12518:18624, 18624:31346] == 1)[0])/(6106*12722), len(np.where(matrix[18624:31346, :339] == 1)[0])/(12722*339), len(np.where(matrix[18624:31346, 339:7587] == 1)[0])/(12722*7248), len(np.where(matrix[18624:31346, 7587:12518] == 1)[0])/(12722*4932), len(np.where(matrix[18624:31346, 12518:18624] == 1)[0])/(12722*6106), len(np.where(matrix[18624:31346, 18624:31346] == 1)[0])/(12722**2)

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
  im = ax.imshow(layers, cmap='Blues')

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


          ax.text(text_x, text_y, label, color='orange', ha='center', va='center')

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
  plt.bar(np.arange(31346), height = sd)
  plt.title('Signed Degree for each Neuron for ' + str(model) + ' Model')
  plt.xlabel('Neuron')
  plt.ylabel('Signed Degree')
  plt.savefig('../output/plots/' + str(folder) + '/' + str(model) + '_sd.png')
  plt.clf()


def overall_degree(matrix, model, folder):
  sd = np.sum(matrix, axis = 1) - np.sum(matrix, axis = 0)
  sd_2 = np.cumsum(sd**2)
  plt.bar(np.arange(31346), sd_2)
  plt.title('Cumulative Directionality for all Neurons for ' + str(model) + ' Model')
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

  # Arrange data to be in form of histogram
  probability, bins = np.histogram(distance_lists, np.arange(0, distance_lists.max(), 75))
  new_prob, new_bin = np.histogram(new_distance_lists, np.arange(0, new_distance_lists.max(), 75))
  ##########################################################################################
  # Plot data
  plt.hist(bins[:-1], bins, weights = probability, label = 'Bio-M', alpha = 0.5)
  plt.hist(new_bin[:-1], new_bin, weights = new_prob, label = str(model), alpha = 0.5)
  plt.xlabel('distance between connected neurons ($\mu m$)')
  plt.ylabel('count')
  plt.title('N2N distance for Bio M Connectome Vs ' + str(model))
  plt.legend()
  plt.savefig('../output/plots/' + str(folder) + '/' + str(model) + '_dist_distr.png')
  plt.clf()


def original_distance_distribution(locations, old_array, new_array, model, folder):
  locations = np.vstack(np.array(populations[i]['locations']) for i in m_type)

  pre = np.where(old_array == 1)[0]
  post = np.where(old_array == 1)[1]

  distance_lists = np.array(np.sqrt(np.sum((locations[post] - locations[pre])**2, axis = 1)))

  # Arrange data to be in form of histogram
  probability, bins = np.histogram(distance_lists, np.arange(0, distance_lists.max(), 75))
  ##########################################################################################
  # Plot data
  plt.hist(bins[:-1], bins, weights = probability, label = 'Bio-M', alpha = 0.5)
  plt.xlabel('distance between connected neurons ($\mu m$)')
  plt.ylabel('count')
  plt.title('N2N distance for Bio M Connectome')
  plt.legend()
  plt.savefig('../output/plots/' + str(folder) + '/' + str(model) + '_dist_distr.png')
  plt.clf()
