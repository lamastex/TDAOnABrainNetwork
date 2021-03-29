import h5py, warnings, bisect, time, datetime, pyflagser, sys
import numpy as np
import matplotlib.pyplot as plt
from m_types import *

warnings.filterwarnings('ignore')
print(datetime.datetime.now())
t = time.process_time()

try:
    mc = sys.argv[1]
    start = int(sys.argv[2])
    stop = int(sys.argv[3])
except IndexError:
    raise SystemExit(f"Usage: {sys.argv[0]} <mc> <start> <stop>")

for z in range(start, stop):
  print(z)
  mc_file = h5py.File('../data/average/cons_locs_pathways_mc' + str(mc) + '_Column.h5', 'r')
  populations = mc_file.get('populations')
  connections = mc_file.get('connectivity')

  ref_list = [[layer1, layer2, layer4, layer5, layer6]]
  ##########################################################################################
  # size_list computes m_type sizes, layer_size computes layer size
  size_list = []
  for i in ref_list:
    for j in i:
      for k in j:
        size_list.append(np.array(connections[k][k]['cMat'], dtype = np.int8).shape[0])

  layer_size = []
  for i in ref_list:
    for j in i:
      layer_size.append(len(j))

  layer_list_end = np.cumsum(layer_size)
  layer_list_start = np.insert(layer_list_end, 0, 0)[:-1]

  the_layers = []
  for i, j in zip(layer_list_start, layer_list_end):
    the_layers.append(np.sum(size_list[i:j]))

  ###########################################################################################
  M = np.zeros((sum(size_list), sum(size_list)), dtype = np.int8)
  ###########################################################################################
  # for M_a, i in zip(m_type, range(len(m_type))):
  #     for M_b, j in zip(m_type, range(len(m_type))):
  #         M_ij = np.array(connections[M_a][M_b]['cMat'], dtype = np.int8)
  #         k_start = np.array(size_list)[:i].sum()
  #         l_start = np.array(size_list)[:j].sum()
  #         for k in range(size_list[i]):
  #             for l in range(size_list[j]):
  #                 M[k_start + k,l_start + l] = M_ij[k,l]
  ###########################################################################################
  # np.save('../output/test_array.npy', M)
  M = np.load('../output/test_array.npy')

  # Have layer arrays split as so
  thing_of_arrays = []
  for i in the_layers:
    for j in the_layers:
      thing_of_arrays.append(M[:i, :j])

  # gather locations for each layer
  locations = np.vstack(populations[i]['locations'] for i in populations.keys())
  location_list = []
  for i in the_layers:
    location_list.append(locations[:i])

  # acquire indices of pre and post for each array
  index_lists = []
  for i in range(len(thing_of_arrays)):
    index_lists.append(np.where(thing_of_arrays[i] == 1))

  # compute the distances of each pairwise neuron set for each array
  arrays_of_distances = []
  for i, j, k in zip(np.arange(25), np.repeat(np.arange(5), 5), np.tile(np.arange(5), 5)):
    arrays_of_distances.append(np.array(np.sqrt(np.sum((location_list[k][index_lists[i][1]] - location_list[j][index_lists[i][0]])**2, axis = 1))))

  # compute probabilistic connections within each array
  prob_list = []
  bin_list =[]
  for i in arrays_of_distances[:]:
    probability, bins = np.histogram(i, np.arange(0, i.max(), 100))
    prob_list.append(probability/max(probability))
    bin_list.append(bins)

  def boundary(num1, breakpoints = [339, 7857, 12518, 18624, 31346], result = (['L1', 'L23', 'L4', 'L5', 'L6'])):
    i = bisect.bisect(breakpoints, num1-1)
    return result[i]

  d = {
    'L1':{
      'L1':(prob_list[0], bin_list[0][:-1]**2), 'L23':(prob_list[1], bin_list[1][:-1]**2), 'L4':(prob_list[2], bin_list[2][:-1]**2), 'L5':(prob_list[3], bin_list[3][:-1]**2), 'L6':(prob_list[4], bin_list[4][:-1]**2)
    },
    'L23':{
      'L1':(prob_list[5][:], bin_list[5][:-1]**2), 'L23':(prob_list[6][:], bin_list[6][:-1]**2), 'L4':(prob_list[7][:-1], bin_list[7][:-2]**2), 'L5':(prob_list[8][:], bin_list[8][:-1]**2), 'L6':(prob_list[9][:9], bin_list[9][:-10]**2)
    },
    'L4':{
      'L1':(prob_list[10][:], bin_list[10][:-1]**2), 'L23':(prob_list[11][:], bin_list[11][:-1]**2), 'L4':(prob_list[12], bin_list[12][:-1]**2), 'L5':(prob_list[13][:-1], bin_list[13][:-2]**2), 'L6':(prob_list[14][:-6], bin_list[14][:-7]**2)
    },
    'L5':{
      'L1':(prob_list[15][:], bin_list[15][:-1]**2), 'L23':(prob_list[16][:], bin_list[16][:-1]**2), 'L4':(prob_list[17][:], bin_list[17][:-1]**2), 'L5':(prob_list[18][:], bin_list[18][:-1]**2), 'L6':(prob_list[19][:-5], bin_list[19][:-6]**2)
    },
    'L6':{
      'L1':(prob_list[20][:-6], bin_list[20][:-7]**2), 'L23':(prob_list[21][:-4], bin_list[21][:-5]**2), 'L4':(prob_list[22][:-4], bin_list[22][:-5]**2), 'L5':(prob_list[23][:-4], bin_list[23][:-5]**2), 'L6':(prob_list[24][:-4], bin_list[24][:-5]**2)
    }
  }

  X = np.where(M == 1)
  pre = X[0]
  post = X[1]
  locations = np.vstack(np.array(populations[i]['locations']) for i in m_type)

  def main_fast(pre, post, locations, d):
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
                  print('tries: ', tries, 'time: {:.2f}'.format(time.process_time() - t), "   r_s_list:", mask_limit)
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

                  if i != post[n] and random() < d[boundary(i)][boundary(post[n])][0][np.searchsorted(d[boundary(i)][boundary(post[n])][1], (loc_post[2]- loc_pre[2])**2 + (loc_post[1] - loc_pre[1])**2 + (loc_post[0] - loc_pre[0])**2)-1]:

                      done = True
                      new_pre[n] = i
                      mask[index_index] = True
                      mask_counter += 1

          print("tries = ", tries)
          return new_pre

  new_pre = main_fast(pre, post, locations, d)

  arr = np.zeros((len(M), len(M)), dtype = np.int8)
  for i, j in zip(new_pre, post):
      arr[i,j] = 1

  simplices = pyflagser.flagser_count_unweighted(arr, directed=True)
  betti_numbers = pyflagser.flagser_unweighted(arr, min_dimension = 5, max_dimension = 6, directed = True)['betti']
  euler_char_simplices = np.sum(pyflagser.flagser_count_unweighted(arr, directed=True)[0::2]) - np.sum(pyflagser.flagser_count_unweighted(arr, directed=True)[1::2])
