import numpy as np
import time
import pyflagser
import sys

try:
    mc = sys.argv[1]
    instance = sys.argv[2]
    mtype = sys.argv[3]
except IndexError:
    raise SystemExit(f"usage: {sys.argv[0]} <mc> <instance> <mtype>")


array =  np.load('../output/array_mc' + str(mc) + '_' + str(mtype) + '.npy')
X = np.where(array == 1)
pre = np.array(np.load('../output/instance_' + str(instance) + '_mc' + str(mc) + '_' + str(mtype) + '.npy'))
post = X[1]

C = np.zeros(len(pre)*2)
C[0::2] = pre
C[1::2] = post
C = C.astype(int)
arr = np.zeros((len(array), len(array)), dtype = np.int8)
idx_list = C.reshape((len(pre), 2))

for i in idx_list:
    if arr[i[0], i[1]] == 1:
        continue
    else:
        arr[i[0], i[1]] = 1

print("Number of simplices: ", pyflagser.flagser_count_unweighted(arr, directed=True))
print("Betti counts: ", pyflagser.flagser_unweighted(arr, min_dimension=5, max_dimension=5, directed=True))
