import numpy as np
import pyflagser
from scipy.sparse import rand

###############################################################
N = 31346
A = np.zeros((N, N), dtype = np.int8)

for i in range(N):
    A[i,:] = np.random.rand(N) < 0.008
###############################################################
np.save("a", A)
###############################################################
M = np.load("a.npy")
print(pyflagser.flagser_count_unweighted(M, min_dimension=0, max_dimension = 7, directed=True))
###############################################################
