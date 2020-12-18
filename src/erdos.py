import numpy as np
import pyflagser
from scipy.sparse import rand

# probability generation function
# def f():
#     if np.random.rand() < 0.008:
#         return 1
#     else:
#         return 0

###############################################################
N = 31346
A = np.zeros((N, N), dtype = np.int8)

# Memory consumption version (generates matrix in one go)
# A[:,:] = np.random.rand(N,N) < 0.1

# row by row, less memory consuimg but faster than element by element
for i in range(N):
    A[i,:] = np.random.rand(N) < 0.008

# Element by element, slowest
# for i in range(N):
#     for j in range(N):
#         A[i,j] = np.random.rand() < 0.008
#         A[i, j] = f() # Alternative to above line


# Builtin function version
# x = rand(N, N, density=0.008, format='csr')
# x.data[:] = 1
###############################################################
np.save("a", A)
###############################################################
M = np.load("a.npy")
print(pyflagser.flagser_count_unweighted(M, min_dimension=0, max_dimension = 7, directed=True))
###############################################################
'''
python3 erdos.py >> erdos.txt
'''









