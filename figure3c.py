import h5py
import pandas as pd
import numpy as np
from pyflagsercontain import flagser_count

#############################################################################################
# difficult to reproduce as it has been mentioned that the Control files are in relation to cloud model connectivity
# general idea is to take all neurons in each layer and take the highest simplicial complex in which it belongs to and take the quartiles for each layer
#############################################################################################
mc0 = h5py.File('../../pathway_average_files/cons_locs_pathways_mc6_Column.h5', 'r')
connectivity = mc0.get('connectivity')
x = np.vstack([np.hstack([connectivity[i][j]['cMat'] for j in connectivity]) for i in connectivity])

simplex_counts = flagser_count(x)
assert len(simplex_counts) == 31346

full_matrix = [i for i in range(31346)]
a_Full =[]

for i in full_matrix:
    a_Full.append(simplex_counts[i])

a_ = pd.DataFrame(a_Full).fillna(0)

L1_indices = a_[0:340]
L23_indices = a_[341:7859]
L4_indices = a_[7860:12519]
L5_indices  =a_[12520:18627]
L6_indices = a_[18628:31346]

# Below computed for each layer
L6_indices.columns = ['Dim0', 'Dim1', 'Dim2', 'Dim3', 'Dim4', 'Dim5', 'Dim6']

d = L6_indices.mask(L6_indices == 0)
tets = L6_indices.assign(
    Start = d.apply(pd.Series.first_valid_index, 1),
    Finish = d.apply(pd.Series.last_valid_index, 1))

e = pd.DataFrame(tets[['Finish']])
f = np.array(e)
f = [[item.replace('Dim', '') for item in lst] for lst in f]
g = pd.DataFrame(f)

g.columns = ['Dim+']
g['Dim+'] = g['Dim+'].astype(int)

numerator = (g['Dim+'].sum(axis = 0)).astype(float)
h = g['Dim+'].value_counts()
denominator = (h.sum(axis = 0)).astype(float)

average = numerator/denominator
print(average)
