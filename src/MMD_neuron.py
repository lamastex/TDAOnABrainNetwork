'''
Script reproduces the mean maximal dimension of all neurons in each layer of the Bio-M model
(Figure 3c of Cliques article)

Connectivity matrix is read in from the downloadable files from BBP website.
The number and dimension of simplices to which each neuron belongs to are detailed in a_.

We then require to know which layer each neuron belongs to, so we collect the m_types from the 
connectivity matrix and from there we can ascertain as to how large each matrix is, thereby 
knowing how many neurons are contained within. Once we knew this, we could append a column
labelling which m_type each neuron belonged to. Further to this, since we require details of
layer only, we stripped everything away that didn't specifically state which layer the neuron
belonged to. We label the layer column 'l_type'.

Finally, since we want to know the column or dimension to which contained the highest non-zero
value, i.e the highest dimension to which the neuron belonged to, we had to strip away the 
columns 'n' and 'l_type' but not before we grouped all neurons into their respective layers.

Since we want to find the mean maximal dimension for all layers, we loop through each layer
to obtain our results and have this printed out at the bottom

NOTE:
run using 'python MMD_neuron.py'
'''
import h5py
import pandas as pd
import numpy as np
from pyflagsercontain import flagser_count
#############################################################################################
mc0 = h5py.File('../../pathway_average_files/cons_locs_pathways_mc6_Column.h5', 'r')
connectivity = mc0.get('connectivity')
x = np.vstack([np.hstack([connectivity[i][j]['cMat'] for j in connectivity]) for i in connectivity])
#############################################################################################
simplex_counts = flagser_count(x)
assert len(simplex_counts) == 31346
a_Full = []

for i in [i for i in range(31346)]:
    a_Full.append(simplex_counts[i])
a_ = pd.DataFrame(a_Full).fillna(0)
#############################################################################################
m_type = [
    'L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 'L23_BP', 
    'L23_BTC', 'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC', 
    'L23_PC', 'L23_SBC', 'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC', 'L4_MC', 
    'L4_NBC', 'L4_NGC', 'L4_PC', 'L4_SBC', 'L4_SP', 'L4_SS', 'L5_BP', 'L5_BTC', 
    'L5_ChC', 'L5_DBC', 'L5_LBC', 'L5_MC', 'L5_NBC', 'L5_NGC', 'L5_SBC', 'L5_STPC', 
    'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC', 'L6_BP', 'L6_BPC', 'L6_BTC', 'L6_ChC', 
    'L6_DBC', 'L6_IPC', 'L6_LBC', 'L6_MC', 'L6_NBC', 'L6_NGC', 'L6_SBC', 
    'L6_TPC_L1', 'L6_TPC_L4', 'L6_UTPC'
    ]

m_type_size = []
for m in m_type:
    m_type_size.append(connectivity[m][m]['cMat'].shape[0])
#############################################################################################
AA = pd.DataFrame(pd.concat([a_, pd.DataFrame(np.repeat(m_type, m_type_size))], axis = 1))
AA.columns = ['Dim0', 'Dim1', 'Dim2', 'Dim3', 'Dim4', 'Dim5', 'Dim6', 'n'] 
AA['l_type'] = AA['n'].str.split('_').str[0]
#############################################################################################
L_list = ['L1', 'L23', 'L4', 'L5', 'L6']
for l in L_list:
    ignore_zeroes = AA.loc[AA['l_type'] == l].drop(['n', 'l_type'], 
        axis = 1).mask(AA.loc[AA['l_type'] == l].drop(['n', 'l_type'], axis = 1) == 0)
    iz = AA.loc[AA['l_type'] == l].drop(['n', 'l_type'], axis = 1).assign(
        Start = ignore_zeroes.apply(pd.Series.first_valid_index, 1),
        Finish = ignore_zeroes.apply(pd.Series.last_valid_index, 1))

    g = pd.DataFrame([[item.replace('Dim', '') for item in lst] for lst in np.array(iz[['Finish']])], 
        columns = ['Dim+'])
    g['Dim+'] = g['Dim+'].astype(int)

    numerator = (g['Dim+'].sum(axis = 0)).astype(float)
    denominator = (g['Dim+'].value_counts().sum(axis = 0)).astype(float)

    average = numerator/denominator
    print(average)
