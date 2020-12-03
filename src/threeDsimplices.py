'''
We read in the original h5 files, specifically the connectivity matrices. This is done for
all instantiations of the Neocortical Model (Bio0-5). From here, we concatenate all 
the submatrices into the large 31,346 x 31,346 matrix. 

Then the number of simplices to which each neuron belongs to is the computed through
simplex counts, and listed in its entirety in a_Full. The column of interest, however,
is the column including the number of 3D simplices in which each neuron belongs to. This
is labelled as L_.

Further to finding the number of 3D simplices belonging to each neuron, we want to see
this for both inhibitory and excitatory neurons. This we successfully split out in 
V_. 

To finally get the data in the form required, we take counts of the number of simplices
each neuron belongs to and list these as such under L_inhib and L_excit. Once we have these
counts, we can then group together these counts into bins.
'''
####################################################################################
import pandas as pd
import h5py
import numpy as np
from pyflagsercontain import flagser_count
####################################################################################
numbers = np.arange(0, 6, 1)
for z in numbers:
    mc0 = h5py.File('../data/cons_locs_pathways_mc' + str(z) + '_Column.h5', 'r')
    connectivity = mc0.get('connectivity')
####################################################################################
    whole_matrix = np.vstack([np.hstack([connectivity[a][b]['cMat'] for b in connectivity]) for a in connectivity])

    simplex_counts = flagser_count(whole_matrix)
    assert len(simplex_counts) == 31346

    a_Full = []
    for i in [i for i in range(31346)]:
        a_Full.append(simplex_counts[i])

    L_ = pd.DataFrame(a_Full).fillna(0)[3]
####################################################################################
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
###################################################################################
    inhib = [
    'L1_DAC','L1_DLAC','L1_HAC','L1_NGC-DA','L1_NGC-SA','L1_SLAC','L23_BTC','L23_DBC',
    'L23_LBC','L23_MC','L23_NBC','L23_NGC','L23_SBC','L4_BTC','L4_DBC','L4_LBC',
    'L4_MC','L4_NBC','L4_NGC','L4_SBC','L5_BTC','L5_DBC','L5_LBC','L5_MC','L5_NBC',
    'L5_NGC','L6_BTC','L6_MC','L6_LBC'
        ]
###################################################################################
    excit = [
    'L23_PC','L4_PC','L4_SP','L4_SS','L5_STPC','L5_TTPC1','L5_TTPC2','L5_UTPC','L6_BPC',
    'L6_IPC','L6_TPC_L1','L6_TPC_L4','L6_UTPC'
        ]
###################################################################################
    m_type_size = []
    for f in m_type:
        size = connectivity[f][f]['cMat'].shape[0]
        m_type_size.append(size)
###################################################################################
    V_ = pd.DataFrame(pd.concat([L_, pd.DataFrame(np.repeat(m_type, m_type_size))], axis = 1))
    V_.columns = ['3D_simplices', 'm_type']
###################################################################################
    inhibitory = V_[V_['m_type'].isin(inhib)]['3D_simplices']
    excitatory = V_[~V_['m_type'].isin(inhib)]['3D_simplices']
###################################################################################
    L_inhib = pd.DataFrame(inhibitory.value_counts())
    L_inhib['simplex'] = L_inhib.index
    L_inhib.columns = ['neuron', 'simplex']
    L_inhib = pd.DataFrame(L_inhib.groupby(pd.cut(L_inhib.simplex, bins = np.linspace(0, 40000, 20)))['neuron'].sum().fillna(0))
    L_inhib.columns = ['neurons']
    L_inhib.to_csv("../output/i" + str(z) + ".csv")
###################################################################################
    L_excit = pd.DataFrame(excitatory.value_counts())
    L_excit['simplex'] = L_excit.index
    L_excit.columns = ['neuron','simplex']
    L_excit = pd.DataFrame(L_excit.groupby(pd.cut(L_excit.simplex, bins = np.linspace(0, 40000, 21)))['neuron'].sum().fillna(0))
    L_excit.columns = ['neurons']
    L_excit.to_csv("../output/e" + str(z) + ".csv")
###################################################################################
