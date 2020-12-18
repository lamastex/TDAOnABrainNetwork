'''
This is the Peters rule implementation of the neural connectivity. The idea is that first we take
the total lengths of the axons and dendrites within each m_type along with the number of afferent and
efferent connection in each m_type and divide the two so that we can an average length of both axons and
dendrites for each m_type. Once we have these, we can then add the two relevant lengths so that we create
an average length synaptic connection between each neuron. We then compute the distances between each 
neuron. With this, we can then subtract the synaptic length from the neuron distances to check if there
could be a potential connection. Having values less than a certain number gave the potential for 
having spatial apositions, i.e. instances where an axon of an outgoing synapse was close enough to the
incoming dendrite of another neuron to form a connection. As mentioned below, some values were used
to experiment with, whereby we settled with 750 in order to garner approximately the same value as 
stated in another paper.

First up, we need to collect the sizes of each matrix which means we run everything 
in the first section up until: np.save('m_type', m_size). Then we can simply load 
this information using the following line: m_typesize = np.load('m_type.npy').

Once completed, we comment this out. The next step is to then compute the average length
of each synapse for each m_type combination. The synapse consists of axons and dendrites,
so details of their absolute length for each m_type can be found on the BPP website, but 
the details have been added below. We then divide these lengths by the counts of afferent 
and efferent connections for each m_type and divide the total length by the count of these.
This is then saved as: # np.save('outters', overlap).

We then load this data in using: overlap = np.load('outters.npy'). The rest above is commented
out and we can then run the for loop afterwards to create the reconstructed matrices. There
has been some experimenting regarding the distance in order to match the number of potential connections
using this type of connectivity, whereupon we have landed at a distance of around 750 or less. 
This give ~600 million potential connections.
'''
import numpy as np
import h5py
import pandas as pd
import itertools
from scipy.spatial import distance
import scipy
import random
import os

# Read in H5 file
mc0_file = h5py.File('../../average/cons_locs_pathways_mc0_Column.h5', 'r')
###########################################################################################
#                                 Neuron Type Locations                                   #
###########################################################################################
populations = mc0_file.get('populations')
connection = mc0_file.get('connectivity')
###########################################################################################
# M_types
m_type = [
        'L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 'L23_BP', 'L23_BTC', 
        'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC', 'L23_PC', 'L23_SBC', 
        'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC', u'L4_MC', 'L4_NBC', 'L4_NGC', 'L4_PC', 
        'L4_SBC', 'L4_SP', 'L4_SS', 'L5_BP', 'L5_BTC', 'L5_ChC', 'L5_DBC', 'L5_LBC', 'L5_MC', 
        'L5_NBC', 'L5_NGC', 'L5_SBC', 'L5_STPC', 'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC', 'L6_BP', 
        'L6_BPC', 'L6_BTC', 'L6_ChC', 'L6_DBC', 'L6_IPC', 'L6_LBC', 'L6_MC', 'L6_NBC', 'L6_NGC', 
        'L6_SBC', 'L6_TPC_L1', 'L6_TPC_L4', 'L6_UTPC'
        ]
###########################################################################################
# matrix_size = []
# for M_a in m_type:
#     ma_size = np.array(pd.DataFrame(connection[M_a][M_a]['cMat']).shape[0])
#     matrix_size.append(ma_size)

# m_size = np.array(matrix_size)
# np.save('m_type', m_size)    
m_typesize = np.load('m_type.npy')
###########################################################################################
# maff = []
# meff = []
# for M_a in m_type:
#     for M_b in m_type:
#         aff = pd.DataFrame(populations[M_a]['nCellAff'])
#         eff = pd.DataFrame(populations[M_b]['nCellEff'])
#         maff.append(aff.sum(axis = 0))
#         meff.append(eff.sum(axis = 0))
###########################################################################################
#                            Axon and Dendrite Total Lengths                              #
# ###########################################################################################
# ax_list_total = np.repeat([
#               332020., 246190., 397900., 340800., 151340., 101790., 134220., 126340., 1082900., 
#               1035100., 3832000., 6671400., 4482100., 381480., 51315000., 2931500., 267620., 
#               1835500., 575870., 666360., 2370500., 2324500., 2230200., 316050., 41523000., 
#               1619300., 3559900., 6749700., 173440., 3785900., 560840., 368770., 5849600., 
#               6248400., 2381800., 87164., 1629000., 5646900., 30428000., 43306000., 8651100., 
#               168130., 16545000., 2085300., 322030., 390450., 16672000., 5096200., 4362000., 
#               4267200., 177210., 1093500., 16353000., 12505000., 15390000.
#               ], 55)

# den_list_total = np.repeat([
#               65834., 35862., 56755., 33329., 34513., 37237., 45289., 667590., 148370., 256080., 
#               845150., 1504200., 1222200., 180580., 13061000., 591530., 32157., 283880., 84597., 
#               123220., 367700., 572000., 355050., 144340., 9763100., 310830., 944030., 14566., 
#               172050., 1100800., 92836., 299520., 1054500., 1158300., 426210., 50163., 319360., 
#               4122500., 12578000., 21216000., 3215500., 166450., 10648000., 611610., 51222., 
#               220350., 7981700., 937160., 1118600., 856760., 103140., 265590., 7776600., 6342100., 
#               7439000.
#               ], 55)
###########################################################################################
#                                 Average Axon + Dendrite                                 #
###########################################################################################
# table_contents = [
#               pd.DataFrame(maff), 
#               pd.DataFrame(meff), 
#               pd.DataFrame(ax_list_total), 
#               pd.DataFrame(den_list_total)
#               ]

# tab = pd.DataFrame(pd.concat(table_contents, axis = 1))

# tab.columns = [
#               'Aff', 
#               'Eff', 
#               'Aff_Total', 
#               'Eff_Total'
#               ]

# tab['aff_average'] = tab['Aff_Total']/tab['Aff']
# tab['eff_average'] = tab['Eff_Total']/tab['Eff']
# tab['synapse_length'] = tab['aff_average'] + tab['eff_average']

# overlap = np.array(tab['synapse_length'])
# np.save('outters', overlap)
###########################################################################################
#                       Matrix for each submatrix of connections                          #
###########################################################################################
overlap = np.load('outters.npy')
chunks = pd.DataFrame(np.array([overlap[x:x+55] for x in range(0, len(overlap), 55)]))
###########################################################################################
#                                        All cases                                        #
###########################################################################################
for i, iz in zip(m_type, chunks):
    for j, jz in zip(m_type, chunks):
        # locations and distances between neurons in each matrix
        L_a = pd.DataFrame(np.matrix(populations[i]['locations']), columns = ['x', 'y', 'z'])
        L_b = pd.DataFrame(np.matrix(populations[j]['locations']), columns = ['x', 'y', 'z']) 
        D_ = scipy.spatial.distance.cdist(L_a, L_b, 'euclidean')
        
        # lengths is an array by array filling of the values from chunks
        lengths = np.full((len(L_a), len(L_b)), chunks[jz][iz])

        # poss_connect takes the distances between actual neurons and the potential synapse length
        poss_connect = np.array(np.array(D_) - np.array(lengths))

        # converts poss_connect to boolean matrix and then remove any 1s in diagonal
        overconnected_peter = ((np.array(poss_connect) < 750).astype(int))
        np.fill_diagonal(overconnected_peter, 0) 
        
        # connect takes the number of 1s in each matrix
        connect = (np.array(overconnected_peter) == 1).sum()

        # rcc reads in real connection matrix
        rcc = np.array(connection[i][j]['cMat'])
        
        # real connect count takes the real number of connections in each matrix
        real_connect_count = (rcc.sum()).sum()

        # idx is used in order to take from the list and replace in the first if statement        
        idx = np.flatnonzero(overconnected_peter)
        
        # num_replaced is used in order to count how many instances of 1 need to be set to 0
        num_replaced = connect - real_connect_count

        if num_replaced > 0:
                peter = overconnected_peter
                np.put(peter,np.random.choice(idx,size=num_replaced,replace=False),0)
        else:
                peter = rcc

        a = pd.DataFrame(peter)
        a.to_csv("reconstruction/" + str(i) + str(j) + ".csv", header = False, index = False)