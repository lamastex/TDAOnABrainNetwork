import numpy as np
import h5py
import pandas as pd
import collections
from mpl_toolkits.mplot3d import Axes3D
import itertools
from scipy.spatial import distance
import scipy
import random
import os

# Read in H5 file
mc0_file = h5py.File('pathway_average_files/cons_locs_pathways_mc0_Column.h5', 'r')
###########################################################################################
#                                 Neuron Type Locations                                   #
###########################################################################################
populations = mc0_file.get('populations')
connection = mc0_file.get('connectivity')

# # M_types
m_type = ['L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 'L23_BP', 'L23_BTC', 'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC', 'L23_PC', 'L23_SBC', 'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC', u'L4_MC', 'L4_NBC', 'L4_NGC', 'L4_PC', 'L4_SBC', 'L4_SP', 'L4_SS', 'L5_BP', 'L5_BTC', 'L5_ChC', 'L5_DBC', 'L5_LBC', 'L5_MC', 'L5_NBC', 'L5_NGC', 'L5_SBC', 'L5_STPC', 'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC', 'L6_BP', 'L6_BPC', 'L6_BTC', 'L6_ChC', 'L6_DBC', 'L6_IPC', 'L6_LBC', 'L6_MC', 'L6_NBC', 'L6_NGC', 'L6_SBC', 'L6_TPC_L1', 'L6_TPC_L4', 'L6_UTPC']

# matrix_size = []
# for M_a in m_type:
#     ma_size = np.array(pd.DataFrame(connection[M_a][M_a]['cMat']).shape[0])
#     matrix_size.append(ma_size)

    

# m_size = np.array(matrix_size)
# np.save('m_type', m_size)    
# m_typesize = np.load('m_type.npy')

# maff = []
# meff = []
# for M_a in m_type:
#     for M_b in m_type:
#         aff = pd.DataFrame(populations[M_a]['nCellAff'])
#         eff = pd.DataFrame(populations[M_b]['nCellEff'])
#         maff.append(aff.sum(axis = 0))
#         meff.append(eff.sum(axis = 0))
        
# # ###########################################################################################
# # #                            Axon and Dendrite Total Lengths                              #
# # ###########################################################################################
# ax_list_total = np.repeat([332020., 246190., 397900., 340800., 151340., 101790., 134220., 126340., 1082900., 1035100., 3832000., 6671400., 4482100., 381480., 51315000., 2931500., 267620., 1835500., 575870., 666360., 2370500., 2324500., 2230200., 316050., 41523000., 1619300., 3559900., 6749700., 173440., 3785900., 560840., 368770., 5849600., 6248400., 2381800., 87164., 1629000., 5646900., 30428000., 43306000., 8651100., 168130., 16545000., 2085300., 322030., 390450., 16672000., 5096200., 4362000., 4267200., 177210., 1093500., 16353000., 12505000., 15390000.], 55)

# den_list_total = np.repeat([65834., 35862., 56755., 33329., 34513., 37237., 45289., 667590., 148370., 256080., 845150., 1504200., 1222200., 180580., 13061000., 591530., 32157., 283880., 84597., 123220., 367700., 572000., 355050., 144340., 9763100., 310830., 944030., 14566., 172050., 1100800., 92836., 299520., 1054500., 1158300., 426210., 50163., 319360., 4122500., 12578000., 21216000., 3215500., 166450., 10648000., 611610., 51222., 220350., 7981700., 937160., 1118600., 856760., 103140., 265590., 7776600., 6342100., 7439000.], 55)

# # ###########################################################################################
# # #                                 Average Axon + Dendrite                                 #
# # ###########################################################################################
# table_contents = [pd.DataFrame(maff), pd.DataFrame(meff), pd.DataFrame(ax_list_total), pd.DataFrame(den_list_total)]
# tab = pd.DataFrame(pd.concat(table_contents, axis = 1))
# tab.columns = ['Aff', 'Eff', 'Aff_Total', 'Eff_Total']

# tab['aff_average'] = tab['Aff_Total']/tab['Aff']
# tab['eff_average'] = tab['Eff_Total']/tab['Eff']
# tab['synapse_length'] = tab['aff_average'] + tab['eff_average']

# overlap = np.array(tab['synapse_length'])
# np.save('outters', overlap)
# ###########################################################################################
# #                       Matrix for each submatrix of connections                          #
# ###########################################################################################
# # Size will need to alter to size of connection matrix
overlap = np.load('outters.npy')
chunks = pd.DataFrame(np.array([overlap[x:x+55] for x in range(0, len(overlap), 55)]))

# ###########################################################################################
# #                                        All cases                                        #
# ###########################################################################################
for i, iz in zip(m_type, chunks):
    for j, jz in zip(m_type, chunks):
        L_a = pd.DataFrame(np.matrix(populations[i]['locations']), columns = ['x', 'y', 'z'])
        L_b = pd.DataFrame(np.matrix(populations[j]['locations']), columns = ['x', 'y', 'z'])
        D_ = scipy.spatial.distance.cdist(L_a, L_b, 'euclidean')
        
        lengths = np.full((len(L_a), len(L_b)), chunks[jz][iz])
        
        # not correct I dont believe
#         poss_connect = np.array(np.array(D_) - np.array(lengths))
#         overconnected_peter = (np.array(poss_connect) < 1).astype(int)
#         np.fill_diagonal(overconnected_peter, 0)
#         overconnected_peter = np.array(overconnected_peter)
        
        # 
        poss_connect = np.array(np.array(D_) - np.array(lengths))
        upper = 1
        lower = -1
#         (np.array(poss_connect)).any(-1 < np.array(poss_connect) < 1).astype(int)
#         overconnected_peter = [i for i in np.array(poss_connect) if lower < i < upper]
#         overconnected_peter = ((np.array(poss_connect) < 1).any() and (np.array(poss_connect)> -1).any()).astype(int)
        overconnected_peter = ((np.array(poss_connect) < 1)&(np.array(poss_connect) > -1).astype(int))

        np.fill_diagonal(overconnected_peter, 0) 
        overconnected_peter = np.array(overconnected_peter)
        
        
        connect = (np.array(overconnected_peter) == 1).sum()
        print(connect)
        rcc = np.array(pd.DataFrame(connection[i][j]['cMat']))
        real_connect_count = (rcc.sum()).sum()
        print(real_connect_count)
        print('#############')
        

        idx = np.flatnonzero(overconnected_peter)
        num_replaced = connect - real_connect_count

#         print(num_replaced)
#         if num_replaced > 0:
#             np.put(overconnected_peter,np.random.choice(idx,size=num_replaced,replace=False),0)
#         else:
#             overconnected_peter = pd.DataFrame(connection['L1_DAC']['L1_DAC']['cMat'], columns = None, index = None)

        a = pd.DataFrame(overconnected_peter)
        a.to_csv("reconstruction/" + str(i) + str(j) + ".csv", header = False, index = False)

# 11m47 to compile

# real	15m6,458s
###########################################################################################
#                                        One case                                         #
###########################################################################################
# L_a = pd.DataFrame(np.matrix(populations['L1_DAC']['locations']), columns = ['x', 'y', 'z'])
# L_b = pd.DataFrame(np.matrix(populations['L1_DAC']['locations']), columns = ['x', 'y', 'z'])
# D_ = scipy.spatial.distance.cdist(L_a, L_b, 'euclidean')
        
# lengths = np.full((len(L_a), len(L_b)), 203)
        
# poss_connect = np.array(np.array(D_) - np.array(lengths))
# overconnected_peter = (np.array(poss_connect) < 1).astype(int)
# print(overconnected_peter)

# np.fill_diagonal(overconnected_peter, 0)
# print(overconnected_peter.sum())

# connect = (np.array(overconnected_peter) == 1).sum()
# print(connect)
# print(pd.DataFrame(connection['L1_DAC']['L1_DAC']['cMat'], columns = None, index = None))
# print(poss_connect.sum())

# L1_DAC to L1_DAC simplices, error from above
# [59, 657, 2118, 2486, 1228, 281, 24]

        
#####################################################################
# WIP. Not sure how to implement this model yet
#####################################################################