import h5py
import pandas as pd
import matplotlib.pyplot as plt

#############################################################################
#                               Graphs                                      #
#############################################################################
#############################################################################
#                             Graph A1                                      #
#############################################################################
# x = [0, 1, 2, 3, 4, 5, 6]
# mc6_full = [31346, 7822274, 76361228, 64064185, 7274386, 156404, 896]
# mc6_excite = [27026, 6827328, 66858753, 56366598, 6490582, 142827, 833]

# plt.plot(x, mc6_full, label = "mc6 Full")
# plt.plot(x, mc6_excite, label = "excitatory subgraph")
# plt.title(label = "Excitatory Sub-Graph (A1)")
# plt.xlabel(xlabel = "Dimension")
# plt.ylabel(ylabel = "Simplices")
# plt.legend()
# # plt.show()
# # plt.savefig('Fig3_A1.png')

# # To get the data, I ran the flagser-count program on the MC6 file in average.
# # To get the excitatory neurons, I checked the Downloads page on BBP website
# # and found which neurons were excitatory and added them to a flag file and
# # performed the flagser-count to find how many simplices belonged to each 
# # dimension.

# #############################################################################
# #                                Graph A2                                   #
# #############################################################################
# x1 = [0, 1, 2, 3, 4, 5]
# mc6_inhibitory = [7502, 313856, 749279, 210627, 9327, 75]

# plt.plot(x1, mc6_inhibitory, label = "Inhibitory subgraph")
# plt.title(label = "Inhibitory Sub-Graph MC6 (A2)")
# plt.xlabel(xlabel = "Dimension")
# plt.ylabel(ylabel = "Simplices")
# plt.legend()
# # plt.show()
# # plt.savefig('Fig3_A2.png')

# # This is similar to above, I just took all the inhibitory neurons, collected
# # them together as a group and calculated the number of simplices in each dimension
# # Not correct (mismatch here, need to find out why) -> with the only difference to the original papers graph being that my count here
# # is approximately an order of magnitude larger. 

# #############################################################################
# #                               Graph A3                                    #
# #############################################################################
# # Excitatory Sub Layers
# x2 = [0, 1, 2, 3, 4, 5, 6]
# excitatory_L1 = [0, 0, 0, 0, 0, 0, 0]
# excitatory_L23 =  [5870, 368777, 1080750, 334125, 15742, 122, 0]
# excitatory_L4 = [4186, 304355, 1198048, 500731, 31864, 344, 0]
# excitatory_L5 = [5039, 616784, 4445800, 2908060, 216084, 2078, 1]
# excitatary_L6 = [11931, 2026213, 17743995, 16586000, 2463924, 74345, 546]

# plt.plot(x2, excitatory_L1, x2, excitatory_L23, x2, excitatory_L4,x2, excitatory_L5,x2, excitatary_L6)
# plt.title(label = "Excitatory Sub-graph MC6 (Average - A3)")
# labels = ["L1", "L23", "L4", "L5", "L6"]
# plt.xlabel(xlabel = "Dimension")
# plt.ylabel(ylabel = "Simplices")
# plt.legend(labels = labels)
# plt.show()
# plt.savefig('Fig3_A3.png')

# For this graphic, I separated the excitatory neurons of MC6 average into the 
# relevant layers and calculated the number of simplices in each dimension. Again, we used the flagser package


#############################################################################
#                               Graph B                                     #
#############################################################################
######################## Figure 3B Excitatory ########################
E_0 = pd.DataFrame(pd.read_csv("figure3b/e0.csv"))['neurons']
E_1 = pd.DataFrame(pd.read_csv("figure3b/e1.csv"))['neurons']
E_2 = pd.DataFrame(pd.read_csv("figure3b/e2.csv"))['neurons']
E_3 = pd.DataFrame(pd.read_csv("figure3b/e3.csv"))['neurons']
E_4 = pd.DataFrame(pd.read_csv("figure3b/e4.csv"))['neurons']
E_5 = pd.DataFrame(pd.read_csv("figure3b/e5.csv"))['neurons']
E_6 = pd.DataFrame(pd.read_csv("figure3b/e6.csv"))['neurons']

e_concat = [E_0, E_1, E_2, E_3, E_4, E_5, E_6]
e_list = pd.DataFrame(pd.concat(e_concat, axis = 1))

e_list1 = pd.DataFrame(e_list.sum(axis = 1))
e_list1.columns = ['neurons']

e_list1['index'] = e_list1.index
ay = plt.bar(x = e_list1['index'], height = e_list1['neurons'], width = 1, label = "Excitatory")
# plt.show()
 
print(e_list1)
# e_list = pd.DataFrame(pd.concat(e_concat))
# e_list.columns = ['neuron', 'Dim3']
# e_list1 = e_list.sort_values(by = 'Dim3', ascending = False)
# e_list1['list'] = range(len(e_list1))
# ay = e_list1.plot(kind = 'line', x = 'Dim3', y = 'list')
# ay.set_xlabel("3D Simplices/ Neuron")
# ay.set_ylabel("Neuron")
# ax.grid('on')

######################## Figure 3B Inhibitory ########################
M_0 = pd.DataFrame(pd.read_csv("figure3b/i0.csv"))['neurons']
M_1 = pd.DataFrame(pd.read_csv("figure3b/i1.csv"))['neurons']
M_2 = pd.DataFrame(pd.read_csv("figure3b/i2.csv"))['neurons']
M_3 = pd.DataFrame(pd.read_csv("figure3b/i3.csv"))['neurons']
M_4 = pd.DataFrame(pd.read_csv("figure3b/i4.csv"))['neurons']
M_5 = pd.DataFrame(pd.read_csv("figure3b/i5.csv"))['neurons']
M_6 = pd.DataFrame(pd.read_csv("figure3b/i6.csv"))['neurons']

concat_list = [M_0, M_1, M_2, M_3, M_4, M_5, M_6]
check_list = pd.DataFrame(pd.concat(concat_list, axis = 1))

check_list1 = pd.DataFrame(check_list.sum(axis = 1))
check_list1.columns = ['neurons']

check_list1['index'] = check_list1.index

ax = plt.bar(x = check_list1['index'], height = check_list1['neurons'], width = 1, label = "Inhibitory")
# plt.show()

# check_list.columns = ['neuron', 'Dim3']
# check_list1 = check_list.sort_values(by = 'Dim3', ascending = False)
# check_list1['list'] = range(len(check_list1))

# ax = check_list1.plot(kind = 'line', x = 'Dim3', y = 'list')
# ax.set_xlabel("3D simplices/neuron")
# ax.set_ylabel("Neuron")
# plt.show()
plt.legend()
plt.show()

# values for the above (3b) were extracted by the "synaptic_code_fig_3b" script
#############################################################################
#                               Graph D                                     #
#############################################################################
# a = pd.read_csv("figure3d/figure3D.csv")
# a.columns = ['x', 'L1', 'L23', 'L4', 'L5', 'L6']

# plt.plot(a.x, a.L1, a.x, a.L23, a.x, a.L4, a.x, a.L5, a.x, a.L6)
# labels = ["L1", "L23", "L4", "L5", "L6"]
# plt.xlabel(xlabel = "Dimension")
# plt.ylabel(ylabel = "Simplices/Neuron")
# plt.title(label = "Number of simplices per Neuron per Layer")
# plt.legend(labels = labels)
# plt.show()
# plt.savefig('Fig3_D.png')

# This represents the number of simplices per neuron for the file 
# average_pathways_6 (Bio-M). So this is the average over all of the tests.
