import matplotlib.pyplot as plt

# Values taken from running either flagser or python on the different filetypes
x = [0, 1, 2, 3, 4, 5, 6] # dimensions
new_erdos = [31000, 7685000.2, 15068000.4, 227500.2, 24000, 0, 0] # need to redo this with reconstruction files
MC0_I0 = [31346, 7648079, 56864241, 28696049, 1673099, 16288, 29] # control file, first rat, instance 0, cloud model
MC6 = [31346, 7822274, 76361228, 64064185, 7274386, 156404, 896] # original reconstruction of neocortical model, average of rats was taken (Bio-M)
MC0_combined = [34525, 5414314, 39712945, 26976971, 2410339, 36825, 47] # possibly controls of MC0 combined
General_Bio = [31346, 7648079, 34103325, 5158826, 43523, 30, 0] # reconstruction of Bio-M computed through python

plt.plot(x, new_erdos, label = "Erdos-Renyi")
plt.plot(x, MC6, label = "Bio-M")
plt.plot(x, MC0_I0, label = "Control Cloud")
plt.plot(x, General_Bio, label = "General_Bio")
plt.xlabel(xlabel = "Dimension")
plt.ylabel(ylabel = "Number of Simplices")
plt.title(label = "Directed Simplex count for each model")
plt.legend()
plt.show()
plt.savefig('Bio_M.png')

#####################################################################
#                       Erdos Renyi generated                       #
#####################################################################
# Not sure how these values were obtained
# Number of nonzero elements in the matrix : 7750
# In dimension 0 we have 310 simplices 
# In dimension 1 we have 7750 simplices 
# In dimension 2 we have 15312 simplices 
# In dimension 3 we have 2208 simplices 
# In dimension 4 we have 0 simplices 
# Euler characteristic : 5664
# 
# 
# Number of nonzero elements in the matrix : 7692
# In dimension 0 we have 310 simplices 
# In dimension 1 we have 7692 simplices 
# In dimension 2 we have 14994 simplices 
# In dimension 3 we have 2400 simplices 
# In dimension 4 we have 0 simplices 
# Euler characteristic : 5212


# In dimension 0 we have 310 simplices 
# In dimension 1 we have 7760 simplices 
# In dimension 2 we have 15786 simplices 
# In dimension 3 we have 2256 simplices 
# In dimension 4 we have 120 simplices 
# In dimension 5 we have 0 simplices 
# Euler characteristic : 6200

# In dimension 0 we have 310 simplices 
# In dimension 1 we have 7528 simplices 
# In dimension 2 we have 14118 simplices 
# In dimension 3 we have 1968 simplices 
# In dimension 4 we have 0 simplices 
# Euler characteristic : 4932


# In dimension 0 we have 310 simplices 
# In dimension 1 we have 7696 simplices 
# In dimension 2 we have 15132 simplices 
# In dimension 3 we have 2544 simplices 
# In dimension 4 we have 0 simplices 
# Euler characteristic : 5202

#####################################################################
#                       Inhibitory and excitatory                   #
#                           MC0 and MC6                             #
#####################################################################
# potentially through flagser after converting h5 files to individual layers
# checking through BBP website which m_types are excitatory and inhibitory
# must be a better way to use plt

# Excitatory Graphs - MC6
x1 = [0, 1, 2, 3, 4, 5, 6]
Excitatory_MC6_L1 = [0, 0, 0, 0, 0, 0, 0]
Excitatory_MC6_L23 = [5870, 368777, 1080750, 334125, 15742, 122, 0]
Excitatory_MC6_L4 = [3089, 180398, 616152, 238916, 14539, 146, 0]
Excitatory_MC6_L5 = [5039, 616784, 4445800, 2908060, 216084, 2078, 1]
Excitatory_MC6_L6 = [6643, 677911, 3357760, 1577303, 101384, 1053, 3]

plt.plot(x1, Excitatory_MC6_L1, label = "Excitatory_MC6_L1")
plt.plot(x1, Excitatory_MC6_L23, label = "Excitatory_MC6_L23")
plt.plot(x1, Excitatory_MC6_L4, label = "Excitatory_MC6_L4")
plt.plot(x1, Excitatory_MC6_L5, label = "Excitatory_MC6_L5")
plt.plot(x1, Excitatory_MC6_L6, label = "Excitatory_MC6_L6")
plt.title(label = "Average Connections MC6")
plt.legend()
plt.show()
plt.savefig('Excitatory_graph_MC6.png')

# Inhibitory Graphs - MC6
Inhibitory_MC6_L1 = [339, 848, 92, 1, 0, 0, 0]
Inhibitory_MC6_L23 = [1648, 21624, 19393, 2070, 26, 0, 0]
Inhibitory_MC6_L4 = [1572, 26104, 17218, 742, 2, 0, 0]
Inhibitory_MC6_L5 = [1067, 6597, 1715, 33, 0, 0, 0]
Inhibitory_MC6_L6 = [9261, 1147117, 7641907, 5769164, 735968, 20030, 184]

plt.plot(x1, Inhibitory_MC6_L1, label = "Inhibitory_MC6_L1")
plt.plot(x1, Inhibitory_MC6_L23, label = "Inhibitory_MC6_L23")
plt.plot(x1, Inhibitory_MC6_L4, label = "Inhibitory_MC6_L4")
plt.plot(x1, Inhibitory_MC6_L5, label = "Inhibitory_MC6_L5")
plt.plot(x1, Inhibitory_MC6_L6, label = "Inhibitory_MC6_L6")
plt.title(label = "Average Connections MC6")
plt.legend()
plt.show()
plt.savefig('Inhibitory_graph_MC6.png')

#####################################################################
