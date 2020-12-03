'''
Beginning with 'Simplicial Counts', the simplical values for the Bio-M model
can be obtained either by running the flagser program or by running the 
'original_bio_m.py' script and then running the 'general_biol_reconstruction.py'
script, whereas the simplicial values for General-Biol can be obtained by 
running the scripts in the general_biol folder.

The values for the excitatory and inhibitory simplicial counts can be obtained
by using the flagser package.

run using: 'python3 Figure2.py'
'''
import matplotlib.pyplot as plt
#####################################################################
#                         Simplicial Counts                         #
#####################################################################
x = [0, 1, 2, 3, 4, 5, 6] # dimensions
plt.plot(x, [31346, 7822274, 76361228, 64064185, 7274386, 156404, 896], label = "Bio-M")
plt.plot(x, [31346, 7648079, 34103325, 5158826, 43523, 30, 0], label = "General_Bio")
plt.xlabel(xlabel = "Dimension")
plt.ylabel(ylabel = "Number of Simplices")
plt.title(label = "Directed Simplex count for each model")
plt.legend()
plt.show()
#####################################################################
#                       Inhibitory and excitatory                   #
#####################################################################
plt.plot(x, [0, 0, 0, 0, 0, 0, 0], label = "Excitatory_MC6_L1")
plt.plot(x, [5870, 368777, 1080750, 334125, 15742, 122, 0], label = "Excitatory_MC6_L23")
plt.plot(x, [3089, 180398, 616152, 238916, 14539, 146, 0], label = "Excitatory_MC6_L4")
plt.plot(x, [5039, 616784, 4445800, 2908060, 216084, 2078, 1], label = "Excitatory_MC6_L5")
plt.plot(x, [6643, 677911, 3357760, 1577303, 101384, 1053, 3], label = "Excitatory_MC6_L6")
plt.title(label = "Average Connections MC6")
plt.legend()
plt.show()

plt.plot(x, [339, 848, 92, 1, 0, 0, 0], label = "Inhibitory_MC6_L1")
plt.plot(x, [1648, 21624, 19393, 2070, 26, 0, 0], label = "Inhibitory_MC6_L23")
plt.plot(x, [1572, 26104, 17218, 742, 2, 0, 0], label = "Inhibitory_MC6_L4")
plt.plot(x, [1067, 6597, 1715, 33, 0, 0, 0], label = "Inhibitory_MC6_L5")
plt.plot(x, [9261, 1147117, 7641907, 5769164, 735968, 20030, 184], label = "Inhibitory_MC6_L6")
plt.title(label = "Average Connections MC6")
plt.legend()
plt.show()
#####################################################################
