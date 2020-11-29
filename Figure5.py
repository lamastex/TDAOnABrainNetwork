import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# 42 variants
# six sets
# seven varying instantiations of each
# so this is the individual pathways
#########################################################################################
#                                      Graph B                                          #
#########################################################################################
Euler_1 = [12252037,11940527,12006995,12321329,11686485,11813241,10625627]
Euler_2 = [9379929,9229543,9639831,9711248,9763679,7912652,8345991]
Euler_3 = [9134068,8079536,7959119,8036974,8725798,7186736,6546731]
Euler_4 = [13160606,12458620,12538837,12069465,12075500,11463811,11217052]
Euler_5 = [12179553,12015541,11993859,11804031,12329492,11941146,11615349] 

Bio_1 = [343,434,346,338,441,417,447]
Bio_2 = [1107,1192,1139,1028,1030,1392,1364]
Bio_3 = [1750,2020,1793,1821,1607,1857,2320]
Bio_4 = [1143,1366,1211,1300,1262,1360,1432]
Bio_5 = [235,288,296,253,267,331,387]

plt.scatter(Euler_1, Bio_1, label = "Bio_1s", color = "red")
plt.scatter(Euler_2, Bio_2, label = "Bio_2s", color = "green")
plt.scatter(Euler_3, Bio_3, label = "Bio_3s", color = "purple")
plt.scatter(Euler_4, Bio_4, label = "Bio_4s", color = "blue")
plt.scatter(Euler_5, Bio_5, label = "Bio_5s", color = "cyan")
plt.title(label = "Reconstruction of Revised Control Dataset")
plt.xlabel(xlabel = "Euler Characterisitc")
plt.ylabel(ylabel = "Betti Number B6")
plt.legend()
plt.show()

#########################################################################################
# These are based on the revised cloud model version, so give a different result to the original figure (fig5b)


