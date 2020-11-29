import pandas as pd
import h5py
import numpy as np
from pyflagsercontain import flagser_count

# file reading and variable setting. Average files (0-5)
# I believe the graphic that comes from this is a summation of all 6 instantiations (noted below)
mc0 = h5py.File('../../pathway_average_files/cons_locs_pathways_mc6_Column.h5', 'r')
connectivity = mc0.get('connectivity')

####################################################################################
#######################           Full Matrix Code           #######################
####################################################################################
# reproduces whole connectivity matrix
whole_matrix = np.vstack([np.hstack([connectivity[a][b]['cMat'] for b in connectivity]) for a in connectivity])

simplex_counts = flagser_count(whole_matrix)
assert len(simplex_counts) == 31346

full_size = [i for i in range(31346)]
a_Full = []
for i in full_size:
    a_Full.append(simplex_counts[i])

a_ = pd.DataFrame(a_Full).fillna(0)
L_ = a_[3]

# better way to write this in figure 2 reconstruction files
L1_DACa = L_[0:59]
L1_DLACa = L_[60:83]
L1_HACa = L_[84:174]
L1_NGCa_DA = L_[175:245]
L1_NGCa_SA = L_[246:297]
L1_SLACa = L_[298:340]
L23_BPa = L_[341:368]
L23_BTCa = L_[369:473]
L23_ChCa = L_[474:534]
L23_DBCa = L_[535:709]
L23_LBCa = L_[710:1162]
L23_MCa = L_[1163:1494]
L23_NBCa = L_[1495:1760]
L23_NGCa = L_[1761:1816]
L23_PCa = L_[1817:7692]
L23_SBCa = L_[7693:7859]
L4_BPa = L_[7860:7867]
L4_BTCa = L_[7868:7887]
L4_ChCa = L_[7888:7895]
L4_DBCa = L_[7896:7932]
L4_LBCa = L_[7933:8053]
L4_MCa = L_[8054:8170]
L4_NBCa = L_[8171:8266]
L4_NGCa = L_[8267:8272]
L4_PCa = L_[8273:10953]
L4_SBCa = L_[10954:11015]
L4_SPa = L_[11016:12109]
L4_SSa = L_[12110:12519]
L5_BPa = L_[12520:12553]
L5_BTCa = L_[12554:12628]
L5_ChCa = L_[12629:12647]
L5_DBCa = L_[12648:12743]
L5_LBCa = L_[12744:12951]
L5_MCa = L_[12952:13349]
L5_NBCa = L_[13350:13551]
L5_NGCa = L_[13552:13559]
L5_SBCa = L_[13560:13584]
L5_STPCa = L_[13585:13885]
L5_TTPC1a = L_[13886:16287]
L5_TTPC2a = L_[16288:18284]
L5_UTPCa = L_[18285:18627]
L6_BPa = L_[18628:18697]
L6_BPCa = L_[18698:21876]
L6_BTCa = L_[21877:21931]
L6_ChCa = L_[21932:21947]
L6_DBCa = L_[21948:21978]
L6_IPCa = L_[21979:25448]
L6_LBCa = L_[25449:25914]
L6_MCa = L_[25915:26247]
L6_NBCa = L_[26248:26445]
L6_NGCa = L_[26446:26462]
L6_SBCa = L_[26463:26528]
L6_TPCa_L1 = L_[26529:28168]
L6_TPCa_L4 = L_[28169:29609]
L6_UTPCa = L_[29610:31346]

# Note of which m_types are inhibitory and excitatory
inhib = [L1_DACa,L1_DLACa,L1_HACa,L1_NGCa_DA,L1_NGCa_SA,L1_SLACa,L23_BTCa,L23_DBCa,L23_LBCa,L23_MCa,L23_NBCa,L23_NGCa,L23_SBCa,L4_BTCa,L4_DBCa,L4_LBCa,L4_MCa,L4_NBCa,L4_NGCa,L4_SBCa,L5_BTCa,L5_DBCa,L5_LBCa,L5_MCa,L5_NBCa,L5_NGCa,L6_BTCa,L6_MCa,L6_LBCa]

excitatory = [L23_PCa,L4_PCa,L4_SPa,L4_SSa,L5_STPCa,L5_TTPC1a,L5_TTPC2a,L5_UTPCa,L6_BPCa,L6_IPCa,L6_TPCa_L1,L6_TPCa_L4,L6_UTPCa]

inhibitory_data = pd.concat(inhib)
excitatory_data = pd.concat(excitatory)


# Inhibitory
L_counts = pd.DataFrame(inhibitory_data.value_counts())#.index.tolist())
L_counts['simplex'] = L_counts.index
L_counts.columns = ['Neuron', 'simplex']
# L_counts.to_csv("qwer.csv")
L_counts_1 = pd.DataFrame(L_counts.groupby(pd.cut(L_counts.simplex, bins = np.linspace(0, 40000, 20)))['Neuron'].sum().fillna(0))
L_counts_1.columns = ['neurons']
L_counts_1.to_csv("i6.csv")

# Excitatory
# L_counts = pd.DataFrame(excitatory_data.value_counts())
# L_counts['simplex'] = L_counts.index
# L_counts.columns = ['neuron','simplex']
# L_counts = pd.DataFrame(L_counts.groupby(pd.cut(L_counts.simplex, bins = np.linspace(0, 40000, 21)))['neuron'].sum().fillna(0))
# L_counts.columns = ['neurons']
# L_counts.to_csv("e6.csv")


####################################################################################
# counts for number of neurons belonging to each of 3 dimensional simplices.
# summed across 6 instantiations of the average files
