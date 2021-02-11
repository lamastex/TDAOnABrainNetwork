# working-manuscript-TopologicalDataAnalysisOnABrainNetwork
working-manuscript-TopologicalDataAnalysisOnABrainNetwork - for Kieran Barber

# Requirements
To run the code within, you will need the following:
* pyflagser: https://pypi.org/project/pyflagser/
* pyflagsercontain: `pip install pyflagsercontain`
* cons_loc_pathway files: https://bbp.epfl.ch/nmc-portal/downloads (average files)

# Empirical Geometric Configuration model (4 steps)
Four scripts to run. Before running any of these scripts, however, you will need to download the pathway files (average files).
Once downloaded, the first script to run is `h5_tocsv.py`. This will convert the h5 files into csv files and will put them in the `reconstruction` folder. 
Second file to run is `concatenate_csv.py`. This will return the array that was in h5, then csv, into a numpy array.
Third file to run is then `reconfig_pre.py`. This will reallocate the pre-synaptic neurons according to a distance-dependence.
Fourth file to run is `compute_stats.py`. This will compute topological statistics for the network.
