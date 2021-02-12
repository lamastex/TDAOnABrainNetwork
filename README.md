# working-manuscript-TopologicalDataAnalysisOnABrainNetwork
working-manuscript-TopologicalDataAnalysisOnABrainNetwork - for Kieran Barber

# Requirements
To run the code within, you will need the following:
* pyflagser: https://pypi.org/project/pyflagser/
* cons_loc_pathway files: https://bbp.epfl.ch/nmc-portal/downloads.html (You want to download the average files which are located at the foor of this page under a file name `average_full.tar`, save the folder into `data` and untar there)

# Empirical Geometric Configuration model (4 steps)
There are four scripts to run, however, fortunately, this can all be run with one run from the shell script. Before running any of these scripts, however, you will need to download the pathway files (average files).
* To run, you have the following from command line: `sh configuration.sh`. There will be a message that informs of computation time and some idea on the input of the arguments to follow.
* Once the average files are downloaded, the first script to run is `h5_tocsv.py`. This will convert the h5 files into csv files and will put them in the `reconstruction` folder. 
* Second file to run is `concatenate_csv.py`. This will return the array that was in h5, then csv, into a numpy array.
* Third file to run is then `reconfig_pre.py`. This will reallocate the pre-synaptic neurons according to a distance-dependence.
* Fourth file to run is `compute_stats.py`. This will compute topological statistics for the network.
