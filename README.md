# working-manuscript-TopologicalDataAnalysisOnABrainNetwork
working-manuscript-TopologicalDataAnalysisOnABrainNetwork - for Kieran Barber

# Requirements
To run the code within, you will need the following:
* cons_loc_pathway files: https://bbp.epfl.ch/nmc-portal/downloads.html (You want to download the average files which are located at the foot of this page under a file named `average_full.tar`, save the folder into `data` and untar there)

# Models
Running of the scripts can be done by locating and running the script `all_scripts.sh` from the `my_model` folder. Options for which model are given upon calling the script. ALternatively, to get a run of all 7 models, you can simply run `auto_run.sh`.

To run some of the computations for the topological statistics there is a system requirement of RAM > 11gb.

A further note, to run any particular statistic, you are just required to uncomment out the relevant statistic in `my_model/src/statistics.py`.
