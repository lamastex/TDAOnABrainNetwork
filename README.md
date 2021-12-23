# Topological Data Analysis on a Brain Network

This is the repository containing the work done by Kieran Barber for his Masters Thesis at The Department of Mathematics, Uppsala University, Uppsala, Sweden, under the supervision of Raazesh Sainudiin.

The work resulted from discussions with several researchers, including, Kathryn Hess, Svante Janson, Wojciech Chachólski, Martina Scolamiero, and Michael Reimann.

# Requirements

To run the code within, you will need the following:
* cons_loc_pathway files: https://bbp.epfl.ch/nmc-portal/downloads.html (You want to download the average files which are located at the foot of this page under a file named `average_full.tar`, save the folder into `data` and untar there)

# Models

Running of the scripts can be done by locating and running the script `all_scripts.sh` from the `my_model` folder. Options for which model are given upon calling the script.

To run some of the computations for the topological statistics there is a system requirement of RAM > 11gb.

A further note, to run any particular statistic, you are just required to uncomment out the relevant statistic in `my_model/src/statistics.py`.

# Citation

```
@misc{Barber_Random_Graph_Models,
author = {Barber, Kieran and Sainudiin, Raazesh},
license = {Unlicense},
title = {{Random graph models of a neocortical column in a rat's brain and their topological statistical distributions}},
howpublished={\url{https://github.com/lamastex/TDAOnABrainNetwork}}
}
```
# Acknowledgements

This research was partially supported by the Wallenberg AI, Autonomous Systems and Software Program funded by Knut and Alice Wallenberg Foundation. 
The distributed computing infrastructure for this pjoject was supported by Databricks University Alliance with AWS credits.
