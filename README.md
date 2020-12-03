# working-manuscript-TopologicalDataAnalysisOnABrainNetwork
working-manuscript-TopologicalDataAnalysisOnABrainNetwork - for Kieran Barber

# Requirements
To run the code within, you will need the following:
* flagser: https://github.com/luetge/flagser
* pyflagser: https://pypi.org/project/pyflagser/
* pyflagsercontain: pip install pyflagsercontain
* cons_loc_pathway files: https://bbp.epfl.ch/nmc-portal/downloads (average files)

## General_biol
From https://bbp.epfl.ch/nmc-portal/downloads, we can download the average file
and store this in the data directory. Then run: `python general_biol_swap.py` (~10 mins runtime)

This saves csv files in the reconstruction folder to be used in the next step.

Once completed, we now run:`python3 general_biol_reconstruction.py` (~15 mins)

This returns output similar to the following: 31346, 7648079, 34103325, 5158826, 43523, 30, 0

## Bio-M
### Using Python
Making sure the reconstruction folder is empty, we can run the following: `python original.py` (~5 mins)

This does exactly as `general_biol_swap.py` apart from the swapping

Rerun the following: `python3 general_biol_reconstruction.py` (~15 mins)

This returns output similar to the following: 31346, 7822274, 76361228, 64064185, 7274386, 156404, 896

### Using Flagser
Alternatively, if we have the flagser package saved in `src` directory, we can 
do the following:

1. Convert the h5 files to .flag files: `./tools/h5toflagser ../data/(h5 file) MC6.flag`

2. Run the following: `./flagser-count MC6.flag --out file.h5`

This will return the same as above.

## Inhibitory and Excitatory neurons in layers (Figure2.py)
For the second and third graphs in `Figure2.py`, we need the relevant 
groups for inhibitory and excitatory neurons and for each layer.

inhib_L1: 'L1_DAC,L1_DLAC,L1_HAC,L1_NGC-DA,L1_NGC-SA,L1_SLAC'
inhib_L23:'L23_BTC,L23_DBC,L23_LBC,L23_MC,L23_NBC,L23_NGC,L23_SBC'
inhib_L4: 'L4_BTC,L4_DBC,L4_LBC,L4_MC,L4_NBC,L4_NGC,L4_SBC'
inhib_L5: 'L5_BTC,L5_DBC,L5_LBC,L5_MC,L5_NBC,L5_NGC'
inhib_L6: 'L6_BTC,L6_MC,L6_LBC'

excit_L23: 'L23_PC'
excit_L4: 'L4_PC,L4_SP,L4_SS'
excit_L5: 'L5_STPC,L5_TTPC1,L5_TTPC2,L5_UTPC'
excit_L6: 'L6_BPC,L6_IPC,L6_TPC_L1,L6_TPC_L4,L6_UTPC'

Again, to run, first convert relevant files to .flag: `./tools/h5toflagser ../data/((h5 file) --groups "L1_DAC,L1_DLAC,L1_HAC,L1_NGC-DA,L1_NGC-SA,L1_SLAC" MC6_L1.flag`

Then we can run: `./flagser-count MC6_L1.flag --out file.h5`

This returns something like the following: 339, 848, 92, 1

## Inhibitory and Excitatory neurons (Figure3.py A graphs)
This is the same as above but for the following groups:

inhib:'L1_DAC,L1_DLAC,L1_HAC,L1_NGC-DA,L1_NGC-SA,L1_SLAC,L23_BTC,L23_DBC,
    L23_LBC,L23_MC,L23_NBC,L23_NGC,L23_SBC,L4_BTC,L4_DBC,L4_LBC,
    L4_MC,L4_NBC,L4_NGC,L4_SBC,L5_BTC,L5_DBC,L5_LBC,L5_MC,L5_NBC,
    L5_NGC,L6_BTC,L6_MC,L6_LBC'

excit: 'L23_PC,L4_PC,L4_SP,L4_SS,L5_STPC,L5_TTPC1,L5_TTPC2,L5_UTPC,L6_BPC,
    L6_IPC,L6_TPC_L1,L6_TPC_L4,L6_UTPC'

## 3D simplices - threeDsimplices and Figure3b
threedsimplices.py runs in conjunction with figure 3b. We first run: `python threesimplices.py`

This will produce csv files in output. Then running: `python3 Figure3.py`

These values will be read in and produce a graph.

## Mean maximal Dimension - MMD_neuron

To run this file: `python MMD_neuron.py`

Output is something similar to:

>3.3392330383480826,
>4.5146315509444,
>4.862046771079168,
>4.9274484113986246,
>5.017607294450558

These values are included in Figure3.py graph 3c.

## Figure3
If all the above has been run, this should be run with: `python3 Figure3.py`

## Figure5
Can be run with: `python3 Figure5.py`
