# working-manuscript-TopologicalDataAnalysisOnABrainNetwork
working-manuscript-TopologicalDataAnalysisOnABrainNetwork - for Kieran Barber

# neuron_swapping (Geneal Biol model)
Script takes the h5 files which can be obtained from BBP website under downloads.
The h5 files contain information about the neurons, locations, connectivity etc.
3D neuronal locations are read in for pre- and post- synaptic neurons
Their Euclidean distances are computed.
These are then binned according to their pairwise distances
Original connection matrix is read in
These values are then shuffled according to which bin they are in.

# reconstructed_file (General Biol model)
