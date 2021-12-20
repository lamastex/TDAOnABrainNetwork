m_type = [
    'L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC', 'L23_BP',
    'L23_BTC', 'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC',
    'L23_PC', 'L23_SBC', 'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC',
    'L4_MC', 'L4_NBC', 'L4_NGC', 'L4_PC', 'L4_SBC', 'L4_SP', 'L4_SS', 'L5_BP',
    'L5_BTC', 'L5_ChC', 'L5_DBC', 'L5_LBC', 'L5_MC', 'L5_NBC', 'L5_NGC', 'L5_SBC',
    'L5_STPC', 'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC', 'L6_BP', 'L6_BPC', 'L6_BTC',
    'L6_ChC', 'L6_DBC', 'L6_IPC', 'L6_LBC', 'L6_MC', 'L6_NBC', 'L6_NGC', 'L6_SBC',
    'L6_TPC_L1', 'L6_TPC_L4', 'L6_UTPC'
    ]

layer1 = ['L1_DAC', 'L1_DLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_NGC-SA', 'L1_SLAC']

layer2 = [
  'L23_BP','L23_BTC', 'L23_ChC', 'L23_DBC', 'L23_LBC', 'L23_MC', 'L23_NBC', 'L23_NGC',
  'L23_PC', 'L23_SBC'
]

layer4 = [
  'L4_BP', 'L4_BTC', 'L4_ChC', 'L4_DBC', 'L4_LBC',
  'L4_MC', 'L4_NBC', 'L4_NGC', 'L4_PC', 'L4_SBC', 'L4_SP', 'L4_SS'
]

layer5 = [
  'L5_BP','L5_BTC', 'L5_ChC', 'L5_DBC', 'L5_LBC', 'L5_MC', 'L5_NBC', 'L5_NGC', 'L5_SBC',
  'L5_STPC', 'L5_TTPC1', 'L5_TTPC2', 'L5_UTPC'
]

layer6 = [
  'L6_BP', 'L6_BPC', 'L6_BTC','L6_ChC', 'L6_DBC', 'L6_IPC', 'L6_LBC', 'L6_MC',
  'L6_NBC', 'L6_NGC', 'L6_SBC','L6_TPC_L1', 'L6_TPC_L4', 'L6_UTPC'
]

layer_name = ['L1', 'L23', 'L4', 'L5', 'L6']


''' ER vars '''
erModel, erModelName  = 'Erdos_Renyi_probability', 'Erdős–Rényi'
erFolder, erModelType = 'ER', 'Proportion of Connections in Block'
erModel2, erModelType2 = 'Erdos_Renyi', 'Connections per Block'
erTest = 'er_bwed'


''' GB vars '''
gbModel, gbFolder = 'layer5dens', 'GB'
gbModelName, gbModelType = 'GB', 'Proportion of Connections in Block'
gbModel2 = 'layer5counts'
gbModel3 = 'morph55'
gbTest = 'gb_bwed'


''' conf vars '''
cFolder = 'configuration'
cModel = 'configuration_densities'
cModelName, cModelType = 'Configuration', 'Proportion of Connections in Block'
cModel2, cModelType2 = 'configuration', 'Connections per Block'
cTest = 'c_bwed'


''' GC vars '''
gcFolder = 'GC'
gcModel = 'GC_densities'
gcModelName, gcModelType = 'GC', 'Proportion of Connections in Block'
gcModel2, gcModelType2 = 'GC', 'Connections per Block'
gcTest = 'gc_bwed'


''' BC vars '''
bcFolder = 'BC'
bcModel = 'Block_Configuration'
bcModelName, bcModelType = 'BC', 'Proportion of Connections in Block'
bcModel2, bcModelType2 = 'Block', 'Connections per Block'
bcTest = 'bc_bwed'


''' BGC vars '''
bgcFolder = 'BGC'
bgcModel = 'GBC_densities'
bgcModelName, bgcModelType = 'BGC', 'Proportion of Connections in Block'
bgcModel2, bgcModelType2 = 'Block', 'Connections per Block'
bgcTest = 'bgc_bwed'


''' Bio-M vars '''
bioFolder = 'Bio-M'
bioModel, bioModel2 = 'BioM_densities', 'BioM'
bioModelName  = 'Bio-M'
bioModelType, bioModelType2 = 'Proportion of Connections in Block', 'Connections per Block'
bioTest = 'BioM_dens_sci'
bioModelType3 = 'Block-wise Edge Density'
