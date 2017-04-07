#!/bin/sh

########################
### SAMPLE LIST ########## 
#######################

declare -a example=("WW" "WZ" "HNMupMup_100" "DYJets")

declare -a tmplist=('WpWp_qcd_madgraph' 'ZG_llG_MCatNLO' 'ZZ_llnunu_powheg' 'ZZ_llqq_MCatNLO' 'ZZ_llll_MCatNLO' 'ZZ_llll_powheg' 'ZZ_pythia8' 'ttHnobb_Powheg' 'ttHtobb_Powheg')

declare -a tmpall_mc=('TTJets_aMC' 'WJets' 'WW'  'WZ' 'ZZ' 'DYJets')

declare -a hn=('DYJets_10to50'  'DYJets' 'WW' 'ZZ' 'WZ' 'TTJets_MG')

declare -a pu_dilepton_list=('DYJets_10to50' 'DYJets' 'WJets' 'TT_powheg'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t'  'SingleTbar_tW' 'SingleTop_tW' 'WGtoLNuG'  'ZGto2LG' 'WZTo3LNu_powheg' 'ZZTo4L_powheg' 'DYJets_MG_10to50' 'DYJets_MG' 'TTJets_aMC' )

#declare -a ch_cb=('TT_powheg' 'CHToCB_M125_madgraph_13TeV_2016')
#declare -a ch_cb=('TT_powheg' 'CHToCB_M090_madgraph_13TeV_2016' 'CHToCB_M100_madgraph_13TeV_2016' 'CHToCB_M110_madgraph_13TeV_2016' 'CHToCB_M120_madgraph_13TeV_2016' 'CHToCB_M125_madgraph_13TeV_2016' 'CHToCB_M130_madgraph_13TeV_2016' 'CHToCB_M140_madgraph_13TeV_2016' 'CHToCB_M150_madgraph_13TeV_2016')
declare -a ch_cb=('TT_powheg')