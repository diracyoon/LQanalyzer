#!/bin/sh

########################
### SAMPLE LIST ########## 
#######################

declare -a trilep_fake_bkg_short=(
'DYJets_10to50'
'SingleTop_s' 'SingleTop_t' 'SingleTbar_t' 'SingleTop_tW' 'SingleTbar_tW'
'WZ' 'ZZ' 'WW'
'ZGto2LG' 'WGtoLNuG'
)
declare -a trilep_fake_bkg_long=('DYJets' 'WJets' 'TT_powheg')

declare -a  qcd_mu_prompt=('QCD_Pt-1000toInf_MuEnriched' 'QCD_Pt-120to170_MuEnriched' 'QCD_Pt-15to20_MuEnriched' 'QCD_Pt-170to300_MuEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_Pt-300to470_MuEnriched' 'QCD_Pt-30to50_MuEnriched' 'QCD_Pt-470to600_MuEnriched' 
'QCD_Pt-50to80_MuEnriched' 'QCD_Pt-600to800_MuEnriched' 'QCD_Pt-800to1000_MuEnriched' 'QCD_Pt-80to120_MuEnriched' 'TT_powheg' 'WJets'  'DYJets')
declare -a hnfail=('HNMumEp_50' 'HNMumEp_500' 'HNMumEp_60' 'HNMumMum_100' 'HNMumMum_1100' 'HNMumMum_1500' 'HNMumMum_200' 'HNMumMum_40' 'HNMumMum_50' 'HNMumMum_500' 'HNMumMum_60' 'HNMoriondLL_Tchannel_MumMum_100' 'HNMoriondLL_Tchannel_MumMum_1100' 'HNMoriondLL_Tchannel_MumMum_200' 'HNMoriondLL_Tchannel_MumMum_500' 'HNMumMup_100' 'HNMumMup_1100' 'HNMumMup_1500' 'HNMumMup_200' 'HNMumMup_40' 'HNMumMup_50' 'HNMumMup_500' 'HNMumMup_60' 'HNMupEm_100' 'HNMupEm_1100' 'HNMupEm_1500' 'HNMupEm_200' 'HNMupEm_40' 'HNMupEm_50' 'HNMupEm_500' 'HNMupEm_60' 'HNMupEp_100')

declare -a dilepton_list2=('DYJets_10to50' 'DYJets' 'WJets' 'WpWpEWK' 'WpWpQCD' 'TT_powheg'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t'  'SingleTbar_tW' 'SingleTop_tW' 'WWW' 'ttW' 'ttZ' 'ttH_nonbb' 'ttH_bb' 'ZZZ' 'WZZ' 'WWZ' 'VBF_HToMuMu' 'WGtoLNuG'  'ZGto2LG' 'WZTo3LNu_powheg' 'ZZTo4L_powheg' 'WWTo2L2Nu' 'WWToLNuQQ' 'TG' 'TTG' 'ggHtoWW' 'ggHtoZZ' 'vbfHtoWW' 'tZq'  'ttZToLL_M-1to10' 'WZG' 'WWG' 'ggZZto4e' 'ttWToLNu' 'ttZToLL_M-10'  'WZto2L2Q_amcatnlo' 'ZZTo2L2Nu_Powheg' 'ZZTo2L2Q_Powheg' 'ggZZto2e2mu' 'ggZZto2mu2nu' 'ggZZto2mu2tau'  'ggZZto4mu' 'ggWWto2L2Nu' 'ww_ds' 'ZZTo4L_amcatnlo' 'WZ' 'WW' 'ZZ' 'TTJets_aMC') 



declare -a qcd_ee=('qcd_15to20_bctoe' 'qcd_170to250_bctoe' 'qcd_20to30_bctoe' 'qcd_250toinf_bctoe' 'qcd_30to80_bctoe' 'qcd_80to170_bctoe' 'QCD_Pt-120to170_EMEnriched' 'QCD_Pt-170to300_EMEnriched' 'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-300toInf_EMEnriched'  'QCD_Pt-30to50_EMEnriched' 'QCD_Pt-50to80_EMEnriched' 'QCD_Pt-80to120_EMEnriched')
declare -a qcd_eemm=('qcd_15to20_bctoe' 'qcd_170to250_bctoe' 'qcd_20to30_bctoe' 'qcd_250toinf_bctoe' 'qcd_30to80_bctoe' 'qcd_80to170_bctoe' 'QCD_Pt-120to170_EMEnriched' 'QCD_Pt-170to300_EMEnriched' 'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-300toInf_EMEnriched'  'QCD_Pt-30to50_EMEnriched' 'QCD_Pt-50to80_EMEnriched' 'QCD_Pt-80to120_EMEnriched' 'WJets' 'TTJets_aMC' 'QCD_Pt-1000toInf_MuEnriched' 'QCD_Pt-120to170_MuEnriched' 'QCD_Pt-15to20_MuEnriched' 'QCD_Pt-170to300_MuEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_Pt-300to470_MuEnriched' 'QCD_Pt-30to50_MuEnriched' 'QCD_Pt-470to600_MuEnriched' 'QCD_Pt-50to80_MuEnriched' 'QCD_Pt-600to800_MuEnriched' 'QCD_Pt-800to1000_MuEnriched' 'QCD_Pt-80to120_MuEnriched' 'TT_powheg')

declare -a hn_ll_mm_1=( 'HNMumMum_100' 'HNMumMum_1000' 'HNMumMum_1100' 'HNMumMum_1200' 'HNMumMum_125' 'HNMumMum_1300' 'HNMumMum_1400' 'HNMumMum_150' 'HNMumMum_1500' 'HNMumMum_200' 'HNMumMum_250' 'HNMumMum_300' 'HNMumMum_40' 'HNMumMum_400' 'HNMumMum_50' 'HNMumMum_500' 'HNMumMum_60' 'HNMumMum_600' 'HNMumMum_70' 'HNMumMum_700' 'HNMumMum_80' 'HNMumMum_800' 'HNMumMum_90' 'HNMumMum_900' 'HNMumMum_Tchannel_1000' 'HNMoriondLL_Tchannel_MumMum_100' 'HNMoriondLL_Tchannel_MumMum_1100' 'HNMumMum_Tchannel_1200' 'HNMumMum_Tchannel_1300' 'HNMumMum_Tchannel_1400' 'HNMumMum_Tchannel_1500' 'HNMoriondLL_Tchannel_MumMum_200' 'HNMumMum_Tchannel_300' 'HNMumMum_Tchannel_400' 'HNMoriondLL_Tchannel_MumMum_500' 'HNMumMum_Tchannel_600' 'HNMumMum_Tchannel_700' 'HNMumMum_Tchannel_800' 'HNMumMum_Tchannel_900'  'HNMupMup_100' 'HNMupMup_1000' 'HNMupMup_1100' 'HNMupMup_1200' 'HNMupMup_125' 'HNMupMup_1300' 'HNMupMup_1400' 'HNMupMup_150' 'HNMupMup_1500' 'HNMupMup_200' 'HNMupMup_250' 'HNMupMup_300' 'HNMupMup_40' 'HNMupMup_400' 'HNMupMup_50' 'HNMupMup_500' 'HNMupMup_60' 'HNMupMup_600' 'HNMupMup_70' 'HNMupMup_700' 'HNMupMup_80' 'HNMupMup_800' 'HNMupMup_90' 'HNMupMup_900' 'HNMupMup_Tchannel_1000' 'HNMoriondLL_Tchannel_MupMup_100' 'HNMoriondLL_Tchannel_MupMup_1100' 'HNMupMup_Tchannel_1200' 'HNMupMup_Tchannel_1300' 'HNMupMup_Tchannel_1400' 'HNMupMup_Tchannel_1500' 'HNMoriondLL_Tchannel_MupMup_200' 'HNMupMup_Tchannel_300' 'HNMupMup_Tchannel_400' 'HNMoriondLL_Tchannel_MupMup_500' 'HNMupMup_Tchannel_600' 'HNMupMup_Tchannel_700' 'HNMupMup_Tchannel_800' 'HNMupMup_Tchannel_900' )


declare -a example=("WW" "WZ" "HNMupMup_100" "DYJets")

declare -a tmplist=('ggZZto4e' 'ttZToLL_M-1to10' 'ttWToLNu' 'ttZToLL_M-10' 'WZG' 'WWG')

declare -a tmpall_mc=('TTJets_aMC' 'LowStat_WJets' 'WW'  'WZ' 'ZZ' 'LowStat_DYJets' 'DYJets_10to50' )

declare -a hn=('DYJets_10to50'  'DYJets' 'TT_powheg' 'WJets' 'WGtoLNuG'  'ZGto2LG' 'qcd_15to20_bctoe' 'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-20to30_MuEnriched')
declare -a hn_fake=('DYJets' 'TT_powheg' 'WJets' 'qcd_15to20_bctoe' 'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_DoubleEMEnriched_30-40_mgg80toinf' 'QCD_DoubleEMEnriched_30-inf_mgg40to80' 'QCD_DoubleEMEnriched_40-inf_mgg80toinf')

declare -a pu_dilepton_list=('DYJets_10to50' 'DYJets' 'WJets' 'TT_powheg'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t'  'SingleTbar_tW' 'SingleTop_tW' 'WGtoLNuG'  'ZGto2LG' 'WZTo3LNu_powheg' 'ZZTo4L_powheg' 'DYJets_MG_10to50' 'DYJets_MG' 'TTJets_aMC' )

declare -a hntmp=('TTTT' 'TG' 'TTG' 'ttWToLNu' 'ttZToLL_M-1to10' 'ttZToLL_M-10'  'tZq' 'ggHtoWW' 'ggHtoZZ' 'WWG' 'WZG' 'WZto2L2Q_amcatnlo' 'ZZTo2L2Nu_Powheg' 'ZZTo2L2Q_Powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau'  'ggZZto4e' 'ggWWto2L2Nu' 'ww_ds'  )

declare -a hn_eetmp=('DYJets_10to50' 'DYJets' 'WJets' 'WpWpEWK' 'WpWpQCD' 'TT_powheg'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t'  'SingleTbar_tW' 'SingleTop_tW' 'WWW' 'ttW' 'ttZ' 'ttH_nonbb' 'ttH_bb' 'ZZZ' 'WZZ'  'VBF_HToMuMu' 'WGtoLNuG'  'ZGto2LG' 'WZTo3LNu_powheg' 'ZZTo4L_powheg'  'WWTo2L2Nu' 'WWToLNuQQ' 'QCD_DoubleEMEnriched_30-40_mgg80toinf' 'QCD_DoubleEMEnriched_30-inf_mgg40to80' 'QCD_DoubleEMEnriched_40-inf_mgg80toinf' 'TG' 'TTG' 'ttWToLNu' 'ttZToLL_M-1to10' 'ttZToLL_M-10'  'tZq' 'ggHtoWW' 'ggHtoZZ' 'WWG' 'WZG' 'WZto2L2Q_amcatnlo' 'ZZTo2L2Nu_Powheg' 'ZZTo2L2Q_Powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau'  'ggZZto4e' 'ggWWto2L2Nu' 'ww_ds'  )
declare -a hn_mm_all=('DYJets_10to50' 'DYJets' 'WJets' 'WpWpEWK' 'WpWpQCD' 'TT_powheg'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t'  'SingleTbar_tW' 'SingleTop_tW' 'WWW' 'ttW' 'ttZ' 'ttH_nonbb' 'ttH_bb' 'ZZZ' 'WZZ'  'VBF_HToMuMu' 'WGtoLNuG'  'ZGto2LG' 'WZTo3LNu_powheg' 'ZZTo4L_powheg'  'WWTo2L2Nu' 'WWToLNuQQ' 'TG' 'TTG' 'ttWToLNu' 'ttZToLL_M-1to10' 'ttZToLL_M-10'  'tZq' 'ggHtoWW' 'ggHtoZZ' 'WWG' 'WZG' 'WZto2L2Q_amcatnlo' 'ZZTo2L2Nu_Powheg' 'ZZTo2L2Q_Powheg' 'ggZZto2e2mu' 'ggZZto2mu2nu' 'ggZZto2mu2tau'  'ggZZto4mu' 'ggWWto2L2Nu' 'ww_ds'  )


declare -a sktmp=('WJets' 'TT_powheg')
declare -a vv=('WZTo3LNu_powheg' 'ZZTo4L_powheg' 'WGtoLNuG'  'ZGto2LG' 'ZZTo2L2Nu_Powheg' 'ZZTo2L2Q_Powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau'  'ggZZto4e' 'ggWWto2L2Nu' 'ww_ds' 'TG' 'TTG' 'ttWToLNu')

declare -a hn_ee_sig=('WpWpEWK' 'WpWpQCD'  'ZZZ' 'WZZ'  'ww_ds'  'ggZZto4e' 'WZTo3LNu_powheg' 'ZZTo4L_powheg'  'ggHtoZZ'  'WWG' 'WZG'    'ttWToLNu' 'ttZToLL_M-1to10' 'ttZToLL_M-10'  'tZq' 'ggHtoWW' )
declare -a hn_ee_sigcf=('DYJets' 'TT_powheg' 'WWTo2L2Nu' )

declare -a hn_ee_type=('DYJets' 'TT_powheg' 'WJets' 'ZGto2LG' 'QCD_Pt-30to50_EMEnriched')


declare -a tmp=(
'DYJets_MG'
'GG_HToMuMu' 
'SingleTop_t'
'TTJets_aMC'
'TTLL_powheg'
'ttZToLL_M-10'
'WWTo2L2Nu_DS'
'ZGto2LG' )

declare -a tmp2=(
'HNEpMup_100' 
'HNEpMup_1100' 
'HNEpMup_1500' 
'HNEpMup_200' 
'HNEpMup_40' 
'HNEpMup_500' 
'HNEpMup_50' 
'HNEpMup_60' 
'HNMoriondLL_Tchannel_EpMup_100' 
'HNMoriondLL_Tchannel_EpMup_1100' 
'HNMoriondLL_Tchannel_EpMup_200' 
'HNMoriondLL_Tchannel_EpMup_500' 
'HN_MuMuMu_1000'
'HN_MuMuMu_10'
'HN_MuMuMu_150'
'HN_MuMuMu_20'
'HN_MuMuMu_200'
 )

declare -a ch_tt=('TT_powheg')
declare -a ch_cb_125=('CHToCB_M090_madgraph_13TeV_2016')
declare -a ch_test=('TT_powheg' 'CHToCB_M125_madgraph_13TeV_2016')
declare -a ch_low=('CHToCB_M090_madgraph_13TeV_2016' 'CHToCB_M100_madgraph_13TeV_2016')
declare -a ch_high=('TT_powheg' 'CHToCB_M110_madgraph_13TeV_2016' 'CHToCB_M120_madgraph_13TeV_2016' 'CHToCB_M125_madgraph_13TeV_2016' 'CHToCB_M130_madgraph_13TeV_2016' 'CHToCB_M140_madgraph_13TeV_2016' 'CHToCB_M150_madgraph_13TeV_2016')
declare -a ch_cb_all=('CHToCB_M090_madgraph_13TeV_2016' 'CHToCB_M100_madgraph_13TeV_2016' 'CHToCB_M110_madgraph_13TeV_2016' 'CHToCB_M120_madgraph_13TeV_2016' 'CHToCB_M125_madgraph_13TeV_2016' 'CHToCB_M130_madgraph_13TeV_2016' 'CHToCB_M140_madgraph_13TeV_2016' 'CHToCB_M150_madgraph_13TeV_2016')
declare -a ch_all=('TT_powheg' 'CHToCB_M090_madgraph_13TeV_2016' 'CHToCB_M100_madgraph_13TeV_2016' 'CHToCB_M110_madgraph_13TeV_2016' 'CHToCB_M120_madgraph_13TeV_2016' 'CHToCB_M125_madgraph_13TeV_2016' 'CHToCB_M130_madgraph_13TeV_2016' 'CHToCB_M140_madgraph_13TeV_2016' 'CHToCB_M150_madgraph_13TeV_2016')

declare -a new_list=('DYJets' 'LowStatDYJets' )
