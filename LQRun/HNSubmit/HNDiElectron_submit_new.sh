sktree -a HNDiElectron  -S DoubleEG -s SKTree_DiLepSkim  -n 30 -F
#sktree -a HNDiElectron -i DYJets -s SKTree_DiLepSkim   -n 100 -tau True
#mv /data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer//data/output/CAT/HNDiElectron/periodBtoH/HNDiElectron_SKDYJets_hndilep_cat_v8-0-7.root /data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer//data/output/CAT/HNDiElectron/periodBtoH/HNDiElectron_SKDYTauTau_dilep_cat_v8-0-7.root
#sktree -a HNDiElectron -i DYJets -s SKTree_DiLepSkim  -n 30 -F
sktree -a HNDiElectron -list hn_eetmp -s SKTree_DiLepSkim
#sktree -a HNDiElectron -list vv -s SKTree_DiLepSkim  -F 200
#sktree -a HNDiElectron  -S DoubleEG -s SKTree_DiLepSkim  -n 30 -fake True -F
#mv  /data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer//data/output/CAT/HNDiElectron/periodBtoH//HNDiElectron_DoubleEG_SKnonprompt_dilep_cat_v8-0-7.root  /data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer//data/output/CAT/HNDiElectron/periodBtoH//HNDiElectron_SKDoubleEG_SKnonprompt_dilep_cat_v8-0-7.root 
#sktree -a HNDiElectron -list qcd_eg -n 20 -F
#sktree -a HNDiElectron -S DoubleEG -fake True -s SKTree_DiLepSkim  -n 100
#sktree -a HNDiElectron  -S DoubleMuon -s SKTree_DiLepSkim  -n 100