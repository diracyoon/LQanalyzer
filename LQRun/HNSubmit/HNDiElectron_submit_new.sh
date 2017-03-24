sktree -a HNDiElectron  -S DoubleEG -s SKTree_DiLepSkim  -n 30 -F
#sktree -a HNDiElectron  -S SingleElectron  -n 50  -F
#sktree -a HNDiElectron -i DYJets -s SKTree_DiLepSkim  -n 30 -F
sktree -a HNDiElectron -list dilepton_list -s SKTree_DiLepSkim
#sktree -a HNDiElectron -list qcd_eg -n 20 -F
#sktree -a HNDiElectron -S DoubleEG -fake True -s SKTree_DiLepSkim  -n 100
#sktree -a HNDiElectron  -S DoubleMuon -s SKTree_DiLepSkim  -n 100