sktree -a PileupValidation -list pu_dilepton_list -s SKTree_DiLepSkim
sktree -a PileupValidation  -S DoubleEG -s SKTree_DiLepSkim  -n 30 -F
#sktree -a PileupValidation  -S SingleElectron  -n 50  -F
#sktree -a PileupValidation -i DYJets -s SKTree_DiLepSkim  -n 30 -F
#sktree -a PileupValidation -list qcd_eg -n 20 -F
#sktree -a PileupValidation -S DoubleEG -fake True -s SKTree_DiLepSkim  -n 100
#sktree -a PileupValidation  -S DoubleMuon -s SKTree_DiLepSkim  -n 100