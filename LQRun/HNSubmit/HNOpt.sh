sktree -a HNDiElectronOptimisation -list hn_ee_sigcf -n 100 -flip True -s SKTree_TriLepSkim  -F
#sktree -a HNDiElectron -list qcd_eg -n 20 -F
sktree -a HNDiElectronOptimisation -s FLATCAT -list hn_ll_ee1 
sktree -a HNDiElectronOptimisation -S DoubleEG -fake True -s SKTree_HNDiLepSkim  -n 200
