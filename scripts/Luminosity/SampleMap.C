#ifndef SampleMap_C
#define SampleMap_C

#include <map>

map<TString, TString>  GetLQMap(TString listpath);
map<TString, TString>  GetLQMap2016(TString listpath, bool sig);
map<TString, Double_t>  GetXSecMap(TString listpath);
map<TString, Double_t>  GetXSecMap2016(TString listpath, bool sig);
map<TString, TString>  GetMissingMap(TString cversion,TString listpath, bool sig);
vector<TString>  GetAvailableMap(TString cversion,TString listpath, bool sig);
map<TString, TString>  GetDatasetNames(TString cversion,TString listpath, bool sig);
map<TString, TString>  GetTriLepMap(TString listpath);
map<TString, TString>  GetTriLepMap2016(TString listpath, bool sig);

bool CheckMaps(TString listpath,bool sig);


bool IsSignal(TString samplename){
  
  if(samplename.Contains("HN")) return true;
  if(samplename.Contains("CH")) return true;
  if(samplename.Contains("TTToHc")) return true;
  return false;
}

bool CheckMaps(TString listpath, bool sig){
  
  // Function checks the initialisation of the maps are consistent
  map<TString,  Double_t> mapxs = GetXSecMap2016(listpath, sig);
  map<TString,  TString>  maplq = GetLQMap2016(listpath, sig);

  bool failcheck=false;
  if (mapxs.size() != maplq.size()) {
    failcheck=true;
    cout << "Maps are not the same size" << endl;
  }
  for(std::map<TString, Double_t>::iterator mit =mapxs.begin(); mit != mapxs.end();++mit){
    std::map<TString, TString>::iterator mit2 =maplq.find(mit->first);
    if(mit2 == maplq.end()){
      failcheck=true;
      cout << "Samples " << mit->first << " not in lqmap " << endl;
    }
  }
  
  for( std::map<TString, TString>::iterator mit = maplq.begin(); mit != maplq.end();++mit){
    std::map<TString, Double_t>::iterator mit2 = mapxs.find(mit->first);
    if(mit2 == mapxs.end()){
      failcheck=true;
      cout << "Samples " << mit->first << " not in xsmap " << endl;
    }
  }
  
  return failcheck;
}

map<TString,  Double_t>  GetXSecMap2016(TString datasetfile, bool sig){
  map<TString,  Double_t> dirmap;

  ifstream datasetname_file(datasetfile.Data());
  if(!datasetname_file){
    cerr << "Did not find " + datasetfile + ", exiting ..." << endl;
    return dirmap ;
  }
  while(!datasetname_file.eof()) {
    string samplealias;
    string datasetname;
    double xsec;
    string fullname;
    datasetname_file>> samplealias;
    datasetname_file >> datasetname;
    datasetname_file >> xsec;
    datasetname_file >> fullname;
    if(samplealias=="END") break;
    if(sig&&IsSignal(samplealias))    {
      if (! datasetname.empty())  dirmap[TString(datasetname)] = xsec;
    }
    if(!sig&&!IsSignal(samplealias))    {
      if (! datasetname.empty())  dirmap[TString(datasetname)] = xsec;
    }

  }
  /*

  dirmap["DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] =18610.;
  dirmap["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"]=5765.4;           
  dirmap["DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]=18610.;
  dirmap["DYJetsToLL_M-5to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]=71310.;
  dirmap["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]=5765.4;

  dirmap["ttbb_4FS_ckm_amcatnlo_madspin_pythia8"]=13.93;
  dirmap["ttHTobb_M125_13TeV_powheg_pythia8"]=0.2934;      
  dirmap["ttHToNonbb_M125_13TeV_powheg_pythia8"]=0.2151;      

  dirmap["WW_TuneCUETP8M1_13TeV-pythia8"]=118.7;
  dirmap["WZ_TuneCUETP8M1_13TeV-pythia8"]=47.13;
  dirmap["ZZ_TuneCUETP8M1_13TeV-pythia8"]=16.523;
  dirmap["WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8"]=0.1651;
  dirmap["WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"]=0.1651;
  dirmap["WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"]=0.05565;
  dirmap["ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"]=0.01398;

  dirmap["TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"]=831.76; 
  dirmap["TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]=831.76; 
  dirmap["TT_TuneCUETP8M1_13TeV-powheg-pythia8"]=831.76;
  dirmap["TT_TuneCUETP8M1_mtop1695_13TeV-powheg-pythia8"]=831.76;
  dirmap["TT_TuneCUETP8M1_mtop1755_13TeV-powheg-pythia8"]=831.76;
  dirmap["TTTo2L2Nu_13TeV-powheg"]=87.31;
  dirmap["TTTo2L2Nu_13TeV_ScaleDown-powheg"]=87.31;
  dirmap["TTTo2L2Nu_13TeV_ScaleUp-powheg"]=87.31;
  dirmap["TTToSemiLeptonic_13TeV-powheg"]=356.4;
  dirmap["TTToSemiLeptonic_13TeV_ScaleDown-powheg"]=356.4;
  dirmap["TTToSemiLeptonic_13TeV_ScaleUp-powheg"]=356.4;

  dirmap["ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1"]=3.36; 
  dirmap["ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1"]=26.38;  
  dirmap["ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1"]=44.33;    
  dirmap["ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1"]=35.85;    
  dirmap["ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1"]=35.85;    
  
  dirmap["ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1"]=35.6;
  dirmap["ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1"]=35.6;

  dirmap["VBF_HToMuMu_M125_13TeV_powheg_pythia8"]=0.000828308;
  dirmap["GluGlu_HToMuMu_M125_13TeV_powheg_pythia8"]=0.00970632;
  
  dirmap["ttWJets_13TeV_madgraphMLM"]=0.61;    
  dirmap["ttZJets_13TeV_madgraphMLM"]=0.78;      

  dirmap["WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]=61526.7;
  dirmap["WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"]=61526.7;


  dirmap["13TeV_HN150_mumumu_width_AUTO_VmuN_0p1"]=0.0004561;
  dirmap["13TeV_HN40_mumumu_width_AUTO_VmuN_0p1"]=8.6;
  dirmap["13TeV_HN60_mumumu_width_AUTO_VmuN_0p1"]=3.416;
  dirmap["13TeV_HN700_mumumu_width_AUTO_VmuN_0p1"]=0.000002644;
  */
  return dirmap;
}

map<TString,  Double_t>  GetXSecMap(TString datasetfile){

   map<TString,  Double_t> dirmap;

   ifstream datasetname_file(datasetfile.Data());
   if(!datasetname_file){
     cerr << "Did not find " + datasetfile + ", exiting ..." << endl;
     return dirmap;
   }
   while(!datasetname_file.eof()) {
     string samplealias;
     string datasetname;
     double xsec;
     string fullname;
     datasetname_file>> samplealias;
     datasetname_file >> datasetname;
     datasetname_file >> xsec;
     datasetname_file >> fullname;
     if(samplealias=="END") break;
     if (! datasetname.empty()) dirmap[TString(datasetname)] = xsec;
   }
   /*
   dirmap["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] = 5765.4;
   dirmap["DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] = 18610. ;
   dirmap["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"] = 5765.4;
   dirmap["DYJetsToLL_M-5to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"] =71310.;
   dirmap["ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1"] = 26.38;
   dirmap["ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1"] =44.33;
   dirmap["ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1"] =35.6;
   dirmap["ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1"] =35.6;
   dirmap["ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1"] = 3.36;

   dirmap["ttHTobb_M125_13TeV_powheg_pythia8"] =0.5058; ////
   dirmap["ttHToNonbb_M125_13TeV_powheg_pythia8"] =0.5058; ////
   dirmap["GluGlu_HToMuMu_M125_13TeV_powheg_pythia8"] =0.00970632; ////
   dirmap["VBF_HToMuMu_M125_13TeV_powheg_pythia8"] =0.000828308; ////

   dirmap["TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] =831.76; 
   dirmap["TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"] =831.76 ;

   dirmap["TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8"] =0.2043;
   dirmap["TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8"] =0.4062;
   dirmap["TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] =0.2529 ;
   dirmap["TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] =0.5297;
   dirmap["TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-scaleup-pythia8"] = 831.76;
   dirmap["TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-scaledown-pythia8"] = 831.76;
   dirmap["TT_TuneCUETP8M1_13TeV-powheg-scaleup-pythia8"]=831.76;
   dirmap["TT_TuneCUETP8M1_13TeV-powheg-scaledown-pythia8"]=831.76;
   dirmap["TT_TuneCUETP8M1_13TeV-powheg-pythia8"]=831.76;
   dirmap["TT_TuneEE5C_13TeV-amcatnlo-herwigpp"]=831.76;
   dirmap["TT_TuneEE5C_13TeV-powheg-herwigpp"]= 831.76;
   dirmap["TT_TuneCUETP8M1_mtop1665_13TeV-powheg-pythia8"]=831.76;
   dirmap["TT_TuneCUETP8M1_mtop1695_13TeV-powheg-pythia8"]=831.76;
   dirmap["TT_TuneCUETP8M1_mtop1715_13TeV-powheg-pythia8"]=831.76;
   dirmap["TT_TuneCUETP8M1_mtop1735_13TeV-powheg-pythia8"]=831.76;
   dirmap["TT_TuneCUETP8M1_mtop1755_13TeV-powheg-pythia8"]=831.76;
   dirmap["TT_TuneCUETP8M1_mtop1785_13TeV-powheg-pythia8"]=831.76;
   dirmap["TT_TuneZ2star_13TeV-powheg-pythia6-tauola"]=831.76;
   dirmap["TTbarDMJets_pseudoscalar_Mchi-10_Mphi-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]=0.1901;
   dirmap["TTbarDMJets_pseudoscalar_Mchi-10_Mphi-10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]=0.01499;

   dirmap["WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] =61526.7;

   dirmap["WW_TuneCUETP8M1_13TeV-pythia8"] =113.826;
   dirmap["WZ_TuneCUETP8M1_13TeV-pythia8"] =47.13;
   dirmap["ZZ_TuneCUETP8M1_13TeV-pythia8"] =16.91;

   dirmap["ZZTo4L_13TeV_powheg_pythia8"] = 1.256;
   dirmap["ZZTo4L_13TeV-amcatnloFXFX-pythia8"] =1.212 ;
   dirmap["ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8"] = 3.22;
   dirmap["ZZTo2L2Nu_13TeV_powheg_pythia8"] = 0.564;
   dirmap["WWTo2L2Nu_13TeV-powheg"] = 12.178;
   dirmap["WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8"] = 5.595;
   dirmap["WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8"] = 4.42965;
   dirmap["WZJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] = 5.26;
   dirmap["WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8"] = 0.02064;
   dirmap["WpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8"] = 0.01538;
   dirmap["WW_DoubleScattering_13TeV-pythia8"] = 1.64;
   dirmap["WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] = 0.05565;
   
   dirmap["ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] = 0.01398 ;
   dirmap["WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] = 0.1651; 
   dirmap["QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] = 3819570.0;
   dirmap["QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] = 558528000*0.0053;
   dirmap["QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] =139803000*0.01182;
   dirmap["QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] =19222500*0.02276;
   dirmap["QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] =2758420*0.03844;
   dirmap["QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] =469797*0.05362;
   dirmap["QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] =117989*0.07335;
   dirmap["QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] =7820.25*0.10196;
   dirmap["QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] =645.528*0.12242;
   dirmap["QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] =187.109*0.13412;
   dirmap["QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] =32.3486*0.14552;
   dirmap["QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] =10.4305*0.15544;
   dirmap["QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8"] =2302200.0;
   dirmap["QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8"] =557600000*0.0096;
   dirmap["QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8"] =136000000*0.073;
   dirmap["QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8"] =19800000*0.146 ;
   dirmap["QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8"] =2800000*0.125;
   dirmap["QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8"] =477000*0.132;
   dirmap["QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8"] =114000*0.165;
   dirmap["QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8"] =9000*0.15;
   dirmap["QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8"] =54120000*0.002;
   dirmap["QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8"] =162060000*0.0016;
   dirmap["QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8"] =108000000*0.000225 ;
   dirmap["QCD_Pt_15to20_bcToE_TuneCUETP8M1_13TeV_pythia8"]=254596.0;
   dirmap["QCD_Pt_20to30_bcToE_TuneCUETP8M1_13TeV_pythia8"] =557627000*0.00059;
   dirmap["QCD_Pt_30to80_bcToE_TuneCUETP8M1_13TeV_pythia8"] =159068000*0.00255;
   dirmap["QCD_Pt_80to170_bcToE_TuneCUETP8M1_13TeV_pythia8"] =3221000*0.01183;
   dirmap["QCD_Pt_170to250_bcToE_TuneCUETP8M1_13TeV_pythia8"] =105771*0.02492;
   dirmap["QCD_Pt_250toInf_bcToE_TuneCUETP8M1_13TeV_pythia8"] =21094.1*0.03375;
   dirmap["GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8"] = 137751*0.001587;
   dirmap["GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8"] = 16792*0.0514;
   dirmap["TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8"] = 3.697;
   dirmap["WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"] = 405.271 ;
   dirmap["ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] = 117.864 ;

   dirmap["GluGluToZZTo2e2mu_BackgroundOnly_13TeV_MCFM"]=0.003194;
   dirmap["GluGluToZZTo2mu2tau_BackgroundOnly_13TeV_MCFM"]=0.003194;
   dirmap["GluGluToZZTo4mu_BackgroundOnly_13TeV_MCFM"]=0.001586;

   dirmap["13TeV_HN150_mumumu_width_AUTO_VmuN_0p1"]=0.0004561;
   dirmap["13TeV_HN40_mumumu_width_AUTO_VmuN_0p1"]=8.6;
   dirmap["13TeV_HN60_mumumu_width_AUTO_VmuN_0p1"]=3.416;
   dirmap["13TeV_HN700_mumumu_width_AUTO_VmuN_0p1"]=0.000002644;

   dirmap["tthwA_3l4j_hc130A30_emu_catcut"] = 0.000373;
   dirmap["tthwA_3l4j_hc130A30_3mu_catcut"] = 0.000373;
   
   dirmap["HN_pair_WR600_Zp400_HN100"] = 1.; 


   dirmap["MajoranaNeutrinoToEE_M-40_TuneZ2star_13TeV-alpgen"]=21.56;
   dirmap["MajoranaNeutrinoToEE_M-100_TuneZ2star_13TeV-alpgen"]=5.678E-02;
   dirmap["MajoranaNeutrinoToEE_M-500_TuneZ2star_13TeV-alpgen"]=7.831E-05;
   dirmap["MajoranaNeutrinoToEE_M-1500_TuneZ2star_13TeV-alpgen"]=5.375E-07;
   dirmap["MajoranaNeutrinoToMM_M-40_TuneZ2star_13TeV-alpgen"]=21.56;
   dirmap["MajoranaNeutrinoToMM_M-100_TuneZ2star_13TeV-alpgen"]=5.678E-02;
   dirmap["MajoranaNeutrinoToMM_M-500_TuneZ2star_13TeV-alpgen"]=7.831E-05;
   dirmap["MajoranaNeutrinoToMM_M-1500_TuneZ2star_13TeV-alpgen"]=5.375E-07;

   dirmap["MajoranaNeutrinoToEE_M-40_schan_lll_ss_TuneZ2star_13TeV-alpgen"]=3.593E+00;
   dirmap["MajoranaNeutrinoToEE_M-100_schan_lll_ss_TuneZ2star_13TeV-alpgen"]=9.460E-03;
   dirmap["MajoranaNeutrinoToEE_M-500_schan_lll_ss_TuneZ2star_13TeV-alpgen"]=1.305E-05;
   dirmap["MajoranaNeutrinoToEE_M-1500_schan_lll_ss_TuneZ2star_13TeV-alpgen"]=9.000E-08;
   dirmap["MajoranaNeutrinoToMuMu_M-40_schan_lll_ss_TuneZ2star_13TeV-alpgen"]=3.593E+00;
   dirmap["MajoranaNeutrinoToMuMu_M-100_schan_lll_ss_TuneZ2star_13TeV-alpgen"]=9.460E-03;
   dirmap["MajoranaNeutrinoToMuMu_M-500_schan_lll_ss_TuneZ2star_13TeV-alpgen"]=1.305E-05;
   dirmap["MajoranaNeutrinoToMuMu_M-1500_schan_lll_ss_TuneZ2star_13TeV-alpgen"]=9.000E-08;
   dirmap["MajoranaNeutrinoToEE_M-40_schan_lll_os_TuneZ2star_13TeV-alpgen"]=3.593E+00;
   dirmap["MajoranaNeutrinoToEE_M-100_schan_lll_os_TuneZ2star_13TeV-alpgen"]=9.460E-03;
   dirmap["MajoranaNeutrinoToEE_M-500_schan_lll_os_TuneZ2star_13TeV-alpgen"]=1.305E-05;
   dirmap["MajoranaNeutrinoToEE_M-1500_schan_lll_os_TuneZ2star_13TeV-alpgen"]=9.000E-08;
   dirmap["MajoranaNeutrinoToMuMu_M-40_schan_lll_os_TuneZ2star_13TeV-alpgen"]=3.593E+00;
   dirmap["MajoranaNeutrinoToMuMu_M-100_schan_lll_os_TuneZ2star_13TeV-alpgen"]=9.460E-03;
   dirmap["MajoranaNeutrinoToMuMu_M-500_schan_lll_os_TuneZ2star_13TeV-alpgen"]=1.305E-05;
   dirmap["MajoranaNeutrinoToMuMu_M-1500_schan_lll_os_TuneZ2star_13TeV-alpgen"]=9.000E-08;
   dirmap["MajoranaNeutrinoToEE_M-40_schan_ll_os_TuneZ2star_13TeV-alpgen"]=1.076E+01;
   dirmap["MajoranaNeutrinoToEE_M-100_schan_ll_os_TuneZ2star_13TeV-alpgen"]=2.8378E-02;
   dirmap["MajoranaNeutrinoToEE_M-500_schan_ll_os_TuneZ2star_13TeV-alpgen"]=3.8151E-05;
   dirmap["MajoranaNeutrinoToEE_M-1500_schan_ll_os_TuneZ2star_13TeV-alpgen"]=2.6865E-07;
   dirmap["MajoranaNeutrinoToMuMu_M-40_schan_ll_os_TuneZ2star_13TeV-alpgen"]=1.076E+01;
   dirmap["MajoranaNeutrinoToMuMu_M-100_schan_ll_os_TuneZ2star_13TeV-alpgen"]=2.8378E-02;
   dirmap["MajoranaNeutrinoToMuMu_M-500_schan_ll_os_TuneZ2star_13TeV-alpgen"]=3.8151E-05;
   dirmap["MajoranaNeutrinoToMuMu_M-1500_schan_ll_os_TuneZ2star_13TeV-alpgen"]=1.076E+01;
   
   
   // v766
   dirmap["13TeV_HN40_mumumu_width_AUTO_VmuN_0p1"] =8.605;
   dirmap["13TeV_HN60_mumumu_width_AUTO_VmuN_0p1"] =3.416;
   dirmap["13TeV_HN150_mumumu_width_AUTO_VmuN_0p1"] =0.0004561;
   dirmap["13TeV_HN700_mumumu_width_AUTO_VmuN_0p1"] =0.000002644;
   dirmap["TTToHcToWA_3l4j_3mu_Mhc130MA30_5fLO_76Xcatcut"] =0.02064;
   dirmap["TTToHcToWA_3l4j_3mu_Mhc130MA5_5fLO_76Xcatcut"] =0.02064;
   dirmap["TTToHcToWA_3l4j_3mu_Mhc155MA30_5fLO_76Xcatcut"] =0.02064;
   dirmap["TTToHcToWA_3l4j_3mu_Mhc90MA5_5fLO_76Xcatcut"] =0.02064;
   dirmap["TTToHcToWA_3l4j_emu_Mhc130MA30_5fLO_76Xcatcut"] =0.02064;
   dirmap["TTToHcToWA_3l4j_emu_Mhc130MA5_5fLO_76Xcatcut"] =0.02064;
   dirmap["TTToHcToWA_3l4j_emu_Mhc155MA30_5fLO_76Xcatcut"] =0.02064;
   dirmap["TTToHcToWA_3l4j_emu_Mhc90MA5_5fLO_76Xcatcut"] =0.02064;
   dirmap["TTToHcToWZp_3l4j_1e2mu_Mhc155MZp5_5fLO_76Xcatcut"] =0.02064;
   dirmap["TTToHcToWZp_3l4j_3mu_Mhc155MZp5_5fLO_76Xcatcut"] =0.02064;
   dirmap["MajoranaNeutrinoToEpEp_Schannel_V0p1_M-100_13TeV-amcatnlo"] =0.02282;
   dirmap["MajoranaNeutrinoToEpEp_Schannel_V0p1_M-200_13TeV-amcatnlo"] =0.001040;
   dirmap["MajoranaNeutrinoToEpEp_Schannel_V0p1_M-500_13TeV-amcatnlo"] =0.00003230;
   dirmap["MajoranaNeutrinoToEpEp_Schannel_V0p1_M-1100_13TeV-amcatnlo"] =0.00000113;
   dirmap["MajoranaNeutrinoToEmEm_Schannel_V0p1_M-100_13TeV-amcatnlo"] =0.01539;
   dirmap["MajoranaNeutrinoToEmEm_Schannel_V0p1_M-200_13TeV-amcatnlo"] =0.0006207;
   dirmap["MajoranaNeutrinoToEmEm_Schannel_V0p1_M-500_13TeV-amcatnlo"] =0.0000006553;
   dirmap["MajoranaNeutrinoToEmEm_Schannel_V0p1_M-1100_13TeV-amcatnlo"] =0.0000004028;
   dirmap["MajoranaNeutrinoToEpMup_Schannel_V0p1_M-100_13TeV-amcatnlo"] =0.02282;
   dirmap["MajoranaNeutrinoToEpMup_Schannel_V0p1_M-200_13TeV-amcatnlo"] =0.001040;
   dirmap["MajoranaNeutrinoToEpMup_Schannel_V0p1_M-500_13TeV-amcatnlo"] =0.00003230;
   dirmap["MajoranaNeutrinoToEpMup_Schannel_V0p1_M-1100_13TeV-amcatnlo"] =0.00000113;
   dirmap["MajoranaNeutrinoToEmMum_Schannel_V0p1_M-100_13TeV-amcatnlo"] =0.01539;
   dirmap["MajoranaNeutrinoToEmMum_Schannel_V0p1_M-200_13TeV-amcatnlo"] =0.0006207;
   dirmap["MajoranaNeutrinoToEmMum_Schannel_V0p1_M-500_13TeV-amcatnlo"] =0.0000006553;
   dirmap["MajoranaNeutrinoToEmMum_Schannel_V0p1_M-1100_13TeV-amcatnlo"] =0.0000004028;
   dirmap["MajoranaNeutrinoToMupMup_Schannel_V0p1_M-100_13TeV-amcatnlo"] =0.02282;
   dirmap["MajoranaNeutrinoToMupMup_Schannel_V0p1_M-200_13TeV-amcatnlo"] =0.001040;
   dirmap["MajoranaNeutrinoToMupMup_Schannel_V0p1_M-500_13TeV-amcatnlo"] =0.00003230;
   dirmap["MajoranaNeutrinoToMupMup_Schannel_V0p1_M-1100_13TeV-amcatnlo"] =0.00000113;
   dirmap["MajoranaNeutrinoToMumMum_Schannel_V0p1_M-100_13TeV-amcatnlo"] =0.01539;
   dirmap["MajoranaNeutrinoToMumMum_Schannel_V0p1_M-200_13TeV-amcatnlo"] =0.0006207;
   dirmap["MajoranaNeutrinoToMumMum_Schannel_V0p1_M-500_13TeV-amcatnlo"] =0.0000006553;
   dirmap["MajoranaNeutrinoToMumMum_Schannel_V0p1_M-1100_13TeV-amcatnlo"] =0.0000004028;
   dirmap["MajoranaNeutrinoToMupEp_Schannel_V0p1_M-100_13TeV-amcatnlo"] =0.02282;
   dirmap["MajoranaNeutrinoToMupEp_Schannel_V0p1_M-200_13TeV-amcatnlo"] =0.001040;
   dirmap["MajoranaNeutrinoToMupEp_Schannel_V0p1_M-500_13TeV-amcatnlo"] =0.00003230;
   dirmap["MajoranaNeutrinoToMupEp_Schannel_V0p1_M-1100_13TeV-amcatnlo"] =0.00000113;
   dirmap["MajoranaNeutrinoToMumEm_Schannel_V0p1_M-100_13TeV-amcatnlo"] =0.01539;
   dirmap["MajoranaNeutrinoToMumEm_Schannel_V0p1_M-200_13TeV-amcatnlo"] =0.0006207;
   dirmap["MajoranaNeutrinoToMumEm_Schannel_V0p1_M-500_13TeV-amcatnlo"] =0.0000006553;
   dirmap["MajoranaNeutrinoToMumEm_Schannel_V0p1_M-1100_13TeV-amcatnlo"] =0.0000004028;
   dirmap["MajoranaNeutrinoToEpEp_Tchannel_V0p1_M-100_13TeV-madgraph"] =0.008149;
   dirmap["MajoranaNeutrinoToEpEp_Tchannel_V0p1_M-500_13TeV-madgraph"] =0.000006136;
   dirmap["MajoranaNeutrinoToEpEp_Tchannel_V0p1_M-1100_13TeV-madgraph"] =0.0000008495;
   dirmap["MajoranaNeutrinoToEpEp_Tchannel_V0p1_M-200_13TeV-madgraph"] =0.00004982;
   dirmap["MajoranaNeutrinoToEmEm_Tchannel_V0p1_M-100_13TeV-madgraph"] =0.05385;
   dirmap["MajoranaNeutrinoToEmEm_Tchannel_V0p1_M-200_13TeV-madgraph"] =0.00002925;
   dirmap["MajoranaNeutrinoToEmEm_Tchannel_V0p1_M-500_13TeV-madgraph"] =0.000002895;
   dirmap["MajoranaNeutrinoToEmEm_Tchannel_V0p1_M-1100_13TeV-madgraph"] =0.0000003091;
   dirmap["MajoranaNeutrinoToMupMup_Tchannel_V0p1_M-100_13TeV-madgraph"] =0.008149;
   dirmap["MajoranaNeutrinoToMupMup_Tchannel_V0p1_M-200_13TeV-madgraph"] =0.00004982;
   dirmap["MajoranaNeutrinoToMupMup_Tchannel_V0p1_M-500_13TeV-madgraph"] =0.000006136;
   dirmap["MajoranaNeutrinoToMupMup_Tchannel_V0p1_M-1100_13TeV-madgraph"] =0.0000008495;
   dirmap["MajoranaNeutrinoToMumMum_Tchannel_V0p1_M-100_13TeV-madgraph"] =0.05385;
   dirmap["MajoranaNeutrinoToMumMum_Tchannel_V0p1_M-200_13TeV-madgraph"] =0.00002925;
   dirmap["MajoranaNeutrinoToMumMum_Tchannel_V0p1_M-500_13TeV-madgraph"] =0.000002895;
   dirmap["MajoranaNeutrinoToMumMum_Tchannel_V0p1_M-1100_13TeV-madgraph"] =0.0000003091;
   */
   
   return dirmap;
}

vector<TString>  GetAvailableMap(TString cversion, TString listpath, bool sig){
  vector<TString> available;

  if(cversion.Contains("v7-4"))   return available;
  if(cversion.Contains("v7-6-2")) return available;

  std::map<TString, TString> mapdir = GetLQMap2016(listpath, sig);

  TString dir = "ls   " + TString(getenv("LQANALYZER_DATASET_DIR"))+"/datasets_" + cversion + "/ > inputlist_map.txt";

  bool cluster = false;
  TString analysisdir = TString(getenv("HOSTNAME"));
  if(analysisdir.Contains("cmscluster.snu.ac.kr")) cluster=true;
  if(cluster) dir = "ls  /data4/LocalNtuples/LQAnalyzer_rootfiles_for_analysis/DataSetLists/datasets_" + cversion + "/ > inputlist_map.txt";

  system(dir.Data());
  
  std::ifstream fin("inputlist_map.txt");
  std::string word;
  vector<std::string> input_datasetlist;
  while ( fin >> word ) {
    input_datasetlist.push_back(word);
  }

  system("rm inputlist_map.txt"); 
  for(unsigned int i=0; i < input_datasetlist.size(); i++){
    std::ifstream fdin( ( string(getenv("LQANALYZER_DATASET_DIR")) +"/datasets_" + cversion + "/" + input_datasetlist.at(i)).Data());

    std::string datasetname="";
    bool missing=true;

    if(TString(input_datasetlist.at(i)).Contains("Run2015")) continue;
    if(TString(input_datasetlist.at(i)).Contains("Run2016")) continue;

    std::string dataword;
    int id=0;
    int id_check=-1;
    while ( fdin >> dataword ) {
      id++;
      if (TString(dataword) == "DataSetName") id_check = id+2;
      if(id==id_check) {
	if(TString(dataword).Contains("/")){
	  if(TString(dataword).Contains("MINIAOD")) datasetname=dataword;
	  else if(TString(dataword).Contains("USER")) datasetname=dataword;
	}
	else datasetname=dataword;
      }

      if(TString(dataword).Contains("/catTuple") ) {
	missing=false;
	if(!(TString(dataword).Contains(cversion))) {
	  cout << "dataword = "<<dataword << endl;
	  cout << "Datasets do not contain catversion in name" << endl;
	  //exit(1);
	}
      }
    }

    TString dir2 = "ls /data2/DATA/cattoflat/MC/" +  cversion +"/ > inputsnu.txt"  ;
    if(cluster) dir2 = "ls /data4/DATA/FlatCatuples/MC/" +  cversion +"/ > inputsnu.txt"  ;

    system(dir2.Data());
    
    std::ifstream fsin("inputsnu.txt");
    std::string sword;
    vector<TString> snu_files;
    while ( fsin >> sword ) {
      snu_files.push_back(TString(sword));
    }
    system("rm inputsnu.txt");
    
    bool atsnu=false;

    for(std::vector<TString>::iterator it = snu_files.begin(); it != snu_files.end();++it){
      if(TString(datasetname).Contains(*it)) atsnu=true;
    }
    if(!atsnu){
      if(!missing) available.push_back(datasetname);
    }
  }
  cout << "Samples that are missing are" << endl; 
  for(unsigned int i=0; i < available.size(); i++){
    cout << available.at(i) << endl;
  }

  return available;
}

map<TString, TString>  GetDatasetNames(TString cversion, TString listpath, bool sig){

  bool cluster = false;
  TString analysisdir = TString(getenv("HOSTNAME"));
  if(analysisdir.Contains("cmscluster.snu.ac.kr")) cluster=true;

  map<TString, TString> datasets;
  std::map<TString, TString> mapdir = GetLQMap2016(listpath, sig);
  TString dir = "ls  "+ TString(getenv("LQANALYZER_DATASET_DIR"))+"/datasets_" + cversion + "/ > inputlist_map.txt";
  if(cluster) dir = "ls  /data4/LocalNtuples/LQAnalyzer_rootfiles_for_analysis/DataSetLists/datasets_" + cversion + "/ > inputlist_map.txt";

  system(dir.Data());
  std::ifstream fin("inputlist_map.txt");
  std::string word;
  vector<std::string> input_datasetlist;
  while ( fin >> word ) {
    input_datasetlist.push_back(word);
  }
  system("rm inputlist_map.txt");
  for(unsigned int i=0; i < input_datasetlist.size(); i++){


    std::ifstream fdin( (string(getenv("LQANALYZER_DATASET_DIR"))+"/datasets_" + cversion + "/" + input_datasetlist.at(i)).Data());

    std::string datasetname="";
    if(TString(input_datasetlist.at(i)).Contains("Run2015")) continue;
    if(TString(input_datasetlist.at(i)).Contains("Run2016")) continue;
    
    std::string dataword;
    int id=0;
    int id_check=-1;
    while ( fdin >> dataword ) {
      id++;
      if (TString(dataword) == "DataSetName") id_check = id+2;
      if(id==id_check) {
	if(TString(dataword).Contains("/")){
	  if(TString(dataword).Contains("MINIAOD")) datasetname=dataword;
	  else if(TString(dataword).Contains("USER")) datasetname=dataword;
	}
	else  datasetname=dataword;
      }
      
      if(TString(dataword).Contains("/catTuple")){
      }
    }
    
    for(std::map<TString, TString>::iterator mit =mapdir.begin(); mit != mapdir.end();++mit){
      if(TString(datasetname).Contains(mit->first)) datasets[mit->second]= TString(datasetname);
    }
  }
  
  return datasets;
}
  

map<TString, TString>  GetMissingMap(TString cversion, TString listpath, bool sig){

  map<TString, TString> map_missing;

  bool cluster = false;
  TString analysisdir = TString(getenv("HOSTNAME"));
  if(analysisdir.Contains("cmscluster.snu.ac.kr")) cluster=true;
  if(cversion.Contains("v7-4"))   return map_missing;
  if(cversion.Contains("v7-6-2")) return map_missing;

  std::map<TString, TString> mapdir = GetLQMap2016(listpath, sig);
  TString dir = "ls  " + TString(getenv("LQANALYZER_DATASET_DIR"))+"/datasets_" + cversion + "/ > inputlist_map.txt";
  if(cluster) dir = "ls  /data4/LocalNtuples/LQAnalyzer_rootfiles_for_analysis/DataSetLists/datasets_" + cversion + "/ > inputlist_map.txt";

  system(dir.Data());

  std::ifstream fin("inputlist_map.txt");
  std::string word;
  vector<std::string> input_datasetlist;
  while ( fin >> word ) {
    input_datasetlist.push_back(word);
  }
  system("rm inputlist_map.txt");
  for(unsigned int i=0; i < input_datasetlist.size(); i++){
    std::ifstream fdin( (string(getenv("LQANALYZER_DATASET_DIR"))+"/datasets_" + cversion + "/" + input_datasetlist.at(i)).Data());


    std::string datasetname="";
    bool missing=true;

    if(TString(input_datasetlist.at(i)).Contains("Run2016")) continue;
    if(TString(input_datasetlist.at(i)).Contains("Run2015")) continue;

    std::string dataword;
    int id_check=-1;
    int id=0;
    while ( fdin >> dataword ) {
      id++;
      if (TString(dataword) == "DataSetName") id_check = id+2;
      if(id==id_check) {
	if(TString(dataword).Contains("/")){
	  if(TString(dataword).Contains("MINIAOD")) datasetname=dataword;
	  else if(TString(dataword).Contains("USER")) datasetname=dataword;
	}
	else datasetname=dataword;
      }
      if(TString(dataword).Contains("/catTuple")){
	missing=false;
      }
    }

    bool inlqmap=false;;

    for(std::map<TString, TString>::iterator mit =mapdir.begin(); mit != mapdir.end();++mit){
      if(TString(datasetname).Contains(mit->first)) inlqmap=true;
    }
    if(TString(input_datasetlist.at(i)).Contains("json")) inlqmap=true;    
    if(TString(input_datasetlist.at(i)).Contains("DYJets_MG")) inlqmap=true;
    if(!inlqmap) {
      cout << "LQMap (GetMissingMap) not complete" << endl;
      cout <<  input_datasetlist.at(i) << endl;
    }
    if(missing){
      for(std::map<TString, TString>::iterator mit =mapdir.begin(); mit != mapdir.end();++mit){
	if(TString(datasetname).Contains(mit->first)) map_missing[mit->second]= TString(datasetname);
      }
    }
  }   
  
  cout << "Samples that are not available at miniAOD are" << endl;
  for(std::map<TString, TString>::iterator mit =map_missing.begin(); mit !=map_missing.end();++mit){
    cout << mit->first << " "  << mit->second << endl;
  }
  return map_missing;
}


map<TString, TString>  GetTriLepMap2016(TString listpath,bool  sig){


  map<TString, TString> trimap;

  ifstream datasetname_file(listpath.Data());
  if(!datasetname_file){
    cerr << "Did not find " + listpath + ", exiting ..." << endl;
    return trimap;
  }
  while(!datasetname_file.eof()) {
    string samplealias;
    string datasetname;
    double xsec;
    string fullname;
    datasetname_file>> samplealias;
    datasetname_file >> datasetname;
    datasetname_file >> xsec;
    datasetname_file >>fullname;

    if(samplealias=="END") break;
    if(sig&&IsSignal(samplealias))    {
      if (! datasetname.empty()) trimap[TString(datasetname)] = samplealias;
    }
    if(!sig&&!IsSignal(samplealias))    {
      if (! datasetname.empty()) trimap[TString(datasetname)] = samplealias;
    }

  }

  return trimap;
}

map<TString, TString>  GetTriLepMap(){

  map<TString, TString> trimap;
  
  trimap["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] = "DY50plus_MCatNLO";
  trimap["DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] = "DY10to50_MCatNLO";
  trimap["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"] = "DY50plus_madgraph";
  trimap["DYJetsToLL_M-5to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"] = "DY5to50_madgraph";

  trimap["ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1"] ="singletop_tbar_Powheg";
  trimap["ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1"] ="singletop_t_Powheg";
  trimap["ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1"] ="singletop_tbarW_Powheg";
  trimap["ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1"] ="singletop_tW_Powheg";
  trimap["ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1"] = "singletop_s_MCatNLO";
  trimap["TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] ="TT_MCatNLO";
  trimap["TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8"] ="ttWJetsToLNu_MCatNLO";
  trimap["TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8"] ="ttWJetsToQQ_MCatNLO";
  trimap["TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] ="ttZToLLNuNu_MCatNLO";
  trimap["TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] ="ttZToQQ_MCatNLO";
  trimap["WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] ="WJets_MCatNLO";
  trimap["WW_TuneCUETP8M1_13TeV-pythia8"] ="WW_pythia8";
  trimap["WZ_TuneCUETP8M1_13TeV-pythia8"] ="WZ_pythia8";
  trimap["ZZ_TuneCUETP8M1_13TeV-pythia8"] ="ZZ_pythia8";
  trimap["ZZTo4L_13TeV_powheg_pythia8"] = "ZZ_llll_powheg";
  trimap["ZZTo4L_13TeV-amcatnloFXFX-pythia8"] = "ZZ_llll_MCatNLO";
  trimap["ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8"] = "ZZ_llqq_MCatNLO";
  trimap["ZZTo2L2Nu_13TeV_powheg_pythia8"] = "ZZ_llnunu_powheg";
  trimap["WWTo2L2Nu_13TeV-powheg"] = "WW_llnn_powheg";
  trimap["WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8"] = "WZ_llqq_MCatNLO";
  trimap["WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8"] = "WZ_lllnu_powheg";
  trimap["WZJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] = "WZ_lllnu_MCatNLO";
  trimap["WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8"] = "WpWp_madgraph";
  trimap["WpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8"] = "WpWp_qcd_madgraph";
  trimap["WW_DoubleScattering_13TeV-pythia8"] = "WW_doublescattering";
  trimap["WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] = "WZZ_MCatNLO";
  trimap["ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] = "ZZZ_MCatNLO";
  trimap["WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] = "WWZ_MCatNLO";

  return trimap;
}
  

map<TString, TString>  GetLQMap2016(TString datasetfile, bool sig){

  map<TString, TString> lqmap;

  ifstream datasetname_file(datasetfile.Data());
  if(!datasetname_file){
    cerr << "Did not find " + datasetfile + ", exiting ..." << endl;
    return lqmap;
  }
  while(!datasetname_file.eof()) {
    string samplealias;
    string datasetname;
    double xsec;
    string fullname;
    datasetname_file>> samplealias;
    datasetname_file >> datasetname;
    datasetname_file >> xsec;
    datasetname_file >>fullname;

    if(samplealias=="END") break;
    if(sig&&IsSignal(samplealias))    {
      if (! datasetname.empty()) lqmap[TString(datasetname)] = samplealias;
    }
    if(!sig&&!IsSignal(samplealias))    {
      if (! datasetname.empty()) lqmap[TString(datasetname)] = samplealias;
    }

  }


  /*
  lqmap["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] = "DY50plus_MCatNLO";
  lqmap["DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] = "DY10to50_MCatNLO";
  lqmap["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"] = "DY50plus_madgraph";
  lqmap["DYJetsToLL_M-5to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"] = "DY5to50_madgraph";
  lqmap["DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]="DY10to50_madgraph";

  lqmap["ttbb_4FS_ckm_amcatnlo_madspin_pythia8"]="ttbb_MCatNLO";

  lqmap["ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1"] ="singletop_tbar_Powheg";
  lqmap["ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1"] ="singletop_t_Powheg";
  lqmap["ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1"] ="singletop_tbarW_Powheg";
  lqmap["ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1"] ="singletop_tW_Powheg";
  lqmap["ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1"] = "singletop_s_MCatNLO";
  lqmap["ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1"]="singletop_tW_nofullhad_Powheg";
  lqmap["ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1"]="singletop_tbarW_nofullhad_Powheg";
  lqmap["ttHTobb_M125_13TeV_powheg_pythia8"] ="ttHtobb_Powheg";
  lqmap["ttHToNonbb_M125_13TeV_powheg_pythia8"] ="ttHnobb_Powheg";

  lqmap["GluGlu_HToMuMu_M125_13TeV_powheg_pythia8"] ="ggHtomm_Powheg";
  lqmap["VBF_HToMuMu_M125_13TeV_powheg_pythia8"] ="vhf_Htomm_Powheg";

  lqmap["WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] ="WJets_MCatNLO";
  lqmap["WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]="WJets_madgraph";
  lqmap["WW_TuneCUETP8M1_13TeV-pythia8"] ="WW_pythia8";
  lqmap["WZ_TuneCUETP8M1_13TeV-pythia8"] ="WZ_pythia8";
  lqmap["ZZ_TuneCUETP8M1_13TeV-pythia8"] ="ZZ_pythia8";

  lqmap["WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] = "WZZ_MCatNLO";
  lqmap["ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] = "ZZZ_MCatNLO";
  lqmap["WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] = "WWZ_MCatNLO";
  lqmap["WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8"]="WWW_MCatNLO";;

  lqmap["TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] ="TT_MCatNLO";
  lqmap["TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]="TT_madgraph";
  lqmap["TT_TuneCUETP8M1_13TeV-powheg-pythia8"]="TT_powheg";
  lqmap["TT_TuneCUETP8M1_mtop1695_13TeV-powheg-pythia8"]="TT_powheg_m169";
  lqmap["TT_TuneCUETP8M1_mtop1755_13TeV-powheg-pythia8"]="TT_powheg_m171";
  lqmap["TTTo2L2Nu_13TeV-powheg"]="TT_2l2nu_powheg";
  lqmap["TTTo2L2Nu_13TeV_ScaleDown-powheg"]="TT_2l2nu_Sup_powheg";
  lqmap["TTTo2L2Nu_13TeV_ScaleUp-powheg"]="TT_2l2nu_Sdownpowheg";
  lqmap["TTToSemiLeptonic_13TeV-powheg"]="TT_semilep_powheg";
  lqmap["TTToSemiLeptonic_13TeV_ScaleDown-powheg"]="TT_semilep_Sup_powheg";
  lqmap["TTToSemiLeptonic_13TeV_ScaleUp-powheg"]="TT_semilep_Sdown_powheg";


  lqmap["ttWJets_13TeV_madgraphMLM"]="ttWJets";
  lqmap["ttZJets_13TeV_madgraphMLM"]="ttZJets";

  
  lqmap["13TeV_HN150_mumumu_width_AUTO_VmuN_0p1"]="HN150_mumumu_VmuN_0p1";
  lqmap["13TeV_HN40_mumumu_width_AUTO_VmuN_0p1"]="HN40_mumumu_VmuN_0p1";
  lqmap["13TeV_HN60_mumumu_width_AUTO_VmuN_0p1"]="HN60_mumumu_VmuN_0p1";
  lqmap["13TeV_HN700_mumumu_width_AUTO_VmuN_0p1"]="HN700_mumumu_VmuN_0p1";

  */
  return lqmap;
}
map<TString, TString>  GetLQMap(TString datasetfile){

  map<TString, TString> lqmap;

  ifstream datasetname_file(datasetfile.Data());
  if(!datasetname_file){
    cerr << "Did not find " + datasetfile + ", exiting ..." << endl;
    return lqmap;
  }
  while(!datasetname_file.eof()) {
    string samplealias;
    string datasetname;
    double xsec;
    string fullname;
    datasetname_file>> samplealias;
    datasetname_file >> datasetname;
    datasetname_file >> xsec;
    datasetname_file >> fullname;
    if(samplealias=="END") break;
    if (! datasetname.empty())  lqmap[TString(datasetname)] = samplealias;
  }
  /*
  lqmap["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] = "DY50plus_MCatNLO";
  lqmap["DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] = "DY10to50_MCatNLO";
  lqmap["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"] = "DY50plus_madgraph";
  lqmap["DYJetsToLL_M-5to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"] = "DY5to50_madgraph";

  lqmap["ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1"] ="singletop_tbar_Powheg";
  lqmap["ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1"] ="singletop_t_Powheg";
  lqmap["ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1"] ="singletop_tbarW_Powheg";
  lqmap["ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1"] ="singletop_tW_Powheg";
  lqmap["ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1"] = "singletop_s_MCatNLO";
  lqmap["ttHTobb_M125_13TeV_powheg_pythia8"] ="ttHtobb_Powheg";
  lqmap["ttHToNonbb_M125_13TeV_powheg_pythia8"] ="ttHnobb_Powheg";

  lqmap["GluGlu_HToMuMu_M125_13TeV_powheg_pythia8"] ="ggHtomm_Powheg";
  lqmap["VBF_HToMuMu_M125_13TeV_powheg_pythia8"] ="vhf_Htomm_Powheg";

  lqmap["TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] ="TT_MCatNLO";
  lqmap["TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-scaleup-pythia8"] = "TT_scaleup_MCatNLO";
  lqmap["TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-scaledown-pythia8"] = "TT_scaledown_MCatNLO";
  lqmap["TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"] ="TT_MG5";
  lqmap["TT_TuneCUETP8M1_13TeV-powheg-scaleup-pythia8"] = "TT_scaleup_powheg";
  lqmap["TT_TuneCUETP8M1_13TeV-powheg-scaledown-pythia8"]= "TT_scaledown_powheg";
  lqmap["TT_TuneCUETP8M1_13TeV-powheg-pythia8"]=  "TT_powheg";
  lqmap["TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8"] ="ttWJetsToLNu_MCatNLO";
  lqmap["TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8"] ="ttWJetsToQQ_MCatNLO";
  lqmap["TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] ="ttZToLLNuNu_MCatNLO";
  lqmap["TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] ="ttZToQQ_MCatNLO";
  lqmap["TT_TuneZ2star_13TeV-powheg-pythia6-tauola"] = "TT_powheg_pythia6";

  lqmap["TT_TuneEE5C_13TeV-amcatnlo-herwigpp"]= "TT_amcatnlo-herwigpp";
  lqmap["TT_TuneEE5C_13TeV-powheg-herwigpp"]="TT_powheg-herwigpp"; 
  lqmap["TT_TuneCUETP8M1_mtop1665_13TeV-powheg-pythia8"]="TT_powheg_mtop1665";
  lqmap["TT_TuneCUETP8M1_mtop1695_13TeV-powheg-pythia8"]="TT_powheg_mtop1695";
  lqmap["TT_TuneCUETP8M1_mtop1715_13TeV-powheg-pythia8"]="TT_powheg_mtop1715";
  lqmap["TT_TuneCUETP8M1_mtop1735_13TeV-powheg-pythia8"]="TT_powheg_mtop1735";
  lqmap["TT_TuneCUETP8M1_mtop1755_13TeV-powheg-pythia8"]="TT_powheg_mtop1755";
  lqmap["TT_TuneCUETP8M1_mtop1785_13TeV-powheg-pythia8"]="TT_powheg_mtop1785"; 
  lqmap["TTbarDMJets_pseudoscalar_Mchi-10_Mphi-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]="TTbarDMJets_pseudoscalar_Mchi-10_Mphi-100";
  lqmap["TTbarDMJets_pseudoscalar_Mchi-10_Mphi-10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"]="TTbarDMJets_pseudoscalar_Mchi-10_Mphi-10";
  lqmap["WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] ="WJets_MCatNLO";
  lqmap["WW_TuneCUETP8M1_13TeV-pythia8"] ="WW_pythia8";
  lqmap["WZ_TuneCUETP8M1_13TeV-pythia8"] ="WZ_pythia8";
  lqmap["ZZ_TuneCUETP8M1_13TeV-pythia8"] ="ZZ_pythia8";
  
  lqmap["QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] = "QCD_mu15to20_pythia8";
  lqmap["QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] ="QCD_mu20to30_pythia8";
  lqmap["QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] ="QCD_mu30to50_pythia8";
  lqmap["QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] ="QCD_mu50to80_pythia8";
  lqmap["QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] ="QCD_mu80to120_pythia8";
  lqmap["QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] ="QCD_mu120to170_pythia8";
  lqmap["QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] ="QCD_mu170to300_pythia8";
  lqmap["QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] ="QCD_mu300to470_pythia8";
  lqmap["QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] ="QCD_mu470to600_pythia8";
  lqmap["QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] ="QCD_mu600to800_pythia8";
  lqmap["QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"] ="QCD_mu800to1000_pythia8";
  lqmap["QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8"]= "QCD_mu1000toINF_pythia8";
  lqmap["QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8"] ="QCD_em15to20_pythia8";
  lqmap["QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8"] ="QCD_em20to30_pythia8";
  lqmap["QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8"] ="QCD_em30to50_pythia8";
  lqmap["QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8"] ="QCD_em50to80_pythia8";
  lqmap["QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8"] ="QCD_em80to120_pythia8";
  lqmap["QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8"] ="QCD_em120to170_pythia8";
  lqmap["QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8"] ="QCD_em170to300_pythia8";
  lqmap["QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8"] ="QCD_em300toINF_pythia8";
  lqmap["QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8"] = "QCD_DoubleEM_40toInf_pythia8";
  lqmap["QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8"] ="QCD_DoubleEM_30toInf_pythia8";
  lqmap["QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8"] ="QCD_DoubleEM_30to40_pythia8";
  lqmap["QCD_Pt_15to20_bcToE_TuneCUETP8M1_13TeV_pythia8"]="QCD_15to20_bcToE_pythia8";
  lqmap["QCD_Pt_20to30_bcToE_TuneCUETP8M1_13TeV_pythia8"] ="QCD_20to30_bcToE_pythia8";
  lqmap["QCD_Pt_80to170_bcToE_TuneCUETP8M1_13TeV_pythia8"] ="QCD_80to170_bcToE_pythia8";
  lqmap["QCD_Pt_170to250_bcToE_TuneCUETP8M1_13TeV_pythia8"] ="QCD_170to250_bcToE_pythia8";
  lqmap["QCD_Pt_250toInf_bcToE_TuneCUETP8M1_13TeV_pythia8"] ="QCD_250toInf_bcToE_pythia8";
  lqmap["QCD_Pt_30to80_bcToE_TuneCUETP8M1_13TeV_pythia8"] = "QCD_30to80_bcToE_pythia8";

  lqmap["ZZTo4L_13TeV_powheg_pythia8"] = "ZZ_llll_powheg";
  lqmap["ZZTo4L_13TeV-amcatnloFXFX-pythia8"] = "ZZ_llll_MCatNLO";
  lqmap["ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8"] = "ZZ_llqq_MCatNLO";
  lqmap["ZZTo2L2Nu_13TeV_powheg_pythia8"] = "ZZ_llnunu_powheg";
  lqmap["WWTo2L2Nu_13TeV-powheg"] = "WW_llnn_powheg";
  lqmap["WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8"] = "WZ_llqq_MCatNLO";
  lqmap["WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8"] = "WZ_lllnu_powheg";
  lqmap["WZJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] = "WZ_lllnu_MCatNLO";
  lqmap["WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8"] = "WpWp_madgraph";
  lqmap["WpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8"] = "WpWp_qcd_madgraph";

  lqmap["WW_DoubleScattering_13TeV-pythia8"] = "WW_doublescattering";
  lqmap["WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] = "WZZ_MCatNLO";
  lqmap["ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] = "ZZZ_MCatNLO";
  lqmap["WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8"] = "WWZ_MCatNLO";
  lqmap["GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8"] = "GJet_20to40_pythia8";
  lqmap["GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8"] = "GJet_40plus_pythia8";
  lqmap["TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8"] = "TTG_MCatNLO";
  lqmap["WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"] = "WG_lnuG_madgraph";
  lqmap["ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"] = "ZG_llG_MCatNLO";
  lqmap["GluGluToZZTo2e2mu_BackgroundOnly_13TeV_MCFM"]="GluGluToZZTo2e2mu";
  lqmap["GluGluToZZTo2mu2tau_BackgroundOnly_13TeV_MCFM"]="GluGluToZZTo2mu2tau";
  lqmap["GluGluToZZTo4mu_BackgroundOnly_13TeV_MCFM"]="GluGluToZZTo4mu";

 
  lqmap["13TeV_HN150_mumumu_width_AUTO_VmuN_0p1"]="HN150_mumumu_VmuN_0p1";
  lqmap["13TeV_HN40_mumumu_width_AUTO_VmuN_0p1"]="HN40_mumumu_VmuN_0p1";
  lqmap["13TeV_HN60_mumumu_width_AUTO_VmuN_0p1"]="HN60_mumumu_VmuN_0p1";
  lqmap["13TeV_HN700_mumumu_width_AUTO_VmuN_0p1"]="HN700_mumumu_VmuN_0p1";
  lqmap["tthwA_3l4j_hc130A30_emu_catcut"]="tthwA_3l4j_hc130A30_emu";
  lqmap["tthwA_3l4j_hc130A30_3mu_catcut"]="tthwA_3l4j_hc130A30_3mu";
  lqmap["MajoranaNeutrinoToEE_M-40_TuneZ2star_13TeV-alpgen"]="HN_EE_M40";
  lqmap["MajoranaNeutrinoToEE_M-100_TuneZ2star_13TeV-alpgen"]="HN_EE_M100";
  lqmap["MajoranaNeutrinoToEE_M-500_TuneZ2star_13TeV-alpgen"]="HN_EE_M500";
  lqmap["MajoranaNeutrinoToEE_M-1500_TuneZ2star_13TeV-alpgen"]="HN_EE_M1500";
  lqmap["MajoranaNeutrinoToMM_M-40_TuneZ2star_13TeV-alpgen"]="HN_MM_M40";
  lqmap["MajoranaNeutrinoToMM_M-100_TuneZ2star_13TeV-alpgen"]="HN_MM_M100";
  lqmap["MajoranaNeutrinoToMM_M-500_TuneZ2star_13TeV-alpgen"]="HN_MM_M500";
  lqmap["MajoranaNeutrinoToMM_M-1500_TuneZ2star_13TeV-alpgen"]="HN_MM_M1500";

  lqmap["HN_pair_WR600_Zp400_HN100"] = "HN_pair_WR600_Zp400_HN100";

  lqmap["MajoranaNeutrinoToEE_M-40_schan_lll_ss_TuneZ2star_13TeV-alpgen"]="HN_schan_lll_ee_ss_40";
  lqmap["MajoranaNeutrinoToEE_M-100_schan_lll_ss_TuneZ2star_13TeV-alpgen"]="HN_schan_lll_ee_ss_100";
  lqmap["MajoranaNeutrinoToEE_M-500_schan_lll_ss_TuneZ2star_13TeV-alpgen"]="HN_schan_lll_ee_ss_500";
  lqmap["MajoranaNeutrinoToEE_M-1500_schan_lll_ss_TuneZ2star_13TeV-alpgen"]="HN_schan_lll_ee_ss_1500";
  lqmap["MajoranaNeutrinoToMuMu_M-40_schan_lll_ss_TuneZ2star_13TeV-alpgen"]="HN_schan_lll_mm_ss_40";
  lqmap["MajoranaNeutrinoToMuMu_M-100_schan_lll_ss_TuneZ2star_13TeV-alpgen"]="HN_schan_lll_mm_ss_100";
  lqmap["MajoranaNeutrinoToMuMu_M-500_schan_lll_ss_TuneZ2star_13TeV-alpgen"]="HN_schan_lll_mm_ss_500";
  lqmap["MajoranaNeutrinoToMuMu_M-1500_schan_lll_ss_TuneZ2star_13TeV-alpgen"]="HN_schan_lll_mm_ss_1500";
  lqmap["MajoranaNeutrinoToEE_M-40_schan_lll_os_TuneZ2star_13TeV-alpgen"]="HN_schan_lll_ee_os_40";
  lqmap["MajoranaNeutrinoToEE_M-100_schan_lll_os_TuneZ2star_13TeV-alpgen"]="HN_schan_lll_ee_os_100";
  lqmap["MajoranaNeutrinoToEE_M-500_schan_lll_os_TuneZ2star_13TeV-alpgen"]="HN_schan_lll_ee_os_500";
  lqmap["MajoranaNeutrinoToEE_M-1500_schan_lll_os_TuneZ2star_13TeV-alpgen"]="HN_schan_lll_ee_os_1500";
  lqmap["MajoranaNeutrinoToMuMu_M-40_schan_lll_os_TuneZ2star_13TeV-alpgen"]="HN_schan_lll_mm_os_40";
  lqmap["MajoranaNeutrinoToMuMu_M-100_schan_lll_os_TuneZ2star_13TeV-alpgen"]="HN_schan_lll_mm_os_100";
  lqmap["MajoranaNeutrinoToMuMu_M-500_schan_lll_os_TuneZ2star_13TeV-alpgen"]="HN_schan_lll_mm_os_500";
  lqmap["MajoranaNeutrinoToMuMu_M-1500_schan_lll_os_TuneZ2star_13TeV-alpgen"]="HN_schan_lll_mm_os_1500";
  lqmap["MajoranaNeutrinoToEE_M-40_schan_ll_os_TuneZ2star_13TeV-alpgen"]="HN_schan_ll_ee_os_40";
  lqmap["MajoranaNeutrinoToEE_M-100_schan_ll_os_TuneZ2star_13TeV-alpgen"]="HN_schan_ll_ee_os_100";
  lqmap["MajoranaNeutrinoToEE_M-500_schan_ll_os_TuneZ2star_13TeV-alpgen"]="HN_schan_ll_ee_os_500";
  lqmap["MajoranaNeutrinoToEE_M-1500_schan_ll_os_TuneZ2star_13TeV-alpgen"]="HN_schan_ll_ee_os_1500";
  lqmap["MajoranaNeutrinoToMuMu_M-40_schan_ll_os_TuneZ2star_13TeV-alpgen"]="HN_schan_ll_mm_os_40";
  lqmap["MajoranaNeutrinoToMuMu_M-100_schan_ll_os_TuneZ2star_13TeV-alpgen"]="HN_schan_ll_mm_os_100";
  lqmap["MajoranaNeutrinoToMuMu_M-500_schan_ll_os_TuneZ2star_13TeV-alpgen"]="HN_schan_ll_mm_os_500";
  lqmap["MajoranaNeutrinoToMuMu_M-1500_schan_ll_os_TuneZ2star_13TeV-alpgen"]="HN_schan_ll_mm_os_1500";

  
  lqmap["13TeV_HN40_mumumu_width_AUTO_VmuN_0p1"] ="HN40_mumumu";
  lqmap["13TeV_HN60_mumumu_width_AUTO_VmuN_0p1"] ="HN60_mumumu";
  lqmap["13TeV_HN150_mumumu_width_AUTO_VmuN_0p1"] ="HN150_mumumu";
  lqmap["13TeV_HN700_mumumu_width_AUTO_VmuN_0p1"] ="HN700_mumumu";
  lqmap["TTToHcToWA_3l4j_3mu_Mhc130MA30_5fLO_76Xcatcut"] ="TTToHcToWA_3l4j_3mu_Mhc130MA3";
  lqmap["TTToHcToWA_3l4j_3mu_Mhc130MA5_5fLO_76Xcatcut"] ="TTToHcToWA_3l4j_3mu_Mhc130MA5";
  lqmap["TTToHcToWA_3l4j_3mu_Mhc155MA30_5fLO_76Xcatcut"] ="TTToHcToWA_3l4j_3mu_Mhc155MA30";
  lqmap["TTToHcToWA_3l4j_3mu_Mhc90MA5_5fLO_76Xcatcut"] ="TTToHcToWA_3l4j_3mu_Mhc90MA5";
  lqmap["TTToHcToWA_3l4j_emu_Mhc130MA30_5fLO_76Xcatcut"] ="TTToHcToWA_3l4j_emu_Mhc130MA30";
  lqmap["TTToHcToWA_3l4j_emu_Mhc130MA5_5fLO_76Xcatcut"] ="TTToHcToWA_3l4j_emu_Mhc130MA5";
  lqmap["TTToHcToWA_3l4j_emu_Mhc155MA30_5fLO_76Xcatcut"] ="TTToHcToWA_3l4j_emu_Mhc155MA30";
  lqmap["TTToHcToWA_3l4j_emu_Mhc90MA5_5fLO_76Xcatcut"] ="TTToHcToWA_3l4j_emu_Mhc90MA5";
  lqmap["TTToHcToWZp_3l4j_1e2mu_Mhc155MZp5_5fLO_76Xcatcut"] ="TTToHcToWZp_3l4j_1e2mu_Mhc155MZp5";
  lqmap["TTToHcToWZp_3l4j_3mu_Mhc155MZp5_5fLO_76Xcatcut"] ="TTToHcToWZp_3l4j_3mu_Mhc155MZp5";
  lqmap["MajoranaNeutrinoToEpEp_Schannel_V0p1_M-100_13TeV-amcatnlo"] ="HN_EpEp_Schannel_100";
  lqmap["MajoranaNeutrinoToEpEp_Schannel_V0p1_M-200_13TeV-amcatnlo"] ="HN_EpEp_Schannel_200";
  lqmap["MajoranaNeutrinoToEpEp_Schannel_V0p1_M-500_13TeV-amcatnlo"] ="HN_EpEp_Schannel_500";
  lqmap["MajoranaNeutrinoToEpEp_Schannel_V0p1_M-1100_13TeV-amcatnlo"] ="HN_EpEp_Schannel_1100";
  lqmap["MajoranaNeutrinoToEmEm_Schannel_V0p1_M-100_13TeV-amcatnlo"] ="HN_EmEm_Schannel_100";
  lqmap["MajoranaNeutrinoToEmEm_Schannel_V0p1_M-200_13TeV-amcatnlo"] ="HN_EmEm_Schannel_200";
  lqmap["MajoranaNeutrinoToEmEm_Schannel_V0p1_M-500_13TeV-amcatnlo"] ="HN_EmEm_Schannel_500";
  lqmap["MajoranaNeutrinoToEmEm_Schannel_V0p1_M-1100_13TeV-amcatnlo"] ="HN_EmEm_Schannel_1100";
  lqmap["MajoranaNeutrinoToEpMup_Schannel_V0p1_M-100_13TeV-amcatnlo"] ="HN_EpMup_Schannel_100";
  lqmap["MajoranaNeutrinoToEpMup_Schannel_V0p1_M-200_13TeV-amcatnlo"] ="HN_EpMup_Schannel_200";
  lqmap["MajoranaNeutrinoToEpMup_Schannel_V0p1_M-500_13TeV-amcatnlo"] ="HN_EpMup_Schannel_500";
  lqmap["MajoranaNeutrinoToEpMup_Schannel_V0p1_M-1100_13TeV-amcatnlo"] ="HN_EpMup_Schannel_1100";
  lqmap["MajoranaNeutrinoToEmMum_Schannel_V0p1_M-100_13TeV-amcatnlo"] ="HN_EmMum_Schannel_100";
  lqmap["MajoranaNeutrinoToEmMum_Schannel_V0p1_M-200_13TeV-amcatnlo"] ="HN_EmMum_Schannel_200";
  lqmap["MajoranaNeutrinoToEmMum_Schannel_V0p1_M-500_13TeV-amcatnlo"] ="HN_EmMum_Schannel_500";
  lqmap["MajoranaNeutrinoToEmMum_Schannel_V0p1_M-1100_13TeV-amcatnlo"] ="HN_EmMum_Schannel_1100";
  lqmap["MajoranaNeutrinoToMupMup_Schannel_V0p1_M-100_13TeV-amcatnlo"] ="HN_EpMup_Schannel_100";
  lqmap["MajoranaNeutrinoToMupMup_Schannel_V0p1_M-200_13TeV-amcatnlo"] ="HN_EpMup_Schannel_200";
  lqmap["MajoranaNeutrinoToMupMup_Schannel_V0p1_M-500_13TeV-amcatnlo"] ="HN_EpMup_Schannel_500";
  lqmap["MajoranaNeutrinoToMupMup_Schannel_V0p1_M-1100_13TeV-amcatnlo"] ="HN_EpMup_Schannel_1100";
  lqmap["MajoranaNeutrinoToMumMum_Schannel_V0p1_M-100_13TeV-amcatnlo"] ="HN_MumMum_Schannel_100";
  lqmap["MajoranaNeutrinoToMumMum_Schannel_V0p1_M-200_13TeV-amcatnlo"] ="HN_MumMum_Schannel_200";
  lqmap["MajoranaNeutrinoToMumMum_Schannel_V0p1_M-500_13TeV-amcatnlo"] ="HN_MumMum_Schannel_500";
  lqmap["MajoranaNeutrinoToMumMum_Schannel_V0p1_M-1100_13TeV-amcatnlo"] ="HN_MumMum_Schannel_1100";
  lqmap["MajoranaNeutrinoToMupEp_Schannel_V0p1_M-100_13TeV-amcatnlo"] ="HN_MupMup_Schannel_100";
  lqmap["MajoranaNeutrinoToMupEp_Schannel_V0p1_M-200_13TeV-amcatnlo"] ="HN_MupMup_Schannel_200";
  lqmap["MajoranaNeutrinoToMupEp_Schannel_V0p1_M-500_13TeV-amcatnlo"] ="HN_MupMup_Schannel_500";
  lqmap["MajoranaNeutrinoToMupEp_Schannel_V0p1_M-1100_13TeV-amcatnlo"] ="HN_MupMup_Schannel_1100";
  lqmap["MajoranaNeutrinoToMumEm_Schannel_V0p1_M-100_13TeV-amcatnlo"] ="HN_MumEm_Schannel_100";
  lqmap["MajoranaNeutrinoToMumEm_Schannel_V0p1_M-200_13TeV-amcatnlo"] ="HN_MumEm_Schannel_200";
  lqmap["MajoranaNeutrinoToMumEm_Schannel_V0p1_M-500_13TeV-amcatnlo"] ="HN_MumEm_Schannel_500";
  lqmap["MajoranaNeutrinoToMumEm_Schannel_V0p1_M-1100_13TeV-amcatnlo"] ="HN_MumEm_Schannel_1100";
  lqmap["MajoranaNeutrinoToEpEp_Tchannel_V0p1_M-100_13TeV-madgraph"] ="HN_EpEp_Tchannel_100";
  lqmap["MajoranaNeutrinoToEpEp_Tchannel_V0p1_M-500_13TeV-madgraph"] ="HN_EpEp_Tchannel_500";
  lqmap["MajoranaNeutrinoToEpEp_Tchannel_V0p1_M-1100_13TeV-madgraph"] ="HN_EpEp_Tchannel_1100";
  lqmap["MajoranaNeutrinoToEpEp_Tchannel_V0p1_M-200_13TeV-madgraph"] ="HN_EpEp_Tchannel_200";
  lqmap["MajoranaNeutrinoToEmEm_Tchannel_V0p1_M-100_13TeV-madgraph"] ="HN_EmEm_Tchannel_100";
  lqmap["MajoranaNeutrinoToEmEm_Tchannel_V0p1_M-200_13TeV-madgraph"] ="HN_EmEm_Tchannel_200";
  lqmap["MajoranaNeutrinoToEmEm_Tchannel_V0p1_M-500_13TeV-madgraph"] ="HN_EmEm_Tchannel_500";
  lqmap["MajoranaNeutrinoToEmEm_Tchannel_V0p1_M-1100_13TeV-madgraph"] ="HN_EmEm_Tchannel_1100";
  lqmap["MajoranaNeutrinoToMupMup_Tchannel_V0p1_M-100_13TeV-madgraph"] ="HN_MupMup_Tchannel_100";
  lqmap["MajoranaNeutrinoToMupMup_Tchannel_V0p1_M-200_13TeV-madgraph"] ="HN_MupMup_Tchannel_200";
  lqmap["MajoranaNeutrinoToMupMup_Tchannel_V0p1_M-500_13TeV-madgraph"] ="HN_MupMup_Tchannel_500";
  lqmap["MajoranaNeutrinoToMupMup_Tchannel_V0p1_M-1100_13TeV-madgraph"] ="HN_MupMup_Tchannel_1100";
  lqmap["MajoranaNeutrinoToMumMum_Tchannel_V0p1_M-100_13TeV-madgraph"] ="HN_MumMum_Tchannel_100";
  lqmap["MajoranaNeutrinoToMumMum_Tchannel_V0p1_M-200_13TeV-madgraph"] ="HN_MumMum_Tchannel_200";
  lqmap["MajoranaNeutrinoToMumMum_Tchannel_V0p1_M-500_13TeV-madgraph"] ="HN_MumMum_Tchannel_500";
  lqmap["MajoranaNeutrinoToMumMum_Tchannel_V0p1_M-1100_13TeV-madgraph"] ="HN_MumMum_Tchannel_1100";
  */
  return lqmap;
}

#endif
