/***************************************************************************
 * @Project: LQFakeRateCalculator_El Frame - ROOT-based analysis framework for Korea SNU
OB * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/


/// Local includes
#include "FakeRateCalculator_El.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (FakeRateCalculator_El);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
FakeRateCalculator_El::FakeRateCalculator_El() :  AnalyzerCore(),  out_electrons(0) {


  // To have the correct name in the log:                                                                                                                            
  SetLogName("FakeRateCalculator_El");

  Message("In FakeRateCalculator_El constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void FakeRateCalculator_El::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
  
  MakeCleverHistograms(sighist_ee, "SingleLooseElJet");
  MakeCleverHistograms(sighist_ee, "SingleTightElJet");

  MakeCleverHistograms(sighist_ee, "DiElectron");
  MakeCleverHistograms(sighist_ee, "DiElectron_SingleLeg");
  MakeCleverHistograms(sighist_ee, "Use40_loose");
  MakeCleverHistograms(sighist_ee, "Use40_tight");

  MakeCleverHistograms(sighist_ee, "SingleElectron_unprescaled");
  return;
}


void FakeRateCalculator_El::ExecuteEvents()throw( LQError ){
    
  //// Initial event cuts
  /// MET FIleters 
  if(!PassMETFilter()) return;     
  
  /// Require good promary vertex 
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex  
  numberVertices = eventbase->GetEvent().nVertices();   

  /// These run on double electron dataset
  if((isData&&k_channel == "DoubleEG") || !isData){

    // GetFakeRateAndPromptRates fills fake rate + prompt rate plots
    // hist names are tagged with second argu
    // dijet method

    
    //// SIGNAL  ID 
    /// ELECTRON16_HN_FAKELOOSE_NOD0 = 0.5 reliso cut
    ///                                dxy / dxysig = N/A
    /// ELECTRON16_HN_TIGHT            iso = Tight iso of POG (checking opt.)
    ///                                dxy cut = safe cuts from https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
    GetFakeRateAndPromptRates("ELECTRON16_HN_FAKELOOSE_NOD0","dijet_nod0",  "ELECTRON16_HN_TIGHT",weight,true,  true);

    ///  ELECTRON16_HN_TIGHT_DXYSIG  = ELECTRON16_HN_TIGHT + dxysig < 3.
    GetFakeRateAndPromptRates("ELECTRON16_HN_FAKELOOSE_NOD0","dijet_nod0_dxysig",  "ELECTRON16_HN_TIGHT_DXYSIG",weight,true,  true);
    GetFakeRateAndPromptRates("ELECTRON16_HN_FAKELOOSE_NOD0","dijet_nod0_dxysig_miniiso",  "ELECTRON16_HN_TIGHT_DXYSIG",weight,true,  true);

    /// Apply D0 cut (in z peak)
    /// ELECTRON16_HN_FAKELOOSE = 0.5 reliso cut   ;  dxy = safe from https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
    GetFakeRateAndPromptRates("ELECTRON16_HN_FAKELOOSE_D0","dijet_d0", "ELECTRON16_HN_TIGHT",weight,true, true);
    GetFakeRateAndPromptRates("ELECTRON16_HN_FAKELOOSE_D0_DXYSIG","dijet_d0_dxysig", "ELECTRON16_HN_TIGHT_DXYSIG",weight,true, true);
    GetFakeRateAndPromptRates("ELECTRON16_HN_FAKELOOSE_D0_DXYSIG","dijet_d0_dxysig_miniiso",  "ELECTRON16_HN_TIGHT_DXYSIG",weight,true,  true);
    
    


    /// ISO SYSTEMATICS
    GetFakeRateAndPromptRates("ELECTRON16_HN_FAKELOOSE_ISO04","dijet_iso04", "ELECTRON16_HN_TIGHT",weight,true, false);
    GetFakeRateAndPromptRates("ELECTRON16_HN_FAKELOOSE_ISO06", "dijet_iso06","ELECTRON16_HN_TIGHT",weight, true,false);
   
    /// POG FAKERATES  iso and d0 cuts are same for tight and medium, only shower id is changed
    GetFakeRateAndPromptRates("ELECTRON16_POG_FAKELOOSE","dijet_pog", "ELECTRON16_FR_POG_TIGHT",weight, true,false);
    GetFakeRateAndPromptRates("ELECTRON16_POG_MEDIUM_FAKELOOSE","dijet_pog", "ELECTRON16_FR_POG_MEDIUM",weight, true,false);
    GetFakeRateAndPromptRates("ELECTRON16_POG_FAKELOOSE_CC","dijet_pog", "ELECTRON16_FR_POG_TIGHT_CC",weight, true,false);
    GetFakeRateAndPromptRates("ELECTRON16_POG_MEDIUM_FAKELOOSE_CC","dijet_pog", "ELECTRON16_FR_POG_MEDIUM_CC",weight, true,false);
    GetFakeRateAndPromptRates("ELECTRON16_MVA_FAKELOOSE_CC","dijet_mva", "ELECTRON16_FR_MVA_TIGHT_CC",weight, true,false);

    GetFakeRateAndPromptRates("ELECTRON16_POG_FAKELOOSE_CC","dijet_pog", "ELECTRON16_FR_POG_TIGHT_DXYCC",weight, true,false);
    GetFakeRateAndPromptRates("ELECTRON16_POG_MEDIUM_FAKELOOSE_CC","dijet_pog", "ELECTRON16_FR_POG_MEDIUM_DXYCC",weight, true,false);
    GetFakeRateAndPromptRates("ELECTRON16_MVA_FAKELOOSE_CC","dijet_mva", "ELECTRON16_FR_MVA_TIGHT_DXYCC",weight, true,false);
    

    GetFakeRateAndPromptRates("ELECTRON16_POG_FAKELOOSE_d0","dijet_pog_d0", "ELECTRON16_FR_POG_TIGHT",weight, true,false);
    GetFakeRateAndPromptRates("ELECTRON16_POG_MEDIUM_FAKELOOSE_d0","dijet_pog_d0", "ELECTRON16_FR_POG_MEDIUM",weight, true,false);
    GetFakeRateAndPromptRates("ELECTRON16_POG_FAKELOOSE_CC_d0","dijet_pog_d0", "ELECTRON16_FR_POG_TIGHT_CC",weight, true,false);
    GetFakeRateAndPromptRates("ELECTRON16_POG_MEDIUM_FAKELOOSE_CC_d0","dijet_pog_d0", "ELECTRON16_FR_POG_MEDIUM_CC",weight, true,false);
    GetFakeRateAndPromptRates("ELECTRON16_MVA_FAKELOOSE_CC_d0","dijet_mva_d0", "ELECTRON16_FR_MVA_TIGHT_CC",weight, true,false);

    GetFakeRateAndPromptRates("ELECTRON16_POG_FAKELOOSE_CC_d0","dijet_pog_d0", "ELECTRON16_FR_POG_TIGHT_DXYCC",weight, true,false);
    GetFakeRateAndPromptRates("ELECTRON16_POG_MEDIUM_FAKELOOSE_CC_d0","dijet_pog_d0", "ELECTRON16_FR_POG_MEDIUM_DXYCC",weight, true,false);
    GetFakeRateAndPromptRates("ELECTRON16_MVA_FAKELOOSE_CC_d0","dijet_mva_d0", "ELECTRON16_FR_MVA_TIGHT_DXYCC",weight, true,false);


    /// use dxy method for systematic
    GetFakeRateAndPromptRates("ELECTRON_HN_HIGHDXY_FAKELOOSE","dxy", "ELECTRON_HN_HIGHDXY_TIGHT",weight, true,false);
    if(!isData)GetFakeRateAndPromptRates("ELECTRON_HN_LOWDXY_FAKELOOSE","dxy", "ELECTRON_HN_LOWDXY_TIGHT",weight, true,false);
    /// GetFakeRateAndPromptRatesPerPeriod fills same fake rate plots as GetFakeRateAndPromptRates but splits MC into 7 data periods
    GetFakeRateAndPromptRatesPerPeriod("ELECTRON16_HN_FAKELOOSE_NOD0","dijet_nod0", "ELECTRON16_HN_TIGHT",weight, true,false);
  }
  if((isData&&k_channel == "SingleElectron") || !isData){
    // MakeSingleElectronCRPlots makes single electron CR plots
    //MakeSingleElectronCRPlots("ELECTRON16_HN_FAKELOOSE_NOD0","dijet_nod0",  "ELECTRON16_HN_TIGHT",weight,true); 
  }
}
void FakeRateCalculator_El::MakeSingleElectronCRPlots(TString looseid, TString eltag, TString tightid, float w, bool usepujetid){
  std::vector<snu::KElectron> electronLooseColl = GetElectrons(false,false,  looseid);
  std::vector<snu::KElectron> electronTightColl = GetElectrons(false,false,  tightid);
  
  if (!k_isdata) {
    w = w * MCweight * eventbase->GetEvent().PileUpWeight();
    
  }

  float id_weight=1.;
  float reco_weight=1.;
  if(!isData){
    for(unsigned int iel=0; iel < electronLooseColl.size(); iel++){
      /// apply tight pog ID, this is not fully correct but as we have no SF and tight id is close to POG tight we use this                                                     
      id_weight*= mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronLooseColl);
      reco_weight*=  mcdata_correction->ElectronRecoScaleFactor(electronLooseColl);
    }
    w*= id_weight;
    w*= reco_weight;
  }
  /// USE HN jets. Add pileup ID? Currently                                                                                                                                     
  std::vector<snu::KJet> jetCollTight = GetJets("JET_HN");
  std::vector<snu::KJet> jetColl           = GetJets("JET_NOLEPTONVETO");
  if(usepujetid){
    jetCollTight = GetJets("JET_HN_PU");
    jetColl           = GetJets("JET_NOLEPTONVETO_PU");
  }

  if(electronLooseColl.size()<1) return;
  std::vector<snu::KMuon> muonColl = GetMuons("MUON_HN_LOOSE");  // loose selection                                                                                             
  if(muonColl.size() > 0) return;
  if(electronTightColl.size() ==1) {

    TString triggerslist_singlelep = "HLT_Ele32_eta2p1_WPTight_Gsf_v";

    /// HLT_Ele32_eta2p1_WPTight_Gsf_v not working beofere 805 cattuples                                                                                                                                                                                                                                                                                   
    float trigger_ps_singlelepweight= WeightByTrigger("HLT_Ele32_eta2p1_WPTight_Gsf_v", TargetLumi);
    
    if(PassTrigger(triggerslist_singlelep) ){
      FillCLHist(sighist_ee, "SingleElectron_unprescaled", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w*trigger_ps_singlelepweight);
    }
  }
}


void FakeRateCalculator_El::GetFakeRateAndPromptRatesPerPeriod(TString looseid, TString eltag, TString tightid, float w, bool usepujetid, bool runall){


  int iperiod = GetMCPeriodRandom();
  if(isData) iperiod = GetDataPeriod();
  TString speriod="";
  if (iperiod == 1)     speriod="B";
  if (iperiod == 2)     speriod="C";
  if (iperiod == 3)     speriod="D";
  if (iperiod == 4)     speriod="E";
  if (iperiod == 5)     speriod="F";
  if (iperiod == 6)     speriod="G";
  if (iperiod == 7)     speriod="H";
  SetupLuminosityMap(false, speriod);
  // setup correct pileup correcion + SFs
  
  GetFakeRateAndPromptRates(looseid, eltag+"_"+speriod, tightid, w, usepujetid, runall);
  
  /// Reset lumi map
  SetupLuminosityMap(false, "None");
  
}

void FakeRateCalculator_El::GetFakeRateAndPromptRates(TString looseid, TString eltag, TString tightid, float w, bool usepujetid, bool runall){

  // PileUpWeight is period dependant 
  // MC events are split into 7 data periods and 
  if (!k_isdata) {
    w = w * MCweight * eventbase->GetEvent().PeriodPileUpWeight(GetMCPeriodRandom());
  }

  /// Four single electron triggers
  TString triggerslist_12="HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  TString triggerslist_18="HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  TString triggerslist_23="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  TString triggerslist_33="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";

  // analysis trigger / diel trigger                                                                                                                                                                                                                                                               
  TString triggerslist="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
  
  /// GetElectrons returns vector of muons. 
  /// false, true arguments means means chargeflip electrons 
  std::vector<snu::KElectron> electronTightColl = GetElectrons(false,false, tightid);
  std::vector<snu::KElectron> electronLooseColltmp = GetElectrons(false,false,looseid);
  
  std::vector<snu::KElectron> electronLooseColl;

  if(eltag.Contains("dijet_d0")){
    for(unsigned int iel=0; iel < electronLooseColltmp.size(); iel++){
      if(fabs(electronLooseColltmp[iel].dxy()) <  0.02){
	electronLooseColl.push_back(electronLooseColltmp[iel]);
      }
    }
  }
  else electronLooseColl=electronLooseColltmp;
 

  /// electronLooseColl has 0.015 dxy cuts applied 
  
  /// electronLooseColltmp has 0.2 dxy cut

  float id_weight=1.;
  float reco_weight=1.;
  if(!isData && tightid==("ELECTRON_POG_TIGHT")){
    for(unsigned int iel=0; iel < electronLooseColl.size(); iel++){
      /// apply tight pog ID, this is not fully correct but as we have no SF and tight id is close to POG tight we use this
      id_weight*= mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronLooseColl);
      reco_weight*=  mcdata_correction->ElectronRecoScaleFactor(electronLooseColl);
    }
    w*= id_weight;
    w*= reco_weight;
  }

  /// USE HN jets. Add pileup ID? Currently 
  std::vector<snu::KJet> jetCollTight = GetJets("JET_HN");
  std::vector<snu::KJet> jetColl           = GetJets("JET_NOLEPTONVETO");
  if(usepujetid){
    jetCollTight = GetJets("JET_HN_PU");
    jetColl      = GetJets("JET_NOLEPTONVETO_PU");
  }

  /// remove events with no loose leptons
  if(electronLooseColl.size()<1) return;
  std::vector<snu::KMuon> muonColl = GetMuons("MUON_HN_LOOSE");  // loose selection
  if(muonColl.size() > 0) return;
    
  /// If HN selection make full set of plots
  if(runall){
    /// Dilepton plots to show DY is well described
    if(PassTrigger(triggerslist) ){
      if(electronTightColl.size() ==2) {
	if(electronTightColl.at(0).Pt() > 20. && electronTightColl.at(1).Pt() > 15. ){
	  float pr_weight= WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi);
	  FillCLHist(sighist_ee, "DiElectron", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w*pr_weight);
	}
      }
    }
  }
  
  /// Get prescale for single el event. Returns 1. or 0. for data
  float prescale_trigger =  GetPrescale(electronLooseColl,  PassTrigger(triggerslist_12), PassTrigger(triggerslist_18), PassTrigger( triggerslist_23), PassTrigger(triggerslist_33), TargetLumi); 
  
  /// Make standard plots for loose and tight collection dyx
  if(eltag.Contains("dxy") && !eltag.Contains("dxysig")){
    float diel_prescale = WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi);
    if(!PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) diel_prescale=0.;
					  
    MakeDXYFakeRatePlots(tightid, eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi), w);
    return;
  }
  
  /// Make standard plots for loose and tight collection dijet                                                                                                                   
  MakeFakeRatePlots(tightid, eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w, false);
  if(runall){

    TString dxysigtag="";
    if(eltag.Contains("dxysig"))dxysigtag+="_dxysig";

    TString dxysigtag_iso="";
    if(eltag.Contains("miniiso"))dxysigtag_iso="_miniiso";
    if(eltag.Contains("dxysig"))dxysigtag_iso+="_dxysig";

    // run this for loose ID only once
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b050_e100", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);    
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b050_e050", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);    
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b050_e040", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);    
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b050_e025", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);    
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b050_e020", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);    
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b050_e010", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);        
    
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b025_e100", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b025_e050", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b025_e040", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b025_e025", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b025_e020", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b025_e010", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);        

    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b015_e100", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b015_e050", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b015_e040", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b015_e025", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b015_e020", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b015_e010", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);        

    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b010_e100", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b010_e050", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b010_e040", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b010_e025", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b010_e020", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b010_e010", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);    
    

    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b017_e017", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b017_e015", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b017_e014", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b017_e013", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b017_e012", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b017_e011", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b017_e010", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);

    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b015_e017", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b015_e015", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b015_e014", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b015_e013", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b015_e012", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b015_e011", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b015_e010", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);

    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b014_e017", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b014_e015", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b014_e014", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b014_e013", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b014_e012", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b014_e011", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b014_e010", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);


    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b013_e017", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b013_e015", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b013_e014", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b013_e013", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b013_e012", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b013_e011", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b013_e010", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b012_e017", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b012_e015", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b012_e014", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b012_e013", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b012_e012", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b012_e011", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b012_e010", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);

    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b011_e017", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b011_e015", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b011_e014", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b011_e013", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b011_e012", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b011_e011", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b011_e010", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);

    
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b010_e017", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b010_e015", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b010_e014", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b010_e013", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b010_e012", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag + "_dxy_b010_e011", eltag, electronLooseColltmp,  jetCollTight, jetColl,  prescale_trigger, w,true);

    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b050_e050", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b050_e0525", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b050_e055", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b050_e060", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b050_e065", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b050_e075", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b050_e100", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b050_e125", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);

    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b0525_e050", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b0525_e0525", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b0525_e055", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b0525_e060", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b0525_e065", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b0525_e075", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b0525_e100", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b0525_e125", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);


    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b055_e050", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b055_e0525", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b055_e055", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b055_e060", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b055_e065", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b055_e075", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b055_e100", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b055_e125", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);


    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b060_e050", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b060_e0525", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b060_e055", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b060_e060", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b060_e065", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b060_e075", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b060_e100", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b060_e125", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);

    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b065_e050", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b065_e0525", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b065_e055", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b065_e060", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b065_e065", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b065_e075", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b065_e100", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b065_e125", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);

    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b075_e050", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b075_e0525", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b075_e055", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b075_e060", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b075_e065", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b075_e075", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b075_e100", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b075_e125", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);

    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b100_e050", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b100_e0525", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b100_e055", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b100_e060", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b100_e065", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b100_e075", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b100_e100", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b100_e125", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);

    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b125_e050", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b125_e0525", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b125_e055", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b125_e060", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b125_e065", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b125_e075", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b125_e100", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    MakeFakeRatePlots("HNTight" + dxysigtag_iso + "_b125_e125", eltag, electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, w,true);
    
  }
  
  if(eltag!="dijet_nod0")  return;

  if(electronLooseColl.size()!= 1) return;
  // now apply prescale to weight (0. or 1. for data)

  if(prescale_trigger == 0.) return;
  w*=prescale_trigger;

  bool truth_match= false;
  if(!k_isdata) {
    if(electronLooseColl.size() > 0){
      /// ONLY INCLUDE ELECTRONS FROM W/Z                                                                                                                                       
      if((electronLooseColl.at(0).MCMatched()))  truth_match=true;
    }
    else truth_match = true;
  }
  else truth_match=true;

  
  if(runall){

    if(electronLooseColl.size() == 1){
      if(jetCollTight.size() >=1) {
	if(jetCollTight.at(0).Pt() > 40.){
	  FillCLHist(sighist_ee, "SingleLooseElJet", eventbase->GetEvent(), muonColl,electronLooseColl,jetCollTight, w);
	}
      }
    }
    
    if(electronTightColl.size() == 1&& electronLooseColl.size() == 1&& muonColl.size()==0){
      if(jetCollTight.size() >=1) {
	if(jetCollTight.at(0).Pt() > 40.){
	  FillCLHist(sighist_ee, "SingleTightElJet", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w);
	}
      }
    }
    if(truth_match){
      if(electronLooseColl.size() == 1){
	if(jetCollTight.size() >=1) {
	  if(jetCollTight.at(0).Pt() > 40.){
	    FillCLHist(sighist_ee, "SingleLooseElJet_tm", eventbase->GetEvent(), muonColl,electronLooseColl,jetCollTight, w);
	  }
	}
      }

      if(electronTightColl.size() == 1&& electronLooseColl.size() == 1&& muonColl.size()==0){
	if(jetCollTight.size() >=1) {
	  if(jetCollTight.at(0).Pt() > 40.){
	    FillCLHist(sighist_ee, "SingleTightElJet_tm", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w);
	  }
	}
      }
    }
    
    
    Double_t MT=0;
    Double_t METdphi=0;
    for(unsigned int iel = 0; iel < electronTightColl.size();  iel++){
      METdphi = TVector2::Phi_mpi_pi(electronTightColl.at(iel).Phi()- eventbase->GetEvent().METPhi());
      MT = sqrt(2.* electronTightColl.at(iel).Et()*eventbase->GetEvent().PFMET() * (1 - cos( METdphi)));
    }
    if(eventbase->GetEvent().PFMET() > 40 && (60. < MT)  &&(MT < 100.) &&truth_match){
      
      if(electronTightColl.size() == 1&& jetCollTight.size() >= 1){
	if(jetCollTight.at(0).Pt() > 30.) {
	  FillCLHist(sighist_ee, "SingleTightElJet_prompt", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w);
	}
	if(jetCollTight.at(0).Pt() > 40.) {
	  FillCLHist(sighist_ee, "SingleTightElJet40_prompt", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w);
	}
      }
      if(jetCollTight.size() >= 1){
	if(jetCollTight.at(0).Pt() > 30.) {
	  FillCLHist(sighist_ee, "SingleLooseElJet_prompt", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w);
	}
	if(jetCollTight.at(0).Pt() > 40.) {
	  FillCLHist(sighist_ee, "SingleLooseElJet40_prompt", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w);
	}
      }
    }
    
    
    bool useevent40 = UseEvent(electronLooseColl , jetColl, 40., prescale_trigger, w);
    
    if(useevent40){
      if ( electronTightColl.size() == 1 && jetCollTight.size() >= 1)FillCLHist(sighist_ee, "TightEl30_notm", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w);
      if (jetCollTight.size() >= 1) FillCLHist(sighist_ee, "LooseEl30_notm", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w);
    
      if(truth_match){
	if ( electronTightColl.size() == 1 && jetCollTight.size() >= 1)FillCLHist(sighist_ee, "TightEl30", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w);
	if (jetCollTight.size() >= 1) FillCLHist(sighist_ee, "LooseEl30", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w);
	
      }
    }
    
  }/// run all plots maker
    
  
  return;
}// End of execute event loop


float FakeRateCalculator_El::GetPrescale( std::vector<snu::KElectron> electrons,  bool pass4, bool pass3, bool pass2, bool pass1, float fake_total_lum ){
  
  float prescale_trigger= 1.;
  if(electrons.size() ==1){

    /// 15 - 20  HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_*
    /// 20 - 25  HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_*
    /// 25 - 35  HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_* 
    /// 35 - inf HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_*

    
    if(electrons.at(0).Pt() >= 35.){
      //HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30
      if(pass1){
	if(isData) return 1.;
	prescale_trigger = WeightByTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v", fake_total_lum); //// 20 + GeV bins
      }
      else {
	if(isData) return 0;
	prescale_trigger = WeightByTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v", fake_total_lum)*0.8;
      }
    }
    else  if(electrons.at(0).Pt() >= 25.){
      //HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30

      if(pass2){
	if(isData) return 1.;
        prescale_trigger =  WeightByTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v", fake_total_lum) ; //// 20 + GeV bins
      }
      else {
	if(isData) return 0;
	prescale_trigger =  WeightByTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v", fake_total_lum) * 0.8; 
      }
    }
    else   if(electrons.at(0).Pt() >= 20.){
      //HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30
      if(pass3){
	if(isData) return 1.;
	prescale_trigger = WeightByTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v", fake_total_lum) ;
      }
      else {
	if(isData) return 0;
	prescale_trigger = WeightByTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v", fake_total_lum)*0.8 ;
      }
    }
    else   if(electrons.at(0).Pt() >= 15.){
      //HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_
      if(pass4){
	if(isData) return 1.;
        prescale_trigger = WeightByTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v", fake_total_lum) ;
      }
      else {
	if(isData) return 0;
        prescale_trigger = WeightByTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v", fake_total_lum)*0.8 ;
      }
    }
    else{
      prescale_trigger = 0.;
    }
  }
  if(prescale_trigger == 0.) return 0.;
  if(k_isdata) return 1.;
  if(electrons.at(0).Pt() < 10.) return 0.;
  

  
  return prescale_trigger;
}


void FakeRateCalculator_El::MakeDXYFakeRatePlots(TString label, TString eltag,  std::vector<snu::KElectron> electrons,  std::vector<snu::KJet> jets, std::vector<snu::KJet> alljets, float prescale_w, float prescale_dielw, float w){

  std::vector<snu::KElectron> electrons_tight =   GetElectrons(false,false,label);
  
  bool truth_match= false;
  /// FOR FAKE RATE SUBTRACTION NEED ONLY PROMPT ELECTRONS                                                                                                                      
  if(!k_isdata) {
    if(electrons.size() > 0){
      /// ONLY INCLUDE ELECTRONS FROM W/Z                                                                                                                                       
      if((electrons.at(0).MCMatched()))  truth_match=true;
    }
    else truth_match = true;
  }
  else truth_match=true;

  //// single el events 
  if(electrons.size()==1 && prescale_w == 0.) return;
  if(electrons.size()==2 && prescale_dielw == 0.) return;
  
  
  if(electrons.size()==1 && truth_match)GetFakeRates(electrons, electrons_tight, label, jets, alljets,  label+"_eldxy", (prescale_w * w),true);
  if(electrons.size()==2 && truth_match)GetFakeRates(electrons, electrons_tight, label, jets, alljets,  label+"_dieldxy", (prescale_dielw * w),true);
  if(!isData){
    if(electrons.size()==1 )GetFakeRates(electrons, electrons_tight, label, jets, alljets,  label+"_eldxy_notm", (prescale_w * w),true);
    if(electrons.size()==2 )GetFakeRates(electrons, electrons_tight, label, jets, alljets,  label+"_dieldxy_notm", (prescale_dielw * w),true);
  }
  
}



void FakeRateCalculator_El::MakeFakeRatePlots(TString label, TString eltag,  std::vector<snu::KElectron> electronstmp,  std::vector<snu::KJet> jets, std::vector<snu::KJet> alljets, float prescale_w, float w, bool makebasicplots){
  
  
  std::vector<snu::KElectron> electrons_tight;
  std::vector<snu::KElectron> electrons;

  if(label.Contains("dxy_b")){
    if(eltag.Contains("dijet_d0")){
      for(unsigned int iel=0; iel < electronstmp.size(); iel++){
	if(fabs(electronstmp[iel].dxy()) < GetDXYCut(electronstmp[iel],label)){
	  electrons.push_back(electronstmp[iel]);
	  float iso= electronstmp[iel].PFRelIso(0.3);
	  if(iso <   0.05) electrons_tight.push_back(electronstmp[iel]);
	}
      }
    }
    else{
      /// no d0 cuts in loose
      electrons= electronstmp;
      
      for(unsigned int iel=0; iel < electronstmp.size(); iel++){
	if(fabs(electronstmp[iel].dxy()) < GetDXYCut(electronstmp[iel],label)){
	  if(eltag.Contains("dxysig")){
	    if(fabs(electronstmp[iel].dxySig()) < 3.){
	      float iso= electronstmp[iel].PFRelIso(0.3);	      
	      if(iso <  0.05)electrons_tight.push_back(electronstmp[iel]);
	    }
	  }
	  else{
	    float iso= electronstmp[iel].PFRelIso(0.3);
	    if(iso <  0.05)electrons_tight.push_back(electronstmp[iel]);
	  }
	}
      }
    }
  }
  else if(label.Contains("_b")){
    if(eltag.Contains("dijet_d0")){
      for(unsigned int iel=0; iel < electronstmp.size(); iel++){
        if(fabs(electronstmp[iel].dxy()) < GetDXYCut(electronstmp[iel],"b020_e020")){
          electrons.push_back(electronstmp[iel]);
          float iso= electronstmp[iel].PFRelIso(0.3);
          if(eltag.Contains("miniiso")) iso=electronstmp[iel].PFRelMiniIso();
	  if(iso <   GetIsoCut(electronstmp[iel],label)) electrons_tight.push_back(electronstmp[iel]);
        }
      }
    }
    else{
      /// no d0 cuts in loose                                                                                                                                              
      electrons= electronstmp;

      for(unsigned int iel=0; iel < electronstmp.size(); iel++){
        if(fabs(electronstmp[iel].dxy()) < GetDXYCut(electronstmp[iel],"b020_e020")){
          if(eltag.Contains("dxysig")){
            if(fabs(electronstmp[iel].dxySig()) < 3.){
              float iso= electronstmp[iel].PFRelIso(0.3);
              if(eltag.Contains("miniiso")) iso=electronstmp[iel].PFRelMiniIso();
	      if(iso <   GetIsoCut(electronstmp[iel],label)) electrons_tight.push_back(electronstmp[iel]);
	      
            }
          }
          else{
            float iso= electronstmp[iel].PFRelIso(0.3);
            if(eltag.Contains("miniiso")) iso=electronstmp[iel].PFRelMiniIso();
	    if(iso <   GetIsoCut(electronstmp[iel],label)) electrons_tight.push_back(electronstmp[iel]);
          }
        }
      }
    }
  }
  
  else {
    electrons=electronstmp;
    electrons_tight =   GetElectrons(false,false ,label);
  }  

  if(electrons.size() ==2){
    if(Zcandidate(electrons, 10., true)){
      Float_t ptbins[10] = { 10., 15.,20.,25.,30.,35.,45.,60.,100., 200.};
      Float_t etabins[5] = { 0.,0.8,  1.479, 2.,  2.5};
      float pr_weight= WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi);
      FillHist(("Prompt_LooseEl_" + label  + "_" + eltag+ "_pt_eta").Data(), electrons[0].Pt(), fabs(electrons[0].Eta()),  weight*pr_weight, ptbins, 9 , etabins, 4);
      FillHist(("Prompt_LooseEl_" + label  + "_" + eltag+ "_pt_eta").Data(), electrons[1].Pt(), fabs(electrons[1].Eta()),  weight*pr_weight, ptbins, 9 , etabins, 4);
      if(PassID(electrons[0],label ))       FillHist(("Prompt_TightEl_" + label + "_" + eltag+ "_pt_eta").Data(), electrons[0].Pt(), fabs(electrons[0].Eta()),  weight*pr_weight, ptbins, 9 , etabins, 4);
      if(PassID(electrons[1],label ))       FillHist(("Prompt_TightEl_" +label  + "_" + eltag+ "_pt_eta").Data(), electrons[1].Pt(), fabs(electrons[1].Eta()),  weight*pr_weight, ptbins, 9 , etabins, 4);
      
    }
  }
  
  if(electrons.size() != 1 ) return;
  
  if(prescale_w==0.) return;
  
  bool useevent20 = UseEvent(electrons , jets, 20., prescale_w, w); 
  bool useevent30 = UseEvent(electrons , jets, 30., prescale_w, w); 
  bool useevent40 = UseEvent(electrons , jets, 40., prescale_w, w); 
  bool useevent60 = UseEvent(electrons , jets, 60., prescale_w, w); 

  bool truth_match= false;
  /// FOR FAKE RATE SUBTRACTION NEED ONLY PROMPT ELECTRONS
  if(!k_isdata) {
    if(electrons.size() > 0){
      /// ONLY INCLUDE ELECTRONS FROM W/Z
      if((electrons.at(0).MCMatched()))  truth_match=true;
    }
    else truth_match = true;
  }
  else truth_match=true;
  
  label=label + "_" + eltag;

  if(truth_match){
    if(jets.size() >= 1){
      if(useevent20)GetFakeRates(electrons, electrons_tight,label, jets, alljets,  label+"_20", (prescale_w * w),makebasicplots);
      if(useevent30)GetFakeRates(electrons, electrons_tight,label, jets, alljets,  label+"_30", (prescale_w * w),makebasicplots);
      if(useevent40)GetFakeRates(electrons, electrons_tight,label, jets, alljets,  label+"_40", (prescale_w * w),makebasicplots);
      if(useevent60)GetFakeRates(electrons, electrons_tight, label,jets, alljets,  label+"_60", (prescale_w * w),makebasicplots);
    }
  }
}

bool FakeRateCalculator_El::UseEvent(std::vector<snu::KElectron> electrons,  std::vector< snu::KJet> jets, float awayjetcut, float precale_weight, float wt){
  
  bool useevent = false;
  if(electrons.size() != 1) return false;
  if(precale_weight != 0.){
    /// Z and W veto
    Double_t MT=0;
    Double_t METdphi=0;
    for(unsigned int w = 0; w < electrons.size();  w++){
      METdphi = TVector2::Phi_mpi_pi(electrons.at(w).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet));
      MT = sqrt(2.* electrons.at(w).Et()*eventbase->GetEvent().MET(snu::KEvent::pfmet) * (1 - cos( METdphi)));

      if(( (eventbase->GetEvent().MET(snu::KEvent::pfmet) < 20) && (MT < 25.)) ) {
	
        for (unsigned int ielT=0; ielT < electrons.size(); ielT++){
          for(unsigned int ij=0; ij < jets.size(); ij++){
            if(jets.at(ij).Pt() < awayjetcut) continue;
            float dphi =fabs(TVector2::Phi_mpi_pi(electrons.at(ielT).Phi()- jets.at(ij).Phi()));
            if( (jets.at(ij).ChargedEMEnergyFraction()) > 0.65)  continue;
	    if(dphi > 2.5) useevent = true;
          }
        }
      }
    }
  }
  return useevent;
}

void FakeRateCalculator_El::GetFakeRates(std::vector<snu::KElectron> loose_el, std::vector<snu::KElectron> tight_el, TString tightlabel,  std::vector<snu::KJet> jets,  std::vector<snu::KJet> alljets, TString tag, double w, bool basicplots){
  
   
  Float_t ptbins[10] = { 10., 15.,20.,25.,30.,35.,45.,60.,100., 200.};
  Float_t ptbinsb[8] = { 10., 15.,20.,30.,45.,60.,100., 200.};
  Float_t etabin[2] = { 0.,  2.5};
  Float_t etabins[3] = { 0., 1.479,  2.5};
  Float_t etabins2[5] = { 0.,0.8,  1.479, 2.,  2.5};

  /// for most cuts just plot pt_eta
  if(basicplots){
    
    if(loose_el.size() == 1 && jets.size() >= 1){
      float el_pt = loose_el.at(0).Pt();
      FillHist(("LooseEl" + tag + "_pt_eta").Data(), el_pt, fabs(loose_el.at(0).Eta()),  w, ptbins, 9 , etabins2, 4);

      if( tight_el.size() == 1){
	FillHist(("TightEl" + tag + "_pt_eta").Data(), el_pt, fabs(tight_el.at(0).Eta()),  w, ptbins, 9 , etabins2, 4);
	
      }
    }
    
    if(loose_el.size() == 2 && jets.size() >= 1){
      float el1_pt = loose_el.at(0).Pt();
      float el2_pt = loose_el.at(1).Pt();
      FillHist(("LooseEl" + tag + "_pt_eta").Data(), el1_pt, fabs(loose_el.at(0).Eta()),  w, ptbins, 9 , etabins2, 4);
      FillHist(("LooseEl" + tag + "_pt_eta").Data(), el2_pt, fabs(loose_el.at(1).Eta()),  w, ptbins, 9 , etabins2, 4);

      if(PassID(loose_el[0], tightlabel)){
	FillHist(("TightEl" + tag + "_pt_eta").Data(), el1_pt, fabs(loose_el.at(0).Eta()),  w, ptbins, 9 , etabins2, 4);
      }
      if(PassID(loose_el[1], tightlabel)){
	FillHist(("TightEl" + tag + "_pt_eta").Data(), el2_pt, fabs(loose_el.at(1).Eta()),  w, ptbins, 9 , etabins2, 4);
      }
    }
    return;
  }
  
  int nbjet(0);
  float tmp_deltaR=1000.;
  for(unsigned int ij =0 ; ij < jets.size() ; ij++){
    if(jets.at(ij).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) nbjet++;
    if(loose_el.at(0).DeltaR(jets.at(ij)) < tmp_deltaR) tmp_deltaR = loose_el.at(0).DeltaR(jets.at(ij));
  }
 
  float awayjetpt = 0.;
  for(unsigned int ij =0 ; ij < alljets.size() ; ij++){
    float dphi =fabs(TVector2::Phi_mpi_pi(loose_el.at(0).Phi()- alljets.at(ij).Phi()));
    if(dphi > 2.5) awayjetpt = alljets.at(ij).Pt();
  }

  if(loose_el.size() == 1 && tight_el.size() == 1 && jets.size() >= 1){
   
    float el_pt = tight_el.at(0).Pt();
   
    // 1D FakeRates
    FillHist(("TightEl" + tag + "_eta").Data(), tight_el.at(0).Eta(), w, -2.5, 2.5,50);
    FillHist(("TightEl" + tag + "_pt").Data(),el_pt, w,  ptbins, 9);
   
    if(fabs(tight_el.at(0).Eta()) < 1.5)     FillHist(("TightEl" + tag + "_pt_barrel").Data(),el_pt, w,  ptbins, 9);
    else FillHist(("TightEl" + tag + "_pt_endcap").Data(),el_pt, w,  ptbins, 9);

    FillHist(("TightEl" + tag + "_njets").Data(), jets.size(), w, 0.,5.,5);
    FillHist(("TightEl" + tag + "_nvertices").Data(), eventbase->GetEvent().nVertices(), w, 0., 30., 30);
    FillHist(("TightEl" + tag + "_nbjet").Data(), nbjet, w, 0., 4.,4); 
    FillHist(("TightEl" + tag + "_pt_eta").Data(), el_pt, fabs(tight_el.at(0).Eta()),  w, ptbins, 9 , etabins2, 4);
   
    if(nbjet > 0){
      FillHist(("TightEl" + tag + "_bjet_eta").Data(), tight_el.at(0).Eta(), w, -2.5, 2.5,50);
      FillHist(("TightEl" + tag + "_bjet_pt").Data(), el_pt, w,  ptbins, 9);
      FillHist(("TightEl" + tag + "_bjet_pt_eta").Data(), el_pt, fabs(tight_el.at(0).Eta()),  w, ptbinsb, 7 , etabins, 2);
    }
    else{
      FillHist(("TightEl" + tag + "_0bjet_eta").Data(), tight_el.at(0).Eta(), w, -2.5, 2.5,50);
      FillHist(("TightEl" + tag + "_0bjet_pt").Data(), el_pt, w,  ptbins,9);
      FillHist(("TightEl" + tag + "_0bjet_pt_eta").Data(), el_pt, fabs(tight_el.at(0).Eta()),  w, ptbins, 9 , etabins2, 4);
    }
  }
   
  if(loose_el.size() == 1 && jets.size() >= 1){
    float el_pt = loose_el.at(0).Pt();
   
    FillHist(("LooseEl" + tag + "_eta").Data(), loose_el.at(0).Eta(), w, -2.5, 2.5,50);
    FillHist(("LooseEl" + tag + "_pt").Data(), el_pt, w,  ptbins, 9);
    FillHist(("LooseEl" + tag + "_njets").Data(), jets.size(), w, 0.,5.,5);
    FillHist(("LooseEl" + tag + "_nvertices").Data(), eventbase->GetEvent().nVertices(), w, 0., 30., 30);
    FillHist(("LooseEl" + tag + "_nbjet").Data(), nbjet, w, 0., 4.,4);
    FillHist(("LooseEl" + tag + "_pt_eta").Data(), el_pt, fabs(loose_el.at(0).Eta()),  w, ptbins, 9 , etabins2, 4);

    if(nbjet > 0){
      FillHist(("LooseEl" + tag + "_bjet_eta").Data(), loose_el.at(0).Eta(), w, -2.5, 2.5,50);
      FillHist(("LooseEl" + tag + "_bjet_pt").Data(), el_pt, w,  ptbins, 9);
      FillHist(("LooseEl" + tag + "_bjet_pt_eta").Data(), el_pt, fabs(loose_el.at(0).Eta()),  w, ptbinsb,7 , etabins, 2);
    }
    else{
      FillHist(("LooseEl" + tag + "_0bjet_eta").Data(), loose_el.at(0).Eta(), w, -2.5, 2.5,50);
      FillHist(("LooseEl" + tag + "_0bjet_pt").Data(), el_pt, w,  ptbins, 9);
      FillHist(("LooseEl" + tag + "_0bjet_pt_eta").Data(), el_pt, fabs(loose_el.at(0).Eta()),  w, ptbins, 9 , etabins2, 4);
    }

  }
  
  
  return;
  
}



void FakeRateCalculator_El::GetHSTRates(std::vector<snu::KElectron> loose_el, std::vector<snu::KElectron> tight_el, std::vector<snu::KJet> jets, TString tag){

  Float_t htbins[14] = { 20.,22.5, 25.,27.5, 30.,35.,40.,45.,50.,60.,80.,100.,200., 1000.};

  float tmp_deltaR=1000.;
  for(unsigned int ij =0 ; ij < jets.size() ; ij++){
      if(loose_el.at(0).DeltaR(jets.at(ij)) < tmp_deltaR) tmp_deltaR = loose_el.at(0).DeltaR(jets.at(ij));
  }
  
  if( tight_el.size() == 1 && jets.size() >= 1){
    FillHist(("MeasuredEl" + tag + "_njets").Data(), jets.size(), weight, 0.,5.,5);
    FillHist(("MeasuredEl" + tag + "_nvertices").Data(), eventbase->GetEvent().nVertices(), weight, 0., 30., 30);
    FillHist(("MeasuredEl" + tag + "_ht").Data(), SumPt(jets), weight, htbins, 13);
    FillHist(("MeasuredEl" + tag + "ejetdr").Data(), tmp_deltaR, weight , 0., 5., 50);
  }
  
  if( loose_el.size() == 1 && jets.size() >= 1){
    FillHist(("HSTLooseEl" + tag + "_njets").Data(), jets.size(), weight, 0.,5.,5);
    FillHist(("HSTLooseEl" + tag + "_nvertices").Data(), eventbase->GetEvent().nVertices(), weight, 0., 30., 30);
    FillHist(("HSTLooseEl" + tag + "_ht").Data(), SumPt(jets), weight, htbins, 13);
    FillHist(("HSTLooseEl" + tag + "ejetdr").Data(), tmp_deltaR, weight , 0., 5., 50);    
    
    float hstweight= m_fakeobj->getFakeRate_electronEta(0, loose_el.at(0).Pt(),loose_el.at(0).Eta(), tag);
    FillHist(("PredictedEl" + tag + "_njets").Data(), jets.size(), hstweight, 0.,5.,5);
    FillHist(("PredictedEl" + tag + "_nvertices").Data(), eventbase->GetEvent().nVertices(), hstweight, 0., 30., 30);
    FillHist(("PredictedEl" + tag + "_ht").Data(), SumPt(jets), hstweight, htbins, 13);
    FillHist(("PredictedEl" + tag + "ejetdr").Data(), tmp_deltaR, weight , 0., 5., 50);

  }  
}




void FakeRateCalculator_El::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger<< INFO << "Number of events that pass 1 7GeV trigger = " << n_17_pass  << LQLogger::endmsg;
  m_logger<< INFO << "Number of events that pass 17 GeV + jet trigger = " << n_17_jet_pass  << LQLogger::endmsg;
  m_logger<< INFO << "Number of events that pass 17 GeV || jet trigger = " << n_17_17_jet_pass  << LQLogger::endmsg;

}

void FakeRateCalculator_El::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //DeclareVariable(out_muons, "Signal_Muons");

  n_17_jet_pass=0;
  n_17_17_jet_pass=0;
  n_17_pass=0;

  
  return;
  
}

FakeRateCalculator_El::~FakeRateCalculator_El() {
  
  Message("In FakeRateCalculator_El Destructor" , INFO);
  
}



void FakeRateCalculator_El::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void FakeRateCalculator_El::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this FakeRateCalculator_ElCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void FakeRateCalculator_El::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


