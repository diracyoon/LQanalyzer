/***************************************************************************
 * @Project: LQFakeRateCalculator_El Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/


/// Local includes
#include "FakeRateCalculator_FinalEl.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (FakeRateCalculator_FinalEl);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
FakeRateCalculator_FinalEl::FakeRateCalculator_FinalEl() :  AnalyzerCore(),  out_electrons(0) {


  // To have the correct name in the log:                                                                                                                            
  SetLogName("FakeRateCalculator_FinalEl");

  Message("In FakeRateCalculator_FinalEl constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void FakeRateCalculator_FinalEl::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
  return;
  MakeCleverHistograms(sighist_ee, "SingleLooseElJet");
  MakeCleverHistograms(sighist_ee, "SingleTightElJet");

  MakeCleverHistograms(sighist_ee, "DiElectron");
  MakeCleverHistograms(sighist_ee, "DiElectron_SingleLeg");
  MakeCleverHistograms(sighist_ee, "Use40_loose");
  MakeCleverHistograms(sighist_ee, "Use40_tight");

  MakeCleverHistograms(sighist_ee, "SingleElectron_unprescaled");
  MakeCleverHistograms(sighist_ee, "SingleElectron_prompt_unprescaled");
  MakeCleverHistograms(sighist_ee, "SingleElectron_prompt_unprescaled2");
  return;
}


void FakeRateCalculator_FinalEl::ExecuteEvents()throw( LQError ){
    

  /// THIS COSE IS ONLY FOR GETTING FAKES USeD IN OPTIMISING ID.........

  Message("In ExecuteEvents() " , DEBUG);

  //// Initial event cuts
  /// MET FIleters 
  if(!PassMETFilter()) return;     
  
  Message("PASS MET FILTER " , DEBUG);

  /// Require good promary vertex 
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex  
  numberVertices = eventbase->GetEvent().nVertices();   

  if (!k_isdata) {
    //weight *= (MCweight * eventbase->GetEvent().PileUpWeight(snu::KEvent::down));
    weight *= (MCweight * mcdata_correction->UserPileupWeight(eventbase->GetEvent()));
  }

  Message("PASS VERTEX " , DEBUG);

  /// These run on double electron dataset
  if((isData&&k_channel == "DoubleEG") || !isData){

 
    TString triggerslist_8="HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v";   /// -> tighter cut in lepton ID form tighter trigger emulation cut                                     
    TString triggerslist_12="HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
    TString triggerslist_18="HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
    TString triggerslist_23="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
    TString triggerslist_33="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"; /// 

    std::vector<snu::KMuon> muonColl = GetMuons("MUON_HN_LOOSE");  // loose selection                                                                                                                                                                                              
    if(muonColl.size() > 0) return;
    
    std::vector<snu::KJet> jetCollTight = GetJets("JET_HN");

    std::vector<snu::KJet> jetColl           = GetJets("JET_NOLEPTONVETO");
    
    std::vector<snu::KElectron> loose_el_hn = GetElectrons(false,false,"ELECTRON_HN_FAKELOOSE");
    std::vector<snu::KElectron> loose_el_hn_iso04 = GetElectrons(false,false,"ELECTRON_HN_FAKELOOSE_ISO04");
    std::vector<snu::KElectron> loose_el_hn_iso08 = GetElectrons(false,false,"ELECTRON_HN_FAKELOOSE_ISO08");
    std::vector<snu::KElectron> tight_el_hn = GetElectrons(false,false,"ELECTRON_HN_TIGHTv4");
    GetFakeRateAndPromptRates(loose_el_hn,"HN_TIGHT", tight_el_hn,0.08,0.08,weight,"ELECTRON_HN_TIGHTv4",false, true);

    GetFakeRateAndPromptRatesPerPeriod(loose_el_hn,"HN_TIGHT", tight_el_hn,0.08,0.08,weight,"ELECTRON_HN_TIGHTv4",false, false);
    //
    GetFakeRateAndPromptRates(loose_el_hn_iso04,"HN_TIGHTISO04", tight_el_hn,0.08,0.08,weight,"ELECTRON_HN_TIGHTv4",false, false);
    GetFakeRateAndPromptRates(loose_el_hn_iso08,"HN_TIGHTISO08", tight_el_hn,0.08,0.08,weight,"ELECTRON_HN_TIGHTv4",false, false);
    GetFakeRateAndPromptRates(loose_el_hn,"HN_TIGHTPUJet", tight_el_hn,0.08,0.08,weight,"ELECTRON_HN_TIGHTv4",true, false);
    
    std::vector<snu::KElectron> loose_el_mva = GetElectrons(false,false,"ELECTRON_MVA_FAKELOOSE");
    std::vector<snu::KElectron> tight_el_mva = GetElectrons(false,false,"ELECTRON_MVA_TIGHT");
    GetFakeRateAndPromptRates(loose_el_mva,"MVA_TIGHT", tight_el_mva,0.05,0.05,weight, "ELECTRON_MVA_TIGHT", false, false);
    
    std::vector<snu::KElectron> loose_el_pogt = GetElectrons(false,false,"ELECTRON_POG_FAKELOOSE");
    std::vector<snu::KElectron> tight_el_pogt = GetElectrons(false,false,"ELECTRON_POG_TIGHT");
    GetFakeRateAndPromptRates(loose_el_pogt,"POG_TIGHT", tight_el_pogt,0.0588, 0.0571,weight, "ELECTRON_POG_TIGHT", false, false);
    
    std::vector<snu::KElectron> loose_el_pogm = GetElectrons(false,false,"ELECTRON_POG_MEDIUM_FAKELOOSE");
    std::vector<snu::KElectron> tight_el_pogm = GetElectrons(false,false,"ELECTRON_POG_MEDIUM");
    GetFakeRateAndPromptRates(loose_el_pogm,"POG_MEDIUM", tight_el_pogm,0.0695,0.0821,weight, "ELECTRON_POG_MEDIUM", false, false);
    
    std::vector<snu::KElectron> loose_el_gent = GetElectrons(false,false,"ELECTRON_GENT_FAKELOOSE");
    std::vector<snu::KElectron> tight_el_gent = GetElectrons(false,false,"ELECTRON_GENT_TIGHT");
    GetFakeRateAndPromptRates(loose_el_gent,"GENT", tight_el_gent,0.1,0.1,weight, "ELECTRON_GENT_TIGHT", false, false);
    

    // HN ID check others
    std::vector<snu::KElectron> loose_el_hn_mva = GetElectrons(false,false,"ELECTRON_HN_MVA_FAKELOOSE");
    std::vector<snu::KElectron> tight_el_hn_mva = GetElectrons(false,false,"ELECTRON_HN_MVA_TIGHT");
    GetFakeRateAndPromptRates(loose_el_hn_mva,"HN_MVA_TIGHT", tight_el_hn_mva,0.08,0.08,weight, "ELECTRON_HN_MVA_TIGHT", false, false);

    std::vector<snu::KElectron> loose_el_hn_pogt = GetElectrons(false,false,"ELECTRON_HN_POG_FAKELOOSE");
    std::vector<snu::KElectron> tight_el_hn_pogt = GetElectrons(false,false,"ELECTRON_HN_POG_TIGHT");
    GetFakeRateAndPromptRates(loose_el_hn_pogt,"HN_POG_TIGHT", tight_el_hn_pogt,0.08, 0.08,weight, "ELECTRON_HN_POG_TIGHT", false, false);

    std::vector<snu::KElectron> loose_el_hn_pogm = GetElectrons(false,false,"ELECTRON_HN_POG_MEDIUM_FAKELOOSE");
    std::vector<snu::KElectron> tight_el_hn_pogm = GetElectrons(false,false,"ELECTRON_HN_POG_MEDIUM");
    GetFakeRateAndPromptRates(loose_el_hn_pogm,"HN_POG_MEDIUM", tight_el_hn_pogm,0.08,0.08,weight, "ELECTRON_HN_POG_MEDIUM", false, false);

    std::vector<snu::KElectron> loose_el_hn_gent = GetElectrons(false,false,"ELECTRON_HN_GENT_FAKELOOSE");
    std::vector<snu::KElectron> tight_el_hn_gent = GetElectrons(false,false,"ELECTRON_HN_GENT_TIGHT");
    GetFakeRateAndPromptRates(loose_el_hn_gent,"HN_GENT", tight_el_hn_gent,0.08,0.08,weight, "ELECTRON_HN_GENT_TIGHT", false, false);


    /// opt

    std::vector<snu::KElectron> tight_el_hn_iso05 = GetElectrons(false,false,"ELECTRON_HN_TIGHTv4_iso05");
    std::vector<snu::KElectron> tight_el_hn_iso06 = GetElectrons(false,false,"ELECTRON_HN_TIGHTv4_iso06");
    std::vector<snu::KElectron> tight_el_hn_iso07 = GetElectrons(false,false,"ELECTRON_HN_TIGHTv4_iso07");
    std::vector<snu::KElectron> tight_el_hn_iso09 = GetElectrons(false,false,"ELECTRON_HN_TIGHTv4_iso09");
    std::vector<snu::KElectron> tight_el_hn_iso10 = GetElectrons(false,false,"ELECTRON_HN_TIGHTv4_iso10");
    GetFakeRateAndPromptRates(loose_el_hn,"HN_TIGHT_iso05", tight_el_hn_iso05,0.05,0.05,weight,"ELECTRON_HN_TIGHTv4_iso05",false, false);
    GetFakeRateAndPromptRates(loose_el_hn,"HN_TIGHT_iso06", tight_el_hn_iso06,0.06,0.06,weight,"ELECTRON_HN_TIGHTv4_iso06",false, false);
    GetFakeRateAndPromptRates(loose_el_hn,"HN_TIGHT_iso07", tight_el_hn_iso07,0.07,0.07,weight,"ELECTRON_HN_TIGHTv4_iso07",false, false);
    GetFakeRateAndPromptRates(loose_el_hn,"HN_TIGHT_iso09", tight_el_hn_iso09,0.09,0.09,weight,"ELECTRON_HN_TIGHTv4_iso09",false, false);
    GetFakeRateAndPromptRates(loose_el_hn,"HN_TIGHT_iso10", tight_el_hn_iso10,0.1,0.1,weight,"ELECTRON_HN_TIGHTv4_iso10",false, false);
    
    std::vector<snu::KElectron> tight_el_hn_dxy005 = GetElectrons(false,false,"ELECTRON_HN_TIGHTv4_dxy_005");
    std::vector<snu::KElectron> tight_el_hn_dxy02 = GetElectrons(false,false,"ELECTRON_HN_TIGHTv4_dxy_02");
    std::vector<snu::KElectron> tight_el_hn_dxy03 = GetElectrons(false,false,"ELECTRON_HN_TIGHTv4_dxy_03");
    std::vector<snu::KElectron> tight_el_hn_dxy04 = GetElectrons(false,false,"ELECTRON_HN_TIGHTv4_dxy_04");
    std::vector<snu::KElectron> tight_el_hn_dxy05 = GetElectrons(false,false,"ELECTRON_HN_TIGHTv4_dxy_05");
    
    std::vector<snu::KElectron> loose_el_hn_dxy005 = GetElectrons(false,false,"ELECTRON_HN_FAKELOOSE_dxy005");
    std::vector<snu::KElectron> loose_el_hn_dxy02 = GetElectrons(false,false,"ELECTRON_HN_FAKELOOSE_dxy02");
    std::vector<snu::KElectron> loose_el_hn_dxy03 = GetElectrons(false,false,"ELECTRON_HN_FAKELOOSE_dxy03");
    std::vector<snu::KElectron> loose_el_hn_dxy04 = GetElectrons(false,false,"ELECTRON_HN_FAKELOOSE_dxy04");
    std::vector<snu::KElectron> loose_el_hn_dxy05 = GetElectrons(false,false,"ELECTRON_HN_FAKELOOSE_dxy05");
    
    GetFakeRateAndPromptRates(loose_el_hn,"HN_TIGHT_dxy005", tight_el_hn_dxy005,0.08,0.08,weight,"ELECTRON_HN_TIGHTv4_dxy_005",false, false);
    GetFakeRateAndPromptRates(loose_el_hn,"HN_TIGHT_dxy02", tight_el_hn_dxy02,0.08,0.08,weight,"ELECTRON_HN_TIGHTv4_dxy_02",false, false);
    GetFakeRateAndPromptRates(loose_el_hn,"HN_TIGHT_dxy03", tight_el_hn_dxy03,0.08,0.08,weight,"ELECTRON_HN_TIGHTv4_dxy_03",false, false);
    GetFakeRateAndPromptRates(loose_el_hn,"HN_TIGHT_dxy04", tight_el_hn_dxy04,0.08,0.08,weight,"ELECTRON_HN_TIGHTv4_dxy_04",false, false);
    GetFakeRateAndPromptRates(loose_el_hn,"HN_TIGHT_dxy05", tight_el_hn_dxy05,0.08,0.08,weight,"ELECTRON_HN_TIGHTv4_dxy_05",false, false);

    

    
  }
  

  if((isData&&k_channel == "SingleElectron") || !isData){

    std::vector<snu::KElectron> loose_el_hn = GetElectrons(false,false,"ELECTRON_HN_FAKELOOSE");
    std::vector<snu::KElectron> tight_el_hn = GetElectrons(false,false,"ELECTRON_HN_TIGHTv4");

    TString triggerslist="HLT_Ele27_WPTight_Gsf_v";

    bool passtrig =  PassTrigger(triggerslist);
    if(!passtrig) return;
    GetFakeRateAndPromptRates(loose_el_hn,"singleel_HN_TIGHT",tight_el_hn,0.08,0.08,weight,"ELECTRON_HN_TIGHTv4",false, true);
    std::vector<snu::KElectron> loose_el_mva = GetElectrons(false,false,"ELECTRON_MVA_FAKELOOSE");
    std::vector<snu::KElectron> tight_el_mva = GetElectrons(false,false,"ELECTRON_MVA_TIGHT");
    GetFakeRateAndPromptRates(loose_el_mva,"singleel_MVA_TIGHT", tight_el_mva,0.05,0.05,weight,"ELECTRON_MVA_TIGHT", false, false);
    

  }
}
void FakeRateCalculator_FinalEl::MakeSingleElectronCRPlots(TString looseid, TString eltag, TString tightid, float w, bool usepujetid){
  std::vector<snu::KElectron> electronLooseColl = GetElectrons(false,false,  looseid);
  std::vector<snu::KElectron> electronTightColl = GetElectrons(false,false,  tightid);
  

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

  std::vector<snu::KJet> jetColl           = GetJets("JET_NOLEPTONVETO");
  std::vector<snu::KJet> jetCollTight = GetJets("JET_HN");
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
      Double_t MT=0;
      Double_t METdphi=0;
      for(unsigned int iel = 0; iel < electronTightColl.size();  iel++){
	METdphi = TVector2::Phi_mpi_pi(electronTightColl.at(iel).Phi()- eventbase->GetEvent().METPhi());
	MT = sqrt(2.* electronTightColl.at(iel).Et()*eventbase->GetEvent().PFMET() * (1 - cos( METdphi)));
      }


      bool truth_match= false;
      if(!k_isdata) {
	if(electronTightColl.size() > 0){
	  if((electronTightColl.at(0).MCMatched()))  truth_match=true;
	}
	else truth_match = true;
      }
      else truth_match=true;
      
      
      if(eventbase->GetEvent().PFMET() > 30 && (60. < MT)  &&(MT < 100.) &&truth_match)       FillCLHist(sighist_ee, "SingleElectron_prompt_unprescaled", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w*trigger_ps_singlelepweight);
      if(eventbase->GetEvent().PFMET() > 20 && (70. < MT)  &&(MT < 120.) &&truth_match)       FillCLHist(sighist_ee," SingleElectron_prompt_unprescaled2", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w*trigger_ps_singlelepweight);
      


    }
  }
}


void FakeRateCalculator_FinalEl::GetFakeRateAndPromptRatesPerPeriod(std::vector<snu::KElectron> electronLooseColl, TString eltag, std::vector<snu::KElectron> electronTightColl, float isocutb, float isocute,  float w, TString tightid,bool usepujetid, bool runall){


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
  
  GetFakeRateAndPromptRates(electronLooseColl, eltag+"_"+speriod,electronTightColl,isocutb,isocute, w, tightid, usepujetid, false);
  
  /// Reset lumi map
  SetupLuminosityMap(false, "None");
  
}

void FakeRateCalculator_FinalEl::GetFakeRateAndPromptRates(std::vector<snu::KElectron> electronLooseColl, TString eltag, std::vector<snu::KElectron> electronTightColl, float isocutb, float isocute,  float w, TString tightid,bool usepujetid, bool runall){

  bool single_lep=false;

  if((isData&&k_channel == "SingleElectron") || !isData) single_lep=true;
  if(!eltag.Contains("singleel")) single_lep=false;
  
  /// Four single electron triggers
  
  /// Loose ID has emul. tighter than CaloIdL_TrackIdL_IsoVL for all pt range
  // 10-15 has WPLoose https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#HLT_safe_selection_for_2016_data
  // 15-inf has CaloIdL_TrackIdL_IsoVL_ID  https://twiki.cern.ch/twiki/bin/view/CMS/ChangesEGMHLTAlgo2014
  TString triggerslist_8="HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v";   /// -> tighter cut in lepton ID form tighter trigger emulation cut
  TString triggerslist_12="HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  TString triggerslist_18="HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  TString triggerslist_23="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  TString triggerslist_33="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"; /// --> no 33 trigger

  TString triggerslist="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
  

  float id_weight=1.;
  float reco_weight=1.;
  if(!isData) w*= mcdata_correction->ElectronRecoScaleFactor(electronLooseColl);

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
  
   /// Get prescale for single el event. Returns 1. or 0. for data

  float prescale_trigger =  GetPrescale(electronLooseColl,  PassTrigger(triggerslist_8), PassTrigger(triggerslist_12), PassTrigger(triggerslist_18), PassTrigger( triggerslist_23), PassTrigger(triggerslist_33), TargetLumi); 
  if(single_lep) {
    prescale_trigger = WeightByTrigger("HLT_Ele27_WPTight_Gsf_v",TargetLumi);
    if(electronLooseColl[0].Pt() < 30) return;
  }
  
  /// Make standard plots for loose and tight collection dijet                                                                                                                  
  MakeFakeRatePlots(tightid, eltag, electronTightColl,electronLooseColl,  jetCollTight, jetColl,  prescale_trigger, isocutb, isocute,w, !runall);

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
    if(eventbase->GetEvent().PFMET() > 30 && (60. < MT)  &&(MT < 100.) &&truth_match){
      
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
      if (electronTightColl.size() == 1 && jetCollTight.size() >= 1)FillCLHist(sighist_ee, "TightEl30_notm", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w);
      if (jetCollTight.size() >= 1) FillCLHist(sighist_ee, "LooseEl30_notm", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w);
    
      if(truth_match){
	if ( electronTightColl.size() == 1 && jetCollTight.size() >= 1)FillCLHist(sighist_ee, "TightEl30", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w);
	if (jetCollTight.size() >= 1) FillCLHist(sighist_ee, "LooseEl30", eventbase->GetEvent(), muonColl,electronTightColl,jetCollTight, w);
	
      }
    }
    
  }/// run all plots maker
    
  
  return;
}// End of execute event loop


float FakeRateCalculator_FinalEl::GetPrescale( std::vector<snu::KElectron> electrons,bool pass5,  bool pass4, bool pass3, bool pass2, bool pass1, float fake_total_lum ){
  
  float prescale_trigger= 1.;
  if(electrons.size() != 1) return 0.;
  if(electrons.size() ==1){

    /// 10 - 15  HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_*
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
	return 0;
	//prescale_trigger = WeightByTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v", fake_total_lum)*0.8;
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
	return 0;
	//prescale_trigger =  WeightByTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v", fake_total_lum) * 0.8; 
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
	return 0;
	//prescale_trigger = WeightByTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v", fake_total_lum)*0.8 ;
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
	return 0;
        //prescale_trigger = WeightByTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v", fake_total_lum)*0.8 ;
      }
    }
    else   if(electrons.at(0).Pt() >= 8.){
      //HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_                                                                                                                               
      if(pass5){
        if(isData) return 1.;
        prescale_trigger = WeightByTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v", fake_total_lum) ;
      }
      else {
        if(isData) return 0;
	return 0;
        //prescale_trigger = WeightByTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v", fake_total_lum)*0.8 ;
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


void FakeRateCalculator_FinalEl::MakeDXYFakeRatePlots(TString label, TString eltag,  std::vector<snu::KElectron> electrons,  std::vector<snu::KJet> jets, std::vector<snu::KJet> alljets, float prescale_w, float prescale_dielw, float w){

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
  if(electrons.size()==2 ){
    if(electrons[0].Pt()  < 35) return;
  }
  
  
}



void FakeRateCalculator_FinalEl::MakeFakeRatePlots(TString label, TString eltag,   std::vector<snu::KElectron> electrons_tight, std::vector<snu::KElectron> electrons,  std::vector<snu::KJet> jets, std::vector<snu::KJet> alljets, float prescale_w, float isocutb, float isocute ,float w, bool makebasicplots){
  


  if(electrons.size() ==2){
    if(Zcandidate(electrons, 10., true)){
      Float_t ptbins[11] = { 10., 15.,20.,25.,30.,40.,45., 50.,70.,100., 200.};
      
      Float_t etabins[5] = { 0.,0.8,  1.479, 2.,  2.5};
      float pr_weight= WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi);
      FillHist(("Prompt_LooseEl_" + label  + "_" + eltag+ "_pt_eta").Data(), electrons[0].Pt(), fabs(electrons[0].Eta()),  weight*pr_weight, ptbins, 10 , etabins, 4);
      FillHist(("Prompt_LooseEl_" + label  + "_" + eltag+ "_pt_eta").Data(), electrons[1].Pt(), fabs(electrons[1].Eta()),  weight*pr_weight, ptbins, 10 , etabins, 4);
      if(PassID(electrons[0],label ))       FillHist(("Prompt_TightEl_" + label + "_" + eltag+ "_pt_eta").Data(), electrons[0].Pt(), fabs(electrons[0].Eta()),  weight*pr_weight, ptbins, 10 , etabins, 4);
      if(PassID(electrons[1],label ))       FillHist(("Prompt_TightEl_" +label  + "_" + eltag+ "_pt_eta").Data(), electrons[1].Pt(), fabs(electrons[1].Eta()),  weight*pr_weight, ptbins, 10, etabins, 4);
      
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
  
  label= eltag;
  if(truth_match){
    if(jets.size() >= 1){
      if(useevent20)GetFakeRates(electrons, electrons_tight,label, jets, alljets,  label+"_20", isocutb, isocute,(prescale_w * w),makebasicplots);
      if(useevent30)GetFakeRates(electrons, electrons_tight,label, jets, alljets,  label+"_30", isocutb, isocute,(prescale_w * w),makebasicplots);
      if(useevent40)GetFakeRates(electrons, electrons_tight,label, jets, alljets,  label+"_40", isocutb, isocute,(prescale_w * w),makebasicplots);
      if(useevent60)GetFakeRates(electrons, electrons_tight, label,jets, alljets,  label+"_60", isocutb, isocute,(prescale_w * w),makebasicplots);
    }
  }
}

bool FakeRateCalculator_FinalEl::UseEvent(std::vector<snu::KElectron> electrons,  std::vector< snu::KJet> jets, float awayjetcut, float precale_weight, float wt){
  
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

void FakeRateCalculator_FinalEl::GetFakeRates(std::vector<snu::KElectron> loose_el, std::vector<snu::KElectron> tight_el, TString tightlabel,  std::vector<snu::KJet> jets,  std::vector<snu::KJet> alljets, TString tag,float isocutb, float isocute, double w, bool basicplots){
  
  Float_t ptbins[11] = { 10., 15.,20.,25.,30.,40.,45., 50.,70.,100., 200.};
  Float_t htbins[7] = { 20.,30., 50.,70.,100., 200.,1000.};
  Float_t ptbinsb[8] = { 10., 15.,20.,30.,45.,60.,100., 200.};
  Float_t etabin[2] = { 0.,  2.5};
  Float_t etabins[4] = { 0., 0.8,1.479,  2.5};
  Float_t etabins2[5] = { 0.,0.8,  1.479, 2.,  2.5};

  /// for most cuts just plot pt_eta
  
  bool closebjet_m = false;
  for(unsigned int ij =0 ; ij < alljets.size() ; ij++){
    if(loose_el.at(0).DeltaR(alljets.at(ij)) < 0.5) {
      if(alljets.at(ij).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) closebjet_m = true;
    }
  }
  
  
  float el_pt = loose_el.at(0).Pt();
  float isocut = isocutb;
  if(fabs(loose_el.at(0).Eta()) > 1.479)isocut = isocute;
  float el_pt_corr = loose_el.at(0).Pt()*(1+max(0.,(loose_el.at(0).PFRelIso(0.3)-isocut))) ; /// will need changing for systematics                                                                                                                                           

  if(loose_el.size() == 1 && jets.size() >= 1){
    
    if(el_pt > 10.)FillHist(("LooseEl" + tag + "_pt_eta").Data(), el_pt, fabs(loose_el.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
    if(el_pt_corr > 10.)FillHist(("LooseEl" + tag + "_ptcorr_eta").Data(), el_pt_corr, fabs(loose_el.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
    if(closebjet_m) {
      if(el_pt > 10.)FillHist(("LooseEl" + tag + "_pt_eta_cbj").Data(), el_pt, fabs(loose_el.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
      if(el_pt_corr > 10.)FillHist(("LooseEl" + tag + "_ptcorr_eta_cbj").Data(), el_pt_corr, fabs(loose_el.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
    }
    else{
      if(el_pt > 10.)FillHist(("LooseEl" + tag + "_pt_eta_nocbj").Data(), el_pt, fabs(loose_el.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
      if(el_pt_corr > 10.)FillHist(("LooseEl" + tag + "_ptcorr_eta_nocbj").Data(), el_pt_corr, fabs(loose_el.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
    }
    if( tight_el.size() == 1){
      if(el_pt > 10.)FillHist(("TightEl" + tag + "_pt_eta").Data(), el_pt, fabs(tight_el.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
      if(el_pt_corr > 10.)FillHist(("TightEl" + tag + "_ptcorr_eta").Data(), el_pt_corr, fabs(tight_el.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
      if(closebjet_m) {
	if(el_pt > 10.)FillHist(("TightEl" + tag + "_pt_eta_cbj").Data(), el_pt, fabs(tight_el.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
	if(el_pt_corr > 10.)FillHist(("TightEl" + tag + "_ptcorr_eta_cbj").Data(), el_pt_corr, fabs(tight_el.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);

      }
      else{
	if(el_pt > 10.)FillHist(("TightEl" + tag + "_pt_eta_cbj").Data(), el_pt, fabs(tight_el.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
	if(el_pt_corr > 10.)FillHist(("TightEl" + tag + "_ptcorr_eta_cbj").Data(), el_pt_corr, fabs(tight_el.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);

      }
    }
  }
  
  if(basicplots) return;
  

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
   
   
    // 1D FakeRates
    FillHist(("TightEl" + tag + "_eta").Data(), tight_el.at(0).Eta(), w, -2.5, 2.5,50);
    FillHist(("TightEl" + tag + "_pt").Data(),el_pt, w,  ptbins, 10);
    FillHist(("TightEl" + tag + "_ptcorr").Data(),el_pt_corr, w,  ptbins, 10);
   
    if(fabs(tight_el.at(0).Eta()) < 1.5)     FillHist(("TightEl" + tag + "_pt_barrel").Data(),el_pt, w,  ptbins, 10);
    else FillHist(("TightEl" + tag + "_pt_endcap").Data(),el_pt, w,  ptbins, 10);

    if(fabs(tight_el.at(0).Eta()) < 1.5)     FillHist(("TightEl" + tag + "_ptcorr_barrel").Data(),el_pt_corr, w,  ptbins, 10);
    else FillHist(("TightEl" + tag + "_ptcorr_endcap").Data(),el_pt_corr ,w,  ptbins,10);


    FillHist(("TightEl" + tag + "_njets").Data(), jets.size(), w, 0.,5.,5);
    FillHist(("TightEl" + tag + "_nvertices").Data(), eventbase->GetEvent().nVertices(), w, 0., 30., 30);
    FillHist(("TightEl" + tag + "_nbjet").Data(), nbjet, w, 0., 4.,4); 
    FillHist(("TightEl" + tag + "_ht").Data(), SumPt(jets), w, htbins, 6); 
    
    if(nbjet > 0){
      FillHist(("TightEl" + tag + "_bjet_eta").Data(), tight_el.at(0).Eta(), w, -2.5, 2.5,50);
      FillHist(("TightEl" + tag + "_bjet_pt").Data(), el_pt, w,  ptbins, 10);
      FillHist(("TightEl" + tag + "_bjet_ptcorr").Data(), el_pt_corr, w,  ptbins, 10);
      FillHist(("TightEl" + tag + "_bjet_pt_eta").Data(), el_pt, fabs(tight_el.at(0).Eta()),  w, ptbinsb, 7 , etabins, 2);
      FillHist(("TightEl" + tag + "_bjet_ptcorr_eta").Data(),el_pt_corr, fabs(tight_el.at(0).Eta()),  w, ptbinsb, 7 , etabins, 2);
    }
    else{
      FillHist(("TightEl" + tag + "_0bjet_eta").Data(), tight_el.at(0).Eta(), w, -2.5, 2.5,50);
      FillHist(("TightEl" + tag + "_0bjet_pt").Data(), el_pt, w,  ptbins,9);
      FillHist(("TightEl" + tag + "_0bjet_ptcorr").Data(), el_pt_corr, w,  ptbins,9);
      FillHist(("TightEl" + tag + "_0bjet_pt_eta").Data(), el_pt, fabs(tight_el.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
      FillHist(("TightEl" + tag + "_0bjet_ptcorr_eta").Data(), el_pt_corr, fabs(tight_el.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
    }
  }
   
  if(loose_el.size() == 1 && jets.size() >= 1){
    
    FillHist(("LooseEl" + tag + "_eta").Data(), loose_el.at(0).Eta(), w, -2.5, 2.5,50);
    FillHist(("LooseEl" + tag + "_pt").Data(), el_pt, w,  ptbins, 10);
    FillHist(("LooseEl" + tag + "_ptcorr").Data(), el_pt_corr, w,  ptbins, 10);
    FillHist(("LooseEl" + tag + "_njets").Data(), jets.size(), w, 0.,5.,5);
    FillHist(("LooseEl" + tag + "_nvertices").Data(), eventbase->GetEvent().nVertices(), w, 0., 30., 30);
    FillHist(("LooseEl" + tag + "_nbjet").Data(), nbjet, w, 0., 4.,4);
    FillHist(("LooseEl" + tag + "_ht").Data(), SumPt(jets), w, htbins, 6);

    if(fabs(loose_el.at(0).Eta()) < 1.5)     FillHist(("LooseEl" + tag + "_pt_barrel").Data(),el_pt, w,  ptbins, 10);
    else FillHist(("LooseEl" + tag + "_pt_endcap").Data(),el_pt, w,  ptbins, 10);

    if(fabs(loose_el.at(0).Eta()) < 1.5)     FillHist(("LooseEl" + tag + "_ptcorr_barrel").Data(),el_pt_corr, w,  ptbins, 10);
    else FillHist(("LooseEl" + tag + "_ptcorr_endcap").Data(),el_pt_corr, w,ptbins, 10);



    if(nbjet > 0){
      FillHist(("LooseEl" + tag + "_bjet_eta").Data(), loose_el.at(0).Eta(), w, -2.5, 2.5,50);
      FillHist(("LooseEl" + tag + "_bjet_pt").Data(), el_pt, w,  ptbins, 10);
      FillHist(("LooseEl" + tag + "_bjet_ptcorr").Data(), el_pt_corr, w,  ptbins, 10);
      FillHist(("LooseEl" + tag + "_bjet_pt_eta").Data(), el_pt, fabs(loose_el.at(0).Eta()),  w, ptbinsb,7 , etabins, 2);
      FillHist(("LooseEl" + tag + "_bjet_ptcorr_eta").Data(), el_pt_corr, fabs(loose_el.at(0).Eta()),  w, ptbinsb,7 , etabins, 2);
    }
    else{
      FillHist(("LooseEl" + tag + "_0bjet_eta").Data(), loose_el.at(0).Eta(), w, -2.5, 2.5,50);
      FillHist(("LooseEl" + tag + "_0bjet_pt").Data(), el_pt, w,  ptbins, 10);
      FillHist(("LooseEl" + tag + "_0bjet_ptcorr").Data(), el_pt_corr, w,  ptbins, 10);
      FillHist(("LooseEl" + tag + "_0bjet_pt_eta").Data(), el_pt, fabs(loose_el.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
      FillHist(("LooseEl" + tag + "_0bjet_ptcorr_eta").Data(), el_pt_corr, fabs(loose_el.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
    }

  }
  
  
  return;
  
}



void FakeRateCalculator_FinalEl::GetHSTRates(std::vector<snu::KElectron> loose_el, std::vector<snu::KElectron> tight_el, std::vector<snu::KJet> jets, TString tag){

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
    
    float hstweight= m_fakeobj->getFakeRate( loose_el.at(0).Pt(),loose_el.at(0).Eta(), tag,15.);
    FillHist(("PredictedEl" + tag + "_njets").Data(), jets.size(), hstweight, 0.,5.,5);
    FillHist(("PredictedEl" + tag + "_nvertices").Data(), eventbase->GetEvent().nVertices(), hstweight, 0., 30., 30);
    FillHist(("PredictedEl" + tag + "_ht").Data(), SumPt(jets), hstweight, htbins, 13);
    FillHist(("PredictedEl" + tag + "ejetdr").Data(), tmp_deltaR, weight , 0., 5., 50);

  }  
}




void FakeRateCalculator_FinalEl::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger<< INFO << "Number of events that pass 1 7GeV trigger = " << n_17_pass  << LQLogger::endmsg;
  m_logger<< INFO << "Number of events that pass 17 GeV + jet trigger = " << n_17_jet_pass  << LQLogger::endmsg;
  m_logger<< INFO << "Number of events that pass 17 GeV || jet trigger = " << n_17_17_jet_pass  << LQLogger::endmsg;

}

void FakeRateCalculator_FinalEl::BeginCycle() throw( LQError ){
  
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

FakeRateCalculator_FinalEl::~FakeRateCalculator_FinalEl() {
  
  Message("In FakeRateCalculator_FinalEl Destructor" , INFO);
  
}



void FakeRateCalculator_FinalEl::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void FakeRateCalculator_FinalEl::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this FakeRateCalculator_FinalElCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void FakeRateCalculator_FinalEl::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


