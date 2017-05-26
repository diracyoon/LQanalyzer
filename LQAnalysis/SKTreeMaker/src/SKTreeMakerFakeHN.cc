
// $Id: SKTreeMakerFakeEl.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQAnalyzer Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "SKTreeMakerFakeHN.h"

//Core includes
#include "EventBase.h"                                                                                                                           


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (SKTreeMakerFakeHN);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
SKTreeMakerFakeHN::SKTreeMakerFakeHN() :  AnalyzerCore(), out_muons(0), out_electrons(0),out_photons(0), out_jets(0), out_fatjets(0), out_genjets(0), out_truth(0), nevents(0),pass_eventcut(0), pass_vertexcut(0) {

  // To have the correct name in the log:                                                                                                                            
  SetLogName("SKTreeMakerFakeHN");
  
  
  
  
}

void SKTreeMakerFakeHN::ExecuteEvents()throw( LQError ){
  
  bool _doubleEG =(k_channel.Contains("DoubleEG"));
  bool _doubleMuon =(k_channel.Contains("DoubleMuon"));
  bool _singleEG =(k_channel.Contains("SingleElectron"));
  bool _singleMuon =(k_channel.Contains("SingleMuon"));

  TString triggerslist_8="HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v";   /// -> tighter cut in lepton ID form tighter trigger emulation cut                                       
  TString triggerslist_12="HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  TString triggerslist_18="HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  TString triggerslist_23="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";

  if(_doubleEG){
    bool passtrig=false;
    if(PassTrigger(triggerslist_8)) passtrig=true;
    if(PassTrigger(triggerslist_12)) passtrig=true;
    if(PassTrigger(triggerslist_18)) passtrig=true;
    if(PassTrigger(triggerslist_23)) passtrig=true;
    if(!passtrig) throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );


    std::vector<snu::KElectron> elColl = GetElectrons("ELECTRON_HN_FAKELOOSEST");
    std::vector<snu::KElectron> elCollpog = GetElectrons("ELECTRON_POG_FAKELOOSEST");


    std::vector<snu::KElectron> elCollTight = GetElectrons("ELECTRON_POG_TIGHT"); 
    std::vector<snu::KElectron> elHNCollTight = GetElectrons("ELECTRON_HN_TIGHT"); 
    if(elCollTight.size() > 1 && elHNCollTight.size() > 1) throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );
    if(elColl.size() < 1 && elCollpog.size() < 1) throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );

    std::vector<snu::KMuon> muonColl = GetMuons("MUON_HN_VETO");  // loose selection                                                                                                                                                               
    if(muonColl.size() > 0) throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );
    
  }
  else if(_doubleMuon){
    bool passtrig=false;
    if(PassTrigger("HLT_Mu3_PFJet40_v")) passtrig=true;
    if(PassTrigger("HLT_Mu8_v")) passtrig=true;
    if(PassTrigger("HLT_Mu17_v")) passtrig=true;
    if(!passtrig) throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );
    std::vector<snu::KElectron> elColl = GetElectrons("ELECTRON_HN_VETO");  // loose selection                                                    
    if(elColl.size() > 0) throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );

    std::vector<snu::KMuon> mucoll = GetMuons("MUON_HN_FAKELOOSEST");
    std::vector<snu::KMuon> mucollTight = GetMuons("MUON_POG_TIGHT");
    if(mucollTight.size() > 1) throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );
    if(mucoll.size() < 1) throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );


  }
  else if(_singleMuon){

    TString analysis_trigger_muon="HLT_IsoMu24_v";
    TString analysis_trigger_tkmuon="HLT_IsoTkMu24_v";
    if(!(PassTrigger(analysis_trigger_muon) || PassTrigger(analysis_trigger_tkmuon)) )throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );
    std::vector<snu::KElectron> elColl = GetElectrons("ELECTRON_HN_VETO");  // loose selection                                              
    
    if(elColl.size() > 0) throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );
    
    std::vector<snu::KMuon> mucoll = GetMuons("MUON_HN_FAKELOOSEST");
    std::vector<snu::KMuon> mucollTight = GetMuons("MUON_POG_TIGHT");
    if(mucollTight.size() > 1) throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );
    if(mucoll.size() < 1) throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );
    if(mucoll.size() >= 1){
      if(mucoll[0].Pt() < 25)  throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );
    }
  }
  else if(_singleEG){

    TString analysis_trigger_1="HLT_Ele27_WPTight_Gsf_v";
    TString analysis_trigger_2="HLT_Ele25_eta2p1_WPTight_Gsf_v";
    if(!(PassTrigger(analysis_trigger_1) || PassTrigger(analysis_trigger_2)) )throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );
    
    std::vector<snu::KElectron> elColl = GetElectrons("ELECTRON_HN_FAKELOOSEST");
    std::vector<snu::KElectron> elCollpog = GetElectrons("ELECTRON_POG_FAKELOOSEST");
    std::vector<snu::KElectron> elCollTight = GetElectrons("ELECTRON_POG_TIGHT");
    if(elCollTight.size() > 1) throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );
    if(elColl.size() < 1 && elCollpog.size() < 1 ) throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );
    if(elColl.size() >= 1){
      if(elColl[0].Pt() < 25)  throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );
    }
    
    std::vector<snu::KMuon> muonColl = GetMuons("MUON_HN_VETO");  // loose selection                                                        
    if(muonColl.size() > 0) throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );

  }
  else {
    bool passtrig=false;
    TString analysis_trigger_1="HLT_Ele27_WPTight_Gsf_v";
    TString analysis_trigger_2="HLT_Ele25_eta2p1_WPTight_Gsf_v";
    TString analysis_trigger_muon="HLT_IsoMu24_v";
    TString analysis_trigger_tkmuon="HLT_IsoTkMu24_v";
    if(PassTrigger(analysis_trigger_muon)) passtrig=true;
    if(PassTrigger(analysis_trigger_tkmuon)) passtrig=true;
    if(PassTrigger(triggerslist_8)) passtrig=true;
    if(PassTrigger(triggerslist_12)) passtrig=true;
    if(PassTrigger(triggerslist_18)) passtrig=true;
    if(PassTrigger(triggerslist_23)) passtrig=true;
    if(PassTrigger("HLT_Mu3_PFJet40_v")) passtrig=true;
    if(PassTrigger("HLT_Mu8_v")) passtrig=true;
    if(PassTrigger("HLT_Mu17_v")) passtrig=true;

    if(PassTrigger(analysis_trigger_1))  passtrig=true;
    if(PassTrigger(analysis_trigger_2))  passtrig=true;

    if(!passtrig) throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );
    
    std::vector<snu::KMuon> muonColl = GetMuons("MUON_HN_VETO");  // loose selection                                                                                                                                                            
    std::vector<snu::KElectron> elColl = GetElectrons("ELECTRON_HN_VETO");  // loose selection                                                                                                                                                  
    if((muonColl.size() > 0) && (elColl.size() > 0)) throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );
   
  }
  
  //////////////////////////////////////////////////////
  //////////// Select objetcs
  //////////////////////////////////////////////////////   
  
  //######   MUON SELECTION ###############                                                                                                                                      
  Message("Selecting Muons", DEBUG);
  std::vector<snu::KMuon> skim_muons;
  /// Apart from eta/pt muons are required to have a global OR tracker track    && be PF                                                                                         
  eventbase->GetMuonSel()->SetPt(5.);
  eventbase->GetMuonSel()->SetEta(3.);
  eventbase->GetMuonSel()->BasicSelection(out_muons, false); /// Muons For SKTree                                                                                                
  SetCorrectedMomentum(out_muons);


  //###### Electron Selection ########                                                                                                                                           
  Message("Selecting electrons", DEBUG);
  std::vector<snu::KElectron> skim_electrons;
  eventbase->GetElectronSel()->SetPt(10.);
  eventbase->GetElectronSel()->SetEta(3.);
  eventbase->GetElectronSel()->BasicSelection(out_electrons);


  //######   MUON SELECTION ###############
  Message("Selecting Muons", DEBUG);

  //###### JET SELECTION  ################
  Message("Selecting jets", DEBUG);
  eventbase->GetJetSel()->SetPt(10);
  eventbase->GetJetSel()->SetEta(5.);
  eventbase->GetJetSel()->BasicSelection(out_jets);
  
  Message("Selecting fat jets", DEBUG);
  eventbase->GetFatJetSel()->SetPt(20);
  eventbase->GetFatJetSel()->SetEta(5.);
  eventbase->GetFatJetSel()->BasicSelection(out_fatjets);

  if(out_jets.size() ==0) throw LQError( "REMOVE TRIGGERED EVENTS for OR",  LQError::SkipEvent );

  out_event   = eventbase->GetEvent();
  out_trigger = eventbase->GetTrigger();
  out_truth   = eventbase->GetTruth();
  
  return;
}// End of execute event loop



void SKTreeMakerFakeHN::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
}


void SKTreeMakerFakeHN::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);

  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  DeclareVariable(out_electrons, "KElectrons", "LQTree");
  DeclareVariable(out_photons, "KPhotons");

  DeclareVariable(out_muons, "KMuons");
  DeclareVariable(out_jets, "KJets");
  DeclareVariable(out_fatjets, "KFatJets");
  DeclareVariable(out_genjets, "KGenJets");
  DeclareVariable(out_trigger, "KTrigger");
  DeclareVariable(out_event, "KEvent");
  DeclareVariable(out_truth, "KTruth");

  //// Set triggers available in sktree
  triggerlist.clear();
  
  
  if(k_isdata){
    AddTriggerToList("HLT_IsoMu1");
    AddTriggerToList("HLT_IsoMu2");
    AddTriggerToList("HLT_IsoTkMu1");
    AddTriggerToList("HLT_IsoTkMu2");
    AddTriggerToList("HLT_Mu8");
    AddTriggerToList("HLT_Mu1");
    AddTriggerToList("HLT_Mu2");
    AddTriggerToList("HLT_Mu3");
    AddTriggerToList("HLT_Mu5");
    AddTriggerToList("HLT_DoubleEle2");
    AddTriggerToList("HLT_DoubleEle8");
    AddTriggerToList("HLT_DoubleEle3");
    AddTriggerToList("HLT_Ele8");
    AddTriggerToList("HLT_Ele1");//// 12-16-18                                                                                                                                                                                                                                                                                                   
    AddTriggerToList("HLT_Ele2");
    AddTriggerToList("HLT_Ele3");
    AddTriggerToList("HLT_DoublePhoton");
    AddTriggerToList("HLT_Photon36_R9Id90");
    AddTriggerToList("HLT_Photon50_R9Id90");
    AddTriggerToList("HLT_Photon90_R9Id90");
    AddTriggerToList("HLT_TripleMu");
    AddTriggerToList("HLT_DiMu");
    
  }
  else {
    AddTriggerToList("HLT_IsoMu1");
    AddTriggerToList("HLT_IsoMu2");
    AddTriggerToList("HLT_IsoTkMu1");
    AddTriggerToList("HLT_IsoTkMu2");
    AddTriggerToList("HLT_Mu8");
    AddTriggerToList("HLT_Mu1");
    AddTriggerToList("HLT_Mu2");
    AddTriggerToList("HLT_Mu3");
    AddTriggerToList("HLT_Mu5");
    AddTriggerToList("HLT_DoubleEle2");
    AddTriggerToList("HLT_DoubleEle8");
    AddTriggerToList("HLT_DoubleEle3");
    AddTriggerToList("HLT_Ele8");
    AddTriggerToList("HLT_Ele1");//// 12-16-18                                                                                                                                                                                                                                                                                                  
    AddTriggerToList("HLT_Ele2");
    AddTriggerToList("HLT_Ele3");
    AddTriggerToList("HLT_DoublePhoton");
    AddTriggerToList("HLT_Photon36_R9Id90");
    AddTriggerToList("HLT_Photon50_R9Id90");
    AddTriggerToList("HLT_Photon90_R9Id90");
    AddTriggerToList("HLT_TripleMu");
    AddTriggerToList("HLT_DiMu");

  }

  
  
    

  return;
  
}

SKTreeMakerFakeHN::~SKTreeMakerFakeHN() {
  
  Message("In Analyzer Destructor" , INFO);

}


void SKTreeMakerFakeHN::FillCutFlow(TString cut, float weight){


  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);

  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 4,0.,4.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"TriLep");
  }
}



void SKTreeMakerFakeHN::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void SKTreeMakerFakeHN::ClearOutputVectors() throw (LQError){
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
  out_photons.clear();
  out_jets.clear();
  out_fatjets.clear();
  out_genjets.clear();
  out_truth.clear();

}



