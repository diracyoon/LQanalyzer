// $Id: SKTreeValidation.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQSKTreeValidation Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "SKTreeValidation.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (SKTreeValidation);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
SKTreeValidation::SKTreeValidation() :  AnalyzerCore(), out_muons(0)  {
  
  
  mapcounter.clear();
  
  SetLogName("SKTreeValidation");
  
  Message("In SKTreeValidation constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

  if(1){
    std::vector<TString> dimuonIDs;
    //dimuonIDs.push_back("POGTight");
    dimuonIDs.push_back("POGTightroch_corrected");
    //    dimuonIDs.push_back("POGTightNoJetSmearroch_corrected");
    //dimuonIDs.push_back("HNTightroch_corrected");
    //dimuonIDs.push_back("POGTightroch_correctedtruthmatch");
    for(unsigned int i=0; i < dimuonIDs.size(); i++){
      MakeCleverHistograms(sighist_mm,"ZMuon"+dimuonIDs.at(i));
      MakeCleverHistograms(sighist_mm,"DiMuon"+dimuonIDs.at(i));
      MakeCleverHistograms(sighist_mm,"DiMuon_puW"+dimuonIDs.at(i));
      MakeCleverHistograms(sighist_mm, "DiMuon_dijet"+dimuonIDs.at(i));
      MakeCleverHistograms(sighist_mm, "DiMuon_01jet"+dimuonIDs.at(i));
      MakeCleverHistograms(sighist_mm, "DiMuon_2jet"+dimuonIDs.at(i));
      MakeCleverHistograms(sighist_mm, "DiMuon_dibjet"+dimuonIDs.at(i));
      MakeCleverHistograms(sighist_mm, "DiMuon_SSPreselection"+dimuonIDs.at(i));
    }
    MakeCleverHistograms(sighist_mm,"DiMuon");
    MakeCleverHistograms(sighist_mm,"DiJet");
    
    MakeCleverHistograms(muhist,"KMuonHists");
    MakeCleverHistograms(muhist,"KMuonHists_POGTIGHT");
    
    /*
      std::vector<TString> electronIDs;
    electronIDs.push_back("POGTight");
    electronIDs.push_back("POGTighttruthmatch");
    for(unsigned int i=0; i < electronIDs.size(); i++){
      MakeCleverHistograms(sighist_ee,"SingleElectron"+electronIDs.at(i));
      MakeCleverHistograms(sighist_ee,"SingleElectron_noW"+electronIDs.at(i));
      MakeCleverHistograms(sighist_ee,"SingleElectron_IDW"+electronIDs.at(i));
      MakeCleverHistograms(sighist_ee,"SingleElectron_puW"+electronIDs.at(i));
      MakeCleverHistograms(sighist_ee,"SingleElectron_Trigger"+electronIDs.at(i));
      MakeCleverHistograms(sighist_ee,"SingleElectron_Wregion"+electronIDs.at(i));
      MakeCleverHistograms(sighist_ee, "SingleElectron_dijet"+electronIDs.at(i));
    }
    */
    std::vector<TString> dielectronIDs;
    dielectronIDs.push_back("POGTight");
    //dielectronIDs.push_back("POGTighttruthmatch");
    
    for(unsigned int i=0; i < dielectronIDs.size(); i++){
      MakeCleverHistograms(sighist_ee,"ZElectron"+dielectronIDs.at(i));
      MakeCleverHistograms(sighist_ee,"DiElectron_SSPreselection"+dielectronIDs.at(i));
      MakeCleverHistograms(sighist_ee,"DiElectron_SS1Jet"+dielectronIDs.at(i));
      MakeCleverHistograms(sighist_ee,"DiElectron_SSMet"+dielectronIDs.at(i));
      MakeCleverHistograms(sighist_ee,"DiElectron"+dielectronIDs.at(i));
      MakeCleverHistograms(sighist_ee,"DiElectron_puW"+dielectronIDs.at(i));
      MakeCleverHistograms(sighist_ee,"DiElectron_Trigger"+dielectronIDs.at(i));
      MakeCleverHistograms(sighist_ee, "DiElectron_dijet"+dielectronIDs.at(i));


    }
    MakeCleverHistograms(elhist,"KElectronHists");
    MakeCleverHistograms(elhist,"KElectronHists_POGTIGHT");
    
    std::vector<TString> muonIDs;
    muonIDs.push_back("POGTightroch_corrected");
    for(unsigned int i=0; i < muonIDs.size(); i++){
      MakeCleverHistograms(sighist_mm,"SingleMuon"+muonIDs.at(i));
      MakeCleverHistograms(sighist_mm,"SingleMuon_noW"+muonIDs.at(i));
      MakeCleverHistograms(sighist_mm,"SingleMuon_IDW"+muonIDs.at(i));
      MakeCleverHistograms(sighist_mm,"SingleMuon_puW"+muonIDs.at(i));
      MakeCleverHistograms(sighist_mm,"SingleMuon_Trigger"+muonIDs.at(i));
      MakeCleverHistograms(sighist_mm,"SingleMuon_Wregion"+muonIDs.at(i));
      MakeCleverHistograms(sighist_mm, "SingleMuon_dijet"+muonIDs.at(i));
      
      MakeCleverHistograms(sighist_em,"EMuon"+muonIDs.at(i));
      MakeCleverHistograms(sighist_em,"EMuon_noW"+muonIDs.at(i));
      MakeCleverHistograms(sighist_em,"EMuon_IDW"+muonIDs.at(i));
      MakeCleverHistograms(sighist_em,"EMuon_puW"+muonIDs.at(i));
      MakeCleverHistograms(sighist_em,"EMuon_Trigger"+muonIDs.at(i));
      MakeCleverHistograms(sighist_em, "EMuon_dijet"+muonIDs.at(i));
      MakeCleverHistograms(sighist_em,"EMuon_SSPreselection"+muonIDs.at(i));
      
    }
  }

}



void SKTreeValidation::FillCutFlow(TString cut, float weight){

  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
   
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 10,0.,10.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"LumiMaskCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"MuMuCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(7,"OSMuMuCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(8,"OSMuMuJJCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(9,"SSMuMuCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(10,"SSMuMuJJCut");   
    
  }
}

void SKTreeValidation::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
  
  Message("Making clever hists for Z ->ll test code", INFO);
  

  return;
}


void SKTreeValidation::ExecuteEvents()throw( LQError ){

  int mcperiod = GetMCPeriod();


  if(!isData) weight*=MCweight;
  //ListTriggersAvailable();
  counter("NoCut",weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));

  FillCutFlow("NoCut",  weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));
  if(TString(getenv("CATVERSION")).Contains("v7")){
    if(isData&& (! eventbase->GetEvent().LumiMask())) return;
    counter("LumiMaskCut", weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));

    FillCutFlow("LumiMaskCut", weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));
  }


  /// Apply the gen weight
  /// Acts on data to remove bad reconstructed event 

  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   

  bool makePUFile=false;
  if(makePUFile){
    //std::vector<snu::KMuon> muons = GetMuons(BaseSelection::MUON_POG_TIGHT");
    std::vector<snu::KMuon> muons = GetMuons("MUON_HN_TRI_TIGHT");
    return;
    
    //if(TriggerEff("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", muons)){
    if(muons.size() ==2) {
      if(!SameCharge(muons)){
	if(muons.at(0).Pt() > 20. && muons.at(1).Pt() > 20.){
	  
	  if(GetDiLepMass(muons) < 120. && GetDiLepMass(muons)  > 60. ){
	      
	    if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 60., 60, "N_{vertex}");
	    else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 60., 60, "N_{vertex}");
	  }
	}
      }
    }
    return;
  }

  ///#### CAT:::PassBasicEventCuts is updated: uses selections as described in https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters: If you <b>see this is out of date please comment

  
  
  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut",  weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));;
  counter("EventCut", weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));

   /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton

  TString dimuon_trigmuon_trig1="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v";
  TString dimuon_trigmuon_trig2="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v";

  
  TString muon_trigmuon_trig1="HLT_IsoMu24_v";
  TString muon_trigmuon_trig2="HLT_IsoTkMu24_v";
  
  TString diel_trig="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
  
  TString em1_trig="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";
  
  TString el_trig="HLT_Ele27_WPTight_Gsf_v";
  if(Is2015Analysis()){

    el_trig="HLT_Ele23_WPLoose_Gsf_v";
    muon_trigmuon_trig1="HLT_IsoMu20_v";
    muon_trigmuon_trig2="HLT_IsoTkMu20_v";
    diel_trig="HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
    em1_trig="HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v";
  }

  if(IsDiEl()){
    vector <snu::KElectron> vecel = GetElectrons("ELECTRON_POG_TIGHT");
    if(vecel.size() > 1) counter("Diel_pass",weight);
    else counter("Diel_fail",weight);
    if(SameCharge(vecel)) counter("SSDiEL",weight);
  }
  std::vector<TString> triggerslist_dimu;
  triggerslist_dimu.push_back(dimuon_trigmuon_trig1);
  //triggerslist_dimu.push_back(dimuon_trigmuon_trig2);
  
  std::vector<TString> triggerslist_mu;
  triggerslist_mu.push_back(muon_trigmuon_trig1);
  triggerslist_mu.push_back(muon_trigmuon_trig2);
  
  std::vector<TString> triggerslist_diel;
  triggerslist_diel.push_back(diel_trig);
  
  std::vector<TString> triggerslist_el;
  triggerslist_el.push_back(el_trig);
  
  std::vector<TString> triggerslist_emu;
   triggerslist_emu.push_back(em1_trig);
   
   
   FillCutFlow("TriggerCut",  weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));;
   counter("TriggerCut", weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));

   if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex

   FillCutFlow("VertexCut",  weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));;
   counter("VertexCut", weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));

   
   std::vector<snu::KMuon>  mus = GetMuons("MUON_POG_TIGHT");
   
   if(mus.size() ==2) {
	   
     counter("DiMu", weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));
     //     m_logger << INFO << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;

     if(!SameCharge(mus)){
       if(mus.at(0).Pt() > 20. && mus.at(1).Pt() > 20.){

	 //for(int i=0; i < GetJets("JET_HN").size() ; i++){
	 //	   cout << "Jet pt = " << GetJets("JET_HN").at(i).Pt() << " eta = " << GetJets("JET_HN").at(i).Eta() << endl;
	 //	 }
	 /*for(int i=0; i < GetMuons("MUON_HN_VETO").size() ; i++){
	   cout << "Muon pt = " << GetMuons("MUON_HN_VETO").at(i).Pt() << " eta = " << GetMuons("MUON_HN_VETO").at(i).Eta()  << endl;
	 }
	 for(int i=0; i < GetElectrons("ELECTRON_HN_VETO").size() ; i++){
           cout << "Electron pt = " << GetElectrons("ELECTRON_HN_VETO").at(i).Pt() << " eta = " << GetElectrons("ELECTRON_HN_VETO").at(i).Eta()  <<endl;
	   }*/

         FillCutFlow("OSMuMuCut",  weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));
	 FillCLHist(sighist_mm,"DiMuon",  eventbase->GetEvent(), GetMuons("MUON_POG_TIGHT"), GetElectrons("ELECTRON_POG_TIGHT"), GetJets("JET_HN"), weight);
	 if(GetJets("JET_HN").size() > 1){
	   FillCLHist(sighist_mm,"DiJet",  eventbase->GetEvent(), GetMuons("MUON_POG_TIGHT"), GetElectrons("ELECTRON_POG_TIGHT"), GetJets("JET_HN"), weight);
	 }
         if(GetJets("JET_HN").size() ==0) counter("0Jet", weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));
	 if(GetJets("JET_HN").size() == 1) counter("1Jet", weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));
	 if(GetJets("JET_HN").size()== 2) counter("2Jet", weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));
	 if(GetJets("JET_HN").size() == 3) counter("3Jet", weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));
	 if(GetJets("JET_HN").size() == 4) counter("4Jet", weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));
	 if(GetJets("JET_HN").size()> 4)   counter("5Jet",  weight*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));
       }
     }
   }



   float pileup_reweight=(1.0);

   if (!k_isdata) {
     pileup_reweight = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());
     //pileup_reweight= eventbase->GetEvent().PileUpWeight();
     FillHist("PUWeightvsNVertex",pileup_reweight, eventbase->GetEvent().nVertices(), weight, 0., 5., 500, 0., 60., 60, "N_{vertex}");
   }
   
   FillHist("PileupWeight" ,  pileup_reweight,weight,  0. , 50., 10, "pu weight");

   ///_______________________________________________________________________________________________________________________________________________________________________________________________________//
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   //// Electrons plots
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   MakeElectronValidationPlots("ELECTRON_POG_TIGHT", weight, pileup_reweight,triggerslist_el, "MUON_POG_TIGHT","JET_HN", "POGTight");
   //   MakeElectronValidationPlots("ELECTRON_POG_TIGHT, weight, pileup_reweight,triggerslist_el, "MUON_POG_TIGHT,"JET_HN", "POGTighttruthmatch");
   FillCLHist(elhist,"KElectronHists", GetElectrons("ELECTRON_NOCUT"), weight);
   FillCLHist(elhist,"KElectronHists_POGTIGHT", GetElectrons("ELECTRON_POG_TIGHT"), weight);
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ///// DiElectron Validation plots
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   MakeDiElectronValidationPlots("ELECTRON16_HN_TIGHT", weight, pileup_reweight,triggerslist_diel, "MUON_POG_TIGHT","JET_HN", "POGTight");
   //MakeDiElectronValidationPlots("ELECTRON_POG_TIGHT", weight, pileup_reweight, triggerslist_diel, "MUON_POG_TIGHT","JET_HN", "POGTighttruthmatch");
   ///_______________________________________________________________________________________________________________________________________________________________________________________________________//
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   /// Single Muon Validation plots
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   //MakeMuonValidationPlots("MUON_POG_TIGHT", weight, pileup_reweight,triggerslist_mu, "ELECTRON_POG_TIGHT","JET_HN", "POGTightroch_corrected");
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   /// DiMuon Validation plots
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   //MakeDiMuonValidationPlots("MUON_POG_TIGHT", weight, pileup_reweight,triggerslist_dimu, "ELECTRON_POG_TIGHT","JET_HN", "POGTight");
   //MakeDiMuonValidationPlots("MUON_POG_TIGHT", weight, pileup_reweight, triggerslist_dimu, "ELECTRON_POG_TIGHT","JET_HN", "POGTightNoJetSmearroch_corrected",false);
   MakeDiMuonValidationPlots("MUON_POG_TIGHT", weight, pileup_reweight, triggerslist_dimu, "ELECTRON_POG_TIGHT","JET_HN", "POGTightroch_corrected", true, mcperiod);
   //MakeDiMuonValidationPlots("MUON_POG_TIGHT", weight, pileup_reweight,triggerslist_dimu, "ELECTRON_POG_TIGHT","JET_HN", "POGTightroch_correctedtruthmatch");
   //MakeDiMuonValidationPlots("MUON_HN_TIGHT", weight, pileup_reweight, triggerslist_dimu, "ELECTRON_POG_TIGHT","JET_HN", "HNTightroch_corrected");
   FillCLHist(muhist,"KMuonHists", GetMuons("MUON_NOCUT"), weight);
   FillCLHist(muhist,"KMuonHists_POGTIGHT", GetMuons("MUON_POG_TIGHT"), weight);
   ///_______________________________________________________________________________________________________________________________________________________________________________________________________//
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   // ElectronMuon Validation
   MakeElMuonValidationPlots("MUON_POG_TIGHT", weight, pileup_reweight,triggerslist_emu, "ELECTRON_POG_TIGHT","JET_HN", "POGTightroch_corrected");
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////!!!!!!////////
   ///_______________________________________________________________________________________________________________________________________________________________________________________________________//


   
   return;
}// End of execute event loop
  
void SKTreeValidation::counter(TString cut, float w){
  w=1.;
  map<TString,float>::iterator itmapcounter = mapcounter.find(cut) ;
  if (itmapcounter == mapcounter.end()){
    mapcounter[cut] = w;
  }
  else{
    float sum = itmapcounter->second;
    mapcounter[cut] = sum+w;
  }

}

void SKTreeValidation::EndCycle()throw( LQError ){
  
  for ( map<TString,float>::iterator itmapcounter = mapcounter.begin(); itmapcounter != mapcounter.end(); itmapcounter++){
    cout << itmapcounter->first << " has count "<< itmapcounter->second << endl;
  }
  Message("In EndCycle" , INFO);

}


void SKTreeValidation::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  

  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //  DeclareVariable(out_muons, "Signal_Muons");

  
  return;
  
}

SKTreeValidation::~SKTreeValidation() {
  
  Message("In SKTreeValidation Destructor" , INFO);
  
}



void SKTreeValidation::MakeMuonValidationPlots(TString muid, float w, float pu_reweight, std::vector<TString> trignames,TString elid, TString jetid, TString tag){

  Message("In MakeMuonValidationPlots " , DEBUG);

  std::vector<snu::KElectron> electrons =  GetElectrons(elid);
  std::vector<snu::KJet> jets =  GetJets(jetid);


  std::vector<snu::KMuon> muons;

  if(k_running_nonprompt){
    muons = GetMuons("MUON_HN_LOOSE",false);
  }
  else if(tag.Contains("truthmatch"))   muons = GetMuons(muid,false);
  else   muons = GetMuons(muid);

  if(tag.Contains("roch"))   CorrectMuonMomentum(muons);

  bool trig_pass= PassTrigger(trignames.at(0)) ;
  if(trignames.size() ==2) trig_pass = trig_pass || PassTrigger(trignames.at(1)) ;


  /// List of all corrections to be applied
  float trigger_sf(1.);
  float id_iso_sf(1.);
  float trigger_ps(1.);

  float ev_weight(1.);

  if(!isData){
    trigger_sf =   mcdata_correction->TriggerScaleFactor(electrons,muons, trignames.at(0));
    id_iso_sf =   mcdata_correction->MuonScaleFactor(muid, muons,0);
    id_iso_sf *=  mcdata_correction->MuonISOScaleFactor(muid, muons,0);

    /// Apply weight to MC to scale to lumniosity
    trigger_ps= WeightByTrigger(trignames, TargetLumi)  ;

    ev_weight = w * trigger_sf * id_iso_sf * pu_reweight*trigger_ps;
  }
  if(k_running_nonprompt){
    ev_weight=1.; /// In case... should not be needed
    ev_weight      *=   m_datadriven_bkg->Get_DataDrivenWeight_M(false,muons,muid, "dijet");

  }

  if (!trig_pass) return;


  if(muons.size() ==1){
    if(muons.at(0).Pt() > 25. ){
      FillCLHist(sighist_mm, "SingleMuon"+tag, eventbase->GetEvent(), muons,electrons,jets, weight);
      if(trig_pass){
	FillHist("nvertex_mu" + tag + "_nopurw", eventbase->GetEvent().nVertices(),  w*trigger_sf*trigger_ps, 0., 40.,40, "N_{vertex}");
        FillCLHist(sighist_mm, "SingleMuon_Trigger"+tag, eventbase->GetEvent(), muons,electrons,jets, w*trigger_sf*trigger_ps);
	FillHist("nvertex_mu" + tag + "_purw", eventbase->GetEvent().nVertices()  ,  w*trigger_sf*trigger_ps*pu_reweight, 0., 40.,40, "N_{vertex}") ;
	
	FillCLHist(sighist_mm, "SingleMuon_puW"+tag, eventbase->GetEvent(), muons,electrons,jets, w*pu_reweight*trigger_ps*trigger_sf);

	FillCLHist(sighist_mm, "SingleMuon_IDW"+tag, eventbase->GetEvent(), muons,electrons,jets, w*pu_reweight*id_iso_sf*trigger_ps*trigger_sf);


	float METdphi = TVector2::Phi_mpi_pi(muons.at(0).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet));
	float MT=(2.* muons.at(0).Et()*eventbase->GetEvent().MET(snu::KEvent::pfmet) * (1 - cos( METdphi)));
	
	if(MT > 40. && eventbase->GetEvent().MET(snu::KEvent::pfmet) > 30.) FillCLHist(sighist_mm, "SingleMuon_Wregion"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
	  
	if(jets.size() >= 2)          FillCLHist(sighist_mm, "SingleMuon_dijet"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
      }
    }
  }
}


void SKTreeValidation::MakeDiMuonValidationPlots(TString muid, float w, float pu_reweight,  std::vector<TString> trignames,TString elid, TString jetid, TString tag, bool smearjets, int mcperiod){

  Message("In MakeDiMuonValidationPlots " , DEBUG);
  

  std::vector<snu::KElectron> electrons =  GetElectrons(elid);
  std::vector<snu::KJet> jets =  GetJets(jetid, smearjets);

  int nbjet(0);
  for(unsigned int ij =0 ; ij < jets.size() ; ij++){
    if(IsBTagged(jets.at(ij),snu::KJet::CSVv2, snu::KJet::Medium, mcperiod)) nbjet++;
  }
  //  if(nbjet == 1) w*= 0.8;
  //  if(nbjet == 2) w*= 0.65;

  if(!isData) w =  w*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",TargetLumi);
  
  
   
  std::vector<snu::KMuon> muons;
  if(k_running_nonprompt){
    muons = GetMuons("MUON_HN_LOOSE",false);
  }
  else if(tag.Contains("truthmatch"))   muons = GetMuons(muid,false);
  else   muons = GetMuons(muid,false);

  FillHist("zpeak_mumu_test" + tag, GetDiLepMass(muons), w, 0., 200.,400);
  
  if(muons.size() != 2) return;
  
  if(tag.Contains("roch"))   CorrectMuonMomentum(muons);
  
  // returns 1. if event passes trigger bit in data or Eff of trigger Eff(pt_muon) if in MC
  //float trig_pass= TriggerEff("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", muons);
  /// returns bool if event passes trigger bit
  bool trig_pass= PassTrigger(trignames.at(0)) ;
  if(trignames.size() ==2) trig_pass = trig_pass || PassTrigger(trignames.at(1)) ;
  
  /// List of all corrections to be applied
  float trigger_sf(1.);
  float id_iso_sf(1.);
  float trigger_ps(1.);
  float ev_weight(1.);
  float puweight(1.);

  if(!isData){
    //trigger_sf = TriggerScaleFactor(electrons,muons, trignames.at(0));
    id_iso_sf =   mcdata_correction->MuonScaleFactor(muid, muons,0);
    id_iso_sf *=  mcdata_correction->MuonISOScaleFactor(muid, muons,0);
    //    id_iso_sf *= MuonTrackingEffScaleFactor(muons);
    //cout << "muon pt/eta = " << muons.at(0).Pt() << " / " << muons.at(0).Eta() << " id scale factor = " <<  MuonScaleFactor(muid, muons,0) << " iso scale factor =" << MuonISOScaleFactor(muid, muons,0) << endl;
 
    // uses cattool 69 mb weight
    //puweight=eventbase->GetEvent().PileUpWeight();
    /// uses own weight
    puweight=mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());
    
    /// scale to lumi of trigger
    //trigger_ps= WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");;
    
    ev_weight = w * id_iso_sf *  pu_reweight;
  }
  //cout << trigger_ps << " " << w << " " << trigger_sf << " " << id_iso_sf << " " << pu_reweight << " " << trig_pass << endl;
  if(k_running_nonprompt){
    ev_weight=1.; /// In case... should not be needed
    //ev_weight      *=  m_datadriven_bkg->Get_DataDrivenWeight_MM(false,muons,muid, "dijet");


  }

  //if (trig_pass > 0.)  w*= trig_pass;
  // cout << "jets.size = " << jets.size() << endl;

  if(muons.size() ==2) {
    if(!SameCharge(muons)){
      if(muons.at(0).Pt() > 20. && muons.at(1).Pt() > 20.){
	FillCutFlow("MuMuCut",  w);
	
	/// Z peak plots
	if(GetDiLepMass(muons) < 120. && GetDiLepMass(muons)  > 60. ){
	  FillHist("zpeak_mumu_notrig" + tag + "_nopurw", GetDiLepMass(muons), w, 0., 200.,400);
	  if(trig_pass){
	    FillHist("zpeak_mumu" + tag + "_nopurw", GetDiLepMass(muons), w*id_iso_sf, 0., 200.,400, "m_Z [GeV])");
	    FillHist("zpeak_mumu" + tag + "_purw", GetDiLepMass(muons),    w*id_iso_sf*pu_reweight, 0., 200.,400, "m_Z [GeV])");
	    FillHist("nvertex_mumu" + tag + "_nopurw", eventbase->GetEvent().nVertices(),  w*id_iso_sf, 0., 40.,40, "n_{vertex}");
	    FillHist("nvertex_mumu" + tag + "_purw", eventbase->GetEvent().nVertices()  , w*id_iso_sf*pu_reweight, 0., 40.,40, "n_{vertex}") ;
	    FillCLHist(sighist_mm, "ZMuon"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
	  }
	}

	if(trig_pass){
	  FillCLHist(sighist_mm, "DiMuon"+tag, eventbase->GetEvent(), muons,electrons,jets, w);
	  FillCLHist(sighist_mm, "DiMuon_puW"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
	  
	  FillHistPerLumi("DiMuon_llmass"+tag,GetDiLepMass(muons), ev_weight, 0., 200.,20, 10);
	  FillHistPerLumi("DiMuon_njet"+tag,jets.size(), ev_weight, 0., 10.,10, 10);
	  FillHistPerLumi("DiMuon_met"+tag,eventbase->GetEvent().MET(snu::KEvent::pfmet), ev_weight, 0., 200.,20, 10);
	  
	  
	  if(jets.size() >= 2) {
	    FillCutFlow("OSMuMuJJCut",  w);
	    FillCLHist(sighist_mm, "DiMuon_dijet"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
	    if(nbjet > 1 ) FillCutFlow("DiMuon_dibjet",  w*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi));
	    if(nbjet > 1 && eventbase->GetEvent().MET(snu::KEvent::pfmet) > 50.) FillCutFlow("DiMuon_met",  w*WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi)); 
	    if(nbjet > 1 && eventbase->GetEvent().MET(snu::KEvent::pfmet) > 50.)             FillCLHist(sighist_mm, "DiMuon_dibjet"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
	    if(muons.at(0).Pt()  < 50.) FillCLHist(sighist_mm, "DiMuon_2jet"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
	  }
	  else{
	    if(muons.at(0).Pt() > 50)FillCLHist(sighist_mm, "DiMuon_01jet"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
	    
	  }
	}
      }
    }
    else{
      if(muons.at(0).Pt() > 20. && muons.at(1).Pt() > 20.){

        /// Z peak plots
        if(GetDiLepMass(muons)  > 15. ){
          if(trig_pass){
	    FillCutFlow("SSMuMuCut",  w);
	    if(jets.size() >= 2){
	      FillCLHist(sighist_mm, "DiMuon_SSPreselection"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
	      FillCutFlow("SSMuMuJJCut",  w);
	    }
	    
	  }
	}
      }
    }
  }
}

void SKTreeValidation::MakeElMuonValidationPlots(TString muid, float w, float pu_reweight,  std::vector<TString> trignames,TString elid, TString jetid, TString tag){

  Message("In MakeElMuonValidationPlots " , DEBUG);

  std::vector<snu::KElectron> electrons ;
  std::vector<snu::KMuon> muons;
  if(k_running_nonprompt){
    muons = GetMuons("MUON_HN_LOOSE",false);
    electrons             = GetElectrons(true, false,"ELECTRON_HN_FAKELOOSE_NOD0");
  }
  else if(tag.Contains("truthmatch")) {
    muons = GetMuons(muid,false);
    electrons = GetElectrons(true, false,elid);
  }
  else {
    muons = GetMuons(muid);
    electrons = GetElectrons(elid);
  }

  if(tag.Contains("roch")){
    CorrectMuonMomentum(muons);
  }

  std::vector<snu::KJet> jets =  GetJets(jetid);

  /// List of all corrections to be applied
  float trigger_sf(1.);
  float id_iso_sf(1.);
  float trigger_ps(1.);
  float ev_weight(1.);

  if(!isData){
    //trigger_sf = TriggerScaleFactor(electrons,muons, trignames.at(0));
    id_iso_sf =  mcdata_correction->MuonScaleFactor(muid, muons,0);
    id_iso_sf *= mcdata_correction->MuonISOScaleFactor(muid, muons,0);
    id_iso_sf *= mcdata_correction->ElectronScaleFactor(elid, electrons,0); ///MUON_POG_TIGHT == MUON_HN_TIGHT
    id_iso_sf *= mcdata_correction->ElectronRecoScaleFactor(electrons);


    /// Tiny effect on unprescaled triggers
    trigger_ps= WeightByTrigger(trignames, TargetLumi)  ;
    ev_weight = w * trigger_sf * id_iso_sf * trigger_ps* pu_reweight;
  }

  bool trig_pass= PassTrigger(trignames.at(0)) ;
  if(trignames.size() ==2) trig_pass = trig_pass || PassTrigger(trignames.at(1)) ;

  if(k_running_nonprompt){
    ev_weight=1.; /// In case... should not be needed
    ev_weight      *=  m_datadriven_bkg->Get_DataDrivenWeight_EM(false,muons,electrons,muid,elid, "dijet");
  }


  if(muons.size() ==1 && electrons.size() ==1) {
    if(muons.at(0).Charge() != electrons.at(0).Charge()){
      if(muons.at(0).Pt() > 20. && electrons.at(0).Pt() > 20.){

	if(trig_pass){
	  FillCLHist(sighist_em, "EMuon_Trigger"+tag, eventbase->GetEvent(), muons,electrons,jets, w*trigger_sf*trigger_ps);
	  snu::KParticle Z = muons.at(0) + electrons.at(0);
          if(Z.M() > 15.){
	    
	    FillCLHist(sighist_em, "EMuon"+tag, eventbase->GetEvent(), muons,electrons,jets, w*trigger_sf*trigger_ps);
	    FillCLHist(sighist_em, "EMuon_puW"+tag, eventbase->GetEvent(), muons,electrons,jets, w*pu_reweight*trigger_sf*trigger_ps);
	    FillCLHist(sighist_em, "EMuon_IDW"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
	  }
	}
      }
    }
    else{
      if(muons.at(0).Pt() > 20. && electrons.at(0).Pt() > 20.){
	if(trig_pass){
          snu::KParticle Z = muons.at(0) + electrons.at(0);
          if(Z.M() > 15.){
	    if(jets.size() >= 2)          FillCLHist(sighist_em, "EMuon_SSPreselection"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
	  }
	}
      }
    }
  }

  }
  


void SKTreeValidation::MakeElectronValidationPlots(TString elid, float w, float pu_reweight,  std::vector<TString> trignames,TString muid, TString jetid, TString tag){

  Message("In MakeElectronValidationPlots " , DEBUG);

  std::vector<snu::KElectron> electrons ;
  if(k_running_nonprompt){
    electrons             = GetElectrons(true, false,"ELECTRON_HN_FAKELOOSE_NOD0");
  }
  else if(tag.Contains("truthmatch"))   electrons = GetElectrons(true, false,elid);
  else   electrons =  GetElectrons(elid);


  std::vector<snu::KElectron> electronVetoColl   = GetElectrons("ELECTRON_HN_VETO");

  std::vector<snu::KMuon> muons =  GetMuons(muid);
  std::vector<snu::KJet> jets =  GetJets(jetid);


  bool trig_pass= PassTrigger(trignames.at(0)) ;
  if(trignames.size() ==2) trig_pass = trig_pass || PassTrigger(trignames.at(1)) ;
  
  /// List of all corrections to be applied
  float trigger_sf(1.);
  float id_iso_sf(1.);
  float trigger_ps(1.);
  float reco_weight=1.;
  
  float ev_weight(1.); 
  if(!isData){
    //trigger_sf = TriggerScaleFactor(electrons,muons, trignames.at(0));
    id_iso_sf=   mcdata_correction->ElectronScaleFactor(elid, electrons,0); ///MUON_POG_TIGHT == MUON_HN_TIGHT
    reco_weight = mcdata_correction->ElectronRecoScaleFactor(electrons);
    /// Tiny effect on unprescaled triggers
    trigger_ps= WeightByTrigger(trignames, TargetLumi)  ;
    ev_weight = w * trigger_sf * id_iso_sf * reco_weight * pu_reweight*trigger_ps;
  }
  if(k_running_nonprompt){
    ev_weight=1.; /// In case... should not be needed
    ev_weight      *=  m_datadriven_bkg->Get_DataDrivenWeight_E(false,electrons,elid, "dijet");

  }


  if(electrons.size() ==1){
    if(electrons.at(0).Pt() > 30. ){
      FillCLHist(sighist_ee, "SingleElectron"+tag, eventbase->GetEvent(), muons,electrons,jets, weight);
      if(trig_pass){
	FillCLHist(sighist_ee, "SingleElectron_Trigger"+tag, eventbase->GetEvent(), muons,electrons,jets, w*trigger_ps*trigger_sf);
	FillHist("nvertex_el" + tag + "_nopurw", eventbase->GetEvent().nVertices(),  w*trigger_ps*trigger_sf, 0., 40.,40);
	FillHist("nvertex_el" + tag + "_purw", eventbase->GetEvent().nVertices()  , w*trigger_ps*trigger_sf*pu_reweight, 0., 40.,40) ;
	
	FillCLHist(sighist_ee, "SingleElectron_puW"+tag, eventbase->GetEvent(), muons,electrons,jets, w*pu_reweight*trigger_ps*trigger_sf);
	
	FillCLHist(sighist_ee, "SingleElectron_IDW"+tag, eventbase->GetEvent(), muons,electrons,jets, w*pu_reweight*trigger_ps*trigger_sf*id_iso_sf);
	
	float METdphi = TVector2::Phi_mpi_pi(electrons.at(0).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet));
        float MT=(2.* electrons.at(0).Et()*eventbase->GetEvent().MET(snu::KEvent::pfmet) * (1 - cos( METdphi)));
	if(MT > 40. && eventbase->GetEvent().MET(snu::KEvent::pfmet) > 30.) FillCLHist(sighist_ee, "SingleElectron_Wregion"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
        if(jets.size() >= 2)          FillCLHist(sighist_ee, "SingleElectron_dijet"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
	
      }
    }
  }
}


void SKTreeValidation::MakeDiElectronValidationPlots(TString elid, float w, float pu_reweight,  std::vector<TString> trignames,TString muid, TString jetid, TString tag){
  counter("MakeDiElectronValidationPlots",w);
  Message("In MakeDiElectronValidationPlots" , DEBUG);

  std::vector<snu::KElectron> electrons ;
  if(k_running_nonprompt){
    electrons             = GetElectrons(true, false,"ELECTRON16_HN_FAKELOOSE_NOD0");
  }
  else if(tag.Contains("truthmatch"))   electrons = GetElectrons(true, false,elid);
  else   electrons =  GetElectrons(elid);
  
  std::vector<snu::KElectron> electronVetoColl   = GetElectrons("ELECTRON_HN_VETO");
  
  std::vector<snu::KMuon> muons =  GetMuons(muid);
  std::vector<snu::KJet> jets =  GetJets(jetid);
  
  bool trig_pass= PassTrigger(trignames.at(0)) ;
  if(trignames.size() ==2) trig_pass = trig_pass || PassTrigger(trignames.at(1)) ;

  float trigger_sf(1.);
  float id_iso_sf(1.);
  float trigger_ps(1.);
  float reco_weight=1.;

  float ev_weight(1.);
  if(!isData){
    //trigger_sf = TriggerScaleFactor(electrons,muons, trignames.at(0)); /// 
    id_iso_sf=   mcdata_correction->ElectronScaleFactor(elid, electrons,0); ///MUON_POG_TIGHT == MUON_HN_TIGHT
    //reco_weight = ElectronRecoScaleFactor(electrons);
    /// Tiny effect on unprescaled triggers
    trigger_ps= WeightByTrigger(trignames.at(0), TargetLumi)  ;
    ev_weight = w * trigger_sf * id_iso_sf * reco_weight * pu_reweight*trigger_ps;
  }
  if(k_running_nonprompt){
    ev_weight=1.; /// In case... should not be needed
    ev_weight      *=  m_datadriven_bkg->Get_DataDrivenWeight_EE(false, electrons);
  }
  weight              *=  m_datadriven_bkg->WeightCFEvent(electrons, k_running_chargeflip);    
  if(!trig_pass) return;
  if( m_datadriven_bkg->WeightCFEvent(electrons, k_running_chargeflip) == 0.) return; 
  
  if(electrons.size() ==2) {
    counter("DiEl",w);

    if(!SameCharge(electrons)){
      if(electrons.at(0).Pt() > 20. && electrons.at(1).Pt() > 20.){

        /// Z peak plots
        if(GetDiLepMass(electrons) < 120. && GetDiLepMass(electrons)  > 60. ){
          if(trig_pass){
            FillHist("zpeak_ee" + tag + "_nopurw", GetDiLepMass(electrons), w*id_iso_sf*trigger_sf, 0., 200.,400);
            FillHist("zpeak_ee" + tag + "_purw", GetDiLepMass(electrons),    w*id_iso_sf*trigger_sf*pu_reweight, 0., 200.,400);
            FillHist("nvertex_ee" + tag + "_nopurw", eventbase->GetEvent().nVertices(),  w*id_iso_sf*trigger_sf, 0., 40.,40);
            FillHist("nvertex_ee" + tag + "_purw", eventbase->GetEvent().nVertices()  , w*id_iso_sf*trigger_sf*pu_reweight, 0., 40.,40) ;
            FillCLHist(sighist_ee, "ZElectron"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
          }
        }

        if(trig_pass){
	  FillCLHist(sighist_ee, "DiElectron"+tag, eventbase->GetEvent(), muons,electrons,jets, w);
	  FillCLHist(sighist_ee, "DiElectron_Trigger"+tag, eventbase->GetEvent(), muons,electrons,jets, w*trigger_sf*trigger_ps);

	  FillCLHist(sighist_ee, "DiElectron_puW"+tag, eventbase->GetEvent(), muons,electrons,jets, w*pu_reweight*trigger_sf*trigger_ps);
	  FillCLHist(sighist_ee, "DiElectron_IDW"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);


          //if(EtaRegion("BB",electrons))  FillCLHist(sighist_ee, "DiElectron_BB"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
          //if(EtaRegion("EB",electrons))  FillCLHist(sighist_ee, "DiElectron_EB"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
          //if(EtaRegion("EE",electrons))  FillCLHist(sighist_ee, "DiElectron_EE"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
	  
          if(jets.size() >= 2)          FillCLHist(sighist_ee, "DiElectron_dijet"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
	  
	  FillHistPerLumi("DiElectron_llmass"+tag,GetDiLepMass(electrons), ev_weight, 0., 200.,20, 10);
          FillHistPerLumi("DiElectron_njet"+tag,jets.size(), ev_weight, 0., 10.,10, 10);
          FillHistPerLumi("DiElectron_met"+tag,eventbase->GetEvent().MET(snu::KEvent::pfmet), ev_weight, 0., 200.,20, 10);
	}
      }
    }
    else{
      counter("SS",w);
      if(electrons.at(0).Pt() > 20. && electrons.at(1).Pt() > 20.){
	if(GetDiLepMass(electrons)  > 15. ){
	  if(trig_pass){
	    counter("SS_trigger",w);
	    if(jets.size() >= 2)   FillCLHist(sighist_ee, "DiElectron_SSPreselection"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
	    if(jets.size() == 1) FillCLHist(sighist_ee, "DiElectron_SS1Jet"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);  
	    if(eventbase->GetEvent().MET(snu::KEvent::pfmet) > 50.) FillCLHist(sighist_ee, "DiElectron_SSMet"+tag, eventbase->GetEvent(), muons,electrons,jets, ev_weight);
	  }
	}
      }
    }
  }
}
  
void SKTreeValidation::BeginEvent( )throw( LQError ){
  
  Message("In BeginEvent() " , DEBUG);

  return;
}



void SKTreeValidation::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this SKTreeValidationCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void SKTreeValidation::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



