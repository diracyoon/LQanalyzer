// $Id: ExampleAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNDiElectron Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *d
 ***************************************************************************/

/// Local includes
#include "HNDiElectron.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to samevallow inheritance for use in LQCore/core classes
ClassImp (HNDiElectron);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
HNDiElectron::HNDiElectron() :  AnalyzerCore(),  out_electrons(0) {

  mapcounter.clear();

  functionality = HNDiElectron::VALIDATION;


  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNDiElectron");

  Message("In HNDiElectron constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void HNDiElectron::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //

   Message("Making clever hists for Z ->ll test code", INFO);

   //// Initialise Plotting class functions
   /// MakeCleverHistograms ( type, "label")  type can be muhist/elhist/jethist/sighist

   
   MakeCleverHistograms(sighist_ee, "LowMassEE_HNV1");
   MakeCleverHistograms(sighist_ee, "LowMassEE_HN");
   MakeCleverHistograms(sighist_ee, "LowMassEE_HN2");
   MakeCleverHistograms(sighist_ee, "LowMassEE_GENT");
   MakeCleverHistograms(sighist_ee, "LowMassEE_HN_LT");
   MakeCleverHistograms(sighist_ee, "LowMassEE_MVA");
   MakeCleverHistograms(sighist_ee, "HighMassEE_HNV1");
   MakeCleverHistograms(sighist_ee, "HighMassEE_HN");
   MakeCleverHistograms(sighist_ee, "HighMassEE_HN2");

   MakeCleverHistograms(sighist_ee, "HighMassEE_GENT");
   MakeCleverHistograms(sighist_ee, "HighMassEE_MVA");

   return;

   if(functionality == HNDiElectron::ANALYSIS){
     vector<TString> labels;
     labels.push_back("hn"); 
     //labels.push_back("medium");   // pog numbers
     labels.push_back("tight");    // pog numbers 
     
     for(unsigned int il = 0; il < labels.size(); il++){
       TString label = labels.at(il);
       continue;
       MakeCleverHistograms(sighist_ee,label + "_SSee_1jet");
       MakeCleverHistograms(sighist_ee,label + "_SSee_DiJet");
       MakeCleverHistograms(sighist_ee,label + "_SSPreselection");
       MakeCleverHistograms(sighist_ee,label + "_OSee_1jet");
       MakeCleverHistograms(sighist_ee,label + "_OSee_DiJet");
       MakeCleverHistograms(sighist_ee,label + "_OSPreselection");
       MakeCleverHistograms(sighist_ee,label + "_LowMass");
       MakeCleverHistograms(sighist_ee,label + "_MediumMass");
       MakeCleverHistograms(sighist_ee,label + "_HighMass");
       MakeCleverHistograms(sighist_ee,label + "_LowMassCR");
       MakeCleverHistograms(sighist_ee,label + "_MediumMassCR");
       MakeCleverHistograms(sighist_ee,label + "_HighMassCR");
     }
   }
   else if(functionality == HNDiElectron::VALIDATION){

     vector<TString> labels;
     labels.push_back("");
     labels.push_back("truthmatched_");

     vector<TString> ids;
     ids.push_back("POGTIGHT");
     ids.push_back("MVATIGHT");
     ids.push_back("HNTIGHT");
     
     for(unsigned int ii = 0; ii <ids.size(); ii++){
       for(unsigned int il = 0; il < labels.size(); il++){
	 TString label = ids[ii]+ labels.at(il);
	 
	 MakeCleverHistograms(sighist_ee,label+"OSNoZDiElectron");
	 MakeCleverHistograms(sighist_ee,label+"OSZDiElectron");
	 MakeCleverHistograms(sighist_ee,label+"OSDiElectronw1");
	 MakeCleverHistograms(sighist_ee,label+"OSDiElectronw2");
	 MakeCleverHistograms(sighist_ee,label+"OSDiElectronw3");
	 MakeCleverHistograms(sighist_ee,label+"OSDiElectronw4");
	 MakeCleverHistograms(sighist_ee,label+"OSDiElectronw5");
	 
	 MakeCleverHistograms(sighist_ee,label+"OSDiElectronDiJet");
	 MakeCleverHistograms(sighist_ee,label+"OSDiElectronDiJet1bw1");
	 MakeCleverHistograms(sighist_ee,label+"OSDiElectronDiJet1bw2");
	 MakeCleverHistograms(sighist_ee,label+"OSDiElectronDiJet1bw4");
	 MakeCleverHistograms(sighist_ee,label+"OSDiElectronDiJet1bw3");
	 MakeCleverHistograms(sighist_ee,label+"OSDiElectronDiJet2b");
	 MakeCleverHistograms(sighist_ee,label+"SSDiElectron");
	 MakeCleverHistograms(sighist_ee,label+"SSZDiElectron");
	 MakeCleverHistograms(sighist_ee,label+"SSnoZDiElectron");
	 MakeCleverHistograms(sighist_ee,label+"SSDiElJet");
	 MakeCleverHistograms(sighist_ee,label+"SS1Jet");
	 MakeCleverHistograms(sighist_e,label+"SingleElectron_Wregion");
	 MakeCleverHistograms(sighist_eeee,label+"ZZElectron1");
	 MakeCleverHistograms(sighist_eeee,label+"ZZElectron2");
	 
	 MakeCleverHistograms(sighist_eee,label+"ZgElectron");
	 MakeCleverHistograms(sighist_eee,label+"ZgElectron2");
	 MakeCleverHistograms(sighist_eee,label+"WZElectron");
	 MakeCleverHistograms(sighist_eee,label+"WZElectron_pu");
	 
	 MakeCleverHistograms(sighist_ee,label+"_LT");
	 
       }
     }

     vector<TString> mmlabels;
     mmlabels.push_back("");
     mmlabels.push_back("truthmatched_");                                                                                                                        
     vector<TString> mmids;
     mmids.push_back("GENT");
     mmids.push_back("MMHNTIGHT");
     mmids.push_back("MMPOGTIGHT");
     
     for(unsigned int ii = 0; ii < mmids.size(); ii++){
       for(unsigned int il = 0; il < mmlabels.size(); il++){
	 TString label = mmids[ii] + mmlabels.at(il);
	 MakeCleverHistograms(sighist_mm,label+"_LT");
       }
     }
   }
   

   /// Validation signalplots
   /*
   MakeCleverHistograms(sighist_ee, "SIGNALVALIDATION_EE" );
   MakeCleverHistograms(sighist_mm, "SIGNALVALIDATION_MM" );
   */

   return;
}


void HNDiElectron::ExecuteEvents()throw( LQError ){
  
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;

  if(!isData)weight*= MCweight;
  
  /// Interested to see effect of MET filters on tail of met
  if(isData) FillHist("MET_PFMet_uncleaned" , eventbase->GetEvent().PFMET(), weight,  0. , 10000., 1000);
  counter("NoCut",weight);
  
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return;
  counter("Vertex",weight);
  weight*=WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi);
  if(!PassMETFilter()) return;     /// Initial event cuts :                                 
  counter("METFiltert",weight);

  if(isData) FillHist("MET_PFMet_cleaned" , eventbase->GetEvent().PFMET(), weight,  0. , 10000., 1000);                                                                                                                                                                      
  std::vector<snu::KMuon> muons_veto = GetMuons("MUON_HN_VETO",false);
  std::vector<snu::KElectron> electrons_veto = GetElectrons(true,true,"ELECTRON_HN_VETO");

  std::vector<snu::KElectron> electrons_hn1=GetElectrons(true,true,"ELECTRON_HN_TIGHT");
  std::vector<snu::KElectron> electrons_hn2=GetElectrons(true,true,"ELECTRON_HN_TIGHTv4");

  if(electrons_hn2.size() ==2){
   
    std::vector<snu::KTruth> truthColl= eventbase->GetTruth();
    int LepType=GetLeptonType(electrons_hn2[0],truthColl);
    cout << electrons_hn2[0].GetType() << " LepType="  << LepType <<endl;
  }
  std::vector<snu::KElectron> electrons_gent=GetElectrons(true,true,"ELECTRON_GENT_TIGHT");
  std::vector<snu::KElectron> electrons_mva=GetElectrons(true,true,"ELECTRON_MVA_TIGHT");
  std::vector<snu::KJet> jets =  GetJets("JET_HN");
  std::vector<snu::KJet> alljets =  GetJets("JET_NOLEPTONVETO");
  std::vector<snu::KFatJet> fatjets = GetFatJets("FATJET_HN");
  std::vector<snu::KJet> hnjets = GetJetsWFT("JET_HN","FATJET_HN");             

 

  std::vector<snu::KElectron> electrons_loose=GetElectrons(true,true , "ELECTRON_HN_FAKELOOSE");

  vector<TString> IDs;
  IDs.push_back("ELECTRON_HN_EFF_PT");
  IDs.push_back("ELECTRON_HN_EFF_DXY");
  IDs.push_back("ELECTRON_HN_EFF_DXY2");
  IDs.push_back("ELECTRON_HN_EFF_DZ");
  IDs.push_back("ELECTRON_HN_EFF_DZ2");
  IDs.push_back("ELECTRON_HN_EFF_charge");
  IDs.push_back("ELECTRON_HN_EFF_conv");
  IDs.push_back("ELECTRON_HN_EFF_dxysig");
  IDs.push_back("ELECTRON_HN_EFF_iso");
  IDs.push_back("ELECTRON_HN_EFF_gentmva");
  IDs.push_back("ELECTRON_HN_TIGHT");
  IDs.push_back("ELECTRON_HN_TIGHTv4");

  bool _diel =   isData ?  (k_channel.Contains("DoubleEG")) : true ;
  bool _singleel =   isData ?  (k_channel.Contains("SingleElectron")) : true ;

  
  std::vector<snu::KElectron> electrons_nocut=GetElectrons(true,true,"ELECTRON_NOCUT");

  FillCLHist(sighist_ee,"NoCut" , eventbase->GetEvent(), muons_veto, electrons_nocut,jets, alljets,weight);

  FillEventCutFlow(0, "NoCut", weight);
  FillEventCutFlow(1, "NoCut", weight);
  for(unsigned int itrig = 0 ; itrig < 2; itrig++){
    TString dataset="";

    if(itrig==0 && !_diel) continue;
    if(itrig==1 && !_singleel) continue;
    if(itrig==0){
      dataset="DoubleEG";
      if(!PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) continue;
    }
    if(itrig==1){
      dataset="SingleElectron";
      if(PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) continue;
      if(!PassTrigger("HLT_Ele27_WPTight_Gsf_v"))  continue;
    }
    FillEventCutFlow(itrig, "Trigger", weight);

    FillCLHist(sighist_ee,dataset+"Trigger" , eventbase->GetEvent(), muons_veto, electrons_nocut,jets, alljets,weight);
  }
  if(muons_veto.size() == 0 && electrons_veto.size()==2){
    for(unsigned int iid=0; iid < IDs.size(); iid++){
      
      TString elid = IDs[iid];
      std::vector<snu::KElectron> electrons_eff=GetElectrons(true,true,elid);
      
      for(unsigned int itrig = 0 ; itrig < 2; itrig++){
	
	if(SameCharge(electrons_eff)){
	  
	  TString dataset="";
	  
	  if(itrig==0 && !_diel) continue;
	  if(itrig==1 && !_singleel) continue;
	  if(itrig==0){
	    dataset="DoubleEG";
	    if(!PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) continue;
	    if(electrons_eff[0].Pt() < 25 || electrons_eff[1].Pt() < 15) continue;
	  }
	  if(itrig==1){
	    dataset="SingleElectron";
	    if(PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) continue;
	    if(!PassTrigger("HLT_Ele27_WPTight_Gsf_v"))  continue;
	    if(electrons_eff[0].Pt() < 30 || electrons_eff[1].Pt() < 10) continue;
	  }
	  
	  FillEventCutFlow(itrig, elid, weight);
	  FillCLHist(sighist_ee, dataset+elid, eventbase->GetEvent(), muons_veto, electrons_eff,jets, alljets,weight);
	  if(IDs[iid] == "ELECTRON_HN_TIGHT"){
	    
	    if(jets.size() > 1) FillEventCutFlow(itrig, "DiJet", weight);
	    if((jets.size() + fatjets.size()) > 1 )  FillEventCutFlow(itrig, "NewDiJet", weight);
	  }
	}
	else if(electrons_eff.size()==2){
	  TString dataset="";

          if(itrig==0 && !_diel) continue;
          if(itrig==1 && !_singleel) continue;
          if(itrig==0){
            dataset="DoubleEG";
            if(!PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) continue;
            if(electrons_eff[0].Pt() < 25 || electrons_eff[1].Pt() < 15) continue;
          }
          if(itrig==1){
            dataset="SingleElectron";
            if(PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) continue;
            if(!PassTrigger("HLT_Ele27_WPTight_Gsf_v"))  continue;
            if(electrons_eff[0].Pt() < 30 || electrons_eff[1].Pt() < 10) continue;
          }

          FillCLHist(sighist_ee, dataset+"OS"+elid, eventbase->GetEvent(), muons_veto, electrons_eff,jets, alljets,weight);
	  
	}
	
	if(IDs[iid] == "ELECTRON_HN_TIGHT"){
	  
	  
	  if(CheckSignalRegion(true,electrons_eff, jets, alljets,"Low", weight)) FillEventCutFlow(itrig, "LowMass", weight);
	  if(CheckSignalRegion(true,electrons_eff, jets, alljets,"", weight)) FillEventCutFlow(itrig, "MidMass", weight);
	  if(CheckSignalRegion(true,electrons_eff, jets, alljets,"High", weight)) FillEventCutFlow(itrig, "HighMass", weight);
	}
	
      }
    }
    
  }
  
  FillHist("effNoCut", 1 ,weight, 0., 5., 5.);


  if(muons_veto.size() == 0 && electrons_veto.size()==2){

     for(unsigned int itrig = 0 ; itrig < 2; itrig++){
      TString dataset="";
      
      if(itrig==0 && !_diel) continue;
      if(itrig==1 && !_singleel) continue;
      if(itrig==0){
	dataset="DoubleEG";
	if(!PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) continue;
	if(electrons_veto[0].Pt() < 25 || electrons_veto[1].Pt() < 15) continue;
      }
      if(itrig==1){
	dataset="SingleElectron";
	if(PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) continue;
	if(!PassTrigger("HLT_Ele27_WPTight_Gsf_v"))  continue;
	if(electrons_veto[0].Pt() < 30 || electrons_veto[1].Pt() < 10) continue;
      }
      
      std::vector<snu::KElectron> hntight_fake  = GetElectrons(true, true, "ELECTRON_HN_FAKELOOSEST");


      if(isData){
	if(CheckSignalRegion(true,hntight_fake , jets, alljets,"Low", weight)){
	  if(fabs(hntight_fake[1].dxySig()) > 4){
	    if(PassID(hntight_fake[0],"ELECTRON_HN_TIGHT")){
	      vector<snu::KElectron> fakeel;
	      fakeel.push_back(hntight_fake[1]);
	      FillCLHist(sighist_ee, dataset+"LowMassEE_HN_LT", eventbase->GetEvent(), muons_veto, fakeel,jets, alljets,weight);
	    }
	  }
	}
	if(CheckSignalRegion(true,hntight_fake , jets, alljets,"High", weight)){
	  if(fabs(hntight_fake[1].dxySig()) > 4){
            if(PassID(hntight_fake[0],"ELECTRON_HN_TIGHT")){
	      
	      vector<snu::KElectron> fakeel;
	      fakeel.push_back(hntight_fake[1]);
	      FillCLHist(sighist_ee, dataset+"HighMassEE_HN_LT", eventbase->GetEvent(), muons_veto,fakeel,jets, alljets,weight);
	    }
	  }
	}
      }
      else{
	if(CheckSignalRegion(true,hntight_fake , jets, alljets,"Low", weight)){
	  if(PassID(hntight_fake[0],"ELECTRON_HN_TIGHT")){
	    if(NonPrompt(hntight_fake[1])){
		vector<snu::KElectron> fakeel;
		fakeel.push_back(hntight_fake[1]);
		FillCLHist(sighist_ee, dataset+"LowMassEE_HN_LT", eventbase->GetEvent(), muons_veto, fakeel,jets, alljets,weight);
	    }
	  }
        }
        if(CheckSignalRegion(true,hntight_fake , jets, alljets,"High", weight)){
	  if(PassID(hntight_fake[0],"ELECTRON_HN_TIGHT")){
            if(NonPrompt(hntight_fake[1])){
	      vector<snu::KElectron> fakeel;
	      fakeel.push_back(hntight_fake[1]);
	      FillCLHist(sighist_ee, dataset+"HighMassEE_HN_LT", eventbase->GetEvent(), muons_veto,fakeel,jets, alljets,weight);
	    }
	  }
	}
      }      
      
      if(CheckSignalRegion(true,electrons_hn1 , jets, alljets,"Low", weight))        FillEventCutFlow(2+itrig, "LowMassEE_HN", weight);
      if(CheckSignalRegion(true,electrons_hn1 , jets, alljets,"High", weight))        FillEventCutFlow(4+itrig, "HighMassEE_HN", weight);
      if(CheckSignalRegion(true,electrons_hn2 , jets, alljets,"Low", weight))	FillEventCutFlow(2+itrig, "LowMassEE_HN2", weight);
      if(CheckSignalRegion(true,electrons_hn2 , jets, alljets,"High", weight))	FillEventCutFlow(4+itrig, "HighMassEE_HN2", weight);


      
      
      if(CheckSignalRegion(true,electrons_gent , jets, alljets,"Low", weight)){
	FillCLHist(sighist_ee, dataset+"LowMassEE_GENT", eventbase->GetEvent(), muons_veto, electrons_gent,jets, alljets,weight);
	FillEventCutFlow(2+itrig, "LowMassEE_GENT", weight);

      }
      if(CheckSignalRegion(true,electrons_gent , jets, alljets,"High", weight)){
	FillCLHist(sighist_ee, dataset+"HighMassEE_GENT", eventbase->GetEvent(), muons_veto, electrons_gent,jets, alljets,weight);
	FillEventCutFlow(4+itrig, "HighMassEE_GENT", weight);

      }
      if(CheckSignalRegion(true,electrons_mva , jets, alljets,"Low", weight)){
	FillCLHist(sighist_ee, dataset+"LowMassEE_MVA", eventbase->GetEvent(), muons_veto, electrons_mva,jets, alljets,weight);
	FillEventCutFlow(2+itrig, "LowMassEE_MVA", weight);

      }
      if(CheckSignalRegion(true,electrons_mva , jets, alljets,"High", weight)){
	FillCLHist(sighist_ee, dataset+"HighMassEE_MVA", eventbase->GetEvent(), muons_veto, electrons_mva,jets, alljets,weight);
	FillEventCutFlow(4+itrig, "HighMassEE_MVA", weight);

      }
    }
  }
  return;
  
  //MakeValidationPlots(weight);
  
  if(functionality == HNDiElectron::VALIDATION) {
    MakeControlPlots(1,"",weight); /// Uses all MC events (no truth matching)
    MakeControlPlots(2,"truthmatched_",weight); //// removes fake leptons but keeps
    MakeMMControlPlots(1,"",weight); 
    MakeMMControlPlots(2,"truthmatched_",weight); 
  }
}
  

void HNDiElectron::FillEventCutFlow(int cf,TString cut,  float weight){

  if(cf==0){
    if(GetHist( "mm_eventcutflow_DoubleEG")) {
      GetHist( "mm_eventcutflow_DoubleEG")->Fill(cut,weight);
      
    }
    else{
      
      vector<TString> IDs;
      IDs.push_back("NoCut");
      IDs.push_back("Trigger");
      IDs.push_back("ELECTRON_HN_EFF_PT");
      IDs.push_back("ELECTRON_HN_EFF_DXY");
      IDs.push_back("ELECTRON_HN_EFF_DXY2");
      IDs.push_back("ELECTRON_HN_EFF_DZ");
      IDs.push_back("ELECTRON_HN_EFF_DZ2");
      IDs.push_back("ELECTRON_HN_EFF_charge");
      IDs.push_back("ELECTRON_HN_EFF_conv");
      IDs.push_back("ELECTRON_HN_EFF_dxysig");
      IDs.push_back("ELECTRON_HN_EFF_iso");
      IDs.push_back("ELECTRON_HN_EFF_gentmva");
      IDs.push_back("ELECTRON_HN_TIGHT");
      IDs.push_back("DiJet");
      IDs.push_back("NewDiJet");
      IDs.push_back("LowMass");
      IDs.push_back("MidMass");
      IDs.push_back("HighMass");
      AnalyzerCore::MakeHistograms( "mm_eventcutflow_DoubleEG",IDs.size(),0.,float(IDs.size()));    
  
      for(unsigned int iid=0; iid < IDs.size(); iid++){

	TString elid = IDs[iid];
	GetHist( "mm_eventcutflow_DoubleEG")->GetXaxis()->SetBinLabel(1+iid,elid);
      }
      GetHist( "mm_eventcutflow_DoubleEG")->Fill(cut,weight);

    }
  }
  if(cf==1){
    if(GetHist( "mm_eventcutflow_SingleElectron")) {
      GetHist( "mm_eventcutflow_SingleElectron")->Fill(cut,weight);

    }
    else{
      
      vector<TString> IDs;
      IDs.push_back("NoCut");
      IDs.push_back("Trigger");
      IDs.push_back("ELECTRON_HN_EFF_PT");
      IDs.push_back("ELECTRON_HN_EFF_DXY");
      IDs.push_back("ELECTRON_HN_EFF_DXY2");
      IDs.push_back("ELECTRON_HN_EFF_DZ");
      IDs.push_back("ELECTRON_HN_EFF_DZ2");
      IDs.push_back("ELECTRON_HN_EFF_charge");
      IDs.push_back("ELECTRON_HN_EFF_conv");
      IDs.push_back("ELECTRON_HN_EFF_dxysig");
      IDs.push_back("ELECTRON_HN_EFF_iso");
      IDs.push_back("ELECTRON_HN_EFF_gentmva");
      IDs.push_back("ELECTRON_HN_TIGHT");
      IDs.push_back("DiJet");
      IDs.push_back("NewDiJet");
      IDs.push_back("LowMass");
      IDs.push_back("MidMass");
      IDs.push_back("HighMass");
      AnalyzerCore::MakeHistograms( "mm_eventcutflow_SingleElectron",IDs.size(),0.,float(IDs.size()));
      
      for(unsigned int iid=0; iid < IDs.size(); iid++){
	
        TString elid = IDs[iid];
	GetHist( "mm_eventcutflow_SingleElectron")->GetXaxis()->SetBinLabel(1+iid,elid);
      }
      GetHist( "mm_eventcutflow_SingleElectron")->Fill(cut,weight);
      
    }
  }
  if(cf==2){
    if(GetHist( "mm_eventcutflow_LowMass_DoubleEG")) {
      GetHist( "mm_eventcutflow_LowMass_DoubleEG")->Fill(cut,weight);
    }
    else{

      vector<TString> IDs;
      IDs.push_back("LowMassEE_HN");
      IDs.push_back("LowMassEE_HN2");
      IDs.push_back("LowMassEE_MVA");
      IDs.push_back("LowMassEE_GENT");
      AnalyzerCore::MakeHistograms( "mm_eventcutflow_LowMass_DoubleEG",IDs.size(),0.,float(IDs.size()));

      for(unsigned int iid=0; iid < IDs.size(); iid++){

        TString elid = IDs[iid];
        GetHist( "mm_eventcutflow_LowMass_DoubleEG")->GetXaxis()->SetBinLabel(1+iid,elid);
      }
      GetHist( "mm_eventcutflow_LowMass_DoubleEG")->Fill(cut,weight);

    }
  }
  
  if(cf==3){
    if(GetHist( "mm_eventcutflow_LowMass_SingleElectron")) {
      GetHist( "mm_eventcutflow_LowMass_SingleElectron")->Fill(cut,weight);
      
    }
    else{

      vector<TString> IDs;
      IDs.push_back("LowMassEE_HN");
      IDs.push_back("LowMassEE_HN2");
      IDs.push_back("LowMassEE_MVA");
      IDs.push_back("LowMassEE_GENT");
      AnalyzerCore::MakeHistograms( "mm_eventcutflow_LowMass_SingleElectron",IDs.size(),0.,float(IDs.size()));

      for(unsigned int iid=0; iid < IDs.size(); iid++){

        TString elid = IDs[iid];
        GetHist( "mm_eventcutflow_LowMass_SingleElectron")->GetXaxis()->SetBinLabel(1+iid,elid);
      }
      GetHist( "mm_eventcutflow_LowMass_SingleElectron")->Fill(cut,weight);
      
    }
  }

  if(cf==4){
    if(GetHist( "mm_eventcutflow_HighMass_DoubleEG")) {
      GetHist( "mm_eventcutflow_HighMass_DoubleEG")->Fill(cut,weight);

    }
    else{

      vector<TString> IDs;
      IDs.push_back("HighMassEE_HN");
      IDs.push_back("HighMassEE_HN2");
      IDs.push_back("HighMassEE_MVA");
      IDs.push_back("HighMassEE_GENT");
      AnalyzerCore::MakeHistograms( "mm_eventcutflow_HighMass_DoubleEG",IDs.size(),0.,float(IDs.size()));

      for(unsigned int iid=0; iid < IDs.size(); iid++){

        TString elid = IDs[iid];
        GetHist( "mm_eventcutflow_HighMass_DoubleEG")->GetXaxis()->SetBinLabel(1+iid,elid);
      }
      GetHist( "mm_eventcutflow_HighMass_DoubleEG")->Fill(cut,weight);

    }
  }
  
  
  
}


void HNDiElectron::MakeValidationPlots(float w){
  
  std::vector<snu::KElectron> electrons=GetElectrons(true,true,"ELECTRON_NOCUT");
  std::vector<snu::KElectron> electrons_nofake=GetElectrons(true,false,"ELECTRON_NOCUT");
  std::vector<snu::KMuon> muons=GetMuons("MUON_NOCUT");
  std::vector<snu::KJet> jets =  GetJets("JET_HN");
  std::vector<snu::KJet> alljets =  GetJets("JET_NOLEPTONVETO");
  
  //if(muons.size()>  0)cout << muons[0].RelIso04() << endl;
  FillCLHist(sighist_eeee,"All", eventbase->GetEvent(),muons,electrons,  jets, alljets, w);
  FillCLHist(sighist_mmmm,"All_MM", eventbase->GetEvent(),muons,electrons,  jets,alljets,  w);

  std::vector<snu::KMuon> muons_veto = GetMuons("MUON_HN_VETO",false);
  if(muons_veto.size() > 2) return;
  std::vector<snu::KMuon> gent=  GetMuons("MUON_HNGENT_TIGHT",false);
  std::vector<snu::KMuon> gent2=  GetMuons("MUON_HNGENT_TIGHT2",false);
  std::vector<snu::KMuon> gent3=  GetMuons("MUON_HNGENT_TIGHT3",false);
  std::vector<snu::KMuon> gent4=  GetMuons("MUON_HNGENT_TIGHT4",false);
  std::vector<snu::KMuon> pogtight_tight  = GetMuons("MUON_POG_TIGHT",false);
  std::vector<snu::KMuon> hntight  = GetMuons("MUON_HN_TIGHT",false);
  std::vector<snu::KMuon> pogmedium_tight  = GetMuons("MUON_POG_MEDIUM",false);
  
  std::vector<snu::KMuon> hntight_fake  = GetMuons("MUON_HN_Tight_FAKELOOSEST",false);
  std::vector<snu::KMuon> hnmedium_fake  = GetMuons("MUON_HN_Medium_FAKELOOSEST",false);

  std::vector<snu::KMuon> pogtight_fake  = GetMuons("MUON_POG_FAKETIGHT",false);
  std::vector<snu::KMuon> pogmedium_fake  = GetMuons("MUON_POG_FAKEMEDIUM",false);

  std::vector<snu::KMuon> pogtight  = GetMuons("MUON_POG_TIGHT",false);
  std::vector<snu::KMuon> pogmedium  = GetMuons("MUON_POG_MEDIUM",false);


  FillHist("effnmuons_fake_tight", hntight_fake.size(),1., 0., 5., 5.);
  FillHist("effnmuons_fake_medium",  hnmedium_fake.size(), 1., 0., 5., 5.);
  
  FillHist("effnmuons_pfake_tight", pogtight_fake.size(), 1., 0., 5., 5.);
  FillHist("effnmuons_pfake_medium", pogmedium_fake.size(), 1., 0., 5., 5.);
  
  FillHist("effnmuons_p_tight", pogtight.size(),1.,  0., 5., 5.);
  FillHist("effnmuons_p_medium", pogmedium.size(), 1., 0., 5., 5.);
  

  
  std::vector<snu::KMuon> pogmedium_tight_hn;
  if(NBJet(alljets, snu::KJet::CSVv2, snu::KJet::Medium) ==0){
    if(electrons.size() ==0)
    if(jets.size() > 1){
      for(unsigned int im=0; im < pogmedium_tight.size(); im++){
	if(fabs(pogmedium_tight[im].dXY()) > 0.005) continue;
	if(fabs(pogmedium_tight[im].dZ()) > 0.04) continue;
	if(fabs(pogmedium_tight[im].dXYSig()) > 3.) continue;
	pogmedium_tight_hn.push_back(pogmedium_tight[im]);
      }
      
      if(SameCharge(pogmedium_tight_hn)) {
	if(GetDiLepMass(pogmedium_tight_hn) > 10){
	  if(GetDiLepMass(pogmedium_tight_hn) < 70)FillCLHist(sighist_mm,"POG_MM_MEDHN_LM", eventbase->GetEvent(),pogmedium_tight_hn,electrons,  jets, alljets, weight);
	  FillCLHist(sighist_mm,"POG_MM_MEDHN_HM", eventbase->GetEvent(),pogmedium_tight_hn,electrons,  jets, alljets, weight);

	}
      }
    }
  }

  

  float wg=weight;
  float wt=weight;
  float wm=weight;

  if(k_running_nonprompt){
    if(jets.size() > 1){
      gent=  GetMuons("MUON_HNGENT_LOOSE",false);
      pogtight_tight = GetMuons("MUON_HN_LOOSE",false);
      pogmedium_tight = GetMuons("MUON_HN_LOOSE",false);
      if(SameCharge(gent))wg = 1.;//m_datadriven_bkg->Get_DataDrivenWeight_MM(false, gent, PassID(gent[0],"MUON_HNGENT_TIGHT"),  PassID(gent[1],"MUON_HNGENT_TIGHT"), "gent","gent",false, false, "ptcorr_eta",0.1,0.1,false, false);
      else wg=0.;
      cout << "wg  = " << wg << endl;
      if(SameCharge(pogmedium_tight))wm = 1.;//m_datadriven_bkg->Get_DataDrivenWeight_MM(false, pogmedium_tight, PassID(pogmedium_tight[0],"MUON_POG_MEDIUM"),  PassID(pogmedium_tight[1],"MUON_POG_MEDIUM"),"pogmedium",  "pogmedium", false, false ,"ptcorr_eta",0.25,0.25,false,false);
      else wm = 0.;
      if(SameCharge(pogtight_tight)) wt =1.;// m_datadriven_bkg->Get_DataDrivenWeight_MM(false, pogtight_tight, PassID(pogtight_tight[0],"MUON_POG_TIGHT"),  PassID(pogtight_tight[1],"MUON_POG_TIGHT"), "pogtight","pogtight", false, false, "ptcorr_eta",0.15,0.15,false,false);
      else wt=0.;
      
      
      if(SameCharge(pogmedium_tight))cout << "gent = " << gent.size() << "  pogmedium_tight=" <<pogmedium_tight.size() << " " << wm << " " << wg << endl;
      if(SameCharge(gent)){
	if(gent[0].Pt() > 25)FillCLHist(sighist_mm,"GENT_MM", eventbase->GetEvent(),gent,electrons,  jets, alljets, wg);
	FillHist("Weights_gent", wg, 1., -5., 5., 100.);
      }
      if(SameCharge(pogmedium_tight)){
	if(pogmedium_tight[0].Pt() > 25.) FillCLHist(sighist_mm,"POG_MED_MM", eventbase->GetEvent(),pogmedium_tight,electrons,  jets,alljets, wm);
	FillHist("Weights_pogm", wm, 1., -5., 5., 100.);
      }
      if(SameCharge(pogtight_tight)){
	if(pogtight_tight[0].Pt() > 25)FillCLHist(sighist_mm,"POG_TIGHT_MM", eventbase->GetEvent(),pogtight_tight,electrons,  jets, alljets, wt);
	FillHist("Weights_pogt", wt, 1., -5., 5., 100.);
      }
    }
  }
  if(SameCharge(gent) && jets.size() > 1)  {
    if(gent[0].Pt() > 25){
      FillCLHist(sighist_mm,"GENT_MM", eventbase->GetEvent(), gent,electrons,  jets, alljets,wg);
      counter("GENT",weight);

    }
  }
  if(SameCharge(gent2) && jets.size() > 1)  {
    if(gent2[0].Pt() > 25){
      counter("GENT2",weight);
      FillCLHist(sighist_mm,"GENT2_MM", eventbase->GetEvent(), gent2,electrons,  jets, alljets,wg);

    }
    cout << PassID(gent2[0],"MUON_HNGENT_TIGHT2") << " " << PassID(gent2[1],"MUON_HNGENT_TIGHT2") << endl;
    cout << PassID(gent2[0],"MUON_HNGENT_TIGHT3") << " " << PassID(gent2[1],"MUON_HNGENT_TIGHT3") << endl;
    cout << gent3.size() << endl;
    cout << gent2[0].dZ() << " " << gent2[1].dZ() << endl;
  }
  
  if(SameCharge(gent3) && jets.size() > 1)  {
    if(gent3[0].Pt() > 25){
      counter("GENT3",weight);
      FillCLHist(sighist_mm,"GENT3_MM", eventbase->GetEvent(), gent3,electrons,  jets, alljets,wg);

    }
  }
  
  if(SameCharge(gent4) && jets.size() > 1)  {
    if(gent4[0].Pt() > 25){
      counter("GENT4",weight);
      FillCLHist(sighist_mm,"GENT4_MM", eventbase->GetEvent(), gent3,electrons,  jets, alljets,wg);

    }
  }

  if(SameCharge(hntight) && jets.size() >1)  {
    if(hntight[0].Pt() > 25){
      counter("NNTIGHT",weight);

    }
  }

  if(SameCharge(pogmedium_tight) && jets.size() > 1) {
    if(pogmedium_tight[0].Pt() > 25.){
      FillCLHist(sighist_mm,"POG_MED_MM", eventbase->GetEvent(),pogmedium_tight,electrons,  jets, alljets,wm);
      if(GetDiLepMass(pogmedium_tight) > 10. && GetDiLepMass(pogmedium_tight) < 70)      FillCLHist(sighist_mm,"POG_MED_MM_LM", eventbase->GetEvent(),pogmedium_tight,electrons,  jets, alljets,wm);
      counter("POGM",weight);
    }

  }
  if(SameCharge(pogtight_tight) && jets.size() > 1)  {
    if(pogtight_tight[0].Pt() > 25){
      FillCLHist(sighist_mm,"POG_TIGHT_MM", eventbase->GetEvent(),pogtight_tight,electrons,  jets, alljets,wt);
      if(GetDiLepMass(pogtight_tight) > 10. &&GetDiLepMass(pogtight_tight) <70)   FillCLHist(sighist_mm,"POG_TIGHT_MM_LM", eventbase->GetEvent(),pogtight_tight,electrons,  jets, alljets,wt);

      counter("POGT",weight);
    }
  }
  
  
  return;

  /// Drop the pt of the lepton to 7 GeV so that ZZ CR has more stats                                                                                                           
  bool keepcf=true;
  bool keepnp=true;
  
  TString pog_veto_elid = "ELECTRON_PTETA"; // POG ID + loose dxy/dz cuts  (no cc or dxysig)          
  TString pog_tight_elid = "ELECTRON_POG_TIGHT"; // POG ID + loose dxy/dz cuts  (no cc or dxysig)          
  std::vector<snu::KElectron> electronVetoColl=GetElectrons(keepcf, keepnp,  pog_veto_elid,15., 2.5);  /// IF k_running_nonprompt loose id                       
  std::vector<snu::KElectron> electronTightColl=GetElectrons(keepcf, keepnp, pog_tight_elid,15., 2.5);  /// IF k_running_nonprompt loose id                       
  std::vector<snu::KElectron> electronTightColl_all=GetElectrons(true, true, pog_tight_elid,15., 2.5);
  std::vector<snu::KMuon> muonVetoColl=GetMuons("MUON_HN_VETO");
  std::vector<snu::KJet> jets_pu =  GetJets("JET_HN_PU");

  std::vector<snu::KElectron> electronColl=GetElectrons(keepcf, keepnp, "ELECTRON_PTETA", 5., 2.5);  /// IF k_running_nonprompt loose id                                                                                  

  
  float tempweight =  w*WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi);
  
  for(unsigned int i=0; i < electronVetoColl.size(); i++){
    

    if(!electronVetoColl[i].PassesConvVeto()){
      cout << "electronVetoColl[i] pt = " << electronVetoColl[i].Pt() << " " << electronVetoColl[i].Eta() << " " << electronVetoColl[i].Phi() << endl;                        
      
      cout << "NonMatched : COnversion " << endl;                                                                                                                             
      cout << "PassesConvVeto " << i+1 << "  = " << electronVetoColl[i].PassesConvVeto() << endl;                                                                             
      cout << "Loose el: " << i+1 << " particle type = " << electronVetoColl[i].GetParticleType() << endl;                                                                    
      cout << "Loose el: " << i+1 << " mother type = " << electronVetoColl[i].GetMotherType() << endl;                                                                        
      cout << "Loose el: " << i+1 << " ismatched = " << electronVetoColl[i].MCMatched() << endl;                                                                              
      cout << "Loose el: " << i+1 << " iscf = " << electronVetoColl[i].MCIsCF()  << endl;                                                                                     
      cout << "Loose el: " << i+1 << " isconv = " << electronVetoColl[i].MCIsFromConversion()  << endl;                                                                       
      cout << "Loose el: " << i+1 << " paretn=tau = " << electronVetoColl[i].MCFromTau()  << endl;                                                                            
      cout << "Loose el: " << i+1 << " matched index = " << electronVetoColl[i].MCTruthIndex()  << endl;                                                                      
      cout << "Loose el: " << i+1 << " pdgid matched = " << electronVetoColl[i].MCMatchedPdgId()  << endl;                                                                    
      cout << "Loose el: " << i+1 << " mother pdgid = " << electronVetoColl[i].MotherPdgId() << " [ " << electronVetoColl[i].MotherTruthIndex()  << " ] " <<  endl ;          
      TruthPrintOut();                                                                                                                                                        
    }

    float ftype=1;
    if(electronVetoColl[i].MCMatched() ) {
      
      if(electronVetoColl[i].GetMotherType() == snu::KElectron::Z || electronVetoColl[i].GetMotherType() == snu::KElectron::W){
        FillHist("Matched_W_Z_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
        FillHist("Matched_W_Z_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
      }

      if(electronVetoColl[i].MCIsCF()){
	FillHist("Matched_CF_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
	FillHist("Matched_CF_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
	if(electronVetoColl[i].MCIsFromConversion()){
	  FillHist("Matched_CF_conv_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
	  FillHist("Matched_CF_conv_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
	  
	}
	else{
	  FillHist("Matched_CF_notconv_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
	  FillHist("Matched_CF_notconv_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
	}
      }
      else   {
        FillHist("Matched_noCF_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
	FillHist("Matched_noCF_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
      }

      if(electronVetoColl[i].MCIsFromConversion()){
	FillHist("Matched_CONV_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
	FillHist("Matched_CONV_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
	if(!electronVetoColl[i].MCIsCF()){
	  FillHist("Matched_CONV_notCF_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
	  FillHist("Matched_CONV_notCFiso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
	}
      }
      else  {
	FillHist("Matched_noCONV_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
	FillHist("Matched_noCONV_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
	
      }
      
      if(electronVetoColl[i].MCFromTau()){
	FillHist("Matched_TAUD_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
	FillHist("Matched_TAUD_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);

	}
      else           {
	FillHist("Matched_noTAUD_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
	FillHist("Matched_noTAUD_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
	
      }
      
      if(electronVetoColl[i].MCIsCF()) ftype+= 10;
      if(electronVetoColl[i].MCIsFromConversion()) ftype+= 100;
      if(electronVetoColl[i].MCFromTau())  ftype+= 1000;

      FillHist("Matched_particleType", electronVetoColl[i].GetParticleType(),tempweight, 0., 6., 6);
      FillHist("Matched_particleType_motherType", electronVetoColl[i].GetParticleType(), electronVetoColl[i].GetMotherType(), tempweight, 0., 6., 6, 0., 5., 5);

      if(electronVetoColl[i].GetParticleType() == snu::KElectron::CONV_CF)         FillHist("Matched_CONV_CF_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
      if(electronVetoColl[i].GetParticleType() == snu::KElectron::CONV_NONECF)         FillHist("Matched_CONV_NONECF_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
      if(!electronVetoColl[i].MCIsFromConversion() && electronVetoColl[i].MCIsCF()) FillHist("Matched_NONECONV_CF_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);

      if(electronVetoColl[i].MCMatchedPdgId() ==22){
	FillHist("Matched_PH_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
        FillHist("Matched_PH_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
      } 
    }
    else{
      ftype*= -1.;

      if(electronVetoColl[i].GetMotherType() == snu::KElectron::Z || electronVetoColl[i].GetMotherType() == snu::KElectron::W){
        FillHist("NonMatched_W_Z_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
        FillHist("NonMatched_W_Z_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
      }
      

      if(electronVetoColl[i].MCIsCF()){
        FillHist("NonMatched_CF_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
        FillHist("NonMatched_CF_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
        if(electronVetoColl[i].MCIsFromConversion()){
          FillHist("NonMatched_CF_conv_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
          FillHist("NonMatched_CF_conv_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);

        }
        else{
          FillHist("NonMatched_CF_notconv_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
          FillHist("NonMatched_CF_notconv_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
        }
      }
      else   {
        FillHist("NonMatched_noCF_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
        FillHist("NonMatched_noCF_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
      }

      if(electronVetoColl[i].MCIsFromConversion()){
        FillHist("NonMatched_CONV_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
        FillHist("NonMatched_CONV_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
        if(!electronVetoColl[i].MCIsCF()){
          FillHist("NonMatched_CONV_notCF_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
          FillHist("NonMatched_CONV_notCFiso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
	}
      }
      else  {
        FillHist("NonMatched_noCONV_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
        FillHist("NonMatched_noCONV_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);

      }
      if(electronVetoColl[i].MCFromTau()){
        FillHist("NonMatched_TAUD_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
	FillHist("NonMatched_TAUD_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
      }
     
      if(electronVetoColl[i].MCMatchedPdgId() ==22){
	FillHist("NonMatched_PH_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
        FillHist("NonMatched_PH_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
      }
      if(electronVetoColl[i].MCIsCF()) ftype-= 10;
      if(electronVetoColl[i].MCIsFromConversion()) ftype-= 100;
      if(electronVetoColl[i].MCFromTau())  ftype-= 1000;
      FillHist("NonMatched_particleType", electronVetoColl[i].GetParticleType(), tempweight,0., 6., 6);
      FillHist("NonMatched_particleType_motherType", electronVetoColl[i].GetParticleType(), electronVetoColl[i].GetMotherType(), tempweight,0., 6., 6, 0., 5., 5);
      if(electronVetoColl[i].GetParticleType() == snu::KElectron::CONV_CF)         FillHist("NonMatched_CONV_CF_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
      if(electronVetoColl[i].GetParticleType() == snu::KElectron::CONV_NONECF)         FillHist("NonMatched_CONV_NONECF_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
      if(!electronVetoColl[i].MCIsFromConversion() && electronVetoColl[i].MCIsCF()) FillHist("NonMatched_NONECONV_CF_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
      if(electronVetoColl[i].GetParticleType() == snu::KElectron::PHOTONFAKE) FillHist("NonMatched_PHFake_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
      if(electronVetoColl[i].GetParticleType() == snu::KElectron::FAKE) FillHist("NonMatched_Fake_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
      
      

    }
    FillHist("EType", ftype, tempweight, -1200, 1200, 480);
  }
  
  vector<bool> vmatched;
  vector<bool> vnonmatched;
  for(unsigned int i=0; i < electronVetoColl.size(); i++){
    if(electronVetoColl[i].MCMatched()){
      vmatched.push_back(true);
      vnonmatched.push_back(false);
    }
    else{
      vnonmatched.push_back(true);
      vmatched.push_back(false);
    }
  }
  
  for(unsigned int i=0; i < electronVetoColl.size(); i++){
    if(electronVetoColl[i].MCMatched()){
      if(electronVetoColl[i].GetParticleType() ==4){
	if(electronVetoColl[i].GetMotherType()){
	  if(!electronVetoColl[i].MCIsCF()){
	    if(electronVetoColl[i].MCIsFromConversion()){
	      if(!electronVetoColl[i].MCFromTau()){
		FillHist("Matched_t1_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
		FillHist("Matched_t1_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
		vmatched[i]=false;
	      }
	    }
	  }
	}
      }
    }
  }  
  
  for(unsigned int i=0; i < electronVetoColl.size(); i++){
    if(electronVetoColl[i].MCMatched()){
      if(electronVetoColl[i].GetParticleType() ==0){
        if(electronVetoColl[i].GetMotherType()==3){
          if(!electronVetoColl[i].MCIsCF()){
            if(!electronVetoColl[i].MCIsFromConversion() ){
              if(!electronVetoColl[i].MCFromTau()){
                FillHist("Matched_t2_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
		FillHist("Matched_t2_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
                vmatched[i]=false;
	      } 
            }
          }
	}
      } 
    }
  }

  for(unsigned int i=0; i < electronVetoColl.size(); i++){
    if(electronVetoColl[i].MCMatched()){
      if(electronVetoColl[i].GetParticleType() ==5){
        if(electronVetoColl[i].GetMotherType()==1){
          if(electronVetoColl[i].MCIsCF()){
            if(!electronVetoColl[i].MCIsFromConversion() ){
              if(!electronVetoColl[i].MCFromTau() ){
                FillHist("Matched_t3_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
		FillHist("Matched_t3_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
		vmatched[i]=false;
	      } 
            }
          }
	}
      } 
    }
  }
  

  bool cout_nonmatched(true);
  for(unsigned int i=0; i < electronVetoColl.size(); i++){
    if(!electronVetoColl[i].MCMatched()){
      if(electronVetoColl[i].GetParticleType() ==3){
        if(electronVetoColl[i].GetMotherType()==1 || electronVetoColl[i].GetMotherType()==4 ){
          if(electronVetoColl[i].MCIsCF()){
            if(electronVetoColl[i].MCIsFromConversion()){
              if(!electronVetoColl[i].MCFromTau()){
                FillHist("NonMatched_t1_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
		FillHist("NonMatched_t1_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
		vnonmatched[i]=false;
	      } 
            }
          }
	}
      } 
    }
  }

  for(unsigned int i=0; i < electronVetoColl.size(); i++){
    if(!electronVetoColl[i].MCMatched()){
      if(electronVetoColl[i].GetParticleType() ==4){
	if(!electronVetoColl[i].MCIsCF()){
	  if(electronVetoColl[i].MCIsFromConversion()){
	    if(!electronVetoColl[i].MCFromTau()){
	      vnonmatched[i]=false;
	      FillHist("NonMatched_t2_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
	      FillHist("NonMatched_t2_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
	    }
	  }
        }
      } 
    }
  }


  if(electronVetoColl.size() == 1){
    if(electronVetoColl[0].MCMatched()){
      if(!electronVetoColl[0].MCIsFromConversion()){
	FillHist("SingleLepMatched_dxy", electronVetoColl[0].dxy(), tempweight, -1., 1., 200);
	FillHist("SingleLepMatched_iso", electronVetoColl[0].PFRelMiniIso(), tempweight, 0., 1., 100);
      }
    }
  }

  if(electronVetoColl.size() ==2){
    
    if(electronVetoColl[0].MCMatched() && electronVetoColl[1].MCMatched()){
      if(!electronVetoColl[1].MCIsFromConversion()){
	if(!electronVetoColl[1].MCIsFromConversion()){
	  if(electronVetoColl[0].GetMotherType() == 1){
	    if(electronVetoColl[1].GetMotherType() == 1){
	      if(electronVetoColl[0].GetParticleType() == snu::KElectron::PROMPT){
		if(electronVetoColl[1].GetParticleType() == snu::KElectron::PROMPT){
		  FillHist("ZMatched_dxy", electronVetoColl[0].dxy(), tempweight, -1., 1., 200);
		  FillHist("ZMatched_iso", electronVetoColl[0].PFRelMiniIso(), tempweight, 0., 1., 100);

		}
	      }
	    }
	  }
	}
      }
    }
    

    for(unsigned int i=0; i < electronVetoColl.size(); i++){
      if(!vnonmatched[i])  continue;
      if(!electronVetoColl[i].MCMatched()){
	/*
	cout << "electronVetoColl[i] pt = " << electronVetoColl[i].Pt() << " " << electronVetoColl[i].Eta() << " " << electronVetoColl[i].Phi() << endl;
	
	cout << "NonMatched : COnversion " << endl;
	cout << "PassesConvVeto " << i+1 << "  = " << electronVetoColl[i].PassesConvVeto() << endl;
	cout << "Loose el: " << i+1 << " particle type = " << electronVetoColl[i].GetParticleType() << endl;
	cout << "Loose el: " << i+1 << " mother type = " << electronVetoColl[i].GetMotherType() << endl;
	cout << "Loose el: " << i+1 << " ismatched = " << electronVetoColl[i].MCMatched() << endl;
	cout << "Loose el: " << i+1 << " iscf = " << electronVetoColl[i].MCIsCF()  << endl;
	cout << "Loose el: " << i+1 << " isconv = " << electronVetoColl[i].MCIsFromConversion()  << endl;
	cout << "Loose el: " << i+1 << " paretn=tau = " << electronVetoColl[i].MCFromTau()  << endl;
	cout << "Loose el: " << i+1 << " matched index = " << electronVetoColl[i].MCTruthIndex()  << endl;
	cout << "Loose el: " << i+1 << " pdgid matched = " << electronVetoColl[i].MCMatchedPdgId()  << endl;
	cout << "Loose el: " << i+1 << " mother pdgid = " << electronVetoColl[i].MotherPdgId() << " [ " << electronVetoColl[i].MotherTruthIndex()  << " ] " <<  endl ;
	TruthPrintOut();
	*/
      }
    }

    
    
    if(electronVetoColl[0].MCMatched() && electronVetoColl[1].MCMatched()){
      for(unsigned int i=0; i < electronVetoColl.size(); i++){
	
	float ftype=1;
	if(electronVetoColl[i].MCMatched() ) {
	  if(electronVetoColl[i].MCIsCF()) ftype+= 10;
	  if(electronVetoColl[i].MCIsFromConversion()) ftype+= 100;
	  if(electronVetoColl[i].MCFromTau())  ftype+= 1000;
	}
	else{
	  ftype*= -1.;
	  if(electronVetoColl[i].MCIsCF()) ftype-= 10;
	  if(electronVetoColl[i].MCIsFromConversion()) ftype-= 100;
	  if(electronVetoColl[i].MCFromTau())  ftype-= 1000;
	  
	}
	FillHist("EType_matched_dilep", ftype, tempweight, -1200, 1200, 480);
	if(SameCharge(electronVetoColl)){
	  FillHist("EType_ss_matched_dilep", ftype, tempweight, -1200, 1200, 480);
	  
	  FillHist("NonMatched_particleType_motherType", electronVetoColl[0].GetParticleType(), electronVetoColl[0].GetMotherType(), tempweight,0., 6., 6, 0., 5., 5);
	  FillHist("NonMatched_particleType_motherType", electronVetoColl[1].GetParticleType(), electronVetoColl[1].GetMotherType(), tempweight,0., 6., 6, 0., 5., 5);
	  
	  FillHist("EType_matched_dilep_mothers", electronVetoColl[0].GetMotherType(),electronVetoColl[1].GetMotherType(), tempweight, 0., 5., 5, 0., 5., 5);
	  FillHist("EType_matched_dilep_partcileType", electronVetoColl[0].GetParticleType(), electronVetoColl[1].GetParticleType(),tempweight, 0, 6., 6 , 0, 6., 6); 
	  if(vmatched[0] && vmatched[1] ){
	    cout << "electronVetoColl[0] pt = " << electronVetoColl[0].Pt() << " " << electronVetoColl[0].Eta() << " " << electronVetoColl[0].Phi() << endl;
	    cout << "electronVetoColl[1] pt = " << electronVetoColl[1].Pt() << " " << electronVetoColl[1].Eta() << " " << electronVetoColl[1].Phi() << endl;

	    for(unsigned int iell = 0 ; iell < electronVetoColl.size(); iell++){
	      cout << "electronVetoColl["<< iell+1<<"] pt = " << electronVetoColl[iell].Pt() << " " << electronVetoColl[iell].Eta() << " " << electronVetoColl[iell].Phi() << endl;
	      cout << "PassesConvVeto " << i+1 << "  = " << electronVetoColl[i].PassesConvVeto() << endl;
	      
	      cout << "Loose el: " << iell+1 << " particle type = " << electronVetoColl[iell].GetParticleType() << endl;
	      cout << "Loose el: " << iell+1 << " mother type = " << electronVetoColl[iell].GetMotherType() << endl;
	      cout << "Loose el: " << iell+1 << " ismatched = " << electronVetoColl[iell].MCMatched() << endl;
	      cout << "Loose el: " << iell+1 << " iscf = " << electronVetoColl[iell].MCIsCF()  << endl;
	      cout << "Loose el: " << iell+1 << " isconv = " << electronVetoColl[iell].MCIsFromConversion()  << endl;
	      cout << "Loose el: " << iell+1 << " paretn=tau = " << electronVetoColl[iell].MCFromTau()  << endl;
	      cout << "Loose el: " << iell+1 << " matched index = " << electronVetoColl[iell].MCTruthIndex()  << endl;
	      cout << "Loose el: " << iell+1 << " pdgid matched = " << electronVetoColl[iell].MCMatchedPdgId()  << endl;
	      cout << "Loose el: " << iell+1 << " mother pdgid = " << electronVetoColl[iell].MotherPdgId() << " [ " << electronVetoColl[iell].MotherTruthIndex()  << " ] " <<  endl;
	      TruthPrintOut();

	      if(electronVetoColl[iell].MCIsCF()  ){
                FillHist("SSMatched_cf_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
                FillHist("SSMatched_cf_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
		if(!electronVetoColl[iell].MCIsFromConversion()){
		  FillHist("SSMatched_cf_noconv_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
		  FillHist("SSMatched_cf_noconv_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);

		}
              }
	      
	      if(electronVetoColl[iell].MCIsFromConversion()  ){
		FillHist("SSMatched_conv_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
		FillHist("SSMatched_conv_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
	      }
	      else if(fabs(electronVetoColl[iell].MotherPdgId()) == 15){
		if(!electronVetoColl[iell].MCIsCF()){
		  FillHist("SSMatched_tau_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
		  FillHist("SSMatched_tau_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
		}
		else{
		  FillHist("SSMatched_cftau_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
                  FillHist("SSMatched_cftau_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);

		}
	      }
	      else {
		if(! ((fabs(electronVetoColl[iell].MotherPdgId())==23) || (fabs(electronVetoColl[iell].MotherPdgId())==24))){
		  FillHist("SSMatched_x_dxy", electronVetoColl[i].dxy(), tempweight, -1., 1., 200);
		  FillHist("SSMatched_x_iso", electronVetoColl[i].PFRelMiniIso(), tempweight, 0., 1., 100);
		}
	      }
	    }
	  }
	}

	
      }
    }
    
  }/// 2 lep
  
  
  if(electronTightColl.size() == 2){
    snu::KParticle Ztt = electronTightColl[0] + electronTightColl[1];
    float tempweight =  weight*WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi);
    if(electronTightColl.at(0).MCFromTau() &&electronTightColl.at(1).MCFromTau()) {
      FillHist("NTauTau_matched", 1, tempweight, 0., 3., 3);
      FillHist("Z_tautau_m_matched", Ztt.M(), tempweight, 0., 200., 40);
    }
    else FillHist("NTauTau_matched", 0, tempweight, 0., 3., 3);
  }


  if(electronTightColl_all.size() == 2){
    snu::KParticle Ztt = electronTightColl_all[0] + electronTightColl_all[1];
    float tempweight =weight*WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi);
    if(electronTightColl_all.at(0).MCFromTau() &&electronTightColl_all.at(1).MCFromTau()) {

      FillHist("NTauTau", 1, tempweight, 0., 3., 3);
      FillHist("Z_tautau_m", Ztt.M(), tempweight, 0., 200., 40);
    }
    else FillHist("NTauTau", 0, tempweight, 0., 3., 3);
    
  }
  

  return;
}


void HNDiElectron::MakeMMControlPlots(int method, TString methodtag, float w)throw( LQError ){
  
  FillByMMTriggerTrigger(2, "MMHNTIGHT", method,methodtag,w);
  FillByMMTriggerTrigger(2, "MMGENT", method,methodtag,w);
  
}

void HNDiElectron::FillByMMTriggerTrigger(int iel_trig, TString ID,int method, TString methodtag, double evw){
  
  
  TString muid_veto="MUON_HN_VETO";
  TString muid_tight="MUON_HN_TIGHT";
  if(k_running_nonprompt)muid_tight="MUON_HN_LOOSE";
  TString muid_tight_fixed="MUON_HN_TIGHT";
  TString muid_loose="MUON_HN_LOOSE";
  methodtag+=ID;
  TString elid_loose="ELECTRON_POG_FAKELOOSE";
  
  if(ID == "MMGENT"){
    muid_veto="MUON_HNGENT_LOOSE";
    muid_tight="MUON_HNGENT_TIGHT";
    if(k_running_nonprompt)muid_tight="MUON_HNGENT_LOOSE";
    muid_loose="MUON_HNGENT_LOOSE";
    
  }
  
  std::vector<snu::KJet> jets =  GetJets("JET_HN");
  int njets = jets.size();
  std::vector<snu::KElectron> electronLoose=GetElectrons(true, true,   elid_loose, 15., 2.5);
  std::vector<snu::KMuon> muonLoose=GetMuons(muid_loose,false);
  
  if(!methodtag.Contains("truth")){
    muonLoose=GetMuons(muid_loose,true);
  }


  Float_t motherbins[18] = { 0., 1., 7.,8., 12., 13., 16., 17., 22., 23., 24., 25., 100., 200., 300., 400., 500., 1000};
  Float_t jetbins[8] = { 0., 1.,2.,3.,4.,5.,6., 7.};

  if(IsFakeEvent(muonLoose) <= 0){
    for(unsigned int imu=0; imu < muonLoose.size(); imu++){
      FillHist(methodtag+"NotFakeJustType_",muonLoose[imu].GetType(),1., -45, 45., 90);
      if(PassID(muonLoose[imu],muid_tight))        FillHist(methodtag+"NotFakeJustType_tight",muonLoose[imu].GetType(),1., -45, 45., 90);
    }
  }
  if(IsFakeEvent(muonLoose) > 0){
    bool LT=false;
    int iFail=-1;
    evw*=WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi);
    

    for(unsigned int imu=0; imu < IsFakeEvent(muonLoose); imu++){
      for(int i=-45; i < 45; i++){
	
	if(muonLoose[imu].Pt() < 10) continue;
	
	if(muonLoose[imu].GetType() == i){
	  stringstream ss;
	  ss << i;
	  string str = ss.str();
	  FillHist(methodtag+"ByType_mother_Full_"+str, fabs(muonLoose[imu].MotherPdgId()),evw, 0, 5000., 5000);
	  if(PassID(muonLoose[imu],muid_tight))             FillHist(methodtag+"ByType_mother_Full_"+str+"_tight", fabs(muonLoose[imu].MotherPdgId()), evw, 0, 5000., 5000);
	  if(muonLoose[imu].MotherTruthIndex() > 0 && muonLoose[imu].MotherTruthIndex() < eventbase->GetTruth().size()){
	    FillHist(methodtag+"ByType_grandmother_Full_"+str, fabs(eventbase->GetTruth().at(muonLoose[imu].MotherTruthIndex()).PdgId()), evw,0, 5000., 5000);
	    if(PassID(muonLoose[imu],muid_tight))            FillHist(methodtag+"ByType_grandmother_Full_"+str+"_tight", fabs(eventbase->GetTruth().at(muonLoose[imu].MotherTruthIndex()).PdgId()), evw,0, 5000., 5000);
	  }
	}
      }
      
      if(muonLoose[imu].GetType() != 1){
	FillHist(methodtag+"JustType_",muonLoose[imu].GetType(),1., -45, 45., 90);
	if(PassID(muonLoose[imu],muid_tight))     FillHist(methodtag+"JustType_tight",muonLoose[imu].GetType(),1., -45, 45., 90);
      }
      if(muonLoose[imu].GetType() >= 8 && muonLoose[imu].GetType() <=10){
	if(muonLoose[imu].GetType() == 9){
	  FillHist(methodtag+"_vmassnu", GetVirtualMass(13,true), 1., 0., 200., 200);
	  if(PassID(muonLoose[imu],muid_tight))    FillHist(methodtag+"_vmassnu_tight", GetVirtualMass(13,true), 1., 0., 200., 200);
	}
	else{
	  FillHist(methodtag+"_vmass", GetVirtualMass(13,false), 1., 0., 200., 200);
          if(PassID(muonLoose[imu],muid_tight))    FillHist(methodtag+"_vmass_tight", GetVirtualMass(13,false), 1., 0., 200., 200);

	}
	if(muonLoose[imu].GetType()  == 2 && fabs(muonLoose[imu].MotherPdgId()) > 5000.) cout << "Mother = " << fabs(muonLoose[imu].MotherPdgId()) << endl;

      }
      
    }
    
    


    if(!TruthMatched(muonLoose[0])) {
      FillHist(methodtag+"TypeByNJet_tm", muonLoose[0].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);
      if(PassID(muonLoose[0],muid_tight))       FillHist(methodtag+"TypeByNJet_tm_tight", muonLoose[0].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);
      FillHist(methodtag+"MotherByNJet_tm", fabs(muonLoose[0].MotherPdgId()), njets, evw, motherbins,17, jetbins,7);
      if(PassID(muonLoose[0],muid_tight))         FillHist(methodtag+"MotherByNJet_tm_tight", fabs(muonLoose[0].MotherPdgId()), njets, evw, motherbins,17, jetbins,7);
    }
    if(!TruthMatched(muonLoose[1])) {
      FillHist(methodtag+"TypeByNJet_tm", muonLoose[1].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);
      if(PassID(muonLoose[1],muid_tight))       FillHist(methodtag+"TypeByNJet_tm_tight", muonLoose[1].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);

      FillHist(methodtag+"MotherByNJet_tm", fabs(muonLoose[1].MotherPdgId()), njets, evw, motherbins,17, jetbins,7);
      if(PassID(muonLoose[1],muid_tight))         FillHist(methodtag+"MotherByNJet_tm_tight", fabs(muonLoose[1].MotherPdgId()), njets, evw, motherbins,17, jetbins,7);

    }
    if(!muonLoose[0].MCMatched()) {
      FillHist(methodtag+"TypeByNJet_mcmatched", muonLoose[0].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);
      if(PassID(muonLoose[0],muid_tight))       FillHist(methodtag+"TypeByNJet_mcmatched_tight", muonLoose[0].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);
      FillHist(methodtag+"MotherByNJet_mcmatched", fabs(muonLoose[0].MotherPdgId()), njets, evw, motherbins,17, jetbins,7);
      if(PassID(muonLoose[0],muid_tight))         FillHist(methodtag+"MotherByNJet_mcmatched_tight", fabs(muonLoose[0].MotherPdgId()), njets, evw, motherbins,17, jetbins,7);
    }
    if(!muonLoose[1].MCMatched()) {
      FillHist(methodtag+"TypeByNJet_mcmatched", muonLoose[1].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);
      if(PassID(muonLoose[1],muid_tight))       FillHist(methodtag+"TypeByNJet_mcmatched_tight", muonLoose[1].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);
      
      FillHist(methodtag+"MotherByNJet_mcmatched", fabs(muonLoose[1].MotherPdgId()), njets, evw, motherbins,17, jetbins,7);
      if(PassID(muonLoose[1],muid_tight))         FillHist(methodtag+"MotherByNJet_mcmatched_tight", fabs(muonLoose[1].MotherPdgId()), njets, evw, motherbins,17, jetbins,7);

    }

    
    if(PassID(muonLoose[0],muid_tight) && !PassID(muonLoose[1],muid_tight)) {LT=true; iFail=1;}
    if(!PassID(muonLoose[0],muid_tight) && PassID(muonLoose[1],muid_tight)) {LT=true; iFail=0;}
    if(LT&& iFail>=0){
      if(PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") || PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")){
        if(muonLoose[0].Pt()>20 && muonLoose[1].Pt() > 10){

          FillHist(methodtag+"TypeByNJet", muonLoose[iFail].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);
          if(!muonLoose[iFail].MCMatched())FillHist(methodtag+"NotMatchedTypeByNJet", muonLoose[iFail].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);
	  
	  float lt_weight =evw;
          if(!isData){
            lt_weight*=mcdata_correction->UserPileupWeight(eventbase->GetEvent());
          }
          FillCLHist(sighist_mm, methodtag+"_LT", eventbase->GetEvent(), muonLoose,electronLoose,jets, lt_weight);

        }
      }
    }

  }

  /*
  TString muid_veto="MUON_HN_VETO";
  TString muid_tight="MUON_HN_TIGHT";
  if(k_running_nonprompt)muid_tight="MUON_HN_LOOSE";
  TString muid_tight_fixed="MUON_HN_TIGHT";
  TString muid_loose="MUON_HN_LOOSE";
  methodtag+=ID;
  TString elid_loose="ELECTRON_POG_FAKELOOSE";
  */

  //std::vector<snu::KMuon> muonZZColl=GetMuons(mu_id,keepnp);  /// IF k_running_nonprompt loose id                                                                                                                                             
  
  //FillCLHist(sighist_mm, methodtag+"SSDiMuon", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);




}

void HNDiElectron::MakeControlPlots(int method, TString methodtag, float w)throw( LQError ){

  
  /// e+e-                                                                                                                                                     
  FillByTriggerTrigger(2, "POGTIGHT", method,methodtag,w);
  FillByTriggerTrigger(2, "MVATIGHT", method,methodtag,w);
  FillByTriggerTrigger(2, "HNTIGHT", method,methodtag,w);
  //if(PassTrigger("HLT_Ele27_WPTight_Gsf_v")) FillByTriggerTrigger(1,method,methodtag,w);
  
}

void HNDiElectron::FillByTriggerTrigger(int iel_trig, TString ID,int method, TString methodtag, double evw){

  TString elid_veto="";
  TString elid_tight="";
  TString elid_tight_fixed="";
  TString elid_loose="";
  methodtag+=ID;

  int id=0;
  TString method_fake="";
  if(ID.Contains("POG")){
    elid_veto="ELECTRON_POG_VETO";
    elid_tight="ELECTRON_POG_TIGHT";
    elid_tight_fixed="ELECTRON_POG_TIGHT";
    elid_loose="ELECTRON_POG_FAKELOOSE";

    id=1;
    method_fake="dijet_ajet40";
  }
  if(ID.Contains("MVA")){
    elid_veto="ELECTRON_POG_VETO";
    elid_tight="ELECTRON_MVA_TIGHT";
    elid_tight_fixed="ELECTRON_MVA_TIGHT";
    elid_loose="ELECTRON_MVA_FAKELOOSE";
    method_fake="dijet_ajet40";
    id=2;
  }
  if(ID.Contains("HN")){
    elid_veto="ELECTRON_HN_VETO";
    elid_tight="ELECTRON_HN_TIGHT";
    elid_tight_fixed="ELECTRON_HN_TIGHT";
    elid_loose="ELECTRON_HN_FAKELOOSE";
    method_fake="mva";
    id=3;
  }

  


  std::vector<snu::KElectron> electronVetoColl=GetElectrons(true, false,   elid_veto, 10., 2.5);

  std::vector<snu::KElectron> electronLoose=GetElectrons(true, true,   elid_loose, 10., 2.5);
  
  TString el_elid = elid_tight;

  if(k_running_nonprompt) el_elid =elid_loose;
  bool keepcf=true;
  bool keepnp=false;
  if(method==1){
    keepcf=false;
    keepnp=false;
  }
  else if(method==2){
    // if sample is QCD sample then no events with have prompt electrons so always keep fakes in this case                                                                                                                                                                                                                                                                  
    if(k_sample_name.Contains("QCD")) keepnp=true;
    if(k_sample_name.Contains("qcd")) keepnp=true;
    
  }
  else if(method==3){
    // if sample is QCD sample then no events with have prompt electrons so always keep fakes in this case                                                                                                                                                                                                                                                                  
    if(k_sample_name.Contains("QCD")) keepnp=true;
    if(k_sample_name.Contains("qcd")) keepnp=true;
    /// electronTightColl = POG ID unless k_running_nonprompt is true                                                                                                                                                                                                                                                                                                       
  }
  
  if(!methodtag.Contains("truth")){
    keepcf=true;
    keepnp=true;
  }
  std::vector<snu::KElectron> electronTightColl=GetElectrons(keepcf, keepnp,  el_elid,10., 2.5);  /// IF k_running_nonprompt loose id                                    
      
  std::vector<snu::KJet> jets =  GetJets("JET_HN");
  std::vector<snu::KMuon> muonVetoColl=GetMuons("MUON_HN_VETO");
  int njets = jets.size();
  
  Float_t motherbins[18] = { 0., 1., 7.,8., 12., 13., 16., 17., 22., 23., 24., 25., 100., 200., 300., 400., 500., 1000};
  Float_t jetbins[8] = { 0., 1.,2.,3.,4.,5.,6., 7.};
  
  if(IsFakeEvent(electronLoose) > 0){
    bool LT=false;
    
    evw*= WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi)  ;
    

    for(unsigned int iel=0; iel < IsFakeEvent(electronLoose); iel++){

      for(int i=-45; i < 45; i++){
	
	if(electronLoose[iel].GetType() == i){
          stringstream ss;
          ss << i;
          string str = ss.str();
          FillHist(methodtag+"ByType_mother_Full_"+str, fabs(electronLoose[iel].MotherPdgId()), evw,0, 5000., 5000);
          if(PassID(electronLoose[iel],elid_tight))             FillHist(methodtag+"ByType_mother_Full_"+str+"_tight", fabs(electronLoose[iel].MotherPdgId()), evw,0, 5000., 5000);
          if(electronLoose[iel].MotherTruthIndex() > 0 && electronLoose[iel].MotherTruthIndex() < eventbase->GetTruth().size()){
	    FillHist(methodtag+"ByType_grandmother_Full_"+str, fabs(eventbase->GetTruth().at(electronLoose[iel].MotherTruthIndex()).PdgId()), evw,0, 5000., 5000);
	    if(PassID(electronLoose[iel],elid_tight))            FillHist(methodtag+"ByType_grandmother_Full_"+str+"_tight", fabs(eventbase->GetTruth().at(electronLoose[iel].MotherTruthIndex()).PdgId()),evw, 0, 5000., 5000);
	  }
        }
      }
      
      if(electronLoose[iel].GetType() != 1){
	FillHist(methodtag+"JustType_", electronLoose[iel].GetType(),1., -45, 45., 90);
	if(PassID(electronLoose[iel],elid_tight))     FillHist(methodtag+"JustType_tight", electronLoose[iel].GetType(),1., -45, 45., 90);
      }
      if(electronLoose[iel].GetType() >= 14 && electronLoose[iel].GetType() <=16){
	if(electronLoose[iel].GetType()==15){
	  FillHist(methodtag+"_vmassnu", GetVirtualMass(11,true), 1., 0., 200., 200);
	  if(PassID(electronLoose[iel],elid_tight)) FillHist(methodtag+"_vmassnu_tight", GetVirtualMass(11,true), 1., 0., 200., 200);
	}
	else{
	  FillHist(methodtag+"_vmass", GetVirtualMass(11,false), 1., 0., 200., 200);
          if(PassID(electronLoose[iel],elid_tight)) FillHist(methodtag+"_vmass_tight", GetVirtualMass(11,false), 1., 0., 200., 200);
	}
      }

    }
    
    if(!TruthMatched(electronLoose[0],false)) {
      FillHist(methodtag+"TypeByNJet_tm", electronLoose[0].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);
      if(PassID(electronLoose[0],elid_tight))       FillHist(methodtag+"TypeByNJet_tm_tight", electronLoose[0].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);
      FillHist(methodtag+"MotherByNJet_tm", fabs(electronLoose[0].MotherPdgId()), njets, evw, motherbins,17, jetbins,7);
      if(PassID(electronLoose[0],elid_tight))         FillHist(methodtag+"MotherByNJet_tm_tight", fabs(electronLoose[0].MotherPdgId()), njets, evw, motherbins,17, jetbins,7);
    }
    if(!TruthMatched(electronLoose[1],false)) {
      FillHist(methodtag+"TypeByNJet_tm", electronLoose[1].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);
      if(PassID(electronLoose[1],elid_tight))       FillHist(methodtag+"TypeByNJet_tm_tight", electronLoose[1].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);

      FillHist(methodtag+"MotherByNJet_tm", fabs(electronLoose[1].MotherPdgId()), njets, evw, motherbins,17, jetbins,7);
      if(PassID(electronLoose[1],elid_tight))         FillHist(methodtag+"MotherByNJet_tm_tight", fabs(electronLoose[1].MotherPdgId()), njets, evw, motherbins,17, jetbins,7);

    } 
    
    if(!electronLoose[0].MCMatched()) {
      FillHist(methodtag+"TypeByNJet_mcmatched", electronLoose[0].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);
      if(PassID(electronLoose[0],elid_tight))       FillHist(methodtag+"TypeByNJet_mcmatched_tight", electronLoose[0].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);
      FillHist(methodtag+"MotherByNJet_mcmatched", fabs(electronLoose[0].MotherPdgId()), njets, evw, motherbins,17, jetbins,7);
      if(PassID(electronLoose[0],elid_tight))         FillHist(methodtag+"MotherByNJet_mcmatched_tight", fabs(electronLoose[0].MotherPdgId()), njets, evw, motherbins,17, jetbins,7);
    }
    if(!electronLoose[1].MCMatched()) {
      FillHist(methodtag+"TypeByNJet_mcmatched", electronLoose[1].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);
      if(PassID(electronLoose[1],elid_tight))       FillHist(methodtag+"TypeByNJet_mcmatched_tight", electronLoose[1].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);

      FillHist(methodtag+"MotherByNJet_mcmatched", fabs(electronLoose[1].MotherPdgId()), njets, evw, motherbins,17, jetbins,7);
      if(PassID(electronLoose[1],elid_tight))         FillHist(methodtag+"MotherByNJet_mcmatched_tight", fabs(electronLoose[1].MotherPdgId()), njets, evw, motherbins,17, jetbins,7);

    }

    
    int iFail=-1;
    if(PassID(electronLoose[0],elid_tight) && !PassID(electronLoose[1],elid_tight)) {LT=true; iFail=1;}
    if(!PassID(electronLoose[0],elid_tight) && PassID(electronLoose[1],elid_tight)) {LT=true; iFail=0;}
    if(LT&&iFail>=0){
      if(GetDiLepMass(electronLoose) > 101.|| GetDiLepMass(electronLoose) < 81.){
	if(PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")){

	  if(electronLoose[0].Pt()>25 && electronLoose[1].Pt() > 15){
	    FillHist(methodtag+"TypeByNJet", electronLoose[iFail].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);
	    if(!electronLoose[iFail].MCMatched())  FillHist(methodtag+"NotMatchedTypeByNJet", electronLoose[iFail].GetType(), njets, evw,-45.,45.,90, 0, 7., 7);
	    float lt_weight =evw;
	    if(!isData){
	      if(id==1)lt_weight*=   mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronLoose,0);
	      if(id>1)lt_weight*=   mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electronLoose,0);
	     
	      lt_weight*= mcdata_correction->ElectronRecoScaleFactor(electronLoose);
	      lt_weight*=mcdata_correction->UserPileupWeight(eventbase->GetEvent());
	    }
	    FillCLHist(sighist_ee, methodtag+"_LT", eventbase->GetEvent(), muonVetoColl,electronLoose,jets, lt_weight);
	  }
	}
      }
    }
    
    for(int ix=0; ix < 3; ix++){
      
      float val=0.;
      TString labelval="";
      if(ix==0){
	val=GetDiLepMass(electronLoose);
	labelval="llmass";
      }
      if(ix>1){
	val=electronLoose[1].Pt();
	labelval="l2pt";
      }

      
      int pass1=false;
      int pass2=false;
      if(PassID(electronLoose[0],elid_tight)) pass1=true;
      if(PassID(electronLoose[1],elid_tight)) pass2=true;
      if(ix==2&&!pass1&&!pass2) continue;
      if(electronLoose[0].GetType() >  22  || electronLoose[1].GetType() >  22) FillHist(labelval+"_check1",val,1.,0, 200., 100);
      else{
	
	if(electronLoose[0].GetType()  < 6  && electronLoose[1].GetType()  < 6  ){
	  FillHist(labelval+"_check2",val,1.,0, 200., 100);
	}
	else{
	  
	  FillHist("pt1_isloose",electronLoose[0].GetType(),1., 0, 30., 30.);
	  FillHist("pt2_isloose",electronLoose[1].GetType(),1., 0,30., 30.);
	  if(pass1)FillHist("pt1_isltight",electronLoose[0].GetType(),1., 0,30., 30.);
	  if(pass2)FillHist("pt2_isltight",electronLoose[1].GetType(),1., 0,30., 30.);
	  
	  if(electronLoose[0].GetType() == 7 || electronLoose[1].GetType() == 7 )                FillHist(labelval+"_check3",val,1.,0, 200., 100);
	  if(electronLoose[0].GetType() == 8 || electronLoose[1].GetType() == 8 )                FillHist(labelval+"_check4",val,1.,0, 200., 100);
	  if(electronLoose[0].GetType() == 9 || electronLoose[1].GetType() == 9 )                FillHist(labelval+"_check5",val,1.,0, 200., 100);
	  if(electronLoose[0].GetType() == 10 || electronLoose[1].GetType() == 10 )                FillHist(labelval+"_check6",val,1.,0, 200., 100);
	  if(electronLoose[0].GetType() == 11 || electronLoose[1].GetType() == 11 )                FillHist(labelval+"_check7",val,1.,0, 200., 100);
	  if(electronLoose[0].GetType() == 12 || electronLoose[1].GetType() == 12 )                FillHist(labelval+"_check8",val,1.,0, 200., 100);
	  if(electronLoose[0].GetType() == 13 || electronLoose[1].GetType() == 13 )                FillHist(labelval+"_check9",val,1.,0, 200., 100);
	  if(electronLoose[0].GetType() == 14 || electronLoose[1].GetType() == 14 )                FillHist(labelval+"_check10",val,1.,0, 200., 100);
	  if(electronLoose[0].GetType() == 15 || electronLoose[1].GetType() ==15)                FillHist(labelval+"_check11",val,1.,0, 200., 100);
	  if(electronLoose[0].GetType() == 16 || electronLoose[1].GetType() == 16 )                FillHist(labelval+"_check12",val,1.,0, 200., 100);
	  if(electronLoose[0].GetType() == 17 || electronLoose[1].GetType() == 17 )                FillHist(labelval+"_check13",val,1.,0, 200., 100);
	  if(electronLoose[0].GetType() == 18 || electronLoose[1].GetType() == 18 )                FillHist(labelval+"_check14",val,1.,0, 200., 100);
	  if(electronLoose[0].GetType() == 19 || electronLoose[1].GetType() == 19 )                FillHist(labelval+"_check15",val,1.,0, 200., 100);
	  if(electronLoose[0].GetType() == 20 || electronLoose[1].GetType() == 20 )                FillHist(labelval+"_check16",val,1.,0, 200., 100);
	  if(electronLoose[0].GetType() == 21 || electronLoose[1].GetType() == 21 )                FillHist(labelval+"_check17",val,1.,0, 200., 100);
	}
      }
    }
  }

  std::vector<snu::KElectron> electronTightColl_all=GetElectrons(false, false, "ELECTRON_POG_TIGHT",10., 2.5);
 
  /// Drop the pt of the lepton to 7 GeV so that ZZ CR has more stats                                                                        
  std::vector<snu::KElectron> electronZZColl=GetElectrons(keepcf, keepnp, el_elid,7., 2.5);  /// IF k_running_nonprompt loose id                                
  std::vector<snu::KElectron> electronNCColl=GetElectrons(keepcf, keepnp, "ELECTRON_NOCUT",0., 2.5);                                                            
  /// GetKFactor returns k factor for any sample. If no k-fact is found returns 1. IF data event returns 1.;                                                                                                                                                                                                                                                                
  evw *= GetKFactor();
  
  /// Function NBJet returns #jets tagged with csv medium wp tagger with SF applied per jet                                                                                                                                                                                                                                                                               
  /// make validation plots for dielectron events + WZ/ZZ CRs                                                                                                                                                                                                                                                                                                               
  
  float trigger_sf(1.);
  float id_iso_sf(1.);
  float trigger_ps(1.);
  float reco_weight=1.;
  float puweight(1.);
  float puweight2(1.);
  float puweight3(1.);
  float puweight4(1.);
  float ev_weight(1.);
  float ev_weight2(1.);
  float ev_weight3(1.);
  float ev_weight4(1.);


  std::vector<snu::KJet> jets_pu =  GetJets("JET_HN_PU");

  int nbjet = NBJet(jets);



  counter("Control",weight);

  if(muonVetoColl.size() > 0) return;
  counter("MuonVeto",weight);

  /// Simply fill nymber of leptons in samples                                                                                                                                                                                                                                                                                                                              
  FillHist(methodtag+"_nelectron_v", electronVetoColl.size() , evw*WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi),0., 5.,5);
  FillHist(methodtag+"_nelectron_tm_v", electronTightColl_all.size() , evw*WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi),0., 5.,5);
  FillHist(methodtag+"_nelectron_v", electronTightColl.size() , evw*WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi),0., 5.,5);
  FillHist(methodtag+"_muon_v", muonVetoColl.size() , evw*WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi),0., 5.,5);

  FillHist(methodtag+"_njets_v", jets.size() , evw*WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi),0., 5.,5);
  FillHist(methodtag+"_njets_withpuid_v", jets_pu.size() , evw*WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi),0., 5.,5);

  TString  analysis_trigger = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";

  
  if(electronTightColl.size() < 2 ) return;
  counter("2Tight",weight);

  /// make it safer for trigger turn on
  if(iel_trig==2){
    if(electronTightColl[0].Pt() < 25) return;
    if(electronTightColl[1].Pt() < 12) return;
  }
  else if(iel_trig==1){
    if(electronTightColl[0].Pt() < 30) return;
    if(electronTightColl[1].Pt() < 10) return;
  }
  counter("2TightPt",weight);

  if(!isData ||  (isData&&PassTrigger(analysis_trigger))  ){
  
    /// apply trigger eff to MC instead of trigger bit * SF

    if(!isData) evw*=  mcdata_correction->GetDoubleEGTriggerEff(electronTightColl);


    
    counter("Trigger_diel",evw);
    
    
    
    // 4 PU weights are checked
    /// Method 1 uses nvtx reweightg for ALL MC events 
    puweight= mcdata_correction->UserPileupWeight(eventbase->GetEvent());
    /// Method 2 uses NVtx reweighting to periods G+H and standard PU weights to B-F
    if(GetMCPeriodRandom() < 6) puweight2=eventbase->GetEvent().PileUpWeight(snu::KEvent::down);
    else puweight2= mcdata_correction->UserPileupWeight(eventbase->GetEvent(), jets.size());
    /// Method 3 uses standard PU weights to B-H
    puweight3= eventbase->GetEvent().PileUpWeight();
    /// Method 4 uses PU weightes using lower xsec for minbias events
    puweight4 = eventbase->GetEvent().PileUpWeight(snu::KEvent::down);
    
    
    if(electronZZColl.size() >= 4){
      
      /// 4 lepton events 
      
      snu::KParticle ZZpart = electronZZColl[0] + electronZZColl[1] + electronZZColl[2] + electronZZColl[3];
      
      //// set ev_weight OF running fake bkg prediction
      if(k_running_nonprompt){
	ev_weight     = 1.;//m_datadriven_bkg->Get_DataDrivenWeight(false, muonVetoColl, "MUON_HN_TRI_TIGHT", muonVetoColl.size(), electronZZColl, elid_tight_fixed, electronZZColl.size(), elid_loose, method_fake);
      }
      
      
      FillHist(methodtag+"_zz_mass_cut0", ZZpart.M() , WeightByTrigger(analysis_trigger, TargetLumi)*evw/mcdata_correction->GetDoubleEGTriggerEff(electronZZColl), 0., 800.,32);
      FillHist(methodtag+"_zz_mass_cut1", ZZpart.M() , WeightByTrigger(analysis_trigger, TargetLumi)*evw, 0., 800.,32);
      
      /// file:///Users/John/Downloads/AN2016_359_v11.pdf selection
      
      if(!(muonVetoColl.size() > 0 )){
	
	FillHist(methodtag+"_zz_mass_cut2", ZZpart.M() , WeightByTrigger(analysis_trigger, TargetLumi)*evw, 0., 800.,32);	    
	
	// IF MC apply corrections to event weight 
	if(!isData){
	  if(id==1)id_iso_sf=   mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronZZColl,0);
	  if(id>1)id_iso_sf=   mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electronZZColl,0);
	  trigger_ps= WeightByTrigger(analysis_trigger, TargetLumi)  ;
	  reco_weight = mcdata_correction->ElectronRecoScaleFactor(electronZZColl);
	  FillHist(methodtag+"_zz_mass_cut3", ZZpart.M() , id_iso_sf*evw*trigger_ps, 0., 800.,32);	    
	  FillHist(methodtag+"_zz_mass_cut4", ZZpart.M() , id_iso_sf*evw*trigger_ps*reco_weight, 0., 800.,32);	    
	  FillHist(methodtag+"_zz_mass_cut5", ZZpart.M() , id_iso_sf*evw*trigger_ps*reco_weight*puweight, 0., 800.,32);	    
	  ev_weight = evw * id_iso_sf * reco_weight * puweight*trigger_ps;
	}
	
	/// Check no OS pairs have m(ll) < 4. geV OR dR(l,l) < 0.02
	bool pass_zz_sel1=true;
	for(unsigned int iel=0; iel < electronZZColl.size() -1; iel++){
	  for(unsigned int iel2= iel+1; iel2 < electronZZColl.size(); iel2++){
	    if( electronZZColl[iel].DeltaR(electronZZColl[iel2]) < 0.02) pass_zz_sel1=false;
	    snu::KParticle zee = electronZZColl[iel] + electronZZColl[iel2];
	    if((zee.Charge() == 0) && (zee.M() < 4.)) pass_zz_sel1=false;
	  }
	}
	if(pass_zz_sel1)       FillHist(methodtag+"_zz_mass_cut6", ZZpart.M() , ev_weight, 0., 800.,32);
	
	
	int ich(0);
	/// fill index of vector that has e+ and e- used for finding os pairs
	vector<int> plus_i;
	vector<int> minus_i;
	for(unsigned int iel=0; iel < electronZZColl.size(); iel++){
	  ich += electronZZColl[iel].Charge();
	  if(electronZZColl[iel].Charge() > 0) plus_i.push_back(iel);
	  else minus_i.push_back(iel);
	}
	
	bool pass_zz_sel2(false);
	bool pass_zz_sel3(false); /// jaesungs 2 Z in Z peak (pt Zmu1 > 20. )
	
	/// ich = sum of charge and must be 0
	/// then require 2e+ and 2e-
	if(ich == 0 && plus_i.size() == 2 && minus_i.size()==2){
	  
	  for(int izp = 0 ; izp< 1;izp++){
	    for(int izm = 0 ; izm < 2;izm++){
	      
	      /// loop finds Z e+e- pairs
	      snu::KParticle Z1 = electronZZColl[plus_i[izp]] + electronZZColl[minus_i[izm]];
	      
	      if(fabs(Z1.M() - 91.19) < 10.) {
		/// Z1 is Zcandidate
		/// now get second Z candidate indices
		int izp2 = izp;
		if(izp2 ==1)izp2=0;
		else izp2=1;
		int izm2 = izm;
		if(izm2 ==1)izm2=0;
		else izm2=1;
		// Z2 uses non z1 electrons
		snu::KParticle Z2= electronZZColl[plus_i[izp2]] + electronZZColl[minus_i[izm2]];
		if(fabs(Z2.M() - 91.19) <10.) {
		  // this event has 2 Z candidates
		  pass_zz_sel3=true;
		  if( !( (electronZZColl[plus_i[izp]].Pt() > 20 ) || (electronZZColl[minus_i[izm]].Pt() > 20)) ) pass_zz_sel3=false;
		  if( !( (electronZZColl[plus_i[izp2]].Pt() > 20 ) || (electronZZColl[minus_i[izm2]].Pt() > 20)) ) pass_zz_sel3=false;
		  if(pass_zz_sel3) break;
		}
	      }
	    }
	  }
	  
	  snu::KParticle Z1cand1 = electronZZColl[plus_i[0]] + electronZZColl[minus_i[0]];
	  snu::KParticle Z2cand1 = electronZZColl[plus_i[1]] + electronZZColl[minus_i[1]];
	  
	  bool cand_found=false;
	  bool cand1=true;
	  if(!(Z1cand1.M() > 4 && Z1cand1.M() < 120) ) cand1=false;
	  if(!(Z2cand1.M() > 4 && Z2cand1.M() < 120) ) cand1=false;
	  
	  snu::KParticle Z1cand2 = electronZZColl[plus_i[0]] + electronZZColl[minus_i[1]];
	  snu::KParticle Z2cand2 = electronZZColl[plus_i[1]] + electronZZColl[minus_i[0]];
	  bool cand2=true;
	  if(!(Z1cand2.M() > 4 && Z1cand2.M() < 120) ) cand2=false;
	  if(!(Z2cand2.M() > 4 && Z2cand2.M() < 120) ) cand2=false;
	  
	  if(cand1 || cand1) cand_found=true;
	  
	  snu::KParticle Z1, Z2;
	  float Zdiff_ev(0.);
	  if(cand1){
	    float zcand1diff1 = fabs(Z1cand1.M() - 91.1);
	    float zcand1diff2 = fabs(Z2cand1.M() - 91.1);
	    if(zcand1diff1 < zcand1diff2) {Z1= Z1cand1; Z2=Z2cand1; Zdiff_ev = zcand1diff1;}
	    else  {Z1= Z2cand1; Z2=Z1cand1; Zdiff_ev = zcand1diff2;}
	  }
	  if(cand2){
	    float zcand2diff1 = fabs(Z1cand2.M() - 91.1);
	    float zcand2diff2 = fabs(Z2cand2.M() - 91.1);
	    if(zcand2diff1 < zcand2diff2) {
	      if(zcand2diff1 < Zdiff_ev) {Z1 = Z1cand2; Z2 = Z2cand2;}
	    }
	    else{
	      if(zcand2diff2 < Zdiff_ev) {Z1 = Z2cand2; Z1 = Z1cand2;}
	    }
	  }
	  
	  // Z1 is closest to Z mass
	  // Z2 is other os pair
	  if(cand_found&& Z1.M() > 40.) pass_zz_sel2=true;
	  
	  if(pass_zz_sel1&&pass_zz_sel3){
	    FillHist(methodtag+"_zz_mass_cut7", ZZpart.M() , id_iso_sf*evw*trigger_ps*reco_weight*puweight, 0., 800.,32);
	    
	    FillCLHist(sighist_eeee, methodtag+"ZZElectron1", eventbase->GetEvent(), muonVetoColl,electronZZColl,jets, ev_weight);
	    
	  }
	  if(pass_zz_sel1&&pass_zz_sel2){
	    FillHist(methodtag+"_zz_mass_cut8", ZZpart.M() , id_iso_sf*evw*trigger_ps*reco_weight*puweight, 0., 800.,32);
	    
	    FillCLHist(sighist_eeee, methodtag+"ZZElectron2", eventbase->GetEvent(), muonVetoColl,electronZZColl,jets, ev_weight);
	    //if(Z1.M() > 60 && Z2.M() > 60) FillCLHist(sighist_eeee, "ZZElectron3", eventbase->GetEvent(), muonVetoColl,electronZZColl,jets, ev_weight);
	  }
	  
	}
      }
    }// END ZZ CR
  
    /// WZ CR / Zg
    if(electronTightColl.size() == 3 && (nbjet == 0) && (electronTightColl[2].Pt() > 10.)){
      
      
      // pt el 1 and 2 are set earlier to comply with trigger (25,12)
      if(k_running_nonprompt){
	ev_weight     = 1.;//m_datadriven_bkg->Get_DataDrivenWeight(false, muonVetoColl, "MUON_HN_TRI_TIGHT", muonVetoColl.size(), electronTightColl,elid_tight_fixed, electronTightColl.size(),elid_loose, method_fake);
	
      }
      if(!(muonVetoColl.size() > 0 ||  electronVetoColl.size() > 3)){
	if(!isData){
	  if(id==1)id_iso_sf=   mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronTightColl,0);
          if(id>1)id_iso_sf=   mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electronTightColl,0);

	  trigger_ps= WeightByTrigger(analysis_trigger, TargetLumi)  ;
	  reco_weight = mcdata_correction->ElectronRecoScaleFactor(electronTightColl);
	  ev_weight = evw * id_iso_sf * reco_weight * puweight*trigger_ps;
	}
	snu::KParticle Z1 = electronTightColl[0] + electronTightColl[1];
	snu::KParticle Z2 = electronTightColl[0] + electronTightColl[2];
	snu::KParticle Z3 = electronTightColl[1] + electronTightColl[2];
	float diff1 = fabs(Z1.M() - 91.19) ;
	float diff2 = fabs(Z2.M() - 91.19) ;
	float diff3 = fabs(Z3.M() - 91.19) ;
	
	bool m_ll4=true;
	
	float diff = 999999999.;
	float el_zpt = 0.;
	float el_wpt = 0.;
	int iZ1(0);
	int iZ2(0);
	if(Z1.Charge() == 0){
	  diff = diff1;
	  el_zpt = electronTightColl[0].Pt();
	  el_wpt = electronTightColl[2].Pt();
	  iZ1=0;
	  iZ1=1;
	  if(Z1.M() < 4.) m_ll4=false;
	}
	if( (diff2 < diff) && (Z2.Charge() == 0)){
	  diff= diff2;
	  el_zpt= electronTightColl[0].Pt();
	  el_wpt = electronTightColl[1].Pt();
	  iZ1=0;
	  iZ1=2;
	  if(Z2.M() < 4.) m_ll4=false;
	}
	if( (diff3 < diff) && (Z3.Charge() == 0)){
	  diff= diff3;
	  el_zpt= electronTightColl[1].Pt();
	  el_wpt = electronTightColl[0].Pt();
	  iZ1=1;
	  iZ1=2;
	  if(Z3.M() < 4.) m_ll4=false;
	}
	
	if(diff < 10.){
	  if(el_zpt > 20.){
	    if(eventbase->GetEvent().MET(snu::KEvent::pfmet) > 30.){
	      if(el_wpt > 20.){
		snu::KParticle WZ = electronTightColl[0] + electronTightColl[1]+ electronTightColl[2];
		if(WZ.M() > 100){
		  FillCLHist(sighist_eee, methodtag+"WZElectron", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);
		  FillCLHist(sighist_eee, methodtag+"WZElectron_pu", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight*(puweight3/puweight));
		}
	      }
	    }
	  }
	}
	else if(diff > 15.){
	  
	  if(el_zpt > 20.&&m_ll4){
	    if(eventbase->GetEvent().MET(snu::KEvent::pfmet) < 50.){
	      snu::KParticle WZ = electronTightColl[0] + electronTightColl[1]+ electronTightColl[2];
	      snu::KParticle Z = electronTightColl[iZ1] + electronTightColl[iZ2];
	      
	      if(fabs(WZ.M() - 91.1) < 10.){
		FillCLHist(sighist_eee, methodtag+"ZgElectron", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);
		if(jets.size() > 0)    FillCLHist(sighist_eee, methodtag+"ZgElectron2", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);
	      }
	    }
	  }
	}
      }
    }
    
    bool externconv=false;
    for(unsigned int x=0; x <  electronVetoColl.size(); x++){
      if( electronVetoColl[x].MCIsExternalConversion())externconv=true;
    }

    if(electronTightColl.size() == 2){
      
      snu::KParticle Zpart = electronTightColl[0] + electronTightColl[1];
      counter("pog_diel",evw);
      if(Zpart.M() < 10.) goto endDiElLoop;
      counter("pog_diel_mll10",evw);
      
      if(!(muonVetoColl.size() > 0 ||  electronVetoColl.size() > 2)){
	counter("vetolepton_diel",evw);
	
	if(!isData){
	  if(id==1)id_iso_sf=   mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronTightColl,0);
	  if(id>1)id_iso_sf=   mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electronTightColl,0);

	  trigger_ps= WeightByTrigger(analysis_trigger, TargetLumi)  ;
	  reco_weight = mcdata_correction->ElectronRecoScaleFactor(electronTightColl);
	  ev_weight = evw * id_iso_sf * reco_weight * puweight*trigger_ps;
	  ev_weight2 = evw *trigger_ps ;
	  ev_weight3 = evw  * id_iso_sf*trigger_ps;
	  ev_weight4 = evw * id_iso_sf * reco_weight *trigger_ps;
	}
	
	
	if(k_running_nonprompt){
	  //ev_weight =  m_datadriven_bkg->Get_DataDrivenWeight_EE( false,electronTightColl, "ELECTRON16_POG_FAKELOOSE","ELECTRON16_FR_POG_TIGHT","dijet_ajet40");     
	  ev_weight     = 1.;//m_datadriven_bkg->Get_DataDrivenWeight(false, muonVetoColl, "MUON_HN_TRI_TIGHT", muonVetoColl.size(), electronTightColl,elid_tight_fixed, electronTightColl.size(),elid_loose, method_fake);
	
	  ev_weight2 = ev_weight;
	  ev_weight3 = ev_weight;
	  ev_weight4 = ev_weight;
	}
	
	//// Dileton [top + DY CR]
	if(OppositeCharge(electronTightColl)){
	  counter("os_diel",evw);
	  if(electronTightColl.at(0).Pt() > 25. && electronTightColl.at(1).Pt() > 20.){
	    TString tightid = el_elid;
	    if(k_running_nonprompt){
	      if(!PassID(electronTightColl.at(0), tightid) && !PassID(electronTightColl.at(1), tightid) )FillHist(methodtag+"_os_ee_ll_pt2", electronTightColl.at(1).Pt() , ev_weight, 0., 100.,100);
	      if(PassID(electronTightColl.at(0), tightid) && !PassID(electronTightColl.at(1), tightid) )FillHist(methodtag+"_os_ee_tl_pt2", electronTightColl.at(1).Pt() , ev_weight, 0., 100.,100);
	      if(!PassID(electronTightColl.at(0), tightid) && PassID(electronTightColl.at(1), tightid) )FillHist(methodtag+"_os_ee_lt_pt2", electronTightColl.at(1).Pt() , ev_weight, 0., 100.,100);
	      if(PassID(electronTightColl.at(0), tightid) && PassID(electronTightColl.at(1), tightid) )FillHist(methodtag+"_os_ee_tt_pt2", electronTightColl.at(1).Pt() , ev_weight, 0., 100.,100);
	      
	      if(!PassID(electronTightColl.at(0), tightid) && !PassID(electronTightColl.at(1), tightid) )FillHist(methodtag+"_osnw_ee_ll_pt2", electronTightColl.at(1).Pt() ,1., 0., 100.,100);
	      if(PassID(electronTightColl.at(0), tightid) && !PassID(electronTightColl.at(1), tightid) )FillHist(methodtag+"_osnw_ee_tl_pt2", electronTightColl.at(1).Pt() , 1., 0., 100.,100);
	      
	      if(PassID(electronTightColl.at(0), tightid) && PassID(electronTightColl.at(1), tightid) )FillHist(methodtag+"_osnw_ee_lt_pt2", electronTightColl.at(1).Pt() , 1., 0., 100.,100);
	      if(PassID(electronTightColl.at(0), tightid) && PassID(electronTightColl.at(1), tightid) )FillHist(methodtag+"_osnw_ee_tt_pt2", electronTightColl.at(1).Pt() , 1., 0., 100.,100);
	      
	      
	    }
	    
	    if(GetDiLepMass(electronTightColl) < 101. && GetDiLepMass(electronTightColl)  > 81.){
	      FillCLHist(sighist_ee, methodtag+"OSZDiElectron", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);

	      if((nbjet==0)){
		FillHist(methodtag+"_zpeak_ee_w1", GetDiLepMass(electronTightColl), evw*trigger_ps, 0., 200.,400);
		FillHist(methodtag+"_zpeak_ee_w2", GetDiLepMass(electronTightColl), evw*trigger_sf*trigger_ps, 0., 200.,400);
		FillHist(methodtag+"_zpeak_ee_w3", GetDiLepMass(electronTightColl),    evw*id_iso_sf*trigger_sf*trigger_ps, 0., 200.,400);
		FillHist(methodtag+"_zpeak_ee_w4", GetDiLepMass(electronTightColl),    evw*id_iso_sf*trigger_sf*puweight*trigger_ps, 0., 200.,400);
		FillHist(methodtag+"_zpeak_ee_w5a", GetDiLepMass(electronTightColl),    evw*id_iso_sf*trigger_sf*puweight*trigger_ps*reco_weight, 0., 200.,400);
		FillHist(methodtag+"_zpeak_ee_w5b", GetDiLepMass(electronTightColl),    evw*id_iso_sf*trigger_sf*puweight2*trigger_ps*reco_weight, 0., 200.,400);
		FillHist(methodtag+"_zpeak_ee_w5c", GetDiLepMass(electronTightColl),    evw*id_iso_sf*trigger_sf*puweight3*trigger_ps*reco_weight, 0., 200.,400);
		FillHist(methodtag+"_zpeak_ee_w5d", GetDiLepMass(electronTightColl),    evw*id_iso_sf*trigger_sf*puweight4*trigger_ps*reco_weight, 0., 200.,400);
	      }
	    }
	    else   FillCLHist(sighist_ee, methodtag+"OSNoZDiElectron", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);
	    
	    
	    if(electronTightColl[1].Pt() > 30.)FillCLHist(sighist_ee,methodtag+ "OSDiElectronw0", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, 1.);
	    
	    counter("os_dielw1",ev_weight);
	    //m_logger << INFO << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
	    //cout << electronTightColl[0].Pt() << " " << electronTightColl[1].Pt() << endl;
	    //cout << electronTightColl[0].Pt()*electronTightColl[0].SmearFactor() << " " << electronTightColl[1].Pt()*electronTightColl[0].SmearFactor() << endl;

	    FillCLHist(sighist_ee, methodtag+"OSDiElectronw1", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);
	    FillHist("llmass_check",GetDiLepMass(electronTightColl),ev_weight,0, 200., 100);

	    FillCLHist(sighist_ee, methodtag+"OSDiElectronw2", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight2);
	    FillCLHist(sighist_ee, methodtag+"OSDiElectronw3", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight3);
	    FillCLHist(sighist_ee, methodtag+"OSDiElectronw4", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight4);
	    if(electronTightColl[0].GetType() ==1 && electronTightColl[1].GetType()==1){
	      FillCLHist(sighist_ee, methodtag+"OSDiElectronw5", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);
	      if(!externconv)FillCLHist(sighist_ee, methodtag+"OSDiElectronw6", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);
	    }
	    if(jets.size() > 1) {
	      FillCLHist(sighist_ee, methodtag+"OSDiElectronDiJet", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);
	      if(eventbase->GetEvent().MET(snu::KEvent::pfmet) > 70.){
		if(nbjet ==1) FillCLHist(sighist_ee, methodtag+"OSDiElectronDiJet1bw1", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);
		if(nbjet ==1) FillCLHist(sighist_ee, methodtag+"OSDiElectronDiJet1bw2", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight2);
		if(nbjet ==1) FillCLHist(sighist_ee, methodtag+"OSDiElectronDiJet1bw3", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight3);
		if(nbjet ==1) FillCLHist(sighist_ee, methodtag+"OSDiElectronDiJet1bw4", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight4);
		if(nbjet > 1) FillCLHist(sighist_ee, methodtag+"OSDiElectronDiJet2b", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);
	      }
	    }
	  } // pt cuts
	}// e+e-
	else{
	  counter("ss_diel",evw);
	  if(GetDiLepMass(electronTightColl) < 101.&& GetDiLepMass(electronTightColl) > 81.)FillCLHist(sighist_ee, methodtag+"SSZDiElectron", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);
	  else {
	    FillCLHist(sighist_ee, methodtag+"SSnoZDiElectron", eventbase->GetEvent(),muonVetoColl,electronTightColl,jets, ev_weight);
	    
	  }
	  FillCLHist(sighist_ee, methodtag+"SSDiElectron", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);
	  FillCLHist(sighist_ee, methodtag+"SSDiElectron_noW", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, 1.);
	  

	  if(jets.size() > 1)FillCLHist(sighist_ee, methodtag+"SSDiElJet", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);
	  if(jets.size()== 1 && GetDiLepMass(electronTightColl) > 100.)FillCLHist(sighist_ee, methodtag+"SS1Jet", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);
	}
      } // veto leptons 
    }// 2 el
  }


 endDiElLoop:
  return;
  /// e+e-                                                                                                                                                                                                                                                               
}


 

void HNDiElectron::CheckJetsCloseToLeptons(std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets){

  if(electrons.size() == 2) {
    if(SameCharge(electrons)) {
      snu::KParticle looseee = electrons.at(0) + electrons.at(1);
      if(fabs(looseee.M() - 90.) > 20. ){
	for(unsigned int ijet=0; ijet < jets.size(); ijet++){
	  for(unsigned int iel=0; iel < electrons.size(); iel++){
	    if( electrons[iel].DeltaR(jets.at(ijet)) < 0.4){
	      /// Jet usually vetoe
	      FillHist(("ptratio_jetel"), electrons[iel].Pt()/jets[ijet].Pt(), weight, 0., 1., 100);
	    }
	  }
	}
      }
    }
  }
  return;
}


void HNDiElectron::OptimiseID(bool isss){
  
  
}

float  HNDiElectron::GetTightWeight(){
  return 1.;
}

float HNDiElectron::GetMediumWeight(){
  
  return 1.;
}

void HNDiElectron::GetSSSignalEfficiency(float w){
  
  // w = lumi weight * MC gen weight
  
  //// Check efficiency of Selecting two same sign electrons
  if(SameCharge(GetElectrons("ELECTRON_NOCUT")))   FillCutFlow("SS_NoCut",w);
  if(SameCharge(GetElectrons("ELECTRON_PTETA")))   FillCutFlow("SS_PtEta",w);
  // TIGHT POG
  if(SameCharge(GetElectrons("ELECTRON16_POG_TIGHT_SC")))    FillCutFlow("SS_Tight",w*GetTightWeight());
  if(SameCharge(GetElectrons("ELECTRON16_POG_TIGHT_D0"))) FillCutFlow("SS_Tight_d0veto",w*GetTightWeight());
  if(SameCharge(GetElectrons("ELECTRON16_POG_TIGHT_RELISO")))     FillCutFlow("SS_Tight_reliso", w*GetTightWeight());  
  if(SameCharge(GetElectrons("ELECTRON16_POG_TIGHT_CHARGECONS"))) FillCutFlow("SS_Tight_chargeconst", w*GetTightWeight());
  if(SameCharge(GetElectrons("ELECTRON16_POG_TIGHT_DXYSIG")))     FillCutFlow("SS_Tight_dxysig",w*GetTightWeight());

  TString analysis_trigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
  /// Trigger List (unprescaled)                                                                                                                                                                                                                                                
  std::vector<TString> triggerslist;
  triggerslist.push_back(analysis_trigger);

  if(PassTrigger(analysis_trigger) && SameCharge(GetElectrons("ELECTRON16_POG_TIGHT_DXYSIG")))  FillCutFlow("SS_Tight_trigger",w*GetTightWeight());



  // MEDIUM POG
  if(SameCharge(GetElectrons("ELECTRON16_POG_MEDIUM")))    FillCutFlow("SS_Medium",w*GetMediumWeight());
  if(SameCharge(GetElectrons("ELECTRON16_POG_MEDIUM_D0"))) FillCutFlow("SS_Medium_d0veto",w*GetMediumWeight());
  if(SameCharge(GetElectrons("ELECTRON16_POG_MEDIUM_RELISO")))     FillCutFlow("SS_Medium_reliso", w*GetMediumWeight());
  if(SameCharge(GetElectrons("ELECTRON16_POG_MEDIUM_CHARGECONS"))) FillCutFlow("SS_Medium_chargeconst", w*GetMediumWeight());
  if(SameCharge(GetElectrons("ELECTRON16_POG_MEDIUM_DXYSIG")))     FillCutFlow("SS_Medium_dxysig",w*GetMediumWeight());

  if(PassTrigger(analysis_trigger) && SameCharge(GetElectrons("ELECTRON16_POG_MEDIUM_DXYSIG")))  FillCutFlow("SS_Medium_trigger",w*GetMediumWeight());

  /// HN ID
  if(SameCharge(GetElectrons("ELECTRON_HN_TIGHT"))) FillCutFlow("SS_anal_el", w);
  
  std::vector<snu::KJet> jets = GetJets("JET_HN");
  std::vector<snu::KJet> jets_pu = GetJets("JET_HN_PU");
  if(SameCharge(GetElectrons("ELECTRON_HN_TIGHT"))){
    if(jets.size() > 0)  FillCutFlow("ss_anal_el_jet", w*GetTightWeight());
    if(jets.size() > 1)  FillCutFlow("ss_anal_el_dijet", w*GetTightWeight());
    if(PassTrigger(analysis_trigger)){
      if(jets.size() > 1)  FillCutFlow("ss_trig_anal_el_dijet", w);
    }
    if(jets_pu.size() > 0)  FillCutFlow("ss_anal_el_pujet", w*GetTightWeight());
    if(jets_pu.size() > 1)  FillCutFlow("ss_anal_el_pudijet", w*GetTightWeight());
    if(PassTrigger(analysis_trigger)){
      if(jets_pu.size() > 1) FillCutFlow("ss_trig_anal_el_pudijet", w*GetTightWeight());
    }
  }
  
}

void HNDiElectron::GetOSSignalEfficiency(float w){
  
  // w = lumi weight * MC gen weight                                                                                                                                                                                                                                            
  //// Check efficiency of Selecting two same sign electrons                                                                                                                                                                                                                    
  
  if(OppositeCharge(GetElectrons("ELECTRON_NOCUT")))   FillCutFlow("OS_NoCut",w);
  if(OppositeCharge(GetElectrons("ELECTRON_PTETA")))   FillCutFlow("OS_PtEta",w);
  // TIGHT POG                                                                                                                                                                                                                                                                  
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_TIGHT_SC")))    FillCutFlow("OS_Tight",w*GetTightWeight());
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_TIGHT_D0"))) FillCutFlow("OS_Tight_d0veto",w*GetTightWeight());
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_TIGHT_RELISO")))     FillCutFlow("OS_Tight_reliso", w*GetTightWeight());
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_TIGHT_CHARGECONS"))) FillCutFlow("OS_Tight_chargeconst", w*GetTightWeight());
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_TIGHT_DXYSIG")))     FillCutFlow("OS_Tight_dxysig",w*GetTightWeight());

  TString analysis_trigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
  /// Trigger List (unprescaled)                                                                                                                                                                                                                                             
  std::vector<TString> triggerslist;
  triggerslist.push_back(analysis_trigger);
  
  if(PassTrigger(analysis_trigger) && OppositeCharge(GetElectrons("ELECTRON16_POG_TIGHT_DXYSIG")))  FillCutFlow("OS_Tight_trigger",w*GetTightWeight());


  // MEDIUM POG                                                                                                                                                                                                                                                                 
  
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_MEDIUM")))    FillCutFlow("OS_Medium",w*GetMediumWeight());
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_MEDIUM_D0"))) FillCutFlow("OS_Medium_d0veto",w*GetMediumWeight());
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_MEDIUM_RELISO")))     FillCutFlow("OS_Medium_reliso", w*GetMediumWeight());
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_MEDIUM_CHARGECONS"))) FillCutFlow("OS_Medium_chargeconst", w*GetMediumWeight());
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_MEDIUM_DXYSIG")))     FillCutFlow("OS_Medium_dxysig",w*GetMediumWeight());

  if(PassTrigger(analysis_trigger) && OppositeCharge(GetElectrons("ELECTRON16_POG_MEDIUM_DXYSIG")))  FillCutFlow("OS_Medium_trigger",w*GetMediumWeight());
  /// HN ID                                                                                                                                                                                                                                                                     
  if(OppositeCharge(GetElectrons("ELECTRON_HN_TIGHT"))) FillCutFlow("OS_anal_el", w);
  
  std::vector<snu::KJet> jets = GetJets("JET_HN");
  std::vector<snu::KJet> jets_pu = GetJets("JET_HN_PU");
  if(OppositeCharge(GetElectrons("ELECTRON_HN_TIGHT"))){
    if(jets.size() > 0)  FillCutFlow("os_anal_el_jet", w);
    if(jets.size() > 1)  FillCutFlow("os_anal_el_dijet", w);
    if(PassTrigger(analysis_trigger)){
	if(jets.size() > 1)  FillCutFlow("os_trig_anal_el_dijet", w);
    }
    if(jets_pu.size() > 0)  FillCutFlow("os_anal_el_pujet", w);
    if(jets_pu.size() > 1)  FillCutFlow("os_anal_el_pudijet", w);
    if(PassTrigger(analysis_trigger)){
      if(jets_pu.size() > 1)  FillCutFlow("os_trig_anal_el_pudijet", w);
    }
  }
}
										   
								   



void HNDiElectron::RunAnalysis(TString plottag, TString tightelid, TString vetoelid, TString looseelid){

  
  FillHist("NoCut" , 1., MCweight,  0. , 2., 2);
  FillHist("NoCut_w" , 1., weight,  0. , 2., 2);
    
  
  return;

  /// FillCutFlow(cut, weight) fills a basic TH1 called cutflow. It is used to check number of events passing different cuts
  /// The string cut must match a bin label in FillCutFlow function
  FillHist("GenWeight" , 1., MCweight,  0. , 2., 2);
  
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  
  //if(!PassBasicEventCuts()) return;     /// Initial event cuts  
  FillCutFlow("EventCut", weight);
  
  TString analysis_trigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
  /// Trigger List (unprescaled)
  std::vector<TString> triggerslist;
  triggerslist.push_back(analysis_trigger);

  std::vector<snu::KMuon> muonColl = GetMuons("MUON_HN_VETO"); // loose selection                                                                                                                                                                              
  
  if(muonColl.size() == 3){
    if(muonColl.at(0).Pt() > 20){
      if(muonColl.at(1).Pt() > 10){
        if(muonColl.at(2).Pt() > 10){
	  
	  bool pass_cut1(true);
          bool pass_cut2(true);
          for(std::vector<snu::KMuon>::iterator it = muonColl.begin(); it != muonColl.end(); it++){
            for(std::vector<snu::KMuon>::iterator it2 = it+1; it2 != muonColl.end(); it2++){
              if(it->Charge() != it2->Charge()){
                if(it->DeltaR(*it2) < 0.5)  pass_cut1 = false;
                float mass_dilep = (*it+*it2).M();
                if(fabs(mass_dilep - 90.) < 10.)  pass_cut2 = false;

              }
            }
          }
          if(pass_cut2){
            FillCLHist(trilephist, "TriMuon_mass", eventbase->GetEvent(), muonColl, GetElectrons("ELECTRON_HN_VETO"),GetJets("JET_HN"), weight);
            if(pass_cut1)          {
              FillCLHist(trilephist, "TriMuon_dr", eventbase->GetEvent(), muonColl, GetElectrons("ELECTRON_HN_VETO"),GetJets("JET_HN"), weight);

              if(muonColl.at(2).Pt() > 15)            FillCLHist(trilephist, "TriMuon_pt_dr", eventbase->GetEvent(), muonColl, GetElectrons("ELECTRON_HN_VETO"),GetJets("JET_HN"), weight);
              if(NBJet(GetJets("JET_HN")) == 0){
                FillCLHist(trilephist, "TriMuon_nobjet", eventbase->GetEvent(), muonColl, GetElectrons("ELECTRON_HN_VETO"),GetJets("JET_HN"), weight);
              }
              if(NBJet(GetJets("JET_HN"))> 1) FillCLHist(trilephist, "TriMuon_bjet", eventbase->GetEvent(), muonColl, GetElectrons("ELECTRON_HN_VETO"),GetJets("JET_HN"), weight);
            }
          }
        }
      }
    }
  }

  //if(!PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", GetElectrons("ELECTRON_HN_VETO"), prescale)) return;

  //return;
    //if(!PassTrigger(triggerslist, prescale)) return;
  

  if(PassTrigger(analysis_trigger)){
    FillCutFlow("TriggerCut", weight);
  }
  
  /// trigger_weight is for MC only: retruns 1 if data.
  /// Checks the luminosity of the trigger and returns weight that applied to 'weight' will correct for difference in luinosity of json file used in data
  float trigger_ps_weight= WeightByTrigger(analysis_trigger, TargetLumi);

  FillHist("PSWeight" , trigger_ps_weight, 1., 0. , 2., 200);

  
  if(PassTrigger(analysis_trigger)){
    FillCutFlow("TriggerWeight", weight);
  }

  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  FillCutFlow("VertexCut", weight);


  /// Use the number of vertices in the event to check effect of pileup reweighting
  numberVertices = eventbase->GetEvent().nVertices();   
  
  //////////////////////////////////////////////////////
  //////////// Select objetcs
  //////////////////////////////////////////////////////

  
  // Get loose muons for veto: Can call  POGSoft/POGLoose/POGMedium/POGTight/HNVeto/HNLoose/HNMedium/HNTight
  //  std::vector<snu::KMuon> muonColl = GetMuons(BaseSelection::MUON_HN_VETO); // loose selection
  
  /// Get tight jets : Can call NoLeptonVeto/Loose/Medium/Tight/HNJets
  std::vector<snu::KJet> jetColl_hn  = GetJets("JET_HN");// pt > 20 ; eta < 2.5; PFlep veto; NO pileup ID
  std::vector<snu::KJet> jetColl_nlv  = GetJets("JET_NOLEPTONVETO");
  std::vector<snu::KJet> jetColl_loose  = GetJets("JET_HN");

  std::vector<snu::KJet> JET_HN_TChannel  = GetJets("JET_HN");

  FillHist("Njets", jetColl_hn.size() ,weight, 0. , 5., 5);

  

  TString elid = "ELECTRON_POG_TIGHT";
  if(k_running_nonprompt) elid="ELECTRON_HN_FAKELOOSE_NOD0";

  std::vector<snu::KElectron> electronColl             = GetElectrons(true, false, elid);
  
  float weight_trigger_sf =  mcdata_correction->TriggerScaleFactor(electronColl, muonColl, analysis_trigger);
  FillHist("TriggerSFWeight" , weight_trigger_sf, 1., 0. , 2., 200);
  
  // Sets weight to weight if not running chargeflip bkg estimate or events are S
  //if(k_running_chargeflip) weight              *= WeightCFEvent(electronColl, k_running_chargeflip);
  

  std::vector<snu::KElectron> electronHNVetoColl   = GetElectrons("ELECTRON_HN_VETO");
  
  

  /// Correct MC for pileup   
  
  float pileup_reweight (1.);
  if (!k_isdata) {
    /// use silver or gold
    pileup_reweight = eventbase->GetEvent().PileUpWeight();
  }
  FillHist("PileupWeight" , pileup_reweight, 1.,  0. , 2., 200);
  
  float id_weight=1.;
  float reco_weight=1.;
  
  if(!isData){
    for(unsigned int iel=0; iel < electronColl.size(); iel++){
      id_weight*= mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronColl);
      
    }
  }
  FillHist("IDWeight" ,  id_weight,1.,  0. , 2., 200);
  FillHist("RecoWeight" ,  reco_weight, 1., 0. , 2., 200);

  if(!isData && !k_running_nonprompt){
    weight*= id_weight;
    weight*=reco_weight;
    weight*=pileup_reweight;
    weight*=trigger_ps_weight;
  }

  if(muonColl.size() > 0) return;
  if(electronHNVetoColl.size() > 2) return;



  if(k_running_nonprompt){
    weight=1.; /// In case... should not be needed
    //weight      *=  m_datadriven_bkg->Get_DataDrivenWeight_EE(false, electronColl,"ELECTRON_POG_TIGHT", "dijet");
  } 

  int nbjet = NBJet(jetColl_hn, snu::KJet::CSVv2, snu::KJet::Medium);
  float event_met = eventbase->GetEvent().PFMET();

  if(electronColl.size() == 2 ) {
   
    if(electronColl.at(0).Pt() > 20. && electronColl.at(1).Pt() > 15. ){
      
      FillCLHist(sighist_ee, "SIGNAL", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight);
      
      if(jetColl_loose.size() > 3 ) {
	FillCLHist(sighist_ee, "SIGNAL_4J", eventbase->GetEvent(), muonColl,electronColl,jetColl_loose, weight);
      } 
      


      if(SameCharge(electronColl))  {
	FillCLHist(sighist_ee, "SS_SIGNAL", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight);
	if(electronColl.at(0).GsfCtfScPixChargeConsistency() && electronColl.at(1).GsfCtfScPixChargeConsistency()){
	  if(electronColl.at(0).PassesConvVeto() && electronColl.at(1).PassesConvVeto()){
	    FillCLHist(sighist_ee, "SS_SIGNAL_CC", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight);
	    if(!Zcandidate(electronColl, 20., false)){

	      if(jetColl_loose.size() > 3 ) {
		bool has_forward_jet(false), has_back_jet(false);
		for(unsigned int ij = 0 ; ij < jetColl_loose.size(); ij++){
		  if(jetColl_loose.at(ij).Eta() > 2.5) has_forward_jet=true;
		  if(jetColl_loose.at(ij).Eta() < -2.5) has_back_jet=true;
		}
		if(has_forward_jet && has_back_jet)
		  FillCLHist(sighist_ee, "SS_SIGNAL_4J", eventbase->GetEvent(), muonColl,electronColl,jetColl_loose, weight);
	      }
	      
	      FillCLHist(sighist_ee, "SS_SIGNAL_noZ", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight);

	      if(electronColl.at(1).Pt() < 20.)  FillCLHist(sighist_ee, "SS_SIGNAL_LowPt", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight);

	      if(electronColl.at(0).IsEBFiducial()   && electronColl.at(1).IsEBFiducial())       FillCLHist(sighist_ee, "SS_SIGNAL_BB", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight);
	      if(electronColl.at(0).IsEEFiducial() && electronColl.at(1).IsEEFiducial())       FillCLHist(sighist_ee, "SS_SIGNAL_EE", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight);
	      
	      if(jetColl_hn.size() == 1 && (GetDiLepMass(electronColl) > 100. )) FillCLHist(sighist_ee, "SS_SIGNAL_1Jet", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight);
	      if(jetColl_hn.size() > 1 ) FillCLHist(sighist_ee, "SS_SIGNAL_Presel", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight);
	      
	      if(jetColl_hn.size() > 3 ) {
		bool has_forward_jet(false), has_back_jet(false);
		for(unsigned int ij = 0 ; ij < jetColl_hn.size(); ij++){
		  if(jetColl_hn.at(ij).Eta() > 1.5) has_forward_jet=true;
		  if(jetColl_hn.at(ij).Eta() < -1.5) has_back_jet=true;
		}
		if(has_forward_jet && has_back_jet){
		  FillCLHist(sighist_ee, "TChannel", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight);
		  FillHist("SigTchannel" , 1., MCweight,  0. , 2., 2);
		}
		
	      }



	      if(NBJet(jetColl_hn, snu::KJet::CSVv2, snu::KJet::Medium) == 0){
		if(electronColl.at(0).IsEBFiducial()   && electronColl.at(1).IsEBFiducial())       FillCLHist(sighist_ee, "SS_SIGNAL_BB_noB", eventbase->GetEvent(),   muonColl,electronColl,jetColl_hn, weight);
		if(electronColl.at(0).IsEEFiducial() && electronColl.at(1).IsEEFiducial())       FillCLHist(sighist_ee, "SS_SIGNAL_EE_noB", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight);
	      }
	      else{
		if(electronColl.at(0).IsEBFiducial()   && electronColl.at(1).IsEBFiducial())       FillCLHist(sighist_ee, "SS_SIGNAL_BB_B", eventbase->GetEvent(),   muonColl,electronColl,jetColl_hn, weight);
                if(electronColl.at(0).IsEEFiducial() && electronColl.at(1).IsEEFiducial())       FillCLHist(sighist_ee, "SS_SIGNAL_EE_B", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight);
              }
	      
	      
	      if(LowMassCheckSignalRegion(electronColl, jetColl_hn, k_running_chargeflip) ) {
		if((nbjet==0) && ( (eventbase->GetEvent().PFMET()  < 30.))){
		  FillCLHist(sighist_ee, "LowMassRegion", eventbase->GetEvent(), muonColl,   electronColl,jetColl_hn, weight);
		} 
		if ((nbjet!=0 && (eventbase->GetEvent().PFMET()  < 30.)) || ((nbjet == 0) && (eventbase->GetEvent().PFMET() > 50.))){
		  FillCLHist(sighist_ee, "LowMassRegionCR", eventbase->GetEvent(), muonColl,   electronColl,jetColl_hn, weight);
		  
		}
	      }
	      if(MidMassCheckSignalRegion(electronColl, jetColl_hn, k_running_chargeflip) ) {
                if((nbjet==0) && ( (eventbase->GetEvent().PFMET()  < 30.))){
                  FillCLHist(sighist_ee, "MediumMassRegion", eventbase->GetEvent(), muonColl,   electronColl,jetColl_hn, weight);
                }
		if ((nbjet!=0 && (eventbase->GetEvent().PFMET()  < 30.)) || ((nbjet == 0) && (eventbase->GetEvent().PFMET() > 50.))){
                  FillCLHist(sighist_ee, "MediumMassRegionCR", eventbase->GetEvent(), muonColl,   electronColl,jetColl_hn, weight);
		}
	      } 
	      if(HighMassCheckSignalRegion(electronColl, jetColl_hn, k_running_chargeflip) ) {
		if((nbjet==0) && ( (eventbase->GetEvent().PFMET()  < 30.))){
                  FillCLHist(sighist_ee, "HighMassRegion", eventbase->GetEvent(), muonColl,   electronColl,jetColl_hn, weight);
                }
              }
	      
	        
	      std::vector<float> masscuts40;
	      masscuts40.push_back(20.);  masscuts40.push_back(15.);   /// pt1 / pt2
	      masscuts40.push_back(10.);  masscuts40.push_back(60.);   /// m(ee)
	      masscuts40.push_back(0.);   masscuts40.push_back(120.);  /// m(jj)
	      masscuts40.push_back(80.);  masscuts40.push_back(160.);  /// m(eejj)
	      masscuts40.push_back(30.);                               /// MET
	      masscuts40.push_back(20.);                               /// jet 1 pt
	      masscuts40.push_back(0.);  masscuts40.push_back(120.);  /// m(e2jj)

	      
	      std::vector<float> masscuts100;
	      masscuts100.push_back(20.);  masscuts100.push_back(15.);
	      masscuts100.push_back(15.);  masscuts100.push_back(10000.);
	      masscuts100.push_back(50.);  masscuts100.push_back(110.);
	      masscuts100.push_back(120.); masscuts100.push_back(10000.);
	      masscuts100.push_back(35.);
	      masscuts100.push_back(30.);
	      masscuts100.push_back(80.);  masscuts100.push_back(120.);


	      std::vector<float> masscuts500;
	      masscuts500.push_back(120);  masscuts500.push_back(40.);
	      masscuts500.push_back(15.);  masscuts500.push_back(10000.);
	      masscuts500.push_back(50.);  masscuts500.push_back(110.);
	      masscuts500.push_back(360.); masscuts500.push_back(10000.);
	      masscuts500.push_back(35.);
	      masscuts500.push_back(40.);
	      masscuts500.push_back(0.);  masscuts500.push_back(10000.);


	      std::vector<float> masscuts1500;
	      masscuts1500.push_back(120);  masscuts1500.push_back(40.);
	      masscuts1500.push_back(15.);  masscuts1500.push_back(10000.);
	      masscuts1500.push_back(50.);  masscuts1500.push_back(110.);
	      masscuts1500.push_back(360.); masscuts1500.push_back(10000.);
	      masscuts1500.push_back(35.);
	      masscuts1500.push_back(40.);
	      masscuts1500.push_back(0.);  masscuts1500.push_back(10000.);
	      
	      std::map<TString, std::vector<float> > sysymap;
	      sysymap["40"] =  masscuts40; 
	      sysymap["100"] =  masscuts100;
	      sysymap["500"] =  masscuts500;
	      sysymap["1500"] =  masscuts1500;
	      std::map<TString, TString > sysymap2;
	      sysymap2["_default"] = "default";



	      std::map< TString,  std::pair< std::vector<snu::KElectron> , float> > leptonmap;
	      leptonmap[""]  = make_pair(electronColl, weight);
  

	      for( std::map< TString,  std::pair< std::vector<snu::KElectron> , float> >::iterator lepmapit = leptonmap.begin(); lepmapit != leptonmap.end(); lepmapit++){
		for( std::map<TString, std::vector<float> >::iterator mapit = sysymap.begin(); mapit != sysymap.end(); mapit++){
		  for( std::map<TString, TString >::iterator mapit2 = sysymap2.begin(); mapit2 != sysymap2.end(); mapit2++){
		    TString label = lepmapit->first + mapit->first+ mapit2->first + "MassRegion_limithist";
		    TString jetlabel = "";
		    if(mapit->first.Contains("fg") && (mapit2->second != "default")) continue;

		    if(OptMassCheckSignalRegion(lepmapit->second.first, jetColl_hn,event_met, k_running_chargeflip , mapit->second, mapit2->second + jetlabel) ){
		      if(mapit2->second == "default"){
			if(!mapit->first.Contains("fg")){
			  FillCLHist(sighist_ee, (mapit->first + lepmapit->first+"MassRegion").Data(), eventbase->GetEvent(), muonColl,lepmapit->second.first,jetColl_hn, lepmapit->second.second);

			  
			}
		      }
		    }
		  }
		}
	      }

	      
	    }
	  }
	}
      }
      else     FillCLHist(sighist_ee, "OS_SIGNAL", eventbase->GetEvent(), muonColl,electronColl,jetColl_hn, weight);
    }
  }



    

  return;
}// End of execute event loop
  

void HNDiElectron::SignalValidation(){
  
  // Check jet properties                                                                                                                                                                                                                                                                                                                                                         
  FillCLHist(sighist_ee, "SIGNALVALIDATION_EE", eventbase->GetEvent(),  GetMuons("MUON_NOCUT"), GetElectrons("ELECTRON_NOCUT"),GetJets("JET_HN"), weight);
  FillCLHist(sighist_ee, "SIGNALVALIDATION_TIGHT_EE", eventbase->GetEvent(),  GetMuons("MUON_NOCUT"), GetElectrons("ELECTRON_POG_TIGHT"),GetJets("JET_HN"), weight);
  FillCLHist(sighist_mm, "SIGNALVALIDATION_MM", eventbase->GetEvent(),  GetMuons("MUON_NOCUT"), GetElectrons("ELECTRON_NOCUT"),GetJets("JET_HN"), weight);
  if(GetJets("JET_HN").size() > 1){
    if(GetElectrons("ELECTRON_POG_TIGHT").size() == 2){
      FillCLHist(sighist_ee, "SIGNALVALIDATION_EE_PRESEL", eventbase->GetEvent(),  GetMuons("MUON_NOCUT"), GetElectrons("ELECTRON_POG_TIGHT"),GetJets("JET_HN"), weight);
    }
    if(GetMuons("MUON_POG_TIGHT").size() == 2){
      FillCLHist(sighist_mm, "SIGNALVALIDATION_MM_PRESEL", eventbase->GetEvent(),  GetMuons("MUON_POG_TIGHT"), GetElectrons("ELECTRON_NOCUT"),GetJets("JET_HN"), weight);
    }
  }
  

}
 

void HNDiElectron::GetTriggEfficiency(){
  //ListTriggersAvailable();                                                                                                                                                                                                                                                                                                                                                      
  vector<int> pt1;
  pt1.push_back(18);
  pt1.push_back(35);
  pt1.push_back(25);
  pt1.push_back(25);
  pt1.push_back(25);
  pt1.push_back(30);
  pt1.push_back(30);
  pt1.push_back(20);
  pt1.push_back(25);
  pt1.push_back(30);
  pt1.push_back(120);
  pt1.push_back(20);
  pt1.push_back(25);
  pt1.push_back(25);
  vector<int>pt2;
  pt2.push_back(10);
  pt2.push_back(35);
  pt2.push_back(25);
  pt2.push_back(10);
  pt2.push_back(10);
  pt2.push_back(10);
  pt2.push_back(10);
  pt2.push_back(15);
  pt2.push_back(15);
  pt2.push_back(10);
  pt2.push_back(10);
  pt2.push_back(15);
  pt2.push_back(15);
  pt2.push_back(15);


  std::vector<TString> lists_triggers;
  lists_triggers.push_back("HLT_Ele15_IsoVVVL_PFHT600_v");
  lists_triggers.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v");
  lists_triggers.push_back("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v");
  lists_triggers.push_back("HLT_Ele22_eta2p1_WPLoose_Gsf_v");
  lists_triggers.push_back("HLT_Ele23_WPLoose_Gsf_v");
  lists_triggers.push_back("HLT_Ele27_WPLoose_Gsf_v");
  lists_triggers.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_v");
  lists_triggers.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v");
  lists_triggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
  lists_triggers.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_HT200");
  lists_triggers.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");
  lists_triggers.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  lists_triggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  lists_triggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded_v1");
  for(unsigned int i=0; i < lists_triggers.size(); i++){
    FillTriggerEfficiency(lists_triggers.at(i), weight, "denominator_nojet", lists_triggers );
  }

  std::vector<snu::KJet> jets=GetJets("JET_HN");
  std::vector<snu::KElectron> electrons= GetElectrons("ELECTRON_POG_TIGHT");
  if(jets.size() > 1){
    for(unsigned int i=0; i < lists_triggers.size(); i++){
      TString trig=lists_triggers.at(i);
      FillTriggerEfficiency(lists_triggers.at(i), weight, "denominator", lists_triggers );
      if(PassTrigger(trig))  {
	FillTriggerEfficiency(lists_triggers.at(i), weight, "numerator",lists_triggers );
	
	if(electrons.size() ==2) {
	  
	  FillTriggerEfficiency(lists_triggers.at(i), weight, "numerator_dimuon",lists_triggers );
	  if(electrons.at(0).Pt() > pt1.at(i) && electrons.at(1).Pt() > pt2.at(i))  FillTriggerEfficiency(lists_triggers.at(i), weight, "numerator_dimuon_pt",lists_triggers );
	}
      }
    }
  }
}

void HNDiElectron::FillTriggerEfficiency(TString cut, float weight, TString label, std::vector<TString> list){

  if(GetHist("TriggerEfficiency_" + label)) {
    GetHist("TriggerEfficiency_"+label)->Fill(cut,weight);

  }
  else{
    int ntrig = list.size();
    AnalyzerCore::MakeHistograms("TriggerEfficiency_"+label,ntrig,0.,float(ntrig));

    for(unsigned int it=0; it < list.size(); it++){
      GetHist("TriggerEfficiency_"+label)->GetXaxis()->SetBinLabel(it+1,list.at(it));
    }
  }
  
}


 bool HNDiElectron::CheckSignalRegion( bool isss,  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, std::vector<snu::KJet> alljets,TString name, float w){

  if(electrons.size() != 2 ) return false ;
  // Set by trigger
  if(electrons.at(0).Pt() < 25.) return false;
  if(electrons.at(1).Pt() < 10.) return false;
  if(isss&&!SameCharge(electrons)) return false;

  if(!isss&&SameCharge(electrons)) return false;
  if(jets.size() < 2) return false;
  snu::KParticle ee = electrons.at(0) + electrons.at(1);
  if(ee.M()  < 10.) return false;

  bool debug=false;
  if(name.Contains("Low")){
    if(ee.M()  > 70.) {if(debug)cout << "Fail mee  " << endl; return false;}
  }
  
  
  float dijetmass_tmp=999.;
  float dijetmass=9990000.;
  int m=-999;
  int n=-999;
  for(UInt_t emme=0; emme<jets.size(); emme++){
    for(UInt_t enne=1; enne<jets.size(); enne++) {
      if(emme == enne) continue;
      dijetmass_tmp = (jets[emme]+jets[enne]).M();
      if ( fabs(dijetmass_tmp-80.4) < fabs(dijetmass-80.4) ) {
        dijetmass = dijetmass_tmp;
        m = emme;
        n = enne;
      }
    }
  }
  snu::KParticle jj = jets.at(m) + jets.at(n) ;

  float dPhi = fabs(TVector2::Phi_mpi_pi(jets[m].Phi() - jets[n].Phi()));
  float contramass=2*jets[m].Pt()*jets[n].Pt()*(1+cos(dPhi));
  contramass=sqrt(contramass);

  if(name.Contains("Low")){
    if(contramass > 100) return false;
    if((jets[0]+jets[1]).M() > 250.) return false;
    //    if((jets[0] + electrons[0] + electrons[1]).M() > 250.)  return false;
  }
  else{
    //if(jj.M() > 120.) {if(debug)cout << "Fail mjj  " << endl; return false;}
    //if(jj.M() < 50.) {if(debug)cout << "Fail mjj  " << endl; return false;}

  }


  if(name.Contains("Low")){
    if(electrons.at(0).DeltaR(electrons.at(1)) > 3.5) return false;

  }

  if((ee.M() > 80.) && (ee.M() < 100.) ) return false;



  float ST = electrons[0].Pt() + electrons[1].Pt();
  float looseST = electrons[0].Pt() + electrons[1].Pt();
  float HT=0.;
  for(unsigned int ij=0; ij <alljets.size(); ij++){
    ST+= alljets[ij].Pt();

  }
  for(unsigned int ij=0; ij <jets.size(); ij++){

    HT+= jets[ij].Pt();

  }

  float LT = electrons[0].Pt() + electrons[1].Pt();
  if(name.Contains("Low")){
    if((pow(eventbase->GetEvent().PFMET(),2.)/ ST) > 15.) {if(debug)cout << "Fail met  " << endl; return false;}
    if(eventbase->GetEvent().PFMET() > 80) return false;
  }
  else if(name.Contains("High")){
    if((pow(eventbase->GetEvent().PFMET(),2.)/ ST) > 8.) {if(debug)cout << "Fail met  " << endl; return false;}
  }
  else{
    if((pow(eventbase->GetEvent().PFMET(),2.)/ ST) > 12.5) {if(debug)cout << "Fail met  " << endl; return false;}
    if(eventbase->GetEvent().PFMET() > 80) return false;

  }

  if(name.Contains("Low")){
    if(ST > 450.)  return false;
  }
  else   if(name.Contains("High")){
    snu::KParticle lljj = jets.at(m) + jets.at(n) + electrons.at(0) + electrons.at(1);

    if(ST < 400.)  return false;
    if(lljj.M() < 200.) return false;
    if(electrons.at(0).Pt() < 50.) return false;
    if(electrons.at(1).Pt() < 25.) return false;
  }

  if(name.Contains("High")){
    //if((electrons.at(0).Pt() + electrons.at(1).Pt() ) < 50.) return false;
    if (ST < 400. ) return false;

    float mindRjj=999.;
    for(unsigned int ij=0; ij <jets.size()-1; ij++){
      for(unsigned int ij2=ij+1; ij2 <jets.size(); ij2++){
        if(jets.at(ij).DeltaR(jets.at(ij2)) < mindRjj) mindRjj=jets.at(ij).DeltaR(jets.at(ij2)) ;
      }
    }
  }

  int nbjet=0;
  for(unsigned int ij=0; ij <alljets.size(); ij++){
    if( alljets[ij].IsBTagged(  snu::KJet::CSVv2, snu::KJet::Medium))  nbjet++;
  }
  if(nbjet > 0) {if(debug)cout << "Fail nbjet  " << endl; return false;}


  if(debug)cout << "PASSES ID" << endl;
  return true;


  
  

  return true;

}

/*
bool HNDiElectron::CheckSignalRegion( bool isss,  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, TString name, float w){

  if(electrons.size() != 2 ) return false ;
  // Set by trigger                                                                                                                                                             
  if(electrons.at(0).Pt() < 25.) return false;
  if(electrons.at(1).Pt() < 15.) return false;
  if(isss&&!SameCharge(electrons)) return false;

  return true;
  if(!isss&&SameCharge(electrons)) return false;
  if(jets.size() < 2) return false;
  snu::KParticle ee = electrons.at(0) + electrons.at(1);
  if(ee.M()  < 10.) return false;

  snu::KParticle jj = jets.at(0) + jets.at(1) ;
  if(jj.M() > 120.) return false;
  if(jj.M() < 40.) return false;

  if(electrons.at(0).DeltaR(electrons.at(1)) > 3.5) return false;
  if(eventbase->GetEvent().SumET() < 200.) return false;;

  if((ee.M() > 80.) && (ee.M() < 100.) ) return false;
  float ST = electrons[0].Pt() + electrons[1].Pt();
  for(unsigned int ij=0; ij <jets.size(); ij++){
    ST+= jets[ij].Pt();
  }
  if((eventbase->GetEvent().PFMET()/pow(2., ST))  > 20.) return false;

  int nbjet=0;
  for(unsigned int ij=0; ij <jets.size(); ij++){
    if(IsBTagged(jets.at(ij), snu::KJet::CSVv2, snu::KJet::Medium, GetPeriod())) nbjet++;
  }
  if(nbjet > 0) return false;


  return true;

}

 */



bool HNDiElectron::OptMassCheckSignalRegion(std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, float evmet, bool runchargeflip , std::vector<float> cuts , TString additional_option){
  
  return OptMassCheckSignalRegion(electrons, jets, 1, evmet, runchargeflip, cuts, additional_option);
}

bool HNDiElectron::OptMassCheckSignalRegion(std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, int btagwp,  float evmet, bool runchargeflip , std::vector<float> cuts, TString additional_option){
  
  if(additional_option.Contains("btag_t")) btagwp = 2; 
  if(additional_option.Contains("btag_l")) btagwp = 0; 

  bool removeclosejet = false;
  if(additional_option.Contains("closejet")) removeclosejet = true;
  
  bool cut_on_me2jj=true;
  if(additional_option.Contains("noMe2jj")) cut_on_me2jj= false;
  
  bool cut_on_ee=true;
  if(additional_option.Contains("noeeupper")) cut_on_ee= false;

  bool cut_on_eejj=true;
  if(additional_option.Contains("noeejjupper")) cut_on_eejj= false;

  bool cut_on_pt=true;
  if(additional_option.Contains("pt10")) cut_on_pt= false;

  
  int nbjet=0;
  snu::KJet::WORKING_POINT wp;
  if(btagwp==0) wp=snu::KJet::Loose;
  if(btagwp==1) wp=snu::KJet::Medium;
  if(btagwp==2) wp=snu::KJet::Tight;

  for(unsigned int ij=0; ij <jets.size(); ij++){
    if(IsBTagged(jets.at(ij), snu::KJet::CSVv2, snu::KJet::Medium, GetPeriod())) nbjet++;
  }

  if(additional_option.Contains("lowmass")) return (nbjet == 0 && (evmet < 30.) && LowMassCheckSignalRegion( electrons,  jets, runchargeflip));
  if(additional_option.Contains("midmass"))  return (nbjet == 0 && (evmet < 35.) && MidMassCheckSignalRegion( electrons,  jets, runchargeflip));
  if(additional_option.Contains("highmass")) return (nbjet == 0 && (evmet < 35.) && HighMassCheckSignalRegion( electrons,  jets, runchargeflip));
  

  if(additional_option.Contains("basic")){
    if(electrons.size() != 2) return false;
    if(electrons.at(0).Pt() < 20.) return false;
    if(electrons.at(1).Pt() < 15.) return false;
    if(jets.size() < 2) return false;
    if(jets.at(0).Pt() < cuts.at(9)) return false;
    if(!runchargeflip){
      if(!SameCharge(electrons)) return false;
    }
    float wmassjj= 100000.;
    int indexj1=0;
    int indexj2=0;

    for(unsigned int ij=0; ij < jets.size()-1; ij++){
      for(unsigned int ij2=ij+1; ij2 < jets.size(); ij2++){
	snu::KParticle jjtmp = jets.at(ij) + jets.at(ij2) ;
	if(fabs(jjtmp.M() - 80.4) < wmassjj) {
	  wmassjj = fabs(jjtmp.M() - 80.4);
	  indexj1=ij;
	  indexj2=ij2;
	}
      }
    }

    snu::KParticle jj = jets.at(indexj1) + jets.at(indexj2);
    if(jj.M() > 110.) return false;
    if(jj.M() < 50.) return false;
    if(evmet > 35.) return false;
    
    return true;
  }
  
  if(electrons.size() != 2 ) return false;
  
  if(electrons.at(0).Pt() < cuts.at(0)) return false;
  if(!cut_on_pt){
    if(electrons.at(1).Pt() < 10.) return false;
  }
  else if(electrons.at(1).Pt() < cuts.at(1)) return false;
  
  if(!runchargeflip){
    if(!SameCharge(electrons)) return false;
  }
  if(jets.size() < 2) return false;
  
  
  snu::KParticle ee = electrons.at(0) + electrons.at(1);
  float mass=ee.M();
  
  if(mass < cuts.at(2)) return false;

  if(!cut_on_ee) {
    if(mass > 100000000.) return false;
  }
  else if(mass > cuts.at(3)) return false;
  if(fabs(mass-90.) < 10.) return false;

  if(jets.at(0).Pt() < cuts.at(9)) return false;
   
  float wmassjj= 100000.;
  int indexj1=0;
  int indexj2=0;

  for(unsigned int ij=0; ij < jets.size()-1; ij++){
    for(unsigned int ij2=ij+1; ij2 < jets.size(); ij2++){
      snu::KParticle jjtmp = jets.at(ij) + jets.at(ij2) ;
      if(fabs(jjtmp.M() - 80.4) < wmassjj) {
        wmassjj = fabs(jjtmp.M() - 80.4);
        indexj1=ij;
        indexj2=ij2;
      }
    }
  }
  snu::KParticle jj = jets.at(indexj1) + jets.at(indexj2);
  if(jj.M() > cuts.at(5)) return false;
  if(jj.M() < cuts.at(4)) return false;

  snu::KParticle eejj = electrons.at(0) + electrons.at(1)+ jets.at(indexj1) + jets.at(indexj2) ;

  if(!cut_on_eejj) {
    if(eejj.M()  > 200.) return false;
  }
  else  if(eejj.M()  > cuts.at(7)) return false;
  if(eejj.M()  < cuts.at(6)) return false;

  snu::KParticle e1jj = electrons.at(0)+ jets.at(indexj1) + jets.at(indexj2) ;
  //if(e1jj.M() > cuts.at(13)) return false;
  //if(e1jj.M() < cuts.at(12)) return false;

  snu::KParticle e2jj = electrons.at(1)+ jets.at(indexj1) + jets.at(indexj2) ;
  
  if(cut_on_me2jj){
    if(e2jj.M() > cuts.at(11)) return false;
    if(e2jj.M() < cuts.at(10)) return false;
  }

  if(removeclosejet){
    if(HasCloseBJet( electrons.at(0)) ) return false;
    if(HasCloseBJet( electrons.at(1)) ) return false;
  }

  if(evmet > cuts.at(8)) return false;

  if(nbjet > 0) return false;

  return true;

}


bool HNDiElectron::LowMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runchargeflip){
  if(electrons.size() != 2 ) return false;
  if(electrons.at(0).Pt() < 20.) return false;
  if(electrons.at(1).Pt() < 15.) return false;
  if(!runchargeflip){
    if(!SameCharge(electrons)) return false;
  }
  if(jets.size() < 2) return false;
  snu::KParticle ee = electrons.at(0) + electrons.at(1);
  if(ee.M()  < 10.) return false;
  if(fabs(ee.M()-90.)  < 10.) return false;
  
  float wmassjj= 100000.;
  int indexj1=0;
  int indexj2=0;
  for(unsigned int ij=0; ij < jets.size()-1; ij++){
    for(unsigned int ij2=ij+1; ij2 < jets.size(); ij2++){
      snu::KParticle jjtmp = jets.at(ij) + jets.at(ij2);
      if(fabs(jjtmp.M() - 80.4) < wmassjj) {
        wmassjj = fabs(jjtmp.M() - 80.4);
        indexj1=ij;
        indexj2=ij2;
      }
    }
  }
  snu::KParticle jj = jets.at(indexj1) + jets.at(indexj2);
  if(jj.M() > 120.) return false;

  snu::KParticle eejj = electrons.at(0) + electrons.at(1)+ jets.at(indexj1) + jets.at(indexj2) ;
  if(eejj.M() > 200.)return false;

  //if(eventbase->GetEvent().PFMET() > 30.) return false;
  
  int nbjet=0;
  for(unsigned int ij=0; ij <jets.size(); ij++){
    if(IsBTagged(jets.at(ij), snu::KJet::CSVv2, snu::KJet::Medium, GetPeriod())) nbjet++;
  }
  //if(nbjet > 0) return false;

  return true;
  
}

bool HNDiElectron::MidMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runchargeflip){
  
  if(electrons.size() != 2 ) return false;
  if(electrons.at(0).Pt() < 20.) return false;
  if(electrons.at(1).Pt() < 15.) return false;
  if(!runchargeflip){
    if(!SameCharge(electrons)) return false;
  }
  if(jets.size() < 2) return false;
  snu::KParticle ee = electrons.at(0) + electrons.at(1);
  if(ee.M()  < 15.) return false;
  if(fabs(ee.M()-90.)  < 10.) return false;

  if(jets.at(0).Pt() < 30.) return false;

  float wmassjj= 100000.;
  int indexj1=0;
  int indexj2=0;
  for(unsigned int ij=0; ij < jets.size()-1; ij++){
    for(unsigned int ij2=ij+1; ij2 < jets.size(); ij2++){
      snu::KParticle jjtmp = jets.at(ij) + jets.at(ij2) ;
      if(fabs(jjtmp.M() - 80.4) < wmassjj) {
        wmassjj = fabs(jjtmp.M() - 80.4);
	indexj1=ij;
        indexj2=ij2;
      }
    }
  }
  snu::KParticle jj = jets.at(indexj1) + jets.at(indexj2);
  if(jj.M() > 110.) return false;
  if(jj.M() < 50.) return false;

  
  //if(eventbase->GetEvent().PFMET() > 35.) return false;

  int nbjet=0;
  for(unsigned int ij=0; ij <jets.size(); ij++){
    if(IsBTagged(jets.at(ij), snu::KJet::CSVv2, snu::KJet::Medium, GetPeriod())) nbjet++;
  }
  //if(nbjet > 0) return false;

  return true;
}



	
	

bool HNDiElectron::HighMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runchargeflip){
  if(electrons.size() != 2 ) return false;
  if(electrons.at(0).Pt() < 20.) return false;
  if(electrons.at(1).Pt() < 15.) return false;
  if(!runchargeflip){
    if(!SameCharge(electrons)) return false;
  }
  if(jets.size() < 2) return false;
  snu::KParticle ee = electrons.at(0) + electrons.at(1);
  if(ee.M()  < 15.) return false;
  if(fabs(ee.M()-90.)  < 10.) return false;

  float wmassjj= 100000.;
  int indexj1=0;
  int indexj2=0;
  for(unsigned int ij=0; ij < jets.size()-1; ij++){
    for(unsigned int ij2=ij+1; ij2 < jets.size(); ij2++){
      snu::KParticle jjtmp = jets.at(ij) + jets.at(ij2) ;
      if(fabs(jjtmp.M() - 80.4) < wmassjj) {
        wmassjj = fabs(jjtmp.M() - 80.4);
        indexj1=ij;
        indexj2=ij2;
      }
    }
  }
  snu::KParticle jj = jets.at(indexj1) + jets.at(indexj2);
  if(jj.M() > 110.) return false;
  if(jj.M() < 50.) return false;

  //  if(eventbase->GetEvent().PFMET() > 35.) return false;
  if(jets.at(0).Pt() < 40.) return false;
  
  int nbjet=0;
  for(unsigned int ij=0; ij <jets.size(); ij++){
    if(IsBTagged(jets.at(ij), snu::KJet::CSVv2, snu::KJet::Medium, GetPeriod())) nbjet++;
  }
  //  if(nbjet > 0) return false;

  return true;
}




void HNDiElectron::counter(TString cut, float w){

  w=1;
  map<TString,float>::iterator itmapcounter = mapcounter.find(cut) ;
  if (itmapcounter == mapcounter.end()){
    mapcounter[cut] = w;
  }
  else{
    float sum = itmapcounter->second;
    mapcounter[cut] = sum+w;
  }

  FillCutFlow(cut,w);

}


void HNDiElectron::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger << DEBUG << "END OF CYCLE: isdata=" << isData <<  LQLogger::endmsg;

  for ( map<TString,float>::iterator itmapcounter = mapcounter.begin(); itmapcounter != mapcounter.end(); itmapcounter++){
    cout << itmapcounter->first << " has count "<< itmapcounter->second << endl;
  }

}


void HNDiElectron::BeginCycle() throw( LQError ){
  
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

HNDiElectron::~HNDiElectron() {
  
  Message("In HNDiElectron Destructor" , INFO);
  
}



void HNDiElectron::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void HNDiElectron::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNDiElectronCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HNDiElectron::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



