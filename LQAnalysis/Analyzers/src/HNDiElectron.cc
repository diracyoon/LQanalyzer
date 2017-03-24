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

//// Needed to allow inheritance for use in LQCore/core classes
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
     MakeCleverHistograms(sighist_ee,"OSDiElectron");
     MakeCleverHistograms(sighist_ee,"OSDiElectronDiJet");
     MakeCleverHistograms(sighist_ee,"SSDiElectron");
     MakeCleverHistograms(sighist_ee,"SingleElectron_Wregion");
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
  
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return;
  
  
  if(!PassMETFilter()) return;     /// Initial event cuts :                                                                                                                                                                                                                 
  if(true){
    TString analysis_trigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
    std::vector<TString> triggerslist;
    triggerslist.push_back(analysis_trigger);
    if(PassTrigger(analysis_trigger)){
      std::vector<snu::KElectron> electronTightColl=GetElectrons(TString("ELECTRON_POG_TIGHT"),15., 2.5);
          
      if(OppositeCharge(electronTightColl)){
	if(electronTightColl.at(0).Pt() > 30. && electronTightColl.at(1).Pt() > 20.){
	  if(GetDiLepMass(electronTightColl) < 120. && GetDiLepMass(electronTightColl)  > 60. ){
	    if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 75., 75);
	    else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 75., 75);
	    if(isData) FillHist("Nint_nocut_data", eventbase->GetEvent().PileUpInteractionsTrue(),weight, 0. , 75., 75);
	    if(isData) FillHist("Nint_nocut_mc", eventbase->GetEvent().PileUpInteractionsTrue(),weight, 0. , 75., 75);
	  }
	}
      }
    }
    return;
  }
  FillHist("PFMet_uncleaned" , eventbase->GetEvent().PFMET(), MCweight,  0. , 1000., 1000);

  FillHist("PFMet" , eventbase->GetEvent().PFMET(), MCweight,  0. , 1000., 1000);
  FillHist("PFMetPhi" , eventbase->GetEvent().METPhi(snu::KEvent::pfmet), MCweight,  -3. , 3., 100);
  FillHist("PFMet_t1" , eventbase->GetEvent().PFMETType1(), MCweight,  0. , 1000., 1000);
  FillHist("PFMetPhi_t1" , eventbase->GetEvent().PFMETType1Phi(), MCweight,  -3. , 3., 100);
  FillHist("PFMet_unsmeared" , eventbase->GetEvent().PFMETUnSmearedType1Pt(), MCweight,  0. , 1000., 1000);
  FillHist("PFMetPhi_unsmeared" , eventbase->GetEvent().PFMETUnSmearedType1Phi(), MCweight,  -3. , 3., 100);

  if(eventbase->GetEvent().PFMET() >  70.){
    FillHist("PFMet_highmet" , eventbase->GetEvent().PFMET(), MCweight,  0. , 1000., 1000);
    FillHist("PFMetPhi_highmet" , eventbase->GetEvent().METPhi(snu::KEvent::pfmet), MCweight,  -3. , 3., 100);
    FillHist("PFMet_t1_highmet" , eventbase->GetEvent().PFMETType1(), MCweight,  0. , 1000., 1000);
    FillHist("PFMetPhi_t1_highmet" , eventbase->GetEvent().PFMETType1Phi(), MCweight,  -3. , 3., 100);
    FillHist("PFMet_unsmeared_highmet" , eventbase->GetEvent().PFMETUnSmearedType1Pt(), MCweight,  0. , 1000., 1000);
    FillHist("PFMetPhi_unsmeared_highmet" , eventbase->GetEvent().PFMETUnSmearedType1Phi(), MCweight,  -3. , 3., 100);
  }


  std::vector<snu::KElectron> electronVetoColl=GetElectrons("ELECTRON_HN_VETO");
  std::vector<snu::KElectron> electronTightColl=GetElectrons(TString("ELECTRON_POG_TIGHT"),15., 2.5);
  std::vector<snu::KMuon> muonVetoColl=GetMuons("MUON_HN_VETO");
  std::vector<snu::KJet> jets =  GetJets("JET_HN");

  if(functionality == HNDiElectron::VALIDATION){

    /// make validation plots for dielectron events

    float trigger_sf(1.);
    float id_iso_sf(1.);
    float trigger_ps(1.);
    float reco_weight=1.;
    float puweight(1.);
    float ev_weight(1.);
    
    bool _diel =   isData ?  (k_channel.Contains("DoubleEG")) : true ;
    bool _singleel =   isData ?  (k_channel.Contains("SingleElectron")) : true ;
    
    if(_diel){
      /// e+e-
      //      TString analysis_trigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
      TString analysis_trigger="HLT_Ele27_WPTight_Gsf_v";
      /// Trigger List (unprescaled)
      
      
      std::vector<TString> triggerslist;
      triggerslist.push_back(analysis_trigger);
      
      if(PassTrigger(analysis_trigger)){
	counter("Trigger_diel",weight);

	trigger_sf =  mcdata_correction->TriggerScaleFactor(electronTightColl, muonVetoColl, analysis_trigger);

	puweight=mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());

	if(electronTightColl.size() == 2){
	  counter("pog_diel",weight);
	  
	  if(!(muonVetoColl.size() > 0 ||  electronVetoColl.size() > 2)){
	    counter("vetolepton_diel",weight);
	    
	    if(!isData){
	      id_iso_sf=   mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronTightColl,0);
	      trigger_ps= WeightByTrigger(analysis_trigger, TargetLumi)  ;
	      reco_weight = mcdata_correction->ElectronRecoScaleFactor(electronTightColl);
	      ev_weight = weight * trigger_sf * id_iso_sf * reco_weight * puweight*trigger_ps;
	    }
	    
	    if(OppositeCharge(electronTightColl)){
	      counter("os_diel",weight);
	      
	      if(electronTightColl.at(0).Pt() > 30. && electronTightColl.at(1).Pt() > 20.){
		if(GetDiLepMass(electronTightColl) < 120. && GetDiLepMass(electronTightColl)  > 60. ){
		  
		  FillHist("zpeak_ee_w1", GetDiLepMass(electronTightColl), weight*trigger_ps, 0., 200.,400);
		  FillHist("zpeak_ee_w2", GetDiLepMass(electronTightColl), weight*trigger_sf*trigger_ps, 0., 200.,400);
		  FillHist("zpeak_ee_w3", GetDiLepMass(electronTightColl),    weight*id_iso_sf*trigger_sf*trigger_ps, 0., 200.,400);
		  FillHist("zpeak_ee_w4", GetDiLepMass(electronTightColl),    weight*id_iso_sf*trigger_sf*puweight*trigger_ps, 0., 200.,400);
		  FillHist("zpeak_ee_w5", GetDiLepMass(electronTightColl),    weight*id_iso_sf*trigger_sf*puweight*trigger_ps*reco_weight, 0., 200.,400);
		  
		  
		}
		FillCLHist(sighist_ee, "OSDiElectron", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);
		if(jets.size() > 1)                 FillCLHist(sighist_ee, "OSDiElectronDiJet", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);

	      } // pt cuts
	    }// e+e-
	    else{
	      counter("ss_diel",weight);
              if(electronTightColl.at(0).Pt() > 30. && electronTightColl.at(1).Pt() > 15.){
	      
		FillCLHist(sighist_ee, "SSDiElectron", eventbase->GetEvent(), muonVetoColl,electronTightColl,jets, ev_weight);
	      }
	    }
	  } // veto leptons 
	}// 2 el
      } /// DiEl Trigger
      return;
    }
    
    if(_singleel){
      /// e+e-                                                                                                                                                                                                                                                               
      TString analysis_trigger="HLT_Ele27_WPTight_Gsf_v";
      /// Trigger List (unprescaled)                                                                                                                                                                                                                                         
      

      std::vector<TString> triggerslist;
      triggerslist.push_back(analysis_trigger);

      if(PassTrigger(analysis_trigger)){
        counter("Trigger_singleel",weight);
	puweight=mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());
	
        if(electronTightColl.size() == 1){
          counter("pog_singleel",weight);
	  
          if(!(muonVetoColl.size() > 0 && electronVetoColl.size() > 1)){
            counter("vetolepton_singleel",weight);
	    
            if(!isData){
              id_iso_sf=   mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronTightColl,0);
              trigger_ps= WeightByTrigger(analysis_trigger, TargetLumi)  ;
              reco_weight = mcdata_correction->ElectronRecoScaleFactor(electronTightColl);
              ev_weight = weight * trigger_sf * id_iso_sf * reco_weight * puweight*trigger_ps;
	      
            }
	    
	    if(electronTightColl.at(0).Pt() > 30.){
	      
	      float METdphi = TVector2::Phi_mpi_pi(electronTightColl.at(0).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet));
	      float MT=(2.* electronTightColl.at(0).Et()*eventbase->GetEvent().MET(snu::KEvent::pfmet) * (1 - cos( METdphi)));
	      
	      if(MT > 40. && eventbase->GetEvent().MET(snu::KEvent::pfmet) > 30.) FillCLHist(sighist_ee, "SingleElectron_Wregion", eventbase->GetEvent(),  muonVetoColl,electronTightColl,jets, ev_weight);
	      
	    } // pt cuts                                                                                                                                                                                                                                                   
	  } // veto leptons                                                                                                                                                                                                                                                  
	}// 2 el                                                                                                                                                                                                                                                             
      } /// El Trigger                                                                                                                                                                                                                                                     
    }
  
  
    
      
    return;
  }
  
  
  /// make plots of POG + AN ID Efficiency
  
  GetSSSignalEfficiency(weight);
  GetOSSignalEfficiency(weight);
    /// Validation of signal MC
  if(IsSignal()) {  
    SignalValidation();
  }
  
  /// Get Efficiency of signal + Bkg using multiple triggers
  GetTriggEfficiency();



  if ((electronVetoColl.size() + muonVetoColl.size()) >2) return;
  
  CheckJetsCloseToLeptons(GetElectrons("ELECTRON_HN_VETO"), GetJets("JET_NOLEPTONVETO"));
  std::vector<snu::KFatJet> fatjetcoll = GetFatJets("FATJET_HN");

  FillHist("NFatJets", fatjetcoll.size(), weight, 0., 10., 10);
  FillHist("2DJets", fatjetcoll.size(), jets.size(),  weight, 0., 10., 10,  0., 10., 10);
  for(unsigned int ifjet=0; ifjet < fatjetcoll.size(); ifjet++){
    FillHist(("tau21"), fatjetcoll[ifjet].Tau2()/fatjetcoll[ifjet].Tau1(), weight, 0., 1., 100);
    FillHist(("PrunedMass"), fatjetcoll[ifjet].PrunedMass(),  weight, 0., 200., 100);
    FillHist(("SoftDropMass"), fatjetcoll[ifjet].SoftDropMass(),  weight, 0., 200., 100);
  }
  
  OptimiseID(true);

  //RunAnalysis("hn", "ELECTRON_HN_TIGHT","ELECTRON_HN_VETO","ELECTRON_HN_FAKELOOSE_NOD0");
  //RunAnalysis("pog_medium","ELECTRON_POG_MEDIUM","ELECTRON_POG_VETO","ELECTRON_POG_FAKELOOSE");
  //RunAnalysis("pog_tight", "ELECTRON_POG_TIGHT","ELECTRON_POG_VETO","ELECTRON_POG_FAKELOOSE");

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
  
  std::vector<snu::KElectron> electronLooseColl_tight = GetElectrons("ELECTRON16_POG_FAKELOOSE_CC");
  std::vector<snu::KElectron> electronLooseColl_medium = GetElectrons("ELECTRON16_POG_MEDIUM_FAKELOOSE_CC");
  
  std::vector<snu::KElectron> electronMedium_chargeconst  = GetElectrons("ELECTRON16_FR_POG_MEDIUM_CC");
  std::vector<snu::KElectron> electronTight_chargeconst  = GetElectrons("ELECTRON16_FR_POG_TIGHT_CC");
  std::vector<snu::KJet> jets = GetJets("JET_HN");

  if(!k_running_nonprompt){
    
    if(CheckSignalRegion(isss,electronLooseColl_tight, jets,"Signal_Mediumlooseiso_d0", weight)) FillHist("IDREF",0.  , weight, 0.,1.,1);
    if(CheckSignalRegion(isss,electronMedium_chargeconst, jets,"Signal_Mediumlooseiso_d0", weight)) FillHist("IDcutflow",0.  , weight, 0.,2.,2);
    if(CheckSignalRegion(isss,electronTight_chargeconst, jets,"Signal_Tightlooseiso_d0", weight))   FillHist("IDcutflow",1.  , weight, 0.,2.,2);
  }
  else{
    float ee_weight_medium = m_datadriven_bkg->Get_DataDrivenWeight_EE(false, electronLooseColl_medium,"ELECTRON16_POG_MEDIUM_FAKELOOSE_CC","ELECTRON16_FR_POG_MEDIUM_CC","dijet_ajet40");
    float ee_weight_tight =   m_datadriven_bkg->Get_DataDrivenWeight_EE(false, electronLooseColl_tight,"ELECTRON16_POG_FAKELOOSE_CC","ELECTRON16_POG_FAKELOOSE_CC","dijet_ajet40");
    if(CheckSignalRegion(isss, electronLooseColl_medium, jets,"Signal_Mediumlooseiso_d0", weight)) FillHist("IDcutflow",0.  , weight*ee_weight_medium, 0.,2.,2);
    if(CheckSignalRegion(isss,electronLooseColl_tight, jets,"Signal_Tightlooseiso_d0", weight))   FillHist("IDcutflow",1.  , weight*ee_weight_tight, 0.,2.,2);
  }
  
  if(CheckSignalRegion(isss, GetElectrons("HNTight_isoref"),jets,"", weight)) FillHist("ISOREF", 0., weight, 0.,1.,1);
  
  if(!k_running_nonprompt){
    std::vector <TString> isocuts;
    isocuts.push_back("045");
    isocuts.push_back("050");
    isocuts.push_back("055");
    isocuts.push_back("060");
    isocuts.push_back("065");
    isocuts.push_back("075");
    isocuts.push_back("100");
    isocuts.push_back("125");
    int iisocut(0);
    for(unsigned int iiso = 0 ; iiso < isocuts.size(); iiso++){
      for(unsigned int iiso2 = 0 ; iiso2 < isocuts.size(); iiso2++, iisocut++){
	if(CheckSignalRegion(isss, GetElectrons(("HNTight_b"+isocuts[iiso]+"_e"+isocuts[iiso2]).Data()),jets,"", weight)) FillHist("ISOcutflow", iisocut , weight,0., 64., 64);
	
      }
    }
  }
  else{

    std::vector <TString> isocuts;
    isocuts.push_back("045");
    isocuts.push_back("050");
    isocuts.push_back("055");
    isocuts.push_back("060");
    isocuts.push_back("065");
    isocuts.push_back("075");
    isocuts.push_back("100");
    isocuts.push_back("125");

    int iisocut(0);
    for(unsigned int iiso = 0 ; iiso < isocuts.size(); iiso++){
      for(unsigned int iiso2 = 0 ; iiso2 < isocuts.size(); iiso2++, iisocut++){
	std::vector<snu::KElectron>  eltight = GetElectrons(("HNTight_b"+isocuts[iiso]+"_e"+isocuts[iiso2]).Data());
	std::vector<snu::KElectron>  elloose=  GetElectrons("ELECTRON16_HN_FAKELOOSE_NOD0");
        float evweight = m_datadriven_bkg->Get_DataDrivenWeight_EE(false,  elloose, "ELECTRON16_HN_FAKELOOSE_NOD0",("HNTight_b"+isocuts[iiso]+"_e"+isocuts[iiso2]).Data(), "opt_dijet_ajet40");
	if(CheckSignalRegion(isss, eltight,jets,"", weight)) FillHist("ISOcutflow", iisocut , evweight,0., 64., 64);
      }
    }
  }
  

  if(CheckSignalRegion(isss,GetElectrons("HNTight_dxy_ref"),jets,"", weight)) FillHist("d0cutflow", 0., weight, 0.,20.,20);
  
  if(!k_running_nonprompt){
    
    std::vector<TString> dxyb;
    dxyb.push_back("050");
    dxyb.push_back("025");
    dxyb.push_back("015");
    dxyb.push_back("010");
    std::vector<TString> dxye;
    dxye.push_back("100");
    dxye.push_back("050");
    dxye.push_back("040");
    dxye.push_back("025");
    dxye.push_back("020");

    int idxy(1);
    for(unsigned int ib=0; ib <dxyb.size(); ib++){
      for(unsigned int ie=0; ie <dxye.size(); ie++, idxy++){
	if(CheckSignalRegion(isss, GetElectrons(("HNTight_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data()) ,jets,"", weight)) FillHist("d0cutflow", idxy, weight,  0.,20.,20);
	
      }
    }
  }
  else{
    std::vector<TString> dxyb;
    dxyb.push_back("050");
    dxyb.push_back("025");
    dxyb.push_back("015");
    dxyb.push_back("010");
    std::vector<TString> dxye;
    dxye.push_back("100");
    dxye.push_back("050");
    dxye.push_back("040");
    dxye.push_back("025");
    dxye.push_back("020");

    int idxy(1);
    for(unsigned int ib=0; ib <dxyb.size(); ib++){
      for(unsigned int ie=0; ie <dxye.size(); ie++, idxy++){
	std::vector<snu::KElectron>  eltight = GetElectrons(("HNTight_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data());
	std::vector<snu::KElectron>  elloose  =GetElectrons("ELECTRON16_HN_FAKELOOSE_NOD0");
	float evweight = m_datadriven_bkg->Get_DataDrivenWeight_EE(false,  elloose, "ELECTRON16_HN_FAKELOOSE_NOD0",("HNTight_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data(), "opt_dijet_ajet40");
	FillHist("d0cutflow", idxy, evweight,  0.,20.,20);
      }
    }
  }
  
  
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
	  
	  FillCLHist(trilephist, "TriMuon", eventbase->GetEvent(), muonColl, GetElectrons("ELECTRON_HN_VETO"),GetJets("JET_HN"), weight);

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

  std::vector<snu::KElectron> electronColl             = GetElectrons(true, false,elid);
  
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
		  for(unsigned int ij1=0; ij1 < jetColl_hn.size(); ij1++){
		    cout << jetColl_hn.at(ij1).Eta() << endl;
		  }
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



