// $Id: ExampleAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQGetSignalEff Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *d
 ***************************************************************************/

/// Local includes
#include "GetSignalEff.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (GetSignalEff);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
GetSignalEff::GetSignalEff() :  AnalyzerCore() {


  // To have the correct name in the log:                                                                                                                            
  SetLogName("GetSignalEff");

  Message("In GetSignalEff constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void GetSignalEff::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //

   Message("Making clever hists for Z ->ll test code", INFO);

   //// Initialise Plotting class functions
   /// MakeCleverHistograms ( type, "label")  type can be muhist/elhist/jethist/sighist


   return;
}


void GetSignalEff::ExecuteEvents()throw( LQError ){
  
  if(!isData){
    std::vector<snu::KTruth> truthColl= eventbase->GetTruth();
    //TruthPrintOut();
    
    for(unsigned int i = 0; i < truthColl.size(); i++){
    if(fabs(eventbase->GetTruth().at(i).PdgId()) == 11) {
      
      if(eventbase->GetTruth().at(i).GenStatus() ==23){
	if(fabs(eventbase->GetTruth().at(eventbase->GetTruth().at(i).IndexMother()).PdgId()) != 9900012){
	  FillHist("Electron_FromW_LowMass_Pt", eventbase->GetTruth().at(i).Pt(), 1., 0., 500., 250);
	  FillHist("Electron_FromW_LowMass_Pt_lowpt", eventbase->GetTruth().at(i).Pt(), 1., 0., 100., 100);
	  FillHist("Electron_FromW_LowMass_Eta", eventbase->GetTruth().at(i).Eta(), 1., -3., 3., 60);
	}	  
      }
      
      if(fabs(eventbase->GetTruth().at(eventbase->GetTruth().at(i).IndexMother()).PdgId()) ==24) {
	FillHist("Electron_FromW_LowMass_Pt", eventbase->GetTruth().at(i).Pt(), 1., 0., 500., 250);
	FillHist("Electron_FromW_LowMass_Pt_lowpt", eventbase->GetTruth().at(i).Pt(), 1., 0., 100., 100);
	FillHist("Electron_FromW_LowMass_Eta", eventbase->GetTruth().at(i).Eta(), 1., -3., 3., 60);

      }
      if(fabs(eventbase->GetTruth().at(eventbase->GetTruth().at(i).IndexMother()).PdgId()) == 9900012){
	FillHist("Electron_FromN_LowMass_Pt", eventbase->GetTruth().at(i).Pt(), 1., 0., 500., 250);
        FillHist("Electron_FromN_LowMass_Pt_lowpt", eventbase->GetTruth().at(i).Pt(), 1., 0., 100., 100);
	FillHist("Electron_FromN_LowMass_Eta", eventbase->GetTruth().at(i).Eta(), 1., -3., 3., 60);
      }
      
      if(fabs(eventbase->GetTruth().at(eventbase->GetTruth().at(i).IndexMother()).PdgId()) == 9900012){
        FillHist("Electron_FromN_LowMass_Pt", eventbase->GetTruth().at(i).Pt(), 1., 0., 500., 250);
        FillHist("Electron_FromN_LowMass_Pt_lowpt", eventbase->GetTruth().at(i).Pt(), 1., 0., 100., 100);
        FillHist("Electron_FromN_LowMass_Eta", eventbase->GetTruth().at(i).Eta(), 1., -3., 3., 60);
      }

      
      if(fabs(eventbase->GetTruth().at(eventbase->GetTruth().at(i).IndexMother()).PdgId()) ==24) {
        FillHist("Electron_FromW_LowMass_Pt", eventbase->GetTruth().at(i).Pt(), 1., 0., 500., 250);
        FillHist("Electron_FromW_LowMass_Pt_lowpt", eventbase->GetTruth().at(i).Pt(), 1., 0., 100., 100);
        FillHist("Electron_FromW_LowMass_Eta", eventbase->GetTruth().at(i).Eta(), 1., -3., 3., 60);

      }
      
    }
    if(fabs(eventbase->GetTruth().at(i).PdgId()) < 6) {
      
      int mother_i = eventbase->GetTruth().at(i).IndexMother();
      
      while(fabs(eventbase->GetTruth().at(mother_i).PdgId()) ==24){
	mother_i = eventbase->GetTruth().at(mother_i).IndexMother();
      }
      if(fabs(eventbase->GetTruth().at(mother_i).PdgId()) == 9900012){
	FillHist("Jets_FromW_LowMass_Pt", eventbase->GetTruth().at(i).Pt(), 1., 0., 500., 250);
	FillHist("Jets_FromW_LowMass_Eta", eventbase->GetTruth().at(i).Eta(), 1., -5., 5., 100);
      }
    }
  }
  }
  vector<TString> muonIDs;
  muonIDs.push_back("MUON_HN_EFF_PT");
  muonIDs.push_back("MUON_HN_EFF_POG");
  muonIDs.push_back("MUON_HN_EFF_DXY");
  muonIDs.push_back("MUON_HN_EFF_DZ");
  muonIDs.push_back("MUON_HN_EFF_dxysig");
  muonIDs.push_back("MUON_HN_EFF_iso");
  muonIDs.push_back("MUON_HN_TIGHT");

  
  if(!isData)weight*= MCweight;
  weight*=WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi);

  FillEventCutFlow(0, "NoCut", weight);


  std::vector<snu::KFatJet> fatjets = GetFatJets("FATJET_HN");
  std::vector<snu::KJet> jets = GetJets("JET_HN");
  std::vector<snu::KJet> hnjets = GetJetsWFT("JET_HN","FATJET_HN");

  FillHist("hnjets", hnjets.size(), weight, 0., 10., 10);
  FillHist("NFatJets", fatjets.size(), weight, 0., 10., 10);
  FillHist("2DJets", fatjets.size(), jets.size(),  weight, 0., 10., 10,  0., 10., 10);
  for(unsigned int ifjet=0; ifjet < fatjets.size(); ifjet++){
    FillHist(("tau21"), fatjets[ifjet].Tau2()/fatjets[ifjet].Tau1(), weight, 0., 1., 100);
    FillHist(("PrunedMass"), fatjets[ifjet].PrunedMass(),  weight, 0., 200., 100);
    FillHist(("SoftDropMass"), fatjets[ifjet].SoftDropMass(),  weight, 0., 200., 100);
  }

  
  if(PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"))   {
    FillEventCutFlow(0, "HLT_Mu17_",weight);
    TString analysis_trigger_muon="HLT_IsoMu24_v";
    TString analysis_trigger_tkmuon="HLT_IsoTkMu24_v";
    if((PassTrigger(analysis_trigger_muon) || PassTrigger(analysis_trigger_tkmuon))){
      FillEventCutFlow(0, "HLT_IsoMu24",weight);
    }
  }
  else {
    TString analysis_trigger_muon="HLT_IsoMu24_v";
    TString analysis_trigger_tkmuon="HLT_IsoTkMu24_v";
    if((PassTrigger(analysis_trigger_muon) || PassTrigger(analysis_trigger_tkmuon))){
      FillEventCutFlow(0, "HLT_IsoMu24",weight);
    }
    else return;
  }
  

  std::vector<snu::KJet> alljets =  GetJets("JET_NOLEPTONVETO");

  

  std::vector<snu::KMuon> muons_veto = GetMuons("MUON_HN_VETO",false);
  std::vector<snu::KElectron> electrons_veto = GetElectrons(true,true,"ELECTRON_HN_VETO");

  if(muons_veto.size() == 2 && electrons_veto.size()==0){
    for(unsigned int iid=0; iid < muonIDs.size(); iid++){

      TString muid = muonIDs[iid];
      std::vector<snu::KMuon> muons_eff=GetMuons(muid);

      if(SameCharge(muons_eff)){

	
	if(PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")) {
	  if(muons_eff[0].Pt() < 20 || muons_eff[1].Pt() < 10) continue;
	}
	else{
	  TString analysis_trigger_muon="HLT_IsoMu24_v";
	  TString analysis_trigger_tkmuon="HLT_IsoTkMu24_v";
	  if((PassTrigger(analysis_trigger_muon) || PassTrigger(analysis_trigger_tkmuon))){
	    if(muons_eff[0].Pt() < 25 || muons_eff[1].Pt() < 5) continue;
          }
	  
	}
	FillEventCutFlow(0, muid, weight);
        
	if(muonIDs[iid] == "MUON_HN_TIGHT"){
	  
	  if(jets.size() > 1) FillEventCutFlow(0, "DiJet", weight);
	  if((jets.size() + fatjets.size()) > 1 )  FillEventCutFlow(0, "NewDiJet", weight);
	  
	  
          if(CheckSignalRegion(true,muons_eff, electrons_veto,jets, alljets,"Low", weight)) FillEventCutFlow(0, "LowMass", weight);
          if(CheckSignalRegion(true,muons_eff,electrons_veto, jets, alljets,"", weight)) FillEventCutFlow(0, "MidMass", weight);
          if(CheckSignalRegion(true,muons_eff,electrons_veto, jets, alljets,"High", weight)) FillEventCutFlow(0, "HighMass", weight);
	}
	
      }
    }
    
  }
  
  TString pogt="MUON_POG_TIGHT";
  TString pogm="MUON_POG_MEDIUM";
  TString gentid="MUON_HNGENT_TIGHT";
  TString hnid = "MUON_HN_TIGHT";
  TString hnid_iso = "MUON_HN_TIGHT_iso";
  TString hnid_dxy = "MUON_HN_TIGHT_dxy";
  TString hnidm = "MUON_HN_MEDIUM";

  if(k_running_nonprompt){
    pogt="MUON_POG_FAKETIGHT";
    pogm="MUON_POG_FAKEMEDIUM";
    gentid="MUON_HNGENT_LOOSE";
    hnid="MUON_HN_LOOSE";
    hnid_iso="MUON_HN_LOOSE";
    hnid_dxy="MUON_HN_LOOSE_dxy";
    hnidm="MUON_HN_LOOSEMEDIUM";
  }
  bool cb_1=false;
  bool cb_2=false;
 
   
  std::vector<snu::KMuon> pogmedium  = GetMuons(pogm,false);
  std::vector<snu::KMuon> pogtight  = GetMuons(pogt,false);
  std::vector<snu::KMuon> gent= GetMuons(gentid,false);
  std::vector<snu::KMuon>  hn=  GetMuons(hnid,false);
  std::vector<snu::KMuon>  hn_med=  GetMuons(hnidm,false);
  std::vector<snu::KMuon>  hn_iso=  GetMuons(hnid_iso,false);
  std::vector<snu::KMuon>  hn_dxy=  GetMuons(hnid_dxy,false);

  
  float weight_pm=weight;
  float weight_pt=weight;
  float weight_g=weight;
  float weight_hn=weight;
  float weight_hnm=weight;
  float weight_hniso=weight;
  float weight_hndxy=weight;


  if(k_running_nonprompt){

    if(SameCharge(pogmedium))weight_pm=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, pogmedium, PassID(pogmedium[0],"MUON_POG_MEDIUM"),  PassID(pogmedium[1],"MUON_POG_MEDIUM"), "pogmedium","pogmedium", cb_1, cb_2,"ptcorr_eta",0.25,0.25,false,true);

    
    if(SameCharge(pogtight))weight_pt =m_datadriven_bkg->Get_DataDrivenWeight_MM(false, pogtight, PassID(pogtight[0],"MUON_POG_TIGHT"),  PassID(pogtight[1],"MUON_POG_TIGHT"), "pogtight","pogtight",cb_1, cb_2, "ptcorr_eta",0.15,0.15,false,true);

    
    if(SameCharge(gent))weight_g=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, gent, PassID(gent[0],"MUON_HNGENT_TIGHT"),  PassID(gent[1],"MUON_HNGENT_TIGHT"),"gent","gent",cb_1, cb_2, "ptcorr_eta",0.1,0.1,false,true);
     
    if(SameCharge(hn))weight_hn=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, hn, PassID(hn[0],"MUON_HN_TIGHT"),  PassID(hn[1],"MUON_HN_TIGHT"), "Tight0.07_0.005_3_0.04","Tight0.07_0.005_3_0.04",cb_1, cb_2, "ptcorr_eta",0.07,0.07,false,true);

    if(SameCharge(hn_med))weight_hnm=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, hn_med, PassID(hn_med[0],"MUON_HN_MEDIUM"),  PassID(hn_med[1],"MUON_HN_MEDIUM"), "Medium0.07_0.005_3_0.04","Medium0.07_0.005_3_0.04",cb_1, cb_2, "ptcorr_eta",0.07,0.07,false,true);

    if(SameCharge(hn_iso))weight_hniso=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, hn_iso, PassID(hn_iso[0],"MUON_HN_TIGHT_iso"),  PassID(hn_iso[1],"MUON_HN_TIGHT_iso"), "Tight0.09_0.005_3_0.04","Tight0.09_0.005_3_0.04",cb_1, cb_2, "ptcorr_eta",0.09,0.09,false,true);

    if(SameCharge(hn_dxy))weight_hndxy=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, hn_dxy, PassID(hn_dxy[0],"MUON_HN_TIGHT_dxy"),  PassID(hn_dxy[1],"MUON_HN_TIGHT_dxy"), "Tight0.07_0.05_3_0.04","Tight0.07_0.05_3_0.04",cb_1, cb_2, "ptcorr_eta",0.07,0.07,false,true);
     
  }
  
  
  if(CheckSignalRegion(true,pogmedium,electrons_veto, jets, alljets,"Low", weight)) FillEventCutFlow(2, "LowMass_POGMEDIUM", weight_pm);
  if(CheckSignalRegion(true,pogtight,electrons_veto, jets, alljets,"Low", weight)) FillEventCutFlow(2, "LowMass_POGTIGHT", weight_pt);
  if(CheckSignalRegion(true,gent,electrons_veto, jets, alljets,"Low", weight)) FillEventCutFlow(2, "LowMass_GENT", weight_g);
  if(CheckSignalRegion(true,hn,electrons_veto, jets, alljets,"Low", weight)) FillEventCutFlow(2, "LowMass_HN", weight_hn);
  if(CheckSignalRegion(true,hn_med,electrons_veto, jets, alljets,"Low", weight)) FillEventCutFlow(2, "LowMass_HN_MEDIUM", weight_hnm);
  if(CheckSignalRegion(true,hn_iso,electrons_veto, jets, alljets,"Low", weight)) FillEventCutFlow(2, "LowMass_HN_ISO", weight_hniso);
  if(CheckSignalRegion(true,hn_dxy,electrons_veto, jets, alljets,"Low", weight)) FillEventCutFlow(2, "LowMass_HN_DXY", weight_hndxy);
  
  if(CheckSignalRegion(true,pogmedium,electrons_veto, jets, alljets,"High", weight)) FillEventCutFlow(4, "HighMass_POGMEDIUM", weight_pm);
  if(CheckSignalRegion(true,pogtight,electrons_veto, jets, alljets,"High", weight)) FillEventCutFlow(4, "HighMass_POGTIGHT", weight_pt);
  if(CheckSignalRegion(true,gent,electrons_veto, jets, alljets,"High", weight)) FillEventCutFlow(4, "HighMass_GENT", weight_g);
  if(CheckSignalRegion(true,hn,electrons_veto, jets, alljets,"High", weight)) FillEventCutFlow(4, "HighMass_HN", weight_hn);
  if(CheckSignalRegion(true,hn_med,electrons_veto, jets, alljets,"High", weight)) FillEventCutFlow(4, "HighMass_HN_MEDIUM", weight_hnm);
  if(CheckSignalRegion(true,hn_iso,electrons_veto, jets, alljets,"High", weight)) FillEventCutFlow(4, "HighMass_HN_ISO", weight_hniso);
  if(CheckSignalRegion(true,hn_dxy,electrons_veto, jets, alljets,"High", weight)) FillEventCutFlow(4, "HighMass_HN_DXY", weight_hndxy);
  
  


  return;
}// End of execute event loop
  



void GetSignalEff::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger << DEBUG << "END OF CYCLE: isdata=" << isData <<  LQLogger::endmsg;

}


void GetSignalEff::BeginCycle() throw( LQError ){
  
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




bool GetSignalEff::CheckSignalRegion( bool isss,  std::vector<snu::KMuon> muons, std::vector<snu::KElectron> el,  std::vector<snu::KJet> jets,  std::vector<snu::KJet> alljets, TString name, float w){
				    
  
  bool debug=false;


  if(el.size() > 0) return false ;
  if(muons.size() != 2 ) {if(debug)cout << "Fail el size" << endl; return false ;}
                                                                                                                                                                                 
  if(muons.at(0).Pt() < 20.) {if(debug)cout << "Fail pt1 " << endl; return false;}
  if(muons.at(1).Pt() < 10.)  {if(debug)cout << "Fail pt2  " << endl; return false;}

  if(isss&&!SameCharge(muons)) {if(debug)cout << "Fail ss " << endl; return false;}

  if(!isss&&SameCharge(muons)) {if(debug)cout << "Fail os  " << endl; return false;}
  if(jets.size() < 2) {if(debug)cout << "Fail jets " << endl; return false;}

  snu::KParticle mm = muons.at(0) + muons.at(1);
  if(mm.M()  < 10.) {if(debug)cout << "Fail mee  " << endl; return false;}

  if(name.Contains("Low")){
    if(mm.M()  > 70.) {if(debug)cout << "Fail mee  " << endl; return false;}
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
    if((jets[0]+jets[1]).M() > 200.) return false;
    if((jets[0] + muons[0] + muons[1]).M() > 250.)  return false;
  }
  else{
    if(jj.M() > 120.) {if(debug)cout << "Fail mjj  " << endl; return false;}
    if(jj.M() < 50.) {if(debug)cout << "Fail mjj  " << endl; return false;}

  }


  float ST = muons[0].Pt() + muons[1].Pt();
  float looseST = muons[0].Pt() + muons[1].Pt();

  std::vector<snu::KJet> loosejets=GetJets("JET_NOCUT");


  float HT=0.;
  for(unsigned int ij=0; ij <alljets.size(); ij++){
    ST+= alljets[ij].Pt();

  }
  for(unsigned int ij=0; ij <jets.size(); ij++){

    HT+= jets[ij].Pt();

  }

  float LT = muons[0].Pt() + muons[1].Pt();
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



  float dphi_1 = fabs(TVector2::Phi_mpi_pi(muons.at(0).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet)));
  float MT_1 = sqrt(2.* muons.at(0).Et()*eventbase->GetEvent().PFMET() * (1 - cos( dphi_1)));

  float dphi_2 = fabs(TVector2::Phi_mpi_pi(muons.at(0).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet)));
  float MT_2 = sqrt(2.* muons.at(1).Et()*eventbase->GetEvent().PFMET() * (1 - cos( dphi_2)));

  if(name.Contains("Low")){
    if(MT_1 > 100.) return false;
    if(MT_2 > 100.) return false;
    if(ST > 450.)  return false;
    snu::KParticle lljj = jets.at(m) + jets.at(n) + muons.at(0) + muons.at(1);

  }
  else   if(name.Contains("High")){
    snu::KParticle lljj = jets.at(m) + jets.at(n) + muons.at(0) + muons.at(1);

    if(ST < 400.)  return false;
    if(lljj.M() < 200.) return false;
    if(muons.at(0).Pt() < 50.) return false;
    if(muons.at(1).Pt() < 25.) return false;
  }
  else{


    if(LT/HT < 0.2)  return false;
    if(LT < 35)  return false;
    if(HT < 35)  return false;

    float dRmj=-999.;
    float ratiomj=0.;
    for(unsigned int im=0; im <muons.size(); im++){
      for(unsigned int ij2=0; ij2 <jets.size(); ij2++){
        if(muons.at(im).DeltaR(jets.at(ij2)) >  0.5){
          if(muons.at(im).DeltaR(jets.at(ij2)) > dRmj){
            dRmj=muons.at(im).DeltaR(jets.at(ij2)) ;
            if(im==0) ratiomj = jets.at(ij2).Pt() / muons[im].Pt();
          }
        }
      }
    }
    if(ratiomj > 3.) return false;
    if(dRmj > 4.5) return false;

    if(muons[0].DeltaR(muons[1]) > 4.) return false;

  }
  if(name.Contains("High")){
    if((muons.at(0).Pt() + muons.at(1).Pt() ) < 50.) return false;
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
}


GetSignalEff::~GetSignalEff() {
  
  Message("In GetSignalEff Destructor" , INFO);
  
}


void GetSignalEff::FillEventCutFlow(int cf,TString cut,  float weight){

  if(cf==0){
    if(GetHist( "mm_eventcutflow_DoubleMuon")) {
      GetHist( "mm_eventcutflow_DoubleMuon")->Fill(cut,weight);

    }
    else{

      vector<TString> muonIDs;
      muonIDs.push_back("NoCut");
      muonIDs.push_back("HLT_Mu17_");      muonIDs.push_back("HLT_IsoMu24");
      muonIDs.push_back("MUON_HN_EFF_PT");
      muonIDs.push_back("MUON_HN_EFF_POG");
      muonIDs.push_back("MUON_HN_EFF_DXY");
      muonIDs.push_back("MUON_HN_EFF_DZ");
      muonIDs.push_back("MUON_HN_EFF_dxysig");
      muonIDs.push_back("MUON_HN_EFF_iso");
      muonIDs.push_back("MUON_HN_TIGHT");

      muonIDs.push_back("DiJet");
      muonIDs.push_back("NewDiJet");
      muonIDs.push_back("LowMass");
      muonIDs.push_back("MidMass");
      muonIDs.push_back("HighMass");
      AnalyzerCore::MakeHistograms( "mm_eventcutflow_DoubleMuon",muonIDs.size(),0.,float(muonIDs.size()));

      for(unsigned int iid=0; iid < muonIDs.size(); iid++){

        TString elid = muonIDs[iid];
        GetHist( "mm_eventcutflow_DoubleMuon")->GetXaxis()->SetBinLabel(1+iid,elid);
      }
      GetHist( "mm_eventcutflow_DoubleMuon")->Fill(cut,weight);
      
    }
  }
  
  if(cf==2){
    if(GetHist( "mm_eventcutflow_LowMass_DoubleMuon")) {
      GetHist( "mm_eventcutflow_LowMass_DoubleMuon")->Fill(cut,weight);
    }
  else{

    vector<TString> IDs;
    IDs.push_back("LowMass_POGMEDIUM");
    IDs.push_back("LowMass_POGTIGHT");
    IDs.push_back("LowMass_GENT");
    IDs.push_back("LowMass_HN");
    IDs.push_back("LowMass_HN_MEDIUM");
    IDs.push_back("LowMass_HN_ISO");
    IDs.push_back("LowMass_HN_DXY");



    AnalyzerCore::MakeHistograms( "mm_eventcutflow_LowMass_DoubleMuon",IDs.size(),0.,float(IDs.size()));

    for(unsigned int iid=0; iid < IDs.size(); iid++){

      TString elid = IDs[iid];
      GetHist( "mm_eventcutflow_LowMass_DoubleMuon")->GetXaxis()->SetBinLabel(1+iid,elid);
    }
    GetHist( "mm_eventcutflow_LowMass_DoubleMuon")->Fill(cut,weight);

  }
 }

if(cf==4){
  if(GetHist( "mm_eventcutflow_HighMass_DoubleMuon")) {
    GetHist( "mm_eventcutflow_HighMass_DoubleMuon")->Fill(cut,weight);

  }
  else{

    vector<TString> IDs;
    IDs.push_back("HighMass_POGMEDIUM");
    IDs.push_back("HighMass_POGTIGHT");
    IDs.push_back("HighMass_GENT");
    IDs.push_back("HighMass_HN");
    IDs.push_back("HighMass_HN_MEDIUM");
    IDs.push_back("HighMass_HN_ISO");
    IDs.push_back("HighMass_HN_DXY");


    AnalyzerCore::MakeHistograms( "mm_eventcutflow_HighMass_DoubleMuon",IDs.size(),0.,float(IDs.size()));

    for(unsigned int iid=0; iid < IDs.size(); iid++){

      TString elid = IDs[iid];
      GetHist( "mm_eventcutflow_HighMass_DoubleMuon")->GetXaxis()->SetBinLabel(1+iid,elid);
    }
    GetHist( "mm_eventcutflow_HighMass_DoubleMuon")->Fill(cut,weight);

  }
 }

}





void GetSignalEff::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void GetSignalEff::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this GetSignalEffCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void GetSignalEff::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
}



