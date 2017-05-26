/***************************************************************************
 * @Project: LQFakeRateMC Frame - ROOT-based analysis framework for Korea SNU
OB * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/


/// Local includes
#include "FakeRateMC.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (FakeRateMC);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
FakeRateMC::FakeRateMC() :  AnalyzerCore(),  out_electrons(0) {


  // To have the correct name in the log:                                                                                                                            
  SetLogName("FakeRateMC");

  Message("In FakeRateMC constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


TDirectory* FakeRateMC::getTemporaryDirectory(void) const
{

  
  gROOT->cd();
  TDirectory* tempDir = 0;
  int counter = 0;
  while (not tempDir) {
    // First, let's find a directory name that doesn't exist yet:               
    std::stringstream dirname;
    dirname << "AnalyzerCore_%i" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    // Let's try to make this directory:                                        
    tempDir = gROOT->mkdir((dirname.str()).c_str());

  }

  return tempDir;
}

void FakeRateMC::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
  
  TDirectory* origDir = gDirectory;

  string analysisdir = getenv("LQANALYZER_DIR");

  TFile* file_fake_elmc  = TFile::Open( (analysisdir + "/data/Fake/"+getenv("yeartag")+"/FakeRate13TeV_El_2016_mcdilep.root").c_str());
  TFile* file_fake_mumc  = TFile::Open( (analysisdir + "/data/Fake/"+getenv("yeartag")+"/FakeRate13TeV_Mu_2016_mcdilep.root").c_str());

  TDirectory* tempDir = getTemporaryDirectory();
  tempDir->cd();

  MuonFR =  dynamic_cast<TH2D*> (( file_fake_mumc->Get("FakeRate_mu_pt"))->Clone());
  MuonFRcorr =  dynamic_cast<TH2D*> (( file_fake_mumc->Get("FakeRate_mu_ptcorr"))->Clone());
  MuonFRcbj =  dynamic_cast<TH2D*> (( file_fake_mumc->Get("FakeRate_mu_cj_pt"))->Clone());
  MuonFRcbjcorr =  dynamic_cast<TH2D*> (( file_fake_mumc->Get("FakeRate_mu_cj_ptcorr"))->Clone());
  MuonFRncbj =  dynamic_cast<TH2D*> (( file_fake_mumc->Get("FakeRate_mu_ncj_pt"))->Clone());
  MuonFRncbjcorr =  dynamic_cast<TH2D*> (( file_fake_mumc->Get("FakeRate_mu_ncj_ptcorr"))->Clone());

  ElFR =  dynamic_cast<TH2D*> (( file_fake_elmc->Get("FakeRate_el_pt"))->Clone());
  ElFRcorr =  dynamic_cast<TH2D*> (( file_fake_elmc->Get("FakeRate_el_ptcorr"))->Clone());
  ElFRcbj =  dynamic_cast<TH2D*> (( file_fake_elmc->Get("FakeRate_el_cj_pt"))->Clone());
  ElFRcbjcorr =  dynamic_cast<TH2D*> (( file_fake_elmc->Get("FakeRate_el_cj_ptcorr"))->Clone());
  ElFRncbj =  dynamic_cast<TH2D*> (( file_fake_elmc->Get("FakeRate_el_ncj_pt"))->Clone());
  ElFRncbjcorr =  dynamic_cast<TH2D*> (( file_fake_elmc->Get("FakeRate_el_ncj_ptcorr"))->Clone());


  file_fake_elmc->Close();
  delete file_fake_elmc;
  file_fake_mumc->Close();
  delete file_fake_mumc;

  origDir->cd();



  return;
}


void FakeRateMC::ExecuteEvents()throw( LQError ){
    
  //// Initial event cuts
  /// MET FIleters 
  if(!PassMETFilter()) return;     
  
  /// Require good promary vertex 
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex  
  numberVertices = eventbase->GetEvent().nVertices();   
  
  weight = weight*WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi)*MCweight;
  ExecuteEventsMuon("MUON_HN_LOOSE","MUON_HN_TIGHT","HN",0.05);
  ExecuteEventsMuon("MUON_HNGENT_LOOSE","MUON_HNGENT_TIGHT","GENT", 0.1);
  ExecuteEventsElectron("ELECTRON_HN_FAKELOOSE","ELECTRON_HN_TIGHT","HN",0.05);
  ExecuteEventsElectron("ELECTRON_MVA_FAKELOOSE","ELECTRON_MVA_TIGHT","MVA",0.05);
  ExecuteEventsElectron("ELECTRON_GENT_FAKELOOSE","ELECTRON_GENT_TIGHT","GENT",0.1);
}


void FakeRateMC::ExecuteEventsMuon(TString looseid, TString tightid, TString tag, float tightiso){
  
  std::vector<snu::KMuon> loosemuons = GetMuons(looseid,true);
  Float_t ptbins[11] = { 5., 10., 15.,20.,25.,30.,35.,45.,60.,100., 200.};
  Float_t etabins2[5] = { 0.,0.8,  1.479, 2.,  2.5};
  std::vector<snu::KJet> alljets = GetJets("JET_NOLEPTONVETO");
  std::vector<snu::KJet> jets = GetJets("JET_HN");

  if(SameCharge(loosemuons)){
    FillHist(tag+"SSmu_mass", (loosemuons[0]+loosemuons[1]).M(), weight , 0., 200., 200);
    for(int x=0; x < loosemuons.size(); x++){
      FillHist(tag+"_MuonType_L_PF",loosemuons[x].GetType(),weight, 0., 41., 41);
      FillHist(tag+"_MuonType_L_PF_mother",fabs(loosemuons[x].MotherPdgId()),weight, 0., 600., 600);
      if(PassID(loosemuons[x], tightid) )  {
	FillHist(tag+"_MuonType_L_PF_tight",loosemuons[x].GetType(),weight, 0., 41., 41);
	FillHist(tag+"_MuonType_L_PF_mother_tight",fabs(loosemuons[x].MotherPdgId()),weight, 0., 600., 600);
      }
      
      if(!loosemuons[x].MCMatched())  {
	FillHist(tag+"_MuonType_LJ_PF",loosemuons[x].GetType(),weight, 0., 41., 41);
	if(PassID(loosemuons[x], tightid) )             FillHist(tag+"_MuonType_LJ_PF_tight",loosemuons[x].GetType(),weight, 0., 41., 41);
	
	if(loosemuons[x].GetType()==2){
	  FillHist(tag+"_MuonType_LJ_PF_mother",fabs(loosemuons[x].MotherPdgId()),weight, 0., 600., 600);
	  if(PassID(loosemuons[x], tightid))  FillHist(tag+"_MuonType_LJ_PF_mother_tight",fabs(loosemuons[x].MotherPdgId()),weight, 0., 600., 600);
	}
	FillHist(tag+"_MuonType_mother_LJ_PF",loosemuons[x].GetType(),fabs(loosemuons[x].MotherPdgId()), weight, 0., 41., 41,  0., 600., 600);
	if(PassID(loosemuons[x], tightid))FillHist(tag+"_MuonType_mother_LJ_PF_tight",loosemuons[x].GetType(),fabs(loosemuons[x].MotherPdgId()), weight, 0., 41., 41,  0., 600., 600);
      }
    }
  }

  if(loosemuons.size()==3){
    FillHist(tag+"mumumu_mass", (loosemuons[0]+loosemuons[1]+loosemuons[2]).M(), weight , 0., 200.,200);

    for(int x=0; x < loosemuons.size(); x++){
      if(!loosemuons[x].MCMatched())  {
	FillHist(tag+"_MuonType_LLJ_PF",loosemuons[x].GetType(),weight, 0., 41., 41);
	if(PassID(loosemuons[x], tightid) )             FillHist(tag+"_MuonType_LLJ_PF_tight",loosemuons[x].GetType(),weight, 0., 41., 41);
	
	FillHist(tag+"_MuonType_LLJ_PF_mother",fabs(loosemuons[x].MotherPdgId()),weight, 0., 600., 600);
	if(PassID(loosemuons[x], tightid))  FillHist(tag+"_MuonType_LLJ_PF_mother_tight",fabs(loosemuons[x].MotherPdgId()),weight, 0., 600., 600);
	
	FillHist(tag+"_MuonType_mother_LLJ_PF",loosemuons[x].GetType(),fabs(loosemuons[x].MotherPdgId()), weight, 0., 41., 41,  0., 600., 600);
	if(PassID(loosemuons[1], tightid))FillHist(tag+"_MuonType_mother_LLJ_PF_tight",loosemuons[x].GetType(),fabs(loosemuons[x].MotherPdgId()), weight, 0., 41., 41,  0., 600., 600);
      }
    }
  }
  
  if( ((k_sample_name.Contains("TT"))&&SameCharge(loosemuons)) || (!k_sample_name.Contains("TT")&&loosemuons.size()==2)){
    if(PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") || PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")){
      if(loosemuons[0].Pt() > 20.){
	int nT=0;
	if(TruthMatched(loosemuons[0]))nT++;
	if(TruthMatched(loosemuons[1]))nT++;
	
	  
	//// Check both are tight
	  if(PassID(loosemuons[0], tightid) && (PassID(loosemuons[1], tightid))){
	    FillHist(tag+"_LJ_TT_dimu", loosemuons[1].Pt() , weight, 0., 100., 20);
	    FillHist(tag+"_LJ_TT_dimu_nw", loosemuons[1].Pt() , 1., 0., 100., 20);
	    float mu_pt_corr = loosemuons.at(1).Pt()*(1+max(0.,(loosemuons.at(1).RelIso04()-tightiso))) ;
	    
	    FillHist(tag+"_LJ_TT_dimu_corr", mu_pt_corr  , weight, 0., 100., 20);
	    FillHist(tag+"_LJ_TT_dimu_met",eventbase->GetEvent().MET(snu::KEvent::pfmet) ,  weight, 0., 100., 20);
	    FillHist(tag+"_LJ_TT_dimu_njets",jets.size() ,  weight, 0., 10., 10);
	    
	    
	  }
	  else if(PassID(loosemuons[0], tightid)){
	    
	    //// TL events

	    float mu_pt_corr = loosemuons.at(1).Pt()*(1+max(0.,(loosemuons.at(1).RelIso04()-tightiso))) ;
	    
	    int bin= MuonFR->FindBin(loosemuons[1].Pt(), fabs(loosemuons[1].Eta()));
	    int bincorr=MuonFRcorr->FindBin(mu_pt_corr, fabs(loosemuons[1].Eta()));
	    float fw = MuonFR->GetBinContent(bin);
	    float fwcorr =MuonFRcorr->GetBinContent(bincorr);
	    fw=fw/(1.-fw);
	    fwcorr=fwcorr/(1.-fwcorr);
	    

	    FillHist(tag+"_LJ_TL_dimu_met",eventbase->GetEvent().MET(snu::KEvent::pfmet) ,  fw*weight, 0., 100., 20);
	    
	    FillHist(tag+"_LJ_TL_dimu", loosemuons[1].Pt() , fw*weight, 0., 100., 20);
	    FillHist(tag+"_LJ_TL_dimu_nw", loosemuons[1].Pt() , 1., 0., 100., 20);
	    FillHist(tag+"_LJ_TL_dimu_corr", mu_pt_corr, fwcorr*weight, 0., 100., 20);
            FillHist(tag+"_LJ_TL_dimu_njets",jets.size() ,  fw*weight, 0., 10., 10);
	    FillHist(tag+"_LJ_TL_dimu_corr_njets",jets.size() ,  fwcorr*weight, 0., 10., 10);
	    

	  }
	  else if(PassID(loosemuons[1], tightid)){
	    float mu_pt_corr = loosemuons.at(0).Pt()*(1+max(0.,(loosemuons.at(0).RelIso04()-tightiso))) ;
	    float mu_pt_corr1 = loosemuons.at(1).Pt()*(1+max(0.,(loosemuons.at(1).RelIso04()-tightiso))) ;
	    
	    int bin= MuonFR->FindBin(loosemuons[0].Pt(), fabs(loosemuons[0].Eta()));
	    int bincorr=MuonFRcorr->FindBin(mu_pt_corr, fabs(loosemuons[0].Eta()));
	    float fw = MuonFR->GetBinContent(bin);
	    float fwcorr =MuonFRcorr->GetBinContent(bincorr);
	    fw=fw/(1.-fw);
	    fwcorr=fwcorr/(1.-fwcorr);
	    
	    FillHist(tag+"_LJ_TL_dimu", loosemuons[1].Pt() , fw*weight, 0., 100., 20);
	    FillHist(tag+"_LJ_TL_dimu_corr", mu_pt_corr1 , fwcorr*weight, 0., 100., 20);
	    FillHist(tag+"_LJ_TL_dimu_nw", loosemuons[1].Pt() , 1., 0., 100., 20);
	    FillHist(tag+"_LJ_TL_dimu_met",eventbase->GetEvent().MET(snu::KEvent::pfmet)	,  fw*weight, 0., 100., 20);
	    FillHist(tag+"_LJ_TL_dimu_njets",jets.size() ,  fw*weight, 0., 10., 10);
            FillHist(tag+"_LJ_TL_dimu_corr_njets",jets.size() ,  fwcorr*weight, 0., 10., 10);
	    
	  }
	  else{
	    FillHist(tag+"_LJ_LL_dimu", loosemuons[1].Pt() , weight, 0., 100., 20);
	    FillHist(tag+"_LJ_LL_dimu_nw", loosemuons[1].Pt() , 1., 0., 100., 20);
	    
	  }
      }
    }
  }
  
  
    /*
	  bool closebjet = false;
	  for(unsigned int ij =0 ; ij < alljets.size() ; ij++){
	    if(loosemuons.at(1).DeltaR(alljets.at(ij)) < 0.5) {
	      if(alljets.at(ij).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) closebjet = true;
	    }
	  }

	  if(closebjet){
	    fw=MuonFRcbj->GetBinContent(bin);
	    fw=fw/(1.-fw);
	    fwcorr =MuonFRcbjcorr->GetBinContent(bincorr);
	    fwcorr=fwcorr/(1.-fwcorr);
	    cout << "close " << fw << " " << fwcorr << endl;
	    FillHist(tag+"_TL_dimu_cj", loosemuons[1].Pt() , weight*fw, 0., 100., 20);
	   FillHist(tag+"_TL_dimu_ccj_orr",mu_pt_corr, weight*fwcorr, 0., 100., 20);
	 
	  else{
	    fw=MuonFRncbj->GetBinContent(bin);
            fw=fw/(1.-fw);
            fwcorr =MuonFRncbjcorr->GetBinContent(bincorr);
            fwcorr=fwcorr/(1.-fwcorr);
	    
            FillHist(tag+"_TL_dimu_cj", loosemuons[1].Pt() , weight*fw, 0., 100., 20);
            FillHist(tag+"_TL_dimu_ccj_corr",mu_pt_corr, weight*fwcorr, 0., 100., 20);
            FillHist(tag+"_TL_dimu_ccj_corr_pt",loosemuons[1].Pt(), weight*fwcorr, 0., 100., 20);

	  }
	  if(!loosemuons[0].MCMatched()) FillHist(tag+"_MuonType_TL",loosemuons[1].GetType(),weight, 0., 41., 41);
	  if(!TruthMatched(loosemuons[1]))FillHist(tag+"_MuonType_TL",loosemuons[1].GetType(),weight, 0., 41., 41);
	  
	}
      }
    }
    */
  

  
  /// This is region used to measure fakes
  
  //if(loosemuons.size()!=1) return;
  if(loosemuons.size()==0) return;
  
  //  if(TruthMatch(loosemuons[0]))return;
  
  TString triggerslist_3="HLT_Mu3_PFJet40_v";
  TString triggerslist_8="HLT_Mu8_TrkIsoVVL_v";
  TString triggerslist_17="HLT_Mu17_TrkIsoVVL_v";
  
  float prescale_trigger =  GetPrescaleMu(loosemuons,   PassTrigger(triggerslist_3),PassTrigger(triggerslist_8), PassTrigger(triggerslist_17),TargetLumi);
  weight=prescale_trigger;

  for(unsigned int imu=0; imu < loosemuons.size(); imu++){

    bool closebjet = false;
    for(unsigned int ij =0 ; ij < alljets.size() ; ij++){
      if(loosemuons.at(imu).DeltaR(alljets.at(ij)) < 0.5) {
	if(alljets.at(ij).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) closebjet = true;
      }
    }
    FillHist(tag+"_MuonType",loosemuons[imu].GetType(),weight, 0., 41., 41);
    FillHist(tag+"_MuonMother",loosemuons[imu].MotherPdgId(),weight, 0., 600., 600);
    if(PassID(loosemuons[imu], tightid)) {
      FillHist(tag+"_MuonType_tight",loosemuons[imu].GetType(),weight, 0., 41., 41);
      FillHist(tag+"_MuonMother_tight",loosemuons[imu].MotherPdgId(),weight, 0., 600., 600);
    }
    float mu_pt_corr = loosemuons.at(imu).Pt()*(1+max(0.,(loosemuons.at(imu).RelIso04()-tightiso))) ; /// will need changing for systematics
    float mu_pt =  loosemuons.at(imu).Pt();

    bool useevent=false;
    float METdphi = TVector2::Phi_mpi_pi(loosemuons.at(0).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet));
    float MT = sqrt(2.* loosemuons.at(0).Et()*eventbase->GetEvent().MET(snu::KEvent::pfmet) * (1 - cos( METdphi)));
    if(( (eventbase->GetEvent().MET(snu::KEvent::pfmet) < 20) && (MT < 25.)) ) {

      for(unsigned int ij=0; ij < jets.size(); ij++){
	if(jets.at(ij).Pt() < 40.) continue;
	float dphi =fabs(TVector2::Phi_mpi_pi(loosemuons.at(imu).Phi()- jets.at(ij).Phi()));
	if(dphi > 2.5) useevent = true;
      }
    }
    if(useevent) {
      FillHist(tag+"_JetMuonType",loosemuons[imu].GetType(),weight, 0., 41., 41);
      FillHist(tag+"_JetMother",loosemuons[imu].MotherPdgId(),weight, 0., 600., 600);

      if(PassID(loosemuons[imu], tightid)) {
	FillHist(tag+"_JetMuonType_tight",loosemuons[imu].GetType(),weight, 0., 41., 41);
	FillHist(tag+"_JetMother_tight",loosemuons[imu].MotherPdgId(),weight, 0., 600., 600);
      }
    }
    

    if(fabs(loosemuons.at(imu).Eta() < 0.8)){
      FillHist(tag+"_Loose_mu_eb1_pt1D",mu_pt, 1., ptbins,10);
      FillHist(tag+"_Loose_mu_eb1_ptcorr1D",mu_pt_corr, 1., ptbins,10);
      if(useevent)       FillHist(tag+"_Loose_mu_eb1_ptcorr1D_full",mu_pt_corr, 1., ptbins,10);

      if(closebjet){
	FillHist(tag+"_Loose_mu_cj_eb1_pt1D",mu_pt, 1., ptbins,10);
      }
      else{
	FillHist(tag+"_Loose_mu_ncj_eb1_pt1D",mu_pt, 1., ptbins,10);
      }
      if(PassID(loosemuons[imu], tightid))   {
	FillHist(tag+"_Tight_mu_eb1_pt1D",mu_pt, 1., ptbins,10);
	FillHist(tag+"_Tight_mu_eb1_ptcorr1D",mu_pt_corr, 1., ptbins,10);
	if(useevent) FillHist(tag+"_Tight_mu_eb1_ptcorr1D_full",mu_pt_corr, 1., ptbins,10);
	if(closebjet){
	  FillHist(tag+"_Tight_mu_cj_eb1_pt1D",mu_pt, 1., ptbins,10);
	}
	else{
	  FillHist(tag+"_Tight_mu_ncj_eb1_pt1D",mu_pt, 1., ptbins,10);
	}

      }
    }
    else  if(fabs(loosemuons.at(imu).Eta() < 1.5)){
      FillHist(tag+"_Loose_mu_eb2_pt1D",mu_pt, 1., ptbins,10);
      FillHist(tag+"_Loose_mu_eb2_ptcorr1D",mu_pt_corr, 1., ptbins,10);
      if(useevent)FillHist(tag+"_Loose_mu_eb2_ptcorr1D_full",mu_pt_corr, 1., ptbins,10);
      if(closebjet){
        FillHist(tag+"_Loose_mu_cj_eb2_pt1D",mu_pt, 1., ptbins,10);
      } 
      else{
        FillHist(tag+"_Loose_mu_ncj_eb2_pt1D",mu_pt, 1., ptbins,10);
      } 

      if(PassID(loosemuons[imu], tightid))   {
	FillHist(tag+"_Tight_mu_eb2_pt1D",mu_pt, 1., ptbins,10);
        FillHist(tag+"_Tight_mu_eb2_ptcorr1D",mu_pt_corr, 1., ptbins,10);
        if(useevent)FillHist(tag+"_Tight_mu_eb2_ptcorr1D_full",mu_pt_corr, 1., ptbins,10);
	
        if(closebjet){
          FillHist(tag+"_Tight_mu_cj_eb2_pt1D",mu_pt, 1., ptbins,10);
        }
        else{
          FillHist(tag+"_Tight_mu_ncj_eb2_pt1D",mu_pt, 1., ptbins,10);
        }
      } 
    } 
    else{
      FillHist(tag+"_Loose_mu_ee_pt1D",mu_pt, 1., ptbins,10);
      FillHist(tag+"_Loose_mu_ee_ptcorr1D",mu_pt_corr, 1., ptbins,10);
      if(useevent)FillHist(tag+"_Loose_mu_ee_ptcorr1D_full",mu_pt_corr, 1., ptbins,10);
      if(closebjet){
        FillHist(tag+"_Loose_mu_cj_ee_pt1D",mu_pt, 1., ptbins,10);
      } 
      else{
        FillHist(tag+"_Loose_mu_ncj_ee_pt1D",mu_pt, 1., ptbins,10);
      } 

      if(PassID(loosemuons[imu], tightid))   {
	FillHist(tag+"_Tight_mu_ee_pt1D",mu_pt, 1., ptbins,10);
        FillHist(tag+"_Tight_mu_ee_ptcorr1D",mu_pt_corr, 1., ptbins,10);
        if(useevent)FillHist(tag+"_Tight_mu_ee_ptcorr1D_full",mu_pt_corr, 1., ptbins,10);
        if(closebjet){
          FillHist(tag+"_Tight_mu_cj_ee_pt1D",mu_pt, 1., ptbins,10);
        }
        else{
          FillHist(tag+"_Tight_mu_ncj_ee_pt1D",mu_pt, 1., ptbins,10);
        }
      } 
    }

    FillHist(tag+"_Loose_mu_pt",mu_pt,fabs(loosemuons.at(imu).Eta()), 1., ptbins,10, etabins2, 4);
    FillHist(tag+"_Loose_mu_ptcorr",mu_pt_corr,fabs(loosemuons.at(imu).Eta()), 1.,ptbins,10, etabins2, 4);
    if(closebjet){
      FillHist(tag+"_Loose_mu_cj_pt",mu_pt,fabs(loosemuons.at(imu).Eta()), 1., ptbins,10, etabins2, 4);
      FillHist(tag+"_Loose_mu_cj_ptcorr",mu_pt_corr,fabs(loosemuons.at(imu).Eta()), 1.,ptbins,10, etabins2, 4);
    }
    else{
      FillHist(tag+"_Loose_mu_ncj_pt",mu_pt,fabs(loosemuons.at(imu).Eta()), 1., ptbins,10, etabins2, 4);
      FillHist(tag+"_Loose_mu_ncj_ptcorr",mu_pt_corr,fabs(loosemuons.at(imu).Eta()), 1.,ptbins,10, etabins2, 4);

    }
    if(PassID(loosemuons[imu], tightid))   {
      FillHist(tag+"_Tight_mu_pt",mu_pt,fabs(loosemuons.at(imu).Eta()), 1., ptbins,10, etabins2, 4);
      FillHist(tag+"_Tight_mu_ptcorr",mu_pt_corr,fabs(loosemuons.at(imu).Eta()), 1., ptbins,10, etabins2, 4);
      if(closebjet){
	FillHist(tag+"_Tight_mu_cj_pt",mu_pt,fabs(loosemuons.at(imu).Eta()), 1., ptbins,10, etabins2, 4);
	FillHist(tag+"_Tight_mu_cj_ptcorr",mu_pt_corr,fabs(loosemuons.at(imu).Eta()), 1.,ptbins,10, etabins2, 4);
      }
      else{
	FillHist(tag+"_Tight_mu_ncj_pt",mu_pt,fabs(loosemuons.at(imu).Eta()), 1., ptbins,10, etabins2, 4);
	FillHist(tag+"_Tight_mu_ncj_ptcorr",mu_pt_corr,fabs(loosemuons.at(imu).Eta()), 1.,ptbins,10, etabins2, 4);

      }
    }

  }

  
  if(loosemuons.size()==1){
    float METdphi = TVector2::Phi_mpi_pi(loosemuons.at(0).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet));
    float MT = sqrt(2.* loosemuons.at(0).Et()*eventbase->GetEvent().MET(snu::KEvent::pfmet) * (1 - cos( METdphi)));
    
  }
}





float FakeRateMC::GetPrescaleMu( std::vector<snu::KMuon> muons,bool pass3, bool pass2, bool pass1, float fake_total_lum ){
  
  float prescale_trigger= 1.;
  if(muons.size() != 1) return 0.;
  if(muons.size() ==1){

    if(muons.at(0).Pt() >= 20.){
      if(pass1){
	return 1.;
      }
      else {
	return 0;
      }
    }
    else  if(muons.at(0).Pt() >= 10.){
      if(pass2){
	return 1.;
      }
      else {
	return 0;
      }
    }
    else  if(muons.at(0).Pt() >= 5.){
      if(pass3){
	return 1.;
      }
      else {
        if(isData) return 0;
	return 0;
      }
    }
  }
  if(prescale_trigger == 0.) return 0.;
  
  return prescale_trigger;
}





float FakeRateMC::GetPrescaleEl( std::vector<snu::KElectron> electrons,bool pass5,  bool pass4, bool pass3, bool pass2, bool pass1, float fake_total_lum ){
  
  float prescale_trigger= 1.;
  if(electrons.at(0).Pt() >= 35.){
    //HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30
    if(pass1){
      return 1;
    }
    else {
      return 0;
      //prescale_trigger = WeightByTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v", fake_total_lum)*0.8;
    }
  }
  else  if(electrons.at(0).Pt() >= 25.){
    //HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30

    if(pass2){
      return 1;
    }
    else {
      return 0;
    }
  }
  else   if(electrons.at(0).Pt() >= 20.){
    if(pass3){
      return 1;

    }
    else {
      return 0;
    }
  }
  else   if(electrons.at(0).Pt() >= 15.){
    //HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_
    if(pass4){
      return 1;
    }
    else {
      return 0;
    }
  }
  else   if(electrons.at(0).Pt() >= 8.){
    if(pass5){
      return 1;
    }
    else {
      return 0;
      //prescale_trigger = WeightByTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v", fake_total_lum)*0.8 ;
    }
  }
  else{
    prescale_trigger = 0.;
  }
  
  return 0;

}
void FakeRateMC::ExecuteEventsElectron(TString looseid, TString tightid, TString tag, float tightiso){

  
  std::vector<snu::KElectron> looseelectrons = GetElectrons(true,true,looseid);
  Float_t ptbins[11] = {5., 10., 15.,20.,25.,30.,35.,45.,60.,100., 200.};
  Float_t etabins2[5] = { 0.,0.8,  1.479, 2.,  2.5};

  std::vector<snu::KJet> alljets = GetJets("JET_NOLEPTONVETO");
  std::vector<snu::KJet> jets = GetJets("JET_HN");
  

  if(SameCharge(looseelectrons)){
    FillHist(tag+"SSee_mass", (looseelectrons[0]+looseelectrons[1]).M(), weight , 0., 200.,200);
    for(int x=0; x < looseelectrons.size(); x++){

      FillHist(tag+"_ElectronType_L_PF",looseelectrons[x].GetType(),weight, 0., 41., 41);
      FillHist(tag+"_ElectronType_L_PF_mother",fabs(looseelectrons[x].MotherPdgId()),weight, 0., 600., 600);
      if(PassID(looseelectrons[x], tightid) )  {
        FillHist(tag+"_ElectronType_L_PF_tight",looseelectrons[x].GetType(),weight, 0., 41., 41);
        FillHist(tag+"_ElectronType_L_PF_mother_tight",fabs(looseelectrons[x].MotherPdgId()),weight, 0., 600., 600);
      }
      if(!looseelectrons[x].MCMatched())  {
        FillHist(tag+"_ElectronType_LJ_PF",looseelectrons[x].GetType(),weight, 0., 41., 41);
        if(PassID(looseelectrons[x], tightid) )             FillHist(tag+"_ElectronType_LJ_PF_tight",looseelectrons[x].GetType(),weight, 0., 41., 41);

        if(looseelectrons[x].GetType()==2){
          FillHist(tag+"_ElectronType_LJ_PF_mother",fabs(looseelectrons[x].MotherPdgId()),weight, 0., 600., 600);
          if(PassID(looseelectrons[x], tightid))  FillHist(tag+"_ElectronType_LJ_PF_mother_tight",fabs(looseelectrons[x].MotherPdgId()),weight, 0., 600., 600);
        }
        FillHist(tag+"_ElectronType_mother_LJ_PF",looseelectrons[x].GetType(),fabs(looseelectrons[x].MotherPdgId()), weight, 0., 41., 41,  0., 600., 600);
        if(PassID(looseelectrons[x], tightid))FillHist(tag+"_ElectronType_mother_LJ_PF_tight",looseelectrons[x].GetType(),fabs(looseelectrons[x].MotherPdgId()), weight, 0., 41., 41,  0., 600.,600);
						       
      }
    }
  }

  if(looseelectrons.size()==3){
    FillHist(tag+"eee_mass", (looseelectrons[0]+looseelectrons[1]+looseelectrons[2]).M(), weight , 0., 200.,200);

    for(int x=0; x < looseelectrons.size(); x++){
      if(!looseelectrons[x].MCMatched())  {
        FillHist(tag+"_ElectronType_LLJ_PF",looseelectrons[x].GetType(),weight, 0., 41., 41);
        if(PassID(looseelectrons[x], tightid) )             FillHist(tag+"_ElectronType_LLJ_PF_tight",looseelectrons[x].GetType(),weight, 0., 41., 41);

        FillHist(tag+"_ElectronType_LLJ_PF_mother",fabs(looseelectrons[x].MotherPdgId()),weight, 0., 600., 600);
        if(PassID(looseelectrons[x], tightid))  FillHist(tag+"_ElectronType_LLJ_PF_mother_tight",fabs(looseelectrons[x].MotherPdgId()),weight, 0., 600., 600);

        FillHist(tag+"_ElectronType_mother_LLJ_PF",looseelectrons[x].GetType(),fabs(looseelectrons[x].MotherPdgId()), weight, 0., 41., 41,  0., 600., 600);
        if(PassID(looseelectrons[1], tightid))FillHist(tag+"_ElectronType_mother_LLJ_PF_tight",looseelectrons[x].GetType(),fabs(looseelectrons[x].MotherPdgId()), weight, 0., 41., 41,  0., 600., 600);
						       
      }
    }
  }


  if( ((k_sample_name.Contains("TT"))&&SameCharge(looseelectrons)) || (!k_sample_name.Contains("TT")&&looseelectrons.size()==2)){

    if(PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")){
      if(looseelectrons[0].Pt() > 25.){
	int nT=0;
	if(looseelectrons[0].MCMatched() ) nT++;
	if(looseelectrons[1].MCMatched() ) nT++;
	

	FillHist(tag+"_GetVirtualMass", GetVirtualMass() , weight, 0., 100., 100);
	if(PassID(looseelectrons[0], tightid))  {
	  if(PassID(looseelectrons[0], tightid)) {
	    FillHist(tag+"_GetVirtualMassTight", GetVirtualMass() , weight, 0., 100., 100);
	    
	  }
	}
	if(PassID(looseelectrons[0], tightid))  {
	  if(PassID(looseelectrons[1], tightid))  {
	    FillHist(tag+"_TT_diel_all", looseelectrons[1].Pt() , weight, 0., 100., 20);

	  }
	}
	if(PassID(looseelectrons[0], tightid)) {
	  if(!PassID(looseelectrons[1], tightid))  {
	    if(nT==1){
	      float el_pt_corr = looseelectrons.at(1).Pt()*(1+max(0.,(looseelectrons.at(1).PFRelIso(0.3)-tightiso))) ; 
	      int bin= ElFR->FindBin(looseelectrons[1].Pt(), fabs(looseelectrons[1].Eta()));
	      int bincorr=ElFRcorr->FindBin(el_pt_corr, fabs(looseelectrons[1].Eta()));
	      float fw = ElFR->GetBinContent(bin);
	      float fwcorr =ElFRcorr->GetBinContent(bincorr); 
	      fw=fw/(1.-fw);
	      fwcorr=fwcorr/(1.-fwcorr);
	      FillHist(tag+"_TL_diel", looseelectrons[1].Pt() , weight*fw, 0., 100., 20);
	      FillHist(tag+"_TL_diel_corr",el_pt_corr, weight*fwcorr, 0., 100., 20);
	      
	      
	      bool closebjet = false;
	      for(unsigned int ij =0 ; ij < alljets.size() ; ij++){
		if(looseelectrons.at(1).DeltaR(alljets.at(ij)) < 0.5) {
		  if(alljets.at(ij).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) closebjet = true;
		}
	      }
	      if(closebjet){
		fw=ElFRcbj->GetBinContent(bin);
		fw=fw/(1.-fw);
		fwcorr =ElFRcbjcorr->GetBinContent(bincorr);
		fwcorr=fwcorr/(1.-fwcorr);
		
		FillHist(tag+"_TL_diel_cj", looseelectrons[1].Pt() , weight*fw, 0., 100., 20);
		FillHist(tag+"_TL_diel_ccj_orr",el_pt_corr, weight*fwcorr, 0., 100., 20);
		
	      }
	      else{
		fw=ElFRncbj->GetBinContent(bin);
		fw=fw/(1.-fw);
		fwcorr =ElFRncbjcorr->GetBinContent(bincorr);
		fwcorr=fwcorr/(1.-fwcorr);
		
		FillHist(tag+"_TL_diel_cj", looseelectrons[1].Pt() , weight*fw, 0., 100., 20);
		FillHist(tag+"_TL_diel_ccj_orr",el_pt_corr, weight*fwcorr, 0., 100., 20);
		
	      }
	      
	      
	      
	    }
	  }
	}
	if(!PassID(looseelectrons[0], tightid)) {
	  if(PassID(looseelectrons[1], tightid))  {
	    if(nT==1){
	      
	      float el_pt_corr = looseelectrons.at(0).Pt()*(1+max(0.,(looseelectrons.at(0).PFRelIso(0.3)-tightiso))) ; 
	      float el_pt_corr1 = looseelectrons.at(1).Pt()*(1+max(0.,(looseelectrons.at(1).PFRelIso(0.3)-tightiso))) ; 
	      int bin= ElFR->FindBin(looseelectrons[0].Pt(), fabs(looseelectrons[0].Eta()));
	      int bincorr=ElFRcorr->FindBin(el_pt_corr, fabs(looseelectrons[0].Eta()));
	      float fw = ElFR->GetBinContent(bin);
	      float fwcorr =ElFRcorr->GetBinContent(bincorr);
	      
	      fw=fw/(1.-fw);
	      fwcorr=fwcorr/(1.-fwcorr);
	  
	      FillHist(tag+"_TL_diel", looseelectrons[1].Pt() , weight*fw, 0., 100., 20);
	      FillHist(tag+"_TL_diel_corr",el_pt_corr1, weight*fwcorr, 0., 100., 20);
	  

	      bool closebjet = false;
	      for(unsigned int ij =0 ; ij < alljets.size() ; ij++){
		if(looseelectrons.at(0).DeltaR(alljets.at(ij)) < 0.5) {
		  if(alljets.at(ij).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) closebjet = true;
            }
	      }
	      if(closebjet){
		fw=ElFRcbj->GetBinContent(bin);
		fw=fw/(1.-fw);
		fwcorr =ElFRcbjcorr->GetBinContent(bincorr);
		fwcorr=fwcorr/(1.-fwcorr);
		
		FillHist(tag+"_TL_diel_cj", looseelectrons[1].Pt() , weight*fw, 0., 100., 20);
		FillHist(tag+"_TL_diel_ccj_orr",el_pt_corr1, weight*fwcorr, 0., 100., 20);
		
	      }
	      else{
		fw=ElFRncbj->GetBinContent(bin);
		fw=fw/(1.-fw);
		fwcorr =ElFRncbjcorr->GetBinContent(bincorr);
		fwcorr=fwcorr/(1.-fwcorr);
		
		FillHist(tag+"_TL_diel_cj", looseelectrons[1].Pt() , weight*fw, 0., 100., 20);
		FillHist(tag+"_TL_diel_ccj_orr",el_pt_corr1, weight*fwcorr, 0., 100., 20);
		
	      }
	      
	    } 
	  }
	}
      }
    }
  }

  

  if(looseelectrons.size()!=1) return;


  TString triggerslist_8="HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v";   /// -> tighter cut in lepton ID form tighter trigger emulation cut                        
  TString triggerslist_12="HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  TString triggerslist_18="HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  TString triggerslist_23="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  TString triggerslist_33="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"; ///       

  float prescale_trigger =  GetPrescaleEl(looseelectrons,    PassTrigger(triggerslist_8), PassTrigger(triggerslist_12), PassTrigger(triggerslist_18), PassTrigger( triggerslist_23), PassTrigger(triggerslist_33), TargetLumi);								
  weight=prescale_trigger;
  


  for(unsigned int iel=0; iel < looseelectrons.size(); iel++){
    
    bool closebjet = false;
    for(unsigned int ij =0 ; ij < alljets.size() ; ij++){
      if(looseelectrons.at(iel).DeltaR(alljets.at(ij)) < 0.5) {
	if(alljets.at(ij).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) closebjet = true;
      }
    }
    
    FillHist(tag+"_ElectronType",looseelectrons[iel].GetType(),weight, 0., 41., 41);
    FillHist(tag+"_ElectronMother",looseelectrons[iel].MotherPdgId(),weight, 0., 600., 600);
    
    if(PassID(looseelectrons[iel], tightid))     FillHist(tag+"_ElectronType_tight",looseelectrons[iel].GetType(),weight, 0., 41., 41);
    float el_pt_corr = looseelectrons.at(iel).Pt()*(1+max(0.,(looseelectrons.at(iel).PFRelIso(0.3)-tightiso))) ; /// will need changing for systematics
    float el_pt =looseelectrons.at(iel).Pt();
    
    bool useevent=false;
    float METdphi = TVector2::Phi_mpi_pi(looseelectrons.at(0).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet));
    float MT = sqrt(2.* looseelectrons.at(0).Et()*eventbase->GetEvent().MET(snu::KEvent::pfmet) * (1 - cos( METdphi)));
    if(( (eventbase->GetEvent().MET(snu::KEvent::pfmet) < 20) && (MT < 25.)) ) {

      for(unsigned int ij=0; ij < jets.size(); ij++){
	if(jets.at(ij).Pt() < 40.) continue;
        float dphi =fabs(TVector2::Phi_mpi_pi(looseelectrons.at(iel).Phi()- jets.at(ij).Phi()));
        if(dphi > 2.5) useevent = true;
      }
    }
    if(useevent){
      FillHist(tag+"_JetElectronType",looseelectrons[iel].GetType(),weight, 0., 41., 41);
      if(PassID(looseelectrons[iel], tightid))     FillHist(tag+"_JetElectronType_tight",looseelectrons[iel].GetType(),weight, 0., 41., 41);

    }
    

    if(fabs(looseelectrons.at(iel).Eta() < 0.8)){
      FillHist(tag+"_Loose_el_eb1_pt1D",el_pt, 1., ptbins,10);
      FillHist(tag+"_Loose_el_eb1_ptcorr1D",el_pt_corr, 1., ptbins,10);
      if(useevent)       FillHist(tag+"_Loose_el_eb1_ptcorr1D_full",el_pt_corr, 1., ptbins,10);

      if(closebjet){
        FillHist(tag+"_Loose_el_cj_eb1_pt1D",el_pt, 1., ptbins,10);
      }
      else{
        FillHist(tag+"_Loose_el_ncj_eb1_pt1D",el_pt, 1., ptbins,10);
      }
      if(PassID(looseelectrons[iel], tightid))   {
        FillHist(tag+"_Tight_el_eb1_pt1D",el_pt, 1., ptbins,10);
        FillHist(tag+"_Tight_el_eb1_ptcorr1D",el_pt_corr, 1., ptbins,10);
        if(useevent) FillHist(tag+"_Tight_el_eb1_ptcorr1D_full",el_pt_corr, 1., ptbins,10);
        if(closebjet){
          FillHist(tag+"_Tight_el_cj_eb1_pt1D",el_pt, 1., ptbins,10);
        }
        else{
          FillHist(tag+"_Tight_el_ncj_eb1_pt1D",el_pt, 1., ptbins,10);
        }

      }
    }
    else  if(fabs(looseelectrons.at(iel).Eta() < 1.5)){
      FillHist(tag+"_Loose_el_eb2_pt1D",el_pt, 1., ptbins,10);
      FillHist(tag+"_Loose_el_eb2_ptcorr1D",el_pt_corr, 1., ptbins,10);
      if(useevent)FillHist(tag+"_Loose_el_eb2_ptcorr1D_full",el_pt_corr, 1., ptbins,10);
      if(closebjet){
        FillHist(tag+"_Loose_el_cj_eb2_pt1D",el_pt, 1., ptbins,10);
      }
      else{
        FillHist(tag+"_Loose_el_ncj_eb2_pt1D",el_pt, 1., ptbins,10);
      }

      if(PassID(looseelectrons[iel], tightid))   {
        FillHist(tag+"_Tight_el_eb2_pt1D",el_pt, 1., ptbins,10);
        FillHist(tag+"_Tight_el_eb2_ptcorr1D",el_pt_corr, 1., ptbins,10);
        if(useevent)FillHist(tag+"_Tight_el_eb2_ptcorr1D_full",el_pt_corr, 1., ptbins,10);

        if(closebjet){
          FillHist(tag+"_Tight_el_cj_eb2_pt1D",el_pt, 1., ptbins,10);
        }
	else{
          FillHist(tag+"_Tight_el_ncj_eb2_pt1D",el_pt, 1., ptbins,10);
        }
      }
    }
    else{
      FillHist(tag+"_Loose_el_ee_pt1D",el_pt, 1., ptbins,10);
      FillHist(tag+"_Loose_el_ee_ptcorr1D",el_pt_corr, 1., ptbins,10);
      if(useevent)FillHist(tag+"_Loose_el_ee_ptcorr1D_full",el_pt_corr, 1., ptbins,10);
      if(closebjet){
        FillHist(tag+"_Loose_el_cj_ee_pt1D",el_pt, 1., ptbins,10);
      }
      else{
        FillHist(tag+"_Loose_el_ncj_ee_pt1D",el_pt, 1., ptbins,10);
      }

      if(PassID(looseelectrons[iel], tightid))   {
        FillHist(tag+"_Tight_el_ee_pt1D",el_pt, 1., ptbins,10);
        FillHist(tag+"_Tight_el_ee_ptcorr1D",el_pt_corr, 1., ptbins,10);
        if(useevent)FillHist(tag+"_Tight_el_ee_ptcorr1D_full",el_pt_corr, 1., ptbins,10);
        if(closebjet){
          FillHist(tag+"_Tight_el_cj_ee_pt1D",el_pt, 1., ptbins,10);
        }
        else{
          FillHist(tag+"_Tight_el_ncj_ee_pt1D",el_pt, 1., ptbins,10);
        }
      }
    }


    FillHist(tag+"_Loose_el_pt",el_pt,fabs(looseelectrons.at(iel).Eta()), 1, ptbins,10 , etabins2, 4);
    FillHist(tag+"_Loose_el_ptcorr",el_pt_corr,fabs(looseelectrons.at(iel).Eta()), 1, ptbins,10 , etabins2, 4);
    if(closebjet) {
      FillHist(tag+"_Loose_el_cj_pt",el_pt,fabs(looseelectrons.at(iel).Eta()), 1, ptbins,10 , etabins2, 4);
      FillHist(tag+"_Loose_el_cj_ptcorr",el_pt_corr,fabs(looseelectrons.at(iel).Eta()), 1, ptbins,10 , etabins2, 4);
    }
    else{
      FillHist(tag+"_Loose_el_ncj_pt",el_pt,fabs(looseelectrons.at(iel).Eta()), 1, ptbins,10 , etabins2, 4);
      FillHist(tag+"_Loose_el_ncj_ptcorr",el_pt_corr,fabs(looseelectrons.at(iel).Eta()), 1, ptbins,10 , etabins2, 4);
      
    }
    if(PassID(looseelectrons[iel], tightid))   {
      FillHist(tag+"_Tight_el_pt",el_pt,fabs(looseelectrons.at(iel).Eta()), 1, ptbins,10 , etabins2, 4);
      FillHist(tag+"_Tight_el_ptcorr",el_pt_corr,fabs(looseelectrons.at(iel).Eta()), 1, ptbins,10 , etabins2, 4);
      if(closebjet) {
	FillHist(tag+"_Tight_el_cj_pt",el_pt,fabs(looseelectrons.at(iel).Eta()), 1, ptbins,10 , etabins2, 4);
	FillHist(tag+"_Tight_el_cj_ptcorr",el_pt_corr,fabs(looseelectrons.at(iel).Eta()), 1, ptbins,10 , etabins2, 4);
      }
      else{
	FillHist(tag+"_Tight_el_ncj_pt",el_pt,fabs(looseelectrons.at(iel).Eta()), 1, ptbins,10 , etabins2, 4);
	FillHist(tag+"_Tight_el_ncj_ptcorr",el_pt_corr,fabs(looseelectrons.at(iel).Eta()), 1, ptbins,10 , etabins2, 4);
	
      }
    }
    
  }
  
  FillHist(tag+"_JetElectronType",looseelectrons[0].GetType(),weight, 0., 41., 41);
  if(PassID(looseelectrons[0], tightid))     FillHist(tag+"_JetElectronType_tight",looseelectrons[0].GetType(),weight, 0., 41., 41);


}



void FakeRateMC::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger<< INFO << "Number of events that pass 1 7GeV trigger = " << n_17_pass  << LQLogger::endmsg;
  m_logger<< INFO << "Number of events that pass 17 GeV + jet trigger = " << n_17_jet_pass  << LQLogger::endmsg;
  m_logger<< INFO << "Number of events that pass 17 GeV || jet trigger = " << n_17_17_jet_pass  << LQLogger::endmsg;

}

void FakeRateMC::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //DeclareVariable(out_muons, "Signal_Muons");


  n_17_17_jet_pass=0;
  n_17_pass=0;

  
  return;
  
}

FakeRateMC::~FakeRateMC() {
  
  Message("In FakeRateMC Destructor" , INFO);
  
}



void FakeRateMC::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void FakeRateMC::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this FakeRateMCCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void FakeRateMC::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


