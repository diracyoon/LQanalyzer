/***************************************************************************
ak * @Project: LQFakeRateMC Frame - ROOT-based analysis framework for Korea SNU
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

  MuonFR =  dynamic_cast<TH2D*> (( file_fake_mumc->Get("FakeRateGENT__mu_pt"))->Clone());
  MuonFRcorr =  dynamic_cast<TH2D*> (( file_fake_mumc->Get("FakeRateGENT__mu_ptcorr"))->Clone());
  MuonFRcbj =  dynamic_cast<TH2D*> (( file_fake_mumc->Get("FakeRateGENT__mu_cj_pt"))->Clone());
  MuonFRcbjcorr =  dynamic_cast<TH2D*> (( file_fake_mumc->Get("FakeRateGENT__mu_cj_ptcorr"))->Clone());
  MuonFRncbj =  dynamic_cast<TH2D*> (( file_fake_mumc->Get("FakeRateGENT__mu_ncj_pt"))->Clone());
  MuonFRncbjcorr =  dynamic_cast<TH2D*> (( file_fake_mumc->Get("FakeRateGENT__mu_ncj_ptcorr"))->Clone());

  ElFR =  dynamic_cast<TH2D*> (( file_fake_elmc->Get("FakeRateGENT__el_pt"))->Clone());
  ElFRcorr =  dynamic_cast<TH2D*> (( file_fake_elmc->Get("FakeRateGENT__el_ptcorr"))->Clone());
  ElFRcbj =  dynamic_cast<TH2D*> (( file_fake_elmc->Get("FakeRateGENT__el_cj_pt"))->Clone());
  ElFRcbjcorr =  dynamic_cast<TH2D*> (( file_fake_elmc->Get("FakeRateGENT__el_cj_ptcorr"))->Clone());
  ElFRncbj =  dynamic_cast<TH2D*> (( file_fake_elmc->Get("FakeRateGENT__el_ncj_pt"))->Clone());
  ElFRncbjcorr =  dynamic_cast<TH2D*> (( file_fake_elmc->Get("FakeRateGENT__el_ncj_ptcorr"))->Clone());


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
  ExecuteEventsMuon("MUON_HN_LOOSE","MUON_HN_TIGHT","HN",0.07,weight);
  ExecuteEventsMuon("MUON_HN_LOOSE","MUON_HN_TIGHT","HNNoWeight",0.07,1.);

  //ExecuteEventsMuon("MUON_HNGENT_LOOSE","MUON_HNGENT_TIGHT","GENT", 0.1,weight);
  ExecuteEventsElectron("ELECTRON_HN_FAKELOOSE_mediumMVA","ELECTRON_HN_TIGHTv4","HNmediummva",0.08,weight);
  ExecuteEventsElectron("ELECTRON_HN_FAKELOOSE_noMVA","ELECTRON_HN_TIGHTv4","HNnomva",0.08,weight);
  ExecuteEventsElectron("ELECTRON_HN_FAKELOOSE","ELECTRON_HN_TIGHTv4","HN",0.08,weight);
  ExecuteEventsElectron("ELECTRON_MVA_FAKELOOSE","ELECTRON_MVA_TIGHTJB","MVA",0.1,weight);
  //ExecuteEventsElectron("ELECTRON_GENT_FAKELOOSE","ELECTRON_GENT_TIGHT","GENT",0.1,weight);
}


void FakeRateMC::ExecuteEventsMuon(TString looseid, TString tightid, TString tag, float tightiso, double wmu){
  
  std::vector<snu::KMuon> loosemuons = GetMuons(looseid,true);
  Float_t ptbins[11] = { 5., 10., 15.,20.,25.,30.,35.,45.,60.,100., 200.};
  Float_t etabins2[5] = { 0.,0.8,  1.479, 2.,  2.5};
  std::vector<snu::KJet> alljets = GetJets("JET_NOLEPTONVETO", 10., 2.5);
  std::vector<snu::KJet> jets = GetJets("JET_HN");

  if(SameCharge(loosemuons)){
    FillHist(tag+"SSmu_mass", (loosemuons[0]+loosemuons[1]).M(), wmu , 0., 200., 200);
    
    for(int x=0; x < loosemuons.size(); x++){
      if(!loosemuons[x].MCMatched())  {
	int jettype =  CloseJetType(loosemuons.at(x), alljets);
	FillHist(tag + "_fake_loose_ss_muon_jetflavour", jettype, wmu , -4., 4., 8);
      }
    }
  }
  for(int x=0; x < loosemuons.size(); x++){
    FillHist(tag+"_MuonType_L_PF",loosemuons[x].GetType(),wmu, 0., 41., 41);
    FillHist(tag+"_MuonType_L_PF_mother",fabs(loosemuons[x].MotherPdgId()),wmu, 0., 2000., 2000);
    if(PassID(loosemuons[x], tightid) )  {
      FillHist(tag+"_MuonType_L_PF_tight",loosemuons[x].GetType(),wmu, 0., 41., 41);
      FillHist(tag+"_MuonType_L_PF_mother_tight",fabs(loosemuons[x].MotherPdgId()),wmu, 0., 2000., 2000);
    }
    
    if(!loosemuons[x].MCMatched())  {

      
      FillHist(tag+"_MuonType_LJ_PF",loosemuons[x].GetType(),wmu, 0., 41., 41);
      if(PassID(loosemuons[x], tightid) )             FillHist(tag+"_MuonType_LJ_PF_tight",loosemuons[x].GetType(),wmu, 0., 41., 41);
      
      if(loosemuons[x].GetType()==2){
	FillHist(tag+"_MuonType_LJ_PF_mother",fabs(loosemuons[x].MotherPdgId()),wmu, 0., 2000., 2000);
	if(PassID(loosemuons[x], tightid))  FillHist(tag+"_MuonType_LJ_PF_mother_tight",fabs(loosemuons[x].MotherPdgId()),wmu, 0., 2000., 2000);
      }
      FillHist(tag+"_MuonType_mother_LJ_PF",loosemuons[x].GetType(),fabs(loosemuons[x].MotherPdgId()), wmu, 0., 41., 41,  0., 2000., 2000);
      if(PassID(loosemuons[x], tightid))FillHist(tag+"_MuonType_mother_LJ_PF_tight",loosemuons[x].GetType(),fabs(loosemuons[x].MotherPdgId()), wmu, 0., 41., 41,  0., 2000., 2000);
    }
  }
  
  if( ((k_sample_name.Contains("TT"))&&SameCharge(loosemuons)) || (!k_sample_name.Contains("TT")&&loosemuons.size()==2)){
    if(PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") || PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")){
      if(loosemuons[0].Pt() > 20.&&loosemuons[1].Pt() > 10){
	int nT=0;
	if(TruthMatched(loosemuons[0]))nT++;
	if(TruthMatched(loosemuons[1]))nT++;
	if(!(NonPrompt(loosemuons[0]) || NonPrompt(loosemuons[1])))nT=0;
	if(nT==1){
	  //// Check both are tight
	  if(PassID(loosemuons[0], tightid) && (PassID(loosemuons[1], tightid))){
	    FillHist(tag+"_LJ_TT_dimu", loosemuons[1].Pt() , wmu, 0., 100., 20);
	    FillHist(tag+"_LJ_TT_dimu_nw", loosemuons[1].Pt() , 1., 0., 100., 20);
	    float mu_pt_corr = loosemuons.at(1).Pt()*(1+max(0.,(loosemuons.at(1).RelIso04()-tightiso))) ;
	    
	    FillHist(tag+"_LJ_TT_dimu_corr", mu_pt_corr  , wmu, 0., 100., 20);
	    FillHist(tag+"_LJ_TT_dimu_met",eventbase->GetEvent().MET(snu::KEvent::pfmet) ,  wmu, 0., 100., 20);
	    FillHist(tag+"_LJ_TT_dimu_njets",jets.size() ,  wmu, 0., 10., 10);
	    
	    
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
	    

	    FillHist(tag+"_LJ_TL_dimu_met",eventbase->GetEvent().MET(snu::KEvent::pfmet) ,  fw*wmu, 0., 100., 20);
	    
	    FillHist(tag+"_LJ_TL_dimu", loosemuons[1].Pt() , fw*wmu, 0., 100., 20);
	    FillHist(tag+"_LJ_TL_dimu_nw", loosemuons[1].Pt() , 1., 0., 100., 20);
	    FillHist(tag+"_LJ_TL_dimu_corr", mu_pt_corr, fwcorr*wmu, 0., 100., 20);
            FillHist(tag+"_LJ_TL_dimu_njets",jets.size() ,  fw*wmu, 0., 10., 10);
	    FillHist(tag+"_LJ_TL_dimu_corr_njets",jets.size() ,  fwcorr*wmu, 0., 10., 10);
	    
            if(!TruthMatched(loosemuons[1]))FillHist(tag+"_MuonType_TL",loosemuons[1].GetType(),wmu, 0., 41., 41);                                                    

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
	    
	    FillHist(tag+"_LJ_TL_dimu", loosemuons[1].Pt() , fw*wmu, 0., 100., 20);
	    FillHist(tag+"_LJ_TL_dimu_corr", mu_pt_corr1 , fwcorr*wmu, 0., 100., 20);
	    FillHist(tag+"_LJ_TL_dimu_nw", loosemuons[1].Pt() , 1., 0., 100., 20);
	    FillHist(tag+"_LJ_TL_dimu_met",eventbase->GetEvent().MET(snu::KEvent::pfmet)	,  fw*wmu, 0., 100., 20);
	    FillHist(tag+"_LJ_TL_dimu_njets",jets.size() ,  fw*wmu, 0., 10., 10);
            FillHist(tag+"_LJ_TL_dimu_corr_njets",jets.size() ,  fwcorr*wmu, 0., 10., 10);
	    if(!TruthMatched(loosemuons[0])){
	      FillHist(tag+"_MuonType_TL",loosemuons[0].GetType(),wmu, 0., 41., 41);                                                                 
	    }
	  }
	  else{
	    FillHist(tag+"_LJ_LL_dimu", loosemuons[1].Pt() , wmu, 0., 100., 20);
	    FillHist(tag+"_LJ_LL_dimu_nw", loosemuons[1].Pt() , 1., 0., 100., 20);
	    
	  }
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
	    FillHist(tag+"_TL_dimu_cj", loosemuons[1].Pt() , wmu*fw, 0., 100., 20);
	   FillHist(tag+"_TL_dimu_ccj_orr",mu_pt_corr, wmu*fwcorr, 0., 100., 20);
	 
	  else{
	    fw=MuonFRncbj->GetBinContent(bin);
            fw=fw/(1.-fw);
            fwcorr =MuonFRncbjcorr->GetBinContent(bincorr);
            fwcorr=fwcorr/(1.-fwcorr);
	    
            FillHist(tag+"_TL_dimu_cj", loosemuons[1].Pt() , wmu*fw, 0., 100., 20);
            FillHist(tag+"_TL_dimu_ccj_corr",mu_pt_corr, wmu*fwcorr, 0., 100., 20);
            FillHist(tag+"_TL_dimu_ccj_corr_pt",loosemuons[1].Pt(), wmu*fwcorr, 0., 100., 20);

	  }
	  if(!loosemuons[0].MCMatched()) FillHist(tag+"_MuonType_TL",loosemuons[1].GetType(),wmu, 0., 41., 41);
	  if(!TruthMatched(loosemuons[1]))FillHist(tag+"_MuonType_TL",loosemuons[1].GetType(),wmu, 0., 41., 41);
	  
	}
      }
    }
    */
  

  
  /// This is region used to measure fakes
  
  //if(loosemuons.size()!=1) return;
  if(loosemuons.size()==0) return;
  vector<snu::KMuon> fake_muons;
  for(unsigned int imu=0; imu < loosemuons.size(); imu++){
    if(loosemuons[imu].GetType()==2)fake_muons.push_back(loosemuons[imu]);
  }
  if(fake_muons.size() != 1) return;
  
  TString triggerslist_3="HLT_Mu3_PFJet40_v";
  TString triggerslist_8="HLT_Mu8_TrkIsoVVL_v";
  TString triggerslist_17="HLT_Mu17_TrkIsoVVL_v";
  
  for(unsigned int imu=0; imu < fake_muons.size(); imu++){

    bool closebjet = false;
    for(unsigned int ij =0 ; ij < alljets.size() ; ij++){
      if(fake_muons.at(imu).DeltaR(alljets.at(ij)) < 0.5) {
	if(alljets.at(ij).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) closebjet = true;
      }
    }
    FillHist(tag+"_MuonType",fake_muons[imu].GetType(),wmu, 0., 41., 41);
    FillHist(tag+"_MuonMother",fake_muons[imu].MotherPdgId(),wmu, 0., 2000., 2000);
    if(PassID(fake_muons[imu], tightid)) {
      FillHist(tag+"_MuonType_tight",fake_muons[imu].GetType(),wmu, 0., 41., 41);
      FillHist(tag+"_MuonMother_tight",fake_muons[imu].MotherPdgId(),wmu, 0., 2000., 2000);
    }
    float mu_pt_corr = fake_muons.at(imu).Pt()*(1+max(0.,(fake_muons.at(imu).RelIso04()-tightiso))) ; /// will need changing for systematics
    float mu_pt =  fake_muons.at(imu).Pt();

    bool useevent=false;
    float METdphi = TVector2::Phi_mpi_pi(fake_muons.at(0).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet));
    float MT = sqrt(2.* fake_muons.at(0).Et()*eventbase->GetEvent().MET(snu::KEvent::pfmet) * (1 - cos( METdphi)));
    if(( (eventbase->GetEvent().MET(snu::KEvent::pfmet) < 20) && (MT < 25.)) ) {

      for(unsigned int ij=0; ij < jets.size(); ij++){
	if(jets.at(ij).Pt() < 40.) continue;
	float dphi =fabs(TVector2::Phi_mpi_pi(fake_muons.at(imu).Phi()- jets.at(ij).Phi()));
	if(dphi > 2.5) useevent = true;
      }
    }
    if(useevent) {
      FillHist(tag+"_JetMuonType",fake_muons[imu].GetType(),wmu, 0., 41., 41);
      FillHist(tag+"_JetMother",fake_muons[imu].MotherPdgId(),wmu, 0., 2000., 2000);

      if(PassID(fake_muons[imu], tightid)) {
	FillHist(tag+"_JetMuonType_tight",fake_muons[imu].GetType(),wmu, 0., 41., 41);
	FillHist(tag+"_JetMother_tight",fake_muons[imu].MotherPdgId(),wmu, 0., 2000., 2000);
      }
    }
    

    if(fabs(fake_muons.at(imu).Eta() < 0.8)){
      FillHist(tag+"_Loose_mu_eb1_pt1D",mu_pt, wmu, ptbins,10);
      FillHist(tag+"_Loose_mu_eb1_ptcorr1D",mu_pt_corr, wmu, ptbins,10);
      if(useevent)       FillHist(tag+"_Loose_mu_eb1_ptcorr1D_full",mu_pt_corr, wmu, ptbins,10);

      if(closebjet){
	FillHist(tag+"_Loose_mu_cj_eb1_pt1D",mu_pt, wmu, ptbins,10);
      }
      else{
	FillHist(tag+"_Loose_mu_ncj_eb1_pt1D",mu_pt, wmu, ptbins,10);
      }
      if(PassID(fake_muons[imu], tightid))   {
	FillHist(tag+"_Tight_mu_eb1_pt1D",mu_pt, wmu, ptbins,10);
	FillHist(tag+"_Tight_mu_eb1_ptcorr1D",mu_pt_corr, wmu, ptbins,10);
	if(useevent) FillHist(tag+"_Tight_mu_eb1_ptcorr1D_full",mu_pt_corr, wmu, ptbins,10);
	if(closebjet){
	  FillHist(tag+"_Tight_mu_cj_eb1_pt1D",mu_pt, wmu, ptbins,10);
	}
	else{
	  FillHist(tag+"_Tight_mu_ncj_eb1_pt1D",mu_pt, wmu, ptbins,10);
	}

      }
    }

    else  if(fabs(fake_muons.at(imu).Eta() < 1.5)){
      FillHist(tag+"_Loose_mu_eb2_pt1D",mu_pt, wmu, ptbins,10);
      FillHist(tag+"_Loose_mu_eb2_ptcorr1D",mu_pt_corr, wmu, ptbins,10);
      if(useevent)FillHist(tag+"_Loose_mu_eb2_ptcorr1D_full",mu_pt_corr, wmu, ptbins,10);
      if(closebjet){
        FillHist(tag+"_Loose_mu_cj_eb2_pt1D",mu_pt, wmu, ptbins,10);
      } 
      else{
        FillHist(tag+"_Loose_mu_ncj_eb2_pt1D",mu_pt, wmu, ptbins,10);
      } 

      if(PassID(fake_muons[imu], tightid))   {
	FillHist(tag+"_Tight_mu_eb2_pt1D",mu_pt, wmu, ptbins,10);
        FillHist(tag+"_Tight_mu_eb2_ptcorr1D",mu_pt_corr, wmu, ptbins,10);
        if(useevent)FillHist(tag+"_Tight_mu_eb2_ptcorr1D_full",mu_pt_corr, wmu, ptbins,10);
	
        if(closebjet){
          FillHist(tag+"_Tight_mu_cj_eb2_pt1D",mu_pt, wmu, ptbins,10);
        }
        else{
          FillHist(tag+"_Tight_mu_ncj_eb2_pt1D",mu_pt, wmu, ptbins,10);
        }
      } 
    } 
    else{
      FillHist(tag+"_Loose_mu_ee_pt1D",mu_pt, wmu, ptbins,10);
      FillHist(tag+"_Loose_mu_ee_ptcorr1D",mu_pt_corr, wmu, ptbins,10);
      if(useevent)FillHist(tag+"_Loose_mu_ee_ptcorr1D_full",mu_pt_corr, wmu, ptbins,10);
      if(closebjet){
        FillHist(tag+"_Loose_mu_cj_ee_pt1D",mu_pt, wmu, ptbins,10);
      } 
      else{
        FillHist(tag+"_Loose_mu_ncj_ee_pt1D",mu_pt, wmu, ptbins,10);
      } 

      if(PassID(fake_muons[imu], tightid))   {
	FillHist(tag+"_Tight_mu_ee_pt1D",mu_pt, wmu, ptbins,10);
        FillHist(tag+"_Tight_mu_ee_ptcorr1D",mu_pt_corr, wmu, ptbins,10);
        if(useevent)FillHist(tag+"_Tight_mu_ee_ptcorr1D_full",mu_pt_corr, wmu, ptbins,10);
        if(closebjet){
          FillHist(tag+"_Tight_mu_cj_ee_pt1D",mu_pt, wmu, ptbins,10);
        }
        else{
          FillHist(tag+"_Tight_mu_ncj_ee_pt1D",mu_pt, wmu, ptbins,10);
        }
      } 
    }

    
    int jettype =  CloseJetType(fake_muons.at(imu), alljets);

    FillHist(tag + "_fake_loose_muon_jetflavour", jettype, wmu , -4., 4., 8);
    if(jettype==1) {
      FillHist(tag+"_Loose_mu_b_pt",mu_pt,fabs(fake_muons.at(imu).Eta()), wmu, ptbins,10, etabins2, 4);
      FillHist(tag+"_Loose_mu_b_ptcorr",mu_pt_corr,fabs(fake_muons.at(imu).Eta()), wmu,ptbins,10, etabins2, 4);
      if(PassID(fake_muons[imu], tightid))   {
	FillHist(tag+"_Tight_mu_b_pt",mu_pt,fabs(fake_muons.at(imu).Eta()), wmu, ptbins,10, etabins2, 4);
	FillHist(tag+"_Tight_mu_b_ptcorr",mu_pt_corr,fabs(fake_muons.at(imu).Eta()), wmu, ptbins,10, etabins2, 4);
      }
    }
    if(jettype==2) {
      FillHist(tag+"_Loose_mu_c_pt",mu_pt,fabs(fake_muons.at(imu).Eta()), wmu, ptbins,10, etabins2, 4);
      FillHist(tag+"_Loose_mu_c_ptcorr",mu_pt_corr,fabs(fake_muons.at(imu).Eta()), wmu,ptbins,10, etabins2, 4);
      if(PassID(fake_muons[imu], tightid))   {
        FillHist(tag+"_Tight_mu_c_pt",mu_pt,fabs(fake_muons.at(imu).Eta()), wmu, ptbins,10, etabins2, 4);
        FillHist(tag+"_Tight_mu_c_ptcorr",mu_pt_corr,fabs(fake_muons.at(imu).Eta()), wmu, ptbins,10, etabins2, 4);
      }
    } 
    if(jettype==3) {
      FillHist(tag+"_Loose_mu_l_pt",mu_pt,fabs(fake_muons.at(imu).Eta()), wmu, ptbins,10, etabins2, 4);
      FillHist(tag+"_Loose_mu_l_ptcorr",mu_pt_corr,fabs(fake_muons.at(imu).Eta()), wmu,ptbins,10, etabins2, 4);
      if(PassID(fake_muons[imu], tightid))   {
        FillHist(tag+"_Tight_mu_l_pt",mu_pt,fabs(fake_muons.at(imu).Eta()), wmu, ptbins,10, etabins2, 4);
        FillHist(tag+"_Tight_mu_l_ptcorr",mu_pt_corr,fabs(fake_muons.at(imu).Eta()), wmu, ptbins,10, etabins2, 4);
      }
    }
    



    FillHist(tag+"_Loose_mu_pt",mu_pt,fabs(fake_muons.at(imu).Eta()), wmu, ptbins,10, etabins2, 4);
    FillHist(tag+"_Loose_mu_ptcorr",mu_pt_corr,fabs(fake_muons.at(imu).Eta()), wmu,ptbins,10, etabins2, 4);
    if(closebjet){
      FillHist(tag+"_Loose_mu_cj_pt",mu_pt,fabs(fake_muons.at(imu).Eta()), wmu, ptbins,10, etabins2, 4);
      FillHist(tag+"_Loose_mu_cj_ptcorr",mu_pt_corr,fabs(fake_muons.at(imu).Eta()), wmu,ptbins,10, etabins2, 4);
    }
    else{
      FillHist(tag+"_Loose_mu_ncj_pt",mu_pt,fabs(fake_muons.at(imu).Eta()), wmu, ptbins,10, etabins2, 4);
      FillHist(tag+"_Loose_mu_ncj_ptcorr",mu_pt_corr,fabs(fake_muons.at(imu).Eta()), wmu,ptbins,10, etabins2, 4);

    }
    if(PassID(fake_muons[imu], tightid))   {
      FillHist(tag+"_Tight_mu_pt",mu_pt,fabs(fake_muons.at(imu).Eta()), wmu, ptbins,10, etabins2, 4);
      FillHist(tag+"_Tight_mu_ptcorr",mu_pt_corr,fabs(fake_muons.at(imu).Eta()), wmu, ptbins,10, etabins2, 4);
      if(closebjet){
	FillHist(tag+"_Tight_mu_cj_pt",mu_pt,fabs(fake_muons.at(imu).Eta()), wmu, ptbins,10, etabins2, 4);
	FillHist(tag+"_Tight_mu_cj_ptcorr",mu_pt_corr,fabs(fake_muons.at(imu).Eta()), wmu,ptbins,10, etabins2, 4);
      }
      else{
	FillHist(tag+"_Tight_mu_ncj_pt",mu_pt,fabs(fake_muons.at(imu).Eta()), wmu, ptbins,10, etabins2, 4);
	FillHist(tag+"_Tight_mu_ncj_ptcorr",mu_pt_corr,fabs(fake_muons.at(imu).Eta()), wmu,ptbins,10, etabins2, 4);

      }
    }

  }

  
  if(fake_muons.size()==1){
    float METdphi = TVector2::Phi_mpi_pi(fake_muons.at(0).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet));
    float MT = sqrt(2.* fake_muons.at(0).Et()*eventbase->GetEvent().MET(snu::KEvent::pfmet) * (1 - cos( METdphi)));
    
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
      //prescale_trigger = WmuByTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v", fake_total_lum)*0.8;
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
      //prescale_trigger = WmuByTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v", fake_total_lum)*0.8 ;
    }
  }
  else{
    prescale_trigger = 0.;
  }
  
  return 0;

}
void FakeRateMC::ExecuteEventsElectron(TString looseid, TString tightid, TString tag, float tightiso, double wel){

  
  std::vector<snu::KElectron> looseelectrons ;
  std::vector<snu::KElectron> looseelectronstmp = GetElectrons(true,true,looseid);
 
  for(unsigned int iel = 0; iel < looseelectronstmp.size(); iel++){
    
    bool pass_trigger_emulation=true;
    
    if(looseelectronstmp[iel].Pt() < 15.){
      if(!looseelectronstmp[iel].PassHLTID()) pass_trigger_emulation=false;
    }
    else {
      if(!looseelectronstmp[iel].IsTrigMVAValid()) pass_trigger_emulation=false;
    }
    if(pass_trigger_emulation) looseelectrons.push_back(looseelectronstmp[iel]);
  }


  Float_t ptbins[11] = {5., 10., 15.,20.,25.,30.,35.,45.,60.,100., 200.};
  Float_t etabins2[5] = { 0.,0.8,  1.479, 2.,  2.5};


  std::vector<snu::KJet> alljets = GetJets("JET_NOLEPTONVETO", 10., 2.5);
  std::vector<snu::KJet> jets = GetJets("JET_HN");
  

  if(IsFakeEvent(looseelectrons)){
    if(looseelectrons.size() > 1)FillHist(tag+"SSee_mass", (looseelectrons[0]+looseelectrons[1]).M(), wel , 0., 200.,200);
    
    for(int x=0; x < looseelectrons.size(); x++){
      
      int jettype =  CloseJetTypeHAD(looseelectrons.at(x), alljets);
      
      TString tmpstring="";
      snu::KElectron el = looseelectrons.at(x);
      
      std::vector<snu::KTruth> truthColl= eventbase->GetTruth();
      int LepType=GetLeptonType(looseelectrons.at(x),truthColl);
      
      
      FillHist(tag+"_ElectronType_L_PF",looseelectrons[x].GetType(),wel, 0., 41., 41);
      FillHist(tag+"_ElectronType_L_PF_mother",fabs(looseelectrons[x].MotherPdgId()),wel, 0., 2000., 2000);
      
      if(PassID(looseelectrons[x], tightid) )  {
	FillHist(tag+"_ElectronType_L_PF_tight",looseelectrons[x].GetType(),wel, 0., 41., 41);
	FillHist(tag+"_ElectronType_L_PF_mother_tight",fabs(looseelectrons[x].MotherPdgId()),wel, 0., 2000., 2000);
      }
      if(!looseelectrons[x].MCMatched())  {
	FillHist(tag + "_fake_loose_jetflavour", jettype, wel , -4., 4., 8); 
	if(PassID(looseelectrons[x], tightid) )          FillHist(tag + "_fake_tight_jetflavour", jettype, wel , -4., 4., 8);

	bool closebjet = false;
	for(unsigned int ij =0 ; ij < alljets.size() ; ij++){
	  if(looseelectrons[x].DeltaR(alljets.at(ij)) < 0.5) {
	    if(alljets.at(ij).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) closebjet = true;
	  }
	}
	if(closebjet) FillHist(tag + "_close_bjet", 1, wel, 2, 0., 2.);
	else  FillHist(tag + "_close_bjet", 0, wel, 2, 0., 2.);
	
	if(closebjet){
	  FillHist(tag + "_fake_closebjet_jetflavour", jettype, wel , -4., 4., 8);
	}
	else           FillHist(tag + "_fake_noclosebjet_jetflavour", jettype, wel , -4., 4., 8);


	if(LepType < 0){
	  if(jettype==1)          FillHist(tag+"_ElectronType_LJ_PF_b",looseelectrons[x].GetType(),wel, 0., 41., 41);
	  if(jettype==2)          FillHist(tag+"_ElectronType_LJ_PF_c",looseelectrons[x].GetType(),wel, 0., 41., 41);	
	  if(jettype==3)          FillHist(tag+"_ElectronType_LJ_PF_l",looseelectrons[x].GetType(),wel, 0., 41., 41);
	  FillHist(tag+"_ElectronType_LJ_PF",looseelectrons[x].GetType(),wel, 0., 41., 41);
	  if(PassID(looseelectrons[x], tightid) )             FillHist(tag+"_ElectronType_LJ_PF_tight",looseelectrons[x].GetType(),wel, 0., 41., 41);
	  
	  if(looseelectrons[x].GetType()==7){
	    FillHist(tag+"_ElectronType_LJ_PF_mother",fabs(looseelectrons[x].MotherPdgId()),wel, 0., 2000., 2000);
	    if(PassID(looseelectrons[x], tightid))  FillHist(tag+"_ElectronType_LJ_PF_mother_tight",fabs(looseelectrons[x].MotherPdgId()),wel, 0., 2000., 2000);
	  }
	  FillHist(tag+"_ElectronType_mother_LJ_PF",looseelectrons[x].GetType(),fabs(looseelectrons[x].MotherPdgId()), wel, 0., 41., 41,  0., 2000., 2000);
	  if(PassID(looseelectrons[x], tightid))FillHist(tag+"_ElectronType_mother_LJ_PF_tight",looseelectrons[x].GetType(),fabs(looseelectrons[x].MotherPdgId()), wel, 0., 41., 41,  0., 600.,600);
	  
	}
      }
    }
  }
  
  
  if( ((k_sample_name.Contains("TT"))&&SameCharge(looseelectrons)) || (!k_sample_name.Contains("TT")&&looseelectrons.size()==2)){

    if(PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")){
      if(looseelectrons[0].Pt() > 25. && looseelectrons[1].Pt() > 15.){
	int nT=0;
	if(TruthMatched(looseelectrons[0],false) ) nT++;
	if(TruthMatched(looseelectrons[1],false) ) nT++;
        if(!(NonPrompt(looseelectrons[0]) || NonPrompt(looseelectrons[1])))nT=0;

	if(nT==1){
	  
	  FillHist(tag+"_GetVirtualMass", GetVirtualMass() , wel, 0., 100., 100);
	  if(PassID(looseelectrons[0], tightid))  {
	    if(PassID(looseelectrons[0], tightid)) {
	      FillHist(tag+"_GetVirtualMassTight", GetVirtualMass() , wel, 0., 100., 100);
	      
	    }
	  }
	  if(PassID(looseelectrons[0], tightid))  {
	    if(PassID(looseelectrons[1], tightid))  {
	      FillHist(tag+"_TT_diel_all", looseelectrons[1].Pt() , wel, 0., 100., 20);
	      
	    }
	  }
	  if(PassID(looseelectrons[0], tightid)) {
	    if(!PassID(looseelectrons[1], tightid))  {
	      float el_pt_corr = looseelectrons.at(1).Pt()*(1+max(0.,(looseelectrons.at(1).PFRelIso(0.3)-tightiso))) ; 
	      int bin= ElFR->FindBin(looseelectrons[1].Pt(), fabs(looseelectrons[1].Eta()));
	      int bincorr=ElFRcorr->FindBin(el_pt_corr, fabs(looseelectrons[1].Eta()));
	      float fw = ElFR->GetBinContent(bin);
	      float fwcorr =ElFRcorr->GetBinContent(bincorr); 
	      fw=fw/(1.-fw);
	      fwcorr=fwcorr/(1.-fwcorr);
	      FillHist(tag+"_TL_diel", looseelectrons[1].Pt() , wel*fw, 0., 100., 20);
	      FillHist(tag+"_TL_diel_corr",el_pt_corr, wel*fwcorr, 0., 100., 20);
	      
	      
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
		
		FillHist(tag+"_TL_diel_cj", looseelectrons[1].Pt() , wel*fw, 0., 100., 20);
		FillHist(tag+"_TL_diel_ccj_orr",el_pt_corr, wel*fwcorr, 0., 100., 20);
		
	      }
	      else{
		fw=ElFRncbj->GetBinContent(bin);
		fw=fw/(1.-fw);
		fwcorr =ElFRncbjcorr->GetBinContent(bincorr);
		fwcorr=fwcorr/(1.-fwcorr);
		
		FillHist(tag+"_TL_diel_cj", looseelectrons[1].Pt() , wel*fw, 0., 100., 20);
		FillHist(tag+"_TL_diel_ccj_orr",el_pt_corr, wel*fwcorr, 0., 100., 20);
		
	      }
	      
	      
	      
	    }
	  }
	  
	  if(!PassID(looseelectrons[0], tightid)) {
	    if(PassID(looseelectrons[1], tightid))  {
	      
	      float el_pt_corr = looseelectrons.at(0).Pt()*(1+max(0.,(looseelectrons.at(0).PFRelIso(0.3)-tightiso))) ; 
	      float el_pt_corr1 = looseelectrons.at(1).Pt()*(1+max(0.,(looseelectrons.at(1).PFRelIso(0.3)-tightiso))) ; 
	      int bin= ElFR->FindBin(looseelectrons[0].Pt(), fabs(looseelectrons[0].Eta()));
	      int bincorr=ElFRcorr->FindBin(el_pt_corr, fabs(looseelectrons[0].Eta()));
	      float fw = ElFR->GetBinContent(bin);
	      float fwcorr =ElFRcorr->GetBinContent(bincorr);
	      
	      fw=fw/(1.-fw);
	      fwcorr=fwcorr/(1.-fwcorr);
	      
	      FillHist(tag+"_TL_diel", looseelectrons[1].Pt() , wel*fw, 0., 100., 20);
	      FillHist(tag+"_TL_diel_corr",el_pt_corr1, wel*fwcorr, 0., 100., 20);
	      
	      
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
		
		FillHist(tag+"_TL_diel_cj", looseelectrons[1].Pt() , wel*fw, 0., 100., 20);
		FillHist(tag+"_TL_diel_ccj_orr",el_pt_corr1, wel*fwcorr, 0., 100., 20);
		
	      }
	      else{
		fw=ElFRncbj->GetBinContent(bin);
		fw=fw/(1.-fw);
		fwcorr =ElFRncbjcorr->GetBinContent(bincorr);
		fwcorr=fwcorr/(1.-fwcorr);
		
		FillHist(tag+"_TL_diel_cj", looseelectrons[1].Pt() , wel*fw, 0., 100., 20);
		FillHist(tag+"_TL_diel_ccj_orr",el_pt_corr1, wel*fwcorr, 0., 100., 20);
		
	      }
	      
	    } 
	  }
	}
      }
    }
  }
  
  

  if(looseelectrons.size()==0) return;
  

  vector<snu::KElectron> fake_electrons;
  vector<snu::KElectron> fake_electrons_all;

  for(unsigned int iel=0; iel < looseelectrons.size(); iel++){
    if(looseelectrons[iel].GetType()==7) fake_electrons.push_back(looseelectrons[iel]);
    if(!looseelectrons[iel].MCMatched())fake_electrons_all.push_back(looseelectrons[iel]);
  }
  
  MakeMCFakes(fake_electrons, tag , jets, alljets, tightid, wel);
  MakeMCFakes(fake_electrons_all, tag+"all", jets, alljets, tightid, wel);
  MakeMCFakes(fake_electrons, tag+"noWeight" , jets, alljets, tightid, 1.);
  MakeMCFakes(fake_electrons_all, tag+"noWeightall", jets, alljets, tightid, 1.);

}

void FakeRateMC::MakeMCFakes(std::vector<snu::KElectron> fake_electrons, TString tag, std::vector<snu::KJet> jets,  std::vector<snu::KJet> alljets, TString tightid, float wel_fake){

  Float_t ptbins[11] = {5., 10., 15.,20.,25.,30.,35.,45.,60.,100., 200.};
  Float_t etabins2[5] = { 0.,0.8,  1.479, 2.,  2.5};

  if(fake_electrons.size() != 1) return;
  
  TString triggerslist_8="HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v";   /// -> tighter cut in lepton ID form tighter trigger emulation cut                        
  TString triggerslist_12="HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  TString triggerslist_18="HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  TString triggerslist_23="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
  TString triggerslist_33="HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"; ///       

  
  for(unsigned int iel=0; iel < fake_electrons.size(); iel++){
    
    bool closebjet = false;
    for(unsigned int ij =0 ; ij < alljets.size() ; ij++){
      if(fake_electrons.at(iel).DeltaR(alljets.at(ij)) < 0.5) {
	if(alljets.at(ij).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) closebjet = true;
      }
    }
    

    int jettype =  CloseJetTypeHAD(fake_electrons.at(iel), alljets);
    int jettypeHAD =  CloseJetType(fake_electrons.at(iel), alljets);
    
    FillHist(tag+"_ElectronType",fake_electrons[iel].GetType(),wel_fake, 0., 41., 41);
    FillHist(tag+"_ElectronMother",fake_electrons[iel].MotherPdgId(),wel_fake, 0., 2000., 2000);
    
    if(jettype==1)     FillHist(tag+"_ElectronType_b",fake_electrons[iel].GetType(),wel_fake, 0., 41., 41);
    if(jettype==2)     FillHist(tag+"_ElectronType_c",fake_electrons[iel].GetType(),wel_fake, 0., 41., 41);
    if(jettype==3)     FillHist(tag+"_ElectronType_l",fake_electrons[iel].GetType(),wel_fake, 0., 41., 41);
    
    if(jettype==1)    FillHist(tag+"_ElectronMother_b",fake_electrons[iel].MotherPdgId(),wel_fake, 0., 2000., 2000);
    if(jettype==2)    FillHist(tag+"_ElectronMother_c",fake_electrons[iel].MotherPdgId(),wel_fake, 0., 2000., 2000);
    if(jettype==3)    FillHist(tag+"_ElectronMother_l",fake_electrons[iel].MotherPdgId(),wel_fake, 0., 2000., 2000);

    


    if(PassID(fake_electrons[iel], tightid))     FillHist(tag+"_ElectronType_tight",fake_electrons[iel].GetType(),wel_fake, 0., 41., 41);
    float el_pt_corr = fake_electrons.at(iel).Pt()*(1+max(0.,(fake_electrons.at(iel).PFRelIso(0.3)- 0.08))) ; /// will need changing for systematics
    float el_pt =fake_electrons.at(iel).Pt();
    
    bool useevent=false;
    float METdphi = TVector2::Phi_mpi_pi(fake_electrons.at(0).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet));
    float MT = sqrt(2.* fake_electrons.at(0).Et()*eventbase->GetEvent().MET(snu::KEvent::pfmet) * (1 - cos( METdphi)));
    if(( (eventbase->GetEvent().MET(snu::KEvent::pfmet) < 20) && (MT < 25.)) ) {

      for(unsigned int ij=0; ij < jets.size(); ij++){
	if(jets.at(ij).Pt() < 40.) continue;
        float dphi =fabs(TVector2::Phi_mpi_pi(fake_electrons.at(iel).Phi()- jets.at(ij).Phi()));
        if(dphi > 2.5) useevent = true;
      }
    }
    

    
    if(useevent){
      FillHist(tag+"_JetElectronType",fake_electrons[iel].GetType(),wel_fake, 0., 41., 41);
      if(PassID(fake_electrons[iel], tightid))     FillHist(tag+"_JetElectronType_tight",fake_electrons[iel].GetType(),wel_fake, 0., 41., 41);
      
    }


    FillHist(tag+"_Loose_el_ptcorr_allmc",el_pt_corr, wel_fake, ptbins,10);
    if(jettype==1) FillHist(tag+"_Loose_el_ptcorr_b_allmc",el_pt_corr, wel_fake, ptbins,10);
    if(jettype==2) FillHist(tag+"_Loose_el_ptcorr_c_allmc",el_pt_corr, wel_fake, ptbins,10);
    if(jettype==3) FillHist(tag+"_Loose_el_ptcorr_l_allmc",el_pt_corr, wel_fake, ptbins,10);
    

    if(jettypeHAD==1) FillHist(tag+"_Loose_el_ptcorr_bh_allmc",el_pt_corr, wel_fake, ptbins,10);
    if(jettypeHAD==2) FillHist(tag+"_Loose_el_ptcorr_ch_allmc",el_pt_corr, wel_fake, ptbins,10);
    if(jettypeHAD==3) FillHist(tag+"_Loose_el_ptcorr_lh_allmc",el_pt_corr, wel_fake, ptbins,10);
    
    
    if(PassID(fake_electrons[iel], tightid))   {
      FillHist(tag+"_Tight_el_ptcorr_allmc",el_pt_corr, wel_fake, ptbins,10);
      if(jettype==1) FillHist(tag+"_Tight_el_ptcorr_b_allmc",el_pt_corr, wel_fake, ptbins,10);
      if(jettype==2) FillHist(tag+"_Tight_el_ptcorr_c_allmc",el_pt_corr, wel_fake, ptbins,10);
      if(jettype==3) FillHist(tag+"_Tight_el_ptcorr_l_allmc",el_pt_corr, wel_fake, ptbins,10);
      if(jettypeHAD==1) FillHist(tag+"_Tight_el_ptcorr_bh_allmc",el_pt_corr, wel_fake, ptbins,10);
      if(jettypeHAD==2) FillHist(tag+"_Tight_el_ptcorr_ch_allmc",el_pt_corr, wel_fake, ptbins,10);
      if(jettypeHAD==3) FillHist(tag+"_Tight_el_ptcorr_lh_allmc",el_pt_corr, wel_fake, ptbins,10);

    }
    

    FillHist(tag + "_fake_onlyfakes_jetflavour", jettype, wel_fake , -4., 4., 8);
    if(useevent)FillHist(tag + "_fake_onlyfakes_jetflavour_qcdCR", jettype, wel_fake , -4., 4., 8);
    

    FillHist(tag+"_Loose_el_eta_allmc",fake_electrons.at(iel).SCEta(), wel_fake, -3., 3., 60);
    FillHist(tag+"_Loose_el_pt_allmc",el_pt, wel_fake, ptbins,10);

    FillHist(tag+"_Loose_el_njet_allmc",jets.size(), wel_fake, 0., 6.,6);
    FillHist(tag+"_Loose_el_nbjet_allmc",NBJet(jets), wel_fake, 0., 6.,6);
    

    if(jettype==1) FillHist(tag+"_Loose_el_pt_b_allmc",el_pt, wel_fake, ptbins,10);
    if(jettype==2) FillHist(tag+"_Loose_el_pt_c_allmc",el_pt, wel_fake, ptbins,10);
    if(jettype==3) FillHist(tag+"_Loose_el_pt_l_allmc",el_pt, wel_fake, ptbins,10);

    if(PassID(fake_electrons[iel], tightid))   {
      FillHist(tag+"_Tight_el_eta_allmc",fake_electrons.at(iel).SCEta(), wel_fake, -3., 3., 60);
      FillHist(tag+"_Tight_el_pt_allmc",el_pt, wel_fake, ptbins,10);
      if(jettype==1) FillHist(tag+"_Tight_el_pt_b_allmc",el_pt, wel_fake, ptbins,10);
      if(jettype==2) FillHist(tag+"_Tight_el_pt_c_allmc",el_pt, wel_fake, ptbins,10);
      if(jettype==3) FillHist(tag+"_Tight_el_pt_l_allmc",el_pt_corr, wel_fake, ptbins,10);
      FillHist(tag+"_Tight_el_njet_allmc",jets.size(), wel_fake, 0., 6.,6);
      FillHist(tag+"_Tight_el_nbjet_allmc",NBJet(jets), wel_fake, 0., 6.,6);

    }
  

    if( (k_sample_name.Contains("QCD") ||  k_sample_name.Contains("qcd"))&& !useevent) return;
    
    if(fabs(fake_electrons.at(iel).SCEta()) < 0.8){
      FillHist(tag+"_Loose_el_eb1_pt1D",el_pt, wel_fake, ptbins,10);
      FillHist(tag+"_Loose_el_eb1_ptcorr1D",el_pt_corr, wel_fake, ptbins,10);
      if(jettype==1)     FillHist(tag+"_Loose_el_eb1_ptcorr1D_heavyb",el_pt_corr, wel_fake, ptbins,10);
      if(jettype==2)     FillHist(tag+"_Loose_el_eb1_ptcorr1D_heavyc",el_pt_corr, wel_fake, ptbins,10);
      if(jettype==3)     FillHist(tag+"_Loose_el_eb1_ptcorr1D_light",el_pt_corr, wel_fake, ptbins,10);
      
      FillHist(tag+"_Loose_el_eta1_allmc",fake_electrons.at(iel).SCEta(), wel_fake, -3., 3., 60);

      
      for(unsigned int ibinX = 0 ; ibinX < 13; ibinX++){
	for(unsigned int ibinY = 0 ; ibinY < 20; ibinY++){
	  float isocut = 0. + 0.05*float(ibinX)+0.0001;
	  float mvacut = -1. + 0.1*float(ibinY)+0.0001;
	  std::vector<snu::KElectron> tmpel;
	  if(fake_electrons.at(iel).PFRelIso(0.3) < isocut){
	    if(fake_electrons.at(iel).MVA() > mvacut){
	      tmpel.push_back(fake_electrons.at(iel));
	    }
	  }
	  if(tmpel.size()==1){
	    if(jettype==1)     FillHist(tag+"_opt_loosemva_Loose_el_eb1_ptcorr1D_heavyb", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20); 
	    if(jettype==2)     FillHist(tag+"_opt_loosemva_Loose_el_eb1_ptcorr1D_heavyc", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20); 
	    if(jettype==3)     FillHist(tag+"_opt_loosemva_Loose_el_eb1_ptcorr1D_light", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20); 
	    if(PassID(tmpel[0], tightid))   {
	      if(jettype==1)     FillHist(tag+"_opt_loosemva_Tight_el_eb1_ptcorr1D_heavyb", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20);
	      if(jettype==2)     FillHist(tag+"_opt_loosemva_Tight_el_eb1_ptcorr1D_heavyc", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20);
	      if(jettype==3)     FillHist(tag+"_opt_loosemva_Tight_el_eb1_ptcorr1D_light", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20);
	    }
	  }
	}
      }

      for(unsigned int ibinX = 0 ; ibinX < 13; ibinX++){
        for(unsigned int ibinY = 0 ; ibinY < 100; ibinY++){
          float isocut = 0. + 0.05*float(ibinX)+0.0001;
          float mvacut = -1. + 0.02*float(ibinY)+0.0001;
	  std::vector<snu::KElectron> tmpel;
          if(fake_electrons.at(iel).PFRelIso(0.3) < isocut){
            if(fake_electrons.at(iel).MVA() > mvacut){
              tmpel.push_back(fake_electrons.at(iel));
            }
          }
          if(tmpel.size()==1){

            if(jettype==1)     FillHist(tag+"_opt_finer_loosemva_Loose_el_eb1_ptcorr1D_heavyb", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
            if(jettype==2)     FillHist(tag+"_opt_finer_loosemva_Loose_el_eb1_ptcorr1D_heavyc", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
            if(jettype==3)     FillHist(tag+"_opt_finer_loosemva_Loose_el_eb1_ptcorr1D_light", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
            if(PassID(tmpel[0], tightid))   {
              if(jettype==1)     FillHist(tag+"_opt_finer_loosemva_Tight_el_eb1_ptcorr1D_heavyb", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
              if(jettype==2)     FillHist(tag+"_opt_finer_loosemva_Tight_el_eb1_ptcorr1D_heavyc", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
              if(jettype==3)     FillHist(tag+"_opt_finer_loosemva_Tight_el_eb1_ptcorr1D_light", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
            }
          }
        }
      }

      
      if(closebjet){
        FillHist(tag+"_Loose_el_cj_eb1_pt1D",el_pt, wel_fake, ptbins,10);
      }
      else{
        FillHist(tag+"_Loose_el_ncj_eb1_pt1D",el_pt, wel_fake, ptbins,10);
      }
      if(PassID(fake_electrons.at(iel), tightid))   {
	FillHist(tag+"_Tight_el_eta1_allmc",fake_electrons.at(iel).SCEta(), wel_fake, -3., 3., 60);
        FillHist(tag+"_Tight_el_eb1_pt1D",el_pt, wel_fake, ptbins,10);
        FillHist(tag+"_Tight_el_eb1_ptcorr1D",el_pt_corr, wel_fake, ptbins,10);


	if(jettype==1)     FillHist(tag+"_Tight_el_eb1_ptcorr1D_heavyb",el_pt_corr, wel_fake, ptbins,10);
	if(jettype==2)     FillHist(tag+"_Tight_el_eb1_ptcorr1D_heavyc",el_pt_corr, wel_fake, ptbins,10);
	if(jettype==3)     FillHist(tag+"_Tight_el_eb1_ptcorr1D_light",el_pt_corr, wel_fake, ptbins,10);

        if(closebjet){
          FillHist(tag+"_Tight_el_cj_eb1_pt1D",el_pt, wel_fake, ptbins,10);
        }
        else{
          FillHist(tag+"_Tight_el_ncj_eb1_pt1D",el_pt, wel_fake, ptbins,10);
        }

      }
    }
    else  if(fabs(fake_electrons.at(iel).Eta()) < 1.5){
      FillHist(tag+"_Loose_el_eb2_pt1D",el_pt, wel_fake, ptbins,10);
      FillHist(tag+"_Loose_el_eb2_ptcorr1D",el_pt_corr, wel_fake, ptbins,10);

      if(jettype==1)     FillHist(tag+"_Loose_el_eb2_ptcorr1D_heavyb",el_pt_corr, wel_fake, ptbins,10);
      if(jettype==2)     FillHist(tag+"_Loose_el_eb2_ptcorr1D_heavyc",el_pt_corr, wel_fake, ptbins,10);
      if(jettype==3)     FillHist(tag+"_Loose_el_eb2_ptcorr1D_light",el_pt_corr, wel_fake, ptbins,10);

      FillHist(tag+"_Loose_el_eta2_allmc",fake_electrons.at(iel).SCEta(), wel_fake, -3., 3., 60);

      for(unsigned int ibinX = 0 ; ibinX < 13; ibinX++ ){
        for(unsigned int ibinY = 0 ; ibinY < 20;ibinY++){
          float isocut = 0. + 0.05*float(ibinX)+0.0001;;
          float mvacut = -1. + 0.1*float(ibinY)+0.0001;;

	  std::vector<snu::KElectron> tmpel;
          if(fake_electrons.at(iel).PFRelIso(0.3) < isocut){
            if(fake_electrons.at(iel).MVA() > mvacut){
	      tmpel.push_back(fake_electrons.at(iel));
            }
          }
          if(tmpel.size()==1){
            if(jettype==1)     FillHist(tag+"_opt_loosemva_Loose_el_eb2_ptcorr1D_heavyb", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20);
            if(jettype==2)     FillHist(tag+"_opt_loosemva_Loose_el_eb2_ptcorr1D_heavyc", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20);
            if(jettype==3)     FillHist(tag+"_opt_loosemva_Loose_el_eb2_ptcorr1D_light", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20);
            if(PassID(tmpel[0], tightid))   {
	      if(jettype==1)     FillHist(tag+"_opt_loosemva_Tight_el_eb2_ptcorr1D_heavyb", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20);
              if(jettype==2)     FillHist(tag+"_opt_loosemva_Tight_el_eb2_ptcorr1D_heavyc", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20);
              if(jettype==3)     FillHist(tag+"_opt_loosemva_Tight_el_eb2_ptcorr1D_light", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20);
            }
          }

        }
      } 
      
      for(unsigned int ibinX = 0 ; ibinX < 13; ibinX++ ){
        for(unsigned int ibinY = 0 ; ibinY < 100;ibinY++){
          float isocut = 0. + 0.05*float(ibinX)+0.0001;;
          float mvacut = -1. + 0.02*float(ibinY)+0.0001;;
	  std::vector<snu::KElectron> tmpel;
          if(fake_electrons.at(iel).PFRelIso(0.3) < isocut){
            if(fake_electrons.at(iel).MVA() > mvacut){
              tmpel.push_back(fake_electrons.at(iel));
            }
          }
          if(tmpel.size()==1){
            if(jettype==1)     FillHist(tag+"_opt_finer_loosemva_Loose_el_eb2_ptcorr1D_heavyb", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
            if(jettype==2)     FillHist(tag+"_opt_finer_loosemva_Loose_el_eb2_ptcorr1D_heavyc", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
            if(jettype==3)     FillHist(tag+"_opt_finer_loosemva_Loose_el_eb2_ptcorr1D_light", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
            if(PassID(fake_electrons[iel], tightid))   {
              if(jettype==1)     FillHist(tag+"_opt_finer_loosemva_Tight_el_eb2_ptcorr1D_heavyb", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
              if(jettype==2)     FillHist(tag+"_opt_finer_loosemva_Tight_el_eb2_ptcorr1D_heavyc", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
              if(jettype==3)     FillHist(tag+"_opt_finer_loosemva_Tight_el_eb2_ptcorr1D_light", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
            }
          }

        }
      }


      if(closebjet){
        FillHist(tag+"_Loose_el_cj_eb2_pt1D",el_pt, wel_fake, ptbins,10);
      }
      else{
        FillHist(tag+"_Loose_el_ncj_eb2_pt1D",el_pt, wel_fake, ptbins,10);
      }

      if(PassID(fake_electrons[iel], tightid))   {
	FillHist(tag+"_Tight_el_eta2_allmc",fake_electrons.at(iel).SCEta(), wel_fake, -3., 3., 60);

        FillHist(tag+"_Tight_el_eb2_pt1D",el_pt, wel_fake, ptbins,10);
        FillHist(tag+"_Tight_el_eb2_ptcorr1D",el_pt_corr, wel_fake, ptbins,10);

	if(jettype==1)     FillHist(tag+"_Tight_el_eb2_ptcorr1D_heavyb",el_pt_corr, wel_fake, ptbins,10);
        if(jettype==2)     FillHist(tag+"_Tight_el_eb2_ptcorr1D_heavyc",el_pt_corr, wel_fake, ptbins,10);
        if(jettype==3)     FillHist(tag+"_Tight_el_eb2_ptcorr1D_light",el_pt_corr, wel_fake, ptbins,10);



        if(closebjet){
          FillHist(tag+"_Tight_el_cj_eb2_pt1D",el_pt, wel_fake, ptbins,10);
        }
	else{
          FillHist(tag+"_Tight_el_ncj_eb2_pt1D",el_pt, wel_fake, ptbins,10);
        }
      }
    }
    

    else{
      /// EE REGION 

      FillHist(tag+"_Loose_el_ee_pt1D",el_pt, wel_fake, ptbins,10);
      FillHist(tag+"_Loose_el_ee_ptcorr1D",el_pt_corr, wel_fake, ptbins,10);


      if(jettype==1)     FillHist(tag+"_Loose_el_ee_ptcorr1D_heavyb",el_pt_corr, wel_fake, ptbins,10);
      if(jettype==2)     FillHist(tag+"_Loose_el_ee_ptcorr1D_heavyc",el_pt_corr, wel_fake, ptbins,10);
      if(jettype==3)     FillHist(tag+"_Loose_el_ee_ptcorr1D_light",el_pt_corr, wel_fake, ptbins,10);
      if(closebjet){
        FillHist(tag+"_Loose_el_cj_ee_pt1D",el_pt, wel_fake, ptbins,10);
      }
      else{
        FillHist(tag+"_Loose_el_ncj_ee_pt1D",el_pt, wel_fake, ptbins,10);
      }

      if(PassID(fake_electrons[iel], tightid))   {
        FillHist(tag+"_Tight_el_ee_pt1D",el_pt, wel_fake, ptbins,10);
        FillHist(tag+"_Tight_el_ee_ptcorr1D",el_pt_corr, wel_fake, ptbins,10);

	for(unsigned int ibinX = 0 ; ibinX < 13;ibinX++){
          for(unsigned int ibinY = 0 ; ibinY < 20;ibinY++){
            float isocut = 0. + 0.05*float(ibinX)+0.0001;;
            float mvacut = -1. + 0.1*float(ibinY)+0.0001;;
	    std::vector<snu::KElectron> tmpel;
            if(fake_electrons.at(iel).PFRelIso(0.3) < isocut){
              if(fake_electrons.at(iel).MVA() > mvacut){
                tmpel.push_back(fake_electrons.at(iel));
              }
            }
            if(tmpel.size()==1){
              if(jettype==1)     FillHist(tag+"_opt_loosemva_Loose_el_ee_ptcorr1D_heavyb", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20);
              if(jettype==2)     FillHist(tag+"_opt_loosemva_Loose_el_ee_ptcorr1D_heavyc", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20);
              if(jettype==3)     FillHist(tag+"_opt_loosemva_Loose_el_ee_ptcorr1D_light", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20);
              if(PassID(fake_electrons[iel], tightid))   {
                if(jettype==1)     FillHist(tag+"_opt_loosemva_Tight_el_ee_ptcorr1D_heavyb", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20);
                if(jettype==2)     FillHist(tag+"_opt_loosemva_Tight_el_ee_ptcorr1D_heavyc", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20);
                if(jettype==3)     FillHist(tag+"_opt_loosemva_Tight_el_ee_ptcorr1D_light", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 20);
              }
            }
          }
        }

	for(unsigned int ibinX = 0 ; ibinX < 13;ibinX++){
          for(unsigned int ibinY = 0 ; ibinY < 100;ibinY++){
            float isocut = 0. + 0.05*float(ibinX)+0.0001;;
            float mvacut = -1. + 0.02*float(ibinY)+0.0001;;
	    std::vector<snu::KElectron> tmpel;
            if(fake_electrons.at(iel).PFRelIso(0.3) < isocut){
              if(fake_electrons.at(iel).MVA() > mvacut){
                tmpel.push_back(fake_electrons.at(iel));
              }
            }
            if(tmpel.size()==1){
              if(jettype==1)     FillHist(tag+"_opt_finer_loosemva_Loose_el_ee_ptcorr1D_heavyb", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
              if(jettype==2)     FillHist(tag+"_opt_finer_loosemva_Loose_el_ee_ptcorr1D_heavyc", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
              if(jettype==3)     FillHist(tag+"_opt_finer_loosemva_Loose_el_ee_ptcorr1D_light", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
              if(PassID(fake_electrons[iel], tightid))   {
                if(jettype==1)     FillHist(tag+"_opt_finer_loosemva_Tight_el_ee_ptcorr1D_heavyb", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
                if(jettype==2)     FillHist(tag+"_opt_finer_loosemva_Tight_el_ee_ptcorr1D_heavyc", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
                if(jettype==3)     FillHist(tag+"_opt_finer_loosemva_Tight_el_ee_ptcorr1D_light", isocut, mvacut,  wel_fake, 0., 0.65, 13, -1., 1., 100);
              }
            }
          }
        }

	

	if(jettype==1)     FillHist(tag+"_Tight_el_ee_ptcorr1D_heavyb",el_pt_corr, wel_fake, ptbins,10);
        if(jettype==2)     FillHist(tag+"_Tight_el_ee_ptcorr1D_heavyc",el_pt_corr, wel_fake, ptbins,10);
        if(jettype==3)     FillHist(tag+"_Tight_el_ee_ptcorr1D_light",el_pt_corr, wel_fake, ptbins,10);
        if(closebjet){
          FillHist(tag+"_Tight_el_cj_ee_pt1D",el_pt, wel_fake, ptbins,10);
        }
        else{
          FillHist(tag+"_Tight_el_ncj_ee_pt1D",el_pt, wel_fake, ptbins,10);
        }
      }
    }

    FillHist(tag+"_Loose_el_pt",el_pt,fabs(fake_electrons.at(iel).Eta()), wel_fake, ptbins,10 , etabins2, 4);
    FillHist(tag+"_Loose_el_ptcorr",el_pt_corr,fabs(fake_electrons.at(iel).Eta()), wel_fake, ptbins,10 , etabins2, 4);
    if(closebjet) {
      FillHist(tag+"_Loose_el_cj_pt",el_pt,fabs(fake_electrons.at(iel).Eta()), wel_fake, ptbins,10 , etabins2, 4);
      FillHist(tag+"_Loose_el_cj_ptcorr",el_pt_corr,fabs(fake_electrons.at(iel).Eta()), wel_fake, ptbins,10 , etabins2, 4);
    }
    else{
      FillHist(tag+"_Loose_el_ncj_pt",el_pt,fabs(fake_electrons.at(iel).Eta()), wel_fake, ptbins,10 , etabins2, 4);
      FillHist(tag+"_Loose_el_ncj_ptcorr",el_pt_corr,fabs(fake_electrons.at(iel).Eta()), wel_fake, ptbins,10 , etabins2, 4);
      
    }
    if(PassID(fake_electrons[iel], tightid))   {
      FillHist(tag+"_Tight_el_pt",el_pt,fabs(fake_electrons.at(iel).Eta()), wel_fake, ptbins,10 , etabins2, 4);
      FillHist(tag+"_Tight_el_ptcorr",el_pt_corr,fabs(fake_electrons.at(iel).Eta()), wel_fake, ptbins,10 , etabins2, 4);
      if(closebjet) {
	FillHist(tag+"_Tight_el_cj_pt",el_pt,fabs(fake_electrons.at(iel).Eta()), wel_fake, ptbins,10 , etabins2, 4);
	FillHist(tag+"_Tight_el_cj_ptcorr",el_pt_corr,fabs(fake_electrons.at(iel).Eta()), wel_fake, ptbins,10 , etabins2, 4);
      }
      else{
	FillHist(tag+"_Tight_el_ncj_pt",el_pt,fabs(fake_electrons.at(iel).Eta()), wel_fake, ptbins,10 , etabins2, 4);
	FillHist(tag+"_Tight_el_ncj_ptcorr",el_pt_corr,fabs(fake_electrons.at(iel).Eta()), wel_fake, ptbins,10 , etabins2, 4);
	
      }
    }
    
  }
  
  FillHist(tag+"_JetElectronType",fake_electrons[0].GetType(),wel_fake, 0., 41., 41);
  if(PassID(fake_electrons[0], tightid))     FillHist(tag+"_JetElectronType_tight",fake_electrons[0].GetType(),wel_fake, 0., 41., 41);


}




int FakeRateMC::CloseJetType(snu::KMuon mu, std::vector<snu::KJet> jets){

  //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools#Jet_flavour_in_PAT                                                                                                                                                                                        
  int jetFlavour=-999999;

  snu::KJet closejet;
  float mindR=0.4;
  bool jetfound=false;
  for(unsigned int ij=0; ij < jets.size(); ij++){
    if(jets.at(ij).DeltaR(mu) < mindR) {
      closejet=jets.at(ij);
      mindR=jets.at(ij).DeltaR(mu) ;
      jetfound=true;
    }
  }

  if(jetfound){

    jetFlavour = fabs(int(closejet.HadronFlavour()));

    if(jetFlavour == 5 ) return 1;
    if(jetFlavour == 4 ) return 2;
    if(jetFlavour==1 || jetFlavour==2 || jetFlavour==3 || jetFlavour==21)  return 3;
    if(jetFlavour == 0) return 4;
  }
  return -1.;
}
int FakeRateMC::CloseJetType(snu::KElectron el, std::vector<snu::KJet> jets){
  
  //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools#Jet_flavour_in_PAT
  int jetFlavour=-999999;
  
  snu::KJet closejet;
  float mindR=0.4;
  bool jetfound=false;
  for(unsigned int ij=0; ij < jets.size(); ij++){
    if(jets.at(ij).DeltaR(el) < mindR) {
      closejet=jets.at(ij);
      mindR=jets.at(ij).DeltaR(el) ;
      jetfound=true;
    }
  }
  
  if(jetfound){
    
    jetFlavour = fabs(int(closejet.PartonFlavour()));
    
    if(jetFlavour == 5 ) return 1;
    if(jetFlavour == 4 ) return 2;
    if(jetFlavour==1 || jetFlavour==2 || jetFlavour==3 || jetFlavour==21)  return 3;
    if(jetFlavour == 0) return 4;
  }
  return -1.;
}

int FakeRateMC::CloseJetTypeHAD(snu::KElectron el, std::vector<snu::KJet> jets){

  int jetFlavour=-999999;

  snu::KJet closejet;
  float mindR=0.4;
  bool jetfound=false;
  for(unsigned int ij=0; ij < jets.size(); ij++){
    if(jets.at(ij).DeltaR(el) < mindR) {
      closejet=jets.at(ij);
      mindR=jets.at(ij).DeltaR(el) ;
      jetfound=true;
    }
  }

  if(jetfound){

    jetFlavour = fabs(int(closejet.HadronFlavour()));

    if(jetFlavour == 5 ) return 1;
    if(jetFlavour == 4 ) return 2;
    if(jetFlavour == 0 ) return 3;

  }
  return -1.;

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


