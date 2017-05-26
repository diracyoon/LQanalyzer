/***************************************************************************
 * @Project: LQFakeRateCalculator_Mu Frame - ROOT-based analysis framework for Korea SNU
OB * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/


/// Local includes
#include "FakeRateCalculator_Mu.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (FakeRateCalculator_Mu);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
FakeRateCalculator_Mu::FakeRateCalculator_Mu() :  AnalyzerCore(),  out_electrons(0) {


  // To have the correct name in the log:                                                                                                                            
  SetLogName("FakeRateCalculator_Mu");

  Message("In FakeRateCalculator_Mu constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void FakeRateCalculator_Mu::InitialiseAnalysis() throw( LQError ) {
  
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
  MakeCleverHistograms(sighist_ee, "SingleElectron_prompt_unprescaled");
  MakeCleverHistograms(sighist_ee, "SingleElectron_prompt_unprescaled2");
  return;
}


void FakeRateCalculator_Mu::ExecuteEvents()throw( LQError ){
    
  //// Initial event cuts
  /// MET FIleters 
  if(!PassMETFilter()) return;     
  
  /// Require good promary vertex 
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex  
  numberVertices = eventbase->GetEvent().nVertices();   

  if (eventbase->GetEvent().MET(snu::KEvent::pfmet) > 20) return;

  

  /// These run on double electron dataset
  if(isData || !isData){

    // GetFakeRateAndPromptRates fills fake rate + prompt rate plots
    // hist names are tagged with second argu
    // dijet method

    
    /// ELECTRON_HN_FAKELOOSEST has no ID cuts
    /// only pt/eta/chargeconst/looseIP/

    std::vector<snu::KMuon> tmploose_mu = GetMuons("MUON_HN_FAKELOOSEST",false);
    
    if(tmploose_mu.size() != 1) return;
    vector<float> vcut_dxy_b;
    vector<TString> vcut_dxy_b_s;
    
    for(unsigned int dxy_b=0;dxy_b < 10; dxy_b++){
      float cut_dxy_b =  float(dxy_b)*0.005 + 0.01;
      vcut_dxy_b.push_back(cut_dxy_b);
      stringstream ss;
      ss <<cut_dxy_b;
      vcut_dxy_b_s.push_back(TString(ss.str()));
    }
    vector<float> vcut_dxysig_b;
    vcut_dxysig_b.push_back(3.);
    vcut_dxysig_b.push_back(4.);
    vector<TString> vcut_dxysig_b_s;
    vcut_dxysig_b_s.push_back("3");
    vcut_dxysig_b_s.push_back("4");
    


    vector<float> vcut_dz_b;
    vector<TString> vcut_dz_b_s;
    
    for(unsigned int dz_b=0;dz_b < 4; dz_b++){
      float cut_dz_b =  float(dz_b)*0.02 + 0.04;
      vcut_dz_b.push_back(cut_dz_b);
      stringstream ss;
      ss <<cut_dz_b;
      vcut_dz_b_s.push_back(TString(ss.str()));
    }   

    

    vector<float> vcut_iso_b;
    vector<TString> vcut_iso_b_s;
    for(unsigned int iso_b=0;iso_b < 5; iso_b++){
      float cut_iso_b = float(iso_b)*0.01 + 0.05;
      vcut_iso_b.push_back(cut_iso_b);
      stringstream ss;
      ss <<cut_iso_b;
      vcut_iso_b_s.push_back(TString(ss.str()));
    }
    
    //// Loop over cuts and fill loose and tight el and get fake rates for ID
    for(unsigned int dxy_b=0; dxy_b < vcut_dxy_b.size(); dxy_b++){
      for(unsigned int dxysig_b=0; dxysig_b < vcut_dxysig_b.size(); dxysig_b++){
	for(unsigned int dz_b=0; dz_b < vcut_dz_b.size(); dz_b++){
	  for(unsigned int iso_b=0; iso_b < vcut_iso_b.size(); iso_b++){
	    
	    std::vector<snu::KMuon> loose_mu;
	    std::vector<snu::KMuon> tight_mu;
	    for(unsigned int imu=0; imu<tmploose_mu.size(); imu++){
	      float reliso = tmploose_mu[imu].RelIso04();
	      if(fabs(tmploose_mu[imu].dZ()) > vcut_dz_b[dz_b]) continue;
	      
	      loose_mu.push_back(tmploose_mu[imu]);
	      if(fabs(tmploose_mu[imu].dXYSig()) > vcut_dxysig_b[dxysig_b]) continue;
	      if(tmploose_mu[imu].GlobalChi2()  > 10.) continue;
	      if(fabs(tmploose_mu[imu].dXY()) > vcut_dxy_b[dxy_b]) continue;
	      if(reliso > vcut_iso_b[iso_b]) continue;
	      tight_mu.push_back(tmploose_mu[imu]);
	      
	      GetFakeRateAndPromptRates(loose_mu,"dijet_"+vcut_iso_b_s[iso_b]+"_"+vcut_dxy_b_s[dxy_b]+"_"+vcut_dxysig_b_s[dxysig_b]+"_"+vcut_dz_b_s[dz_b],tight_mu,weight, vcut_iso_b[iso_b],true,  false);
	    }
	  }
	}//iso
      }//dz
    }//dxy
  }
  return;
  
  
  if((isData&&k_channel == "SingleElectron") || !isData){

  }
}

void FakeRateCalculator_Mu::MakeSingleMuonCRPlots(TString looseid, TString eltag, TString tightid, float w, bool usepujetid){

}


void FakeRateCalculator_Mu::GetFakeRateAndPromptRates(std::vector<snu::KMuon> muonLooseColl, TString mutag, std::vector<snu::KMuon> muonTightColl, float w, float isocut, bool usepujetid, bool runall){

  // PileUpWeight is period dependant 
  // MC events are split into 7 data periods and 
  if (!k_isdata) {
    w = w * MCweight * eventbase->GetEvent().PeriodPileUpWeight(GetMCPeriodRandom());
  }

  TString triggerslist_3="HLT_Mu3_PFJet40_v";
  TString triggerslist_8="HLT_Mu8_TrkIsoVVL_v";
  TString triggerslist_17="HLT_Mu17_TrkIsoVVL_v";

  float id_weight=1.;
  float reco_weight=1.;

  /// USE HN jets. Add pileup ID? Currently 
  std::vector<snu::KJet> jetCollTight = GetJets("JET_HN");
  std::vector<snu::KJet> jetColl           = GetJets("JET_NOLEPTONVETO");
  if(usepujetid){
    jetCollTight = GetJets("JET_HN_PU");
    jetColl      = GetJets("JET_NOLEPTONVETO_PU");
  }

  /// remove events with no loose leptons
  if(muonLooseColl.size()<1) return;
  std::vector<snu::KElectron> elColl = GetElectrons("ELECTRON_HN_VETO");  // loose selection
  if(elColl.size() > 0) return;
  
   /// Get prescale for single el event. Returns 1. or 0. for data
  float prescale_trigger =  GetPrescale(muonLooseColl,   PassTrigger(triggerslist_3),PassTrigger(triggerslist_8), PassTrigger(triggerslist_17),TargetLumi);
  
  /// Make standard plots for loose and tight collection dijet                                                                                                                   
  if(muonLooseColl.size()==1)MakeFakeRatePlots("", mutag, muonTightColl ,muonLooseColl,  jetCollTight, jetColl,  prescale_trigger, isocut, w, true);
  
  if(!runall) return;


  return;
}// End of execute event loop


float FakeRateCalculator_Mu::GetPrescale( std::vector<snu::KMuon> muons,bool pass3, bool pass2, bool pass1, float fake_total_lum ){
  
  float prescale_trigger= 1.;
  if(muons.size() != 1) return 0.;
  if(muons.size() ==1){

    /// 10 - 20  HLT_Mu7
    /// 20 - INF  HLT_Mu17
    
    if(muons.at(0).Pt() >= 20.){
      if((isData&&k_channel != "DoubleMuon")) return 0.;

      if(pass1){
	if(isData) return 1.;
	prescale_trigger = WeightByTrigger("HLT_Mu17_TrkIsoVVL_v", fake_total_lum); //// 20 + GeV bins
      }
      else {
	if(isData) return 0;
	return 0;
      }
    }
    else  if(muons.at(0).Pt() >= 10.){
      if((isData&&k_channel != "DoubleMuon")) return 0.;

      if(pass2){
	if(isData) return 1.;
        prescale_trigger =  WeightByTrigger("HLT_Mu8_TrkIsoVVL_v", fake_total_lum) ; //// 20 + GeV bins
      }
      else {
	if(isData) return 0;
	return 0;
      }
    }
    else  if(muons.at(0).Pt() >= 5.){
      if((isData&&k_channel != "SingleMuon")) return 0.;
      if(pass3){
        if(isData) return 1.;
        prescale_trigger =  WeightByTrigger("HLT_Mu3_PFJet40_v", fake_total_lum) ; //// 20 + GeV bins                                                                         
      }
      else {
        if(isData) return 0;
	return 0;
      }
    }

  }
  if(prescale_trigger == 0.) return 0.;
  if(k_isdata) return 1.;
  if(muons.at(0).Pt() < 10.) return 0.;
  
  
  return prescale_trigger;
}





void FakeRateCalculator_Mu::MakeFakeRatePlots(TString label, TString mutag,   std::vector<snu::KMuon> muons_tight, std::vector<snu::KMuon> muons,  std::vector<snu::KJet> jets, std::vector<snu::KJet> alljets, float prescale_w, float w, float isocut,bool makebasicplots){
  
  
  if(muons.size() != 1 ) return;
  
  if(prescale_w==0.) return;
  
  bool useevent20 = UseEvent(muons , jets, 20., prescale_w, w); 
  bool useevent30 = UseEvent(muons , jets, 30., prescale_w, w); 
  bool useevent40 = UseEvent(muons , jets, 40., prescale_w, w); 
  bool useevent60 = UseEvent(muons , jets, 60., prescale_w, w); 

  bool truth_match= false;
  /// FOR FAKE RATE SUBTRACTION NEED ONLY PROMPT ELECTRONS
  if(!k_isdata) {
    if(muons.size() > 0){
      /// ONLY INCLUDE ELECTRONS FROM W/Z
      if((muons.at(0).MCMatched()))  truth_match=true;
    }
    else truth_match = true;
  }
  else truth_match=true;
  
  label= mutag;

  if(truth_match){
    if(jets.size() >= 1){
      if(makebasicplots){
	if(useevent40)GetFakeRates(muons, muons_tight,label, jets, alljets,  label+"_40", (prescale_w * w),isocut,makebasicplots);
      }
      else{
	//if(useevent20)GetFakeRates(muons, muons_tight,label, jets, alljets,  label+"_20", (prescale_w * w),makebasicplots);
	//if(useevent30)GetFakeRates(muons, muons_tight,label, jets, alljets,  label+"_30", (prescale_w * w),makebasicplots);
	if(useevent40)GetFakeRates(muons, muons_tight,label, jets, alljets,  label+"_40", (prescale_w * w),isocut,makebasicplots);
	//if(useevent60)GetFakeRates(muons, muons_tight, label,jets, alljets,  label+"_60", (prescale_w * w),makebasicplots);
      }
    }
  }
}

bool FakeRateCalculator_Mu::UseEvent(std::vector<snu::KMuon> muons,  std::vector< snu::KJet> jets, float awayjetcut, float precale_weight, float wt){
  
  bool useevent = false;
  if(muons.size() != 1) return false;
  if(precale_weight != 0.){
    /// Z and W veto
    Double_t MT=0;
    Double_t METdphi=0;
    for(unsigned int w = 0; w < muons.size();  w++){
      METdphi = TVector2::Phi_mpi_pi(muons.at(w).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet));
      MT = sqrt(2.* muons.at(w).Et()*eventbase->GetEvent().MET(snu::KEvent::pfmet) * (1 - cos( METdphi)));
      if(( (eventbase->GetEvent().MET(snu::KEvent::pfmet) < 20) && (MT < 25.)) ) {
	
        for (unsigned int ielT=0; ielT < muons.size(); ielT++){
          for(unsigned int ij=0; ij < jets.size(); ij++){
            if(jets.at(ij).Pt() < awayjetcut) continue;
            float dphi =fabs(TVector2::Phi_mpi_pi(muons.at(ielT).Phi()- jets.at(ij).Phi()));
	    
	    if(dphi > 2.5){
	      if((jets.at(ij).Pt() /  muons.at(ielT).Pt() ) < 1.) continue;
	      useevent = true;
	    }
          }
        }
      }
    }
  }
  return useevent;
}

void FakeRateCalculator_Mu::GetFakeRates(std::vector<snu::KMuon> loose_mu, std::vector<snu::KMuon> tight_mu, TString tightlabel,  std::vector<snu::KJet> jets,  std::vector<snu::KJet> alljets, TString tag, double w, float isocut, bool basicplots){
  
   
  Float_t ptbins[11] = { 5.,10., 15.,20.,25.,30.,35.,45.,60.,100., 200.};
  Float_t ptbinsb[9] = { 5., 10., 15.,20.,30.,45.,60.,100., 200.};
  Float_t etabin[2] = { 0.,  2.5};
  Float_t etabins[4] = { 0., 0.8,1.479,  2.5};
  Float_t etabins2[5] = { 0.,0.8,  1.479, 2.,  2.5};

  /// for most cuts just plot pt_eta
  if(basicplots){
    
    if(loose_mu.size() == 1 && jets.size() >= 1){
      float el_pt = loose_mu.at(0).Pt();
      float el_pt_corr = loose_mu.at(0).Pt()*(1+max(0.,(loose_mu.at(0).RelIso04()-isocut))) ; /// will need changing for systematics
      
      FillHist(("LooseMu" + tag + "_pt_eta").Data(), el_pt, fabs(loose_mu.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
      FillHist(("LooseMu" + tag + "_ptcorr_eta").Data(), el_pt_corr, fabs(loose_mu.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
      

      if( tight_mu.size() == 1){
	FillHist(("TightMu" + tag + "_pt_eta").Data(), el_pt, fabs(tight_mu.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
	FillHist(("TightMu" + tag + "_ptcorr_eta").Data(), el_pt_corr, fabs(tight_mu.at(0).Eta()),  w, ptbins, 10 , etabins2, 4);
	
      }
    }

    return;
  }
  
  int nbjet(0);
  float tmp_deltaR=1000.;
  for(unsigned int ij =0 ; ij < jets.size() ; ij++){
    if(jets.at(ij).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) nbjet++;
    if(loose_mu.at(0).DeltaR(jets.at(ij)) < tmp_deltaR) tmp_deltaR = loose_mu.at(0).DeltaR(jets.at(ij));
  }
 
  float awayjetpt = 0.;
  for(unsigned int ij =0 ; ij < alljets.size() ; ij++){
    float dphi =fabs(TVector2::Phi_mpi_pi(loose_mu.at(0).Phi()- alljets.at(ij).Phi()));
    if(dphi > 2.5) awayjetpt = alljets.at(ij).Pt();
  }

  if(loose_mu.size() == 1 && tight_mu.size() == 1 && jets.size() >= 1){
   
    float el_pt = tight_mu.at(0).Pt();
   
    // 1D FakeRates
    FillHist(("TightMu" + tag + "_eta").Data(), tight_mu.at(0).Eta(), w, -2.5, 2.5,50);
    FillHist(("TightMu" + tag + "_pt").Data(),el_pt, w,  ptbins, 9);
   
    if(fabs(tight_mu.at(0).Eta()) < 1.5)     FillHist(("TightMu" + tag + "_pt_barrel").Data(),el_pt, w,  ptbins, 9);
    else FillHist(("TightMu" + tag + "_pt_endcap").Data(),el_pt, w,  ptbins, 9);

    FillHist(("TightMu" + tag + "_njets").Data(), jets.size(), w, 0.,5.,5);
    FillHist(("TightMu" + tag + "_nvertices").Data(), eventbase->GetEvent().nVertices(), w, 0., 30., 30);
    FillHist(("TightMu" + tag + "_nbjet").Data(), nbjet, w, 0., 4.,4); 
    FillHist(("TightMu" + tag + "_pt_eta").Data(), el_pt, fabs(tight_mu.at(0).Eta()),  w, ptbins, 9 , etabins2, 4);
   
    if(nbjet > 0){
      FillHist(("TightMu" + tag + "_bjet_eta").Data(), tight_mu.at(0).Eta(), w, -2.5, 2.5,50);
      FillHist(("TightMu" + tag + "_bjet_pt").Data(), el_pt, w,  ptbins, 9);
      FillHist(("TightMu" + tag + "_bjet_pt_eta").Data(), el_pt, fabs(tight_mu.at(0).Eta()),  w, ptbinsb, 7 , etabins, 2);
    }
    else{
      FillHist(("TightMu" + tag + "_0bjet_eta").Data(), tight_mu.at(0).Eta(), w, -2.5, 2.5,50);
      FillHist(("TightMu" + tag + "_0bjet_pt").Data(), el_pt, w,  ptbins,9);
      FillHist(("TightMu" + tag + "_0bjet_pt_eta").Data(), el_pt, fabs(tight_mu.at(0).Eta()),  w, ptbins, 9 , etabins2, 4);
    }
  }
   
  if(loose_mu.size() == 1 && jets.size() >= 1){
    float el_pt = loose_mu.at(0).Pt();
   
    FillHist(("LooseMu" + tag + "_eta").Data(), loose_mu.at(0).Eta(), w, -2.5, 2.5,50);
    FillHist(("LooseMu" + tag + "_pt").Data(), el_pt, w,  ptbins, 9);
    FillHist(("LooseMu" + tag + "_njets").Data(), jets.size(), w, 0.,5.,5);
    FillHist(("LooseMu" + tag + "_nvertices").Data(), eventbase->GetEvent().nVertices(), w, 0., 30., 30);
    FillHist(("LooseMu" + tag + "_nbjet").Data(), nbjet, w, 0., 4.,4);
    FillHist(("LooseMu" + tag + "_pt_eta").Data(), el_pt, fabs(loose_mu.at(0).Eta()),  w, ptbins, 9 , etabins2, 4);

    if(nbjet > 0){
      FillHist(("LooseMu" + tag + "_bjet_eta").Data(), loose_mu.at(0).Eta(), w, -2.5, 2.5,50);
      FillHist(("LooseMu" + tag + "_bjet_pt").Data(), el_pt, w,  ptbins, 9);
      FillHist(("LooseMu" + tag + "_bjet_pt_eta").Data(), el_pt, fabs(loose_mu.at(0).Eta()),  w, ptbinsb,7 , etabins, 2);
    }
    else{
      FillHist(("LooseMu" + tag + "_0bjet_eta").Data(), loose_mu.at(0).Eta(), w, -2.5, 2.5,50);
      FillHist(("LooseMu" + tag + "_0bjet_pt").Data(), el_pt, w,  ptbins, 9);
      FillHist(("LooseMu" + tag + "_0bjet_pt_eta").Data(), el_pt, fabs(loose_mu.at(0).Eta()),  w, ptbins, 9 , etabins2, 4);
    }

  }
  
  
  return;
  
}



void FakeRateCalculator_Mu::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger<< INFO << "Number of events that pass 1 7GeV trigger = " << n_17_pass  << LQLogger::endmsg;
  m_logger<< INFO << "Number of events that pass 17 GeV + jet trigger = " << n_17_jet_pass  << LQLogger::endmsg;
  m_logger<< INFO << "Number of events that pass 17 GeV || jet trigger = " << n_17_17_jet_pass  << LQLogger::endmsg;

}

void FakeRateCalculator_Mu::BeginCycle() throw( LQError ){
  
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

FakeRateCalculator_Mu::~FakeRateCalculator_Mu() {
  
  Message("In FakeRateCalculator_Mu Destructor" , INFO);
  
}



void FakeRateCalculator_Mu::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void FakeRateCalculator_Mu::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this FakeRateCalculator_MuCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void FakeRateCalculator_Mu::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


