//Core includes
#include "EventBase.h" 
#include "BaseSelection.h"

//Local includes
#include "Jet_Selection_Test_0.h"

//Needed to allow inheritance for use in LQCore/core classes
ClassImp(Jet_Selection_Test_0);
 
//////////

Jet_Selection_Test_0::Jet_Selection_Test_0() : AnalyzerCore()   
{
  //To have the correct name in the log:                                                          
  SetLogName("Jet_Selection_Test_0");
  
  m_logger << INFO << "Construct Jet_Selection_Test_0." << LQLogger::endmsg;

  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

  Bool_t high_mass_fitter = kFALSE;
  chk_debug = kFALSE;

  fitter = new Kinematic_Fitter_1(high_mass_fitter, chk_debug);
}//Jet_Selection_Test_0::Jet_Selection_Test_0()

//////////

Jet_Selection_Test_0::~Jet_Selection_Test_0() 
{
  m_logger << INFO << "Destruct Jet_Selection_Test_0." << LQLogger::endmsg;

  delete fitter;
}//Jet_Selection_Test_0::~Jet_Selection_Test_0()

//////////

void Jet_Selection_Test_0::BeginCycle() throw(LQError)
{
  m_logger << INFO << "BeginCycle is called." << LQLogger::endmsg;

  return;
}//void Jet_Selection_Test_0::BeginCycle()

//////////

void Jet_Selection_Test_0::BeginEvent() throw(LQError)
{
  m_logger << DEBUG << "BeginEvent is called." << LQLogger::endmsg;

  return;
}//void Jet_Selection_Test_0::BeginEvent()

//////////

void Jet_Selection_Test_0::ExecuteEvents() throw(LQError)
{
  m_logger << DEBUG << "RunNumber = " << eventbase->GetEvent().RunNumber() << ", EventNumber = " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "IsData = " << isData << LQLogger::endmsg;

  //Apply the gen weight 
  if(!isData) weight*=MCweight;

  FillCutFlow("NoCut", weight);

  /*////////////////////*/
  /*Gen Truth*/
  /*////////////////////*/
    
  std::vector<snu::KTruth> gen_truth_coll;
  eventbase->GetTruthSel()->Selection(gen_truth_coll);

  snu::KTruth gen_quark[4];
  snu::KTruth gen_neutrino;
  snu::KTruth gen_lepton;
  Bool_t chk_semi_leptonic = Search_Truth_Value(gen_truth_coll, gen_quark, gen_neutrino, gen_lepton);

  if(chk_semi_leptonic==kFALSE) throw LQError("No semi-leptonic decay", LQError::SkipEvent);
  FillCutFlow("SemiLeptonic", weight);
  
  //Vertex cut
  if(!eventbase->GetEvent().HasGoodPrimaryVertex()) throw LQError("Fails vertex cuts", LQError::SkipEvent);
  FillCutFlow("VertexCut", weight);
  
  /*MET*/
  //Pass basic MET filter
  if(!PassMETFilter()) throw LQError("Fail basic MET filter.", LQError::SkipEvent);
  FillCutFlow("METCut", weight);
  
  //Mininum MET requirement
  Double_t met = eventbase->GetEvent().MET();
  Double_t met_phi = eventbase->GetEvent().METPhi();

  if(met<20) throw LQError("Fail minimum MET requirement", LQError::SkipEvent);
  FillCutFlow("MininumMET", weight);
  
  //construct met vector
  TLorentzVector met_vector;
  met_vector.SetPxPyPzE(met*TMath::Cos(met_phi), met*TMath::Sin(met_phi), 0, met);

  /*Electron*/
  std::vector<snu::KElectron> electron_tight_coll;
  electron_tight_coll = GetElectrons(BaseSelection::ELECTRON_POG_TIGHT);

  std::vector<snu::KElectron> electron_veto_coll;
  electron_veto_coll = GetElectrons(BaseSelection::ELECTRON_POG_VETO);

  //electron veto cuts
  if(electron_veto_coll.size()>0) throw LQError("Fails electron veto cuts",  LQError::SkipEvent);
  FillCutFlow("ElectronVeto", weight);
  
  /*Muon*/
  std::vector<snu::KMuon> muon_tight_coll;
  muon_tight_coll = GetMuons(BaseSelection::MUON_POG_TIGHT);
  
  std::vector<snu::KMuon> muon_veto_coll;
  muon_veto_coll = GetMuons(BaseSelection::MUON_POG_LOOSE);
  
  //tight muon==1 cut
  if(muon_tight_coll.size()!=1) throw LQError("Fails tight muon==1 cuts", LQError::SkipEvent);
  FillCutFlow("OneTightMuon", weight);

  //extra muon veto cut
  if(muon_veto_coll.size()>1) throw LQError("Fails extra electron veto cuts", LQError::SkipEvent);
  FillCutFlow("ExtraMuonVeto", weight);

  //construct muon vector
  TLorentzVector muon_vector(muon_tight_coll.at(0));

  /*Jet*/
  std::vector<snu::KJet> jet_hard_coll;
  jet_hard_coll = GetJets("JET_NOCUT", 30, 2.5);

  std::vector<snu::KJet> jet_soft_coll;
  jet_soft_coll = GetJets("JET_NOCUT", 20, 2.5);

  //reordering jet by pt
  sort(jet_hard_coll.begin(), jet_hard_coll.end(), Compare_Jet_Pt);
  sort(jet_soft_coll.begin(), jet_soft_coll.end(), Compare_Jet_Pt);

  Int_t n_jet_soft = jet_soft_coll.size();
  
  vector<TLorentzVector> jet_soft_vec;
  for(Int_t i=0; i<n_jet_soft; i++){ jet_soft_vec.push_back(jet_soft_coll.at(i)); }
  
  //at least four hard jets cut
  if(jet_hard_coll.size()<4) throw LQError("Fails at least four hard jets cuts.", LQError::SkipEvent);
  FillCutFlow("FourJets", weight);
  
  Bool_t* chk_btag = new Bool_t[n_jet_soft];
  Bool_t* target_jet = new Bool_t[n_jet_soft];
  for(Int_t i=0; i<n_jet_soft; i++)
    {
      chk_btag[i] = kFALSE;
      target_jet[i] = kFALSE; 
    }

  //find b tag jet
  Int_t n_bjet_soft = 0;
  for(Int_t i=0; i<n_jet_soft; i++)
    {
      Double_t csv = jet_soft_coll.at(i).BJetTaggerValue(snu::KJet::CSVv2);

      if(csv>CSV_THRESHOLD_MEDIUM)
	{
	  chk_btag[i] = kTRUE;
	  n_bjet_soft++;
	}
    }
 
  //at least two b tagged jets
  Int_t n_bjet_target = 0;
  for(Int_t i=0; i<4; i++){ if(chk_btag[i] == kTRUE) n_bjet_target++; }
  
  //parton-jet matching for all jets
  Int_t permutation_truth[4] = {999, 999, 999, 999};
  Bool_t chk_match_all = Truth_Jet_Match(gen_quark, jet_soft_coll, permutation_truth, 0.3);

  //if(chk_match_all==kFALSE) throw LQError("Fails jet match", LQError::SkipEvent);
  //FillCutFlow("JetMatchAll", weight);
  
  /////////////////////////
  /*target jets selection*/
  /////////////////////////
  
  //leading four jets
  for(Int_t i=0; i<4; i++){ target_jet[i] = kTRUE; }
    
  //////////
  /*Fitter*/
  //////////
  
  fitter->Set(met_vector, muon_vector, jet_soft_vec, target_jet, chk_btag);
  fitter->Fit();

  ////////////////
  /*Save results*/
  ////////////////

  Bool_t chk_convergence = fitter->Get_Convergence_Checker();
  if(chk_convergence==kFALSE) throw LQError("Fitter Fail", LQError::SkipEvent);

  Double_t chi2 = fitter->Get_Chi2();

  Int_t permutation_fitter[4];
  fitter->Get_Permutation(permutation_fitter);
   
  //check true-jet insertion
  Bool_t chk_true_jet_input = Chk_True_Jet_Input(permutation_truth, target_jet, n_jet_soft);
 
  //check jet permutation 
  Bool_t chk_permutation_match = Chk_Permutation_Match(permutation_truth, permutation_fitter);
 
  Double_t chi2_piece[N_CHI2_PIECE];
  fitter->Get_Chi2_Piece(chi2_piece);
    
  TLorentzVector fitted_object[6];
  for(Int_t i=0; i<6; i++){ fitted_object[i] = fitter->Get_Fitted_Object(i); }
  
  FillHist("MET_Truth_Vs_Measured", gen_neutrino.Pt(), met-gen_neutrino.Pt(), weight, 0, 0, 0, 0, 0, 0);
  FillHist("METPhi_Truth_Vs_Measured", gen_neutrino.Phi(), met_phi-gen_neutrino.Phi(), weight, 0, 0, 0, 0, 0, 0);

  FillHist("MET_Truth_Vs_Fitted", gen_neutrino.Pt(), fitted_object[5].Pt()-gen_neutrino.Pt(), weight, 0, 0, 0, 0, 0, 0);
  FillHist("METPhi_Truth_Vs_Fitted", gen_neutrino.Phi(), fitted_object[5].Phi()-gen_neutrino.Phi(), weight, 0, 0, 0, 0, 0, 0);
  FillHist("Pz_Truth_Vs_Fitted", gen_neutrino.Pz(), fitted_object[5].Pz()-gen_neutrino.Pz(), weight, 0, 0, 0, 0, 0, 0);
  FillHist("Eta_Truth_Vs_Fitted", gen_neutrino.Eta(), fitted_object[5].Eta()-gen_neutrino.Eta(), weight, 0, 0, 0, 0, 0, 0);

  TLorentzVector leptonic_w = fitted_object[4] + fitted_object[5];
  Double_t leptonic_w_mass = leptonic_w.M();
    
  TLorentzVector hadronic_top = fitted_object[1] + fitted_object[2] + fitted_object[3];
  Double_t hadronic_t_mass = hadronic_top.M();

  TLorentzVector leptonic_top = fitted_object[0] + fitted_object[4] + fitted_object[5];
  Double_t leptonic_t_mass = leptonic_top.M();
  
  FillHist("Leptonic_W_Mass", leptonic_w_mass, weight);
  FillHist("Hadronic_Top_Mass", hadronic_t_mass, weight);
  FillHist("Leptonic_Top_Mass", leptonic_t_mass, weight);
  
  TLorentzVector hadronic_w_ch = fitted_object[2] + fitted_object[3];
  Double_t hadronic_w_ch_mass = hadronic_w_ch.M();

  if(n_bjet_target==2)
    {
      if(chk_true_jet_input==kTRUE) FillHist("True_Jet_Input_2B", 1, weight);
      else FillHist("True_Jet_Input_2B", 0, weight);

      if(chk_true_jet_input==kTRUE && chk_permutation_match==kTRUE) FillHist("Jet_Permutation_Match_2B", 1, weight);
      else FillHist("Jet_Permutation_Match_2B", 0, weight);

      for(Int_t i=0; i<N_CHI2_PIECE; i++)
	{
	  TString hname = "Chi2_Piece_2B_";
	  hname += i;

	  FillHist(hname, chi2_piece[i], weight);
	  if(chk_true_jet_input==kFALSE) FillHist(hname + "_TF", chi2_piece[i], weight);
	  else if(chk_permutation_match==kFALSE) FillHist(hname + "_TS_PF", chi2_piece[i], weight);
	  else FillHist(hname + "_PS", chi2_piece[i], weight);
	}
      
      FillHist("Chi2_2B", chi2, weight);
      if(chk_true_jet_input==kFALSE) FillHist("Chi2_2B_TF", chi2, weight);
      else if(chk_permutation_match==kFALSE) FillHist("Chi2_2B_TS_PF", chi2, weight);
      else FillHist("Chi2_2B_PS", chi2, weight); 
	
      FillHist("DiJetMass_2B", hadronic_w_ch_mass, weight);
      if(chk_true_jet_input==kFALSE) FillHist("DiJetMass_2B_TF", hadronic_w_ch_mass, weight);
      else if(chk_permutation_match==kFALSE) FillHist("DiJetMass_2B_TS_PF", hadronic_w_ch_mass, weight);
      else FillHist("DiJetMass_2B_PS", hadronic_w_ch_mass, weight);
      
      FillHist("DiJetMass_Chi2_2B", hadronic_w_ch_mass, chi2, weight, 0, 0, 0, 0, 0, 0);
      if(chk_true_jet_input==kFALSE) FillHist("DiJetMass_Chi2_2B_TF", hadronic_w_ch_mass, chi2, weight, 0, 0, 0, 0, 0, 0);
      else if(chk_permutation_match==kFALSE) FillHist("DiJetMass_Chi2_2B_TS_PF", hadronic_w_ch_mass, chi2, weight, 0, 0, 0, 0, 0, 0); 
      else FillHist("DiJetMass_Chi2_2B_PS", hadronic_w_ch_mass, chi2, weight, 0, 0, 0, 0, 0, 0);
    }//2b tagged events
  else
    {
      if(chk_true_jet_input==kTRUE) FillHist("True_Jet_Input_3B", 1, weight);
      else FillHist("True_Jet_Input_3B", 0, weight);
      
      if(chk_true_jet_input==kTRUE && chk_permutation_match==kTRUE) FillHist("Jet_Permutation_Match_3B", 1, weight);
      else FillHist("Jet_Permutation_Match_3B", 0, weight);
      
      for(Int_t i=0; i<N_CHI2_PIECE; i++)
	{
	  TString hname = "Chi2_Piece_3B_";
	  hname += i;

	  FillHist(hname, chi2_piece[i], weight);
	  if(chk_true_jet_input==kFALSE) FillHist(hname + "_TF", chi2_piece[i], weight);
	  else if(chk_permutation_match==kFALSE) FillHist(hname + "_TS_PF", chi2_piece[i], weight);
	  else FillHist(hname + "_PS", chi2_piece[i], weight);
	}

      FillHist("Chi2_3B", chi2, weight);
      if(chk_true_jet_input==kFALSE) FillHist("Chi2_3B_TF", chi2, weight);
      else if(chk_permutation_match==kFALSE) FillHist("Chi2_3B_TS_PF", chi2, weight);
      else FillHist("Chi2_3B_PS", chi2, weight);
      
      FillHist("DiJetMass_3B", hadronic_w_ch_mass, weight);
      if(chk_true_jet_input==kFALSE) FillHist("DiJetMass_3B_TF", hadronic_w_ch_mass, weight);
      else if(chk_permutation_match==kFALSE) FillHist("DiJetMass_3B_TS_PF", hadronic_w_ch_mass, weight);
      else FillHist("DiJetMass_3B_PS", hadronic_w_ch_mass, weight);
      
      FillHist("DiJetMass_Chi2_3B", hadronic_w_ch_mass, chi2, weight, 0, 0, 0, 0, 0, 0);
      if(chk_true_jet_input==kFALSE) FillHist("DiJetMass_Chi2_3B_TF", hadronic_w_ch_mass, chi2, weight, 0, 0, 0, 0, 0, 0);
      else if(chk_permutation_match==kFALSE) FillHist("DiJetMass_Chi2_3B_TS_PF", hadronic_w_ch_mass, chi2, weight, 0, 0, 0, 0, 0, 0);
      else FillHist("DiJetMass_Chi2_3B_PS", hadronic_w_ch_mass, chi2, weight, 0, 0, 0, 0, 0, 0);
    }//3b tagged events
  
  
  FillCLHist(muhist, "Muon", muon_tight_coll, weight);
  FillCLHist(jethist, "Jet", jet_hard_coll, weight);
  
  /*delete*/

  delete[] chk_btag;
  delete[] target_jet;


  return;
}//void Jet_Selection_Test_0::ExecuteEvents()

//////////

void Jet_Selection_Test_0::EndCycle() throw(LQError)
{
  m_logger << INFO << "EndCyle is called." << LQLogger::endmsg;

  return;
}//void Jet_Selection_Test_0::EndCycle()

//////////

void Jet_Selection_Test_0::InitialiseAnalysis() throw(LQError)
{
  m_logger << INFO << "Initialise Jet_Selection_Test_0 analysis." << LQLogger::endmsg;

  //Initialise histograms                                                                             
  MakeHistograms();
  
  //true jet input prob.
  MakeHistograms("True_Jet_Input_2B", 2, -0.5, 1.5);
  MakeHistograms("True_Jet_Input_3B", 2, -0.5, 1.5);
  
  //jet permutation prob.
  MakeHistograms("Jet_Permutation_Match_2B", 2, -0.5, 1.5);
  MakeHistograms("Jet_Permutation_Match_3B", 2, -0.5, 1.5);
  
  //lepton

  //neutrino 
  MakeHistograms2D("MET_Truth_Vs_Measured", 100, 0, 200, 100, -100, 100);
  MakeHistograms2D("METPhi_Truth_Vs_Measured", 100, -4, 4, 100, -4, 4);

  MakeHistograms2D("MET_Truth_Vs_Fitted", 100, 0, 200, 100, -100, 100);
  MakeHistograms2D("METPhi_Truth_Vs_Fitted", 100, -4, 4, 100, -4, 4);
  MakeHistograms2D("Pz_Truth_Vs_Fitted", 100, -200, 200, 100, -200, 200);
  MakeHistograms2D("Eta_Truth_Vs_Fitted", 100, -5, 5, 100, -5, 5);

  //chi2
  for(Int_t i=0; i<N_CHI2_PIECE; i++)
    {
      //2b
      TString hname = "Chi2_Piece_2B_";
      hname += i;

      MakeHistograms(hname, 100, 0, 50);
      MakeHistograms(hname + "_TF", 100, 0, 50);
      MakeHistograms(hname + "_TS_PF", 100, 0, 50);
      MakeHistograms(hname + "_PS", 100, 0, 50);
    
      //3b
      hname = "Chi2_Piece_3B_";
      hname += i;
      
      MakeHistograms(hname, 100, 0, 50);
      MakeHistograms(hname + "_TF", 100, 0, 50);
      MakeHistograms(hname + "_TS_PF", 100, 0, 50);
      MakeHistograms(hname + "_PS", 100, 0, 50);
    }

  MakeHistograms("Chi2_2B", 100, 0, 50);
  MakeHistograms("Chi2_2B_TF", 100, 0, 50);
  MakeHistograms("Chi2_2B_TS_PF", 100, 0, 50);
  MakeHistograms("Chi2_2B_PS", 100, 0, 50);

  MakeHistograms("Chi2_3B", 100, 0, 50);
  MakeHistograms("Chi2_3B_TF", 100, 0, 50);
  MakeHistograms("Chi2_3B_TS_PF", 100, 0, 50);
  MakeHistograms("Chi2_3B_PS", 100, 0, 50);

  //dijet mass
  MakeHistograms("DiJetMass_2B", 50, 0, 200);
  MakeHistograms("DiJetMass_2B_TF", 50, 0, 200);
  MakeHistograms("DiJetMass_2B_TS_PF", 50, 0, 200);
  MakeHistograms("DiJetMass_2B_PS", 50, 0, 200);

  MakeHistograms("DiJetMass_3B", 50, 0, 200);
  MakeHistograms("DiJetMass_3B_TF", 50, 0, 200);
  MakeHistograms("DiJetMass_3B_TS_PF", 50, 0, 200);
  MakeHistograms("DiJetMass_3B_PS", 50, 0, 200);
  
  //chi2 & dijet mass 2d
  MakeHistograms2D("DiJetMass_Chi2_2B", 50, 0, 200, 100, 0, 50);
  MakeHistograms2D("DiJetMass_Chi2_2B_TF", 50, 0, 200, 100, 0, 50);
  MakeHistograms2D("DiJetMass_Chi2_2B_TS_PF", 50, 0, 200, 100, 0, 50);
  MakeHistograms2D("DiJetMass_Chi2_2B_PS", 50, 0, 200, 100, 0, 50);
  
  MakeHistograms2D("DiJetMass_Chi2_3B", 50, 0, 200, 100, 0, 50);
  MakeHistograms2D("DiJetMass_Chi2_3B_TF", 50, 0, 200, 100, 0, 50);
  MakeHistograms2D("DiJetMass_Chi2_3B_TS_PF", 50, 0, 200, 100, 0, 50);
  MakeHistograms2D("DiJetMass_Chi2_3B_PS", 50, 0, 200, 100, 0, 50);
  
  //others  
  MakeHistograms("Leptonic_W_Mass", 50, 55, 105);
  MakeHistograms("Hadronic_Top_Mass", 200, 150, 190);
  MakeHistograms("Leptonic_Top_Mass", 200, 150, 190);

  MakeCleverHistograms(muhist, "Muon");
  MakeCleverHistograms(jethist, "Jet");

  return;
}//void Jet_Selection_Test_0::InitialiseAnalysis()

//////////
