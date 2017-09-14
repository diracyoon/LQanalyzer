//Core includes
#include "EventBase.h" 
#include "BaseSelection.h"

//Local includes
#include "Fitter_Test.h"

//Needed to allow inheritance for use in LQCore/core classes
ClassImp(Fitter_Test);
 
//////////

Fitter_Test::Fitter_Test() : AnalyzerCore()   
{
  //To have the correct name in the log:                                                          
  SetLogName("Fitter_Test");
  
  m_logger << INFO << "Construct Fitter_Test." << LQLogger::endmsg;

  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

  Bool_t high_mass_fitter = kFALSE;
  chk_debug = kFALSE;
  Int_t ts_correction_type = 4;
  
  fitter = new Kinematic_Fitter_1(high_mass_fitter, ts_correction_type, chk_debug);
 
  ts_correction = new TS_Correction(ts_correction_type);
}//Fitter_Test::Fitter_Test()

//////////

Fitter_Test::~Fitter_Test() 
{
  m_logger << INFO << "Destruct Fitter_Test." << LQLogger::endmsg;

  delete fitter;
}//Fitter_Test::~Fitter_Test()

//////////

void Fitter_Test::BeginCycle() throw(LQError)
{
  m_logger << INFO << "BeginCycle is called." << LQLogger::endmsg;

  return;
}//void Fitter_Test::BeginCycle()

//////////

void Fitter_Test::BeginEvent() throw(LQError)
{
  m_logger << DEBUG << "BeginEvent is called." << LQLogger::endmsg;

  return;
}//void Fitter_Test::BeginEvent()

//////////

void Fitter_Test::ExecuteEvents() throw(LQError)
{
  m_logger << DEBUG << "RunNumber = " << eventbase->GetEvent().RunNumber() << ", EventNumber = " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "IsData = " << isData << LQLogger::endmsg;

  //Apply the gen weight 
  if(!isData) weight*=MCweight;

  FillCutFlow("NoCut", weight);
  
  /*//////////*/
  /*Vertex cut*/
  /*//////////*/
  if(!eventbase->GetEvent().HasGoodPrimaryVertex()) throw LQError("Fails vertex cuts", LQError::SkipEvent);
  FillCutFlow("VertexCut", weight);
  
  /*///*/
  /*MET*/
  /*///*/
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

  // electron veto cuts
  if(electron_veto_coll.size()>0) throw LQError("Fails electron veto cuts",  LQError::SkipEvent);
  FillCutFlow("ElectronVeto", weight);
  
  /*Muon*/
  std::vector<snu::KMuon> muon_tight_coll;
  muon_tight_coll = GetMuons(BaseSelection::MUON_POG_TIGHT);
  
  std::vector<snu::KMuon> muon_loose_coll;
  muon_loose_coll = GetMuons(BaseSelection::MUON_POG_LOOSE);
  
  //tight muon==1 cut
  if(muon_tight_coll.size()!=1) throw LQError("Fails tight muon==1 cuts", LQError::SkipEvent);
  FillCutFlow("OneTightMuon", weight);

  //extra muon veto cut
  if(muon_loose_coll.size()>1) throw LQError("Fails extra electron veto cuts", LQError::SkipEvent);
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
  for(Int_t i=0; i<n_jet_soft; i++)
    {
      Double_t cvs = jet_soft_coll.at(i).BJetTaggerValue(snu::KJet::CSVv2);

      if(cvs>CSV_THRESHOLD_MEDIUM)  chk_btag[i] = kTRUE;
      //if(cvs>CSV_THRESHOLD_TIGHT)  chk_btag[i] = kTRUE;
    }

  //at least two b jets in leading four jet
  Int_t n_bjet_target = 0;
  for(Int_t i=0; i<4; i++){ if(chk_btag[i]==kTRUE) n_bjet_target++; }
 
  if(n_bjet_target<2) throw LQError("Fails at least two b tagged jet in leading four hard jets", LQError::SkipEvent);
  FillCutFlow("TwoBJets", weight);

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
  
  //parton-jet matching for all jets
  Int_t permutation_truth[4] = {0, 0, 0, 0};
  Bool_t chk_match_all = Truth_Jet_Match(gen_quark, jet_soft_coll, permutation_truth);
  
  if(chk_match_all==kFALSE) throw LQError("Fail parton-jet match for all jets", LQError::SkipEvent);
  FillCutFlow("JetMatchAll", weight);
  
  ////////////////////////
  /*target jet selection*/
  ////////////////////////

  //pt leading four jet
  for(Int_t i=0; i<4; i++){ target_jet[i] = kTRUE; }
  
  //force true-jet
  
  //use only true-jet for fitter test analysis
  Bool_t chk_true_jet_input = Chk_True_Jet_Input(permutation_truth, target_jet, n_jet_soft);
  if(chk_true_jet_input==kFALSE) throw LQError("Fail true-jet input", LQError::SkipEvent);
  FillCutFlow("True-jet", weight);
  
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

  //neutrino pz solution. Is it real or complex
  Bool_t chk_pz_real = fitter->Get_Neutrino_Pz_Sol_Checker();
  if(chk_pz_real==kTRUE) FillHist("Neutrino_Pz_Sol", 1, weight);
  else FillHist("Neutrino_Pz_Sol", 0, weight);

  //get chi2
  Double_t chi2 = fitter->Get_Chi2();
  Double_t chi2_piece[N_CHI2_PIECE];
  fitter->Get_Chi2_Piece(chi2_piece);
  
  //get jet permutation
  Int_t permutation_fitter[4];
  fitter->Get_Permutation(permutation_fitter);
  
  //check jet permutation match
  Bool_t chk_permutation_match = Chk_Permutation_Match(permutation_truth, permutation_fitter);
  
  //get fitted objects
  TLorentzVector fitted_object[6];
  TLorentzVector unfitted_object[6];
  for(Int_t i=0; i<6; i++)
    { 
      fitted_object[i] = fitter->Get_Fitted_Object(i); 
      unfitted_object[i] = fitter->Get_Unfitted_Object(i);
    }
 
  FillHist("MET_Truth_Vs_Measured", gen_neutrino.Pt(), met-gen_neutrino.Pt(), weight, 0, 0, 0, 0, 0, 0);
  FillHist("METPhi_Truth_Vs_Measured", gen_neutrino.Phi(), met_phi-gen_neutrino.Phi(), weight, 0, 0, 0, 0, 0, 0);
  
  FillHist("MET_Truth_Vs_Fitted", gen_neutrino.Pt(), fitted_object[5].Pt()-gen_neutrino.Pt(), weight, 0, 0, 0, 0, 0, 0);
  FillHist("METPhi_Truth_Vs_Fitted", gen_neutrino.Phi(), fitted_object[5].Phi()-gen_neutrino.Phi(), weight, 0, 0, 0, 0, 0, 0);
  FillHist("Pz_Truth_Vs_Fitted", gen_neutrino.Pz(), fitted_object[5].Pz()-gen_neutrino.Pz(), weight, 0, 0, 0, 0, 0, 0);
  FillHist("Eta_Truth_Vs_Fitted", gen_neutrino.Eta(), fitted_object[5].Eta()-gen_neutrino.Eta(), weight, 0, 0, 0, 0, 0, 0);
  
  if(chk_permutation_match==kTRUE)
    {
      FillHist("MET_Truth_Vs_Fitted_OFS", gen_neutrino.Pt(), fitted_object[5].Pt()-gen_neutrino.Pt(), weight, 0, 0, 0, 0, 0, 0);
      FillHist("METPhi_Truth_Vs_Fitted_OFS", gen_neutrino.Phi(), fitted_object[5].Phi()-gen_neutrino.Phi(), weight, 0, 0, 0, 0, 0, 0);
      FillHist("Pz_Truth_Vs_Fitted_OFS", gen_neutrino.Pz(), fitted_object[5].Pz()-gen_neutrino.Pz(), weight, 0, 0, 0, 0, 0, 0);
      FillHist("Eta_Truth_Vs_Fitted_OFS", gen_neutrino.Eta(), fitted_object[5].Eta()-gen_neutrino.Eta(), weight, 0, 0, 0, 0, 0, 0);
    }
  else
    {
      FillHist("MET_Truth_Vs_Fitted_OFF", gen_neutrino.Pt(), fitted_object[5].Pt()-gen_neutrino.Pt(), weight, 0, 0, 0, 0, 0, 0);
      FillHist("METPhi_Truth_Vs_Fitted_OFF", gen_neutrino.Phi(), fitted_object[5].Phi()-gen_neutrino.Phi(), weight, 0, 0, 0, 0, 0, 0);
      FillHist("Pz_Truth_Vs_Fitted_OFF", gen_neutrino.Pz(), fitted_object[5].Pz()-gen_neutrino.Pz(), weight, 0, 0, 0, 0, 0, 0);
      FillHist("Eta_Truth_Vs_Fitted_OFF", gen_neutrino.Eta(), fitted_object[5].Eta()-gen_neutrino.Eta(), weight, 0, 0, 0, 0, 0, 0);
    }

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

  TLorentzVector hadronic_w_ch_unfitted = unfitted_object[2] + unfitted_object[3];
  Double_t hadronic_w_ch_unfitted_mass = hadronic_w_ch_unfitted.M();

  //if jet permutation match fail
  Bool_t chk_permutation_fail_reason;
  if(chk_permutation_match==kFALSE)
    {
      Double_t chi2_truth = 1000;
      Int_t permutation_fitter_truth[4];
      for(Int_t i=0; i<48; i++)
  	{	  
  	  fitter->Get_Permutation(permutation_fitter_truth, "ALL", i);
  	  Double_t chi2_all = fitter->Get_Chi2("ALL", i);
	  
	  //for true permutation
	  if(permutation_truth[0]==permutation_fitter_truth[0] && permutation_truth[1]==permutation_fitter_truth[1])
	    {
	      if(chi2_all<chi2_truth) chi2_truth = chi2_all;
	    }	  
	}//for loop
      
      if(chi2_truth>900) chk_permutation_fail_reason = kFALSE;
      else chk_permutation_fail_reason = kTRUE;
    }
  
  if(n_bjet_target==2)
    {
      if(chk_permutation_match==kTRUE) FillHist("Jet_Permutation_Match_2B", 1, weight);
      else
	{
	  FillHist("Jet_Permutation_Match_2B", 0, weight);
	  
	  if(chk_permutation_fail_reason==kFALSE) FillHist("Jet_Permutation_Match_Fail_Reason_2B", 0, weight);
	  else FillHist("Jet_Permutation_Match_Fail_Reason_2B", 1, weight);
	}
      
      for(Int_t i=0; i<N_CHI2_PIECE; i++)
	{
	  TString hname = "Chi2_Piece_2B_";
	  hname += i;

	  FillHist(hname, chi2_piece[i], weight);
	  if(chk_permutation_match==kFALSE) FillHist(hname + "_OFF", chi2_piece[i], weight);
	  else FillHist(hname + "_OFS", chi2_piece[i], weight);
	}
     
      FillHist("Chi2_2B", chi2, weight);
      if(chk_permutation_match==kFALSE) FillHist("Chi2_2B_OFF", chi2, weight);
      else FillHist("Chi2_2B_OFS", chi2, weight);
      
      FillHist("DiJetMass_2B", hadronic_w_ch_mass, weight);
      if(chk_permutation_match==kFALSE) FillHist("DiJetMass_2B_OFF", hadronic_w_ch_mass, weight);
      else FillHist("DiJetMass_2B_OFS", hadronic_w_ch_mass, weight);

      FillHist("DiJetMass_Chi2_2B", hadronic_w_ch_mass, chi2, weight, 0, 0, 0, 0, 0, 0);
      if(chk_permutation_match==kFALSE) FillHist("DiJetMass_Chi2_2B_OFF", hadronic_w_ch_mass, chi2, weight, 0, 0, 0, 0, 0, 0);
      else FillHist("DiJetMass_Chi2_2B_OFS", hadronic_w_ch_mass, chi2, weight, 0, 0, 0, 0, 0, 0);

      FillHist("Unfitted_DiJetMass_2B", hadronic_w_ch_unfitted_mass, weight);
      
      //goodness cut study
      for(Int_t cut_level=0; cut_level<20; cut_level++)
        {
          Double_t d_cut_level = 0.05*(cut_level+1);
          Bool_t chk_goodness_cut = fitter->Pass_Goodness_Cut(d_cut_level);

          if(chk_goodness_cut==kTRUE)
            {
              TString hname = "DiJetMass_Chi2_2B_";
              hname += (Int_t)(100*d_cut_level);

              FillHist(hname, hadronic_w_ch_mass, chi2, weight, 0, 0, 0, 0, 0, 0);
            }
        }
    }//two b tagged
  else
    {
      if(chk_permutation_match==kTRUE) FillHist("Jet_Permutation_Match_3B", 1, weight);
      else
	{
	  FillHist("Jet_Permutation_Match_3B", 0, weight);
	  
	  if(chk_permutation_fail_reason==kFALSE) FillHist("Jet_Permutation_Match_Fail_Reason_3B", 0, weight);
          else FillHist("Jet_Permutation_Match_Fail_Reason_3B", 1, weight);
	}
	  
      for(Int_t i=0; i<N_CHI2_PIECE; i++)
	{
	  TString hname = "Chi2_Piece_3B_";
	  hname += i;

	  FillHist(hname, chi2_piece[i], weight);
	  if(chk_permutation_match==kFALSE) FillHist(hname + "_OFF", chi2_piece[i], weight);
	  else FillHist(hname + "_OFS", chi2_piece[i], weight);
	}

      FillHist("Chi2_3B", chi2, weight);
      if(chk_permutation_match==kFALSE) FillHist("Chi2_3B_OFF", chi2, weight);
      else FillHist("Chi2_3B_OFS", chi2, weight);

      FillHist("DiJetMass_3B", hadronic_w_ch_mass, weight);
      if(chk_permutation_match==kFALSE) FillHist("DiJetMass_3B_OFF", hadronic_w_ch_mass, weight);
      else FillHist("DiJetMass_3B_OFS", hadronic_w_ch_mass, weight);
      
      FillHist("DiJetMass_Chi2_3B", hadronic_w_ch_mass, chi2, weight, 0, 0, 0, 0, 0, 0);
      if(chk_permutation_match==kFALSE) FillHist("DiJetMass_Chi2_3B_OFF", hadronic_w_ch_mass, chi2, weight, 0, 0, 0, 0, 0, 0);
      else FillHist("DiJetMass_Chi2_3B_OFS", hadronic_w_ch_mass, chi2, weight, 0, 0, 0, 0, 0, 0);
      
      FillHist("Unfitted_DiJetMass_3B", hadronic_w_ch_unfitted_mass, weight);
    
      TLorentzVector jet_1 = ts_correction->Get_Corrected_Jet(jet_hard_coll.at(permutation_truth[1]), 2);
      TLorentzVector jet_2 = ts_correction->Get_Corrected_Jet(jet_hard_coll.at(permutation_truth[2]), 2);

      if(jet_1.Pt()>jet_2.Pt()) FillHist("Two_B_Pt_Comparison", 0, weight);
      else FillHist("Two_B_Pt_Comparison", 1, weight);
      
      if(jet_1.P()>jet_2.P()) FillHist("Two_B_P_Comparison", 0, weight);
      else FillHist("Two_B_P_Comparison", 1, weight);
    }//three b tagged

  FillCLHist(muhist, "Muon", muon_tight_coll, weight);
  FillCLHist(jethist, "Jet", jet_hard_coll, weight);


  /*delete*/
  delete[] chk_btag;
  delete[] target_jet;

  return;
}//void Fitter_Test::ExecuteEvents()

//////////

void Fitter_Test::EndCycle() throw(LQError)
{
  m_logger << INFO << "EndCyle is called." << LQLogger::endmsg;

  return;
}//void Fitter_Test::EndCycle()

//////////

void Fitter_Test::InitialiseAnalysis() throw(LQError)
{
  m_logger << INFO << "Initialise Fitter_Test analysis." << LQLogger::endmsg;

  //Initialise histograms                                                                             
  MakeHistograms();
  
  //jet permutation matching efficiency
  MakeHistograms("Jet_Permutation_Match_2B", 2, -0.5, 1.5);
  MakeHistograms("Jet_Permutation_Match_3B", 2, -0.5, 1.5);
  
  MakeHistograms("Jet_Permutation_Match_Fail_Reason_2B", 2, -0.5, 1.5);
  MakeHistograms("Jet_Permutation_Match_Fail_Reason_3B", 2, -0.5, 1.5);

  //lepton

  //neutrino
  MakeHistograms("Neutrino_Pz_Sol", 2, -0.5, 1.5);
  MakeHistograms("Neutrino_Pz_Sol_OFS", 2, -0.5, 1.5);
  MakeHistograms("Neutrino_Pz_Sol_OFF", 2, -0.5, 1.5);

  MakeHistograms2D("MET_Truth_Vs_Measured", 100, 0, 200, 100, -100, 100);
  MakeHistograms2D("MET_Truth_Vs_Fitted", 100, 0, 200, 100, -100, 100);
  MakeHistograms2D("MET_Truth_Vs_Fitted_OFS", 100, 0, 200, 100, -100, 100);
  MakeHistograms2D("MET_Truth_Vs_Fitted_OFF", 100, 0, 200, 100, -100, 100);

  MakeHistograms2D("METPhi_Truth_Vs_Measured", 100, -4, 4, 100, -4, 4);
  MakeHistograms2D("METPhi_Truth_Vs_Fitted", 100, -4, 4, 100, -4, 4);
  MakeHistograms2D("METPhi_Truth_Vs_Fitted_OFS", 100, -4, 4, 100, -4, 4);
  MakeHistograms2D("METPhi_Truth_Vs_Fitted_OFF", 100, -4, 4, 100, -4, 4);

  MakeHistograms2D("Pz_Truth_Vs_Fitted", 100, -200, 200, 100, -200, 200);
  MakeHistograms2D("Pz_Truth_Vs_Fitted_OFS", 100, -200, 200, 100, -200, 200);
  MakeHistograms2D("Pz_Truth_Vs_Fitted_OFF", 100, -200, 200, 100, -200, 200);
  
  MakeHistograms2D("Eta_Truth_Vs_Fitted", 100, -5, 5, 100, -5, 5);
  MakeHistograms2D("Eta_Truth_Vs_Fitted_OFS", 100, -5, 5, 100, -5, 5);
  MakeHistograms2D("Eta_Truth_Vs_Fitted_OFF", 100, -5, 5, 100, -5, 5);

  
  //chi2
  for(Int_t i=0; i<N_CHI2_PIECE; i++)
    {
      //2b
      TString hname = "Chi2_Piece_2B_";
      hname += i;

      MakeHistograms(hname, 100, 0, 50);
      MakeHistograms(hname + "_OFF", 100, 0, 50);
      MakeHistograms(hname + "_OFS", 100, 0, 50);
   
      //3b
      hname = "Chi2_Piece_3B_";
      hname += i;

      MakeHistograms(hname, 100, 0, 50);
      MakeHistograms(hname + "_OFF", 100, 0, 50);
      MakeHistograms(hname + "_OFS", 100, 0, 50);
    }

  MakeHistograms("Chi2_2B", 100, 0, 50);
  MakeHistograms("Chi2_2B_OFS", 100, 0, 50);
  MakeHistograms("Chi2_2B_OFF", 100, 0, 50);

  MakeHistograms("Chi2_3B", 100, 0, 50);
  MakeHistograms("Chi2_3B_OFS", 100, 0, 50);
  MakeHistograms("Chi2_3B_OFF", 100, 0, 50);

  //dijet mass
  MakeHistograms("DiJetMass_2B", 50, 0, 200); 
  MakeHistograms("DiJetMass_2B_OFS", 50, 0, 200);
  MakeHistograms("DiJetMass_2B_OFF", 50, 0, 200);
  
  MakeHistograms("DiJetMass_3B", 50, 0, 200);
  MakeHistograms("DiJetMass_3B_OFS", 50, 0, 200);
  MakeHistograms("DiJetMass_3B_OFF", 50, 0, 200);

  MakeHistograms("Unfitted_DiJetMass_2B", 50, 0, 200);
  MakeHistograms("Unfitted_DiJetMass_3B", 50, 0, 200);
  
  //chi2 & dijet mass 2d
  MakeHistograms2D("DiJetMass_Chi2_2B", 50, 0, 200, 100, 0, 50);
  MakeHistograms2D("DiJetMass_Chi2_2B_OFS", 50, 0, 200, 100, 0, 50);
  MakeHistograms2D("DiJetMass_Chi2_2B_OFF", 50, 0, 200, 100, 0, 50);

  MakeHistograms2D("DiJetMass_Chi2_3B", 50, 0, 200, 100, 0, 50);
  MakeHistograms2D("DiJetMass_Chi2_3B_OFS", 50, 0, 200, 100, 0, 50);
  MakeHistograms2D("DiJetMass_Chi2_3B_OFF", 50, 0, 200, 100, 0, 50);

  for(Int_t cut_level=0; cut_level<20; cut_level++)
    {
      Double_t d_cut_level = 0.05*(cut_level+1);
      TString hname = "DiJetMass_Chi2_2B_";
      hname += (Int_t)(100*d_cut_level);

      MakeHistograms2D(hname, 50, 0, 200, 100, 0, 50);
    }
  
  //others
  MakeHistograms("Two_B_Pt_Comparison", 2, -0.5, 1.5);
  MakeHistograms("Two_B_P_Comparison", 2, -0.5, 1.5);
  
  MakeHistograms("Leptonic_W_Mass", 50, 55, 105);
  MakeHistograms("Hadronic_Top_Mass", 200, 150, 190);
  MakeHistograms("Leptonic_Top_Mass", 200, 150, 190);
  
  MakeCleverHistograms(muhist, "Muon");
  MakeCleverHistograms(jethist, "Jet");

  return;
}//void Fitter_Test::InitialiseAnalysis()

//////////
