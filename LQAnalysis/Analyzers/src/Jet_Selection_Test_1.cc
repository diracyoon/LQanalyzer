//Core includes
#include "EventBase.h" 
#include "BaseSelection.h"

//Local includes
#include "Jet_Selection_Test_1.h"

//Needed to allow inheritance for use in LQCore/core classes
ClassImp(Jet_Selection_Test_1);
 
//////////

Jet_Selection_Test_1::Jet_Selection_Test_1() : AnalyzerCore()   
{
  //To have the correct name in the log:                                                          
  SetLogName("Jet_Selection_Test_1");
  
  m_logger << INFO << "Construct Jet_Selection_Test_1." << LQLogger::endmsg;

  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

  chk_debug = kFALSE;

  fitter = new Kinematic_Fitter(chk_debug);
  //fitter = new Kinematic_Fitter_Old(chk_debug);
}//Jet_Selection_Test_1::Jet_Selection_Test_1()

//////////

Jet_Selection_Test_1::~Jet_Selection_Test_1() 
{
  m_logger << INFO << "Destruct Jet_Selection_Test_1." << LQLogger::endmsg;

  delete fitter;
}//Jet_Selection_Test_1::~Jet_Selection_Test_1()

//////////

void Jet_Selection_Test_1::BeginCycle() throw(LQError)
{
  m_logger << INFO << "BeginCycle is called." << LQLogger::endmsg;

  return;
}//void Jet_Selection_Test_1::BeginCycle()

//////////

void Jet_Selection_Test_1::BeginEvent() throw(LQError)
{
  m_logger << DEBUG << "BeginEvent is called." << LQLogger::endmsg;

  return;
}//void Jet_Selection_Test_1::BeginEvent()

//////////

void Jet_Selection_Test_1::ExecuteEvents() throw(LQError)
{
  m_logger << DEBUG << "RunNumber = " << eventbase->GetEvent().RunNumber() << ", EventNumber = " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "IsData = " << isData << LQLogger::endmsg;

  //Apply the gen weight 
  if(!isData) weight*=MCweight;

  FillCutFlow("NoCut", weight);

  /*////////////////////*/
  /*////////////////////*/
  /*Gen Truth*/
  /*////////////////////*/
  /*////////////////////*/
    
  std::vector<snu::KTruth> gen_parton_coll;
  eventbase->GetTruthSel()->Selection(gen_parton_coll);

  Int_t n_parton = gen_parton_coll.size();
 
  Int_t gen_truth_index[10];
  for(Int_t i=0; i<10; i++){ gen_truth_index[i] = BLANK; }
  
  snu::KTruth gen_truth[10];
  
  for(Int_t i=0; i<n_parton; i++)
    {
      snu::KTruth truth = gen_parton_coll.at(i);
      
      Int_t gen_index = i; 
      Int_t pdg_id =  truth.PdgId();
      Int_t index_mother = truth.IndexMother();

      //cout << gen_index << "\t" << pdg_id << "\t" << index_mother << endl;
      
      //top observed
      if(pdg_id==6 && gen_truth_index[TOP]==BLANK)
	{
	  gen_truth[TOP] = truth;
	  gen_truth_index[TOP]= gen_index;
	}
      
      //anti_top observed
      if(pdg_id==-6 && gen_truth_index[A_TOP]==BLANK)
	{
	  gen_truth[A_TOP] = truth;
          gen_truth_index[A_TOP]= gen_index;
	}
      
      //t -> b w+ or t -> b h+
      if(index_mother==gen_truth_index[TOP] && pdg_id!=gen_truth[TOP].PdgId())
	{
	  //bottom from top decay observed
	  if(pdg_id==5)
	    {
	      gen_truth[BOTTOM] = truth;
	      gen_truth_index[BOTTOM]= gen_index;
	    }

	  //w+ or h+ from top decay observed
	  else
	    {
	      gen_truth[D_0] = truth;
	      gen_truth_index[D_0]= gen_index;
	    }
	}

      //D_0 -> D_0_0 D_0_1
      if(index_mother==gen_truth_index[D_0] && pdg_id!=gen_truth[D_0].PdgId())
	{
	  if(gen_truth_index[D_0_0]==BLANK)
	    {
	      gen_truth[D_0_0] = truth;
              gen_truth_index[D_0_0]= gen_index;
	    }
	  else
	    {
	      gen_truth[D_0_1] = truth;
              gen_truth_index[D_0_1]= gen_index;
	    }
	}

      //anti_t -> anti_b w- or anti_t -> anti_b h-
      if(index_mother==gen_truth_index[A_TOP] && pdg_id!=gen_truth[A_TOP].PdgId())
        {
          //anti_bottom from top decay observed
          if(pdg_id==-5)
            {
	      gen_truth[A_BOTTOM]= truth;
	      gen_truth_index[A_BOTTOM]= gen_index;
            }

          //w- or h- from top decay observed 
          else
            {
	      gen_truth[D_1] = truth;
              gen_truth_index[D_1]= gen_index;
            }
        }
      
      //D_1 -> D_1_0 D_1_1
      if(index_mother==gen_truth_index[D_1] && pdg_id!=gen_truth[D_1].PdgId())
	{
          if(gen_truth_index[D_1_0]==BLANK)
            {
	      gen_truth[D_1_0] = truth;
              gen_truth_index[D_1_0]= gen_index;
            }
          else
	    {
	      gen_truth[D_1_1] = truth;
              gen_truth_index[D_1_1]= gen_index;
            }
        }

      //tracing parton shower
      for(Int_t j=0; j<10; j++)
      {
        if(index_mother==gen_truth_index[j] && pdg_id==gen_truth[j].PdgId())
          {
	      gen_truth[j] = truth;
	      gen_truth_index[j] = gen_index;
	  }
      }      
    }//for loop over ntruth
  
  //cout << endl;
  //for(Int_t i=0; i<10; i++){ cout << gen_truth_index[i] << "\t" << gen_truth[i].PdgId() << endl; }
  //cout << endl;
  
  //Check one lepton decay and one hadronic decay
  Bool_t chk_leptonic_decay_0 = kFALSE;
  Bool_t chk_leptonic_decay_1 = kFALSE;
  
  if(gen_truth[D_0_0].PdgId()==-11 || gen_truth[D_0_0].PdgId()==-13 || gen_truth[D_0_1].PdgId()==-11 || gen_truth[D_0_1].PdgId()==-13) chk_leptonic_decay_0 = kTRUE;
  if(gen_truth[D_1_0].PdgId()==11 || gen_truth[D_1_0].PdgId()==13 || gen_truth[D_1_1].PdgId()==11 || gen_truth[D_1_1].PdgId()==13) chk_leptonic_decay_1 = kTRUE;
  
  if((chk_leptonic_decay_0==kTRUE && chk_leptonic_decay_1==kTRUE) || (chk_leptonic_decay_0==kFALSE && chk_leptonic_decay_1==kFALSE)) throw LQError("Fail truth cuts", LQError::SkipEvent);
  FillCutFlow("Truth", weight);
 
  //constuct gen parton array for easy handling
  snu::KTruth four_gen_truth[4];
  snu::KTruth gen_neutrino;
  snu::KTruth gen_lepton;
  if(chk_leptonic_decay_0==kFALSE)
    {
      four_gen_truth[0] = gen_truth[A_BOTTOM];
      four_gen_truth[1] = gen_truth[BOTTOM];
      four_gen_truth[2] = gen_truth[D_0_0];
      four_gen_truth[3] = gen_truth[D_0_1];
      
      if(gen_truth[D_1_0].PdgId()==11 || gen_truth[D_1_0].PdgId()==13)
	{
	  gen_neutrino = gen_truth[D_1_1];
	  gen_lepton = gen_truth[D_1_0];
	}
      else
	{
	  gen_neutrino = gen_truth[D_1_0];
	  gen_lepton = gen_truth[D_1_1];
	}
    }
  else
    {
      four_gen_truth[0] = gen_truth[BOTTOM];
      four_gen_truth[1] = gen_truth[A_BOTTOM];
      four_gen_truth[2] = gen_truth[D_1_0];
      four_gen_truth[3] = gen_truth[D_1_1];

      if(gen_truth[D_0_0].PdgId()==-11 || gen_truth[D_0_0].PdgId()==-13)
	{
	  gen_neutrino = gen_truth[D_0_1];
	  gen_lepton = gen_truth[D_0_0];
	}
      else
	{
	  gen_neutrino = gen_truth[D_0_0];
	  gen_lepton = gen_truth[D_0_1];
	}
    }
  
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

  // electron veto cuts
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

  //at least four hard jets cut
  if(jet_hard_coll.size()<4) throw LQError("Fails at least four hard jets cuts.", LQError::SkipEvent);
  FillCutFlow("FourJets", weight);

  Int_t n_jet_soft = jet_soft_coll.size();
  
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
      Double_t cvs = jet_soft_coll.at(i).BJetTaggerValue(snu::KJet::CSVv2);

      if(cvs>CSV_THRESHOLD_MEDIUM)
	{
	  chk_btag[i] = kTRUE;
	  n_bjet_soft++;
	}
    }
  
  if(n_bjet_soft<2) throw LQError("Fails at least two b tagged jets", LQError::SkipEvent);
  FillCutFlow("TwoBJets", weight);
  
  //parton-jet matching for all jets
  Int_t permutation_real[4] = {0, 0, 0, 0};
  Bool_t chk_match_all = Parton_Jet_Match(four_gen_truth, jet_soft_coll, permutation_real);

  if(chk_match_all==kFALSE) throw LQError("Fails jet match", LQError::SkipEvent);
  FillCutFlow("JetMatchAll", weight);
  
  //ordering permutation_real
  Int_t permutation_real_temp[4];
  Int_t permutation_real_order[4];
  for(Int_t i=0; i<4; i++){ permutation_real_temp[i] = permutation_real[i]; }

  sort(permutation_real_temp, permutation_real_temp+4);

  for(Int_t i=0; i<4; i++)
    {
      for(Int_t j=0; j<4; j++)
        {
          if(permutation_real[i] == permutation_real_temp[j]) permutation_real_order[i] = j;
        }
    }
    
  /////////////////////////
  /*target jets selection*/
  /////////////////////////
  
  cout << endl;
  if(n_bjet_soft==4)
    {
      //in case of 4 b tagged jets events, just use the 4 b tagged jets
      
      for(Int_t i=0; i<n_jet_soft; i++){ if(chk_btag[i]==kTRUE) target_jet[i] = kTRUE; }

      //fit
      fitter->Set(met_vector, muon_vector, jet_soft_coll, target_jet, chk_btag);
      fitter->Fit();
    }//if(n_bjet_soft==4)
  
 else if(n_bjet_soft==3)
    {
      //in case of 3 b tagged jets events, use the 3 b tagged jets, and iterate the other non-b tagged jets.
      //then, choose jet with lowest chi^2 
      
      throw LQError("No b", LQError::SkipEvent);
    }//
 
  else if(n_bjet_soft==2)
    {
      //in case of 2 b tagged jets events, user the 2 b tagged jets, and iterate all the combitation of choosing 2 non-b tagged jets
      //then, choose 2 non b tagged jets combination with lowest chi^2

      //let's choose 2 b tagged jet first
      for(Int_t i=0; i<n_jet_soft; i++){ if(chk_btag[i]==kTRUE) target_jet[i] = kTRUE; 

      for(Int_t i=0; i<n_jet_soft; i++)
	{
	  //clear first
	  for(Int_t j=0; j<n_jet_soft; j++){ if(chk_btag[i]==kFALSE) target_jet[j] = kFALSE; }
	  
	  //let's iterate 
	  if(chk_btag[i]==kFALSE) target_jet[i] = kTRUE;
	  
	  for(Int_t j=0; j<n_jet_soft; j++)
	    {
	      if(i<j) continue;
	      
	      cout << i << "\t" << j << endl;
	      
	      //
	      if(chk_btag[j]==kFALSE) target_jet[j] = kTRUE;

	      
	    }//for loop over j
	}//for loop over i
      
      }//if(n_bjet_soft==2)
    
  delete[] chk_btag;
  delete[] target_jet;

  Bool_t chk_convergence = fitter->Get_Convergence_Checker();
  if(chk_convergence==kFALSE) throw LQError("Fitter Fail", LQError::SkipEvent);

  Double_t chi2 = fitter->Get_Chi2();

  Int_t permutation_fitter[4];
  fitter->Get_Permutation(permutation_fitter);
  
  //check jet permutation match
  Int_t jet_permutation_match = 0;
  if(permutation_real_order[0]==permutation_fitter[0] && permutation_real_order[1]==permutation_fitter[1]) jet_permutation_match = 1;
  else jet_permutation_match = 0;
    
  Double_t para_result[9];
  fitter->Get_Parameters(para_result);

  TLorentzVector fitted_object[6];
  for(Int_t i=0; i<6; i++){ fitted_object[i] = fitter->Get_Fitted_Object(i); }
  
  TLorentzVector leptonic_w = fitted_object[4] + fitted_object[5];
  Double_t leptonic_w_mass = leptonic_w.M();
    
  TLorentzVector hadronic_top = fitted_object[1] + fitted_object[2] + fitted_object[3];
  Double_t hadronic_t_mass = hadronic_top.M();

  TLorentzVector leptonic_top = fitted_object[0] + fitted_object[4] + fitted_object[5];
  Double_t leptonic_t_mass = leptonic_top.M();
  
  FillHist("MET_Truth_Vs_Measured", gen_neutrino.Pt(), met-gen_neutrino.Pt(), weight, 0, 0, 0, 0, 0, 0);
  FillHist("METPhi_Truth_Vs_Measured", gen_neutrino.Phi(), met_phi-gen_neutrino.Phi(), weight, 0, 0, 0, 0, 0, 0);
  
  FillHist("MET_Truth_Vs_Fitted", gen_neutrino.Pt(), fitted_object[5].Pt()-gen_neutrino.Pt(), weight, 0, 0, 0, 0, 0, 0);
  FillHist("METPhi_Truth_Vs_Fitted", gen_neutrino.Phi(), fitted_object[5].Phi()-gen_neutrino.Phi(), weight, 0, 0, 0, 0, 0, 0);
  FillHist("Pz_Truth_Vs_Fitted", gen_neutrino.Pz(), fitted_object[5].Pz()-gen_neutrino.Pz(), weight, 0, 0, 0, 0, 0, 0);
  FillHist("Eta_Truth_Vs_Fitted", gen_neutrino.Eta(), fitted_object[5].Eta()-gen_neutrino.Eta(), weight, 0, 0, 0, 0, 0, 0);

  FillHist("Leptonic_W_Mass", leptonic_w_mass, weight);
  FillHist("Hadronic_Top_Mass", hadronic_t_mass, weight);
  FillHist("Leptonic_Top_Mass", leptonic_t_mass, weight);
  
  Double_t chi2_piece[11];
  fitter->Get_Chi2_Piece(chi2_piece);
  
  FillHist("Chi2_Piece_0", chi2_piece[0], weight);
  FillHist("Chi2_Piece_1", chi2_piece[1], weight);
  FillHist("Chi2_Piece_2", chi2_piece[2], weight);
  FillHist("Chi2_Piece_3", chi2_piece[3], weight);
  FillHist("Chi2_Piece_4", chi2_piece[4], weight);
  FillHist("Chi2_Piece_5", chi2_piece[5], weight);
  FillHist("Chi2_Piece_6", chi2_piece[6], weight);
  FillHist("Chi2_Piece_7", chi2_piece[7], weight);
  FillHist("Chi2_Piece_8", chi2_piece[8], weight);
  FillHist("Chi2_Piece_9", chi2_piece[9], weight);
  FillHist("Chi2_Piece_10", chi2_piece[10], weight);
  
  if(n_bjet_soft==2)
    {
      FillHist("Chi2_2B", chi2, weight);
      FillHist("DiJetMass_2B", para_result[8], weight);
      FillHist("DiJetMass_Chi2_2B", para_result[8], chi2, weight, 0, 0, 0, 0, 0, 0);

      //check jet permutation match
      FillHist("Jet_Permutation_Match_2B", jet_permutation_match, weight);
      
      //goodness cut study
      for(Int_t cut_level=0; cut_level<20; cut_level++)
        {
          Double_t d_cut_level = 0.05*(cut_level+1);
          Bool_t chk_goodness_cut = fitter->Pass_Goodness_Cut(d_cut_level);

          if(chk_goodness_cut==kTRUE)
            {
              TString hname = "DiJetMass_Chi2_2B_";
              hname += (Int_t)(100*d_cut_level);

              FillHist(hname, para_result[8], chi2, weight, 0, 0, 0, 0, 0, 0);
            }
        }
    }
  else
    {
      FillHist("Chi2_3B", chi2, weight);
      FillHist("DiJetMass_3B", para_result[8], weight);
      FillHist("DiJetMass_Chi2_3B", para_result[8], chi2, weight, 0, 0, 0, 0, 0, 0);
   
      //check jet permutation match
      FillHist("Jet_Permutation_Match_3B", jet_permutation_match, weight);
    }

  FillCLHist(muhist, "Muon", muon_tight_coll, weight);
  FillCLHist(jethist, "Jet", jet_hard_coll, weight);

  return;
}//void Jet_Selection_Test_1::ExecuteEvents()

//////////

void Jet_Selection_Test_1::EndCycle() throw(LQError)
{
  m_logger << INFO << "EndCyle is called." << LQLogger::endmsg;

  return;
}//void Jet_Selection_Test_1::EndCycle()

//////////

void Jet_Selection_Test_1::InitialiseAnalysis() throw(LQError)
{
  m_logger << INFO << "Initialise Jet_Selection_Test_1 analysis." << LQLogger::endmsg;

  //Initialise histograms                                                                             
  MakeHistograms();
  
  //jet permutation efficiency
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

  MakeHistograms("DiJetMass_2B", 50, 0, 200); 
  MakeHistograms("DiJetMass_3B", 50, 0, 200);
  
  MakeHistograms("Chi2_2B", 100, 0, 50);
  MakeHistograms("Chi2_3B", 100, 0, 50);
  
  MakeHistograms("Chi2_Piece_0", 100, 0, 50);
  MakeHistograms("Chi2_Piece_1", 100, 0,50);
  MakeHistograms("Chi2_Piece_2", 100, 0,50);
  MakeHistograms("Chi2_Piece_3", 100, 0,50);
  MakeHistograms("Chi2_Piece_4", 100, 0,50);
  MakeHistograms("Chi2_Piece_5", 100, 0,50);
  MakeHistograms("Chi2_Piece_6", 100, 0,50);
  MakeHistograms("Chi2_Piece_7", 100, 0,50);
  MakeHistograms("Chi2_Piece_8", 100, 0,50);
  MakeHistograms("Chi2_Piece_9", 100, 0,50);
  MakeHistograms("Chi2_Piece_10", 100, 0,50);

  MakeHistograms("Leptonic_W_Mass", 50, 55, 105);
  MakeHistograms("Hadronic_Top_Mass", 200, 150, 190);
  MakeHistograms("Leptonic_Top_Mass", 200, 150, 190);

  MakeHistograms2D("DiJetMass_Chi2_2B", 50, 0, 200, 100, 0, 50);
  MakeHistograms2D("DiJetMass_Chi2_3B", 50, 0, 200, 100, 0, 50);

  for(Int_t cut_level=0; cut_level<20; cut_level++)
    {
      Double_t d_cut_level = 0.05*(cut_level+1);
      TString hname = "DiJetMass_Chi2_2B_";
      hname += (Int_t)(100*d_cut_level);
      
      MakeHistograms2D(hname, 50, 0, 200, 100, 0, 50);
    }

  MakeCleverHistograms(muhist, "Muon");
  MakeCleverHistograms(jethist, "Jet");

  return;
}//void Jet_Selection_Test_1::InitialiseAnalysis()

//////////

Double_t Jet_Selection_Test_1::Distance(const snu::KTruth& truth, const snu::KJet& jet)
{
  Double_t truth_phi = truth.Phi();
  Double_t truth_eta = truth.Eta();
  
  Double_t jet_phi = jet.Phi();
  Double_t jet_eta = jet.Eta();

  Double_t distance = TMath::Sqrt(TMath::Power(truth_phi-jet_phi, 2) + TMath::Power(truth_eta-jet_eta, 2));
  
  //cout << truth_phi << "\t" << truth_eta << "\t" << jet_phi << "\t" << jet_eta << "\t" << distance << endl;
  
  return distance;
}//Double_t Jet_Selection_Test_1::Distance(const snu::KTruth& truth, const snu::KJet& jet)

//////////

Bool_t Jet_Selection_Test_1::Parton_Jet_Match(const snu::KTruth gen_truth[], const vector<snu::KJet>& jet_vector, Int_t permutation_real[4])
{
  //Bool_t chk_print = kTRUE;
  Bool_t chk_print = kFALSE;
  
  if(chk_print) cout << endl;
  
  Int_t njet = jet_vector.size();

  Bool_t* match = new Bool_t[njet];
  Bool_t* chk_btag = new Bool_t[njet];
  Int_t nbjet = 0;
    for(Int_t i=0; i<njet; i++)
    {
      match[i] = kFALSE; 

      snu::KJet jet = jet_vector.at(i);
      Double_t cvs = jet.BJetTaggerValue(snu::KJet::CSVv2);
      if(cvs>CSV_THRESHOLD_MEDIUM)
	{
	  chk_btag[i] = kTRUE;
	  nbjet++;
	}
      else chk_btag[i] = kFALSE;
    }

    if(chk_print) cout << "njet = " << njet << ", nbjet = " << nbjet << endl;
    
  for(Int_t i=0; i<4; i++)
    {
      for(Int_t j=0; j<njet; j++)
	{
	  if(match[j]==kTRUE) continue;

	  snu::KJet jet = jet_vector.at(j);
	  Double_t distance = Distance(gen_truth[i], jet);
	  
	  if(chk_print) cout << i << "\t" << gen_truth[i].PdgId() << "\t" << gen_truth[i].Pt() << "\t" << j << "\t" << jet.Pt() << "\t" << distance << "\t" << chk_btag[j] << endl;
	  
	  if(distance<DISTANCE_MATCH)
	    { 
	      if(chk_print) cout << "Distance Match" << endl;
	      
	      if(chk_btag[j]==kTRUE && (gen_truth[i].PdgId()==5||gen_truth[i].PdgId()==-5) )
	      	{
	      	  match[j] = kTRUE;
	      	  permutation_real[i] = j;
		  
	      	  if(chk_print) cout << "Found b tag" << endl;
		  
	      	  break;
	      	}//b tag jet
	      else
	      	{
	      	  match[j] = kTRUE; 
	      	  permutation_real[i] = j;
		  
	      	  if(chk_print) cout << "Found No b tag" << endl;

	      	  break;
	      	}//no b tag jet
	      
	      // if((i==0||i==1)&&chk_btag[j]==kTRUE)
	      // 	{
	      // 	  match[j] = kTRUE;
	      // 	  permutation_real[i] = j;
		  
	      // 	  if(chk_print) cout << "Found" << endl;

	      // 	  break;
	      // 	}
	      // else
	      // 	{
	      // 	  match[j] = kTRUE;
	      // 	  permutation_real[i] = j;

	      // 	  if(chk_print) cout << "Found" << endl;
		  
	      // 	  break;
	      // 	}
	    
	    }//If distance
	}//for loop over leading four jet
    }//for loop over four parton
  
  Int_t count = 0;
  for(Int_t i=0; i<njet; i++){ if(match[i]==kTRUE) count++; }
  
  Bool_t chk_match = kFALSE;
  if(count==4) chk_match = kTRUE;
  
  if(chk_print)
    {
      if(chk_match==kTRUE) cout << "GOOD" << endl;
      else cout << "BAD" << endl;
    }
  
  delete[] match;
  delete[] chk_btag;

  return chk_match;
}//Bool_t Jet_Selection_Test_1::Jet_Match()

//////////
