//Core includes
#include "EventBase.h" 
#include "BaseSelection.h"

//Local includes
#include "Jet_Inclusion_Prob.h"

//Needed to allow inheritance for use in LQCore/core classes
ClassImp(Jet_Inclusion_Prob);
 
//////////

Jet_Inclusion_Prob::Jet_Inclusion_Prob() : AnalyzerCore()   
{
  //To have the correct name in the log:                                                          
  SetLogName("Jet_Inclusion_Prob");
  
  m_logger << INFO << "Construct Jet_Inclusion_Prob." << LQLogger::endmsg;

  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

  chk_debug = kFALSE;

  fitter = new Kinematic_Fitter(chk_debug);
  //fitter = new Kinematic_Fitter_Old(chk_debug);
}//Jet_Inclusion_Prob::Jet_Inclusion_Prob()

//////////

Jet_Inclusion_Prob::~Jet_Inclusion_Prob() 
{
  m_logger << INFO << "Destruct Jet_Inclusion_Prob." << LQLogger::endmsg;

  delete fitter;
}//Jet_Inclusion_Prob::~Jet_Inclusion_Prob()

//////////

void Jet_Inclusion_Prob::BeginCycle() throw(LQError)
{
  m_logger << INFO << "BeginCycle is called." << LQLogger::endmsg;

  return;
}//void Jet_Inclusion_Prob::BeginCycle()

//////////

void Jet_Inclusion_Prob::BeginEvent() throw(LQError)
{
  m_logger << DEBUG << "BeginEvent is called." << LQLogger::endmsg;

  return;
}//void Jet_Inclusion_Prob::BeginEvent()

//////////

void Jet_Inclusion_Prob::ExecuteEvents() throw(LQError)
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

  //at least four hard jets cut
  if(jet_hard_coll.size()<4) throw LQError("Fails at least four hard jets cuts.", LQError::SkipEvent);
  FillCutFlow("FourJets", weight);

  Int_t n_jet_soft = jet_soft_coll.size();
  Int_t n_bjet_soft = 0;
  for(Int_t i=0; i<n_jet_soft; i++)
    {
      Double_t csv = jet_soft_coll.at(i).BJetTaggerValue(snu::KJet::CSVv2);
      
      if(csv>CSV_THRESHOLD_MEDIUM) n_bjet_soft++;
    }
  
  //at least two b tagged jet cut
  if(n_bjet_soft<2) throw LQError("Fails at least two b tagged jets cuts.", LQError::SkipEvent);
  FillCutFlow("TwoBJets", weight);
  
  //reordering jet by pt
  sort(jet_soft_coll.begin(), jet_soft_coll.end(), Compare_Jet_Pt);
 
  //parton-jet matching for all jets
  Int_t permutation_pt[4] = {0, 0, 0, 0};
  Bool_t chk_match_all = Truth_Jet_Match(gen_quark, jet_soft_coll, permutation_pt);

  if(chk_match_all==kFALSE) throw LQError("Fails jet match", LQError::SkipEvent);
  FillCutFlow("JetMatchAll", weight);
   
  
  Bool_t chk_inclusion_pt_leading_4_jet = kTRUE;
  Bool_t chk_inclusion_pt_leading_5_jet = kTRUE;
  Bool_t chk_inclusion_pt_leading_6_jet = kTRUE;
  Bool_t chk_inclusion_pt_leading_7_jet = kTRUE;
  Bool_t chk_inclusion_pt_leading_8_jet = kTRUE;
  
  for(Int_t i=0; i<4; i++)
    {
      if(permutation_pt[i]>3) chk_inclusion_pt_leading_4_jet = kFALSE;
      if(permutation_pt[i]>4) chk_inclusion_pt_leading_5_jet = kFALSE;
      if(permutation_pt[i]>5) chk_inclusion_pt_leading_6_jet = kFALSE;
      if(permutation_pt[i]>6) chk_inclusion_pt_leading_7_jet = kFALSE;
      if(permutation_pt[i]>7) chk_inclusion_pt_leading_8_jet = kFALSE;
    }
  
  //pt leading 4 jets
  if(chk_inclusion_pt_leading_4_jet==kTRUE) FillHist("Jets_Inclusion_Prob_Pt_Leading_4_Jets", 1, weight);
  else FillHist("Jets_Inclusion_Prob_Pt_Leading_4_Jets", 0, weight);
    
  //pt leading 5 jets
  if(chk_inclusion_pt_leading_5_jet==kTRUE) FillHist("Jets_Inclusion_Prob_Pt_Leading_5_Jets", 1, weight);
  else FillHist("Jets_Inclusion_Prob_Pt_Leading_5_Jets", 0, weight);

  //pt leading 6 jets
  if(chk_inclusion_pt_leading_6_jet==kTRUE) FillHist("Jets_Inclusion_Prob_Pt_Leading_6_Jets", 1, weight);
  else FillHist("Jets_Inclusion_Prob_Pt_Leading_6_Jets", 0, weight);
  
  //pt leading 7 jets
  if(chk_inclusion_pt_leading_7_jet==kTRUE) FillHist("Jets_Inclusion_Prob_Pt_Leading_7_Jets", 1, weight);
  else FillHist("Jets_Inclusion_Prob_Pt_Leading_7_Jets", 0, weight);
  
  //pt leading 8 jets
  if(chk_inclusion_pt_leading_8_jet==kTRUE) FillHist("Jets_Inclusion_Prob_Pt_Leading_8_Jets", 1, weight);
  else FillHist("Jets_Inclusion_Prob_Pt_Leading_8_Jets", 0, weight);

  //reorderind jet by csv
  sort(jet_soft_coll.begin(), jet_soft_coll.end(), Compare_Jet_CSV);
  
  //parton-jet matching
  Int_t permutation_csv[4] = {0, 0, 0, 0};
  Truth_Jet_Match(gen_quark, jet_soft_coll, permutation_csv);
  
  Bool_t chk_inclusion_csv_leading_4_jet = kTRUE;
  Bool_t chk_inclusion_csv_leading_5_jet = kTRUE;
  Bool_t chk_inclusion_csv_leading_6_jet = kTRUE;
  Bool_t chk_inclusion_csv_leading_7_jet = kTRUE;
  Bool_t chk_inclusion_csv_leading_8_jet = kTRUE;
  for(Int_t i=0; i<4; i++)
    {
      if(permutation_csv[i]>3) chk_inclusion_csv_leading_4_jet = kFALSE;
      if(permutation_csv[i]>4) chk_inclusion_csv_leading_5_jet = kFALSE;
      if(permutation_csv[i]>5) chk_inclusion_csv_leading_6_jet = kFALSE;
      if(permutation_csv[i]>6) chk_inclusion_csv_leading_7_jet = kFALSE;
      if(permutation_csv[i]>7) chk_inclusion_csv_leading_8_jet = kFALSE;
    }

  //csv leading 4 jets
  if(chk_inclusion_csv_leading_4_jet==kTRUE) FillHist("Jets_Inclusion_Prob_CSV_Leading_4_Jets", 1, weight);
  else FillHist("Jets_Inclusion_Prob_CSV_Leading_4_Jets", 0, weight);

  //csv leading 5 jets
  if(chk_inclusion_csv_leading_5_jet==kTRUE) FillHist("Jets_Inclusion_Prob_CSV_Leading_5_Jets", 1, weight);
  else FillHist("Jets_Inclusion_Prob_CSV_Leading_5_Jets", 0, weight);

  //csv leading 6 jets
  if(chk_inclusion_csv_leading_6_jet==kTRUE) FillHist("Jets_Inclusion_Prob_CSV_Leading_6_Jets", 1, weight);
  else FillHist("Jets_Inclusion_Prob_CSV_Leading_6_Jets", 0, weight);

  //csv leading 7 jets
  if(chk_inclusion_csv_leading_7_jet==kTRUE) FillHist("Jets_Inclusion_Prob_CSV_Leading_7_Jets", 1, weight);
  else FillHist("Jets_Inclusion_Prob_CSV_Leading_7_Jets", 0, weight);

  //csv leading 8 jets
  if(chk_inclusion_csv_leading_8_jet==kTRUE) FillHist("Jets_Inclusion_Prob_CSV_Leading_8_Jets", 1, weight);
  else FillHist("Jets_Inclusion_Prob_CSV_Leading_8_Jets", 0, weight);
  
  return;
}//void Jet_Inclusion_Prob::ExecuteEvents()

//////////

void Jet_Inclusion_Prob::EndCycle() throw(LQError)
{
  m_logger << INFO << "EndCyle is called." << LQLogger::endmsg;

  return;
}//void Jet_Inclusion_Prob::EndCycle()

//////////

void Jet_Inclusion_Prob::InitialiseAnalysis() throw(LQError)
{
  m_logger << INFO << "Initialise Jet_Inclusion_Prob analysis." << LQLogger::endmsg;

  //Initialise histograms                                                                             
  MakeHistograms();
  
  //jet inclusion probability for pt leading jets
  MakeHistograms("Jets_Inclusion_Prob_Pt_Leading_4_Jets", 2, -0.5, 1.5);
  MakeHistograms("Jets_Inclusion_Prob_Pt_Leading_5_Jets", 2, -0.5, 1.5);
  MakeHistograms("Jets_Inclusion_Prob_Pt_Leading_6_Jets", 2, -0.5, 1.5);
  MakeHistograms("Jets_Inclusion_Prob_Pt_Leading_7_Jets", 2, -0.5, 1.5);
  MakeHistograms("Jets_Inclusion_Prob_Pt_Leading_8_Jets", 2, -0.5, 1.5);

  //jet inclusion probability for csv leading jets
  MakeHistograms("Jets_Inclusion_Prob_CSV_Leading_4_Jets", 2, -0.5, 1.5);
  MakeHistograms("Jets_Inclusion_Prob_CSV_Leading_5_Jets", 2, -0.5, 1.5);
  MakeHistograms("Jets_Inclusion_Prob_CSV_Leading_6_Jets", 2, -0.5, 1.5);
  MakeHistograms("Jets_Inclusion_Prob_CSV_Leading_7_Jets", 2, -0.5, 1.5);
  MakeHistograms("Jets_Inclusion_Prob_CSV_Leading_8_Jets", 2, -0.5, 1.5);

  return;
}//void Jet_Inclusion_Prob::InitialiseAnalysis()

//////////
