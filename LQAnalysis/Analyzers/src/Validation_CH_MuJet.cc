#include "Validation_CH_MuJet.h"

ClassImp(Validation_CH_MuJet);

//////////

Validation_CH_MuJet::Validation_CH_MuJet() : AnalyzerCore()
{
  //to have the correct name in the log:
  SetLogName("Validation_CH_MuJet");
  
  //this function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

  m_logger << INFO << "Construct Validation_CH_MuJet" << LQLogger::endmsg;
}//Validation_CH_MuJet::Validation_CH_MuJet()

//////////

Validation_CH_MuJet::~Validation_CH_MuJet()
{
  m_logger << INFO << "Destruct Validation_CH_MuJet" << LQLogger::endmsg;
}//Validation_CH_MuJet::~Validation_CH_MuJet()

//////////

void Validation_CH_MuJet::BeginCycle() throw(LQError)
{
  m_logger << INFO << "BeginCycle is called." << LQLogger::endmsg;
  
  return;
}//void Validation_CH_MuJet::BeginCycle()

//////////

void Validation_CH_MuJet::BeginEvent() throw(LQError)
{
  m_logger << DEBUG << "BeginEvent is called." << LQLogger::endmsg;

  return;
}//void Validation_CH_MuJet::BeginEvent()

//////////

void Validation_CH_MuJet::ExecuteEvents() throw(LQError)
{
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
  
  FillCutFlow("NoCut", weight);

  Int_t runNumber = eventbase->GetEvent().RunNumber();

  /*///////*/
  /*Trigger*/
  /*///////*/
  std::vector<TString> list_triggers;
  list_triggers.push_back("HLT_IsoMu24_v");
  list_triggers.push_back("HLT_IsoTkMu24_v");
  
  //or logic
  Bool_t chk_trigger_pass = kFALSE;
  for(UInt_t i=0; i<list_triggers.size(); i++)
    {
      if(PassTrigger(list_triggers.at(i))==kTRUE)
	{
	  chk_trigger_pass = kTRUE;
	  break;
	}
    }
  if(chk_trigger_pass==kFALSE) throw LQError("Fail trigger requirement cut", LQError::SkipEvent);
  FillCutFlow("Trigger", weight);

  //Vertex cut
  if(!eventbase->GetEvent().HasGoodPrimaryVertex()) throw LQError("Fails vertex cuts", LQError::SkipEvent);;

  ///////
  /*MET*/
  ///////
  
  //Basic MET filter
  if(!PassMETFilter()) throw LQError("Fail basic MET filter.", LQError::SkipEvent);
  FillCutFlow("METFilter", weight);

  Double_t met = eventbase->GetEvent().MET(KEvent::pfmet);
  Double_t met_phi = eventbase->GetEvent().METPhi(KEvent::pfmet);

  /*Electron*/
  vector<KElectron> vector_electron_tight = GetElectrons(BaseSelection::ELECTRON_POG_TIGHT, (Float_t)LEPTON_PT_CUT, (Float_t)LEPTON_ETA_CUT);
  vector<KElectron> vector_electron_veto = GetElectrons(BaseSelection::ELECTRON_POG_VETO, (Float_t)LEPTON_PT_CUT, (Float_t)LEPTON_ETA_CUT);
   
  /*Muon*/
  vector<KMuon> vector_muon_tight = GetMuons(BaseSelection::MUON_POG_TIGHT, (Float_t)LEPTON_PT_CUT, (Float_t)LEPTON_ETA_CUT);
  vector<KMuon> vector_muon_loose = GetMuons(BaseSelection::MUON_POG_LOOSE, (Float_t)LEPTON_PT_CUT, (Float_t)LEPTON_ETA_CUT);

  /*JET*/
  vector<KJet> vector_jet_hard = GetJets("JET_NOCUT", 30, 2.4);
  vector<KJet> vector_jet_soft = GetJets("JET_NOCUT", 20, 2.4);

  //reordering by pt
  sort(vector_jet_hard.begin(), vector_jet_hard.end(), Compare_Jet_Pt);  
  sort(vector_jet_soft.begin(), vector_jet_soft.end(), Compare_Jet_Pt);
    
  /*TRUTH*/
  vector<KTruth> vector_truth = eventbase->GetTruth();

  /////////////////
  /*Pre-selection*/
  /////////////////

  //minimum met requirement
  //if(met<20) throw LQError("Fail minimum met requirement", LQError::SkipEvent);
  //FillCutFlow("Minimum_MET", weight);

  //eletron veto
  if(vector_electron_veto.size()!=0) throw LQError("Fail eletron veto", LQError::SkipEvent);
  FillCutFlow("Electron_veto", weight);

  //exact one muon
  if(vector_muon_loose.size()!=1) throw LQError("Fail only one loose muon", LQError::SkipEvent);
  FillCutFlow("One_loose_muon", weight);
  
  //muon pt>30
  if(vector_muon_loose.at(0).Pt()<30) throw LQError("Fail muon pt requirement", LQError::SkipEvent);
  FillCutFlow("Muon_Pt", weight);

  //at least four hard jet
  if(vector_jet_hard.size()<4) throw LQError("Fail at least four hard jets", LQError::SkipEvent);
  FillCutFlow("Four_hard_jets", weight);

  //at least one b-tagged jet
  Bool_t chk_b_tag[4] = {kFALSE, kFALSE, kFALSE, kFALSE};
  Int_t n_b_tag = 0;
  for(Int_t i=0; i<4; i++)
    {
      Double_t cvs = vector_jet_hard.at(i).BJetTaggerValue(snu::KJet::CSVv2);
      if(cvs>CSV_THRESHOLD_MEDIUM)
	{
	  chk_b_tag[i] = kTRUE;
	  n_b_tag++;
	}
    }
  
  if(n_b_tag<1) throw LQError("Fail at least one b-tagged jet", LQError::SkipEvent);
  FillCutFlow("B_Tag", weight);

  /////////////////
  /*Scale factors*/
  /////////////////

  //Top pair pt reweighting
  Double_t top_pair_reweight = 1;
  if(k_isdata==kFALSE)
    {
      Bool_t chk_top_found = kFALSE;
      Bool_t chk_a_top_found = kFALSE;
      Double_t top_pt;
      Double_t a_top_pt;
      for(UInt_t i=0; i<vector_truth.size(); i++)
	{
	  if(vector_truth.at(i).PdgId()==PDG_TOP && vector_truth.at(i).GenStatus()==22)
	    {
	      top_pt = vector_truth.at(i).Pt();
	      chk_top_found = kTRUE;
	    }
	  
	  if(vector_truth.at(i).PdgId()==PDG_A_TOP && vector_truth.at(i).GenStatus()==22)
	    {
	      a_top_pt = vector_truth.at(i).Pt();
	      chk_a_top_found= kTRUE;
	    }
	}

      if(chk_top_found==kTRUE && chk_a_top_found==kTRUE && top_pt<400 && a_top_pt<400) top_pair_reweight = Top_Pair_Reweight(top_pt, a_top_pt);
    }

  //Weight by trigger scale factor
  Double_t weight_by_trigger = 1;
  if(k_isdata==kFALSE) weight_by_trigger =  WeightByTrigger("HLT_IsoMu24_v", TOTAL_LUMINOSITY);
  
  //Trigger scale factor
  Double_t trigger_sf[3] = {1, 1, 1};
  if(k_isdata==kFALSE)
    {
      trigger_sf[0] = mcdata_correction->TriggerScaleFactor(vector_electron_tight, vector_muon_loose, "HLT_IsoMu24_v", -1);
      trigger_sf[1] = mcdata_correction->TriggerScaleFactor(vector_electron_tight, vector_muon_loose, "HLT_IsoMu24_v", 0);
      trigger_sf[2] = mcdata_correction->TriggerScaleFactor(vector_electron_tight, vector_muon_loose, "HLT_IsoMu24_v", 1);
    }

  //Pileup reweight
  Double_t pileup_reweight[3] = {1, 1, 1}; 
  if(MC_pu&&k_isdata==kFALSE)
    {
      pileup_reweight[0] = mcdata_correction->CatPileupWeight(eventbase->GetEvent(), -1);
      pileup_reweight[1] = mcdata_correction->CatPileupWeight(eventbase->GetEvent(), 0);
      pileup_reweight[2] = mcdata_correction->CatPileupWeight(eventbase->GetEvent(), 1);
    }
  
  //Muoon ID scale factor
  Double_t muon_id_sf[3] = {1, 1, 1};
  if(k_isdata==kFALSE)
    {
      muon_id_sf[0] = mcdata_correction->MuonScaleFactor("MUON_POG_LOOSE", vector_muon_loose, -1);
      muon_id_sf[1] = mcdata_correction->MuonScaleFactor("MUON_POG_LOOSE", vector_muon_loose, 0);
      muon_id_sf[2] = mcdata_correction->MuonScaleFactor("MUON_POG_LOOSE", vector_muon_loose, 1);
    }
  
  //Muon Iso. scale factor
  Double_t muon_iso_sf[3] = {1, 1, 1};
  if(k_isdata==kFALSE)
    {
      muon_iso_sf[0] = mcdata_correction->MuonISOScaleFactor("MUON_POG_LOOSE", vector_muon_loose, -1);
      muon_iso_sf[1] = mcdata_correction->MuonISOScaleFactor("MUON_POG_LOOSE", vector_muon_loose, 0);
      muon_iso_sf[2] = mcdata_correction->MuonISOScaleFactor("MUON_POG_LOOSE", vector_muon_loose, 1);
    }

  //Muon tracking effi. scale factor
  Double_t muon_tracking_effi_sf = 1;
  if(k_isdata==kFALSE) muon_tracking_effi_sf = mcdata_correction->MuonTrackingEffScaleFactor(vector_muon_loose);

  //B-tagging scale factor
  Int_t mc_period = GetMCPeriod();
  Double_t b_tag_reweight = 1;
  if(k_isdata==kFALSE) b_tag_reweight = BTagScaleFactor_1a(vector_jet_hard, snu::KJet::CSVv2, snu::KJet::Medium , mc_period);
  
  ///////////////////////////
  /*saving ntuple variables*/
  ///////////////////////////
  Double_t variables[18];
  
  //met
  variables[0] = met;
  
  //muon
  KMuon muon = vector_muon_loose.at(0);
  variables[1] = muon.Eta();
  variables[2] = muon.Pt();
  variables[3] = muon.RelIso04();

  //jets
  KJet jet = vector_jet_hard.at(0);
  variables[4] = jet.Eta();
  variables[5] = jet.Pt();

  jet = vector_jet_hard.at(1);
  variables[6] = jet.Eta();
  variables[7] = jet.Pt();
  
  jet = vector_jet_hard.at(2);
  variables[8] = jet.Eta();
  variables[9] = jet.Pt();

  jet = vector_jet_hard.at(3);
  variables[10] = jet.Eta();
  variables[11] = jet.Pt();
  
  //other event realated variables
  variables[12] = eventbase->GetEvent().nVertices();
  variables[13] = n_b_tag;
  
  //weight and scale factors
  variables[14] = weight;
  variables[15] = MCweight;
  variables[16] = top_pair_reweight;
  variables[17] = weight_by_trigger;
  //variables[18] = trigger_sf[0];
  //variables[19] = trigger_sf[1];
  //variables[20] = trigger_sf[2];
  //variables[18] = pileup_reweight[0];
  //variables[19] = pileup_reweight[1];
  //variables[20] = pileup_reweight[2];
  // variables[24] = muon_id_sf[0];
  // variables[25] = muon_id_sf[1];
  // variables[26] = muon_id_sf[2];
  // variables[27] = muon_iso_sf[0];
  // variables[28] = muon_iso_sf[1];
  // variables[29] = muon_iso_sf[2];
  // variables[30] = muon_tracking_effi_sf; 
  // variables[31] = b_tag_reweight;

  //FillNtp("tuple_variables", variables);

  return;
}//void Validation_CH_MuJet::ExecuteEvents()

//////////

void Validation_CH_MuJet::EndCycle() throw(LQError)
{
  m_logger << INFO << "EndCyle is called." << LQLogger::endmsg;

  return;
}//void Validation_CH_MuJet::EndCycle()

//////////

void Validation_CH_MuJet::InitialiseAnalysis() throw(LQError)
{
  m_logger << INFO << "Initialise Validation_CH_MuJet analysis" << LQLogger::endmsg;

  //Initialise histograms
  MakeHistograms();

  MakeNtp("tuple_variables", "met:muon_eta:muon_pt:muon_reliso04:jet0_eta:jet0_pt:jet1_eta:jet1_pt:jet2_eta:jet2_pt:jet3_eta:jet3_pt:n_vertices:n_b_tag:weight:mc_weight:top_pair_reweight:weight_by_trigger");//:pileup_reweight_down");//:pileup_reweight");//:pileup_reweight_up");

//trigger_sf_down:trigger_sf:trigger_sf_up");//:pileup_reweight_down:pileup_reweight:pileup_reweight_up:muon_id_sf_down:muon_id_sf:muon_id_sf_up:muon_iso_sf_down:muon_iso_sf:muon_iso_sf_up:muon_tracking_eff_sf:b_tag_reweight");

  return;
}//void Validation_CH_MuJet::InitialiseAnalysis()

//////////
