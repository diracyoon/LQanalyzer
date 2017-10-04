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

  //eletron veto
  if(vector_electron_veto.size()!=0) throw LQError("Fail eletron veto", LQError::SkipEvent);
  FillCutFlow("Electron_veto", weight);

  //exact one muon
  if(vector_muon_tight.size()!=1 || vector_muon_loose.size()>1) throw LQError("Fail only one tight muon", LQError::SkipEvent);
  FillCutFlow("One_tight_muon", weight);
  
  //muon pt>20
  if(vector_muon_tight.at(0).Pt()<20) throw LQError("Fail muon pt requirement", LQError::SkipEvent);
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

  cout << MCweight << endl;

  /////////////////
  /*Scale factors*/
  /////////////////

  //Top pair pt reweighting
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

  Double_t top_pair_reweight = 1;
  if(k_isdata==kFALSE && chk_top_found==kTRUE && chk_a_top_found==kTRUE && top_pt<400 && a_top_pt<400) top_pair_reweight = Top_Pair_Reweight(top_pt, a_top_pt);
  
  //Trigger scale factor
  Double_t trigger_sf[3] = { 1, 1, 1};
  if(k_isdata==kFALSE)
    {
      trigger_sf[0] = mcdata_correction->TriggerScaleFactor(vector_electron_tight, vector_muon_tight, "HLT_IsoMu24_v", -1);
      trigger_sf[1] = mcdata_correction->TriggerScaleFactor(vector_electron_tight, vector_muon_tight, "HLT_IsoMu24_v", 0);
      trigger_sf[2] = mcdata_correction->TriggerScaleFactor(vector_electron_tight, vector_muon_tight, "HLT_IsoMu24_v", 1);
    }

  //Weight by trigger scale factor
  Double_t weight_by_trigger_sf = 1;
  if(k_isdata==kFALSE) weight_by_trigger_sf =  WeightByTrigger("HLT_IsoMu24_v", TOTAL_LUMINOSITY);

  //Pile up reweight
  
  //Muoon ID scale factor
  Double_t muon_id_sf[3] = {1, 1, 1};
  if(k_isdata==kFALSE)
    {
      muon_id_sf[0] = mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", vector_muon_tight, -1);
      muon_id_sf[1] = mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", vector_muon_tight, 0);
      muon_id_sf[2] = mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", vector_muon_tight, 1);
    }
  
  //Muon Iso. scale factor
  Double_t muon_iso_sf[3] = {1, 1, 1};
  if(k_isdata==kFALSE)
    {
      muon_iso_sf[0] = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", vector_muon_tight, -1);
      muon_iso_sf[1] = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", vector_muon_tight, 0);
      muon_iso_sf[2] = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", vector_muon_tight, 1);
    }

  //Muon tracking effi. scale factor
  Double_t muon_tracking_effi_sf = 1;
  if(k_isdata==kFALSE) muon_tracking_effi_sf = mcdata_correction->MuonTrackingEffScaleFactor(vector_muon_tight);

  ///////////////////////////
  /*saving ntuple variables*/
  ///////////////////////////
  Double_t variables[29];
  
  //met
  variables[0] = met;
  
  //muon
  KMuon muon = vector_muon_tight.at(0);
  variables[1] = muon.Eta();
  variables[2] = muon.Pt();
  
  //jets
  for(Int_t i=0; i<4; i++)
    {
      KJet jet = vector_jet_hard.at(i);
      variables[2*i+3] = jet.Eta();
      variables[2*i+4] = jet.Pt();
    }
  
  //
  variables[11] = eventbase->GetEvent().nVertices();
  variables[12] = n_b_tag;

  //weight and scale factors
  variables[13] = weight;
  variables[14] = top_pair_reweight;
  variables[15] = trigger_sf[0];
  variables[16] = trigger_sf[1];
  variables[17] = trigger_sf[2];
  variables[18] = weight_by_trigger_sf;
  variables[19] = 1;
  variables[20] = 1;
  variables[21] = 1;
  variables[22] = muon_id_sf[0];
  variables[23] = muon_id_sf[1];
  variables[24] = muon_id_sf[2];
  variables[25] = muon_iso_sf[0];
  variables[26] = muon_iso_sf[1];
  variables[27] = muon_iso_sf[2];
  variables[28] = muon_tracking_effi_sf; 

  FillNtp("tuple_variables", variables);

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

  MakeNtp("tuple_variables", "met:muon_eta:muon_pt:jet0_eta:jet0_pt:jet1_eta:jet1_pt:jet2_eta:jet2_pt:jet3_eta:jet3_pt:n_vertices:n_b_tag:weight:top_pair_reweight:trigger_sf_down:trigger_sf:trigger_sf_up:weight_by_trigger:muon_id_sf_down:muon_id_sf:muon_id_sf_up:muon_iso_sf_down:muon_iso_sf:muon_iso_sf_up:muon_tracking_eff_sf");

  return;
}//void Validation_CH_MuJet::InitialiseAnalysis()

//////////
