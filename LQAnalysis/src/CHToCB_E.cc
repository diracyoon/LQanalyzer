//Core includes
#include "EventBase.h" 
#include "BaseSelection.h"

//Local includes
#include "CHToCB_E.h"

//Needed to allow inheritance for use in LQCore/core classes
ClassImp(CHToCB_E);
 
//////////

CHToCB_E::CHToCB_E() : AnalyzerCore(), out_electrons(0), out_muons(0)
{
  //To have the correct name in the log:                                                          
  SetLogName("CHToCB_E");
  
  m_logger << INFO << "Construct CHToCB_E." << LQLogger::endmsg;

  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

  //MakeCleverHistograms(sighist_mm, "DiMuon");
}//CHToCB_E::CHToCB_E()

//////////

CHToCB_E::~CHToCB_E() 
{
  m_logger << INFO << "Destruct CHToCB_E." << LQLogger::endmsg;
}//CHToCB_E::~CHToCB_E()

//////////

void CHToCB_E::BeginCycle() throw( LQError )
{
  m_logger << INFO << "BeginCycle is called." << LQLogger::endmsg;

  return;
}//void CHToCB_E::BeginCycle()

//////////

void CHToCB_E::BeginEvent() throw(LQError)
{
  m_logger << DEBUG << "BeginEvent is called." << LQLogger::endmsg;
  
  return;
}//void CHToCB_E::BeginEvent()

//////////

void CHToCB_E::ClearOutputVectors() throw(LQError) 
{
  out_muons.clear();
  out_electrons.clear();

  return;
}//void CHToCB_E::ClearOutputVectors()

//////////

void CHToCB_E::ExecuteEvents() throw(LQError)
{
  

  m_logger << DEBUG << "RunNumber = " << eventbase->GetEvent().RunNumber() << ", EventNumber = " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "IsData = " << isData << LQLogger::endmsg;

  //Apply the gen weight 
  if(!isData) weight*=MCweight;

  /*////////////////////*/
  /*////////////////////*/
  /*Event Selection*/
  /*////////////////////*/
  /*////////////////////*/
  
  FillCutFlow("NoCut", weight);

  /*Event*/
  //Trigger 
  std::vector<TString> list_triggers;
  list_triggers.push_back("HLT_Ele32_eta2p1_WPTight_Gsf_v");//HLT_Ele27_eta2p1_WPLoose_Gsf_v");
  for(unsigned int i=0; i<list_triggers.size(); i++)
    {
      TString trigger = list_triggers.at(i);
      // cout << "test 1" << endl;
      // std::vector<std::string> trigger_in_data = eventbase->GetTrigger().GetHLTInsideDatasetTriggerNames();
      // for(unsigned int j=0; j<trigger_in_data.size(); j++)
      // 	{
      // 	  cout << trigger_in_data.at(j) << endl;
      // 	}
      // cout << endl;
      if(!PassTrigger(trigger)) throw LQError("Fails Trigger requirement cuts", LQError::SkipEvent);
    }
  FillCutFlow("Trigger requirement", weight);
   
  //Vertex cut
  if(!eventbase->GetEvent().HasGoodPrimaryVertex()) throw LQError("Fails vertex cuts", LQError::SkipEvent);
  FillCutFlow("VertexCut", weight);

  /*MET*/
  //Pass basic MET filter
  if(!PassMETFilter()) throw LQError("Fail basic MET filter.", LQError::SkipEvent);
  FillCutFlow("METCut", weight);

  //Mininum MET requirement
  double met = eventbase->GetEvent().PFMET();
  if(met<20) throw LQError("Fail minimum MET requirement", LQError::SkipEvent);
  FillCutFlow("MininumMET", weight);
  
  //double met_phi = eventbase->GetEvent().PFMETphi();

  /*Electron*/
  std::vector<snu::KElectron> electron_tight_coll;
  eventbase->GetElectronSel()->SelectElectrons(electron_tight_coll, BaseSelection::ELECTRON_TOP_TIGHT);

  std::vector<snu::KElectron> electron_veto_coll;
  eventbase->GetElectronSel()->SelectElectrons(electron_veto_coll, BaseSelection::ELECTRON_TOP_VETO);
  
  //tight electron==1 cut
  if(electron_tight_coll.size()!=1) throw LQError("Fails tight electron==1 cuts", LQError::SkipEvent);
  FillCutFlow("OneTightElectron", weight);

  //extra electron veto cut
  if(electron_veto_coll.size()>1) throw LQError("Fails extra electron veto cuts", LQError::SkipEvent);
  FillCutFlow("ExtraElectronVeto", weight);

  /*Muon*/
  std::vector<snu::KMuon> muon_tight_coll;
  eventbase->GetMuonSel()->SelectMuons(muon_tight_coll, BaseSelection::MUON_TOP_TIGHT);

  std::vector<snu::KMuon> muon_veto_coll;
  eventbase->GetMuonSel()->SelectMuons(muon_veto_coll, BaseSelection::MUON_TOP_VETO);
  
  //muon veto cuts
  if(muon_veto_coll.size()>0) throw LQError("Fails muon veto cuts",  LQError::SkipEvent);
  FillCutFlow("MuonVeto", weight);

  /*Jet*/
  std::vector<snu::KJet> jet_hard_coll;
  eventbase->GetJetSel()->SelectJets(isData, jet_hard_coll, muon_tight_coll, electron_tight_coll, "JET_NOCUT", 30., 2.5);//What is a role of isData, in here?

  std::vector<snu::KJet> jet_soft_coll;
  eventbase->GetJetSel()->SelectJets(isData, jet_soft_coll, muon_tight_coll, electron_tight_coll, "JET_NOCUT", 10., 2.5);

  //at least four hard jets cut
  if(jet_hard_coll.size()<4) throw LQError("Fails at least four hard jets cuts.", LQError::SkipEvent);
  FillCutFlow("FourJets", weight);
  
  //reordering jet by pt
  sort(jet_hard_coll.begin(), jet_hard_coll.end(), Compare_Jet_Pt);
  sort(jet_soft_coll.begin(), jet_soft_coll.end(), Compare_Jet_Pt);

  //at leat two b jets cut
  int nbjet_hard = 0;
  int njet_hard = jet_hard_coll.size();
  for(int i=0; i<njet_hard; i++)
    {
      double cvs = jet_hard_coll.at(i).BJetTaggerValue(snu::KJet::CSVv2);
      
      if(cvs>CSV_THRESHOLD) nbjet_hard++;
    }

  if(nbjet_hard<2) throw LQError("Fails at least four hard jets cuts.", LQError::SkipEvent);
  FillCutFlow("TwoBJets", weight);

  /*////////////////////*/
  /*////////////////////*/
  /*Kinematic Fitter*/
  /*////////////////////*/
  /*////////////////////*/

  //contructing Electron four vector
  TLorentzVector electron_vector;
  //electron_vector.SetPtEtaPhiE();
  
  //kinematic_fitter.Set();
  //kinematic_fitter.Fit();
  
  return;
}//void CHToCB_E::ExecuteEvents()

//////////

void CHToCB_E::EndCycle() throw(LQError)
{
  m_logger << INFO << "EndCyle is called." << LQLogger::endmsg;

  return;
}//void CHToCB_E::EndCycle()

//////////

void CHToCB_E::EndEvent() throw(LQError)
{
  m_logger << DEBUG << "EndEvent is called." << LQLogger::endmsg;

  return;
}//void CHToCB_E::EndEvent() throw(LQError)

//////////

void CHToCB_E::InitialiseAnalysis() throw(LQError)
{
  m_logger << INFO << "Initialise CHToCB_E analysis." << LQLogger::endmsg;

  //Initialise histograms                                                                             
  MakeHistograms();

  return;
}//void CHToCB_E::InitialiseAnalysis()

//////////
