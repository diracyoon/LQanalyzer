//Core includes
#include "EventBase.h" 
#include "BaseSelection.h"

//Local includes
#include "CHToCB_M.h"

//Needed to allow inheritance for use in LQCore/core classes
ClassImp(CHToCB_M);
 
//////////

CHToCB_M::CHToCB_M() : AnalyzerCore()   
{
  //To have the correct name in the log:                                                          
  SetLogName("CHToCB_M");
  
  m_logger << INFO << "Construct CHToCB_M." << LQLogger::endmsg;

  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

  chk_debug = kFALSE;

  fitter = new Kinematic_Fitter_Old(chk_debug);
}//CHToCB_M::CHToCB_M()

//////////

CHToCB_M::~CHToCB_M() 
{
  m_logger << INFO << "Destruct CHToCB_M." << LQLogger::endmsg;

  delete fitter;
}//CHToCB_M::~CHToCB_M()

//////////

void CHToCB_M::BeginCycle() throw(LQError)
{
  m_logger << INFO << "BeginCycle is called." << LQLogger::endmsg;

  return;
}//void CHToCB_M::BeginCycle()

//////////

void CHToCB_M::BeginEvent() throw(LQError)
{
  m_logger << DEBUG << "BeginEvent is called." << LQLogger::endmsg;

  return;
}//void CHToCB_M::BeginEvent()

//////////

void CHToCB_M::ExecuteEvents() throw(LQError)
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
  list_triggers.push_back("HLT_IsoMu24_v");

  for(unsigned int i=0; i<list_triggers.size(); i++)
    {
      TString trigger = list_triggers.at(i);
      // std::vector<std::string> trigger_in_data = eventbase->GetTrigger().GetHLTInsideDatasetTriggerNames();                                                                                               
      // for(unsigned int j=0; j<trigger_in_data.size(); j++)                                          
      //  	{                                                                                      
      // 	  cout << trigger_in_data.at(j) << endl;                                               
      //  	}                                                                                      
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
  Double_t met = eventbase->GetEvent().MET();
  Double_t met_phi = eventbase->GetEvent().METPhi();
  
  if(met<20) throw LQError("Fail minimum MET requirement", LQError::SkipEvent);
  FillCutFlow("MininumMET", weight);

  //construct met vector
  TLorentzVector met_vector;
  met_vector.SetPxPyPzE(met*TMath::Cos(met_phi), met*TMath::Sin(met_phi), 0, met);
  
  /*Electron*/
  std::vector<snu::KElectron> electron_tight_coll;
  eventbase->GetElectronSel()->SelectElectrons(electron_tight_coll, BaseSelection::ELECTRON_POG_TIGHT);
  

  std::vector<snu::KElectron> electron_veto_coll;
  eventbase->GetElectronSel()->SelectElectrons(electron_veto_coll, BaseSelection::ELECTRON_POG_VETO);

  //electron veto cuts
  if(electron_veto_coll.size()>0) throw LQError("Fails electron veto cuts",  LQError::SkipEvent);
  FillCutFlow("ElectronVeto", weight);

  /*Muon*/
  std::vector<snu::KMuon> muon_tight_coll;
  eventbase->GetMuonSel()->SelectMuons(muon_tight_coll, BaseSelection::MUON_POG_TIGHT);

  std::vector<snu::KMuon> muon_veto_coll;
  eventbase->GetMuonSel()->SelectMuons(muon_veto_coll, BaseSelection::MUON_POG_LOOSE);
  
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
  eventbase->GetJetSel()->SelectJets(isData, jet_hard_coll, muon_tight_coll, electron_tight_coll, "JET_NOCUT", 30., 2.5);

  std::vector<snu::KJet> jet_soft_coll;  
  eventbase->GetJetSel()->SelectJets(isData, jet_soft_coll, muon_tight_coll, electron_tight_coll, "JET_NOCUT", 20., 2.5);

  //reordering jet by pt
  sort(jet_hard_coll.begin(), jet_hard_coll.end(), Compare_Jet_Pt);
  sort(jet_soft_coll.begin(), jet_soft_coll.end(), Compare_Jet_Pt);
  
  //at least four hard jets cut
  if(jet_hard_coll.size()<4) throw LQError("Fails at least four hard jets cuts.", LQError::SkipEvent);
  FillCutFlow("FourJets", weight);

  //at leat two b jets cut
  Bool_t chk_btag[4] = {kFALSE};
  Int_t nbjet_hard = 0;
  for(Int_t i=0; i<4; i++)
    {
      Double_t cvs = jet_hard_coll.at(i).BJetTaggerValue(snu::KJet::CSVv2);

      if(cvs>CSV_THRESHOLD)
	{
	  chk_btag[i] = kTRUE;
	  nbjet_hard++;
	}
    }

  if(nbjet_hard<2) throw LQError("Fails at least four hard jets cuts.", LQError::SkipEvent);
  FillCutFlow("TwoBJets", weight);

  //construct jet vector for easy handling
  vector<TLorentzVector> jet_vector;
  Int_t njet_soft = jet_soft_coll.size();
  for(Int_t i=0; i<njet_soft; i++){ jet_vector.push_back(jet_soft_coll.at(i)); }

  
  if(chk_debug) cout << "START " << eventbase->GetEvent().EventNumber() << endl;
  
  fitter->Set(met_vector, muon_vector, jet_vector, chk_btag);
  fitter->Fit();  
  
  Bool_t chk_convergence = fitter->Get_Convergence();
  if(chk_convergence==kFALSE) throw LQError("Fitter Fail", LQError::SkipEvent);  

  Double_t chi2 = fitter->Get_Chi2();
  
  Int_t permutation[4];
  fitter->Get_Permutation(permutation);
  
  Double_t para_result[9];
  fitter->Get_Parameters(para_result);
  
  Double_t t_mass = 0;//fitter->Get_Top_Mass();

  if(chk_debug)
    {
      cout << "Result " << eventbase->GetEvent().EventNumber() << endl;
      cout << "nbtag= " <<  nbjet_hard << ", Mass=" << para_result[8] << " " << permutation[0] << " " << permutation[1] << " " << permutation[2] << " " << permutation[3] << endl << endl;
    }
  
  if(nbjet_hard==2)
    {
      FillHist("Chi2_2B", chi2, weight);
      FillHist("DiJetMass_2B", para_result[8], weight);
      FillHist("DiJetMass_Chi2_2B", para_result[8], chi2, weight, 0, 0, 0, 0, 0, 0);

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
   
      FillHist("Hadronic_Top_Mass", t_mass, weight);
    }
  else
    {
      FillHist("Chi2_3B", chi2, weight);
      FillHist("DiJetMass_3B", para_result[8], weight);
      FillHist("DiJetMass_Chi2_3B", para_result[8], chi2, weight, 0, 0, 0, 0, 0, 0);
    }
  
  FillCLHist(muhist, "Muon", muon_tight_coll, weight);
  FillCLHist(jethist, "Jet", jet_hard_coll, weight);
  
  return;
}//void CHToCB_M::ExecuteEvents()

//////////

void CHToCB_M::EndCycle() throw(LQError)
{
  m_logger << INFO << "EndCyle is called." << LQLogger::endmsg;

  return;
}//void CHToCB_M::EndCycle()

//////////

void CHToCB_M::InitialiseAnalysis() throw(LQError)
{
  m_logger << INFO << "Initialise CHToCB_E analysis." << LQLogger::endmsg;

  //Initialise histograms                                                                             
  MakeHistograms();
  MakeHistograms("DiJetMass_2B", 50, 0, 200); 
  MakeHistograms("DiJetMass_3B", 50, 0, 200);
  
  MakeHistograms("Chi2_2B", 25, 0, 50);
  MakeHistograms("Chi2_3B", 25, 0, 50);
  
  MakeHistograms("Hadronic_Top_Mass", 200, 120, 220);

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
}//void CHToCB_M::InitialiseAnalysis()

//////////
