// $Id: ExampleAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNDiLepton Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *d
 ***************************************************************************/

/// Local includes
#include "HNDiLepton.h"


//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to samevallow inheritance for use in LQCore/core classes
ClassImp (HNDiLepton);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
HNDiLepton::HNDiLepton() :  AnalyzerCore(),  out_electrons(0) {
  
  mapcounter.clear();

  _mm_channel=false;
  _m_channel=false;
  _ee_channel=false;
  _e_channel=false;
  _em_channel=false;

  _mm_met_presel_cut=80.;
  _ee_met_presel_cut=80.;
  _em_met_presel_cut=80.;

  _mm_mll_presel_cut=10.;
  _ee_mll_presel_cut=10.;
  _em_mll_presel_cut=10.;

  _m_tightid = "MUON_HN_TIGHT";
  _e_tightid = "ELECTRON_HN_TIGHTv4";
  _m_looseid = "MUON_HN_LOOSEv7_SIP3";
  _e_looseid = "ELECTRON_HN_FAKELOOSEv7";



  configure_fakes=true;

  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNDiLepton");

  Message("In HNDiLepton constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  //InitialiseAnalysis();
  
}


void HNDiLepton::InitialiseAnalysis() throw( LQError ) {
  
  if(k_flags.size() > 0){

    if(k_flags[0] == "CutFlow")  functionality = HNDiLepton::CUTFLOW;
    else if(k_flags[0] == "OPT")  functionality = HNDiLepton::OPT;
    else functionality = HNDiLepton::ANALYSIS;
  }

  
  if(  functionality == HNDiLepton::ANALYSIS){
    /// Initialise histograms
    MakeHistograms();  

    //
    // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
    // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
    //
    
    Message("Making clever hists for Z ->ll test code", INFO);
    
    vector<TString> labels;
    labels.push_back("DiMuon");
    labels.push_back("DiMuon_isomu");
    labels.push_back("DiMuon_allFJ");
    //labels.push_back("PogDiMuon");
    
    for(unsigned int il = 0; il < labels.size(); il++){
      TString label = labels.at(il);
      MakeCleverHistograms(sighist_mmm,label+"WZMuon");
      MakeCleverHistograms(sighist_mmm,label+"ZZMuon1");
      MakeCleverHistograms(sighist_mmm,label+"ZZMuon2");
      MakeCleverHistograms(sighist_mm, label + "SSMM_nolepveto");
      MakeCleverHistograms(sighist_mm, label+"SSMM_Preselection");
      MakeCleverHistograms(sighist_mm, label+"SSMM_Preselection_aj20");
      MakeCleverHistograms(sighist_mm, label+"SSMM_Preselection_aj60");
      MakeCleverHistograms(sighist_mm, label+"SSMM_Preselection_sf");
      MakeCleverHistograms(sighist_mm, label+"SSMM_Preselection_df");
      MakeCleverHistograms(sighist_mm, label+"SSMM_1JetCR");
      MakeCleverHistograms(sighist_mm, label+"SSMM_1JetOpt");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_fix");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_veto");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_muonveto");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_noW");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_xc");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_xc2");
      MakeCleverHistograms(sighist_mm, label+"SSMM_0Jet");
      MakeCleverHistograms(sighist_mm, label+"SSMM_1Jet");
      MakeCleverHistograms(sighist_mm, label+"SSMM_2Jet");
      MakeCleverHistograms(sighist_mm, label+"SSMM_3Jet");
      MakeCleverHistograms(sighist_mm, label+"SSMM_4Jet");
      MakeCleverHistograms(sighist_mm, label+"SSMM_0Jet_xc");
      MakeCleverHistograms(sighist_mm, label+"SSMM_1Jet_xc");
      MakeCleverHistograms(sighist_mm, label+"SSMM_2Jet_xc");
      MakeCleverHistograms(sighist_mm, label+"SSMM_2Jet_rmbin");
      MakeCleverHistograms(sighist_mm, label+"SSMM_2Jet_xc_lt");
      MakeCleverHistograms(sighist_mm, label+"SSMM_2Jet_xc_met");
      MakeCleverHistograms(sighist_mm, label+"SSMM_2Jet_xc_met_cfj");
      MakeCleverHistograms(sighist_mm, label+"SSMM_2Jet_xc_met_cfj_l1mass");
      MakeCleverHistograms(sighist_mm, label+"SSMM_2Jet_xc_met_cfj_l1masshigh");
      MakeCleverHistograms(sighist_mm, label+"SSMM_2Jet_xc_minus");
      MakeCleverHistograms(sighist_mm, label+"SSMM_2Jet_xc_plus");

      MakeCleverHistograms(sighist_mm, label+"SSMM_DiTJet");
      MakeCleverHistograms(sighist_mm, label+"SSMM_CR");
      MakeCleverHistograms(sighist_mm, label+"SSMM_LowMass_SR");
      MakeCleverHistograms(sighist_mm, label+"SSMM_LowMass_CR");
      MakeCleverHistograms(sighist_mm, label+"SSMM_HighMass_SR");
      MakeCleverHistograms(sighist_mm, label+"SSMM_HighMass_CR");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_noB");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_FailnoB");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_passMET");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_passHT");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_passDR");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_passDR2");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_passMT");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_passLJ");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_passLTHT");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_passmll");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_passpt2");
      MakeCleverHistograms(sighist_mm, label+"SSMM_DiJet_passpt2mll");
      MakeCleverHistograms(sighist_mm, label+"SSMM_FatJet");
      MakeCleverHistograms(sighist_mm, label+"SSMM_NoFatJet");
      
      MakeCleverHistograms(sighist_mm, label+"OSMM_Preselection");
      MakeCleverHistograms(sighist_mm, label+"OSMM_1Jet");
      MakeCleverHistograms(sighist_mm, label+"OSMM_DiJet");
      MakeCleverHistograms(sighist_mm, label+"OSMM_CR");
      MakeCleverHistograms(sighist_mm, label+"OSMM_LowMass");
      MakeCleverHistograms(sighist_mm, label+"OSMM_MidMass");
      MakeCleverHistograms(sighist_mm, label+"OSMM_HighMass");
    }
    
    MakeCleverHistograms(sighist_mm, "TriggerTop");
    MakeCleverHistograms(sighist_mm, "TriggerOSZMM");
    MakeCleverHistograms(sighist_mm, "DiMuTriggerTop");
    MakeCleverHistograms(sighist_mm, "DiMuTriggerOSZMM");
    MakeCleverHistograms(sighist_mm, "DiMuTriggerOSMM");
    MakeCleverHistograms(sighist_mm, "DiJetTriggerOSMM");
    
    MakeCleverHistograms(sighist_mm, "TriggerTop_p1");
    MakeCleverHistograms(sighist_mm, "TriggerOSZMM_p1");
    MakeCleverHistograms(sighist_mm, "DiMuTriggerTop_p1");
    MakeCleverHistograms(sighist_mm, "DiMuTriggerOSZMM_p1");
    
    MakeCleverHistograms(sighist_mm, "TriggerTop_p2");
    MakeCleverHistograms(sighist_mm, "TriggerOSZMM_p2");
    MakeCleverHistograms(sighist_mm, "DiMuTriggerTop_p2");
    MakeCleverHistograms(sighist_mm, "DiMuTriggerOSZMM_p2");
    
    MakeCleverHistograms(sighist_mm, "OSZMM");
    
  }
  

  return;
}

void HNDiLepton::DoCutFlow(float w){

  std::vector<snu::KElectron> electrons_veto = GetElectrons(false,false,"ELECTRON_HN_VETO");
  std::vector<snu::KJet> alljets = GetJets("JET_NOLEPTONVETO", 20., 2.5);
  std::vector<snu::KJet> jets = GetJets("JET_HN");


  counter("1_NoCut",w);
  std::vector<snu::KMuon> muons = GetMuons("MUON_HN_TIGHT",false);
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return;
  if(!PassMETFilter()) return;     /// Initial event cuts :                                                                                                                                                                                                                    
  if((PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") || PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"))){
    counter( "2_EventClean", w);

    if(muons.size()==2){
      if(muons[1].Pt() > 10){
	if(muons[0].Pt() > 20){
	  counter( "3_DiMu", w);

	  std::vector<snu::KMuon> muons_veto         = GetMuons("MUON_HN_VETO",false);
	  if(muons_veto.size()!=2) {

	    m_logger << INFO  << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
	  }
	    
          if(muons_veto.size() ==2){
            counter("4_MuVeto", w);
            if(electrons_veto.size() ==0){
              counter("5_ElVeto", w);
	      if(GetDiLepMass(muons) > 10.){
		counter("6_mll", w);
		if(jets.size() > 1){
		  counter("7_DiJet", w);
		  CorrectedMETRochester(muons);
		  if(eventbase->GetEvent().PFMET()  < 50.){
		    counter("8_MET", w);
		    
		    float dijetmass_tmp=999.;
		    float dijetmass=9990000.;
		    int m=-999;
		    int n=-999;
		    for(UInt_t emme=0; emme<jets.size(); emme++){
		      for(UInt_t enne=1; enne<jets.size(); enne++) {
			if(emme == enne) continue;
			dijetmass_tmp = (jets[emme]+jets[enne]).M();
			if ( fabs(dijetmass_tmp-80.4) < fabs(dijetmass-80.4) ) {
			  dijetmass = dijetmass_tmp;
			  m = emme;
			  n = enne;
			}
		      }
		    }
		    snu::KParticle jj = jets.at(m) + jets.at(n) ;
		    if(jj.M() < 200) {
		      counter("9_mjj", w);
		      
		      int nbjet = NBJet(jets);
		      if(nbjet ==0){
			counter("91_bjetv1", w);
			int nbjetv2 = NBJet(alljets);
			if(nbjetv2 ==0){
			  counter("92_bjetv2", w);
			}
		      }
		    }// mjj
		  }// MET
		}
	      }
	    }
	  }
	}
      }
    }
  }
  

}


void HNDiLepton::ExecuteEvents()throw( LQError ){
  
  std::vector<snu::KMuon> loose_mu = GetMuons(_m_looseid,true);

  //  for(unsigned int ipdf=0; ipdf < eventbase->GetEvent().PdfWeights().size() ; ipdf++){
  //  cout << ipdf << " " << eventbase->GetEvent().PdfWeights().at(ipdf) << " " <<   eventbase->GetEvent().LHEWeight()<< " "<< eventbase->GetEvent().MCWeight() <<  endl;
    
  //}
  
  
  
  if(k_running_nonprompt){
    while(!fake_configured){
      /// SET UP OWN FAKES HISTS --> KEY ////                                                                                                                                                                                                                                    
      std::map<TString, std::pair<std::pair<TString,TString>  ,std::pair<float,TString> > >fake_hists;
      /// ONLY double (TH2D*) are configured at the minute
      fake_hists["fr_muon_central"] = std::make_pair(std::make_pair("Muon_Data_v7_SIP3_FR.root","Muon_Data_v7_SIP3_FR_Awayjet40"), std::make_pair(70., "TH2D"));
      fake_hists["fr_muon_awayjet20"] = std::make_pair(std::make_pair("Muon_Data_v7_SIP3_FR.root","Muon_Data_v7_SIP3_FR_Awayjet20"), std::make_pair(70., "TH2D"));
      fake_hists["fr_muon_awayjet30"] = std::make_pair(std::make_pair("Muon_Data_v7_SIP3_FR.root","Muon_Data_v7_SIP3_FR_Awayjet30"), std::make_pair(70., "TH2D"));
      fake_hists["fr_muon_awayjet60"] = std::make_pair(std::make_pair("Muon_Data_v7_SIP3_FR.root","Muon_Data_v7_SIP3_FR_Awayjet60"), std::make_pair(70., "TH2D"));
      fake_hists["fr_electron_central"] = std::make_pair(std::make_pair("Electron_Data_v7_FR.root","Electron_Data_v7_FR_Awayjet40") , std::make_pair(70., "TH2D"));
      /// END SET UP OWN FAKE HISTS ////                                                                                                                                                                                                                                           
      ConfigureFakeHists("/data1/LQAnalyzer_rootfiles_for_analysis/CATAnalysis2016/Fake/DiLep/", fake_hists);
    }
  }


  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;

  
  //// Set weight for MC (nevents / xsec (ish))
  if(!isData){
    weight*= MCweight;

  }
  if(!isData)weight*= GetKFactor();
  if(!isData)weight*= MC_CR_Correction(0);
  

  if(k_flags.size() > 0){
    
    if(k_flags[0] == "CutFlow")  functionality = HNDiLepton::CUTFLOW;
    else if(k_flags[0] == "OPT")  functionality = HNDiLepton::OPT;
    else functionality = HNDiLepton::ANALYSIS;
  }

  
  if(functionality==HNDiLepton::CUTFLOW){
    DoCutFlow(1.); return;
  }
  
  _mm_channel =   isData ?  (k_channel.Contains("DoubleMuon")) : true ;
  _m_channel =   isData ?  (k_channel.Contains("SingleMuon")) : true ;
  _ee_channel =   isData ?  (k_channel.Contains("DoubleEG")) : true ;
  _e_channel =   isData ?  (k_channel.Contains("SingleElectron")) : true ;
  _em_channel =   isData ?  (k_channel.Contains("MuonEG")) : true ;

  TString muid=_m_tightid;
  if(k_running_nonprompt) muid=_m_looseid;
  
  TString elid="ELECTRON_HN_TIGHTv4";
  if(k_running_nonprompt) elid="ELECTRON_HN_FAKELOOSEv7";

  
  std::vector<snu::KMuon> muons = GetMuons(muid,true);
  std::vector<snu::KMuon> muons_veto = GetMuons("MUON_HN_VETO",true);
  CorrectedMETRochester(muons);
  //double  met = eventbase->GetEvent().PFMET();
  //double  met_phi = eventbase->GetEvent().PFMETPhi();  

  //CorrectedMETMuon(1, muons, met, met_phi);

  /// Set up electrons                                                                                                                                        
  std::vector<snu::KElectron> electrons = GetElectrons(elid);
  std::vector<snu::KElectron> electrons_veto = GetElectrons("ELECTRON_HN_VETO");

  std::vector<snu::KJet> alljets = GetJets("JET_NOLEPTONVETO", 20., 2.5);
  std::vector<snu::KJet> tchanjets = GetJets("JET_HN_TChannel",20., 5.);
  
  std::vector<snu::KJet>  jets = GetJetsWFT("JET_HN", "FATJET_HN");
  std::vector<snu::KFatJet> fatjetcoll = GetFatJets("FATJET_HN");
  std::vector<snu::KFatJet> fatjetcoll_all = GetFatJets("FATJET_NOCUT");
  cout << fatjetcoll.size() << " " << fatjetcoll_all.size() << endl;

  vector<int> ijets;  

  if(functionality==HNDiLepton::OPT){
    return;
  }

  if(functionality==HNDiLepton::ANALYSIS){
    
    float mm_weight=weight*MMWeight(muons,muid);
    float ee_weight=weight*EEWeight(electrons,elid);
    float em_weight=weight*EMWeight(electrons, muons,elid, muid);
    
    std::vector<TString> triggerlist_DiMuon, triggerlist_DiElectron;
    std::vector<TString> triggerlist_DiMuon_singleLep, triggerlist_DiElectron_singleLep;

    std::vector<TString> triggerlist_EMu, triggerlist_EMu_PeriodBtoG, triggerlist_EMu_PeriodH, triggerlist_EMu_Mu8Ele23, triggerlist_EMu_Mu23Ele8;
    triggerlist_DiMuon.clear();
    triggerlist_DiElectron.clear();
    triggerlist_EMu.clear();

    triggerlist_DiMuon.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
    triggerlist_DiMuon.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

    triggerlist_DiMuon_singleLep.push_back("HLT_IsoMu24_v");
    triggerlist_DiMuon_singleLep.push_back("HLT_IsoTkMu24_v");

    triggerlist_DiElectron.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    triggerlist_DiElectron_singleLep.push_back("HLT_Ele27_WPTight_Gsf_v");

    triggerlist_EMu.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    triggerlist_EMu.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
    triggerlist_EMu.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    triggerlist_EMu.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");
    triggerlist_EMu_PeriodBtoG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    triggerlist_EMu_PeriodBtoG.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");

    triggerlist_EMu_PeriodH.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
    triggerlist_EMu_PeriodH.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");

    triggerlist_EMu_Mu8Ele23.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
    triggerlist_EMu_Mu8Ele23.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");

    triggerlist_EMu_Mu23Ele8.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
    triggerlist_EMu_Mu23Ele8.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");

    
    std::vector<snu::KMuon> pogmuons = GetMuons("MUON_POG_TIGHT",true);
    if(_mm_channel)RunMM("DiMuon",muons,muons_veto,electrons_veto,alljets,jets, fatjetcoll,tchanjets,mm_weight, triggerlist_DiMuon,20., 10.);
    if(_m_channel)RunMM("DiMuon_isomu",muons,muons_veto,electrons_veto,alljets,jets, fatjetcoll,tchanjets,mm_weight, triggerlist_DiMuon_singleLep,25., 5.);
    if(_mm_channel)RunMM("DiMuon_allFJ",muons,muons_veto,electrons_veto,alljets,jets, fatjetcoll_all,tchanjets,mm_weight, triggerlist_DiMuon,20., 10.);
    //RunMM("PogDiMuon",pogmuons,muons_veto,electrons_veto,alljets,jets, fatjetcoll,tchanjets,mm_weight, triggerlist_DiMuon,20., 10.);

    //RunEE();
    //    RunEM();
  }
  //cout << "PT cone " << muons_test[0].PTCone(0.4,eventbase->GetMuonSel()->IsoCutValue(muons_test[0],"MUON_HN_TIGHT")) << " " << muons_test.at(0).Pt()*(1+max(0.,(muons_test.at(0).RelIso04()-0.07))) << endl;
}


float HNDiLepton::EEWeight(std::vector<snu::KElectron> electrons,TString id){
  
  if(isData) return 1.;
  
  double mc_weight = mcdata_correction->ElectronScaleFactor(id, electrons, 0);
  double electron_RecoSF =  mcdata_correction->ElectronRecoScaleFactor(electrons);
  mc_weight*= electron_RecoSF;

  float puweight=1.;
  puweight=eventbase->GetEvent().PileUpWeight(snu::KEvent::down);
  mc_weight*= puweight;

  std::vector<snu::KMuon> muons;
  double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, id, muons, "muid", 1, 0, 0);
  double trigger_eff_MC = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, id, muons, "muid", 1, 1, 0);
  float trigger_sf = trigger_eff_Data/trigger_eff_MC;
  mc_weight*= trigger_sf;


  
  return mc_weight;
}
float HNDiLepton::EMWeight(std::vector<snu::KElectron> electrons,std::vector<snu::KMuon> muons,TString idel, TString idmu){
  
  if(isData) return 1.;
  double mc_weight = mcdata_correction->ElectronScaleFactor(idel, electrons, 0);
  double electron_RecoSF =  mcdata_correction->ElectronRecoScaleFactor(electrons);
  mc_weight*= electron_RecoSF;

  float puweight=1.;
  puweight=eventbase->GetEvent().PileUpWeight(snu::KEvent::down);
  mc_weight*= puweight;
  mc_weight*= mcdata_correction->MuonTrackingEffScaleFactor(muons);
  mc_weight*= mcdata_correction->MuonScaleFactor(idel, muons, 0);
  return mc_weight;
}


float HNDiLepton::MMWeight(std::vector<snu::KMuon> muons, TString id){

  if(isData) return 1.;
  
  float mc_weight=mcdata_correction->MuonTrackingEffScaleFactor(muons);

  float puweight=1.;
  //if(GetMCPeriodRandom() < 6) puweight=eventbase->GetEvent().PileUpWeight(snu::KEvent::down);
  //else puweight= mcdata_correction->UserPileupWeight(eventbase->GetEvent(), jets.size());
  puweight=eventbase->GetEvent().PileUpWeight(snu::KEvent::down);                                                                                      
  mc_weight*= puweight;

  double muon_id_iso_sf = mcdata_correction->MuonScaleFactor(id, muons, 0);
  mc_weight*= muon_id_iso_sf;
  
  std::vector<snu::KElectron> el;
  double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(el, "elid", muons, id, 0, 0, 0);
  double trigger_eff_MC = mcdata_correction->TriggerEfficiencyLegByLeg(el, "elid", muons, id, 0, 1, 0);
  float trigger_sf = trigger_eff_Data/trigger_eff_MC;
  mc_weight*= trigger_sf;
  
  return mc_weight;

}

void HNDiLepton::RunMM(TString label, vector<snu::KMuon> muons, vector<snu::KMuon> muons_veto, vector<snu::KElectron> electrons, vector<snu::KJet> alljets, vector<snu::KJet> jets, vector<snu::KFatJet> fatjets, vector<snu::KJet> tjets,float mm_weight ,vector<TString> mm_trig, float pt1, float pt2){
  
  
  FillEventCutFlow(4,"NoCut",mm_weight*WeightByTrigger(mm_trig[0],TargetLumi),label );


  /// efine event variables                                                                                                                                                  
  snu::KParticle mm;
  if(muons.size() ==2) mm= muons[0] + muons[1];
  int nbjet_l=0;
  int nbjet_m=0;
  int nbjet_t=0;
  float ST (0.);
  if(muons.size() ==2)ST= muons[0].Pt() + muons[1].Pt();

  for(unsigned int ij=0; ij <alljets.size(); ij++){
    if( alljets[ij].IsBTagged(  snu::KJet::CSVv2, snu::KJet::Loose))  nbjet_l++;
    if( alljets[ij].IsBTagged(  snu::KJet::CSVv2, snu::KJet::Medium)) nbjet_m++;
    if( alljets[ij].IsBTagged(  snu::KJet::CSVv2, snu::KJet::Tight))  nbjet_t++;
  }


  float HT=0.;
  for(unsigned int ij=0; ij <jets.size(); ij++){HT+= jets[ij].Pt();       ST+= jets[ij].Pt();}
  for(unsigned int ij=0; ij < fatjets.size(); ij++){ HT+= fatjets[ij].Pt(); ST+= fatjets[ij].Pt();}
  ST+=eventbase->GetEvent().PFMET();
  
  
  float met = eventbase->GetEvent().PFMET();
  float met2_st = pow(eventbase->GetEvent().PFMET(),2.)/ ST;
  
  bool closejet=false;
  bool closejet2=false;
  float mindR=9999.;
  for(unsigned int im=0; im <muons.size(); im++){
    for(unsigned int ij2=0; ij2 <alljets.size(); ij2++){
      if(muons.at(im).DeltaR(alljets.at(ij2)) <  0.5) closejet=true;
      else{
	if(muons.at(im).DeltaR(alljets.at(ij2)) <  mindR){
	  mindR = muons.at(im).DeltaR(alljets.at(ij2));
	}
      }
    }
  }

  float min_eleadawayjet_Dr= -10000.;
  for(UInt_t i=0; i<muons.size(); i++){
    for(UInt_t emme=0; emme<jets.size(); emme++){
      float dR =muons[i].DeltaR(jets[emme]);
      if(dR > 0.5) {
	if(dR > min_eleadawayjet_Dr) {
	  min_eleadawayjet_Dr=  dR;
	  
	}
      }
    }
  }
  
  float mindRjj=999.;
  if(jets.size() > 1){
    for(unsigned int ij=0; ij <jets.size()-1; ij++){
      for(unsigned int ij2=ij+1; ij2 <jets.size(); ij2++){
	if(jets.at(ij).DeltaR(jets.at(ij2)) < mindRjj) mindRjj=jets.at(ij).DeltaR(jets.at(ij2)) ;
      }
    }
  }

  
  if(muons.size()==2){
    if(muons[1].Pt() > pt2 && muons[0].Pt() > pt1 ){

      FillHist("Nveto_1",(muons_veto.size() +electrons.size()), mm_weight*WeightByTrigger(mm_trig[0], TargetLumi));
      
      if(!PassTriggerOR(mm_trig)) return;
      if(!isData)mm_weight*= WeightByTrigger(mm_trig[0], TargetLumi) ;

      FillHist("Nveto_2",(muons_veto.size() +electrons.size()), mm_weight);
      FillEventCutFlow(4,"Trigger",mm_weight,label);
	
      counter("Trigger",mm_weight);
	
      if(electrons.size()>0) return;
	
      counter("ElVeto",mm_weight);
	
      FillEventCutFlow(4,"ElVeto",mm_weight,label);
	
      /// Apply weights to MC
          
      
      FillEventCutFlow(4,"DiLooseMM",mm_weight,label);	
      
      //if(CheckEventComparison("jalmond","test19581_periodB_SKDoubleMuon_dilep_cat_v8-0-7_HNDiLepton","jalmond", "test19581_periodB_SKDoubleMuon_dilep_cat_v8-0-7_HNDiLepton2",false)){

      vector<int> ijets;  
      
      if(SameCharge(muons)){
	
	FillHist("Nveto_3",(muons_veto.size() +electrons.size()), mm_weight);

	FillEventCutFlow(4,"SSMMLoose",mm_weight,label);
	
	if((muons[0] + muons[1]).M()  > _mm_mll_presel_cut) {
	    
	  FillHist("Nveto_4",(muons_veto.size() +electrons.size()), mm_weight);


	  FillEventCutFlow(4,"MLL",mm_weight,label);
	    
	    FillCLHist(sighist_mm, label+"SSMM_nolepveto", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets,  mm_weight);
	    if(muons_veto.size()==2 && electrons.size()==0){
	      if(jets.size() >= 2 )   FillEventCutFlow(0,"DiJet",mm_weight,label);
	      
	      FillEventCutFlow(4,"LepVeto",mm_weight,label);
	      
	      //FillEventComparisonFile("test");
	      
	      float mm_weight_fix=mm_weight;
	      float mm_weight_sf=mm_weight;
	      float mm_weight_df=mm_weight;
	      float mm_weight_awayjet60=mm_weight;
	      float mm_weight_awayjet20=mm_weight;
	      if(k_running_nonprompt){
		
		mm_weight=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, muons, "MUON_HN_TIGHT", "ptcone", "fr_muon_central",0);
		mm_weight_awayjet20=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, muons, "MUON_HN_TIGHT", "ptcone", "fr_muon_awayjet20");
		mm_weight_awayjet60=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, muons, "MUON_HN_TIGHT", "ptcone", "fr_muon_awayjet60");
		mm_weight_sf=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, muons, "MUON_HN_TIGHT", "ptcone", "fr_muon_central",1);
		mm_weight_df=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, muons, "MUON_HN_TIGHT", "ptcone", "fr_muon_central",2);
		
		if(muons[1].Pt()> 30) mm_weight_fix=mm_weight*2.;

		if(PassID(muons[0], "MUON_HN_TIGHT") && PassID(muons[1] , "MUON_HN_TIGHT"))counter("SSMM_TT_Preselection", mm_weight);
		if(PassID(muons[0], "MUON_HN_TIGHT") && PassID(muons[1] , "MUON_HN_TIGHT"))counter("SSMM_TT_Preselection_now", 1.);
		if(PassID(muons[0], "MUON_HN_TIGHT") && !PassID(muons[1] , "MUON_HN_TIGHT"))counter("SSMM_TL_Preselection", mm_weight);
		if(PassID(muons[0], "MUON_HN_TIGHT") && !PassID(muons[1] , "MUON_HN_TIGHT"))counter("SSMM_TL_Preselection_now", 1.);
		if(!PassID(muons[0], "MUON_HN_TIGHT") && PassID(muons[1] , "MUON_HN_TIGHT"))counter("SSMM_LT_Preselection", mm_weight);
		if(!PassID(muons[0], "MUON_HN_TIGHT") && PassID(muons[1] , "MUON_HN_TIGHT"))counter("SSMM_LT_Preselection_now", 1.);
		if(!PassID(muons[0], "MUON_HN_TIGHT") && !PassID(muons[1] , "MUON_HN_TIGHT"))counter("SSMM_LL_Preselection", mm_weight);
		if(!PassID(muons[0], "MUON_HN_TIGHT") && !PassID(muons[1] , "MUON_HN_TIGHT"))counter("SSMM_LL_Preselection_now", 1.);
		counter("SSMM_Preselection_sf", mm_weight_sf);
		counter("SSMM_Preselection_df", mm_weight_df);

	      }
	      
	      counter("SSMM_Preselection", mm_weight);
	      
	      FillCLHist(sighist_mm, label+"SSMM_Preselection", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets,  mm_weight);
	      FillCLHist(sighist_mm, label+"SSMM_Preselection_aj60", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets,  mm_weight_awayjet60);
	      FillCLHist(sighist_mm, label+"SSMM_Preselection_aj20", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets,  mm_weight_awayjet20);
	      FillCLHist(sighist_mm, label+"SSMM_Preselection_sf",  eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets, mm_weight_sf);
	      FillCLHist(sighist_mm, label+"SSMM_Preselection_df",  eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets, mm_weight_df);
	      if((muons[0].DeltaR(muons[1])) > 0.4){
		FillEventCutFlow(4,"mmdR",mm_weight,label);

		FillHist("Nveto_5",(muons_veto.size() +electrons.size()), mm_weight);
		
		if(jets.size()==1) {
		  if(GetDiLepMass(muons) > 110){
		    FillCLHist(sighist_mm, label+"SSMM_1JetCR", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets,  mm_weight);
		    FillEventCutFlow(4,"SSMM_1Jet",mm_weight,label);
		  }
		  else{
		    if(GetDiLepMass(muons) > 15){
		      if(eventbase->GetEvent().PFMET() < 60.) {
			if(NBJet(alljets) ==0){
			  if(muons[0].DeltaR(muons[1])  < 3.4) {
			    if( (GetLT(muons)/(SumPt(jets)+SumPt(fatjets)) ) > 0.8){
			      FillEventCutFlow(4,"SSMM_1JetOpt",mm_weight,label);
			      FillCLHist(sighist_mm, label+"SSMM_1JetOpt", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets,  mm_weight);
			      if(GetDiLepMass(muons) < 60){
				FillCLHist(sighist_mm, label+"SSMM_1JetOpt2", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets,  mm_weight);
				FillEventCutFlow(4,"SSMM_1JetOpt2",mm_weight,label);
			      }
			      if(muons[1].Pt() > 15){
				FillCLHist(sighist_mm, label+"SSMM_1JetOpt_pt2", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets,  mm_weight);
				if(GetDiLepMass(muons) < 60)FillCLHist(sighist_mm, label+"SSMM_1JetOpt2_pt2", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets,  mm_weight);
				FillEventCutFlow(4,"SSMM_1JetOpt_pt2",mm_weight,label);
				if(GetDiLepMass(muons) < 60)FillEventCutFlow(4,"SSMM_1JetOpt2_pt2",mm_weight);
				
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
		
		if(jets.size() == 0)                   FillCLHist(sighist_mm, label+"SSMM_0Jet", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
		if(jets.size() == 1)                   FillCLHist(sighist_mm, label+"SSMM_1Jet", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
		if(jets.size() == 2)                   FillCLHist(sighist_mm, label+"SSMM_2Jet", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
		if(jets.size() == 3)                   FillCLHist(sighist_mm, label+"SSMM_3Jet", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
		if(jets.size() >= 4)                   FillCLHist(sighist_mm, label+"SSMM_4Jet", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
		
		if(muons[1].Pt() > 35){
		  if(NBJet(alljets)==0){
                    if(eventbase->GetEvent().PFMET() < 60.) {
		      if(eventbase->GetEvent().PFMET() > 20.) {
			
			if(jets.size() == 0)                   FillCLHist(sighist_mm, label+"SSMM_0Jet_xc", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
			if(jets.size() == 1)                   FillCLHist(sighist_mm, label+"SSMM_1Jet_xc", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
			if(jets.size() >= 2) {
			  FillCLHist(sighist_mm, label+"SSMM_2Jet_xc", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
			  if(k_running_nonprompt){
			    if(PassID(muons[0], "MUON_HN_TIGHT") && !PassID(muons[1] , "MUON_HN_TIGHT")) FillCLHist(sighist_mm, label+"SSMM_2Jet_TL_xc", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,1.);
			    if(!PassID(muons[0], "MUON_HN_TIGHT") && PassID(muons[1] , "MUON_HN_TIGHT")) FillCLHist(sighist_mm, label+"SSMM_2Jet_LT_xc", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,1.);
			  }
			  
			  if(muons[0].Charge() < 0)                           FillCLHist(sighist_mm, label+"SSMM_2Jet_xc_minus", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
			  if(muons[0].Charge() > 0)                           FillCLHist(sighist_mm, label+"SSMM_2Jet_xc_plus", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
			  if(eventbase->GetEvent().PFMET() < 60.) {
			    if(eventbase->GetEvent().PFMET() > 20.) {
			      FillCLHist(sighist_mm, label+"SSMM_2Jet_xc_met", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);

			      FillCLHist(sighist_mm, label+"SSMM_2Jet_xc_met_cfj", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
			      if(GetMasses("l1jj", muons, jets,  ijets,false) > 110. && GetMasses("l1jj", muons, jets,  ijets,false) < 140) FillCLHist(sighist_mm, label+"SSMM_2Jet_xc_met_cfj_l1mass", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
			      if(GetMasses("l1jj", muons, jets,  ijets,false) > 320.)FillCLHist(sighist_mm, label+"SSMM_2Jet_xc_met_cfj_l1masshigh", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
			      
			    }
			  }
			  if((muons[0].Pt() + muons[1].Pt()) > 125){

			    FillCLHist(sighist_mm, label+"SSMM_2Jet_xc_lt", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);

			  }
			  
			  if(!(GetMasses("l1jj", muons, jets,  ijets,false) > 110. && GetMasses("l1jj", muons, jets,  ijets,false) < 140)){
			    FillCLHist(sighist_mm, label+"SSMM_2Jet_rmbin", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);

			  }
			  
			}


		      }
		    }
		  }
		}

	      	
		if(jets.size() > 1){
		  FillHist("Nveto_6",(muons_veto.size() +electrons.size()), mm_weight);

		  FillCLHist(sighist_mm, label+"SSMM_DiJet", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
		  FillCLHist(sighist_mm, label+"SSMM_DiJet_noW", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,WeightByTrigger(mm_trig[0],TargetLumi)* GetKFactor()*weight);
		  if((muons[0].DeltaR(muons[1])) <  1.)                   FillCLHist(sighist_mm, label+"SSMM_DiJet_xc", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
		  
		  FillCLHist(sighist_mm, label+"SSMM_DiJet_fix", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight_fix);
		  
		  
		  if(GetMasses("l2jj", muons, jets,  ijets,true) > 250 && GetMasses("l2jj", muons, jets,  ijets,true) < 400) FillCLHist(sighist_mm, label+"SSMM_DiJet_xc2", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
		  
		  if(NBJet(alljets) ==0)                   FillCLHist(sighist_mm, label+"SSMM_DiJet_noB", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
		  else FillCLHist(sighist_mm, label+"SSMM_DiJet_FailnoB", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
		  if(NBJet(alljets) ==0)   {
		    FillEventCutFlow(4,"SSMM_DiJet_noB",mm_weight,label);
		    if(eventbase->GetEvent().PFMET() < 80.) {
		      FillCLHist(sighist_mm, label+"SSMM_DiJet_passMET", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
		    }
		    
		    if((eventbase->GetEvent().PFMET() < 60.)) {
		      if(HT < 300){
			FillEventCutFlow(4,"SSMM_DiJet_passHT300",mm_weight,label);
			FillCLHist(sighist_mm, label+"SSMM_DiJet_passHT3", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
			bool closejet=false;
			bool closejet2=false;
			float mindR=9999.;
			for(unsigned int im=0; im <muons.size(); im++){
			  for(unsigned int ij2=0; ij2 <jets.size(); ij2++){
			    if(muons.at(im).DeltaR(jets.at(ij2)) <  0.5) closejet=true;
			    if(muons.at(im).DeltaR(jets.at(ij2)) <  mindR){
			      mindR = muons.at(im).DeltaR(jets.at(ij2));
			    }
			  }
			}
			if(mindR > 3.5) closejet2=true;
			if(!closejet){
			  
			  if(muons[0].DeltaR(muons[1])  < 4.) FillCLHist(sighist_mm, label+"SSMM_DiJet_passDR", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
			  if(muons[0].DeltaR(muons[1])  < 3.4) {
			    FillEventCutFlow(4,"SSMM_DiJet_passDR2",mm_weight,label);

			    FillCLHist(sighist_mm, label+"SSMM_DiJet_passDR2", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);

			    float dphi_1 = fabs(TVector2::Phi_mpi_pi(muons.at(0).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet)));
			    float MT_1 = sqrt(2.* muons.at(0).Et()*eventbase->GetEvent().PFMET() * (1 - cos( dphi_1)));

			    float dphi_2 = fabs(TVector2::Phi_mpi_pi(muons.at(1).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet)));
			    float MT_2 = sqrt(2.* muons.at(1).Et()*eventbase->GetEvent().PFMET() * (1 - cos( dphi_2)));

			    //if( ( (muons[0].Pt() + muons[1].Pt()) / HT ) > 0.25){
			    //FillCLHist(sighist_mm, label+"SSMM_DiJet_passLTHT", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);	


			    if(!closejet2){
			      FillEventCutFlow(4,"SSMM_DiJet_passLJ",mm_weight,label);

			      FillCLHist(sighist_mm, label+"SSMM_DiJet_passLJ", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
			      if(MT_1 < 100 && MT_2 < 100)                                 FillEventCutFlow(4,"SSMM_DiJet_passMT100",mm_weight,label);
			      if(MT_1 < 90 && MT_2 < 90)                                 FillEventCutFlow(4,"SSMM_DiJet_passMT90",mm_weight,label);
			      
			      if(MT_1 < 80 || MT_2 < 80){
				FillEventCutFlow(4,"SSMM_DiJet_passMT80",mm_weight,label);

				FillCLHist(sighist_mm, label+"SSMM_DiJet_passMT", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
				if(muons[1].Pt() > 15) {
				  FillEventCutFlow(4,"SSMM_DiJet_passpt2",mm_weight,label);

				  FillCLHist(sighist_mm, label+"SSMM_DiJet_passpt2", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
				}

				int s_awayjetdr(0), s_mll(0), s_pt(0), s_mllj(0), s_met(0), s_l2jj(0), s_lljj(0), s_contraMass(0), s_st(0);
				int nsig_reg=0;
				if(GetDiLepMass(muons) < 70) s_mll=1;
				else s_mll=2;
				
				if(muons[1].Pt() < 15) s_pt=1;
				else s_pt=2;
				
				if((muons[0]+muons[1]+jets[0]).M() < 100) s_mllj=1;
				else if((muons[0]+muons[1]+jets[0]).M() < 150)  s_mllj=2;
				else if((muons[0]+muons[1]+jets[0]).M() < 175)  s_mllj=3;
				else if((muons[0]+muons[1]+jets[0]).M() < 200)  s_mllj=4;
				else if((muons[0]+muons[1]+jets[0]).M() < 225)s_mllj=5;
				else if((muons[0]+muons[1]+jets[0]).M() < 250)s_mllj=6;
				else if((muons[0]+muons[1]+jets[0]).M() < 300)s_mllj=7;

				if(eventbase->GetEvent().PFMET() < 20)    s_met=1;
				else if(eventbase->GetEvent().PFMET() < 40)    s_met=2;
				else if(eventbase->GetEvent().PFMET() < 60)    s_met=3;
				
				if(ST < 150) s_st=1;
				else if(ST < 175) s_st=2;
				else if(ST < 200) s_st=3;
				else if(ST < 250) s_st=4;
				else s_st=5;
				
				if(GetMasses("l2jj", muons, jets, ijets, true) < 80) s_l2jj=1;
				else if(GetMasses("l2jj", muons, jets, ijets, true) < 120) s_l2jj=2;
				else if(GetMasses("l2jj", muons, jets,  ijets,true) < 200) s_l2jj=3;
				else if(GetMasses("l2jj", muons, jets, ijets, true) < 300) s_l2jj=4;
				
				if(GetMasses("lljj", muons, jets,  ijets,true) < 100) s_lljj=1;
				else if(GetMasses("lljj", muons, jets,  ijets,true) < 150) s_lljj=2;
				else if(GetMasses("lljj", muons, jets,  ijets,true) < 175) s_lljj=3;
				else if(GetMasses("lljj", muons, jets,  ijets,true) < 200) s_lljj=4;
				else if(GetMasses("lljj", muons, jets,  ijets,true) < 250) s_lljj=5;
				else if(GetMasses("lljj", muons, jets,  ijets,true) < 300) s_lljj=6;
				else if(GetMasses("lljj", muons, jets,  ijets,true) < 400) s_lljj=7;
				
				if(GetMasses("contMT", muons, jets,  ijets,true) < 25) s_contraMass=1;
				else 				if(GetMasses("contMT", muons, jets,  ijets,true) < 50) s_contraMass=2;
				else 				if(GetMasses("contMT", muons, jets,  ijets,true) < 100) s_contraMass=3;

				
				float min_eleadawayjet_Dr= -10000.;
				for(UInt_t i=0; i<muons.size(); i++){
				  for(UInt_t emme=0; emme<jets.size(); emme++){
				    float dR =muons[i].DeltaR(jets[emme]);
				    if(dR > 0.5) {
				      if(dR > min_eleadawayjet_Dr) {
					min_eleadawayjet_Dr=  dR;
					
				      }
				    }
				  }
				}
				
				if(min_eleadawayjet_Dr < 3.5)  s_awayjetdr=1;
				else  s_awayjetdr=2;
				

				for(unsigned int ipt =1; ipt < 3; ipt++){
                                  for(unsigned int imll =1; imll < 3; imll++){
				    for(unsigned int  illjj =1; illjj < 8; illjj++){
				      for(unsigned int imet=1; imet< 4; imet++,nsig_reg++){
					if(s_mll ==imll){
					  if(s_pt==ipt){
					    if(s_lljj==illjj){
					      if(s_met==imet){
						FillHist(label+"Sig_Regions1",nsig_reg, mm_weight, 0., 90., 90);
					      }
					    }
					  }
					}
				      }
				    }
				  }
				}

				nsig_reg=0;
				for(unsigned int ipt =1; ipt < 3; ipt++){
                                  for(unsigned int imll =1; imll < 3; imll++){
				    for(unsigned int i_cm=1; i_cm < 4; i_cm++){
				      for(unsigned int i_st=1; i_st < 6; i_st++,nsig_reg++){
					if(s_mll ==imll){
                                          if(s_pt==ipt){
					    if(s_st==i_st){
					      if(s_contraMass==i_cm){
                                                FillHist(label+"Sig_Regions2",nsig_reg, mm_weight, 0., 90., 90);
					      }
					    }
					  }
					}
				      }
				    }
				  }
				}
				nsig_reg=0;
                                for(unsigned int ipt =1; ipt < 3; ipt++){
                                  for(unsigned int imll =1; imll < 3; imll++){
				    for(unsigned int iaj =1; iaj < 3; iaj++){
				      for(unsigned int imet=1; imet< 4; imet++, nsig_reg++){
					for(unsigned int i_cm=1; i_cm < 4; i_cm++){
					  
					  //for(unsigned int i_llj=0; i_llj < 4; i_llj++, nsig_reg++){
					  if(s_mll ==imll){
					    if(s_awayjetdr==iaj){
					      if(s_pt==ipt){
						if(s_met==imet){
						  if(s_contraMass==i_cm){
						    
						    FillHist(label+"Sig_Regions3",nsig_reg, mm_weight, 0., 90., 90);
						  }
						}
					      }
					    }
					  }
					}
				      }
				    }
				  }
				}
				
				if(GetDiLepMass(muons) < 70){
				  FillEventCutFlow(4,"SSMM_DiJet_passmll",mm_weight,label);
				  
				  FillCLHist(sighist_mm, label+"SSMM_DiJet_passmll", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);   
				  if(muons[1].Pt() > 15) {
				    FillEventCutFlow(4,"SSMM_DiJet_passpt2mll",mm_weight,label);
				    
				    FillCLHist(sighist_mm, label+"SSMM_DiJet_passpt2mll", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
				  }
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
		
		
		if(fatjets.size()==0) {
		  FillCLHist(sighist_mm, label+"SSMM_NoFatJet", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
		  FillEventCutFlow(4,"SSMM_NoFatJet",mm_weight,label);
		  
		}
		else {
		  FillCLHist(sighist_mm, label+"SSMM_FatJet", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
		  FillEventCutFlow(4,"SSMM_FatJet",mm_weight,label);
		}
		
		if(tjets.size()> 1)FillCLHist(sighist_mm, label+"SSMM_DiTJet", eventbase->GetEvent(), muons, electrons,tjets, alljets,  fatjets,mm_weight);
		

		if(CheckSignalRegion(true,muons,electrons , jets, alljets,"Low", mm_weight)){
		  if(met < 80)FillCLHist(sighist_mm, label+"SSMM_LowMass_SR", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
		  if(met > 100 || nbjet_m > 0)FillCLHist(sighist_mm, label+"SSMM_LowMass_CR", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
		}

		if(CheckSignalRegion(true,muons,electrons , jets, alljets,"High", mm_weight)){
		  if(met2_st < 15) FillCLHist(sighist_mm, label+"SSMM_HighMass_SR", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
		  if(met2_st > 20 || nbjet_m > 0)FillCLHist(sighist_mm, label+"SSMM_HighMass_CR", eventbase->GetEvent(), muons, electrons,jets, alljets,  fatjets,mm_weight);
		}
	      }
	    }
	}
      }
      
  
    

      else{
	if(muons[1].Pt() > 20 && muons[0].Pt() > 25 ){
	  
	  if(PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")) {
	    if(_mm_channel){
	      FillCLHist(sighist_mm,"DiMuTriggerOSMM", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets, mm_weight);
	      if(jets.size() > 1) 		FillCLHist(sighist_mm,"DiJetTriggerOSMM", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets, mm_weight);
	    }
	    
	  }
	  
	  
	  if(GetDiLepMass(muons)< 100. && GetDiLepMass(muons) > 80.) {
	    FillCLHist(sighist_mm, label+"OSZMM", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets,  mm_weight);
	    TString analysis_trigger_muon="HLT_IsoMu24_v";
	    TString analysis_trigger_tkmuon="HLT_IsoTkMu24_v";
	    if((PassTrigger(analysis_trigger_muon) || PassTrigger(analysis_trigger_tkmuon))){
	      float trigw=1.;
	      if(_m_channel){
		FillCLHist(sighist_mm,"TriggerOSZMM", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets,  trigw*mm_weight);
		if(AllPrompt(muons,0))FillCLHist(sighist_mm,"TriggerOSZMM_p1", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets,  trigw*mm_weight);
		if(AllPrompt(muons,1))FillCLHist(sighist_mm,"TriggerOSZMM_p2", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets,  trigw*mm_weight);
	      }
	    }
	    if(PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")) {
	      if(_mm_channel){
		FillCLHist(sighist_mm,"DiMuTriggerOSZMM", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets, mm_weight);
		if(AllPrompt(muons,0)) FillCLHist(sighist_mm,"DiMuTriggerOSZMM_p1", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets, mm_weight);
		if(AllPrompt(muons,1)) FillCLHist(sighist_mm,"DiMuTriggerOSZMM_p2", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets, mm_weight);
	      }	    
	      
	    }
	  }
	  else{
	    if(jets.size() > 3){
	      if(NBJet(jets) > 1){
		
		TString analysis_trigger_muon="HLT_IsoMu24_v";
		TString analysis_trigger_tkmuon="HLT_IsoTkMu24_v";
		if((PassTrigger(analysis_trigger_muon) || PassTrigger(analysis_trigger_tkmuon))){
		  float trigw=1.;

		  if(_m_channel){
		    FillCLHist(sighist_mm, label+"TriggerTop", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets,  trigw*mm_weight);
		    if(AllPrompt(muons,0) )                 FillCLHist(sighist_mm, label+"TriggerTop_p1", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets,  trigw*mm_weight);
		    if(AllPrompt(muons,1))                  FillCLHist(sighist_mm, label+"TriggerTop_p2", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets,  trigw*mm_weight);
		    
		  }
		}
		if(PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")) {
		  if(_mm_channel){
		    if(muons[1].Pt() > 20 && muons[0].Pt() > 40 ){
		      
		      FillCLHist(sighist_mm, label+"DiMuTriggerTop", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets, mm_weight);
		      if(AllPrompt(muons,0)) 		FillCLHist(sighist_mm, label+"DiMuTriggerTop_p1", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets, mm_weight);
		      if(AllPrompt(muons,1) )		FillCLHist(sighist_mm, label+"DiMuTriggerTop_p2", eventbase->GetEvent(), muons, electrons,jets, alljets, fatjets, mm_weight);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
    

  
  return;
}


void HNDiLepton::FillEventCutFlow(int cf,TString cut,  float weight, TString label){
  
  if(cf==0){
    if(GetHist( label+"mm_eventcutflow")) {
      GetHist( label+"mm_eventcutflow")->Fill(cut,weight);
      
    }
    else{
      AnalyzerCore::MakeHistograms( label+"mm_eventcutflow",15,0.,15.);

      GetHist( label+"mm_eventcutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
      GetHist( label+"mm_eventcutflow")->GetXaxis()->SetBinLabel(2,"DiLooseMM");
      GetHist( label+"mm_eventcutflow")->GetXaxis()->SetBinLabel(3,"ElVeto");
      GetHist( label+"mm_eventcutflow")->GetXaxis()->SetBinLabel(4,"L1Pt");
      GetHist( label+"mm_eventcutflow")->GetXaxis()->SetBinLabel(5,"L2Pt5");
      GetHist( label+"mm_eventcutflow")->GetXaxis()->SetBinLabel(6,"L2Pt10");
      GetHist( label+"mm_eventcutflow")->GetXaxis()->SetBinLabel(7,"L2Pt15");
      GetHist( label+"mm_eventcutflow")->GetXaxis()->SetBinLabel(8,"L153");
      GetHist( label+"mm_eventcutflow")->GetXaxis()->SetBinLabel(9,"L153L2");
      GetHist( label+"mm_eventcutflow")->GetXaxis()->SetBinLabel(10,"SSMMLoose");
      GetHist( label+"mm_eventcutflow")->GetXaxis()->SetBinLabel(11,"MLL");
      GetHist( label+"mm_eventcutflow")->GetXaxis()->SetBinLabel(12,"mmdR");
      GetHist( label+"mm_eventcutflow")->GetXaxis()->SetBinLabel(13,"DiJet");
      GetHist( label+"mm_eventcutflow")->GetXaxis()->SetBinLabel(14,"DiNewJet");
      GetHist( label+"mm_eventcutflow")->GetXaxis()->SetBinLabel(15,"SSMMTight");
      GetHist( label+"mm_eventcutflow")->Fill(cut,weight);

    }
  }
  if(cf==1){
    if(GetHist( label+"cutflow_challenge")) {
      GetHist( label+"cutflow_challenge")->Fill(cut,weight);
      
    }
    else{
      AnalyzerCore::MakeHistograms( label+"cutflow_challenge",11,0.,11.);

      GetHist( label+"cutflow_challenge")->GetXaxis()->SetBinLabel(1,"1_NoCut");
      GetHist( label+"cutflow_challenge")->GetXaxis()->SetBinLabel(2,"2_EventClean");
      GetHist( label+"cutflow_challenge")->GetXaxis()->SetBinLabel(3,"3_DiMu");
      GetHist( label+"cutflow_challenge")->GetXaxis()->SetBinLabel(4,"4_MuVeto");
      GetHist( label+"cutflow_challenge")->GetXaxis()->SetBinLabel(5,"5_ElVeto");
      GetHist( label+"cutflow_challenge")->GetXaxis()->SetBinLabel(6,"6_mll");
      GetHist( label+"cutflow_challenge")->GetXaxis()->SetBinLabel(7,"7_DiJet");
      GetHist( label+"cutflow_challenge")->GetXaxis()->SetBinLabel(8,"8_MET");
      GetHist( label+"cutflow_challenge")->GetXaxis()->SetBinLabel(9,"9_mjj");
      GetHist( label+"cutflow_challenge")->GetXaxis()->SetBinLabel(10,"91_bjetv1");
      GetHist( label+"cutflow_challenge")->GetXaxis()->SetBinLabel(11,"92_bjetv2");

      GetHist( label+"cutflow_challenge")->Fill(cut,weight);
      
    }
  }
  if(cf==-1){
    if(GetHist( label+"noWeightcutflow_challenge")) {
      GetHist( label+"noWeightcutflow_challenge")->Fill(cut,weight);

    }
    else{
      
      AnalyzerCore::MakeHistograms( label+"noWeightcutflow_challenge",11,0.,11.);

      GetHist( label+"noWeightcutflow_challenge")->GetXaxis()->SetBinLabel(1,"NoWeight_1_NoCut");
      GetHist( label+"noWeightcutflow_challenge")->GetXaxis()->SetBinLabel(2,"NoWeight_2_EventClean");
      GetHist( label+"noWeightcutflow_challenge")->GetXaxis()->SetBinLabel(3,"NoWeight_3_DiMu");
      GetHist( label+"noWeightcutflow_challenge")->GetXaxis()->SetBinLabel(4,"NoWeight_4_MuVeto");
      GetHist( label+"noWeightcutflow_challenge")->GetXaxis()->SetBinLabel(5,"NoWeight_5_ElVeto");
      GetHist( label+"noWeightcutflow_challenge")->GetXaxis()->SetBinLabel(6,"NoWeight_6_mll");
      GetHist( label+"noWeightcutflow_challenge")->GetXaxis()->SetBinLabel(7,"NoWeight_7_DiJet");
      GetHist( label+"noWeightcutflow_challenge")->GetXaxis()->SetBinLabel(8,"NoWeight_8_MET");
      GetHist( label+"noWeightcutflow_challenge")->GetXaxis()->SetBinLabel(9,"NoWeight_9_mjj");
      GetHist( label+"noWeightcutflow_challenge")->GetXaxis()->SetBinLabel(10,"NoWeight_10_bjetv1");
      GetHist( label+"noWeightcutflow_challenge")->GetXaxis()->SetBinLabel(11,"NoWeight_11_bjetv2");
      
      
      GetHist( label+"noWeightcutflow_challenge")->Fill(cut,weight);
      
    }
  }
  
  if(cf==2){
    if(GetHist(label +  "mm_eventcutflow_LowMassDoubleMuon")) {
      GetHist(label +  "mm_eventcutflow_LowMassDoubleMuon")->Fill(cut,weight);

    }
    else{
      
      vector<TString> IDs;
      IDs.push_back("LowMassMM_HNTIGHT");
      IDs.push_back("LowMassMM_HNTIGHT2");
      IDs.push_back("LowMassMM_HNTIGHT3");
      IDs.push_back("LowMassMM_HNTIGHT4");
      IDs.push_back("LowMassMM_HNTIGHT5");
      IDs.push_back("LowMassMM_HNTIGHT6");
      AnalyzerCore::MakeHistograms(label +  "mm_eventcutflow_LowMassDoubleMuon",IDs.size(),0.,float(IDs.size()));

      for(unsigned int iid=0; iid < IDs.size(); iid++){

        TString elid = IDs[iid];
        GetHist(label +  "mm_eventcutflow_LowMassDoubleMuon")->GetXaxis()->SetBinLabel(1+iid,elid);
      }
      GetHist(label +  "mm_eventcutflow_LowMassDoubleMuon")->Fill(cut,weight);

    }
  }
  
  if(cf==3){
    if(GetHist( label + "mm_eventcutflow_HighMassDoubleMuon")) {
      GetHist( label + "mm_eventcutflow_HighMassDoubleMuon")->Fill(cut,weight);

    }
    else{

      vector<TString> IDs;
      IDs.push_back("HighMassMM_HNTIGHT");
      IDs.push_back("HighMassMM_HNTIGHT2");
      IDs.push_back("HighMassMM_HNTIGHT3");
      IDs.push_back("HighMassMM_HNTIGHT4");
      IDs.push_back("HighMassMM_HNTIGHT5");
      IDs.push_back("HighMassMM_HNTIGHT6");
      AnalyzerCore::MakeHistograms( label + "mm_eventcutflow_HighMassDoubleMuon",IDs.size(),0.,float(IDs.size()));

      for(unsigned int iid=0; iid < IDs.size(); iid++){

	TString elid = IDs[iid];
	GetHist( label + "mm_eventcutflow_HighMassDoubleMuon")->GetXaxis()->SetBinLabel(1+iid,elid);
      }
      GetHist( label + "mm_eventcutflow_HighMassDoubleMuon")->Fill(cut,weight);

    }
  }
  if(cf==4){
    if(GetHist( label + "mm_lowmassopt_cutflow")) {
      GetHist( label + "mm_lowmassopt_cutflow")->Fill(cut,weight);

    }
    else{

      vector<TString> IDs;
      IDs.push_back("NoCut");
      IDs.push_back("ElVeto");
      IDs.push_back("DiLooseMM");
      IDs.push_back("SSMMLoose");
      IDs.push_back("MLL");
      IDs.push_back("mmdR");
      IDs.push_back("LepVeto");
      IDs.push_back("SSMM_DiJet_noB");
      IDs.push_back("SSMM_FatJet");
      IDs.push_back("SSMM_NoFatJet");
      IDs.push_back("SSMM_1Jet");
      IDs.push_back("SSMM_1JetOpt");
      IDs.push_back("SSMM_1JetOpt2");
      IDs.push_back("SSMM_1JetOpt_pt2");
      IDs.push_back("SSMM_1JetOpt2_pt2");
      IDs.push_back("SSMM_DiJet");      
      IDs.push_back("SSMM_DiJet_passMET45");
      IDs.push_back("SSMM_DiJet_passMET50");
      IDs.push_back("SSMM_DiJet_passMET55");
      IDs.push_back("SSMM_DiJet_passMET60");
      IDs.push_back("SSMM_DiJet_passMET65");
      IDs.push_back("SSMM_DiJet_passMET70");
      IDs.push_back("SSMM_DiJet_passMET80");
      IDs.push_back("SSMM_DiJet_passMETST7_5");
      IDs.push_back("SSMM_DiJet_passMETST10");
      IDs.push_back("SSMM_DiJet_passMETST12_5");
      IDs.push_back("SSMM_DiJet_passMETST15");
      IDs.push_back("SSMM_DiJet_passMETST17_5");
      IDs.push_back("SSMM_DiJet_passMETST20");
      IDs.push_back("SSMM_DiJet_passHT240");
      IDs.push_back("SSMM_DiJet_passHT260");
      IDs.push_back("SSMM_DiJet_passHT280");
      IDs.push_back("SSMM_DiJet_passHT300");
      IDs.push_back("SSMM_DiJet_passDR2");
      IDs.push_back("SSMM_DiJet_passMT80");
      IDs.push_back("SSMM_DiJet_passMT90");
      IDs.push_back("SSMM_DiJet_passMT100");
      IDs.push_back("SSMM_DiJet_passLJ");
      IDs.push_back("SSMM_DiJet_passmll");
      IDs.push_back("SSMM_DiJet_passpt2");
      IDs.push_back("SSMM_DiJet_passpt2mll");


      AnalyzerCore::MakeHistograms( label + "mm_lowmassopt_cutflow",IDs.size(),0.,float(IDs.size()));

      for(unsigned int iid=0; iid < IDs.size(); iid++){

        TString elid = IDs[iid];
        GetHist( label + "mm_lowmassopt_cutflow")->GetXaxis()->SetBinLabel(1+iid,elid);
      }
      GetHist( label + "mm_lowmassopt_cutflow")->Fill(cut,weight);

    }
  }




}

void HNDiLepton::FillEventCutFlow(TString cut, TString label, float weight){

  if(GetHist(label + "_eventcutflow")) {
    GetHist(label + "_eventcutflow")->Fill(cut,weight);

  }
  else{
    AnalyzerCore::MakeHistograms(label + "_eventcutflow",10,0.,10.);

    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(1,"SSMM_loose");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(2,"SSMM_tightiso");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(3,"SSMM_tightiso2");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(4,"SSMM_tightmdiso1");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(5,"SSMM_tightmdiso2");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(6,"SSMM_tightmdiso3");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(7,"SSMM_tightmdiso4");
  }
  
}

void HNDiLepton::GetTriggEfficiency(){
  //ListTriggersAvailable();                                                                                                                                                                                                                                                                                                                                                      
  vector<int> pt1;
  pt1.push_back(18);
  pt1.push_back(35);
  pt1.push_back(25);
  pt1.push_back(25);
  pt1.push_back(25);
  pt1.push_back(30);
  pt1.push_back(30);
  pt1.push_back(20);
  pt1.push_back(25);
  pt1.push_back(30);
  pt1.push_back(120);
  pt1.push_back(20);
  pt1.push_back(25);
  pt1.push_back(25);
  vector<int>pt2;
  pt2.push_back(10);
  pt2.push_back(35);
  pt2.push_back(25);
  pt2.push_back(10);
  pt2.push_back(10);
  pt2.push_back(10);
  pt2.push_back(10);
  pt2.push_back(15);
  pt2.push_back(15);
  pt2.push_back(10);
  pt2.push_back(10);
  pt2.push_back(15);
  pt2.push_back(15);
  pt2.push_back(15);


  std::vector<TString> lists_triggers;
  lists_triggers.push_back("HLT_Ele15_IsoVVVL_PFHT600_v");
  lists_triggers.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v");
  lists_triggers.push_back("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v");
  lists_triggers.push_back("HLT_Ele22_eta2p1_WPLoose_Gsf_v");
  lists_triggers.push_back("HLT_Ele23_WPLoose_Gsf_v");
  lists_triggers.push_back("HLT_Ele27_WPLoose_Gsf_v");
  lists_triggers.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_v");
  lists_triggers.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v");
  lists_triggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
  lists_triggers.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_HT200");
  lists_triggers.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");
  lists_triggers.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  lists_triggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  lists_triggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded_v1");
  for(unsigned int i=0; i < lists_triggers.size(); i++){
    FillTriggerEfficiency(lists_triggers.at(i), weight, "denominator_nojet", lists_triggers );
  }

  std::vector<snu::KJet> jets=GetJets("JET_HN");
  std::vector<snu::KElectron> electrons= GetElectrons("ELECTRON_POG_TIGHT");
  if(jets.size() > 1){
    for(unsigned int i=0; i < lists_triggers.size(); i++){
      TString trig=lists_triggers.at(i);
      FillTriggerEfficiency(lists_triggers.at(i), weight, "denominator", lists_triggers );
      if(PassTrigger(trig))  {
	FillTriggerEfficiency(lists_triggers.at(i), weight, "numerator",lists_triggers );
	
	if(electrons.size() ==2) {
	  
	  FillTriggerEfficiency(lists_triggers.at(i), weight, "numerator_dimuon",lists_triggers );
	  if(electrons.at(0).Pt() > pt1.at(i) && electrons.at(1).Pt() > pt2.at(i))  FillTriggerEfficiency(lists_triggers.at(i), weight, "numerator_dimuon_pt",lists_triggers );
	}
      }
    }
  }
}

void HNDiLepton::FillTriggerEfficiency(TString cut, float weight, TString label, std::vector<TString> list){

  if(GetHist("TriggerEfficiency_" + label)) {
    GetHist("TriggerEfficiency_"+label)->Fill(cut,weight);

  }
  else{
    int ntrig = list.size();
    AnalyzerCore::MakeHistograms("TriggerEfficiency_"+label,ntrig,0.,float(ntrig));

    for(unsigned int it=0; it < list.size(); it++){
      GetHist("TriggerEfficiency_"+label)->GetXaxis()->SetBinLabel(it+1,list.at(it));
    }
  }
  
}


bool HNDiLepton::CheckSignalRegion( bool isss,  std::vector<snu::KMuon> muons, std::vector<snu::KElectron> el,  std::vector<snu::KJet> jets,  std::vector<snu::KJet> alljets, TString name, float w){

  bool debug=false;
  if(el.size() > 0) return false ;
  if(muons.size() != 2 ) {if(debug)cout << "Fail el size" << endl; return false ;}
  // Set by trigger                                                                                                                                                                                                                                                
  if(muons.at(0).Pt() < 20.) {if(debug)cout << "Fail pt1 " << endl; return false;}
  if(muons.at(1).Pt() < 5.)  {if(debug)cout << "Fail pt2  " << endl; return false;}

  if(isss&&!SameCharge(muons)) {if(debug)cout << "Fail ss " << endl; return false;}

  if(!isss&&SameCharge(muons)) {if(debug)cout << "Fail os  " << endl; return false;}
  if(jets.size() < 2) {if(debug)cout << "Fail jets " << endl; return false;}

  snu::KParticle mm = muons.at(0) + muons.at(1);
  if(mm.M()  < _mm_mll_presel_cut) {if(debug)cout << "Fail mee  " << endl; return false;}

  
  float dijetmass_tmp=999.;
  float dijetmass=9990000.;
  int m=-999;
  int n=-999;
  for(UInt_t emme=0; emme<jets.size(); emme++){
    for(UInt_t enne=1; enne<jets.size(); enne++) {
      if(emme == enne) continue;
      if(name.Contains("Low")){
	dijetmass_tmp = (jets[emme]+jets[enne]+muons[0] + muons[1]).M();
      }
      else  dijetmass_tmp = (jets[emme]+jets[enne]).M();
      if ( fabs(dijetmass_tmp-80.4) < fabs(dijetmass-80.4) ) {
        dijetmass = dijetmass_tmp;
        m = emme;
        n = enne;
      }
    }
  }
  

  snu::KParticle jj = jets.at(m) + jets.at(n) ;

  float dPhi = fabs(TVector2::Phi_mpi_pi(jets[m].Phi() - jets[n].Phi()));
  float contramass=2*jets[m].Pt()*jets[n].Pt()*(1+cos(dPhi));
  contramass=sqrt(contramass);

  if(name.Contains("Low")){
    if((jets[m] + jets[n] + muons[0] + muons[1]).M() > 300.)  return false;
  }
  else{
    if(jj.M() > 150.) {if(debug)cout << "Fail mjj  " << endl; return false;}

  }
  if(debug)cout << "PASSES ID" << endl;
  return true;
}




bool HNDiLepton::OptMassCheckSignalRegion(std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, float evmet, bool runchargeflip , std::vector<float> cuts , TString additional_option){
  
  return OptMassCheckSignalRegion(electrons, jets, 1, evmet, runchargeflip, cuts, additional_option);
}

bool HNDiLepton::OptMassCheckSignalRegion(std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, int btagwp,  float evmet, bool runchargeflip , std::vector<float> cuts, TString additional_option){
  
  return true;

}


bool HNDiLepton::LowMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runchargeflip){
 

  return true;
  
}

bool HNDiLepton::MidMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runchargeflip){
  
 
  return true;
}



	
	

bool HNDiLepton::HighMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runchargeflip){
 

  return true;
}




void HNDiLepton::counter(TString cut, float w){

  map<TString,float>::iterator itmapcounter = mapcounter.find(cut) ;
  if (itmapcounter == mapcounter.end()){
    mapcounter[cut] = w;
  }
  else{
    float sum = itmapcounter->second;
    mapcounter[cut] = sum+w;
  }
  FillEventCutFlow(1,cut,w);
  //FillEventCutFlow(-1,cut,1);

}


void HNDiLepton::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger << DEBUG << "END OF CYCLE: isdata=" << isData <<  LQLogger::endmsg;

  for ( map<TString,float>::iterator itmapcounter = mapcounter.begin(); itmapcounter != mapcounter.end(); itmapcounter++){
    cout << itmapcounter->first << " has count "<< itmapcounter->second << endl;
  }

}


void HNDiLepton::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  

  InitialiseAnalysis();
  
  /// If ConfigureFake is used then no fake histograms are setup in HNCommonFake code so you MUST setup in your anlayis code if you wish to use this code

  if(configure_fakes)ConfigureFake();
  

  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //  DeclareVariable(out_muons, "Signal_Muons");
  
  return;
  
}

HNDiLepton::~HNDiLepton() {
  
  Message("In HNDiLepton Destructor" , INFO);
  
}



void HNDiLepton::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void HNDiLepton::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNDiLeptonCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HNDiLepton::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



