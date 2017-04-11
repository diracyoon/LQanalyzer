#include "SKTreeFiller.h"
#include <stdio.h>  

#include <stdlib.h>
#include <iostream>

using namespace snu;
using namespace std;


SKTreeFiller::SKTreeFiller() :Data() {
  
  TString fitParametersFile = "";
};


SKTreeFiller::~SKTreeFiller() {};


bool SKTreeFiller::SkipTrigger(TString tname){
  
  m_logger << DEBUG << "Trigger: " << tname << LQLogger::endmsg;  
  /// Remove extra unnecisary  triggers (from v7-6-4+ this will not be needed))
  if((tname.Contains("Jpsi")
       || tname.Contains("NoFilters")
       || tname.Contains("Upsilon")
       || tname.Contains("7p5")
       || tname.Contains("Save")
       || tname.Contains("R9Id")
       || tname.Contains("PFMET")
       || tname.Contains("PFHT")
       || tname.Contains("NoHE")
       || tname.Contains("HE10")
       || tname.Contains("PFJet50")
       || tname.Contains("Boost")
       || tname.Contains("LooseIso")
       || tname.Contains("MediumIso")
       || tname.Contains("Mass")
       || tname.Contains("Central")
       || tname.Contains("MW")
       || tname.Contains("EBOnly_VBF")
       || tname.Contains("dEta18"))) return true;
  
  return false;
}


snu::KTrigger SKTreeFiller::GetTriggerInfo(std::vector<TString> trignames){
  snu::KTrigger ktrigger;

  if(!LQinput){
    ktrigger = *k_inputtrigger;
    return ktrigger;
  }
  m_logger << DEBUG << "Filling trigger Info" << LQLogger::endmsg;


  std::vector<std::string> vHLTInsideDatasetTriggerNames;
  std::vector<bool> vHLTInsideDatasetTriggerDecisions;
  std::vector<int> vHLTInsideDatasetTriggerPrescales;
  

  /// trignames should only be empty id user is running on Catuples and not SKTreeMaker. In this case all triggers are used 
  if(trignames.size() == 0 ){
    for (UInt_t i=0; i< vtrignames->size(); i++) {
      std::string tgname = vtrignames->at(i);
      Int_t ps = vtrigps->at(i);
      vHLTInsideDatasetTriggerNames.push_back(tgname);
      if(ps > 0) vHLTInsideDatasetTriggerDecisions.push_back(true);
      else vHLTInsideDatasetTriggerDecisions.push_back(false);
      vHLTInsideDatasetTriggerPrescales.push_back(ps);
    }
  }

  
  /// vtrigname is vector of ALL triggers in Catuples
  for (UInt_t i=0 ; i< vtrignames->size(); i++) {
    // trignames is vector of trigger names that we want to store in SKTrees
    // trigname contains names substrings X (where X is for example "HLT_mu") and we store all triggers that start with X
    
    
    std::string tgname = vtrignames->at(i);
    if(TString(CatVersion).Contains("v7-6-2")) {
      if(SkipTrigger(TString(tgname)))continue;
    }

    Int_t ps = vtrigps->at(i);

    for (std::vector<TString>::reverse_iterator it (trignames.end());
	 it != std::vector<TString>::reverse_iterator (trignames.begin());
	 ++it) {

      TString tmpHLT = vtrignames->at(i);
      if ( tmpHLT.BeginsWith(*it)){
	
	vHLTInsideDatasetTriggerNames.push_back(tgname);
	if(ps > 0) vHLTInsideDatasetTriggerDecisions.push_back(true);
	else vHLTInsideDatasetTriggerDecisions.push_back(false);
	vHLTInsideDatasetTriggerPrescales.push_back(ps);
	
	// if trigger is accepted break from loop
	break;
      }
    } // end of trignames loop
  }// loop of all triggers  
  
  ktrigger.SetHLTInsideDatasetTriggerNames(vHLTInsideDatasetTriggerNames);
  ktrigger.SetHLTInsideDatasetTriggerDecisions(vHLTInsideDatasetTriggerDecisions);
  ktrigger.SetHLTInsideDatasetTriggerPrescales(vHLTInsideDatasetTriggerPrescales);
    
  return ktrigger;
  
}

snu::KEvent SKTreeFiller::GetEventInfo(){
 
  snu::KEvent kevent;

  if(!LQinput){
    kevent = *k_inputevent;
    if(k_cat_version < 3){
      if(!TString(kevent.CatVersion()).Contains("v7-4"))kevent.SetCatVersion(CatVersion);
    }
    return kevent;
  }
  //  lumimask = snu::KEvent::gold

  m_logger << DEBUG << "Filling Event Info" << LQLogger::endmsg;
  
  // New variable to set catversion. Add this to flat ntuples for next iteration
  kevent.SetCatVersion(CatVersion);

  if(k_cat_version > 7)  {
    
    /// type 1
    double met_type1 =  sqrt(met_jetRes_Px_up->at(0)*met_jetRes_Px_up->at(0) + met_jetRes_Py_up->at(0)*met_jetRes_Py_up->at(0));
    double phi_type1 =  TMath::ATan2(met_jetRes_Py_up->at(0),met_jetRes_Px_up->at(0)); 
    // type 1 + ohi corrections
    double met_type1xy = sqrt(met_xyshift_px->at(0)*met_xyshift_px->at(0) + met_xyshift_py->at(0)*met_xyshift_py->at(0));
    double phi_type1xy =  TMath::ATan2(met_xyshift_py->at(0), met_xyshift_px->at(0));


    /// Default MET is now xy shifted typ1
    if(isData)  {
      kevent.SetMET(snu::KEvent::pfmet, met_type1xy, phi_type1xy, met_xyshift_sumet->at(0));
      kevent.SetPFMETx(met_xyshift_px->at(0));
      kevent.SetPFMETy(met_xyshift_py->at(0));

      /// Also for completness store type1 without phi corrections
      kevent.SetPFMETType1x(met_jetRes_Px_up->at(0));
      kevent.SetPFMETType1y(met_jetRes_Py_up->at(0));
      kevent.SetPFMETType1SumEt(met_sumet->at(0));
    }
    /// set unsmeared met variables
    kevent.SetPFMETType1Unsmearedx(met_jetRes_Px_up->at(0));
    kevent.SetPFMETType1Unsmearedy(met_jetRes_Py_up->at(0));
    kevent.SetPFMETType1xyUnsmearedx(met_xyshift_px->at(0));
    kevent.SetPFMETType1xyUnsmearedy(met_xyshift_py->at(0));
    
  }
  m_logger << DEBUG << "Filling Event Info [2]" << LQLogger::endmsg;
  /// Since some versions of catuples have no metNoHF due to bug in met code 


  if(k_cat_version > 3){
    /// k_cat_version > 3 == v765+
    if(PDFWeights){
      if(PDFWeights->size() > 0){
	kevent.SetPDFWeights(*PDFWeights);
      }
    }
    if(ScaleWeights){
      if(ScaleWeights->size() > 0){
      kevent.SetScaleWeights(*ScaleWeights);
      }
    }

  }

  

  if(!isData){
    float jpx(0.), jpy(0.), sjpx(0.), sjpy(0.), sjpxup(0.), sjpxdown(0.),sjpyup(0.), sjpydown(0.) ;

    /// only smear jets not close to leptons (use top projection id)
    for(unsigned int ij = 0 ; ij < jets_pt->size(); ij++){
      bool close_to_lepton(false);
      if(jets_pt->at(ij) < 10.) continue;
      for(unsigned int im=0; im < muon_pt->size(); im++){
	if(muon_pt->at(im) < 10.) continue;
	if(fabs(muon_eta->at(im)) > 2.5) continue;
	// find full definition for 13 TeV
	//if(muon_relIso04->at(im) > 0.2)  continue;
        double dr = sqrt( pow(fabs( jets_eta->at(ij) - muon_eta->at(im)),2.0) +  pow( fabs(TVector2::Phi_mpi_pi( jets_phi->at(ij) - muon_phi->at(im))),2.0));
	if(dr < 0.4){
	  close_to_lepton=true;
	}
      }
      for(unsigned int iel=0; iel < electrons_pt->size(); iel++){
	if(electrons_pt->at(iel) < 10.) continue;
        if(fabs(electrons_eta->at(iel)) > 2.5) continue;
	// find full definition for 13 TeV                                                                                                                                          if(electrons_relIso03->at(ilep) > 0.15)  continue;
        double dr = sqrt( pow(fabs( jets_eta->at(ij) - electrons_eta->at(iel)),2.0) +  pow( fabs(TVector2::Phi_mpi_pi( jets_phi->at(ij) - electrons_phi->at(iel))),2.0));
        if(dr < 0.4){
          close_to_lepton=true;
        }
      }
      
      if(close_to_lepton) continue;
      
      float jets_px = jets_pt->at(ij) *TMath::Cos(jets_phi->at(ij)); 
      float jets_py = jets_pt->at(ij) *TMath::Sin(jets_phi->at(ij));
      jpx +=  jets_px;
      jpy +=  jets_py;
      
      if(!isData){
	sjpx +=  jets_smearedRes->at(ij) *jets_px;
	sjpy +=  jets_smearedRes->at(ij) *jets_py;
      }
      else{
	sjpx +=  jets_px;
        sjpy +=  jets_py;
      }
      sjpxup +=  jets_smearedResUp->at(ij) *jets_px;
      sjpyup +=  jets_smearedResUp->at(ij) *jets_py;
      
      sjpxdown +=  jets_smearedResDown->at(ij) *jets_px;
      sjpydown +=  jets_smearedResDown->at(ij) *jets_py;

    }

    // met_jetRes_Px_up ==met_Px since no smearing is applied in miniaods -> cattools
    float met_x  = met_xyshift_px->at(0)  +  jpx - sjpx;
    float met_y  = met_xyshift_py->at(0)  +  jpy - sjpy;
    float met_newpt = sqrt(met_x*met_x+ met_y*met_y);
    float met_newphi = TMath::ATan2(met_y,met_x);
    
    //cout << met_newphi << "  phi " << met_xyshift_phi->(0) << endl;
    kevent.SetMET(snu::KEvent::pfmet,  met_newpt,met_newphi,  met_xyshift_sumet->at(0));  
    kevent.SetPFMETx(met_x);
    kevent.SetPFMETy(met_y);

    /// correct MET for jets smearing
    float type1_met_x  = met_jetRes_Px_up->at(0)  +  jpx - sjpx;
    float type1_met_y  = met_jetRes_Py_up->at(0)  +  jpy - sjpy;
    float type1_met_newpt = sqrt(type1_met_x*type1_met_x+ type1_met_y*type1_met_y);
    float type1_met_newphi = TMath::ATan2(type1_met_y,type1_met_x);
    
    kevent.SetPFMETType1x(type1_met_x);
    kevent.SetPFMETType1y(type1_met_y);		  
    kevent.SetPFMETType1SumEt(met_sumet->at(0));           

    /// Fix met phi 
    float met_x_jer_up  = met_xyshift_px->at(0) +  jpx - sjpxup;
    float met_y_jer_up   = met_xyshift_py->at(0)  +  jpy - sjpyup;
    float met_newpt_jerup = sqrt(met_x_jer_up*met_x_jer_up+ met_y_jer_up*met_y_jer_up);
    float met_x_jer_down   = met_xyshift_px->at(0)  +  jpx - sjpxdown;
    float met_y_jer_down  = met_xyshift_py->at(0)  +  jpy -sjpydown;
    float met_newpt_jerdown = sqrt(met_x_jer_down*met_x_jer_down+ met_y_jer_down*met_y_jer_down);

      
    kevent.SetPFMETShift  (snu::KEvent::up,     snu::KEvent::JetRes,     met_newpt_jerup);
    kevent.SetPFMETShift  (snu::KEvent::down,   snu::KEvent::JetRes,     met_newpt_jerdown);
    

  }
  

  m_logger << DEBUG << "Filling Event Info [3]" << LQLogger::endmsg;
  
  if(k_cat_version > 2){
    if(met_unclusteredEn_Px_up){
      if(met_unclusteredEn_Px_up->at(0)){
	
	// catools use slimmedMETs in MiniAOD
	// this uses type 1 corrected MET
	// as explained in here https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#ETmiss
	// The type1 corrections is computed from ak4PFJetsCHS jets with pT > 15 GeV,
	// smearing mc jets
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
	// change in jet pt needs propagating into MET
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JERCReference
	/// For details see hete
	/// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETRun2Corrections 
	kevent.SetPFMETShift  (snu::KEvent::up,     snu::KEvent::MuonEn,     sqrt(met_muonEn_Px_up*met_muonEn_Px_up + met_muonEn_Py_up*met_muonEn_Py_up));
	kevent.SetPFMETShift  (snu::KEvent::down,   snu::KEvent::MuonEn,     sqrt(met_muonEn_Px_down*met_muonEn_Px_down + met_muonEn_Py_down*met_muonEn_Py_up));
	kevent.SetPFMETShift  (snu::KEvent::up,     snu::KEvent::ElectronEn, sqrt(met_electronEn_Px_up*met_electronEn_Px_up + met_electronEn_Py_up*met_electronEn_Py_up));
	kevent.SetPFMETShift  (snu::KEvent::down,   snu::KEvent::ElectronEn, sqrt(met_electronEn_Px_down*met_electronEn_Px_down + met_electronEn_Py_down*met_electronEn_Py_down));
	kevent.SetPFMETShift  (snu::KEvent::up,     snu::KEvent::Unclustered,sqrt(met_unclusteredEn_Px_up->at(0)*met_unclusteredEn_Px_up->at(0) + met_unclusteredEn_Py_up->at(0)*met_unclusteredEn_Py_up->at(0)));
	kevent.SetPFMETShift  (snu::KEvent::down,   snu::KEvent::Unclustered,sqrt(met_unclusteredEn_Px_down->at(0)*met_unclusteredEn_Px_down->at(0) + met_unclusteredEn_Py_down->at(0)*met_unclusteredEn_Py_down->at(0)));
	kevent.SetPFSumETShift(snu::KEvent::up,     snu::KEvent::Unclustered,met_unclusteredEn_SumEt_up->at(0));
	kevent.SetPFSumETShift(snu::KEvent::down,   snu::KEvent::Unclustered,met_unclusteredEn_SumEt_down->at(0));
	kevent.SetPFMETShift  (snu::KEvent::up,     snu::KEvent::JetEn,      sqrt(met_jetEn_Px_up->at(0)*met_jetEn_Px_up->at(0) + met_jetEn_Py_up->at(0)*met_jetEn_Py_up->at(0)));
	kevent.SetPFMETShift  (snu::KEvent::down,   snu::KEvent::JetEn,      sqrt(met_jetEn_Px_down->at(0)*met_jetEn_Px_down->at(0) + met_jetEn_Py_down->at(0)*met_jetEn_Py_down->at(0)));
	kevent.SetPFSumETShift(snu::KEvent::up,     snu::KEvent::JetEn,      met_jetEn_SumEt_up->at(0));
	kevent.SetPFSumETShift(snu::KEvent::down,   snu::KEvent::JetEn,      met_jetEn_SumEt_down->at(0));

	/// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/python/patPFMETCorrections_cff.py
	/// jets > 15 GeV in mc smeared. This is not done in cattolls so branches have no change,
	/// Apply this here
	/// 

      }
    }
  }
  m_logger << DEBUG << "Filling Event Info [4]" << LQLogger::endmsg;
  
  /// Filling event variables
    
  kevent.SetIsData(isData);
  kevent.SetRunNumber(run);
  kevent.SetEventNumber(event);
  kevent.SetLumiSection(lumi);
  
  if(!isData){
  
    kevent.SetPUWeight(snu::KEvent::central,double(puWeightGold));
    kevent.SetPUWeight(snu::KEvent::down,double(puWeightGoldDn));
    kevent.SetPUWeight(snu::KEvent::up,  double(puWeightGoldUp));
    if(k_cat_version == 4){
      if(puWeightGold_xs71000){
	kevent.SetAltPUWeight(snu::KEvent::central,double(puWeightGold_xs71000));
	kevent.SetAltPUWeight(snu::KEvent::down,double(puWeightGoldDn_xs71000));
	kevent.SetAltPUWeight(snu::KEvent::up,  double(puWeightGoldUp_xs71000));
      }
    }
    else{
      kevent.SetAltPUWeight(snu::KEvent::central,double(puWeightGold));
      kevent.SetAltPUWeight(snu::KEvent::down,double(puWeightGoldDn));
      kevent.SetAltPUWeight(snu::KEvent::up,  double(puWeightGoldUp));

      if(k_cat_version > 7){
	kevent.SetPeriodPileupWeight(double(puWeightGoldB),double(puWeightGoldC),double(puWeightGoldD),double(puWeightGoldE),double(puWeightGoldF),double(puWeightGoldG),double(puWeightGoldH));
      }
    }
  }
  if(isData){
    if(k_cat_version > 2&&k_cat_version < 5){
      kevent.SetLumiMask(snu::KEvent::silver, lumiMaskSilver);
      kevent.SetLumiMask(snu::KEvent::gold,   lumiMaskGold);
    }
    else{
      kevent.SetLumiMask(snu::KEvent::silver, 1);
      kevent.SetLumiMask(snu::KEvent::gold,   1);
    }
  }
  kevent.SetGenId(genWeight_id1, genWeight_id2);
  kevent.SetLHEWeight(lheWeight);
  kevent.SetGenX(genWeightX1, genWeightX2);
  kevent.SetGenQ(genWeightQ);
  if(genWeight > 0.) kevent.SetWeight(1.);
  else kevent.SetWeight(-1.);
  
  kevent.SetVertexInfo(vertex_X, vertex_Y, vertex_Z,0. );
  
  /// MET filter cuts/checks

  
  /// 
  kevent.SetPileUpInteractionsTrue(nTrueInteraction);
    
  kevent.SetNVertices(nPV);
  kevent.SetNGoodVertices(nGoodPV);
  
  kevent.SetIsGoodEvent(nGoodPV);

  /// MET filter cuts/checks
  if(k_cat_version > 4){
    kevent.SetPassEcalDeadCellTriggerPrimitiveFilter(ecalDCTRFilter);
    kevent.SetPassHBHENoiseFilter(HBHENoiseFilter);
    kevent.SetPassHBHENoiseIsoFilter(HBHENoiseIsoFilter);
    kevent.SetPassCSCHaloFilterTight(csctighthaloFilter);
    kevent.SetPassBadEESupercrystalFilter(eeBadScFilter);
    kevent.SetPassTightHalo2016Filter(Flag_globalTightHalo2016Filter);
  }
  else{
    kevent.SetPassEcalDeadCellTriggerPrimitiveFilter(ecalDCTRFilter);
    kevent.SetPassHBHENoiseFilter(HBHENoiseFilter);
    kevent.SetPassCSCHaloFilterTight(csctighthaloFilter);
    kevent.SetPassBadEESupercrystalFilter(eeBadScFilter);
  }
  if(k_cat_version > 6){
    kevent.SetPassBadChargedCandidateFilter(BadChargedCandidateFilter);
    kevent.SetPassBadPFMuonFilter(BadPFMuonFilter);
  }
  return kevent;
}


std::vector<KPhoton> SKTreeFiller::GetAllPhotons(){

  std::vector<KPhoton> photons;

  if(k_cat_version < 3) return photons;
  if(k_cat_version > 4) return photons;
  
  if(!LQinput){
    for(std::vector<KPhoton>::iterator kit  = k_inputphotons->begin(); kit != k_inputphotons->end(); kit++){
      photons.push_back(*kit);
    }
    return photons;
  }
  for (UInt_t iph=0; iph< photons_eta->size(); iph++) {
    if(photons_pt->at(iph) != photons_pt->at(iph)) continue;
    KPhoton ph;
    
    ph.SetPtEtaPhiE(photons_pt->at(iph),photons_eta->at(iph), photons_phi->at(iph),photons_energy->at(iph));

    ph.SetIsLoose(photons_photonID_loose->at(iph));
    ph.SetIsMedium(photons_photonID_medium->at(iph));
    ph.SetIsTight(photons_photonID_tight->at(iph));
    ph.SetPassMVA(photons_photonID_mva->at(iph));
    ph.SetMCMatched(photons_mcMatched->at(iph));
    ph.SetHasPixSeed(photons_haspixseed->at(iph));
    ph.SetPassElVeto(photons_passelectronveto->at(iph));

    ph.SetChargedHadIsoNoEA(photons_chargedHadronIso->at(iph));
    ph.SetpuChargedHadIsoNoEA(photons_puChargedHadronIso->at(iph));
    ph.SetNeutalHadIsoNoEA(photons_neutralHadronIso->at(iph));
    ph.SetPhotonIsoNoEA(photons_photonIso->at(iph));
    ph.SetRhoIso(photons_rhoIso->at(iph));
    ph.SetChargedHadIso(photons_chargedHadronIsoWithEA->at(iph));
    ph.SetPhotonIso(photons_photonIsoWithEA->at(iph));
    ph.SetNeutalHadIso(photons_neutralHadronIsoWithEA->at(iph));
    ph.SetSigmaIetaIeta(photons_sigmaietaieta->at(iph));
    ph.SetR9(photons_r9->at(iph));
    ph.SetHoverE(photons_hovere->at(iph));
    ph.SetSCEta(photons_sceta->at(iph));
    ph.SetSCPhi(photons_scphi->at(iph));
    ph.SetSCRawE(photons_scrawenergy->at(iph));
    ph.SetSCPreShowerE(photons_scpreshowerenergy->at(iph));
    
    photons.push_back(ph);
  }
  std::sort( photons.begin(), photons.end(), isHigherPt );

  return photons;

}

std::vector<KElectron> SKTreeFiller::GetAllElectrons(){

  std::vector<KElectron> electrons;

  if(!LQinput){
    for(std::vector<KElectron>::iterator kit  = k_inputelectrons->begin(); kit != k_inputelectrons->end(); kit++){
      electrons.push_back(*kit);
    }
    return electrons;
  }

  m_logger << DEBUG << "Filling electron Info " << electrons_eta->size() << LQLogger::endmsg;
  
  vector<int> matched_truth;
  for (UInt_t iel=0; iel< electrons_eta->size(); iel++) {
    
    if(electrons_pt->at(iel) != electrons_pt->at(iel))    continue;
    
    KElectron el;

    /// Kinematic Variables
    el.SetPtEtaPhiE(electrons_pt->at(iel),electrons_eta->at(iel), electrons_phi->at(iel),electrons_energy->at(iel));

    el.SetTrigMatch(electron_trigmatch->at(iel));
    el.SetSCEta(electrons_scEta->at(iel));
   
    el.Setdz( electrons_dz->at(iel));
    el.Setdxy(electrons_dxy->at(iel) );
    if(electrons_sigdxy){
      if(electrons_sigdxy->size() > 0 )el.Setdxy_sig(electrons_sigdxy->at(iel) );
    }
    el.SetPFChargedHadronIso(0.3, electrons_puChIso03->at(iel));
    el.SetPFPhotonIso(0.3,electrons_phIso03->at(iel));
    el.SetPFNeutralHadronIso(0.3,electrons_nhIso03->at(iel));
    el.SetPFRelIso(0.3,electrons_relIso03->at(iel));



    m_logger << DEBUG << "Filling electron_minirelIso " << LQLogger::endmsg;
    if(electrons_minirelIso) el.SetPFRelMiniIso(electrons_minirelIso->at(iel));
    
    m_logger << DEBUG << "Filling electron Info 2" << LQLogger::endmsg;
    
    el.SetPFChargedHadronIso(0.4,electrons_puChIso04->at(iel));
    el.SetPFPhotonIso(0.4,electrons_phIso04->at(iel));
    el.SetPFNeutralHadronIso(0.4,electrons_nhIso04->at(iel));
    el.SetPFRelIso(0.4,electrons_relIso04->at(iel));
    
    el.SetPFAbsIso(0.3,electrons_absIso03->at(iel));
    el.SetPFAbsIso(0.4,electrons_absIso04->at(iel));


    /// set Charge variables
    el.SetCharge(electrons_q->at(iel));
    el.SetGsfCtfScPixCharge(electrons_isGsfCtfScPixChargeConsistent->at(iel));
    
    m_logger << DEBUG << "Filling electron Info 3" << LQLogger::endmsg;
    /// set conversion variables
    
    if(electrons_shiftedEnDown){
      el.SetShiftedEUp(electrons_shiftedEnUp->at(iel));
      el.SetShiftedEDown(electrons_shiftedEnDown->at(iel));
    }

    el.SetSNUID(electrons_electronID_snu->at(iel));
    el.SetPassVeto(electrons_electronID_veto->at(iel));
    el.SetPassLoose(electrons_electronID_loose->at(iel));
    el.SetPassMedium(electrons_electronID_medium->at(iel));
    el.SetPassTight(electrons_electronID_tight->at(iel));
    
    /// HEEP
    //el.SetPassHEEP(electrons_electronID_heep->at(iel));

    // MVA
    el.SetPassMVATrigMedium(electrons_electronID_mva_trig_medium->at(iel));
    el.SetPassMVATrigTight(electrons_electronID_mva_trig_tight->at(iel));
    el.SetPassMVANoTrigMedium(electrons_electronID_mva_medium->at(iel));
    el.SetPassMVANoTrigTight(electrons_electronID_mva_tight->at(iel));

    el.SetIsPF(electrons_isPF->at(iel));
    if(electrons_isTrigMVAValid) el.SetIsTrigMVAValid(electrons_isTrigMVAValid->at(iel));
    //el.SetIsMCMatched(electrons_mcMatched->at(iel));
    el.SetHasMatchedConvPhot(electrons_passConversionVeto->at(iel));
    
    el.SetTrkVx(electrons_x->at(iel));
    el.SetTrkVy(electrons_y->at(iel));
    el.SetTrkVz(electrons_z->at(iel));
    m_logger << DEBUG << "Filling electron Info 4" << LQLogger::endmsg;    

    //// Set Is ChargeFlip
    bool isprompt= false;
    bool from_tau = false;
    
    int mother_index=-1;
    int mother_pdgid=-1;
    int matched_index=-1;
    int mc_pdgid=-1;
    bool matched_in_Dr=false;

    int           eltype=0;
    
    if(k_cat_version  > 3){
      
      if(gen_pt){
	// Default deltaR setting for matching
	float min_Dr=0.1;
	/// Loop over all gen particles
	for (UInt_t it=0; it< gen_pt->size(); it++ ){
	  
	  
	  /// Requirements to make sure no crash or warnings with pt=0
	if(gen_motherindex->at(it) <= 0)continue;
	if(gen_motherindex->at(it) >= int(gen_pt->size()))continue;
	if(gen_pt->at(it) < 0.1) continue;
	
	
	double match_eta =electrons_eta->at(iel);
	double match_phi =electrons_phi->at(iel);
	double dr = sqrt( pow(fabs( match_eta - gen_eta->at(it)),2.0) +  pow( fabs(TVector2::Phi_mpi_pi( match_phi - gen_phi->at(it))),2.0));
	
	/// Matching using instructions on
	/// https://indico.cern.ch/event/292928/contributions/1650088/attachments/547844/755123/talk_electron_contribution.pdf
	/// 
	
	/// Match required to status 1 electron
	if(gen_status->at(it) != 1) continue;
	if(fabs(gen_pdgid->at(it)) != 11) continue;
	if(gen_pt->at(it) < 5.) continue;

	/// Check status 1 electron is not matched already to areco electron
	bool already_matched=false;
	for(unsigned int im=0; im < matched_truth.size();im++){
          if(it == unsigned(matched_truth.at(im))) already_matched=true;
        }
        if(already_matched) continue;
	  
	if(matched_in_Dr){
	  /// This is for multiple matched status 1 el.
	  /// In case multiple status 1 electrons are matched with same mother check pt
	  if(gen_motherindex->at(it) == gen_motherindex->at(matched_index)){
	    if (dr < 0.1){
	      // in case 2+ electrons from same mother electron (conversion) also match in pt
	      if( fabs(gen_pt->at(it)-electrons_pt->at(iel)) < fabs(gen_pt->at(matched_index)-electrons_pt->at(iel))) matched_index=it;
	    }
	  }
	  else if ((dr < min_Dr) ){
	    
	    /// find closest match in dR to status 1
	    matched_in_Dr=true;
	    min_Dr= dr;
	    
	    /// set index of matched status 1 electron
	    matched_index=it;
	  }
	}
	else{
	  /// first match status 1 electron
	  if ((dr < min_Dr) ){
	    
	    /// find closest match in dR to status 1
	    matched_in_Dr=true;
	    min_Dr= dr;
	    
	    /// set index of matched status 1 electron
	    matched_index=it;
	    
	  }
	}
      }// end of gen loop to find status 1 electron


      ///// treat case where there is a matched status 1 electron:
      //// classify into prompt:Fake:FromTau

      if(matched_in_Dr){
	/// Find closest non electron ancesteror
	float pdgid = gen_pdgid->at(matched_index);

	// mc_pdgid = closest matched status 1 pdgid
	mc_pdgid= int(pdgid);

	// mindex = mother index: will loop to find first non el mother
	int mindex= matched_index;

	while ( (fabs(gen_pdgid->at(mindex)) == 11)) {
	  pdgid = gen_pdgid->at(mindex);
	  mindex=gen_motherindex->at(mindex);
	}

	/// pdgid is now of electron from non electron mother
	//  mindex = index for mother of non electron ancestor
	
	if( (fabs(gen_pdgid->at(mindex)) == 23) || (fabs(gen_pdgid->at(mindex)) == 24)) {
	  /// Check if el from Z/W is CF and if it is from a photon conversion
	  
	  eltype=1;
	  int n_el_from_el=0;
	  float charge_sum=0.;
	  /// Loop over electrons: Find mother of matched status 1 and see what other daughters there are:
	  /// In case of a conversion i.e  Z->ee->eephoton->eeee the status 23 electorn decays to 3 electrons e+e+e- or e-e-e+
	  bool isthirdel_fromconv(false);
	  for (UInt_t itx=0; itx< gen_pt->size(); itx++ ){
	    if(gen_motherindex->at(itx) <= 0)continue;
	    if(gen_motherindex->at(itx) >= int(gen_pt->size()))continue;
	    if(gen_pt->at(itx) < 0.1) continue;
	    if(fabs(gen_pdgid->at(itx)) ==11) {
	      if(gen_motherindex->at(itx) == gen_motherindex->at(matched_index)) { 
		charge_sum+= gen_pdgid->at(itx); n_el_from_el++;
		if(n_el_from_el==3){
		  if(itx == matched_index) isthirdel_fromconv=true;
		}
	      }

	    }
	  }

	  /// Set if conversion i.ei e->eee
	  /// Two methods: 
	  /// 1) check pdgid of status 1 el vs mother. if < 0 it is a converison
	  /// 2) In case closest status 1 el is not opposite charge truth check number of electrons from mother if 3 it is a conversion
	  if((gen_pdgid->at(matched_index)  * pdgid) < 0 )  {el.SetIsPhotonConversion(true);           eltype=2;} 
	  else  el.SetIsPhotonConversion(false);
	  
	  if(!isthirdel_fromconv){
	    if(n_el_from_el ==3&& (fabs(charge_sum) == 11)) { eltype=3; el.SetIsPhotonConversion(true); }
	  }
	  else{
	    if(pdgid * electrons_q->at(iel) > 0 )  {
	      if(n_el_from_el ==3&& (fabs(charge_sum) == 11)) { eltype=3; el.SetIsPhotonConversion(true);}
	    }
	  }

	  /// Check if it is a chargeflip.
	  /// Either from a conversion or just reconstructed charge is wrong
	  if(pdgid * electrons_q->at(iel) > 0 )    
	    { el.SetIsChargeFlip(true); 
	      if(eltype == 2 || eltype == 3){
		if(eltype == 2) eltype=4;
		if(eltype == 3) eltype=5;
	      }
	      else eltype=6;
	    }

	  else     el.SetIsChargeFlip(false);
	  
	  mother_index=mindex;
	  mother_pdgid=gen_pdgid->at(mindex);
	  isprompt=true; /// means is prompt
	  
	}/// end of Z/W
	else {
	  if(gen_status->at(mindex) == 2){
	    if(fabs(gen_pdgid->at(mindex)) > 50) {isprompt=false; mother_pdgid=gen_pdgid->at(mindex); mother_index=mindex; from_tau=false;
	      eltype=7;
	      
	      if(gen_isprompt->at(matched_index)){
		cout << "matched FAKE, but isPrompt flag??" << endl;
		cout << "------------------CF "<< endl;
		cout << "gen_isprompt = " << gen_isprompt->at(matched_index)  << endl;
		cout << "gen_isdecayedleptonhadron = " <<gen_isdecayedleptonhadron->at(matched_index)  << endl;
		cout << "gen_isdirecthadrondecayproduct  = " <<gen_isdirecthadrondecayproduct->at(matched_index)  << endl;
		cout << "gen_ishardprocess  = " << gen_ishardprocess->at(matched_index)  << endl;
		cout << "gen_istaudecayproduct =  " << gen_istaudecayproduct->at(matched_index)  << endl;
		cout << "gen_isprompttaudecayproduct =  " <<  gen_isprompttaudecayproduct->at(matched_index)  << endl;
	      }
	      
	    }
	    else {
	      isprompt=true;
	      mother_pdgid=gen_pdgid->at(mindex); mother_index=mindex; from_tau=false;
	      eltype=8;
	      
	      
	      if(fabs(gen_pdgid->at(mindex)) == 22){
		if(fabs(gen_pdgid->at(gen_motherindex->at(mindex))) > 50){
		  eltype=9;
		}
		else {
		  eltype=10;
		}
	      }
	    }
	    
	    if(fabs(gen_pdgid->at(mindex)) == 15){
	      eltype=11;

	      isprompt=true; mother_pdgid=gen_pdgid->at(mindex);  mother_index=mindex; from_tau=true;
	      // Check if el from tau  is CF
	      
	      int n_el_from_el=0;
	      float charge_sum=0.;
	      for (UInt_t itx=0; itx< gen_pt->size(); itx++ ){
		if(itx == matched_index) continue;
		if(gen_motherindex->at(itx) <= 0)continue;
		if(gen_motherindex->at(itx) >= int(gen_pt->size()))continue;
		if(gen_pt->at(itx) < 0.1) continue;
		if(fabs(gen_pdgid->at(itx)) ==11) {
		  if(gen_motherindex->at(itx) == gen_motherindex->at(matched_index)) { charge_sum+= gen_pdgid->at(itx); n_el_from_el++;
		  }
		}
	      }// end of truth loop to check Conv
	      if((gen_pdgid->at(matched_index)  * pdgid) < 0 )  el.SetIsPhotonConversion(true);
	      else  el.SetIsPhotonConversion(false);
	      if(n_el_from_el ==3&& (fabs(charge_sum) == 11))  el.SetIsPhotonConversion(true);

	      if(fabs(gen_pdgid->at(gen_motherindex->at(mother_index))) > 50) {isprompt=false; eltype=12;}
	      
	      if(pdgid * electrons_q->at(iel) > 0 )     {el.SetIsChargeFlip(true); eltype=13; }
	      else     el.SetIsChargeFlip(false);
	    }
	    
	   
	    
	  }/// end of status 2 check
	  else {
	    /// using new method for matching: These events are set as prompt 
	    isprompt=true;mother_pdgid=-99999; mother_index=mindex; from_tau=false; 
	    eltype=14;
	    int n_el_from_eg=0;  
	    vector<KTruth> vel_tmp;
	    bool isthirdel_fromconv(false);
	    bool neutrino_invertex(false);
	    for (UInt_t itx=0; itx< gen_pt->size(); itx++ ){
	      if(gen_motherindex->at(itx) <= 0)continue;
	      if(gen_motherindex->at(itx) >= int(gen_pt->size()))continue;
	      if(gen_pt->at(itx) < 0.1) continue;
	      if(fabs(gen_pdgid->at(itx)) ==11) {
		if(gen_motherindex->at(itx) == gen_motherindex->at(matched_index)) {  n_el_from_eg++;
		  if(n_el_from_eg==3){isthirdel_fromconv=true; }
		  if(gen_status->at(itx) ==1){
		    KTruth truthe;
		    truthe.SetPtEtaPhiE(gen_pt->at(itx), gen_eta->at(itx), gen_phi->at(itx), gen_energy->at(itx));
		    vel_tmp.push_back(truthe);
		  }
		}
	      }
	      
	      if(fabs(gen_pdgid->at(itx)) ==12) {

		int index_mother_nu=gen_motherindex->at(itx);
		while (fabs(index_mother_nu) ==12){
		  index_mother_nu=gen_motherindex->at(index_mother_nu);
		}
		
		if(index_mother_nu == mindex) {
		  neutrino_invertex=true;
		}
	      }
	    } // end of truth loop to check Conv
	    
	    if(neutrino_invertex) eltype=15;

	    if(vel_tmp.size() ==2) {
	      KParticle ll = vel_tmp[0] + vel_tmp[1];
	      if(fabs(ll.M()) < 5.) eltype=16;
	    }


	    if((gen_pdgid->at(matched_index)  * pdgid) < 0 )  {el.SetIsPhotonConversion(true);  eltype=17;}
	    else el.SetIsPhotonConversion(false);
	    
	    if(n_el_from_eg ==3&&!isthirdel_fromconv)  {el.SetIsPhotonConversion(true); eltype=18;}
	    if(isthirdel_fromconv&&n_el_from_eg ==3){
	      if(pdgid * electrons_q->at(iel) > 0 )   {
		el.SetIsPhotonConversion(true);
		eltype=18;
	      }
	    }
	    if(pdgid * electrons_q->at(iel) > 0 )  {
	      el.SetIsChargeFlip(true);
	      if(eltype==17  || eltype == 18){
		if(eltype==17 ) eltype=19;
		if(eltype==18 ) eltype=20;
	      }
	      else eltype=21;
	    }
	    else     el.SetIsChargeFlip(false);
	    
	    
	  }  
	}
      }      /// In case no status 1 electron is found : classify electron fake
      else{
	if(gen_pt){
	  for (UInt_t it=0; it< gen_pt->size(); it++ ){
	    if(gen_motherindex->at(it) <= 0)continue;
	    if(gen_motherindex->at(it) >= int(gen_pt->size()))continue;
	    if(gen_pt->at(it) < 0.1) continue;
	    
	    double match_eta =electrons_eta->at(iel);
	    double match_phi =electrons_phi->at(iel);
	    double dr = sqrt( pow(fabs( match_eta - gen_eta->at(it)),2.0) +  pow( fabs(TVector2::Phi_mpi_pi( match_phi - gen_phi->at(it))),2.0));
	  
	    bool already_matched=false;
	    for(unsigned int im=0; im < matched_truth.size();im++){
	      if(it == unsigned(matched_truth.at(im))) already_matched=true;
	    }
	    
	    // in coversion case  2 matched electrons to one photon
	    if(fabs(gen_pdgid->at(it)) != 22 && already_matched) continue;
	    
	    if (dr <0.1){
	      matched_in_Dr=true;
	      int mindex= gen_motherindex->at(it);
	      float pdgid = gen_pdgid->at(it);

	      
	      /// Unlikely to have mother as electron but just in case
	      while ( (fabs(gen_pdgid->at((mindex))) == 11)) {
		mindex=gen_motherindex->at(mindex);
	      }
	      // isprompt = false since it failed status 1 matching
	      isprompt=false;
	      /// mother index of first non electron
	      mother_pdgid=gen_pdgid->at(mindex);
	      mother_index=mindex;
	      matched_index = it;
	      mc_pdgid= int(gen_pdgid->at(it));
	      if(fabs(pdgid) == 22) {

		if(fabs(mother_pdgid) > 50) eltype=22;
		else eltype=23;		
		
		//// This case is not a conversion
		//el.SetIsPhotonConversion(true);
		
		//if(gen_pdgid->at(gen_motherindex->at(it)) * electrons_q->at(iel) > 0 )     el.SetIsChargeFlip(true);
		//else     el.SetIsChargeFlip(false);
		
		from_tau=false;
		break;
	      }
	      
	      else if(fabs(pdgid) == 15){
		from_tau=true;
		if(fabs(mother_pdgid) > 50) eltype=24;
		else eltype=25;
	      }
	      else if(fabs(pdgid) == 1){
		eltype=26;
	      }
	      else if(fabs(pdgid) == 2){
                eltype=27;
              }
	      else if(fabs(pdgid) == 3){
                eltype=28;
              }
	      else if(fabs(pdgid) == 4){
                eltype=29;
              }
	      else if(fabs(pdgid) == 5){
                eltype=30;
              }
	      else if(fabs(pdgid) == 21){
                eltype=31;
              }
	      else eltype=32;
	      
	    }// dr req
	  }// loop over gen vector
	}// require gen info
      }// no status 1 match
      }/// END OF TRUTH MATCHING
      
      matched_truth.push_back(matched_index);
      ///matched_index is index which matches reco muon with smallest dR
    ///- If multiple status 1 muons are matched look at closest in pt
    ///- In no status 1 is matched set as not prompt butlook for closest particle in dR
    /// - In noparticles within dR < 0.1 matched_in_Dr= false
      
      el.SetType(eltype);
    if(!matched_in_Dr){
      el.SetIsMCMatched(false);
      el.SetIsFromTau(false);
      el.SetMCMatchedPdgId(-1);
      el.SetMotherPdgId(0);
      el.SetMotherTruthIndex(-1);
      el.SetMCTruthIndex(-1);
    }
    else{
      
      if(!isprompt){
	if((gen_isprompt->at(matched_index) ==1 )&& (gen_status->at(matched_index) == 1)){
	  
	  //cout << "gen_istaudecayproduct =  " << gen_istaudecayproduct->at(matched_index)  << endl;
	  //cout << "gen_isprompttaudecayproduct =  " <<  gen_isprompttaudecayproduct->at(matched_index)  << endl;
	  if(!(gen_istaudecayproduct->at(matched_index)   || gen_isprompttaudecayproduct->at(matched_index))){
	    //cout << "matched as prompt yet status flag is not prompt" << endl;
	    //cout << "matched_index = " << matched_index << endl;
	    //cout << "reco "<< electrons_pt->at(iel)<< " " << electrons_eta->at(iel)  << " " << electrons_phi->at(iel) << endl;;
	    //for (UInt_t it=0; it< gen_pt->size(); it++ ){
	    //	      if(gen_motherindex->at(it) <= 0)continue;
	    //if(gen_motherindex->at(it) >= int(gen_pt->size()))continue;
	    //if(gen_pt->at(it) < 0.1) continue;
	    //cout << it << " " << gen_pt->at(it)  << " " << gen_eta->at(it) << " " << gen_phi->at(it)<< " " << gen_pdgid->at(it) << "  " << gen_status->at(it) << " " << gen_pdgid->at(gen_motherindex->at(it)) <<" "  <<  gen_motherindex->at(it) << " " << gen_isprompt->at(it)  <<endl;
	      
	    //}
	  }
	}
      }
      
      //if(gen_isprompt->at(matched_index) 
      el.SetIsMCMatched(isprompt);
      el.SetIsFromTau(from_tau);
      el.SetMotherPdgId(mother_pdgid);
      el.SetMCMatchedPdgId(mc_pdgid);
      el.SetMotherTruthIndex(mother_index);
      el.SetMCTruthIndex(matched_index);
    }
    }
    electrons.push_back(el);
  }
  m_logger << DEBUG << "END electrons " << LQLogger::endmsg;
  std::sort( electrons.begin(), electrons.end(), isHigherPt );
  
  return electrons;
}


void SKTreeFiller::ERRORMessage(TString comment){
  
  m_logger << ERROR << "SKTreeFiller had a probleming filling " << comment << ". This variable is not present in the current LQntuples." << LQLogger::endmsg;   
}



std::vector<KGenJet> SKTreeFiller::GetAllGenJets(){

  std::vector<KGenJet> genjets;
  if(isData) return genjets;
  if(!LQinput){
    if(k_inputgenjets){
      for(std::vector<KGenJet>::iterator kit  = k_inputgenjets->begin(); kit != k_inputgenjets->end(); kit++){
	genjets.push_back(*kit);
      }
    }
    return genjets;
  }
  if(k_cat_version < 3){
    for (UInt_t ijet=0; ijet< slimmedGenJets_pt->size(); ijet++) {
      KGenJet jet;
      jet.SetPtEtaPhiE(slimmedGenJets_pt->at(ijet), slimmedGenJets_eta->at(ijet), slimmedGenJets_phi->at(ijet), slimmedGenJets_energy->at(ijet));
      genjets.push_back(jet);
    }
    return genjets;
  }
  
  for (UInt_t ijet=0; ijet< genjet_pt->size(); ijet++) {
    KGenJet jet;
    jet.SetPtEtaPhiE(genjet_pt->at(ijet), genjet_eta->at(ijet), genjet_phi->at(ijet), genjet_energy->at(ijet));
    jet.SetGenJetEMF(genjet_emf->at(ijet));
    jet.SetGenJetHADF(genjet_hadf->at(ijet));
    jet.SetGenJetPDGID(int(genjet_hadf->at(ijet)));
    
    genjets.push_back(jet);
  }
  return genjets;
}


std::vector<KJet> SKTreeFiller::GetAllJets(){

  std::vector<KJet> jets;
  if(!LQinput){

    for(std::vector<KJet>::iterator kit  = k_inputjets->begin(); kit != k_inputjets->end(); kit++){
      jets.push_back(*kit);
    }
    return jets;
  }

  for (UInt_t ijet=0; ijet< jets_eta->size(); ijet++) {
    KJet jet;
    if(jets_pt->at(ijet) != jets_pt->at(ijet)) continue;

    if(isData){
      jet.SetPtEtaPhiE(jets_pt->at(ijet), jets_eta->at(ijet), jets_phi->at(ijet), jets_energy->at(ijet));
    }
    else{
      jet.SetPtEtaPhiE(jets_pt->at(ijet), jets_eta->at(ijet), jets_phi->at(ijet), jets_energy->at(ijet));

      // https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution      
      /// Measurements show that the jet energy resolution (JER) in data is worse than in the simulation and the jets in MC need to be smeared to describe the data.

      jet*= jets_smearedRes->at(ijet);
      jet.SetIsMCSmeared(true);
    }
    jet.SetJetPassLooseID(jets_isLoose->at(ijet));
    jet.SetJetPassTightID(jets_isTight->at(ijet));
    jet.SetJetPassTightLepVetoID(jets_isTightLepVetoJetID->at(ijet));
    
    jet.SetJetPileupIDMVA(jets_PileupJetId->at(ijet));

    if(jets_PileupJetId){ 
      if(std::abs(jets_eta->at(ijet)) < 2.6){
	if(jets_PileupJetId->at(ijet) > 0.3) jet.SetJetPileupIDLooseWP(true);
	else jet.SetJetPileupIDLooseWP(false);
	if(jets_PileupJetId->at(ijet) > 0.7) jet.SetJetPileupIDMediumWP(true);
	else jet.SetJetPileupIDMediumWP(false);
	if(jets_PileupJetId->at(ijet) > 0.9)jet.SetJetPileupIDTightWP(true);
	else jet.SetJetPileupIDTightWP(false);
      }
      else{
	if(jets_PileupJetId->at(ijet) > -0.55) jet.SetJetPileupIDLooseWP(true);
        else jet.SetJetPileupIDLooseWP(false);
        if(jets_PileupJetId->at(ijet) > -0.3) jet.SetJetPileupIDMediumWP(true);
	else jet.SetJetPileupIDMediumWP(false);
        if(jets_PileupJetId->at(ijet) > -0.1)jet.SetJetPileupIDTightWP(true);
	else jet.SetJetPileupIDTightWP(false);
      }
    }
    
    
    /// BTAG variables
    if(jets_CSVInclV2) jet.SetBTagInfo(snu::KJet::CSVv2, jets_CSVInclV2->at(ijet));
    if(jets_CMVAV2)    jet.SetBTagInfo(snu::KJet::cMVAv2, jets_CMVAV2->at(ijet));
    if(jets_JetProbBJet)  jet.SetBTagInfo(snu::KJet::JETPROB, jets_JetProbBJet->at(ijet)); 
    
    //if(jets_iCSVCvsL) {
    //      if(jets_iCSVCvsL->size() > 0)jet.SetCTagInfo(snu::KJet::iCSVCvsL, jets_iCSVCvsL->at(ijet));
    //    }
    if(jets_CCvsLT){
      if(jets_CCvsLT->size() > 0) jet.SetCTagInfo(snu::KJet::CCvsLT, jets_CCvsLT->at(ijet));
    }
    if(jets_CCvsBT){
      if(jets_CCvsBT->size() > 0)jet.SetCTagInfo(snu::KJet::CCvsBT, jets_CCvsBT->at(ijet));
    }
    jet.SetVtxMass(jets_vtxMass->at(ijet));
    jet.SetVtx3DVal(jets_vtx3DVal->at(ijet));
    jet.SetVtx3DSig(jets_vtx3DSig->at(ijet));
    jet.SetVtxNTracks(jets_vtxNtracks->at(ijet));
    
    // flavour
    jet.SetJetPartonFlavour(jets_partonFlavour->at(ijet));
    jet.SetJetHadronFlavour(jets_hadronFlavour->at(ijet));    
    jet.SetJetPartonPdgId(jets_partonPdgId->at(ijet));
    
    jet.SetJetChargedEmEF(jets_chargedEmEnergyFraction->at(ijet));
    /// JEC and uncertainties
    jet.SetJetScaledDownEnergy(jets_shiftedEnDown->at(ijet));
    jet.SetJetScaledUpEnergy(jets_shiftedEnUp->at(ijet));
    jet.SetSmearedResDown(jets_smearedResDown->at(ijet));
    jet.SetSmearedResUp(jets_smearedResUp->at(ijet));
    jet.SetSmearedRes(jets_smearedRes->at(ijet));
    
    jets.push_back(jet);
  }// end of jet 
  
  
  std::sort( jets.begin(), jets.end(), isHigherPt );
  
  m_logger << DEBUG << "PFJet size = " << jets.size() << LQLogger::endmsg;
  return jets;
}



std::vector<KFatJet> SKTreeFiller::GetAllFatJets(){

  std::vector<KFatJet> fatjets;

  if(k_cat_version <  7) return fatjets;

  if(!LQinput){

    for(std::vector<KFatJet>::iterator kit  = k_inputfatjets->begin(); kit != k_inputfatjets->end(); kit++){
      fatjets.push_back(*kit);
    }
    return fatjets;
  }

  for (UInt_t ijet=0; ijet< fatjets_eta->size(); ijet++) {
    KFatJet jet;
    if(fatjets_pt->at(ijet) != fatjets_pt->at(ijet)) continue;
    jet.SetPtEtaPhiE(fatjets_pt->at(ijet), fatjets_eta->at(ijet), fatjets_phi->at(ijet), fatjets_energy->at(ijet));

    jet.SetJetPassLooseID(fatjets_isLoose->at(ijet));
    jet.SetJetPassTightID(fatjets_isTight->at(ijet));
    jet.SetJetPassTightLepVetoID(fatjets_isTightLepVetoJetID->at(ijet));

    jet.SetJetPileupIDMVA(fatjets_PileupJetId->at(ijet));

    if(fatjets_PileupJetId){
      if(std::abs(fatjets_eta->at(ijet)) < 2.6){
        if(fatjets_PileupJetId->at(ijet) > 0.3) jet.SetJetPileupIDLooseWP(true);
        else jet.SetJetPileupIDLooseWP(false);
        if(fatjets_PileupJetId->at(ijet) > 0.7) jet.SetJetPileupIDMediumWP(true);
        else jet.SetJetPileupIDMediumWP(false);
        if(fatjets_PileupJetId->at(ijet) > 0.9)jet.SetJetPileupIDTightWP(true);
        else jet.SetJetPileupIDTightWP(false);
      }
      else{
        if(fatjets_PileupJetId->at(ijet) > -0.55) jet.SetJetPileupIDLooseWP(true);
        else jet.SetJetPileupIDLooseWP(false);
        if(fatjets_PileupJetId->at(ijet) > -0.3) jet.SetJetPileupIDMediumWP(true);
        else jet.SetJetPileupIDMediumWP(false);
        if(fatjets_PileupJetId->at(ijet) > -0.1)jet.SetJetPileupIDTightWP(true);
        else jet.SetJetPileupIDTightWP(false);
      }
    }


    /// BTAG variables                                                                                                                                                                                                                                                                                          
    if(fatjets_CSVInclV2) jet.SetBTagInfo(snu::KFatJet::CSVv2, fatjets_CSVInclV2->at(ijet));
    if(fatjets_CMVAV2)    jet.SetBTagInfo(snu::KFatJet::cMVAv2, fatjets_CMVAV2->at(ijet));
    if(fatjets_JetProbBJet)  jet.SetBTagInfo(snu::KFatJet::JETPROB, fatjets_JetProbBJet->at(ijet));

    if(fatjets_CCvsLT){
      if(fatjets_CCvsLT->size() > 0) jet.SetCTagInfo(snu::KFatJet::CCvsLT, fatjets_CCvsLT->at(ijet));
    }
    if(fatjets_CCvsBT){
      if(fatjets_CCvsBT->size() > 0)jet.SetCTagInfo(snu::KFatJet::CCvsBT, fatjets_CCvsBT->at(ijet));
    }
    jet.SetVtxMass(fatjets_vtxMass->at(ijet));
    jet.SetVtx3DVal(fatjets_vtx3DVal->at(ijet));
    jet.SetVtx3DSig(fatjets_vtx3DSig->at(ijet));
    jet.SetVtxNTracks(fatjets_vtxNtracks->at(ijet));

    // flavour                                                                                                                                                                                                                                                                                                  
    jet.SetJetPartonFlavour(fatjets_partonFlavour->at(ijet));
    jet.SetJetHadronFlavour(fatjets_hadronFlavour->at(ijet));
    jet.SetJetPartonPdgId(fatjets_partonPdgId->at(ijet));

    jet.SetJetChargedEmEF(fatjets_chargedEmEnergyFraction->at(ijet));

    jet.SetJetScaledDownEnergy(fatjets_shiftedEnDown->at(ijet));
    jet.SetJetScaledUpEnergy(fatjets_shiftedEnUp->at(ijet));
    jet.SetSmearedResDown(fatjets_smearedResDown->at(ijet));
    jet.SetSmearedResUp(fatjets_smearedResUp->at(ijet));
    jet.SetSmearedRes(fatjets_smearedRes->at(ijet));


    jet.SetTau1(fatjets_tau1->at(ijet));
    jet.SetTau2(fatjets_tau2->at(ijet));
    jet.SetTau3(fatjets_tau3->at(ijet));

    jet.SetPrunedMass(fatjets_prunedmass->at(ijet));
    jet.SetSoftDropMass(fatjets_softdropmass->at(ijet));

    jet.SetPuppiTau1(fatjets_puppi_tau1->at(ijet));
    jet.SetPuppiTau2(fatjets_puppi_tau2->at(ijet));
    jet.SetPuppiTau3(fatjets_puppi_tau3->at(ijet));
    jet.SetPuppiPt(fatjets_puppi_pt->at(ijet));
    jet.SetPuppiEta(fatjets_puppi_eta->at(ijet));
    jet.SetPuppiPhi(fatjets_puppi_phi->at(ijet));
    jet.SetPuppiM(fatjets_puppi_m->at(ijet));
    
      
    fatjets.push_back(jet);
  }// end of jet                                                   
  std::sort( fatjets.begin(), fatjets.end(), isHigherPt );

  m_logger << DEBUG << "PFJet size = " << fatjets.size() << LQLogger::endmsg;
  return fatjets;
}



std::vector<KMuon> SKTreeFiller::GetAllMuons(){

  std::vector<KMuon> muons ;
  
  if(!LQinput){
    for(std::vector<KMuon>::iterator kit  = k_inputmuons->begin(); kit != k_inputmuons->end(); kit++){
      muons.push_back(*kit);
    }  
    return muons;
  }

  m_logger << DEBUG << "Filling Muons" << LQLogger::endmsg;

  vector<int> matched_truth;
  for (UInt_t ilep=0; ilep< muon_eta->size(); ilep++) {
    KMuon muon;
    if(muon_pt->at(ilep) != muon_pt->at(ilep)) continue;
    m_logger << DEBUG << "Filling global pt/eta ... " << LQLogger::endmsg;
   
    muon.SetTrigMatch(muon_trigmatch->at(ilep));
      
    /// GENERAL
    
    muon.SetISPF(muon_isPF->at(ilep));
    muon.SetIsGlobal(muon_isGlobal->at(ilep));

    muon.SetIsTracker(muon_isTracker->at(ilep));
    muon.SetIsLoose(muon_isLoose->at(ilep));
    muon.SetIsMedium(muon_isMedium->at(ilep));
    muon.SetIsTight(muon_isTight->at(ilep));
    muon.SetIsSoft(muon_isSoft->at(ilep));

    if(muon_shiftedEup){
      muon.SetShiftedEUp(muon_shiftedEup->at(ilep));
      muon.SetShiftedEDown(muon_shiftedEdown->at(ilep));
    }
    
    
    muon.SetPtEtaPhiE(muon_pt->at(ilep), muon_eta->at(ilep),muon_phi->at(ilep), muon_energy->at(ilep));
    if(k_cat_version > 4){
      muon.SetRochPt(muon_pt->at(ilep));
      muon.SetRochEta(muon_roch_eta->at(ilep));
      muon.SetRochPhi(muon_roch_phi->at(ilep));
      muon.SetRochE(muon_roch_energy->at(ilep));
      muon.SetRochM(muon_roch_m->at(ilep));
    }
    else{
      muon.SetRochPt(muon_pt->at(ilep));
      muon.SetRochEta(muon_eta->at(ilep));
      muon.SetRochPhi(muon_phi->at(ilep));
      muon.SetRochE(muon_energy->at(ilep));
      muon.SetRochM(muon_m->at(ilep));
    }
    muon.SetCharge(muon_q->at(ilep));
     
    m_logger << DEBUG << "Filling ms pt/eta ... " << LQLogger::endmsg;
 
    muon.SetRelIso(0.3,muon_relIso03->at(ilep));
    muon.SetRelIso(0.4,muon_relIso04->at(ilep));
    if(muon_minirelIso)muon.SetRelMiniIso(muon_minirelIso->at(ilep));

    if(k_cat_version  > 7){
      muon.SetMiniAODPt(muon_pt->at(ilep));
      muon.SetMiniAODRelIso(0.3,muon_relIso03->at(ilep));
      muon.SetMiniAODRelIso(0.4,muon_relIso04->at(ilep));
      muon.SetIsRochesterCorrected(false);
    }
    muon.Setdz(muon_dz->at(ilep));
    muon.Setdxy(muon_dxy->at(ilep));
    if(muon_sigdxy)muon.Setdxy_sig(muon_sigdxy->at(ilep));
    //// chi2
    muon.SetGlobalchi2( muon_normchi->at(ilep));
        
    /// hits
    muon.SetValidHits( muon_validhits->at(ilep));
    muon.SetPixelValidHits( muon_validpixhits->at(ilep));
    muon.SetValidStations( muon_matchedstations->at(ilep));
    muon.SetLayersWithMeasurement ( muon_trackerlayers->at(ilep));
    
    muon.SetMCMatched(muon_matched->at(ilep));


    muon.SetTrackVx(muon_x->at(ilep));
    muon.SetTrackVy(muon_y->at(ilep));
    muon.SetTrackVz(muon_z->at(ilep));

    //// Set Is ChargeFlip
    bool isprompt= false;
    bool from_tau = false;

    int mother_index=-1;
    int mother_pdgid=-1;
    int matched_index=-1;
    int mc_pdgid=-1;
    bool matched_in_Dr=false;
    
    if(k_cat_version > 3){

    if(gen_pt){
      float min_Dr=0.1;

      for (UInt_t it=0; it< gen_pt->size(); it++ ){
        if(gen_motherindex->at(it) <= 0)continue;
        if(gen_motherindex->at(it) >= int(gen_pt->size()))continue;
        if(gen_pt->at(it) < 0.1) continue;
	
	if(gen_pt->at(it) < 5.) continue;

	double match_eta =muon_eta->at(ilep);
	double match_phi =muon_phi->at(ilep);
	double dr = sqrt( pow(fabs( match_eta - gen_eta->at(it)),2.0) +  pow( fabs(TVector2::Phi_mpi_pi( match_phi - gen_phi->at(it))),2.0));
	/// Matching using instructions on
        /// https://indico.cern.ch/event/292928/contributions/1650088/attachments/547844/755123/talk_electron_contribution.pdf
        ///

	bool already_matched=false;
	for(unsigned int im=0; im < matched_truth.size();im++){
	  if(it == unsigned(matched_truth.at(im))) already_matched=true;
	}
	if(already_matched) continue;

        /// Match requires to status 1 muon
        if(gen_status->at(it) != 1) continue;
        if(fabs(gen_pdgid->at(it)) != 13) continue;
	
	if(matched_in_Dr){
	  if(gen_motherindex->at(it) == gen_motherindex->at(matched_index)){
	    if (dr < 0.1){
	      // in case 2+ electrons from same mother electron (conversion) also match in pt
	      if( fabs(gen_pt->at(it)-muon_pt->at(ilep)) < fabs(gen_pt->at(matched_index)-muon_pt->at(ilep))) matched_index=it;
	    }
	  }
	  else   if (dr < min_Dr){
	    /// find closest match in dR to status 1
	    matched_in_Dr=true;
	    min_Dr= dr;
	    /// set index of matched status 1 muon
	    matched_index=it;
	  }
	}
	else if (dr < min_Dr){
          /// find closest match in dR to status 1
          matched_in_Dr=true;
          min_Dr= dr;
	  
	  /// set index of matched status 1 muon
          matched_index=it;
        }
      }// end of gen loop to find status 1 muon

      ///// treat case where there is a matched status 1 muon:
      //// classify into prompt:Fake:FromTau

      
      if(matched_in_Dr){
        /// Find closest non muon ancesteror
        float pdgid = gen_pdgid->at(matched_index);

        // mc_pdgid = closest matched status 1 pdgid
        mc_pdgid= int(pdgid);

        // mindex = mother index: will loop to find first non el mother
        int mindex= matched_index;

        while ( (fabs(gen_pdgid->at(mindex)) == 13)) {
          pdgid = gen_pdgid->at(mindex);
          mindex=gen_motherindex->at(mindex);
        }
        /// pdgid is now of muon from non muon mother
        //  mindex = index for mother of non muon ancestor

        if( (fabs(gen_pdgid->at(mindex)) == 23) || (fabs(gen_pdgid->at(mindex)) == 24)) {
	  /// Check if it is a chargeflip.
          if(pdgid * muon_q->at(ilep) > 0 )     muon.SetIsChargeFlip(true);
          else     muon.SetIsChargeFlip(false);

	  int n_mu_from_mother=0;
	  for (UInt_t itx=0; itx< gen_pt->size(); itx++ ){
	    if(gen_motherindex->at(itx) <= 0)continue;
	    if(gen_motherindex->at(itx) >= int(gen_pt->size()))continue;
	    if(gen_pt->at(itx) < 0.1) continue;
	    if(fabs(gen_pdgid->at(itx)) ==13) {
	      if(gen_motherindex->at(itx) == gen_motherindex->at(matched_index)) n_mu_from_mother++;
	    }
	  }
	  if(n_mu_from_mother == 3)  muon.SetIsPhotonConversion(true);
          mother_index=mindex;
          mother_pdgid=gen_pdgid->at(mindex);
          isprompt=true; /// means is prompt
        }/// end of Z/W
        else {
          if(gen_status->at(mindex) == 2){
            if(fabs(gen_pdgid->at(mindex)) > 50) {isprompt=false; mother_pdgid=gen_pdgid->at(mindex); mother_index=mindex; from_tau=false;}
	    else isprompt=true;
            if(fabs(gen_pdgid->at(mindex)) == 15){
              isprompt=true; mother_pdgid=gen_pdgid->at(mindex);  mother_index=mindex; from_tau=true;
              // Check if el from tau  is CF

	      if(pdgid * muon_q->at(ilep) > 0 )     muon.SetIsChargeFlip(true);
              else     muon.SetIsChargeFlip(false);

	      int n_mu_from_mother=0;
	      for (UInt_t itx=0; itx< gen_pt->size(); itx++ ){
		if(gen_motherindex->at(itx) <= 0)continue;
		if(gen_motherindex->at(itx) >= int(gen_pt->size()))continue;
		if(gen_pt->at(itx) < 0.1) continue;
		if(fabs(gen_pdgid->at(itx)) ==13) {
		  if(gen_motherindex->at(itx) == gen_motherindex->at(matched_index)) n_mu_from_mother++;
		}
	      }
	      if(n_mu_from_mother == 3)  muon.SetIsPhotonConversion(true);

            }
          }/// end of status 2 check
          else {
            /// using new method for matching: These events are set as prompt
            isprompt=true;mother_pdgid=-99999; mother_index=mindex; from_tau=false;

            if(pdgid * muon_q->at(ilep) > 0 )    muon.SetIsChargeFlip(true);
            else     muon.SetIsChargeFlip(false);
	    int n_mu_from_mother=0;
	    for (UInt_t itx=0; itx< gen_pt->size(); itx++ ){
	      if(gen_motherindex->at(itx) <= 0)continue;
	      if(gen_motherindex->at(itx) >= int(gen_pt->size()))continue;
	      if(gen_pt->at(itx) < 0.1) continue;
	      if(fabs(gen_pdgid->at(itx)) ==13) {
		if(gen_motherindex->at(itx) == gen_motherindex->at(matched_index)) n_mu_from_mother++;
	      }
	    }
	    if(n_mu_from_mother == 3)  muon.SetIsPhotonConversion(true);
          }
        }
      }      /// In case no status 1 muon is found : classify muon fake
      else{
	if(muon_q->size() ==2){
	  for (UInt_t itxx=0; itxx< gen_pt->size(); itxx++ ){
            if(gen_motherindex->at(itxx) <= 0)continue;
            if(gen_motherindex->at(itxx) >= int(gen_pt->size()))continue;
            if(gen_pt->at(itxx) < 0.1) continue;
	    //cout << itxx << " " << gen_pdgid->at(itxx) << " " << gen_pdgid->at(gen_motherindex->at(itxx)) << " " << gen_eta->at(itxx) << " " << gen_phi->at(itxx) << endl;
	  }

	}
        if(gen_pt){
	  for (UInt_t it=0; it< gen_pt->size(); it++ ){
            if(gen_motherindex->at(it) <= 0)continue;
            if(gen_motherindex->at(it) >= int(gen_pt->size()))continue;
            if(gen_pt->at(it) < 0.1) continue;

	    bool already_matched=false;
	    for(unsigned int im=0; im < matched_truth.size();im++){
	      if(it == unsigned(matched_truth.at(im))) already_matched=true;
	    }
	    if(already_matched) continue;

            double match_eta =muon_eta->at(ilep);
            double match_phi =muon_phi->at(ilep);
            double dr = sqrt( pow(fabs( match_eta - gen_eta->at(it)),2.0) +  pow( fabs(TVector2::Phi_mpi_pi( match_phi - gen_phi->at(it))),2.0));

            if (dr <0.1){
              matched_in_Dr=true;
              int mindex= gen_motherindex->at(it);
              float pdgid = gen_pdgid->at(it);
	      
	      /// Unlikely to have mother as muon but just in case
              while ( (fabs(gen_pdgid->at(mindex)) == 13)) {
                mindex=gen_motherindex->at(mindex);
              }
              // isprompt = false since it failed status 1 matching
              isprompt=false;
              /// mother index of first non muon
              mother_pdgid=gen_pdgid->at(mindex);
              mother_index=mindex;
              matched_index = it;
	      mc_pdgid= int(pdgid);
              if(fabs(pdgid) == 22) {
                from_tau=false;
                break;
              }
              if(fabs(pdgid) == 15)from_tau=true;
	      
	      if(fabs(gen_pdgid->at(gen_motherindex->at(it))) == 13){
		int n_mu_from_mother=0;
		for (UInt_t itx=0; itx< gen_pt->size(); itx++ ){
		  if(gen_motherindex->at(itx) <= 0)continue;
		  if(gen_motherindex->at(itx) >= int(gen_pt->size()))continue;
		  if(gen_pt->at(itx) < 0.1) continue;
		  if(fabs(gen_pdgid->at(itx)) ==13) {
		    if(gen_motherindex->at(itx) == gen_motherindex->at(matched_index)) n_mu_from_mother++;
		  }
		}
		if(n_mu_from_mother == 3)  muon.SetIsPhotonConversion(true);
	      }else{
		muon.SetIsPhotonConversion(false);
	      }
	    }// dr req
          }// loop over gen vector
        }// require gen info
      }// no status 1 match
    }
    
    /// matched_index is index which matches reco muon with smallest dR
    /// - If multiple status 1 muons are matched look at closest in pt
    /// - In no status 1 is matched set as not prompt but look for closest particle in dR 
    /// - In no particles within dR < 0.1 matched_in_Dr = false
    matched_truth.push_back(matched_index);

    if(!matched_in_Dr){
      muon.SetMCMatched(false);
      muon.SetIsFromTau(false);
      muon.SetMCMatchedPdgId(-1);
      muon.SetMotherPdgId(0);
      muon.SetMotherTruthIndex(-1);
      muon.SetMCTruthIndex(-1);
    }
    else{

      muon.SetMCMatched(isprompt);
      muon.SetIsFromTau(from_tau);
      muon.SetMotherPdgId(mother_pdgid);
      muon.SetMCMatchedPdgId(mc_pdgid);
      muon.SetMotherTruthIndex(mother_index);
      muon.SetMCTruthIndex(matched_index);
    }
    }

    /// Fill vector
    muons.push_back(muon);
  }
  
  std::sort( muons.begin(), muons.end(), isHigherPt );
  m_logger << DEBUG << "End of Muon Filling" << LQLogger::endmsg;
  return muons;
}

  

std::vector<snu::KTruth>   SKTreeFiller::GetTruthParticles(int np){
  
  m_logger << DEBUG << "Filling Truth" << LQLogger::endmsg;
  std::vector<snu::KTruth> vtruth;

  if(isData) return vtruth;

  int counter=0;

  if(!LQinput){

    for(std::vector<KTruth>::iterator kit  = k_inputtruth->begin(); kit != k_inputtruth->end(); kit++, counter++){
      if(counter == np)  break;
      vtruth.push_back(*kit);
    }
    return vtruth;
  }
  
  m_logger << DEBUG << "Filling truth Info: " << gen_pt->size() << LQLogger::endmsg;

  for (UInt_t it=0; it< gen_pt->size(); it++ , counter++) {
    
    if(counter == np)  break;
    KTruth truthp;
    truthp.SetPtEtaPhiE(double(gen_pt->at(it)), double(gen_eta->at(it)), double(gen_phi->at(it)), double(gen_energy->at(it)));
    truthp.SetParticlePdgId(gen_pdgid->at(it));
    truthp.SetParticleStatus(gen_status->at(it));
    truthp.SetParticleIndexMother(gen_motherindex->at(it));
    
    if(k_cat_version > 3){
      // To save space set a single int as the flag. 
      // 
      int truth_flag = 0;
      if(gen_isprompt->at(it)) truth_flag+=1;
      if(gen_isdecayedleptonhadron->at(it)) truth_flag+=10;
      if(gen_istaudecayproduct->at(it)) truth_flag+=100;
      if(gen_isprompttaudecayproduct->at(it)) truth_flag+=1000;
      if(gen_isdirecthadrondecayproduct->at(it)) truth_flag+=10000;
      if(gen_ishardprocess->at(it)) truth_flag+=100000;
      if(gen_fromhardprocess->at(it)) truth_flag+=1000000;
      if(gen_fromhardprocess_beforeFSR->at(it)) truth_flag+=10000000;
      truthp.SetStatusFlag(truth_flag);
    }
    
    vtruth.push_back(truthp);  
  }
  
  return vtruth;
}
