#include "FatJetSelection.h"
#include <iostream>


using namespace snu;

FatJetSelection::FatJetSelection(LQEvent ev) :BaseSelection() {
  k_lqevent = ev;
}

FatJetSelection::~FatJetSelection() {}

//// This code is used to make selection cuts to vectors of KFatJets


void FatJetSelection::BasicSelection(std::vector<KFatJet>& jetColl) {
  
  //// This is a basic set of cuts on jets

  std::vector<KFatJet> alljets = k_lqevent.GetFatJets();

  for (std::vector<KFatJet>::iterator jit = alljets.begin(); jit!=alljets.end(); jit++){
  
    if ( (jit->Pt() >= pt_cut_min) &&  (fabs(jit->Eta()) < eta_cut)){
      if ( PassUserID(PFJET_LOOSE, *jit) &&    (jit->Pt() >= pt_cut_min) &&  (fabs(jit->Eta()) < eta_cut))  jetColl.push_back(*jit);
    }
       
  }
}

void FatJetSelection::Selection(std::vector<KFatJet>& jetColl){
  
  std::vector<KFatJet> alljets = k_lqevent.GetFatJets();
  
  for (std::vector<KFatJet>::iterator jit = alljets.begin(); jit!=alljets.end(); jit++){

    bool pileupjet=false;
    if(applypileuptool) pileupjet =  ( !jit->PileupJetIDLoose());  ///---> CHECK THIS


    if(apply_ID) {
      if ( jit->Pt() >= pt_cut_min && jit->Pt() < pt_cut_max &&
           fabs(jit->Eta()) < eta_cut
           &&PassUserID(k_id, *jit) && !pileupjet)  jetColl.push_back(*jit);
    }
    else{
      if ( jit->Pt() >= pt_cut_min && jit->Pt() < pt_cut_max &&
           fabs(jit->Eta()) < eta_cut
           && PassUserID(PFJET_LOOSE, *jit)&& !pileupjet)  jetColl.push_back(*jit);
    }
  }

  BaseSelection::reset();
  return;

}  

void FatJetSelection::Selection(std::vector<KFatJet>& jetColl, bool LepVeto, std::vector<KMuon>& muonColl, std::vector<KElectron>& electronColl) {
  
  std::vector<KFatJet> alljets = k_lqevent.GetFatJets();
  
  std::vector<KFatJet> prejetColl; 
  
  for (std::vector<KFatJet>::iterator jit = alljets.begin(); jit!=alljets.end(); jit++){
    
    bool pileupjet=false;
    if(applypileuptool) pileupjet =  ( !jit->PileupJetIDLoose());  ///---> CHECK THIS


    if(apply_ID) {
      if ( jit->Pt() >= pt_cut_min && jit->Pt() < pt_cut_max &&
	   fabs(jit->Eta()) < eta_cut
	   &&PassUserID(k_id, *jit) && !pileupjet)  prejetColl.push_back(*jit);
    }
    else{
      if ( jit->Pt() >= pt_cut_min && jit->Pt() < pt_cut_max && 
	   fabs(jit->Eta()) < eta_cut
	   && PassUserID(PFJET_LOOSE, *jit)&& !pileupjet)  prejetColl.push_back(*jit);
    }
  } 

  for (UInt_t ijet = 0; ijet < prejetColl.size(); ijet++) {
    jetIsOK = true;
    for (UInt_t ilep = 0; ilep < muonColl.size(); ilep++) {
      if (muonColl[ilep].DeltaR( prejetColl[ijet] ) < 1.) {
	jetIsOK = false;  ilep = muonColl.size();
      }
    }/// End of muon loop
    for (UInt_t ilep = 0; ilep < electronColl.size(); ilep++) {
      if (electronColl[ilep].DeltaR( prejetColl[ijet] ) < 1. ) {
	jetIsOK = false;  ilep = electronColl.size();
      }
    }/// End of electron loop
    
    if(LepVeto){
      if (jetIsOK) jetColl.push_back( prejetColl[ijet] );
    }
    else{
      jetColl.push_back( prejetColl[ijet] );
    }
  }
  
  
  BaseSelection::reset();
  return;
  
}

void FatJetSelection::SelectFatJets( std::vector<KFatJet>& jetColl,  vector<pair<TString, TString> > vids, vector<pair<TString, float> > vidf,  float ptcut, float etacut ) {

  std::vector<KFatJet> alljets = k_lqevent.GetFatJets();
  
  int icut(0);
  if (ptcut == -999. || etacut == -999.){
    for(unsigned int iv=0; iv < vidf.size(); iv++){
      if(!Check(vidf[iv].second)) continue;
      if (vidf[iv].first =="ptmin") { icut++; if(ptcut == -999.)ptcut=vidf[iv].second;}
      if (vidf[iv].first =="|etamax|") {icut++;  if (etacut == -999.)etacut=vidf[iv].second;}
      if(icut ==2) break;
    }
  }

  for (std::vector<KFatJet>::iterator jit = alljets.begin(); jit!=alljets.end(); jit++){

    bool pass_selection=true;
    if (!PassUserID(*jit, vids)) pass_selection=false;
    if ( (jit->Pt() >= ptcut)  && fabs(jit->Eta()) < etacut && pass_selection )  jetColl.push_back(*jit);
  }

  
  
}

void FatJetSelection::SelectFatJets(std::vector<KFatJet>& jetColl, std::vector<KMuon> muonColl, std::vector<KElectron> electronColl, vector<pair<TString, TString> > vids, vector<pair<TString, float> > vidf,  float ptcut , float etacut ) {
  
  std::vector<KFatJet> pre_jetColl; 
  std::vector<KFatJet> alljets = k_lqevent.GetFatJets();


  int icut(0);
  float tau21cut(-999.);
  float masscut_min(-999.);
  float masscut_max(-999.);
  if (ptcut == -999. || etacut == -999.){
    for(unsigned int iv=0; iv < vidf.size(); iv++){
      if(!Check(vidf[iv].second)) continue;
      if (vidf[iv].first =="ptmin") { icut++; if(ptcut == -999.)ptcut=vidf[iv].second;}
      if (vidf[iv].first =="|etamax|") {icut++;  if (etacut == -999.)etacut=vidf[iv].second;}
      if (vidf[iv].first =="tau21") {icut++;  tau21cut=vidf[iv].second;}
      if (vidf[iv].first =="mass_min") {icut++;  masscut_min=vidf[iv].second;}
      if (vidf[iv].first =="mass_max") {icut++;  masscut_max=vidf[iv].second;}
      if(icut ==3) break;
    }
  }
  

  for (std::vector<KFatJet>::iterator jit = alljets.begin(); jit!=alljets.end(); jit++){
    
    bool pass_selection=true;
    if (!PassUserID(*jit, vids)) pass_selection=false;

    if(jit->Tau2()/jit->Tau1() > tau21cut) pass_selection=false;
    if(jit->PrunedMass() > masscut_min)  pass_selection=false;
    if(jit->PrunedMass() < masscut_max)  pass_selection=false;
    //cout << jit->Pt() << " " << jit->Eta() << " " << jit->PrunedMass() << " " << jit->Tau2()/jit->Tau1() << endl;
    if ( (jit->Pt() >= ptcut)  && fabs(jit->Eta()) < etacut && pass_selection )  pre_jetColl.push_back(*jit);

  }


  for (UInt_t ijet = 0; ijet < pre_jetColl.size(); ijet++) {
    jetIsOK = true;
    for (UInt_t ilep = 0; ilep < muonColl.size(); ilep++) {
      if (muonColl[ilep].DeltaR( pre_jetColl[ijet] ) < 1.) {
        jetIsOK = false;
	//cout << "Muon eta/phi = " << muonColl[ilep].Eta() << " " << muonColl[ilep].Phi() << endl;
        //cout << "FatJet eta/phi = " <<  pre_jetColl[ijet].Eta() << " " <<  pre_jetColl[ijet].Phi() << endl;

	ilep = muonColl.size();
      }
    }/// End of muon loop
    for (UInt_t ilep = 0; ilep < electronColl.size(); ilep++) {
      if (electronColl[ilep].DeltaR( pre_jetColl[ijet] ) < 1. ) {
        jetIsOK = false;
        ilep = electronColl.size();
      }
    }/// End of electron loop
    
    if (jetIsOK) jetColl.push_back( pre_jetColl[ijet] );
  }/// End of FatJet loop
  
}



bool FatJetSelection::PassUserID (ID id, snu::KFatJet jet){
  if      ( id == PFJET_LOOSE  ) return PassUserID_PFFatJetLoose  (jet);
  else if ( id == PFJET_TIGHT  ) return PassUserID_PFFatJetTight  (jet);
  else return false;
}


bool FatJetSelection::PassUserID (snu::KFatJet jet,vector<pair<TString, TString> > vids ){ 

  for(unsigned int idel =0; idel < vids.size(); idel++){
    if(vids[idel].second == "false") continue;

    if(vids[idel].first == "LooseID") {
      if(!jet.PassLooseID()) return false;
    }
    if(vids[idel].first == "TightID") {
      if(!jet.PassTightID())  return false;
    }
    if(vids[idel].first == "TightIDLepVeto"){
      if(!jet.PassTightLepVetoID()) return false;
    }
  }

  return true;
}



bool FatJetSelection::PassUserID_PFFatJetLoose ( snu::KFatJet jet){
  
  return jet.PassLooseID();
}


bool FatJetSelection::PassUserID_PFFatJetTight ( snu::KFatJet jet)
{
  return jet.PassTightID();
}


FatJetSelection& FatJetSelection::operator= (const FatJetSelection& ms) {
  if(this != &ms){    
    BaseSelection::operator = (ms);
    k_lqevent = ms.k_lqevent;  
  }
  return *this;
};

FatJetSelection::FatJetSelection(const FatJetSelection& ms):
  BaseSelection(ms)
{
  k_lqevent = ms.k_lqevent;  
};



