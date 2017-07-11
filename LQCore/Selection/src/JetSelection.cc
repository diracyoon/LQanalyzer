#include "JetSelection.h"
#include <iostream>


using namespace snu;

JetSelection::JetSelection(LQEvent ev) :BaseSelection() {
  k_lqevent = ev;
}

JetSelection::~JetSelection() {}

//// This code is used to make selection cuts to vectors of KJets


void JetSelection::BasicSelection(std::vector<KJet>& jetColl) {
  
  //// This is a basic set of cuts on jets

  std::vector<KJet> alljets = k_lqevent.GetJets();

  for (std::vector<KJet>::iterator jit = alljets.begin(); jit!=alljets.end(); jit++){
  
    if ( (jit->Pt() >= pt_cut_min) &&  (fabs(jit->Eta()) < eta_cut)){
      if ( PassUserID(PFJET_LOOSE, *jit) &&    (jit->Pt() >= pt_cut_min) &&  (fabs(jit->Eta()) < eta_cut))  jetColl.push_back(*jit);
    }
       
  }
}

void JetSelection::Selection(std::vector<KJet>& jetColl){

  std::vector<KJet> alljets = k_lqevent.GetJets();

  for (std::vector<KJet>::iterator jit = alljets.begin(); jit!=alljets.end(); jit++){

    bool IsNotPileUpJet = true;
    if(applypileuptool){
      IsNotPileUpJet = jit->PassPileUpMVA(PUJetIDWP);
    }


    if(apply_ID) {
      if ( jit->Pt() >= pt_cut_min && jit->Pt() < pt_cut_max &&
           fabs(jit->Eta()) < eta_cut
           && PassUserID(k_id, *jit) && IsNotPileUpJet)  jetColl.push_back(*jit);
    }
    else{
      if ( jit->Pt() >= pt_cut_min && jit->Pt() < pt_cut_max &&
           fabs(jit->Eta()) < eta_cut
           && PassUserID(PFJET_LOOSE, *jit)&& IsNotPileUpJet)  jetColl.push_back(*jit);
    }
  }

  BaseSelection::reset();
  return;

}  

void JetSelection::Selection(std::vector<KJet>& jetColl, bool LepVeto, std::vector<KMuon>& muonColl, std::vector<KElectron>& electronColl) {
 
  std::vector<KJet> alljets = k_lqevent.GetJets();
  
  std::vector<KJet> prejetColl; 

  for (std::vector<KJet>::iterator jit = alljets.begin(); jit!=alljets.end(); jit++){
    
    bool IsNotPileUpJet = true;
    if(applypileuptool){
      IsNotPileUpJet = jit->PassPileUpMVA(PUJetIDWP);
    }

    if(apply_ID) {
      if ( jit->Pt() >= pt_cut_min && jit->Pt() < pt_cut_max &&
	   fabs(jit->Eta()) < eta_cut
	   && PassUserID(k_id, *jit) && IsNotPileUpJet)  prejetColl.push_back(*jit);
    }
    else{
      if ( jit->Pt() >= pt_cut_min && jit->Pt() < pt_cut_max && 
	   fabs(jit->Eta()) < eta_cut
	   && PassUserID(PFJET_LOOSE, *jit) && IsNotPileUpJet)  prejetColl.push_back(*jit);
    }
  } 

  for (UInt_t ijet = 0; ijet < prejetColl.size(); ijet++) {
    jetIsOK = true;
    for (UInt_t ilep = 0; ilep < muonColl.size(); ilep++) {
      if (muonColl[ilep].DeltaR( prejetColl[ijet] ) < 0.4) {
	jetIsOK = false;  ilep = muonColl.size();
      }
    }/// End of muon loop
    for (UInt_t ilep = 0; ilep < electronColl.size(); ilep++) {
      if (electronColl[ilep].DeltaR( prejetColl[ijet] ) < 0.4 ) {
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


void JetSelection::SelectJets(std::vector<KJet>& jetColl,  vector<pair<TString, TString> > vids, vector<pair<TString, float> > vidf,float ptcut, float etacut) {

  std::vector<KJet> alljets = k_lqevent.GetJets();

  int icut(0);
  if (ptcut == -999. || etacut == -999.){
    for(unsigned int iv=0; iv < vidf.size(); iv++){
      if(!Check(vidf[iv].second)) continue;
      if (vidf[iv].first =="ptmin") { icut++; if(ptcut == -999.)ptcut=vidf[iv].second;}
      if (vidf[iv].first =="|etamax|") {icut++;  if (etacut == -999.)etacut=vidf[iv].second;}
      if(icut ==2) break;
    }
  }

  for (std::vector<KJet>::iterator jit = alljets.begin(); jit!=alljets.end(); jit++){

    bool pass_selection=true;
    if (!PassUserID(*jit, vids)) pass_selection=false;
    if ( (jit->Pt() >= ptcut)  && fabs(jit->Eta()) < etacut && pass_selection )  jetColl.push_back(*jit);
  }


} 

void JetSelection::SelectJets(std::vector<KJet>& jetColl, std::vector<KMuon> muonColl, std::vector<KElectron> electronColl,vector<pair<TString, TString> > vids, vector<pair<TString, float> > vidf,  float ptcut, float etacut) {
  
  std::vector<KJet> pre_jetColl; 
  std::vector<KJet> alljets = k_lqevent.GetJets();

  int icut(0);
  if (ptcut == -999. || etacut == -999.){
    for(unsigned int iv=0; iv < vidf.size(); iv++){
      if(!Check(vidf[iv].second)) continue;
      if (vidf[iv].first =="ptmin") { icut++; if(ptcut == -999.)ptcut=vidf[iv].second;}
      if (vidf[iv].first =="|etamax|") {icut++;  if (etacut == -999.)etacut=vidf[iv].second;}
      if(icut ==2) break;
    }
  }

  
  for (std::vector<KJet>::iterator jit = alljets.begin(); jit!=alljets.end(); jit++){
    
    bool pass_selection=true;
    if (!PassUserID(*jit, vids)) pass_selection=false;

    if ( (jit->Pt() >= ptcut)  && fabs(jit->Eta()) < etacut && pass_selection )  pre_jetColl.push_back(*jit);
  }


  for (UInt_t ijet = 0; ijet < pre_jetColl.size(); ijet++) {
    jetIsOK = true;
    for (UInt_t ilep = 0; ilep < muonColl.size(); ilep++) {
      if (muonColl[ilep].DeltaR( pre_jetColl[ijet] ) < 0.4) {
        jetIsOK = false;

	ilep = muonColl.size();
      }
    }/// End of muon loop
    for (UInt_t ilep = 0; ilep < electronColl.size(); ilep++) {
      if (electronColl[ilep].DeltaR( pre_jetColl[ijet] ) < 0.4 ) {
        jetIsOK = false;
        ilep = electronColl.size();
      }
    }/// End of electron loop
    
    if (jetIsOK) jetColl.push_back( pre_jetColl[ijet] );
  }/// End of Jet loop

}



void JetSelection::SelectJets(std::vector<KJet>& jetColl, std::vector<snu::KFatJet> fatjets , std::vector<KMuon> muonColl, std::vector<KElectron> electronColl,vector<pair<TString, TString> > vids, vector<pair<TString, float> > vidf,  float ptcut, float etacut) {

  std::vector<KJet> pre_jetColl;
  std::vector<KJet> alljets = k_lqevent.GetJets();

  int icut(0);
  if (ptcut == -999. || etacut == -999.){
    for(unsigned int iv=0; iv < vidf.size(); iv++){
      if(!Check(vidf[iv].second)) continue;
      if (vidf[iv].first =="ptmin") { icut++; if(ptcut == -999.)ptcut=vidf[iv].second;}
      if (vidf[iv].first =="|etamax|") {icut++;  if (etacut == -999.)etacut=vidf[iv].second;}
      if(icut ==2) break;
    }
  }

  ///https://indico.cern.ch/event/613931/contributions/2476973/attachments/1413452/2181248/update_WTagging_JMAR15022017_Simone.pdf
  /// see slides for delat R cuts between objetst
  for (std::vector<KJet>::iterator jit = alljets.begin(); jit!=alljets.end(); jit++){

    bool pass_selection=true;
    if (!PassUserID(*jit, vids)) pass_selection=false;

    if ( (jit->Pt() >= ptcut)  && fabs(jit->Eta()) < etacut && pass_selection )  pre_jetColl.push_back(*jit);
  }

  for (UInt_t ijet = 0; ijet < pre_jetColl.size(); ijet++) {
    jetIsOK = true;
  
    for (UInt_t ifjet = 0; ifjet < fatjets.size(); ifjet++) {
      if (fatjets[ifjet].DeltaR( pre_jetColl[ijet] ) < 0.8) {
        jetIsOK = false;

        ifjet = fatjets.size();
      }
    }/// End of fatjet loop     
    
    for (UInt_t ilep = 0; ilep < muonColl.size(); ilep++) {
      if (muonColl[ilep].DeltaR( pre_jetColl[ijet] ) < 0.3) {
        jetIsOK = false;

        ilep = muonColl.size();
      }
    }/// End of muon loop                                                                                                                                                    
    for (UInt_t ilep = 0; ilep < electronColl.size(); ilep++) {
      if (electronColl[ilep].DeltaR( pre_jetColl[ijet] ) < 0.3 ) {
        jetIsOK = false;
        ilep = electronColl.size();
      }
    }/// End of electron loop                                                                                                                                                

    if (jetIsOK) jetColl.push_back( pre_jetColl[ijet] );
  }/// End of Jet loop                                                                                                                                                       
}



bool JetSelection::PassUserID (ID id, snu::KJet jet){
  if      ( id == PFJET_LOOSE  ) return PassUserID_PFJetLoose  (jet);
  else if ( id == PFJET_TIGHT  ) return PassUserID_PFJetTight  (jet);
  else return false;
}


bool JetSelection::PassUserID (snu::KJet jet,vector<pair<TString, TString> > vids){ 

  for(unsigned int idel =0; idel < vids.size(); idel++){
    if(vids[idel].second == "false") continue;
    if(vids[idel].first == "pileup") { 
      if (!jet.PassPileUpMVA(vids[idel].second)) return false;
    }
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


bool JetSelection::PassUserID_PFJetLoose ( snu::KJet jet){
  
  return jet.PassLooseID();
}


bool JetSelection::PassUserID_PFJetTight ( snu::KJet jet)
{
  return jet.PassTightID();
}


JetSelection& JetSelection::operator= (const JetSelection& ms) {
  if(this != &ms){    
    BaseSelection::operator = (ms);
    k_lqevent = ms.k_lqevent;  
  }
  return *this;
};

JetSelection::JetSelection(const JetSelection& ms):
  BaseSelection(ms)
{
  k_lqevent = ms.k_lqevent;  
};



