#include "BaseSelection.h"

BaseSelection::BaseSelection() {
 
  dz_cut = 1000000.;
  pt_cut_max = 1000000.;
  pt_cut_min = 0.;
  eta_cut_min = -100.;
  eta_cut = 100.;
  jpt_cut_max = 1000000.;
  jpt_cut_min = 0.;
  jeta_cut_min = -100.;
  jeta_cut = 100.;
  relIso_cut  = 10000000.;
  relIsoMIN_cut = 0.;
  chiNdof_cut = 10000.;
  chiNdofMIN_cut = 0.;
  dxy_cut = 100000.;
  dxyMIN_cut = -1000000.;
  k_id = ELECTRON_HN_TIGHT;
  
  apply_ptcut=false;
  apply_etacut=false;
  apply_jptcut=false;
  applypileuptool= false;
  apply_jetacut=false;
  apply_relisocut=false;
  apply_chi2cut=false;
  apply_dxycut=false;
  apply_dzcut=false;
  apply_general=false;
  apply_deposit=false;
  
  apply_ID= false;
  apply_convcut= false;
  apply_chargeconst = false;

  apply_BETrRegIncl = false;
  apply_BESepCut    = false;
  apply_relisocut   = false;
  apply_HLTSafeCut  = false;     HLTSafeLevel     = "CaloIdL_TrackIdL_IsoVL";
  apply_dxysigmin   = false;     apply_dxysigmax  = false;
  relIsoBarrel_max  = 100000.;   relIsoEndcap_max = 100000.;
  relIsoBarrel_min  = 0.     ;   relIsoEndcap_min = 0.     ;
  dxyBarrel_max     = 100000.;   dxyEndcap_max    = 100000.;
  dzBarrel_max      = 100000.;   dzEndcap_max     = 100000.;
  dxySig_min        = 0.     ;   dxySig_max       = 100000.;
  RelIsoType        = "Default";

  PUJetIDWP         = "Loose";
}

////////// PREDEFINED MUON SELECTIONS                                                                                                        


bool  BaseSelection::Check(float val){

  if(abs(val) == 999.) return false;
  
  return true;
}

TString BaseSelection::GetString(ID id){
  if ( id == MUON_HN_VETO) return "MUON_HN_VETO";
  if ( id == MUON_HN_FAKELOOSE) return "MUON_HN_FAKELOOSE";
  if ( id == MUON_HN_TIGHT) return  "MUON_HN_TIGHT";
  if ( id == MUON_POG_LOOSE) return "MUON_POG_LOOSE";
  if ( id == MUON_POG_MEDIUM) return "MUON_POG_MEDIUM";
  if ( id == MUON_POG_TIGHT) return  "MUON_POG_TIGHT";
  if ( id == ELECTRON_POG_MVA_WP80 ) return "ELECTRON_POG_MVA_WP80";
  if ( id == ELECTRON_POG_MVA_WP90 ) return "ELECTRON_POG_MVA_WP90"; 
  if ( id == ELECTRON_HctoWA_FAKELOOSE ) return "ELECTRON_HctoWA_FAKELOOSE";


  if ( id == MUON_TOP_VETO) return  "MUON_TOP_VETO";
  if ( id == MUON_TOP_LOOSE) return  "MUON_TOP_LOOSE";
  if ( id == MUON_TOP_TIGHT) return "MUON_TOP_TIGHT";
  if ( id == ELECTRON_HN_VETO) return "ELECTRON_HN_VETO";
  if ( id == ELECTRON_HN_FAKELOOSE) return "ELECTRON_HN_FAKELOOSE";
  if ( id == ELECTRON_HN_FAKELOOSE_NOD0) return "ELECTRON_HN_FAKELOOSE_NODO";
  if ( id == ELECTRON_HN_MVA_LOOSE) return "ELECTRON_HN_MVA_LOOSE";
  if ( id == ELECTRON_HN_MVA_TIGHT) return "ELECTRON_HN_MVA_TIGHT";
  if ( id == ELECTRON_HN_TIGHT) return  "ELECTRON_HN_TIGHT";
  if ( id == ELECTRON_TOP_VETO) return  "ELECTRON_TOP_VETO";
  if ( id == ELECTRON_TOP_LOOSE) return  "ELECTRON_TOP_LOOSE";
  if ( id == ELECTRON_TOP_TIGHT) return "ELECTRON_TOP_TIGHT";
  if ( id == ELECTRON_NOCUT) return "ELECTRON_NOCUT";
  if ( id == ELECTRON_POG_VETO) return "ELECTRON_POG_VETO";
  if ( id == ELECTRON_POG_LOOSE) return "ELECTRON_POG_LOOSE";
  if ( id == ELECTRON_POG_MEDIUM) return "ELECTRON_POG_MEDIUM";
  if ( id == ELECTRON_POG_TIGHT) return  "ELECTRON_POG_TIGHT";

  
  cerr << "Did not find ID for object. Please enter TString of ID not enum" << endl;
  exit(EXIT_FAILURE);

}


void BaseSelection::SetIDSMap(std::map<TString, vector<pair<TString, TString> > > smap){


  k_stringmap = smap;

}

void BaseSelection::SetIDFMap(std::map<TString, vector<pair<TString,float> > > fmap){
  k_floatmap = fmap;

}

void BaseSelection::PrintSIDMap(){


  for(std::map<TString, vector<pair<TString,TString> > >::iterator it = k_stringmap.begin() ; it != k_stringmap.end(); it++){

    cout << it->first << " : " << it->second.size() << endl;
  }
}

void BaseSelection::PrintFIDMap(){


  for(std::map<TString, vector<pair<TString,float> > >::iterator it = k_floatmap.begin() ; it != k_floatmap.end();it++){

    cout << it->first << " : " << it->second.size() << endl;
  }
  return;
}


vector<pair<TString,TString> > BaseSelection::GetStringList(TString id){
  std::map<TString, vector<pair<TString,TString> > >::iterator it = k_stringmap.find(id);
  if(it == k_stringmap.end()){
    cerr << "Did not find string " + id +". Please check string in GetBases"<< endl;
    exit(EXIT_FAILURE);
  }
  
  return it->second;

}

vector<pair<TString,float> > BaseSelection::GetFloatList(TString id){
  std::map<TString, vector<pair<TString,float> > >::iterator it = k_floatmap.find(id);
  if(it == k_floatmap.end()){
    cerr << "Did not find string " + id +". Please check string in GetBases"<< endl;
    exit(EXIT_FAILURE);
  }

  return it->second;

}

TString BaseSelection::AccessStringMap(TString label, TString id){

  std::map<TString, vector<pair<TString,TString> > >::iterator it = k_stringmap.find(id);
  if(it == k_stringmap.end()){
    cerr << "Did not find string " + id +". Please check string in GetBases"<< endl;
    exit(EXIT_FAILURE);
  }
  else{
    for (unsigned int i=0; i  < it->second.size(); i++){
      if ( it->second.at(i).first == label) return it->second.at(i).second;
    }
  }
  return "";
}

float BaseSelection::AccessFloatMap(TString label, TString id){

  std::map<TString, vector<pair<TString,float> > >::iterator it = k_floatmap.find(id);
  if(it == k_floatmap.end()){
    cout << " k_floatmap.size() = " <<  k_floatmap.size() << endl;
    cerr << "Did not find ID " + id + " Please check string in GetBases"<< endl;
    exit(EXIT_FAILURE);
  }
  else{
    for(unsigned int i=0; i  <it->second.size(); i++){
      if ( it->second.at(i).first == label) return it->second.at(i).second;
    }
  }
  return -999.;
}


bool BaseSelection::CheckCutString(TString label, TString id){

  if(AccessStringMap(label,id).Contains("true")) return true;
  else return false;

}



TString BaseSelection::GetCutString(TString label, TString id){

  return AccessStringMap(label,id);


}


bool BaseSelection::CheckCutFloat(TString label, TString id){

  if(AccessFloatMap(label,id) == -999.) return false;
  if(AccessFloatMap(label,id) == 999.) return false;
  else return true;

}


BaseSelection& BaseSelection::operator= (const BaseSelection& bs) 
{
  if(this != &bs){
    k_stringmap = bs.k_stringmap;
    k_floatmap = bs.k_floatmap;
    dz_cut = bs.dz_cut;
    pt_cut_max = bs.pt_cut_max;  
    pt_cut_min = bs.pt_cut_min; 
    eta_cut_min = bs.eta_cut_min; 
    eta_cut = bs.eta_cut; 
    jpt_cut_max = bs.jpt_cut_max;
    jpt_cut_min = bs.jpt_cut_min;
    jeta_cut_min = bs.jeta_cut_min; 
    jeta_cut = bs.jeta_cut; 
    relIso_cut  = bs.relIso_cut; 
    relIsoMIN_cut = bs.relIsoMIN_cut;
    chiNdof_cut = bs.chiNdof_cut; 
    chiNdofMIN_cut = bs.chiNdofMIN_cut;
    dxy_cut = bs.dxy_cut; 
    dxyMIN_cut = bs.dxyMIN_cut;
    k_id = bs.k_id;
    applypileuptool=bs.applypileuptool;
    apply_ptcut= bs.apply_ptcut;
    apply_etacut= bs.apply_etacut;
    apply_jptcut= bs.apply_jptcut;
    apply_jetacut= bs.apply_jetacut;
    apply_relisocut= bs.apply_relisocut;
    apply_chi2cut= bs.apply_chi2cut;
    apply_dxycut= bs.apply_dxycut;
    apply_dzcut= bs.apply_dzcut;
    apply_general= bs.apply_general;
    apply_deposit= bs.apply_deposit;
    apply_ID= bs.apply_ID;
    apply_convcut= bs.apply_convcut;
    apply_chargeconst = bs.apply_chargeconst;

    apply_BETrRegIncl = bs.apply_BETrRegIncl;
    apply_BESepCut    = bs.apply_BESepCut;
    apply_relisocut   = bs.apply_relisocut;
    apply_HLTSafeCut  = bs.apply_HLTSafeCut; HLTSafeLevel     = bs.HLTSafeLevel;
    apply_dxysigmin   = bs.apply_dxysigmin;  apply_dxysigmax  = bs.apply_dxysigmax;
    relIsoBarrel_max  = bs.relIsoBarrel_max; relIsoEndcap_max = bs.relIsoEndcap_max;
    relIsoBarrel_min  = bs.relIsoBarrel_min; relIsoEndcap_min = bs.relIsoEndcap_min;
    dxyBarrel_max     = bs.dxyBarrel_max;    dxyEndcap_max    = bs.dxyEndcap_max;
    dzBarrel_max      = bs.dzBarrel_max;     dzEndcap_max     = bs.dzEndcap_max;
    dxySig_min        = bs.dxySig_min;       dxySig_max       = bs.dxySig_max;
    RelIsoType        = bs.RelIsoType;
  
    PUJetIDWP         = bs.PUJetIDWP;
  }
  return *this;
}


BaseSelection::BaseSelection(const BaseSelection& bs) {
  
  k_stringmap = bs.k_stringmap;
  k_floatmap = bs.k_floatmap;
  dz_cut = bs.dz_cut;
  pt_cut_max = bs.pt_cut_max;  
  pt_cut_min = bs.pt_cut_min; 
  eta_cut_min = bs.eta_cut_min; 
  eta_cut = bs.eta_cut; 
  jpt_cut_max = bs.jpt_cut_max;
  jpt_cut_min = bs.jpt_cut_min;
  jeta_cut_min = bs.jeta_cut_min; 
  jeta_cut = bs.jeta_cut; 
  relIso_cut  = bs.relIso_cut; 
  relIsoMIN_cut = bs.relIsoMIN_cut;
  chiNdof_cut = bs.chiNdof_cut; 
  chiNdofMIN_cut = bs.chiNdofMIN_cut;
  dxy_cut = bs.dxy_cut; 
  dxyMIN_cut = bs.dxyMIN_cut;
  k_id = bs.k_id;
  applypileuptool=bs.applypileuptool;
  apply_ptcut= bs.apply_ptcut;
  apply_etacut= bs.apply_etacut;
  apply_jptcut= bs.apply_jptcut;
  apply_jetacut= bs.apply_jetacut;
  apply_relisocut= bs.apply_relisocut;
  apply_chi2cut= bs.apply_chi2cut;
  apply_dxycut= bs.apply_dxycut;
  apply_dzcut= bs.apply_dzcut;
  apply_general= bs.apply_general;
  apply_deposit= bs.apply_deposit;
  apply_ID= bs.apply_ID;
  apply_convcut= bs.apply_convcut;
  apply_chargeconst = bs.apply_chargeconst;
  
  apply_BETrRegIncl = bs.apply_BETrRegIncl;
  apply_BESepCut    = bs.apply_BESepCut;
  apply_relisocut   = bs.apply_relisocut;
  apply_HLTSafeCut  = bs.apply_HLTSafeCut; HLTSafeLevel     = bs.HLTSafeLevel;
  apply_dxysigmin   = bs.apply_dxysigmin;  apply_dxysigmax  = bs.apply_dxysigmax;
  relIsoBarrel_max  = bs.relIsoBarrel_max; relIsoEndcap_max = bs.relIsoEndcap_max;
  relIsoBarrel_min  = bs.relIsoBarrel_min; relIsoEndcap_min = bs.relIsoEndcap_min;
  dxyBarrel_max     = bs.dxyBarrel_max;    dxyEndcap_max    = bs.dxyEndcap_max;
  dzBarrel_max      = bs.dzBarrel_max;     dzEndcap_max     = bs.dzEndcap_max;
  dxySig_min        = bs.dxySig_min;       dxySig_max       = bs.dxySig_max;
  RelIsoType        = bs.RelIsoType;

  PUJetIDWP         = bs.PUJetIDWP;
}



BaseSelection::~BaseSelection() {}

void BaseSelection::reset(){
  dz_cut = 1000000.;
  pt_cut_max = 1000000.;
  pt_cut_min = 0.;
  jpt_cut_max = 1000000.;
  jpt_cut_min = 0.;
  eta_cut_min = -100.;
  eta_cut = 100.;
  jeta_cut_min = -100.;
  jeta_cut = 100.;
  relIso_cut  = 10000000.;
  relIsoMIN_cut = 0.;
  chiNdof_cut = 10000.;
  chiNdofMIN_cut = 0.;
  dxy_cut = 100000.;
  dxyMIN_cut = -1000000.;
  k_id = ELECTRON_HN_TIGHT;
  applypileuptool=false;
  apply_ptcut=false;
  apply_etacut=false;
  apply_jptcut=false;
  apply_jetacut=false;
  apply_relisocut=false;
  apply_chi2cut=false;
  apply_dxycut=false;
  apply_dzcut=false;
  apply_ID= false;
  apply_convcut= false;
  apply_chargeconst = false;

  apply_BETrRegIncl = false;
  apply_BESepCut    = false;
  apply_relisocut   = false;
  apply_HLTSafeCut  = false;     HLTSafeLevel     = "CaloIdL_TrackIdL_IsoVL";
  apply_dxysigmin   = false;     apply_dxysigmax  = false;
  relIsoBarrel_max  = 100000.;   relIsoEndcap_max = 100000.;
  relIsoBarrel_min  = 0.;        relIsoEndcap_min = 0.;
  dxyBarrel_max     = 100000.;   dxyEndcap_max    = 100000.;
  dzBarrel_max      = 100000.;   dzEndcap_max     = 100000.;
  dxySig_min        = 0.     ;   dxySig_max       = 100000.;
  RelIsoType        = "Default";

  PUJetIDWP         = "Loose";
}

void BaseSelection::SetPt(Double_t minPt) {
  apply_ptcut=true;
  minPt ? pt_cut_min=minPt : pt_cut_min=0.0;
  pt_cut_max=10000.0;

}

void BaseSelection::SetPt(Double_t minPt, Double_t maxPt) {
  apply_ptcut=true;
  minPt ? pt_cut_min=minPt : pt_cut_min=0.0;
  maxPt ? pt_cut_max=maxPt : pt_cut_max=10000.0;
}

void BaseSelection::SetID(ID id){
  apply_ID = true;
  k_id = id;
}


void BaseSelection::SetCheckCharge(bool check){
  apply_chargeconst = check;
}

void BaseSelection::SetUseJetPileUp(bool use){
  applypileuptool=use;
}

void BaseSelection::SetPileUpJetID(bool use, TString pujetidwp){
  applypileuptool=use;
  PUJetIDWP=pujetidwp;
}

void BaseSelection::SetApplyConvVeto(bool apply){
  apply_convcut=apply;
}

void BaseSelection::SetHLTSafeCut(TString level){
  apply_HLTSafeCut = true;  HLTSafeLevel = level; 
}

void BaseSelection::SetJetPt(Double_t minPt) {
  apply_jptcut=true;
  minPt ? jpt_cut_min=minPt : jpt_cut_min=0.0;
  jpt_cut_max=10000.0;
}



void BaseSelection::SetEta(Double_t Eta) {
  apply_etacut=true;
  Eta ? eta_cut=Eta : eta_cut=3.0;
  eta_cut_min=0.0;
  
}


void BaseSelection::SetJetEta(Double_t Eta) {
  apply_jetacut=true;
  Eta ? jeta_cut=Eta : jeta_cut=3.0;
  jeta_cut_min=0.0;
}

void BaseSelection::SetEta(Double_t minEta, Double_t Eta) {
  apply_etacut=true;
  minEta ? eta_cut_min=minEta : eta_cut_min=0.0;
  Eta ? eta_cut=Eta : eta_cut=3.0;
}

void BaseSelection::SetRelIso(Double_t RelIso) {
  apply_relisocut=true;
  RelIso ? relIso_cut=RelIso : relIso_cut=10.0;
  relIsoMIN_cut=0.0;
}

void BaseSelection::SetRelIso(Double_t RelIsoMIN, Double_t RelIso) {
  apply_relisocut=true;
  RelIsoMIN ? relIsoMIN_cut=RelIsoMIN : relIsoMIN_cut=0.0;
  RelIso ? relIso_cut=RelIso : relIso_cut=10.0;
}

void BaseSelection::SetChiNdof(Double_t ChiNdof) {
  apply_chi2cut=true;
  ChiNdof ? chiNdof_cut=ChiNdof : chiNdof_cut=100.0;
  chiNdofMIN_cut=0.0;
}

void BaseSelection::SetChiNdof(Double_t ChiNdofMIN, Double_t ChiNdof) {
  apply_chi2cut=true;
  ChiNdofMIN ? chiNdofMIN_cut=ChiNdofMIN : chiNdofMIN_cut=0.0;
  ChiNdof ? chiNdof_cut=ChiNdof : chiNdof_cut=100.0;
}

void BaseSelection::SetBSdxy(Double_t set_dxy) {
  apply_dxycut=true;
  set_dxy ? dxy_cut=set_dxy : dxy_cut=2.0;
  dxyMIN_cut=0.0;
}

void BaseSelection::SetBSdxy(Double_t dxyMIN, Double_t set_dxy) {
  apply_dxycut=true;
   dxyMIN ? dxyMIN_cut=dxyMIN : dxyMIN_cut=0.0;
   set_dxy ? dxy_cut=dxy : dxy_cut=2.0;
}

void BaseSelection::SetBSdz(Double_t set_dz) {
  apply_dzcut=true;
  set_dz ? dz_cut=set_dz : dz_cut=2.0;
}

void BaseSelection::SetdxySigMin(Double_t dxySigMin){
  apply_dxysigmin=true;
  dxySigMin ? dxySig_min=dxySigMin : 0.;
}

void BaseSelection::SetdxySigMax(Double_t dxySigMax){
  apply_dxysigmax=true;
  dxySigMax ? dxySig_max=dxySigMax : 10000.;
}

void BaseSelection::SetBETrRegIncl(bool apply){
  apply_BETrRegIncl=apply;
}

void BaseSelection::SetRelIsoBEMax(Double_t BarrelIso, Double_t EndcapIso){
  apply_BESepCut=true;
  apply_relisocut=true;
  BarrelIso ? relIsoBarrel_max=BarrelIso : 10.0;
  EndcapIso ? relIsoEndcap_max=EndcapIso : 10.0;
}

void BaseSelection::SetRelIsoBEMin(Double_t BarrelIso, Double_t EndcapIso){
  apply_BESepCut=true;
  apply_relisocut=true;
  BarrelIso ? relIsoBarrel_min=BarrelIso : 0.0;
  EndcapIso ? relIsoEndcap_min=EndcapIso : 0.0;
}

void BaseSelection::SetdxyBEMax(Double_t Barreldxy, Double_t Endcapdxy){
  apply_BESepCut=true;
  apply_dxycut=true;
  Barreldxy ? dxyBarrel_max=Barreldxy : 2.0;
  Endcapdxy ? dxyEndcap_max=Endcapdxy : 2.0;
}

void BaseSelection::SetdzBEMax(Double_t Barreldz, Double_t Endcapdz){
  apply_BESepCut=true;
  apply_dzcut=true;
  Barreldz ? dzBarrel_max=Barreldz : 2.0;
  Endcapdz ? dzEndcap_max=Endcapdz : 2.0;
}

void BaseSelection::SetRelIsoType(TString type){
  apply_relisocut=true;
  RelIsoType=type;
}
