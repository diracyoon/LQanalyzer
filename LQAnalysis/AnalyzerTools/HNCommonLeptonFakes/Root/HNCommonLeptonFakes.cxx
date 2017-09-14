// STL include(s):                                                                                                                                                       
#include <stdexcept>
#include <sstream>
#include <iostream>


// ROOT include(s):                                                                                                                                                                
#include "TFile.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TEfficiency.h"

/// local includes
#include "HNCommonLeptonFakes/HNCommonLeptonFakes.h"

using namespace std;


void HNCommonLeptonFakes::InitialiseFake(TString path, std::map<TString, std::pair<std::pair<TString,TString> ,std::pair<float,TString> > > fake_hist_config){
  
  if(_isinitialised) {  exit(1); }

  bool failedinitialisation=false;
  
  cout << "HNCommonLeptonFakes : Initializing" << endl;

  // Remember which TDirectory the system is in at the beginning of the method:                                                                                                      
  TDirectory* origDir = gDirectory;
  
  /// List files with fake rates

  string lqdir = getenv("LQANALYZER_DIR");
  int nhists(0);
  TDirectory* tempDir1 = getTemporaryDirectory();
  
  for(std::map<TString, std::pair<std::pair<TString,TString>,std::pair<float, TString> > >::iterator it=fake_hist_config.begin(); it!= fake_hist_config.end(); it++){
    TFile* file_fake   = TFile::Open( (path + it->second.first.first).Data());
    CheckFile(file_fake);

    TDirectory* tempDir1 = getTemporaryDirectory();
    tempDir1->cd();
    cout << it->second.first.second << endl;
    cout << dynamic_cast<TH2D*>((file_fake->Get(it->second.first.second))->Clone()) << endl;
    if(it->second.second.second=="TH2D") {
      cout << "Setting up " << it->first << endl;
      _2DEfficiencyMap_Double[it->first] = dynamic_cast<TH2D*>((file_fake->Get(it->second.first.second))->Clone());
    }

    cout << "Hist [Fake] name = " << it->second.first.second << endl;
    //cout << "X-axis rage = " << _2DEfficiencyMap_Double[it->first].find(it->first).second->GetXaxis()->GetBinLowEdge(1) << " - " << _2DEfficiencyMap_Double[it->first].find(it->first).second->GetXaxis()->GetBinLowEdge(_2DEfficiencyMap_Double[it->first].find(it->first).second->GetNbins()) << endl;
    _2DPtMaxMap_Double [it->first] = it->second.second.first;
    cout << "Setting  " << it->first << " " << it->second.second.first << endl;
    file_fake->Close();
    origDir->cd();
    nhists++;
  }

  if(nhists> 0)  _isinitialised = true;
  else _isinitialised=false;
  return;
}

bool HNCommonLeptonFakes::IsInitialised(){
  return _isinitialised;
}

// Constructor
HNCommonLeptonFakes::HNCommonLeptonFakes(){
  _isdata=false;
  NullTotals();
  DataPeriod = "B";
  _isinitialised=false;
  NullTotals();
  _minpt_muon=0.;
  _minpt_electron=0.;
  m_debug=false;
  
}

void HNCommonLeptonFakes::SetupFake(){

}

void HNCommonLeptonFakes::SetupFake(TString path, std::map<TString, std::pair<std::pair<TString,TString>  ,std::pair<float,TString> > >fake_hist_config){
  if(_isinitialised) return;
  InitialiseFake(path, fake_hist_config);
  
}

HNCommonLeptonFakes::HNCommonLeptonFakes(TString path, std::map<TString, std::pair<std::pair<TString,TString>  ,std::pair<float,TString> > >fake_hist_config  ){
  _isinitialised=false;
  InitialiseFake(path, fake_hist_config);
  m_debug=false;
  _isdata=false;
  NullTotals();
  _minpt_muon=0.;
  _minpt_electron=0.;

}

void HNCommonLeptonFakes::SetMinPtMuon(float minpt){
  _minpt_muon=minpt;
}

void HNCommonLeptonFakes::SetMinPtElectron(float minpt){
  _minpt_electron=minpt;
}


void HNCommonLeptonFakes::IsData(bool isdata){
  _isdata=isdata;
}

/// Destructor
HNCommonLeptonFakes::~HNCommonLeptonFakes() {
  CutHists.clear();
}
 

float HNCommonLeptonFakes::get_dilepton_mm_eventweight(bool geterr, std::vector<TLorentzVector> muons, bool ismu1tight, bool ismu2tight,TString frkey, float _mu1_val, float _mu2_val,int bkg_sample, TString prkey){

  if(muons.size() != 2) return 0.;

  //// vectors need to be ordered in pT                                                                                                                                                                                                           
  if(muons[0].Pt() < muons[1].Pt()) return -100000000000.;

  float _mu1_eta=fabs(muons.at(0).Eta());
  float _mu2_eta=fabs(muons.at(1).Eta());


  float fr1(0.),fr2(0.),r1(0.),r2(0.);

  // Need to create this function in HNCommonFakes                                                                                                                                                                                                
  r1=getPromptRate(_mu1_val, _mu1_eta, prkey);
  r2=getPromptRate(_mu2_val, _mu2_eta, prkey);
  
  fr1=  getFakeRate(_mu1_val, _mu1_eta,frkey,_minpt_muon,false);
  fr2=  getFakeRate(_mu2_val, _mu2_eta,frkey,_minpt_muon,false);

  if(geterr) {
    float fr1_err=  getFakeRate(_mu1_val, _mu1_eta,frkey,_minpt_muon,geterr);
    float fr2_err=  getFakeRate(_mu2_val, _mu2_eta,frkey,_minpt_muon,geterr);
    float ev_weight = CalculateDiLepMMWeight(fr1_err,fr1,fr2_err,fr2, ismu1tight, ismu2tight, 3);
    return ev_weight;

  }
  // Calculate event weight                                                                                                                                                                                                                      
  float ev_weight = CalculateDiLepMMWeight(r1,fr1,r2,fr2, ismu1tight, ismu2tight, bkg_sample);
  if(bkg_sample < 0) ev_weight = CalculateDiLepMMWeight(r1,fr1,r2,fr2, ismu1tight, ismu2tight);
  
  return ev_weight;
}


float HNCommonLeptonFakes::get_dilepton_ee_eventweight(bool geterr, std::vector<TLorentzVector> electrons, bool isel1tight, bool isel2tight,TString frkey, float value1, float value2, int bkg_sample, TString prkey){

  if(electrons.size() != 2) return 0.;

  float _el1_val=0.;
  float _el2_val=0.;
  if(value1==0.){
    _el1_val=electrons.at(0).Pt();
    _el2_val=electrons.at(1).Pt();
  }
  else {
    _el1_val=value1;
    _el2_val=value2;
  }
  //// vectors need to be ordered in pT                                                                                                                                         
  if(_el1_val < _el2_val) return -100000000000.;

  float _el1_eta=fabs(electrons.at(0).Eta());
  float _el2_eta=fabs(electrons.at(1).Eta());

  float fr1(0.),fr2(0.),r1(0.),r2(0.);

  // Need to create this function in HNCommonFakes                                                                                                                              
  r1=getPromptRate(_el1_val, _el1_eta, prkey);
  r2=getPromptRate(_el2_val, _el2_eta, prkey);
  
  fr1=  getFakeRate(_el1_val, _el1_eta,frkey,_minpt_electron);
  fr2=  getFakeRate(_el2_val, _el2_eta,frkey,_minpt_electron);

  if(geterr) {
    float fr1_err=  getFakeRate(_el1_val, _el1_eta,frkey,_minpt_electron,geterr);
    float fr2_err=  getFakeRate(_el2_val, _el2_eta,frkey,_minpt_electron,geterr);
    float ev_weight = CalculateDiLepMMWeight(fr1_err,fr1,fr2_err,fr2, isel1tight, isel2tight, 3);
    return ev_weight;

  }
  float ev_weight = CalculateDiLepMMWeight(r1,fr1,r2,fr2, isel1tight, isel2tight, bkg_sample);
  if(bkg_sample < 0) ev_weight = CalculateDiLepMMWeight(r1,fr1,r2,fr2, isel1tight, isel2tight);

  return ev_weight;
}

float HNCommonLeptonFakes::get_dilepton_em_eventweight(bool geterr, std::vector<TLorentzVector> electrons, std::vector<TLorentzVector> muons, bool iseltight, bool ismutight,TString frkey_el,TString frkey_mu, float el_value, float mu_value, int bkg_sample){

  TString prkey_el="";
  TString prkey_mu="";
  
  if(muons.size()!=1) {
    return (0.);
  }
  if(electrons.size()!=1) {
    return (0.);
  }

  float _lep1_val=0.;
  float _lep2_val=0.;
  bool _istight1=false;
  bool _istight2=false;
  TString key1="";
  TString key2="";

  float _el_eta=fabs(electrons.at(0).Eta());
  float _mu_eta=fabs(muons.at(0).Eta());

  float fr1(0.),fr2(0.),r1(0.),r2(0.);
  if(electrons.at(0).Pt() > muons.at(0).Pt()){
    _lep1_val=el_value;
    _lep2_val=mu_value;

    _istight1=iseltight;
    _istight2=ismutight;
    fr1=  getFakeRate(_lep1_val, _el_eta,frkey_el,_minpt_electron);
    fr2= getFakeRate(_lep2_val, _mu_eta,frkey_mu,_minpt_muon);
   
    r1 = getPromptRate(_lep1_val, _el_eta, prkey_el);
    r2 = getPromptRate(_lep2_val, _mu_eta, prkey_mu);


    if(geterr) {
      float fr1_err=  getFakeRate(_lep1_val, _el_eta,frkey_el,_minpt_electron,geterr);
      float fr2_err=  getFakeRate(_lep2_val, _mu_eta,frkey_mu,_minpt_muon,geterr);
      float ev_weight = CalculateDiLepMMWeight(fr1_err,fr1,fr2_err,fr2, _istight1, _istight2, 3);
      return ev_weight;
    }
  }
  else{
    
    _lep1_val=mu_value;
    _lep2_val=el_value;

    _istight2=iseltight;
    _istight1=ismutight;

    fr1= getFakeRate(_lep1_val, _mu_eta,frkey_mu,_minpt_muon);
    fr2=  getFakeRate(_lep2_val, _el_eta,frkey_el,_minpt_electron);
    r1 = getPromptRate(_lep1_val, _mu_eta, prkey_mu);
    r2 = getPromptRate(_lep2_val, _el_eta, prkey_el);

    if(geterr) {
      float fr1_err=  getFakeRate(_lep1_val, _mu_eta,frkey_mu,_minpt_muon,geterr);
      float fr2_err=  getFakeRate(_lep2_val, _el_eta,frkey_el,_minpt_electron,geterr);
      float ev_weight = CalculateDiLepMMWeight(fr1_err,fr1,fr2_err,fr2, _istight1, _istight2, 3);
      return ev_weight;
    }
  }

  float ev_weight = CalculateDiLepMMWeight(r1,fr1,r2,fr2, _istight1, _istight2, bkg_sample);
  if(bkg_sample < 0) ev_weight = CalculateDiLepMMWeight(r1,fr1,r2,fr2, _istight1, _istight2);

  return ev_weight;

}



float HNCommonLeptonFakes::getPromptRate(float val, float eta, TString tag){

  float eff_prompt=1.;
  if (tag == "") return 1.;

  if(fabs(eta) >= 2.5) eta = 2.4;

  map<TString,TH2D*>::const_iterator mapit;
  
  mapit = _2DEfficiencyMap_Double.find(tag);
  if(mapit!=_2DEfficiencyMap_Double.end()){

    map<TString,float>::const_iterator ptmapit=_2DPtMaxMap_Double.find(tag);
    if(ptmapit==_2DPtMaxMap_Double.end())  return 1.;
   
    if(val < 15.) val = 16.;
    if(val >= ptmapit->second) val = ptmapit->second-1.;

    int binx =  mapit->second->FindBin(val,fabs(eta));
    eff_prompt =  mapit->second->GetBinContent(binx);
  }
  else if (tag == "") return 1.;
  
  return eff_prompt;
}


float HNCommonLeptonFakes::getFakeRate(float val, float eta , TString tag, float ptmin, bool geterr){
  
  float eff_fake=0.;
  if(fabs(eta) > 2.5) return -9999999.;
  if(val < ptmin) return -999999.;

  
  map<TString,float>::const_iterator ptmapit;
  ptmapit=_2DPtMaxMap_Double.find(tag);

  

  if(val >= ptmapit->second) val = ptmapit->second-1.;

  map<TString,TH2D*>::const_iterator mapit;
  mapit = _2DEfficiencyMap_Double.find(tag);
  if(mapit!=_2DEfficiencyMap_Double.end()){
    int binx =  mapit->second->FindBin(val, fabs(eta));
    eff_fake =  float(mapit->second->GetBinContent(binx));
    if(geterr) eff_fake =  float(mapit->second->GetBinError(binx));
  }
  else NoHist((tag));
  return eff_fake;
  
}




void HNCommonLeptonFakes::NoHist(TString hist, TString function){
  
  cout << function << ": No Histogram found for string " << hist << endl;

  cout << "Possible _EfficiencyMap hists are : " << endl;
  for(map<TString, TH1*>::iterator mit = _EfficiencyMap.begin(); mit != _EfficiencyMap.end(); mit++){
    cout << mit->first << endl;  
  }

  cout << "Possible _2DEfficiencyMap_Double hists are : " << endl;
  for(map<TString, TH2D*>::iterator mit = _2DEfficiencyMap_Double.begin(); mit != _2DEfficiencyMap_Double.end(); mit++){
    cout << mit->first <<" " << mit->second << endl;
  }
  

  exit(1);
}

void HNCommonLeptonFakes::NullTotals(){
  
}


void HNCommonLeptonFakes::CheckFile(TFile* file){

  if(file) cout << "HNCommonLeptonFakes: File " << file->GetName() << " was found." << endl;
  else cout << "HNCommonLeptonFakes  " << file->GetName()  << "  : ERROR Rootfile failed to open." << endl;
  
  if(!file) exit(0);
  return;
}



float HNCommonLeptonFakes::CalculateDiLepMMWeight(float r1, float fr1, float r2, float fr2, bool mu1Tight, bool mu2Tight){

  
  // Function calculates event weight given r/f of two muons in the event
  double alpha = 1./((1.- fr1)*(1.- fr2));

  // Initialise weight
  float w_mm=-999.;
    
  // Terms for RF and FR events from MM
  double termTT = 0.;
  double termTL = alpha*(fr2*(1.-fr1));
  double termLT = alpha*(fr1*(1.-fr2));
  double termLL = -2.*alpha*(fr2*(fr1));
  
  
  // Term for FF events
  termTT += 0.;
  termTL += 0.;
  termLT += 0.;
  termLL += alpha*(fr1*fr2);

  
  // decide which type of event we have so to return correct weight
  bool isLL = (!mu1Tight&&!mu2Tight);
  bool isLT = (!mu1Tight&&mu2Tight);
  bool isTL = (mu1Tight && !mu2Tight);
  bool isTT = (mu1Tight && mu2Tight);

  if(isLL)w_mm = termLL;
  if(isTL)w_mm = termTL;
  if(isLT)w_mm = termLT;
  if(isTT)w_mm=  termTT;  



  return w_mm;   
}


float HNCommonLeptonFakes::CalculateDiLepMMWeight(float fr1_err, float fr1, float fr2_err, float fr2, bool mu1Tight, bool mu2Tight, int eventtype){


  // Function calculates event weight given r/f of two muons in the event
  double alpha = 1./((1.- fr1)*(1.- fr2));
  
  // Initialise weight
  float w_mm=-999.;


  // Terms for RF and FR events from MM
  double termTT = 0.;
  double termTL = alpha*(fr2*(1.-fr1));
  double termLT = alpha*(fr1*(1.-fr2));
  double termLL = -2.*alpha*(fr2*(fr1));

  if(eventtype == 0){
    termTT += 0.;
    termTL += 0.;
    termLT += 0.;
    termLL += alpha*(fr1*fr2);
  }

  if(eventtype == 2){
    // Term for FF events
    termTT = 0.;
    termTL = 0.;
    termLT = 0.;
    termLL = alpha*(fr1*fr2);
  }

  // decide which type of event we have so to return correct weight
  bool isLL = (!mu1Tight&&!mu2Tight);
  bool isLT = (!mu1Tight&&mu2Tight);
  bool isTL = (mu1Tight && !mu2Tight);
  bool isTT = (mu1Tight && mu2Tight);

  
  if(eventtype == 3) {
    // error onfake
    float error = 0.;
    if(isLL) {
      error = sqrt( pow(fr1_err*fr2*(1-fr2),2) + pow(fr2_err*fr1*(1-fr1),2) ) / pow((1-fr1)*(1-fr2),2);
    }
    if(isLT){
      error= fr1_err/pow(fr1-1,2);
    }
    if(isTL){
      error= fr2_err/pow(fr2-1,2);
    }
    return error;
  }
    
  if(isLL)w_mm = termLL;
  if(isTL)w_mm = termTL;
  if(isLT)w_mm = termLT;
  if(isTT)w_mm=  termTT;


  return w_mm;
}



double HNCommonLeptonFakes::lepton_weight(bool loose, double r, double f) {
  // Function to get event weight for a single muon event
  double nl=loose ? 1 :0;
  double nt=!loose ? 1 :0;
  if(r==f) return 0.;
  double a=1./(r-f);
  float nf=a*f*r*nl;
  float nr=a*f*(r - 1.)*nt;
  return nf+nr;
}


/// ====== START TRIMUON FUNCTIONALITY
void HNCommonLeptonFakes::APPLY_scalefactor(bool applysf) {
  apply_scalefactor = applysf;
}



TDirectory* HNCommonLeptonFakes::getTemporaryDirectory(void) const
{
  //                                                                                                                                                                                
  // Create a unique directory in memory to hold the histograms:                                                                                                                    
  //                                                                                                                                                                               
  gROOT->cd();
  TDirectory* tempDir = 0;
  int counter = 0;
  while (not tempDir) {
    // First, let's find a directory name that doesn't exist yet:                                                                                                                    
    std::stringstream dirname;
    dirname << "HNCommonLeptonFakes_%i" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    // Let's try to make this directory:                                                                                                                                             
    tempDir = gROOT->mkdir((dirname.str()).c_str());

  }

  return tempDir;
}

