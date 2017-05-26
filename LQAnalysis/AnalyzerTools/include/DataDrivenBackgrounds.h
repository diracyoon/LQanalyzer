#ifndef DataDrivenBackgrounds_h
#define DataDrivenBackgrounds_h

#include <vector>
#include <string>
#include <map>
#include <sstream>  
#include <fstream>      // std::ifstream
#include <iostream> 

class TString;

/// SNU codes
#include "KMuon.h"
#include "KElectron.h"
#include "KPhoton.h"
#include "KGenJet.h"
#include "KJet.h"
#include "KFatJet.h"
#include "KTruth.h"
#include "KTrigger.h"
#include "KEvent.h"


//forward declarations                                                                                                                                                              

class EventBase;


// ROOT 
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"

// local includes
#include "HNCommonLeptonFakes/HNCommonLeptonFakes/HNCommonLeptonFakes.h"

class DataDrivenBackgrounds{

 public:
    
  DataDrivenBackgrounds();
  ~DataDrivenBackgrounds();

  void SetMCPeriod(int mcperiod);
  void SetIsData(bool isdata);
  void PrintSummary();

  void SetEventBase(EventBase* ev);
  void CheckEventBase();

  //// Charge flip HN
  float WeightCFEvent(std::vector<snu::KElectron> electrons, bool runchargeflip);
  float CFRate(snu::KElectron el);
  float CFRate_Run2(snu::KElectron el, TString el_id);
  std::vector<snu::KElectron>  ShiftElectronEnergy(std::vector<snu::KElectron> el, bool applyshift);


  /// Fake Backgrounds e+m HN

  /// == Dilepton
  float Get_DataDrivenWeight_EE(bool geterr, std::vector<snu::KElectron> k_electrons);
  float Get_DataDrivenWeight_EEmva(bool geterr,std::vector<snu::KElectron> k_electrons, bool tight1, bool tight2, TString key1, TString key2);
  float Get_DataDrivenWeight_EE(bool geterr, std::vector<snu::KElectron> k_electrons, TString IDloose,TString IDtight, TString method);
  float Get_DataDrivenWeight_MM(bool geterr, std::vector<snu::KMuon> k_muons);
  float Get_DataDrivenWeight_MM(bool geterr, std::vector<snu::KMuon> k_muons,  TString ID, TString method);
  float Get_DataDrivenWeight_EM(bool geterr,std::vector<snu::KMuon> k_muons, std::vector<snu::KElectron> k_electrons,  TString IDe, TString IDm,  TString method);

  float GetFakeRateEl(float pt, float eta, TString cut);

  /// ===  trilepton
  float Get_DataDrivenWeight_MMM(bool geterr, std::vector<snu::KMuon> k_muons);

  TString GetElFRKey(TString elidloose, TString elidtight, TString method);
  std::vector<TString> GetElFRKey( TString IDloose,TString IDtight, TString method, std::vector<TString> regs1);
  //==== General lepton fakes
  float Get_DataDrivenWeight(bool geterr, std::vector<snu::KMuon> k_muons, TString muid, int n_muons, std::vector<snu::KElectron> k_electrons, TString elid, int n_electrons, TString elidloose="ELECTRON16_POG_FAKELOOSE", TString method="dijet_ajet40");

  /// = single lepton
  float Get_DataDrivenWeight_E(bool geterr, std::vector<snu::KElectron> k_electrons ,TString ID, TString method);
  float Get_DataDrivenWeight_M(bool geterr, std::vector<snu::KMuon> k_muons,  TString ID, TString method);
  

  void Test();

  // Class object to get event weights for Fakes

  /// GENERAL FUNCTIONS
  std::vector<TLorentzVector> MakeTLorentz( std::vector<snu::KElectron> el);
  std::vector<TLorentzVector> MakeTLorentz( std::vector<snu::KMuon> mu);
  std::vector<TLorentzVector> MakeTLorentz( std::vector<snu::KJet> jet);
  std::vector<TLorentzVector> MakeTLorentz( std::vector<snu::KFatJet> jet);


  inline HNCommonLeptonFakes*   GetFakeObj() const {return m_fakeobj;}
  inline EventBase*  GetEventBase() const {return  dd_eventbase;}


 private:
  bool corr_isdata;
  double k_mcperiod;

  HNCommonLeptonFakes* m_fakeobj;
  EventBase* dd_eventbase;



};
#endif
