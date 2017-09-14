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
#include "HNCommonLeptonFakes/HNCommonLeptonFakesTriLep/HNCommonLeptonFakesTriLep.h"

class DataDrivenBackgrounds{

 public:
    
  DataDrivenBackgrounds();
  DataDrivenBackgrounds(bool setupown);
  ~DataDrivenBackgrounds();

  void SetMCPeriod(int mcperiod);
  void SetIsData(bool isdata);
  void PrintSummary();

  void SetEventBase(EventBase* ev);
  void CheckEventBase();

  bool  SetupFake(TString path, std::map<TString, std::pair<std::pair<TString,TString>  ,std::pair<float,TString> > >fake_hist_config);
  bool  SetupFake();

  //// Charge flip HN
  float WeightCFEvent(std::vector<snu::KElectron> electrons, bool runchargeflip);
  float CFRate(snu::KElectron el);
  float CFRate_Run2(snu::KElectron el, TString el_id);
  std::vector<snu::KElectron>  ShiftElectronEnergy(std::vector<snu::KElectron> el, bool applyshift);


  /// Fake Backgrounds e+m HN

  /// == Dilepton
  float Get_DataDrivenWeight_EE(bool geterr, std::vector<snu::KElectron> k_electrons, TString ID, TString variable,TString tag, int FakeType=0);
  float Get_DataDrivenWeight_MM(bool geterr, std::vector<snu::KMuon> k_muons, TString ID, TString variable, TString tag, int FakeType=0);

  float Get_DataDrivenWeight_EM(bool geterr, std::vector<snu::KElectron> k_electrons,  std::vector<snu::KMuon> k_muons,  TString IDe, TString IDm, TString variable, TString tage,TString tagm, int FakeType=0);


  float Get_DataDrivenWeight_LL(bool geterr, std::vector<snu::KElectron> k_electrons, TString ID, float var1, float var2,TString tag, int FakeType=0);
  float Get_DataDrivenWeight_LL(bool geterr, std::vector<snu::KMuon> k_muons, TString ID, float var1, float var2, TString tag, int FakeType=0);

  float Get_DataDrivenWeight_LL(bool geterr,std::vector<snu::KElectron> k_electrons, std::vector<snu::KMuon> k_muons,  TString IDe, TString IDm, float var1, float var2, TString tage,TString tagm, int FakeType=0);
			


  float GetFakeRateEl(float pt, float eta, TString cut);

  /// ===  trilepton
  float Get_DataDrivenWeight_MMM(bool geterr, std::vector<snu::KMuon> k_muons);

  TString GetElFRKey(TString elidloose, TString elidtight, TString method);
  std::vector<TString> GetElFRKey( TString IDloose,TString IDtight, TString method, std::vector<TString> regs1);
  //==== General lepton fakes
  float Get_DataDrivenWeight(bool geterr, std::vector<snu::KMuon> k_muons, TString muid, int n_muons, std::vector<snu::KElectron> k_electrons, TString elid, int n_electrons, TString elidloose="ELECTRON16_POG_FAKELOOSE", TString method="dijet_ajet40", int HalfSampleErrorDir=0);

  /// = single lepton
  float Get_DataDrivenWeight_E(bool geterr, std::vector<snu::KElectron> k_electrons ,TString ID, TString method);
  float Get_DataDrivenWeight_M(bool geterr, std::vector<snu::KMuon> k_muons,  TString ID, TString method);
  

  void Test();

  // Class object to get event weights for Fakes

  //==== Use pt-cone-corrected
  void SetUsePtCone(bool b);
  double MuonConePt(snu::KMuon muon, double tightiso);

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

  float _electron_isocone;
  float _muon_isocone;
  HNCommonLeptonFakes* m_fakeobj;
  HNCommonLeptonFakesTriLep* m_fakeobjtrilep;

  EventBase* dd_eventbase;
  bool UsePtCone;



};
#endif
