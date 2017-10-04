#ifndef MCDataCorrections_h
#define MCDataCorrections_h

#include <vector>
#include <string>
#include <map>
#include <sstream>  
#include <fstream>      // std::ifstream
#include <iostream> 

class TString;
class Reweight;

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

/// rochester corrections
#include "rochcor2016/RoccoR.h"

// ROOT 
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraphAsymmErrors.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
class MCDataCorrections{

 public:
    
  MCDataCorrections();
  MCDataCorrections(bool isdata);
  ~MCDataCorrections();

  void SetPeriod(int period);
  void SetIsData(bool isdata);
  void PrintSummary();


  TH2F* GetCorrectionHist(TString label);
  bool CheckCorrectionHist(TString label);
  TGraphAsymmErrors* GetCorrectionGraph(TString label);
  bool CheckCorrectionGraph(TString label);
  void FillCorrectionHist(std::string label, std::string dirname, std::string filename, std::string histsname, std::string histtype);
  void FillCorrectionHists();
  void CheckFile(TFile* file);
  TDirectory*   getTemporaryDirectory(void) const;

  void SetupDoubleEGTriggerSF(int ileg, std::string sleg);
  double GetEffDEG1(snu::KElectron el);
  double GetEffDEG2(snu::KElectron el);

  double GetDoubleEGTriggerEff(std::vector<snu::KElectron> el);



  /// Lepton Scale Factors
  double ElectronScaleFactor( TString  elid, std::vector<snu::KElectron> el, int sys=0);
  double ElectronRecoScaleFactor(std::vector<snu::KElectron> el, int sys=0);
  double MuonScaleFactor(TString  muid, std::vector<snu::KMuon> mu, int sys=0);
  double MuonISOScaleFactor(TString muid, std::vector<snu::KMuon> mu,int sys=0);
  double MuonTrackingEffScaleFactor(std::vector<snu::KMuon> mu);  

  double MuonISOScaleFactorPeriodDependant(TString muid, std::vector<snu::KMuon> mu , int cat_period, int sys=0);
  double MuonScaleFactorPeriodDependant(TString  muid, std::vector<snu::KMuon> mu, int cat_period,int sys=0);
  
  /// Trigger 
  double TriggerScaleFactor( std::vector<snu::KElectron> el, std::vector<snu::KMuon> mu, TString trigname, int direction=0);
  double TriggerScaleFactorPeriodDependant( std::vector<snu::KElectron> el, std::vector<snu::KMuon> mu, TString trigname, int catperiod, int direction=0);

  double TriggerEfficiencyLegByLeg(std::vector<snu::KElectron> el, TString elid, std::vector<snu::KMuon> mu, TString muid, int TriggerCategory, int DataOrMC, int direction=0);
  double TriggerEfficiencyLegByLegPeriodDependant(std::vector<snu::KElectron> el, TString elid, std::vector<snu::KMuon> mu, TString muid, int TriggerCategory, int catperiod, int DataOrMC, int direction=0);
  double TriggerEfficiency_DiMuon_passing_DoubleMuonTrigger(snu::KMuon mu1, snu::KMuon mu2, TString leg1, TString leg2, TString muid, int DataOrMC, int catperiod, int direction=0);
  double TriggerEfficiency_DiElectron_passing_DoubleElectronTrigger(snu::KElectron el1, snu::KElectron el2, TString leg1, TString leg2, TString elid, int DataOrMC, int catperiod, int direction=0);
  
  

  /// Other corrections
  float PileupWeightByPeriod(snu::KEvent ev);
  float UserPileupWeight(snu::KEvent ev, int nj = -999);
  float CatPileupWeight(snu::KEvent ev, int syst=0);
  void CorrectMuonMomentum(std::vector<snu::KMuon>& k_muons, std::vector<snu::KTruth> truth);
  float GetCorrectedMuonMomentum(snu::KMuon muon, std::vector<snu::KTruth> truth);
  std::vector<TLorentzVector> MakeTLorentz( std::vector<snu::KElectron> el);
  std::vector<TLorentzVector> MakeTLorentz( std::vector<snu::KMuon> mu);
  std::vector<TLorentzVector> MakeTLorentz( std::vector<snu::KJet> jet);
  std::vector<TLorentzVector> MakeTLorentz( std::vector<snu::KFatJet> jet);

 private:
  
  bool corr_isdata;
  double k_period;

  std::map<TString, TH2F*>  CorrectionMap;
  std::map<TString, TGraphAsymmErrors*>  CorrectionMapGraph;
  std::map<float, std::vector<float>* > deg_etaptmap_leg1;
  std::map<float, std::vector<float>* > deg_etaptmap_leg2;
  Reweight *reweightPU;
  RoccoR *rc;

};
#endif
