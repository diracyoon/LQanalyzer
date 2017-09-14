#ifndef HNEMu_h
#define HNEMu_h

#include "AnalyzerCore.h"


class HNEMu : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HNEMu();
  ~HNEMu();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);
  void FillEventCutFlow(TString cut, TString label, float w);
  void FillIsoCutFlow(TString cut, float w);
  void CheckJetsCloseToLeptons(std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets,  TString name);
  

  float WeightCFEvent(std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons,  bool runchargeflip);  
  float IsDiLep(std::vector<snu::KElectron> electrons);

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;
  float m_os_Z_nw;
  float m_ss_Z_nw;
  float m_os_Z;
  float m_ss_Z;


  float k_met, k_emumass, k_emujjmass, k_l1jjmass, k_l2jjmass,  k_st, k_ht;
  float k_emujjmass_lowmass, k_l1jjmass_lowmass, k_l2jjmass_lowmass;

  int k_njet, k_nbjet_l, k_nbjet_m, k_nbjet_t;
  float  k_l1eta, k_l2eta,k_jjmass_lowmass;
  float k_weight, k_l1pt, k_l2pt, k_j1pt, k_jjmass;


  ClassDef ( HNEMu, 1);
};
#endif
