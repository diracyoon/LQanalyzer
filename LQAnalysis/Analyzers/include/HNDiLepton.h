#ifndef HNDiLepton_h
#define HNDiLepton_h

#include "AnalyzerCore.h"

class HNDiLepton : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HNDiLepton();
  ~HNDiLepton();

  enum FUNC {
    VALIDATION=0,
    ANALYSIS=1, 
    OPT=2,
    CUTFLOW=3
  };

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void DoCutFlow(float w);

  void RunMM(TString label, std::vector<snu::KMuon> muons, std::vector<snu::KMuon> muons_veto,std::vector<snu::KElectron> el, std::vector<snu::KJet> alljets, std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets, std::vector<snu::KJet> tjets,float mm_weight ,std::vector<TString> mm_trig, float pt1, float pt2);
			 

  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillEventCutFlow(int cf,TString cut,  float weight, TString label);
  void FillEventCutFlow(TString cut, TString label , float weight);

  float WeightCFEvent(std::vector<snu::KElectron> electrons, bool runchargeflip);  
  float IsDiLep(std::vector<snu::KElectron> electrons);
  bool LowMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runcf);
  bool OptMassCheckSignalRegion(std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, float evmet, bool runchargeflip , std::vector<float> cuts, TString opt);
  bool OptMassCheckSignalRegion(std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, int bjetwp, float evmet, bool runchargeflip , std::vector<float> cuts,  TString opt);
  bool MidMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runcf);
  bool HighMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runcf);

  bool CheckSignalRegion( bool isss,  std::vector<snu::KMuon> muons,  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, std::vector<snu::KJet> alljets,  TString name, float w);


  void FillTriggerEfficiency(TString cut, float w, TString label,  std::vector<TString> list);
  void CheckJetsCloseToLeptons(std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets);

  float GetTightWeight();
  float  GetMediumWeight();
  void GetSSSignalEfficiency(float w);
  void GetOSSignalEfficiency(float w);
  void GetTriggEfficiency();
  void SignalValidation();
  void RunAnalysis();
  void OptimiseID(bool isss);

  float MMWeight(std::vector<snu::KMuon> muons,TString id);
  float EEWeight( std::vector<snu::KElectron> electrons, TString id);
  float EMWeight( std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons,TString elid, TString muid);

  void counter(TString cut, float w);

 private:

  FUNC functionality ;
  bool _ee_channel;
  bool _mm_channel;
  bool _m_channel;
  bool _e_channel;
  bool _em_channel;

  TString _m_tightid;
  TString _m_looseid;
  TString _e_tightid;
  TString _e_looseid;

  float _mm_mll_presel_cut;
  float _ee_mll_presel_cut;
  float _em_mll_presel_cut;

  float _mm_met_presel_cut;
  float _ee_met_presel_cut;
  float _em_met_presel_cut; 


  std::map<TString, float> mapcounter;
  map<TString, pair<pair<TString,TString>, pair<float,TString> >  > fake_hists;
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;
 

  ClassDef ( HNDiLepton, 1);
};
#endif
