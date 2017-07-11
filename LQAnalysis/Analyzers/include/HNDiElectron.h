#ifndef HNDiElectron_h
#define HNDiElectron_h

#include "AnalyzerCore.h"

class HNDiElectron : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HNDiElectron();
  ~HNDiElectron();

  enum FUNC {
    VALIDATION=0,
    ANALYSIS=1, 
    OPT=2
  };

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillEventCutFlow(int cf,TString cut,  float weight);


  void MakeControlPlots(int method, TString methodtag, float w)throw( LQError );
  void MakeMMControlPlots(int method, TString methodtag, float w)throw( LQError );
  void FillByTriggerTrigger(int iel_trig, TString tag, int method, TString methodtag, double evw);
  void FillByMMTriggerTrigger(int iel_trig, TString tag, int method, TString methodtag, double evw);
  void MakeValidationPlots(float w);

  void RunAnalysis(TString plottag, TString tightelid, TString vetoelid, TString looseelid);

  float WeightCFEvent(std::vector<snu::KElectron> electrons, bool runchargeflip);  
  float IsDiLep(std::vector<snu::KElectron> electrons);
  bool LowMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runcf);
  bool OptMassCheckSignalRegion(std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, float evmet, bool runchargeflip , std::vector<float> cuts, TString opt);
  bool OptMassCheckSignalRegion(std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, int bjetwp, float evmet, bool runchargeflip , std::vector<float> cuts,  TString opt);
  bool MidMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runcf);
  bool HighMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runcf);

  bool CheckSignalRegion( bool isss,  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, std::vector<snu::KJet> alljets,TString name, float w);

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
  

  void counter(TString cut, float w);

 private:

  FUNC functionality ;


  std::map<TString, float> mapcounter;
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;
 

  ClassDef ( HNDiElectron, 1);
};
#endif
