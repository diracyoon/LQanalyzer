#ifndef HNSignalEfficiencies_h
#define HNSignalEfficiencies_h

#include "AnalyzerCore.h"

class HNSignalEfficiencies : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HNSignalEfficiencies();
  ~HNSignalEfficiencies();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillEventCutFlow(TString cut, TString label , float weight);
  void counter(TString cut, float w);

  void RunAnalysis(TString plottag, TString tightelid, TString vetoelid, TString looseelid);

  float WeightCFEvent(std::vector<snu::KElectron> electrons, bool runchargeflip);  
  float IsDiLep(std::vector<snu::KElectron> electrons);

  bool LowMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runcf);
  bool OptMassCheckSignalRegion(std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, float evmet, bool runchargeflip , std::vector<float> cuts, TString opt);
  bool OptMassCheckSignalRegion(std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, int bjetwp, float evmet, bool runchargeflip , std::vector<float> cuts,  TString opt);
  bool MidMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runcf);
  bool HighMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runcf);

  bool CheckSignalRegion( bool isss,  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, TString name, float w);
  bool CheckSignalRegionNN( bool isss,  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, TString name, float w);
  bool CheckSignalRegionTriLep(  std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, TString name, float w);

  void FillTriggerEfficiency(TString cut, float w, TString label,  std::vector<TString> list);

  float GetTightWeight();
  float  GetMediumWeight();
 
  void GetTriggEfficiency();
  void SignalValidation();
  void RunAnalysis();
   

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;
  std::map<TString, float> mapcounter;


  ClassDef ( HNSignalEfficiencies, 1);
};
#endif
