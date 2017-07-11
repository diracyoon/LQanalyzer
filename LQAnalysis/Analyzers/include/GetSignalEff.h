#ifndef GetSignalEff_h
#define GetSignalEff_h

#include "AnalyzerCore.h"

class GetSignalEff : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  GetSignalEff();
  ~GetSignalEff();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillEventCutFlow(int cf,TString cut,  float weight);


  bool CheckSignalRegion( bool isss,  std::vector<snu::KMuon> muons, std::vector<snu::KElectron> el,  std::vector<snu::KJet> jets,  std::vector<snu::KJet> alljets, TString name, float w);
				      
  
 private:

 

  ClassDef ( GetSignalEff, 1);
};
#endif
