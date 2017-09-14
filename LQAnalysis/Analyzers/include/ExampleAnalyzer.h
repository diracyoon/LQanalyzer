#ifndef ExampleAnalyzer_h
#define ExampleAnalyzer_h

#include "AnalyzerCore.h"
class ExampleAnalyzer : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  ExampleAnalyzer();
  ~ExampleAnalyzer();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( ExampleAnalyzer, 1);
};
#endif
