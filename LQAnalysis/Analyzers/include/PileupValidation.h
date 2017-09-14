#ifndef PileupValidation_h
#define PileupValidation_h

#include "AnalyzerCore.h"

class PileupValidation : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  PileupValidation();
  ~PileupValidation();

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


  void counter(TString cut, float w);

 private:

  FUNC functionality ;

  std::map<TString, float> mapcounter;
  
 

  ClassDef ( PileupValidation, 1);
};
#endif
