#ifndef __Jet_Inclusion_Prob_h__
#define __Jet_Inclusion_Prob_h__

#include "AnalyzerCore.h"

#define CSV_THRESHOLD_LOOSE 0.460
#define CSV_THRESHOLD_MEDIUM 0.800
#define CSV_THRESHOLD_TIGHT 0.935

#define DISTANCE_MATCH 0.2

#define BLANK -999

class Jet_Inclusion_Prob : public AnalyzerCore 
{
 public:
  //constructors                                                                                   
  Jet_Inclusion_Prob();
  ~Jet_Inclusion_Prob();

  //Functions from core
  virtual void BeginCycle() throw(LQError);
  virtual void BeginEvent() throw(LQError);
  virtual void ExecuteEvents() throw(LQError);
  virtual void EndCycle() throw(LQError);
   
  void InitialiseAnalysis() throw(LQError);
  
 private:
  Bool_t chk_debug;

  ClassDef(Jet_Inclusion_Prob, 1);
};

#endif /* __Jet_Inclusion_Prob_h__ */
