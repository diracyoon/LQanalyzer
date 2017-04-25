#ifndef __Jet_Selection_Test_1_h__
#define __Jet_Selection_Test_1_h__

#include "AnalyzerCore.h"
#include "Kinematic_Fitter_Base.h"
#include "Kinematic_Fitter.h"
#include "Kinematic_Fitter_Old.h"

#define CSV_THRESHOLD_LOOSE 0.460
#define CSV_THRESHOLD_MEDIUM 0.800
#define CSV_THRESHOLD_TIGHT 0.935

class Jet_Selection_Test_1 : public AnalyzerCore 
{
 public:
  //constructors                                                                                   
  Jet_Selection_Test_1();
  ~Jet_Selection_Test_1();

  //Functions from core
  virtual void BeginCycle() throw(LQError);
  virtual void BeginEvent() throw(LQError);
  virtual void ExecuteEvents() throw(LQError);
  virtual void EndCycle() throw(LQError);
   
  void InitialiseAnalysis() throw(LQError);
  
 private:
  Bool_t chk_debug;

  Kinematic_Fitter_Base* fitter;

  ClassDef(Jet_Selection_Test_1, 1);
};

#endif /* __Jet_Selection_Test_1_h__ */
