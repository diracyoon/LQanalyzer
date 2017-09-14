#ifndef __Fitter_Test_h__
#define __Fitter_Test_h__

#include "AnalyzerCore.h"

#include "My_Defs.h"
#include "Kinematic_Fitter_Base.h"
#include "Kinematic_Fitter_0.h"
#include "Kinematic_Fitter_1.h"
#include "Kinematic_Fitter_2.h"

class Fitter_Test : public AnalyzerCore 
{
 public:
  //constructors                                                                                   
  Fitter_Test();
  ~Fitter_Test();

  //Functions from core
  virtual void BeginCycle() throw(LQError);
  virtual void BeginEvent() throw(LQError);
  virtual void ExecuteEvents() throw(LQError);
  virtual void EndCycle() throw(LQError);
   
  void InitialiseAnalysis() throw(LQError);
  
 private:
  Bool_t chk_debug;

  Kinematic_Fitter_Base* fitter;
  TS_Correction* ts_correction;
  
  ClassDef(Fitter_Test, 1);
};

#endif /* __Fitter_Test_h__ */
