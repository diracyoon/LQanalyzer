#ifndef __Fitter_Test_h__
#define __Fitter_Test_h__

#include "AnalyzerCore.h"
#include "Kinematic_Fitter_Base.h"
#include "Kinematic_Fitter_0.h"
#include "Kinematic_Fitter_1.h"
#include "Kinematic_Fitter_2.h"

#define CSV_THRESHOLD_LOOSE 0.460
#define CSV_THRESHOLD_MEDIUM 0.800
#define CSV_THRESHOLD_TIGHT 0.935

#define DISTANCE_MATCH 0.2

#define BLANK -999

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

  ClassDef(Fitter_Test, 1);
};

#endif /* __Fitter_Test_h__ */
