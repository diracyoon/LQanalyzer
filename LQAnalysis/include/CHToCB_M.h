#ifndef __CHToCB_M_h__
#define __CHToCB_M_h__

#include "AnalyzerCore.h"
#include "Kinematic_Fitter.h"

#define CSV_THRESHOLD 0.679

class CHToCB_M : public AnalyzerCore 
{
 public:
  //constructors                                                                                   
  CHToCB_M();
  ~CHToCB_M();

  //Functions from core
  virtual void BeginCycle() throw(LQError);
  virtual void BeginEvent() throw(LQError);
  virtual void ExecuteEvents() throw(LQError);
  virtual void EndCycle() throw(LQError);
   
  void InitialiseAnalysis() throw(LQError);
  
 private:
  Kinematic_Fitter* fitter;
  
  
  ClassDef(CHToCB_M, 1);
};

#endif /*__CHToCB_M_h__*/
