#ifndef __CHToCB_E_h__
#define __CHToCB_E_h__

#include "AnalyzerCore.h"

#include "Kinematic_Fitter.h"

#define CSV_THRESHOLD 0.679

class CHToCB_E : public AnalyzerCore 
{
 public:
  //constructors                                                                                   
  CHToCB_E();
  ~CHToCB_E();

  //Functions from core
  virtual void BeginCycle() throw(LQError);
  virtual void BeginEvent() throw(LQError);
  virtual void ClearOutputVectors() throw(LQError);
  virtual void ExecuteEvents() throw(LQError);
  virtual void EndCycle() throw(LQError);
  virtual void EndEvent() throw(LQError);

  void InitialiseAnalysis() throw(LQError);
  
 private:
  //The output variables 
  //Vectors for output objetcs
  std::vector<snu::KElectron> out_electrons;
  std::vector<snu::KMuon> out_muons;

  //Kinematic_Fitter kinematic_fitter;
  
  ClassDef(CHToCB_E, 1);
};

#endif /*__CHToCB_E_h__*/
