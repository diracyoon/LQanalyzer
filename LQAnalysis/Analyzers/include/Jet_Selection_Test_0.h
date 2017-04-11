#ifndef __Jet_Selection_Test_0_h__
#define __Jet_Selection_Test_0_h__

#include "AnalyzerCore.h"
#include "Kinematic_Fitter_Base.h"
#include "Kinematic_Fitter.h"
#include "Kinematic_Fitter_Old.h"

#define CSV_THRESHOLD_LOOSE 0.460
#define CSV_THRESHOLD_MEDIUM 0.800
#define CSV_THRESHOLD_TIGHT 0.935

#define DISTANCE_MATCH 0.2

#define BLANK -999

class Jet_Selection_Test_0 : public AnalyzerCore 
{
 public:
  //constructors                                                                                   
  Jet_Selection_Test_0();
  ~Jet_Selection_Test_0();

  //Functions from core
  virtual void BeginCycle() throw(LQError);
  virtual void BeginEvent() throw(LQError);
  virtual void ExecuteEvents() throw(LQError);
  virtual void EndCycle() throw(LQError);
   
  void InitialiseAnalysis() throw(LQError);
  
 private:
  Bool_t chk_debug;

  enum TRUTH_TABLE {TOP, A_TOP, BOTTOM, D_0, D_0_0, D_0_1, A_BOTTOM, D_1, D_1_0, D_1_1} Truth_Table;

  Kinematic_Fitter_Base* fitter;

  Double_t Distance(const snu::KTruth& truth, const snu::KJet& jet);
  Bool_t Parton_Jet_Match(const snu::KTruth gen_truth[], const vector<snu::KJet>& jet_vector, Int_t permutation_real[4]);

  ClassDef(Jet_Selection_Test_0, 1);
};

#endif /* __Jet_Selection_Test_0_h__ */