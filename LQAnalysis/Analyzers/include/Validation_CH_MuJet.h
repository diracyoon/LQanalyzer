#ifndef __Validation_CH_MuJet_h__
#define __Validation_CH_MuJet_h__
                                                                              
#include "EventBase.h"
#include "AnalyzerCore.h"

using namespace std;
using namespace snu;

class Validation_CH_MuJet : public AnalyzerCore
{
 public:
  Validation_CH_MuJet();
  ~Validation_CH_MuJet();

  //functions from core
  virtual void BeginCycle() throw(LQError);
  virtual void BeginEvent() throw(LQError);
  virtual void ExecuteEvents() throw(LQError);
  virtual void EndCycle() throw(LQError);

  void InitialiseAnalysis() throw(LQError);

 protected:

  ClassDef(Validation_CH_MuJet, 1);
};
  
#endif /* __Validation_CH_MuJet_h__ */ 
