#ifndef GetFakeSF_h
#define GetFakeSF_h

#include "AnalyzerCore.h"

class GetFakeSF : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  GetFakeSF();
  ~GetFakeSF();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();



 private:

 

  ClassDef ( GetFakeSF, 1);
};
#endif
