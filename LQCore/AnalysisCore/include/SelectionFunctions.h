#ifndef SelectionFunctions_h
#define SelectionFunctions_h

#include <vector>
#include <string>

class TString;

static const double Mass_Z = 91.1876;
static const double Mass_W = 80.398;

class SelectionFunctions {

 public:
  
  //// constructors
  SelectionFunctions();
  ~SelectionFunctions();
  
  bool TriggerSelector (std::vector<TString> triggernames, std::vector<std::string> inputtriggers, std::vector<bool> triggerdecision, std::vector<int> HLTPrescales, int &prescaler);

  
  double EnergyResScaleFactor (float eta); 
  double EnergyResScaleError  (float eta);
  
};
#endif

