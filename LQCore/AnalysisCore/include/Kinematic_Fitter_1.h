#ifndef __KINEMATIC_FITTER_H__
#define __KINEMATIC_FITTER_H__

#include "Kinematic_Fitter_Base.h"

using namespace std;

class Kinematic_Fitter : public Kinematic_Fitter_Base
{
 public:
  Kinematic_Fitter(const Bool_t& a_chk_debug=kFALSE);
  virtual ~Kinematic_Fitter();
  
  virtual void Fit();
  virtual void Set(const TLorentzVector& a_met, const TLorentzVector& a_lepton, const vector<TLorentzVector>& a_jet_vector, const Bool_t* a_target_jet, const Bool_t* a_b_tag);
  
 protected:
  static Double_t Chi2_Func(const Double_t* par);
  
  ClassDef(Kinematic_Fitter, 1);
};

#endif /*__KINEMATIC_FITTER_H__*/

