#ifndef __Kinematic_Fitter_2_H__
#define __Kinematic_Fitter_2_H__

#include "My_Defs.h"
#include "Kinematic_Fitter_Base.h"

using namespace std;

class Kinematic_Fitter_2 : public Kinematic_Fitter_Base
{
 public:
  Kinematic_Fitter_2(const Bool_t& a_chk_high_mass_fitter=kTRUE, const Bool_t& a_chk_debug=kFALSE);
  virtual ~Kinematic_Fitter_2();
  
  virtual void Fit();
  virtual void Set(const TLorentzVector& a_met, const TLorentzVector& a_lepton, const vector<TLorentzVector>& a_jet_vector, const Bool_t* a_target_jet, const Bool_t* a_b_tag);
  
 protected:
  static Double_t Chi2_Func(const Double_t* par);
  
  ClassDef(Kinematic_Fitter_2, 1);
};

#endif /*__KINEMATIC_FITTER_H__*/

