#ifndef __KINEMATIC_FITTER_H__
#define __KINEMATIC_FITTER_H__

#include "vector"

#include "TLorentzVector.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#define C_MASS 1.27
#define B_MASS 4.8
#define T_MASS 172.5
#define T_WIDTH 1.5
#define W_MASS 80.398
#define W_WIDTH 2.141

#define NULL_VALUE 99999
#define NFIT 9

#include "iostream"

using namespace std;

class Kinematic_Fitter
{
 public:
  Kinematic_Fitter();
  ~Kinematic_Fitter();
  
  void Clear();
  void Fit();
  void Set(const TLorentzVector& a_met, const TLorentzVector& a_lepton, const vector<TLorentzVector>& a_jet_vector, const Bool_t a_chk_b_tag[4], const TLorentzVector& a_ue);

 private:
  vector<TLorentzVector> jet_vector;
  
  static TLorentzVector measured_met;
  static TLorentzVector measured_lepton;
  TLorentzVector measured_jet[4];
  static TLorentzVector reordered_jet[4];
  static TLorentzVector sum_extra_jet;
  static TLorentzVector measured_ue;
  
  Int_t n_b_tag;
  Bool_t chk_b_tag[4];
  Bool_t reordered_b_tag[4];
  
  Double_t error_measured_jet_pt[4];
  static Double_t error_reordered_jet_pt[4];
  static Double_t error_lepton_pt;
  static Double_t error_ue;
  
  Double_t best_chi2;
  Double_t best_parameter_error[NFIT];
  Double_t best_parameter_value[NFIT];

  ROOT::Math::Minimizer* minimizer;
  
  static Double_t Chi2_Func(const Double_t* par);
  
  ClassDef(Kinematic_Fitter, 1);
};

#endif /*__KINEMATIC_FITTER_H__*/

