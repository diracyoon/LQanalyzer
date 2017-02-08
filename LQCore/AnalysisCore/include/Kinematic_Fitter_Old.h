#ifndef __KINEMATIC_FITTER_OLD_H__
#define __KINEMATIC_FITTER_OLD_H__

#include "vector"

#include "TLorentzVector.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "TS_Correction.h"

#define C_MASS 1.27
#define B_MASS 4.8
#define T_MASS 172.5
#define T_WIDTH 1.5
#define W_MASS 80.398
#define W_WIDTH 2.141

#define NFIT 9

#include "iostream"

using namespace std;

class Kinematic_Fitter_Old
{
 public:
  Kinematic_Fitter_Old(Bool_t a_chk_debug=kFALSE);
  ~Kinematic_Fitter_Old();
  
  void Clear();
  void Fit();
  Double_t Get_Chi2(const TString& type="BEST", const Int_t& index=-1);
  void Get_Parameters(Double_t parameter_return[NFIT], const TString& type="BEST", const Int_t& index=-1);
  void Get_Permutation(Int_t permuatation_return[4]);
  void Print();
  void Set(const TLorentzVector& a_met, const TLorentzVector& a_lepton, const vector<TLorentzVector>& a_jet_vector, const Bool_t a_chk_b_tag[4]);

 private:
  Bool_t chk_debug;

  Double_t best_chi2;
  Double_t best_parameter[NFIT];
  Double_t best_permutation[4];  

  Double_t chi2[2][24];
  Double_t parameter[2][24][NFIT];

  //first index: 0 for light, 1 for c, 2 for b
  //second index: four jet
  Double_t ts_corr_value[4][3];
  Double_t ts_corr_error[4][3];
  
  static TLorentzVector measured_met;
  static TLorentzVector measured_lepton;
  TLorentzVector measured_jet[4];
  static TLorentzVector reordered_jet[4];
  static TLorentzVector sum_extra_jet;
  static TLorentzVector measured_ue;
  
  Int_t n_b_tag;
  Bool_t chk_b_tag[4];
  Bool_t reordered_b_tag[4];
  
  static Double_t error_reordered_jet_pt[4];
  static Double_t error_lepton_pt;
  static Double_t error_ue;
  
  ROOT::Math::Minimizer* minimizer;
  
  static Double_t Chi2_Func(const Double_t* par);
  
  TS_Correction* ts_correction;

  ClassDef(Kinematic_Fitter_Old, 1);
};

#endif /*__KINEMATIC_FITTER_OLD_H__*/

