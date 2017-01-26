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

#define NFIT 9

#include "iostream"

using namespace std;

class Kinematic_Fitter 
{
 public:
  Kinematic_Fitter(Bool_t a_chk_debug=kFALSE);
  ~Kinematic_Fitter();
  
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

  //first index: 0 for b, 1 for c, 2 for light quark
  //second index: four jet
  Double_t ts_corr_value[3][4];
  Double_t ts_corr_error[3][4];
  
  vector<TLorentzVector> jet_vector;
  
  static TLorentzVector measured_extra_jet;
  static TLorentzVector measured_met;
  static TLorentzVector measured_lepton;
  TLorentzVector measured_jet[4];
  static TLorentzVector reordered_jet[4];
  
  Int_t n_b_tag;
  Bool_t chk_b_tag[4];
  Bool_t reordered_b_tag[4];
  
  static Double_t error_extra_jet;
  static Double_t error_reordered_jet_pt[4];
  static Double_t error_lepton_pt;
  
  ROOT::Math::Minimizer* minimizer;
  
  static Double_t Chi2_Func(const Double_t* par);
  //void Top_Specific_Correction(const TLorentzVector& jet, const Int_t& jet_type, Double_t par[2]);
  
  ClassDef(Kinematic_Fitter, 1);
};

#endif /*__KINEMATIC_FITTER_H__*/

