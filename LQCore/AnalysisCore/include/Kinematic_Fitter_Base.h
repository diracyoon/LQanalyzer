#ifndef __KINEMATIC_FITTER_BASE_H__
#define __KINEMATIC_FITTER_BASE_H__

#include "TLorentzVector.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "KJet.h"

#include "TS_Correction.h"

#define NFIT 9

#define C_MASS 1.27
#define B_MASS 4.8
#define T_MASS 172.5
#define T_WIDTH 1.5
#define W_MASS 80.398
#define W_WIDTH 2.141

#define GOODNESS_CUT_SLIDING 0.1
#define CHI2_CUT 5

class Kinematic_Fitter_Base
{
 public:
  Kinematic_Fitter_Base(const Bool_t& a_chk_debug=kFALSE);
  virtual ~Kinematic_Fitter_Base();

  virtual void Fit()=0;
  virtual void Set(const TLorentzVector& a_met, const TLorentzVector& a_lepton, const vector<snu::KJet>& a_jet_vector, const Bool_t a_b_tag[4])=0;
    
  Double_t Get_Chi2(const TString& type="BEST", const Int_t& index=-1);
  Bool_t Get_Convergence(){ return chk_convergence; }
  void Get_Parameters(Double_t parameter_return[NFIT], const TString& type="BEST", const Int_t& index=-1);
  void Get_Permutation(Int_t permuatation_return[4]);
  Double_t Get_Top_Mass();

  Bool_t Pass_Goodness_Cut();
  Bool_t Pass_Goodness_Cut(const Double_t& cut_level);

 protected:
  Bool_t chk_debug;

  Int_t reordering_index[4];
  
  Double_t neutrino_pz[2];

  //first index: four jet
  //second index: 0 for light, 1 for c, 2 for b
  Double_t ts_corr_value[4][3];
  Double_t ts_corr_error[4][3];

  string parameter_name[NFIT];
  Double_t parameter_start[NFIT];
  Double_t parameter_step[NFIT] ;

  Bool_t chk_convergence;
  
  vector<snu::KJet> jet_vector;

  static TLorentzVector measured_extra_jet;
  static TLorentzVector measured_met;
  static TLorentzVector measured_lepton;
  TLorentzVector measured_jet[4];
  static TLorentzVector reordered_jet[4];
  static TLorentzVector sum_extra_jet;
  static TLorentzVector measured_ue;

  static Double_t error_extra_jet;
  static Double_t error_reordered_jet_pt[4];
  static Double_t error_lepton_pt;
  static Double_t error_ue;

  TLorentzVector best_fitted_jet[4];
  TLorentzVector fitted_jet[2][24][4];

  Double_t best_chi2;
  Double_t chi2[2][24];
  
  Int_t best_permutation[4];
  Int_t permutation[2][24][4];
  
  Double_t best_parameter[NFIT];
  Double_t parameter[2][24][NFIT];

  Int_t n_b_tag;
  Bool_t measured_b_tag[4];
  Bool_t reordered_b_tag[4];

  TS_Correction* ts_correction;
  
  ROOT::Math::Minimizer* minimizer;
  
  void Apply_TS_Correction();
  void Clear();
  Bool_t Pass_B_Tag_Configuration();
  Bool_t Pass_Native_Top_Mass();
  void Permutation();
  void Reordering_Jets();
  void Set_Minimizer_Parameters(const Int_t& i);
  Bool_t Sol_Neutrino_Pz();
  void Store_Results(const Int_t& i, const Int_t& j);

  ClassDef(Kinematic_Fitter_Base, 1);
};

#endif /* __KINEMATIC_FITTER_BASE_H__ */
