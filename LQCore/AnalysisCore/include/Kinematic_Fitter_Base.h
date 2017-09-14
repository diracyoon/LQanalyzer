#ifndef __Kinematic_Fitter_Base_H__
#define __Kinematic_Fitter_Base_H__

#include "TLorentzVector.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "My_Defs.h"
#include "TS_Correction.h"
#include "Fitter_Result_Container.h"

class Kinematic_Fitter_Base
{
 public:
  Kinematic_Fitter_Base(const Bool_t& a_chk_high_mass_fitter=kTRUE, const Int_t& ts_correction_type=4, const Bool_t& a_chk_debug=kFALSE);
  virtual ~Kinematic_Fitter_Base();

  virtual void Fit()=0;
  virtual void Set(const TLorentzVector& a_met, const TLorentzVector& a_lepton, const vector<TLorentzVector>& a_jet_vector, const Bool_t* a_target_jet,  const Bool_t* a_b_tag)=0;
  void Get_B_Tag(Bool_t b_tag_return[4], const TString& type="BEST", const Int_t& index=-1);
  Double_t Get_Chi2(const TString& type="BEST", const Int_t& index=-1);
  void Get_Chi2_Piece(Double_t chi2_piece_return[],  const TString& type="BEST", const Int_t& index=-1);
  Bool_t Get_Convergence_Checker(){ return chk_convergence; }
  Bool_t Get_Neutrino_Pz_Sol_Checker(){ return chk_neutrino_pz_real; }
  TLorentzVector& Get_Fitted_Object(const Int_t& obj_index, const TString& type="BEST", const Int_t& index=-1);
  TLorentzVector& Get_Unfitted_Object(const Int_t& obj_index, const TString& type="BEST", const Int_t& index=-1);
  Fitter_Result_Container Get_Fitter_Result();
  void Get_Parameters(Double_t parameter_return[NFIT], const TString& type="BEST", const Int_t& index=-1);
  void Get_Permutation(Int_t permuatation_return[4], const TString& type="BEST", const Int_t& index=-1);
  Bool_t Pass_Goodness_Cut();
  Bool_t Pass_Goodness_Cut(const Double_t& cut_level);

 protected:
  Bool_t chk_high_mass_fitter;
  Bool_t chk_debug;
  Bool_t chk_neutrino_pz_real; 
  
  Int_t reordering_index[4];
  
  Int_t n_b_tag;
  Bool_t measured_b_tag[4];
  Bool_t reordered_b_tag[4];
  
  Double_t neutrino_pz[2];

  //first index: four jet
  //second index: 0 for light, 1 for c, 2 for b
  Double_t ts_corr_value[4][3];
  Double_t ts_corr_error[4][3];

  string parameter_name[NFIT];
  Double_t parameter_start[NFIT];
  Double_t parameter_step[NFIT] ;

  Bool_t chk_convergence;
  
  vector<TLorentzVector> jet_vector;

  static TLorentzVector measured_extra_jet;
  static TLorentzVector measured_met;
  static TLorentzVector measured_lepton;
  static TLorentzVector measured_jet[4];
  static TLorentzVector reordered_jet[4];
  static TLorentzVector sum_extra_jet;
  static TLorentzVector measured_ue;

  static Double_t error_extra_jet;
  static Double_t error_reordered_jet_pt[4];
  static Double_t error_lepton_pt;
  static Double_t error_ue;

  static TLorentzVector fitting_jet[4];
  static TLorentzVector fitting_lepton;
  static TLorentzVector fitting_neutrino;

  static Double_t f_chi2_piece[N_CHI2_PIECE];
 
  TLorentzVector best_fitted_jet[4];
  TLorentzVector best_unfitted_jet[4];
  TLorentzVector fitted_jet[2][24][4];
  TLorentzVector unfitted_jet[24][4];
  
  TLorentzVector best_fitted_lepton;
  TLorentzVector best_unfitted_lepton;
  TLorentzVector fitted_lepton[2][24];
  TLorentzVector unfitted_lepton;

  TLorentzVector best_fitted_neutrino;
  TLorentzVector best_unfitted_neutrino;
  TLorentzVector fitted_neutrino[2][24];
  TLorentzVector unfitted_neutrino[2];
  
  Double_t best_chi2;
  Double_t chi2[2][24];
  
  Double_t best_chi2_piece[N_CHI2_PIECE];
  Double_t chi2_piece[2][24][N_CHI2_PIECE];
  
  Int_t best_permutation[4];
  Int_t permutation[2][24][4];
  
  Double_t best_parameter[NFIT];
  Double_t parameter[2][24][NFIT];
  
  Bool_t best_b_tag[4];
  Bool_t b_tag[2][24][4];

  TS_Correction* ts_correction;
 
  ROOT::Math::Minimizer* minimizer;
  
  void Apply_TS_Correction();
  void Clear();
  Bool_t Pass_B_Tag_Configuration();
  Bool_t Pass_Native_Top_Mass();
  void Permutation();
  void Reordering_Jets(const Int_t& i, const Int_t& j);
  void Resol_Neutrino_Pt();
  void Set_Minimizer_Parameters(const Int_t& i);
  Bool_t Sol_Neutrino_Pz();
  void Store_Results(const Int_t& i, const Int_t& j, const Bool_t& b_tag_config_check, const Bool_t& native_top_mass_check);
  void Store_Unfitted_Object(const Int_t& i, const Int_t& j);
  
  ClassDef(Kinematic_Fitter_Base, 1);
};

#endif /*__Kinematic_Fitter_Base_H__*/
