#ifndef __FITTER_RESULT_CONTAINER_H__
#define __FITTER_RESULT_CONTAINER_H__

#include "iostream"

#include "TLorentzVector.h"

#define NFIT 9

using namespace std;

class Fitter_Result_Container
{
 public:
  Fitter_Result_Container();
  ~Fitter_Result_Container();
  
  //
  void Set_Chi2(Double_t a_best_chi2){ best_chi2 = a_best_chi2; }
  void Set_Fitted_Object(const Int_t& obj_index, const TLorentzVector& a_best_fitted_obj);
  void Set_Parameters(const Double_t a_parameters[NFIT]);
  void Set_Permutation(const Int_t a_best_permutation[4]);
  
  //
  Double_t Get_Chi2(){ return best_chi2; }
  TLorentzVector& Get_Fitted_Object(const Int_t& obj_index);
  void Get_Parameters(Double_t parameters_return[NFIT]); 
  void Get_Permutation(Int_t permutation_return[4]);

  Fitter_Result_Container& operator= (const Fitter_Result_Container& a_container);

 private:
  Double_t best_chi2;
  
  Int_t best_permutation[4];
 
  Double_t best_parameter[NFIT]; 
  
  TLorentzVector best_fitted_jet[4];
  TLorentzVector best_fitted_lepton;
  TLorentzVector best_fitted_neutrino;
  
  ClassDef(Fitter_Result_Container, 1);
};

#endif /*__FITTER_RESULT_CONTAINER_H__ */
