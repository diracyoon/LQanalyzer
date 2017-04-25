#include "Fitter_Result_Container.h"

//////////

ClassImp(Fitter_Result_Container);

//////////

Fitter_Result_Container::Fitter_Result_Container()
{
  best_chi2 = 999;
}//Fitter_Result_Container::Fitter_Result_Container()

//////////

Fitter_Result_Container::~Fitter_Result_Container()
{
}//Fitter_Result_Container::~Fitter_Result_Container()

//////////

void Fitter_Result_Container::Set_Fitted_Object(const Int_t& obj_index, const TLorentzVector& a_best_fitted_obj)
{
  if(obj_index<4) best_fitted_jet[obj_index] = a_best_fitted_obj;
  else if(obj_index==4) best_fitted_lepton = a_best_fitted_obj;
  else if(obj_index==5) best_fitted_neutrino = a_best_fitted_obj;
  
  return;
}//void Fitter_Result_Container::Set_Fitted_Object(const Int_t& obj_index, const TLorentzVector& a_best_fitted_obj)

//////////

void Fitter_Result_Container::Set_Parameters(const Double_t a_parameters[NFIT])
{
  for(Int_t i=0; i<NFIT; i++){ best_parameter[i] = a_parameters[i]; }
  
  return;
}//void Fitter_Result_Container::Set_Parameters(const Double_t a_parameters[NFIT])

//////////

void Fitter_Result_Container::Set_Permutation(const Int_t a_best_permutation[4])
{
  for(Int_t i=0; i<4; i++){ best_permutation[i] = a_best_permutation[i]; } 
   
  return;  
}//void Fitter_Result_Container::Set_Permutation(const Int_t a_best_permutation[])

//////////

TLorentzVector& Fitter_Result_Container::Get_Fitted_Object(const Int_t& obj_index)
{
  if(obj_index<4) return best_fitted_jet[obj_index];
  else if(obj_index==4) return best_fitted_lepton;
  else if(obj_index==5) return best_fitted_neutrino;

  TLorentzVector empty;
  return empty;
}//TLorentzVector& Fitter_Result_Container::Get_Fitted_Object(const Int_t&obj_index)

//////////

void Fitter_Result_Container::Get_Parameters(Double_t parameter_return[NFIT])
{
  for(Int_t i=0; i<NFIT; i++){ parameter_return[i] = best_parameter[i]; }
  return;
}//void Fitter_Result_Container::Get_Parameters(Double_t parameter_return[NFIT])

//////////

void Fitter_Result_Container::Get_Permutation(Int_t permutation_return[4])
{
  for(Int_t i=0; i<4; i++){ permutation_return[i] = best_permutation[i]; }

  return;
}//TLorentzVector& Fitter_Result_Container::Get_Permutation(Int_t* permutation_return)

//////////

Fitter_Result_Container& Fitter_Result_Container::operator= (const Fitter_Result_Container& a_container)
{
  //chi2
  best_chi2 = a_container.best_chi2;
  
  //permutation
  for(Int_t i=0; i<4; i++){ best_permutation[i] = a_container.best_permutation[i]; }

  //parameters
  for(Int_t i=0; i<NFIT; i++){ best_parameter[i] = a_container.best_parameter[i]; }
  
  for(Int_t i=0; i<4; i++){ best_fitted_jet[i] = a_container.best_fitted_jet[i]; }
  best_fitted_lepton = a_container.best_fitted_lepton;
  best_fitted_neutrino = a_container.best_fitted_neutrino;
  
  return *this;
}//Fitter_Result_Container& Fitter_Result_Container::operator= (const Fitter_Result_Container& a_container)

//////////
