#ifndef __TS_CORRECTION_H__
#define __TS_CORRECTION_H__

#include <iostream>
#include <fstream>
#include <sstream>

#include <TLorentzVector.h>
#include <TString.h>

#define N_QUARK_TYPE 3

using namespace std;

class TS_Correction
{
 public:
  TS_Correction(const Int_t& a_correction_type);
  ~TS_Correction();

  void Get_Correction(const TLorentzVector& jet, const Int_t jet_type, Double_t corr_val[2]);

 private:
  Int_t correction_type;

  Int_t n_corr_para;
  Int_t n_corr_para_error;
  Int_t n_eta_bin;
 
  Double_t* corr_eta_bin;
  Double_t*** corr_para;
  Double_t*** corr_para_error;

  void (TS_Correction::*corr_func_ptr)(const Double_t& pt, const Int_t& eta_bin, const Int_t& jet_type, Double_t corr_val[2]);
    
  //correction fuction for 8 TeV ana
  void Correction_Func_Type_Old(const Double_t& pt, const Int_t& eta_bin, const Int_t& jet_type, Double_t corr_val[2]);
  
  //correction fuction for 13 TeV ana
  void Correction_Func_Type_New(const Double_t& pt, const Int_t& eta_bin, const Int_t& jet_type, Double_t corr_val[2]);
  
  Int_t Find_Eta_Bin(const Double_t& eta);
  void Parameters_Reader();
  
  ClassDef(TS_Correction, 1);
};

#endif /* __TS_CORRECTION_H__ */
