#ifndef __TS_CORRECTION_H__
#define __TS_CORRECTION_H__

#include <iostream>

#include <TLorentzVector.h>

#define N_QUARK_TYPE 3

using namespace std;

class TS_Correction
{
 public:
  TS_Correction(const Int_t& correction_type);
  ~TS_Correction();

  void Get_Correction(const TLorentzVector& jet, const Int_t jet_type, Double_t corr_val[2]);

 private:
  Int_t n_corr_para;
  Int_t n_eta_bin;
 
  Double_t*** corr_para;
  Double_t*** corr_para_error;

  void Parameters_Reader();
  
  ClassDef(TS_Correction, 1);
};

#endif /* __TS_CORRECTION_H__ */
