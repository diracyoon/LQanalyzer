#include "TS_Correction.h"

//////////

ClassImp(TS_Correction);

//////////

TS_Correction::TS_Correction(const Int_t& correction_type)
{
  //8 TeV correction
  if(correction_type==0)
    {
      n_corr_para = 4;
      n_eta_bin = 9;
    }//8 TeV correction

  Parameters_Reader();

  for(Int_t i=0; i<3; i++)
    {
      for(Int_t j=0; j<n_eta_bin; j++)
	{
	  
	}//n eta bin
    }//quark type
}//TS_Correction::TS_Correction()

//////////

TS_Correction::~TS_Correction()
{
  
  for(Int_t i=0; i<N_QUARK_TYPE; i++)
    {
      for(Int_t j=0; j<n_eta_bin; j++){ delete[] corr_para[i][j]; }
	    
      delete[] corr_para[i]; 
    }
  delete[] corr_para;
}//TS_Correction::~TS_Correction()

//////////

void TS_Correction::Get_Correction(const TLorentzVector& jet, const Int_t jet_type, Double_t corr_val[2])
{
  Double_t pt = jet.Pt();                                                                           
  Double_t eta = TMath::Abs(jet.Eta());
  
  //corr_para[jet_type][][0];
  //Double_t mean = par0 + par1_b*TMath::Sqrt(pt)+par2_b/pt+par3_b*pt

}//void Get_Correction(const TLorentzVector& jet, const Int_t jet_type, Double_t par[2])

//////////

void TS_Correction::Parameters_Reader()
{
  //construct memory space
  corr_para = new Double_t**[N_QUARK_TYPE];
  corr_para_error = new Double_t**[N_QUARK_TYPE];

  for(Int_t i=0; i<N_QUARK_TYPE; i++)
    {
      corr_para[i] = new Double_t* [n_eta_bin];
      corr_para_error[i] = new Double_t* [n_eta_bin];

      for(Int_t j=0; j<n_eta_bin; j++)
        {
          corr_para[i][j] = new Double_t[n_corr_para];
          corr_para_error[i][j] = new Double_t[n_corr_para];
        }//n eta bin
    }//n quark type
  
  string lqdir = getenv("LQANALYZER_DIR");
  string path = lqdir + "df";
  for(Int_t i=0; i<3; i++)
    {
      for(Int_t j=0; j<n_eta_bin; j++)
        {

        }//n eta bin                                                                                                                                                                                           
    }//quark type

}//void TS_Correction::Parameters_Reader()

//////////
