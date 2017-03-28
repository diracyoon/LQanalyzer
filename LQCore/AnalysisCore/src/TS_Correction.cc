#include "TS_Correction.h"

//////////

ClassImp(TS_Correction);

//////////

TS_Correction::TS_Correction(const Int_t& a_correction_type)
{
  correction_type = a_correction_type;
  
  cout << "TS_Correction: corrrection_ type " << correction_type << endl; 

  //8 TeV correction
  if(correction_type==0)
    {
      n_corr_para = 4;
      n_corr_para_error = 4;
      n_eta_bin = 9;

      corr_func_ptr = &TS_Correction::Correction_Func_Type_Old;
    }//8 TeV correction
  else
    {
      n_corr_para = 4;
      n_corr_para_error = 3;
      n_eta_bin = 9;

      corr_func_ptr = &TS_Correction::Correction_Func_Type_New;
    }

  Parameters_Reader();
  
}//TS_Correction::TS_Correction()

//////////

TS_Correction::~TS_Correction()
{
  
  for(Int_t i=0; i<N_QUARK_TYPE; i++)
    {
      for(Int_t j=0; j<n_eta_bin; j++)
	{
	  delete[] corr_para[i][j]; 
	  delete[] corr_para_error[i][j];
	}
	    
      delete[] corr_para[i]; 
      delete[] corr_para_error[i];
    }
  
  delete[] corr_para;
  delete[] corr_para_error;

  delete[] corr_eta_bin;
}//TS_Correction::~TS_Correction()

//////////

void TS_Correction::Get_Correction(const TLorentzVector& jet, const Int_t jet_type, Double_t corr_val[2])
{
  Double_t pt = jet.Pt();     
                                                                      
  Double_t eta = TMath::Abs(jet.Eta());
  Int_t eta_bin = Find_Eta_Bin(eta);
  
  (this->*corr_func_ptr)(pt, eta_bin, jet_type, corr_val);

  return;
}//void Get_Correction(const TLorentzVector& jet, const Int_t jet_type, Double_t par[2])

//////////

void TS_Correction::Correction_Func_Type_Old(const Double_t& pt, const Int_t& eta_bin, const Int_t& jet_type, Double_t corr_val[2])
{
  corr_val[0] = corr_para[jet_type][eta_bin][0];
  corr_val[0] += corr_para[jet_type][eta_bin][1]*TMath::Sqrt(pt);
  corr_val[0] += corr_para[jet_type][eta_bin][2]/pt;
  corr_val[0] += corr_para[jet_type][eta_bin][3]*pt;
  corr_val[0] += 1;

  corr_val[1] = corr_para_error[jet_type][eta_bin][0];
  corr_val[1] += corr_para_error[jet_type][eta_bin][1]*TMath::Sqrt(pt);
  corr_val[1] += corr_para_error[jet_type][eta_bin][2]/pt;
  corr_val[1] += corr_para_error[jet_type][eta_bin][3]*pt;
  corr_val[1] *= corr_val[0]*pt;

  return;
}//void TS_Correction::Correction_Func_Type_Old(const Double_t& pt, const Int_t& eta_bin, const Int_t& jet_type, Double_t corr_val[2])

//////////

void TS_Correction::Correction_Func_Type_New(const Double_t& pt, const Int_t& eta_bin, const Int_t& jet_type, Double_t corr_val[2])
{
  corr_val[0] = corr_para[jet_type][eta_bin][0];
  corr_val[0] += corr_para[jet_type][eta_bin][1]*pt;
  corr_val[0] += corr_para[jet_type][eta_bin][2]*TMath::Sqrt(pt);
  corr_val[0] += corr_para[jet_type][eta_bin][3]/pt;
  corr_val[0] += 1;

  corr_val[1] = corr_para_error[jet_type][eta_bin][0];
  corr_val[1] += corr_para_error[jet_type][eta_bin][1]*pt;
  corr_val[1] += corr_para_error[jet_type][eta_bin][2]/pt;
  corr_val[1] = TMath::Exp(corr_val[1]);
  corr_val[1] = corr_val[0]*corr_val[1]*pt;

  return;
}//void TS_Correction::Correction_Func_Type_New(const Double_t& pt, const Int_t& eta_bin, const Int_t& jet_type, Double_t corr_val[2])

//////////

Int_t TS_Correction::Find_Eta_Bin(const Double_t& eta)
{
  Int_t bin = -1;
  
  for(Int_t i=0; i<n_eta_bin; i++)
    {
      if(eta<corr_eta_bin[i])
	{
	  bin = i;
	  break;
	}
    }

  return bin;
}//Int_t TS_Correction::Find_Eta_Bin(const Double_t& eta)

//////////

void TS_Correction::Parameters_Reader()
{
  //construct memory space
  corr_eta_bin = new Double_t[n_eta_bin];
  
  corr_para = new Double_t**[N_QUARK_TYPE];
  corr_para_error = new Double_t**[N_QUARK_TYPE];

  for(Int_t i=0; i<N_QUARK_TYPE; i++)
    {
      corr_para[i] = new Double_t* [n_eta_bin];
      corr_para_error[i] = new Double_t* [n_eta_bin];

      for(Int_t j=0; j<n_eta_bin; j++)
        {
          corr_para[i][j] = new Double_t[n_corr_para];
          corr_para_error[i][j] = new Double_t[n_corr_para_error];
        }//n eta bin
    }//n quark type
  
  //open TS_Correctino.dat
  string lqdir = getenv("LQANALYZER_DIR");
  string path = lqdir + "/CATConfig/TS_Correction_Config/TS_Correction.dat";
  
  ifstream para_dat;
  para_dat.open(path.c_str());
  
  if(!para_dat.is_open())
    {
      cerr << "Can not find TS_Correction.dat!! Abort the process!!"  << endl;
      exit(EXIT_FAILURE);
    }
     
  TString target("Correction Type : ");
  target += correction_type;
  
  string buf;
  while(1)
    {
      getline(para_dat, buf);
      if(para_dat.eof()==kTRUE)
	{
	  cerr << "Can not find target correction type!!" << endl;
	  exit(EXIT_FAILURE);
	}
      
      TString tbuf(buf);
      
      if(tbuf.Contains(target)==kTRUE)
	{
	  getline(para_dat, buf);
	  break;
	} 
    }
  
  //read parameters
  istringstream iss;
  for(Int_t i=0; i<n_eta_bin; i++)
    {
      for(Int_t j=0; j<N_QUARK_TYPE; j++)
	{
	  getline(para_dat, buf);
	  if(para_dat.eof() || buf.compare("")==0)
	    {
	      cerr << "Wrong data format!!"<< endl;
	      exit(EXIT_FAILURE);
	    }
	  
	  iss.clear();
	  iss.str(buf);
	  
	  for(Int_t k=0; k<n_corr_para; k++){ iss >> corr_para[j][i][k]; } //cout << corr_para[j][i][k] << " "; } 
	  for(Int_t k=0; k<n_corr_para_error; k++){ iss >> corr_para_error[j][i][k]; } //cout << corr_para_error[j][i][k] << " "; }
	  //cout << endl;
	}//n quark type
    
      getline(para_dat, buf);
      if(para_dat.eof() || buf.compare("")==0)
	{
          cerr << "Wrong data format!!" << endl;
          exit(EXIT_FAILURE);
        }

      iss.clear();
      iss.str(buf);

      iss >> corr_eta_bin[i];
            
      getline(para_dat, buf);
    }//n eta bin
  
  return;
}//void TS_Correction::Parameters_Reader()

//////////
