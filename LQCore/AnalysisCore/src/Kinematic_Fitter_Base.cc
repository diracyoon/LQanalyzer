#include "Kinematic_Fitter_Base.h"

//////////

ClassImp(Kinematic_Fitter_Base);

//////////

TLorentzVector Kinematic_Fitter_Base::measured_extra_jet;
TLorentzVector Kinematic_Fitter_Base::measured_met;
TLorentzVector Kinematic_Fitter_Base::measured_lepton;
TLorentzVector Kinematic_Fitter_Base::reordered_jet[4];
TLorentzVector Kinematic_Fitter_Base::sum_extra_jet;
TLorentzVector Kinematic_Fitter_Base::measured_ue;

Double_t Kinematic_Fitter_Base::error_extra_jet;
Double_t Kinematic_Fitter_Base::error_reordered_jet_pt[4];
Double_t Kinematic_Fitter_Base::error_lepton_pt;
Double_t Kinematic_Fitter_Base::error_ue;

/////////

Kinematic_Fitter_Base::Kinematic_Fitter_Base(const Bool_t& a_chk_debug)
{  
  chk_debug = a_chk_debug;

  const char* min_name = "Minuit";
  const char* algo_name = "";

  minimizer = ROOT::Math::Factory::CreateMinimizer(min_name, algo_name);

  minimizer->SetMaxFunctionCalls(10000);
  minimizer->SetMaxIterations(10000);
  minimizer->SetTolerance(0.001);
  minimizer->SetPrintLevel(chk_debug);
  
  ts_correction = new TS_Correction(1);

  //set parameters for minimization
  parameter_name[0] = "B_Leptonic_Side";
  parameter_name[1] = "B_Hadronic_Side";
  parameter_name[2] = "W_CH_Jet_0";
  parameter_name[3] = "W_CH_Jet_1";
  parameter_name[4] = "Lepton";
  parameter_name[5] = "UE_PX";
  parameter_name[6] = "UE_PY";
  parameter_name[7] = "Neutrino_Pz";
  parameter_name[8] = "W_Or_CH_Mass";

  parameter_start[0] = 1.0;
  parameter_start[1] = 1.0;
  parameter_start[2] = 1.0;
  parameter_start[3] = 1.0;
  parameter_start[4] = 1.0;
  parameter_start[5] = 1.0;
  parameter_start[6] = 1.0;
  parameter_start[7] = 1.0;
  parameter_start[8] = W_MASS;

  parameter_step[0] = 0.03;
  parameter_step[1] = 0.03;
  parameter_step[2] = 0.03;
  parameter_step[3] = 0.03;
  parameter_step[4] = 0.01;
  parameter_step[5] = 0.1;
  parameter_step[6] = 0.1;
  parameter_step[7] = 0.1;
  parameter_step[8] = 0.1;
}//Kinematic_Fitter_Base::Kinematic_Fitter_Base(const Bool_t& a_chk_debug)

//////////

Kinematic_Fitter_Base::~Kinematic_Fitter_Base()
{
  delete ts_correction;
}//Kinematic_Fitter_Base::~Kinematic_Fitter_Base()

//////////

Double_t Kinematic_Fitter_Base::Get_Chi2(const TString& type, const Int_t& index)
{
  if(type=="BEST") return best_chi2;
  else
    {
      Int_t i = index/24;
      Int_t j = index%24;

      return chi2[i][j];
    }

  return -1;
}//Double_t Kinematic_Fitter_Base::Get_Chi2(const TString& type, const Int_t& index)

//////////

void Kinematic_Fitter_Base::Get_Parameters(Double_t parameter_return[NFIT], const TString& type, const Int_t& index)
{
  if(type=="BEST")
    {
      for(Int_t i=0; i<NFIT; i++){ parameter_return[i] = best_parameter[i]; }
    }
  else
    {
      Int_t i = index/24;
      Int_t j = index%24;

      for(Int_t k=0; k<NFIT; k++){ parameter_return[k] = parameter[i][j][k]; }
    }

  return;
}//void Kinematic_Fitter_Base::Get_Parameters(Double_t parameter_return[NFIT], const TString& type, const Int_t& index)

//////////

void Kinematic_Fitter_Base::Get_Permutation(Int_t permutation_return[4])
{
  for(Int_t i=0; i<4; i++){ permutation_return[i] = best_permutation[i]; }

  return;
}//void Kinematic_Fitter_Base::Get_Permutation(Double_t permuatation_return[4])

//////////

Double_t Kinematic_Fitter_Base::Get_Top_Mass()
{
  TLorentzVector t_jjj;
  t_jjj.SetPtEtaPhiM(0, 0, 0, 0);

  for(Int_t i=0; i<4; i++)
    {
      if(best_permutation[i]!=0) t_jjj += best_fitted_jet[best_permutation[i]];
    }

  Double_t t_mass = t_jjj.M();

  return t_mass;
}//Double_t Kinematic_Fitter_Base::Get_Top_Mass()

//////////

Bool_t Kinematic_Fitter_Base::Pass_Goodness_Cut()
{
  Bool_t chk_goodness_cut = kTRUE;
  for(Int_t i=0; i<4; i++)
    {
      //sliding cut
      if(GOODNESS_CUT_SLIDING<TMath::Abs(1-best_parameter[i])) chk_goodness_cut = kFALSE;
    }
  
  return chk_goodness_cut;
}//Bool_t Kinematic_Fitter_Base::Pass_Goodness_Cut()

//////////

Bool_t Kinematic_Fitter_Base::Pass_Goodness_Cut(const Double_t& cut_level)
{
  Bool_t chk_goodness_cut = kTRUE;
  for(Int_t i=0; i<4; i++)
    {
      //sliding cut
      if(cut_level<TMath::Abs(1-best_parameter[i])) chk_goodness_cut = kFALSE;
    }

  return chk_goodness_cut;
}//Bool_t Kinematic_Fitter_Base::Pass_Goodness_Cut(const Double_t& cut_level)

//////////

void Kinematic_Fitter_Base::Apply_TS_Correction()
{
  for(Int_t i=0;i<4; i++)
    {
      for(Int_t j=0; j<3; j++)                                                                      
	{
	  Double_t value[2];
	  
	  ts_correction->Get_Correction(measured_jet[i], j, value);
	  
	  ts_corr_value[i][j] = value[0];
	  ts_corr_error[i][j] = value[1];

	}//quark type 
    }//four jet
  
  return;
}//void Kinematic_Fitter_Base::Apply_TS_Correction()

//////////

void Kinematic_Fitter_Base::Clear()
{
  chk_convergence = kFALSE;

  best_chi2 = 999;

  for(Int_t i=0; i<NFIT; i++){ best_parameter[i] = 999; }
  for(Int_t i=0; i<4; i++)
    {
      reordering_index[i] = i;
      best_permutation[i] = 999; 
    }

  for(Int_t i=0; i<2; i++)
    {
      for(Int_t j=0; j<24; j++)
	{
	  chi2[i][j] = -999;

	  for(Int_t k=0; k<NFIT; k++){ parameter[i][j][k] = 999; }
        }
    }

  return;

}//void Kinematic_Fitter_Base::Clear()

//////////



//////////

Bool_t Kinematic_Fitter_Base::Pass_B_Tag_Configuration()
{
  Bool_t b_tag_config_check = kFALSE;

  if(n_b_tag==1)
    {
      //[0] jet in leptonic side.
      //[1] jet in leptonic side.
      //[2] jet from hadronic w decay or charged higgs decay. In case of charged higgs, it should be a b jet.
      if(reordered_b_tag[0]==kTRUE || reordered_b_tag[1]==kTRUE) b_tag_config_check = kTRUE;
    }
  else if(n_b_tag==2)
    {
      if(reordered_b_tag[0]==kTRUE && reordered_b_tag[1]==kTRUE) b_tag_config_check = kTRUE;
    }
  else
    {
      //for three b tag event, let's choose jet[3] as b tagged jet. It doesn't destroy generality
      if(reordered_b_tag[0]==kTRUE && reordered_b_tag[1]==kTRUE && reordered_b_tag[2]==kTRUE) b_tag_config_check = kTRUE;
    }
      
  return b_tag_config_check;
}///Bool_t kinematic_Fitter_Base::Pass_B_Tag_Configuration()

//////////

Bool_t Kinematic_Fitter_Base::Pass_Native_Top_Mass()
{
  TLorentzVector t_jjj = reordered_jet[1] + reordered_jet[2] + reordered_jet[3];
  Double_t t_mass = t_jjj.M();

  Bool_t chk_native_top_mass = kFALSE;
  if(125<t_mass && t_mass<225) chk_native_top_mass = kTRUE;

  return chk_native_top_mass;
}//Bool_t Kinematic_Fitter_Base::Pass_Native_Top_Mass()

//////////

void Kinematic_Fitter_Base::Reordering_Jets()
{
  for(Int_t i=0; i<4; i++)
    {
      
      //apply top specific correction & set quark mass
      Double_t ts_corr = 1;
      Double_t jet_pt_error = 1;
      Double_t q_mass = 1;
      
      if(reordering_index[i]==3)
	{
	  ts_corr = ts_corr_value[reordering_index[i]][1];
	  jet_pt_error = ts_corr_error[reordering_index[i]][1];
	  
	  q_mass = C_MASS;
	}
      else 
	{
	  ts_corr = ts_corr_value[reordering_index[i]][2];                                  
	  jet_pt_error = ts_corr_error[reordering_index[i]][2];                             
	  q_mass = B_MASS;                                                                  
	}   
	  
      Double_t pt = measured_jet[reordering_index[i]].Pt()*ts_corr;                         
      Double_t eta = measured_jet[reordering_index[i]].Eta();                               
      Double_t phi = measured_jet[reordering_index[i]].Phi(); 
	
      //reodered_jet
      reordered_jet[i].SetPtEtaPhiM(pt, eta, phi, q_mass);

      //jet pt error
      error_reordered_jet_pt[i] = jet_pt_error;

      //reordering b tag
      reordered_b_tag[i] = measured_b_tag[reordering_index[i]];
    }

  return;
}//void Kinematic_Fitter_Base::Reordering_Jets()

//////////

void Kinematic_Fitter_Base::Set_Minimizer_Parameters(const Int_t& i)
{
  parameter_start[7] = neutrino_pz[i];
  
  for(Int_t j=0; j<NFIT; j++){ minimizer->SetVariable(j, parameter_name[j], parameter_start[j], parameter_step[j]); }
    
}//void Kinematic_Fitter_Base::Set_Minimizer_Parameters(const Int_t& i)

//////////

Bool_t Kinematic_Fitter_Base::Sol_Neutrino_Pz()
{
  Double_t k = W_MASS*W_MASS/2.0 + measured_lepton.Px()*measured_met.Px() + measured_lepton.Py()*measured_met.Py();
  Double_t a = TMath::Power(measured_lepton.Px(), 2.0) + TMath::Power(measured_lepton.Py(), 2.0);   
  Double_t b = -2*k*measured_lepton.Pz();                                                           
  Double_t c = TMath::Power(measured_lepton.Pt(), 2.0)*TMath::Power(measured_met.Pt(), 2.0) - TMath::Power(k, 2.0);

  Double_t determinant = TMath::Power(b, 2.0) - 4*a*c;
  
  //real soluion
  if(determinant>=0)
    {
      neutrino_pz[0] = (-b + TMath::Sqrt(TMath::Power(b, 2.0) - 4*a*c))/(2*a);                      
      neutrino_pz[1] = (-b - TMath::Sqrt(TMath::Power(b, 2.0) - 4*a*c))/(2*a);                      
      
      return kTRUE;
    }
  //complex solution. Let's take real part only.
  else                                                                                              
    {
      neutrino_pz[0] = -b/(2*a);
      neutrino_pz[1] = neutrino_pz[0];
      
      return kFALSE;
    }      

  return kFALSE;
}//Bool_t Kinematic_Fitter_Base::Sol_Neutrino_Pz()

//////////

void Kinematic_Fitter_Base::Store_Results(const Int_t& i, const Int_t& j)
{
  //store results
  chi2[i][j] = minimizer->MinValue();
  
  const Double_t* parameter_result = minimizer->X(); 
  
  for(Int_t k=0; k<NFIT; k++)
    {
      parameter[i][j][k] = parameter_result[k]; 
      
      if(k<4)
	{
	  Double_t pt = parameter[i][j][k]*reordered_jet[k].Pt();                           
	  Double_t eta = reordered_jet[k].Eta();                                            
	  Double_t phi = reordered_jet[k].Phi();                                            
	  Double_t mass = reordered_jet[k].M();                                             
	  
	  fitted_jet[i][j][k].SetPtEtaPhiM(pt, eta, phi, mass);
	}
    }
  
  //update best results
  if(chi2[i][j]<best_chi2) 
    {
      chk_convergence = kTRUE;
      
      best_chi2 = chi2[i][j];                                                           
      for(Int_t k=0; k<NFIT; k++)
	{
	  best_parameter[k] = parameter[i][j][k]; 
	  
	  if(k<4)
	    {
	      best_permutation[k] = reordering_index[k];                                    
	      best_fitted_jet[k] =  fitted_jet[i][j][k];                                    
	    }                                                                               
	}
    }
  
  return;
}//void Kinematic_Fitter_Base::Store_Results(const Int_t& i, const Int_t& j)

//////////
