#include "Kinematic_Fitter_Base.h"

//////////

ClassImp(Kinematic_Fitter_Base);

//////////

TLorentzVector Kinematic_Fitter_Base::measured_extra_jet;
TLorentzVector Kinematic_Fitter_Base::measured_met;
TLorentzVector Kinematic_Fitter_Base::measured_lepton;
TLorentzVector Kinematic_Fitter_Base::measured_jet[4];
TLorentzVector Kinematic_Fitter_Base::reordered_jet[4];
TLorentzVector Kinematic_Fitter_Base::sum_extra_jet;
TLorentzVector Kinematic_Fitter_Base::measured_ue;

Double_t Kinematic_Fitter_Base::error_extra_jet;
Double_t Kinematic_Fitter_Base::error_reordered_jet_pt[4];
Double_t Kinematic_Fitter_Base::error_lepton_pt;
Double_t Kinematic_Fitter_Base::error_ue;

TLorentzVector Kinematic_Fitter_Base::fitting_jet[4];
TLorentzVector Kinematic_Fitter_Base::fitting_lepton;
TLorentzVector Kinematic_Fitter_Base::fitting_neutrino;

Double_t Kinematic_Fitter_Base::f_chi2_piece[N_CHI2_PIECE];

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
  
  parameter_start[0] = 1.0;
  parameter_start[1] = 1.0;
  parameter_start[2] = 1.0;
  parameter_start[3] = 1.0;
  parameter_start[4] = 1.0;
  parameter_start[5] = 1.0;
  parameter_start[6] = 1.0;
  parameter_start[7] = 1.0;
  
  parameter_step[0] = 0.03;
  parameter_step[1] = 0.03;
  parameter_step[2] = 0.03;
  parameter_step[3] = 0.03;
  parameter_step[4] = 0.01;
  parameter_step[5] = 0.03;
  parameter_step[6] = 0.03;
  parameter_step[7] = 0.03;
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

void Kinematic_Fitter_Base::Get_Chi2_Piece(Double_t chi2_piece_return[11], const TString& type, const Int_t& index)
{
  if(type=="BEST")
    {
      for(Int_t i=0; i<11; i++){ chi2_piece_return[i] = best_chi2_piece[i]; }
      
      return;
    }
  else
    {
      Int_t i = index/24;
      Int_t j = index%24;

      for(Int_t k=0; k<11; k++){ chi2_piece_return[k] = chi2_piece[i][j][k]; }

      return;
    }
    
  return;
}//void Kinematic_Fitter_Base::Get_Chi2_Piece(Double_t chi2_return[11], const TString& type, const Int_t& index)

//////////

TLorentzVector& Kinematic_Fitter_Base::Get_Fitted_Object(const Int_t& obj_index, const TString& type, const Int_t& index)
{
  Int_t i = index/24;  
  Int_t j = index%24;

  //jets
  if(obj_index<4)
    {
      if(type=="BEST") return best_fitted_jet[obj_index];
      else return fitted_jet[i][j][obj_index]; 
    }//jets
    
  //lepton
  else if(obj_index==4)
    {
      if(type=="BEST") return best_fitted_lepton;
      else return fitted_lepton[i][j];
    }

  //neutrino
  else if(obj_index==5)
    {
      if(type=="BEST") return best_fitted_neutrino;
      else return fitted_neutrino[i][j];
    }
  
  TLorentzVector empty;
  return empty; 
}//TLorentzVector Kinematic_Fitter_Base::Get_Fitted_Object(const Int_t& obj_index)

//////////

TLorentzVector& Kinematic_Fitter_Base::Get_Unfitted_Object(const Int_t& obj_index, const TString& type, const Int_t& index)
{
  Int_t i = index/24;
  Int_t j = index%24;
  
  //jets
  if(obj_index<4)
    {
      if(type=="BEST") return best_unfitted_jet[obj_index];
      else return unfitted_jet[j][obj_index];
    }//jet

  //lepton
  else if(obj_index==4)
    {
      if(type=="BEST") return best_unfitted_lepton;
      else return unfitted_lepton;
    }

  //neutrino
  else if(obj_index==5)
    {
      if(type=="BEST") return best_unfitted_neutrino;
      else return unfitted_neutrino[i];
    }
  
  TLorentzVector empty;
  return empty;
}//TLorentzVector& Kinematic_Fitter_Base::Get_Unfitted_Object(const Int_t& obj_index, const TString& type, const Int_t& index)

//////////

Fitter_Result_Container Kinematic_Fitter_Base::Get_Fitter_Result()
{
  Fitter_Result_Container result_container;
  
  //chi2
  result_container.Set_Chi2(best_chi2);
  
  //permutation
  result_container.Set_Permutation(best_permutation);
  
  //parameter
  result_container.Set_Parameters(best_parameter);

  //fitted objects
  for(Int_t i=0; i<4; i++){ result_container.Set_Fitted_Object(i, best_fitted_jet[i]); }
  result_container.Set_Fitted_Object(4, best_fitted_lepton);
  result_container.Set_Fitted_Object(5, best_fitted_neutrino);
    
  return result_container; 
}//Fitter_Result_Container& Kinematic_Fitter_Base::Get_Fitter_Result()

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
  else if(n_b_tag==3)
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
  if(150<t_mass && t_mass<200) chk_native_top_mass = kTRUE;
  
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
      
      //b quark
      if(reordering_index[i]==0 || reordering_index[i]==1)
	{
	  ts_corr = ts_corr_value[reordering_index[i]][2];
          jet_pt_error = ts_corr_error[reordering_index[i]][2];

          q_mass = B_MASS; 
	}
      
      // //c quark in charged higgs case
      // else if(reordering_index[i]==3)
      // 	{
      // 	  ts_corr = ts_corr_value[reordering_index[i]][1];
      // 	  jet_pt_error = ts_corr_error[reordering_index[i]][1];
	  
      // 	  q_mass = C_MASS;
      // 	}
      
      else 
      	{
      	  ts_corr = ts_corr_value[reordering_index[i]][0];                                  
      	  jet_pt_error = ts_corr_error[reordering_index[i]][0];                             
	 
      	  q_mass = 0;                                                                  
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

void Kinematic_Fitter_Base::Resol_Neutrino_Pt()
{
  //cout << "Kinematic_Fitter_Base::Resol_Neutino_Pt" << endl;
    
  //recal MET
  Double_t lepton_mass = measured_lepton.M();
  Double_t cosine = TMath::Cos(measured_met.Phi());
  Double_t sine = TMath::Sin(measured_met.Phi());
  
  Double_t a = measured_lepton.E()*measured_lepton.E() - measured_lepton.Pz()*measured_lepton.Pz() - TMath::Power(measured_lepton.Px()*cosine + measured_lepton.Py()*sine , 2.0);
  Double_t b = (measured_lepton.Px()*cosine + measured_lepton.Py()*sine)*(lepton_mass*lepton_mass - W_MASS*W_MASS);
  Double_t determinant = TMath::Power(lepton_mass*lepton_mass - W_MASS*W_MASS, 2.)*(measured_lepton.E()*measured_lepton.E() - measured_lepton.Pz()*measured_lepton.Pz());
  
  //cout << "a = " << a << endl;
  //cout << "b = " << b << endl;
  //cout << "det = "<< determinant << endl;

  Double_t met_recal[2];
  met_recal[0] = (-b + TMath::Sqrt(determinant))/2./a;
  met_recal[1] = (-b - TMath::Sqrt(determinant))/2./a;
  
  Double_t mass_diff[2];
  TLorentzVector met_recal_vector[2];
  for(Int_t i=0; i<2; i++)
    {
      met_recal_vector[i].SetPxPyPzE(met_recal[i]*cosine, met_recal[i]*sine, 0, met_recal[i]);
      
      TLorentzVector w_lnu;
      w_lnu = measured_lepton + met_recal_vector[i];
      
      Double_t w_lnu_mass = w_lnu.M();
      mass_diff[i] = TMath::Abs(W_MASS - w_lnu_mass);
      
      //cout << measured_met.Pt() << "\t" << met_recal[i] << "\t" << w_lnu.M() << "\t" << mass_diff[i] << endl;
    }
  //cout << endl;
  
  if(mass_diff[0]<mass_diff[1]) measured_met = met_recal_vector[0];
  else measured_met = met_recal_vector[1];
  
  return;
}//void Kinematic_Fitter_Base::Resol_Neutrino_Pt()

//////////

void Kinematic_Fitter_Base::Set_Minimizer_Parameters(const Int_t& i)
{
  parameter_start[7] = neutrino_pz[i];
  
  for(Int_t j=0; j<NFIT; j++){ minimizer->SetVariable(j, parameter_name[j], parameter_start[j], parameter_step[j]); }
    
}//void Kinematic_Fitter_Base::Set_Minimizer_Parameters(const Int_t& i)

//////////

Bool_t Kinematic_Fitter_Base::Sol_Neutrino_Pz()
{
  //cout << "Kinematic_Fitter_Base::Sol_Neutrino_Pz" << endl;
  
  Double_t lepton_mass =  measured_lepton.M();

  Double_t k = W_MASS*W_MASS/2.0 - lepton_mass*lepton_mass/2.0 + measured_lepton.Px()*measured_met.Px() + measured_lepton.Py()*measured_met.Py();
  Double_t a = TMath::Power(measured_lepton.Px(), 2.0) + TMath::Power(measured_lepton.Py(), 2.0);   
  Double_t b = -2*k*measured_lepton.Pz();                                                           
  Double_t c = TMath::Power(measured_lepton.Pt(), 2.0)*TMath::Power(measured_met.Pt(), 2.0) - TMath::Power(k, 2.0);

  Double_t determinant = TMath::Power(b, 2.0) - 4*a*c;
  
  //cout << "determinant = " << determinant << endl;
  
  //real soluion
  if(determinant>=0)
    {
      neutrino_pz[0] = (-b + TMath::Sqrt(determinant))/(2*a);                      
      neutrino_pz[1] = (-b - TMath::Sqrt(determinant))/(2*a);                      
      
      return kTRUE;
    }
  
  //complex solution. Let's take real part and recompute MET.
  else                                                                                              
    {
      neutrino_pz[0] = -b/(2*a);
      
      Resol_Neutrino_Pt();
      
      return kFALSE;
    }      

  return kFALSE;
}//Bool_t Kinematic_Fitter_Base::Sol_Neutrino_Pz()

//////////

void Kinematic_Fitter_Base::Store_Results(const Int_t& i, const Int_t& j)
{
  /*store results*/
  
  //chi2
  chi2[i][j] = minimizer->MinValue();
  
  for(Int_t k=0; k<10; k++){ chi2_piece[i][j][k] = f_chi2_piece[k]; }
  
  //parameters 
  const Double_t* parameter_result = minimizer->X();
  for(Int_t k=0; k<NFIT; k++){ parameter[i][j][k] = parameter_result[k]; }
      
  //fitted jets
  for(Int_t k=0; k<4; k++){ fitted_jet[i][j][k] = fitting_jet[k]; }
  
  //fitted lepton
  fitted_lepton[i][j] = fitting_lepton;

  //fitted neutrino
  fitted_neutrino[i][j] = fitting_neutrino;

  //update best results
  if(chi2[i][j]<best_chi2) 
    {
      chk_convergence = kTRUE;
      
      //chi
      best_chi2 = chi2[i][j];                                                           
      for(Int_t k=0; k<10; k++){ best_chi2_piece[k] = chi2_piece[i][j][k]; } 

      //parameters
      for(Int_t k=0; k<NFIT; k++){ best_parameter[k] = parameter[i][j][k]; } 
	  
      //permutation & jets
      for(Int_t k=0; k<4; k++)
	{
	  best_permutation[k] = reordering_index[k];                                    
	  best_fitted_jet[k] =  fitted_jet[i][j][k];                                    
	  best_unfitted_jet[k] = unfitted_jet[j][k];
	}    
      
      //lepton
      best_fitted_lepton = fitted_lepton[i][j];
      best_unfitted_lepton = unfitted_lepton;
      
      //fitted neutrino
      best_fitted_neutrino = fitted_neutrino[i][j];
      best_unfitted_neutrino = unfitted_neutrino[i];
    }
  
  return;
}//void Kinematic_Fitter_Base::Store_Results(const Int_t& i, const Int_t& j)

//////////

void Kinematic_Fitter_Base::Store_Unfitted_Object(const Int_t& i, const Int_t& j)
{
  //jets
  for(Int_t k=0; k<4; k++){ unfitted_jet[j][k] = reordered_jet[k]; }
  
  //lepton
  unfitted_lepton = measured_lepton;

  //neutrino
  Double_t px = measured_met.Px();
  Double_t py = measured_met.Py();
  Double_t pz = neutrino_pz[i];
  Double_t energy = TMath::Sqrt(px*px + py*py + pz*pz);
  
  unfitted_neutrino[i].SetPxPyPzE(pz, py, pz, energy);

}//void Kinematic_Fitter_Base::Store_Unfitted_Object(const Int_t& i, const Int_t& j)

//////////
