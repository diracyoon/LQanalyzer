#include "Kinematic_Fitter.h"

//////////

ClassImp(Kinematic_Fitter);

//////////

TLorentzVector Kinematic_Fitter::measured_extra_jet;
TLorentzVector Kinematic_Fitter::measured_met;
TLorentzVector Kinematic_Fitter::measured_lepton;
TLorentzVector Kinematic_Fitter::reordered_jet[4];

Double_t Kinematic_Fitter::error_extra_jet;
Double_t Kinematic_Fitter::error_reordered_jet_pt[4];
Double_t Kinematic_Fitter::error_lepton_pt;

//////////

Kinematic_Fitter::Kinematic_Fitter(Bool_t a_chk_debug)
{
  cout << "Kinematic_Fitter type : new" << endl;
  
  chk_debug = a_chk_debug;
  
  const char* min_name = "Minuit";
  const char* algo_name = "";
  
  minimizer = ROOT::Math::Factory::CreateMinimizer(min_name, algo_name);
  
  minimizer->SetMaxFunctionCalls(10000); 
  minimizer->SetMaxIterations(10000);  
  minimizer->SetTolerance(0.001);
  minimizer->SetPrintLevel(chk_debug);

  ts_correction = new TS_Correction(1);
}//Kinematic_Fitter::Kinematic_Fitter()

///////////

Kinematic_Fitter::~Kinematic_Fitter()
{
  delete ts_correction;
}//Kinematic_Fitter::~Kinematic_Fitter()

//////////

void Kinematic_Fitter::Clear()
{
  chk_convergence = kFALSE;

  best_chi2 = 999;

  for(Int_t i=0; i<NFIT; i++){ best_parameter[i] = 999; }
  for(Int_t i=0; i<4; i++){ best_permutation[i] = 999; }

  for(Int_t i=0; i<2; i++)
    {
      for(Int_t j=0; j<24; j++)
	{
	  chi2[i][j] = -999;
	  
	  for(Int_t k=0; k<NFIT; k++){ parameter[i][j][k] = 999; }
	}
    }
  
  return;
}//void Kinematic_Fitter::Clear()

//////////

void Kinematic_Fitter::Fit()
{
  //clear result value storage
  Clear();
  
  /*Calculate neutrion p_z solution*/
  
  //neutrino p_z quadratic equation
  Double_t k = W_MASS*W_MASS/2.0 + measured_lepton.Px()*measured_met.Px() + measured_lepton.Py()*measured_met.Py();
  Double_t a = TMath::Power(measured_lepton.Px(), 2.0) + TMath::Power(measured_lepton.Py(), 2.0);
  Double_t b = -2*k*measured_lepton.Pz();
  Double_t c = TMath::Power(measured_lepton.Pt(), 2.0)*TMath::Power(measured_met.Pt(), 2.0) - TMath::Power(k, 2.0);

  Double_t neutrino_pz[2] = {0};
  Double_t determinant = TMath::Power(b, 2.0) - 4*a*c;

  //real soluion
  if(determinant>=0)
    {
      neutrino_pz[0] = (-b + TMath::Sqrt(TMath::Power(b, 2.0) - 4*a*c))/(2*a);
      neutrino_pz[1] = (-b - TMath::Sqrt(TMath::Power(b, 2.0) - 4*a*c))/(2*a);
    }
  
  //complex solution. Let's take real part only.
  else
    {
      neutrino_pz[0] = -b/(2*a);
      neutrino_pz[1] = neutrino_pz[0];
    }
  
  /*Lepton pt error*/
  error_lepton_pt = 0.01*measured_lepton.Pt();

  /*Top specific correction*/
  for(Int_t i=0;i<4; i++)
    {
      for(Int_t j=0; j<3; j++)
	{
	  Double_t value[2];
	  ts_correction->Get_Correction(measured_jet[i], j, value);

	  ts_corr_value[i][j] = value[0];
	  ts_corr_error[i][j] = value[1];
	  
	  //cout << ts_corr_value[i][j] << " " << ts_corr_error[i][j] << endl;
	}//quark type
    }//four jet
  
  /*Unclustered energy error*/
  Double_t value[2];
  ts_correction->Get_Correction(measured_extra_jet, 0, value);
  error_extra_jet = value[1];
  
  /*Leading four jets & neutrino p_z sol permutation*/
  for(Int_t i=0; i<2; i++)
    {
      Int_t reordering_index[4] = {0, 1, 2, 3};

      for(Int_t j=0; j<24; j++)
	{	  
	  //jet reordering
	  for(Int_t k=0; k<4; k++)
            {
	      //apply top specific correction & set quark mass
	      Double_t ts_corr = 1;
	      Double_t jet_pt_error = 1;
	      Double_t q_mass = 1;
              if(reordering_index[k]==3)
		{
		  ts_corr = ts_corr_value[reordering_index[k]][1];
		  jet_pt_error = ts_corr_error[reordering_index[k]][1];
		  
		  q_mass = C_MASS;
		}
              else
		{
		  ts_corr = ts_corr_value[reordering_index[k]][2];
		  jet_pt_error = ts_corr_error[reordering_index[k]][2];
		
		  q_mass = B_MASS;
		}
	      
              Double_t pt = measured_jet[reordering_index[k]].Pt()*ts_corr;
              Double_t eta = measured_jet[reordering_index[k]].Eta();
              Double_t phi = measured_jet[reordering_index[k]].Phi();

	      reordered_jet[k].SetPtEtaPhiM(pt, eta, phi, q_mass);
              
	      //jet pt error 
	      error_reordered_jet_pt[k] = jet_pt_error;
	      
	      //reordering b tag 
	      reordered_b_tag[k] = chk_b_tag[reordering_index[k]];
            }

	  //b tag configuration check
	  Bool_t b_tag_config_check = Pass_B_Tag_Configuration();
	  
	  //native top mass check
	  Bool_t native_top_mass_check = Pass_Native_Top_Mass();

	  //if b tag configuration is failed, it is not necessary to proceed minimization step.
	  if(b_tag_config_check==kTRUE && native_top_mass_check==kTRUE)
	    {
	      //parameters for fit
	      string parameter_name[NFIT] = {"B_Leptonic_Side", "B_Hadronic_Side", "W_CH_Jet_0", "W_CH_Jet_1", "Lepton", "UE_PX", "UE_PY", "Neutrino_Pz", "W_Or_CH_Mass"};
	      Double_t parameter_start[NFIT] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, neutrino_pz[i], W_MASS};
	      Double_t parameter_step[NFIT] = {0.03, 0.03, 0.03, 0.03, 0.01, 0.1, 0.1, 1.0, 0.1};
	      
	      ROOT::Math::Functor functor(&Chi2_Func, NFIT);
	      minimizer->SetFunction(functor);
	      
	      for(Int_t k=0; k<NFIT; k++){ minimizer->SetVariable(k, parameter_name[k], parameter_start[k], parameter_step[k]); }
	      
	      //do minization
	      minimizer->Minimize();
	      
	      //store results
	      chi2[i][j] = minimizer->MinValue();
	      const Double_t* parameter_result = minimizer->X();
	      for(Int_t k=0; k<NFIT; k++){ parameter[i][j][k] = parameter_result[k]; }
	      
	      
	      for(Int_t k=0; k<4; k++)
		{
		  Double_t pt = parameter[i][j][k]*reordered_jet[k].Pt();
		  Double_t eta = reordered_jet[k].Eta();
		  Double_t phi = reordered_jet[k].Phi();
		  Double_t mass = reordered_jet[k].M();
		  
		  fitted_jet[i][j][k].SetPtEtaPhiM(pt, eta, phi, mass);
		}
	      
	      //update best results
	      if(chi2[i][j]<best_chi2)
		{
		  chk_convergence = kTRUE;
		  
		  best_chi2 = chi2[i][j];
		  
		  for(Int_t k=0; k<NFIT; k++){ best_parameter[k] = parameter[i][j][k]; }
		  for(Int_t k=0; k<4; k++)
		    {
		      best_permutation[k] = reordering_index[k];
		      best_fitted_jet[k] =  fitted_jet[i][j][k];
		    }
		}
	    }//if(b_tag_config_check==kTRUE)
	  	  
	  if(chk_debug) cout << "Test " << i << " " << reordering_index[0] << " " << reordering_index[1] << " " << reordering_index[2] << " " << reordering_index[3] << " " << chi2[i][j] << endl << endl;
	  
	  //jet ordering permutation
	  std::next_permutation(reordering_index, reordering_index+4);
	}//four jet permutation loop 4! = 24 combination
    }//neutrino p_z loop
  
  return;
}//void Kinematic_Fitter::Fit

//////////

Double_t Kinematic_Fitter::Get_Chi2(const TString& type, const Int_t& index)
{
  if(type=="BEST") return best_chi2;
  else
    {
      Int_t i = index/24;
      Int_t j = index%24;
      
      return chi2[i][j];
    }
    
  return -1;
}//Double_t Kinematic_Fitter::Get_Chi2(const TString& type, const Int_t& index)

//////////

void Kinematic_Fitter::Get_Parameters(Double_t parameter_return[NFIT], const TString& type, const Int_t& index)
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
}//void Kinematic_Fitter::Get_Parameter(Double_t parameter_return[NFIT], const TString& type, const Int_t& index)

//////////

void Kinematic_Fitter::Get_Permutation(Int_t permutation_return[4])
{
  for(Int_t i=0; i<4; i++){ permutation_return[i] = best_permutation[i]; }

  return;
}//void Kinematic_Fitter::Get_Permutation(Double_t permuatation_return[4])

//////////

Double_t Kinematic_Fitter::Get_Top_Mass()
{
  TLorentzVector t_jjj;
  t_jjj.SetPtEtaPhiM(0, 0, 0, 0);
  
  for(Int_t i=0; i<4; i++)
    {
      if(best_permutation[i]!=0) t_jjj += best_fitted_jet[best_permutation[i]];
    }
 
  Double_t t_mass = t_jjj.M();
  
  return t_mass;
}//Double_t Kinematic_Fitter::Get_Top_Mass()

//////////

Bool_t Kinematic_Fitter::Pass_Chi2_Cut()
{
  return kTRUE;
}//Bool_t Kinematic_Fitter::Pass_Chi2_Cut()

//////////

Bool_t Kinematic_Fitter::Pass_Goodness_Cut()
{
  Bool_t chk_goodness_cut = kTRUE;
  for(Int_t i=0; i<4; i++)
    {
      //sliding cut
      if(GOODNESS_CUT_SLIDING<TMath::Abs(1-best_parameter[i])) chk_goodness_cut = kFALSE;

      //absolute cut
      //if()
    }
  //cout <<  chk_goodness_cut << endl;

  return chk_goodness_cut;
}//Bool_t Kinematic_Fitter::Pass_Goodness_Cut()

//////////

Bool_t Kinematic_Fitter::Pass_Goodness_Cut(const Double_t& cut_level)
{
  Bool_t chk_goodness_cut = kTRUE;
  for(Int_t i=0; i<4; i++)
    {
      //sliding cut
      if(cut_level<TMath::Abs(1-best_parameter[i])) chk_goodness_cut = kFALSE;
      
    }
  
  return chk_goodness_cut;
}//Bool_t Kinematic_Fitter::Pass_Goodness_Cut(const Double_t& cut_level)

//////////

void Kinematic_Fitter::Print()
{
}//void Kinematic_Fitter::Print()

//////////

void Kinematic_Fitter::Set(const TLorentzVector& a_met, const TLorentzVector& a_lepton, const vector<TLorentzVector>& a_jet_vector, const Bool_t a_chk_b_tag[4])
{
  measured_met = a_met;
  measured_lepton = a_lepton;

  jet_vector = a_jet_vector;
  
  n_b_tag = 0;
  for(Int_t i=0; i<4; i++)
    {
      measured_jet[i] = jet_vector.at(i); 
      chk_b_tag[i] = a_chk_b_tag[i];
      
      if(chk_b_tag[i]==kTRUE) n_b_tag++;
    }
       
  //construct unclustered energy
  measured_extra_jet.SetPtEtaPhiM(0, 0, 0, 0);
  measured_extra_jet -= measured_met;
  measured_extra_jet -= measured_lepton;
  for(Int_t i=0; i<4; i++){ measured_extra_jet -= measured_jet[i]; }

  return;
}//void Kinematic_Fitter::Set

//////////

Bool_t Kinematic_Fitter::Pass_B_Tag_Configuration()
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
      //for three b tag event, let's choose jet[3] as b tagged jet. It doesn't destroy                                                                                                                 
      if(reordered_b_tag[0]==kTRUE && reordered_b_tag[1]==kTRUE && reordered_b_tag[2]==kTRUE) b_tag_config_check = kTRUE;
    }
  
  return b_tag_config_check;
}//Bool_t kinematic_Fitter::Pass_B_Tag_Configuration()

//////////

Bool_t Kinematic_Fitter::Pass_Native_Top_Mass()
{
  TLorentzVector t_jjj = reordered_jet[1] + reordered_jet[2] + reordered_jet[3];
  Double_t t_mass = t_jjj.M();
  
  Bool_t chk_native_top_mass = kFALSE;
  if(125<t_mass && t_mass<225) chk_native_top_mass = kTRUE;
  
  return chk_native_top_mass;
}//Bool_t Kinematic_Fitter::Pass_Native_Top_Mass()

//////////

Double_t Kinematic_Fitter::Chi2_Func(const Double_t* par)
{
  Double_t chi2 = 0;
  
  /*Contructing fit vectors*/

  //contruct vectors for four fit jets
  TLorentzVector fitting_jet[4];
  for(Int_t i=0; i<4; i++)
    {
      Double_t fitting_jet_pt = par[i]*reordered_jet[i].Pt();
      Double_t fitting_jet_eta = reordered_jet[i].Eta();
      Double_t fitting_jet_phi = reordered_jet[i].Phi();
      Double_t fitting_jet_mass = reordered_jet[i].M();
     
      fitting_jet[i].SetPtEtaPhiM(fitting_jet_pt, fitting_jet_eta, fitting_jet_phi, fitting_jet_mass);
    }//loop for four jets

  //contruct a vector for fit lepton
  Double_t fitting_lepton_pt = par[4]*measured_lepton.Pt();
  Double_t fitting_lepton_eta = measured_lepton.Eta();
  Double_t fitting_lepton_phi = measured_lepton.Phi();
  Double_t fitting_lepton_mass = measured_lepton.M();
  
  TLorentzVector fitting_lepton;
  fitting_lepton.SetPtEtaPhiM(fitting_lepton_pt, fitting_lepton_eta, fitting_lepton_phi, fitting_lepton_mass);
  
  //contruct a vector for fit ue
  Double_t fitting_extra_jet_px = par[5]*measured_extra_jet.Px();
  Double_t fitting_extra_jet_py = par[6]*measured_extra_jet.Py();
  Double_t fitting_extra_jet_et = TMath::Sqrt(TMath::Power(fitting_extra_jet_px, 2.0)+TMath::Power(fitting_extra_jet_py, 2.0));
  
  TLorentzVector fitting_extra_jet;
  fitting_extra_jet.SetPxPyPzE(fitting_extra_jet_px, fitting_extra_jet_py, 0, fitting_extra_jet_et);
  
  //contruct a vector for fit neutrino
  Double_t fitting_neutrino_px = 0;
  fitting_neutrino_px -= fitting_lepton.Px();
  for(Int_t i=0; i<4; i++){ fitting_neutrino_px -= fitting_jet[i].Px(); }
  fitting_neutrino_px -= fitting_extra_jet.Px();
 
  Double_t fitting_neutrino_py = 0;
  fitting_neutrino_py -= fitting_lepton.Py();
  for(Int_t i=0; i<4; i++){ fitting_neutrino_py -= fitting_jet[i].Py(); }
  fitting_neutrino_py -= fitting_extra_jet.Py();
   
  Double_t fitting_neutrino_pz = par[7];

  Double_t fitting_neutrino_e = TMath::Sqrt(TMath::Power(fitting_neutrino_px, 2.0) + TMath::Power(fitting_neutrino_py, 2.0) + TMath::Power(fitting_neutrino_pz, 2.0));

  TLorentzVector fitting_neutrino;
  fitting_neutrino.SetPxPyPzE(fitting_neutrino_px, fitting_neutrino_py, fitting_neutrino_pz, fitting_neutrino_e);

  /*Add chi2*/

  //add chi2 of four jets
  for(Int_t i=0; i<4; i++){ chi2 += TMath::Power(fitting_jet[i].Pt()-reordered_jet[i].Pt(), 2.0)/TMath::Power(error_reordered_jet_pt[i], 2.0); }
  
  //add chi2 of lepton
  chi2 += TMath::Power(fitting_lepton.Pt()-measured_lepton.Pt(), 2.0)/TMath::Power(error_lepton_pt, 2.0);
  
  //add chi2 of unclustered energy
  chi2 += TMath::Power(fitting_extra_jet.Px()-measured_extra_jet.Px(), 2.0)/TMath::Power(error_extra_jet, 2.0);
  chi2 += TMath::Power(fitting_extra_jet.Py()-measured_extra_jet.Py(), 2.0)/TMath::Power(error_extra_jet, 2.0);

  //add chi2 of w in leptonic side
  TLorentzVector fit_w_lnu = fitting_lepton + fitting_neutrino;
  chi2 += TMath::Power(fit_w_lnu.M()-W_MASS, 2.0)/TMath::Power(W_WIDTH, 2.0);
  
  //add chi2 of t in leptonic side
  TLorentzVector fit_t_lnuj = fitting_lepton + fitting_neutrino + fitting_jet[0];
  chi2 += TMath::Power(fit_t_lnuj.M()-T_MASS, 2.0)/TMath::Power(T_WIDTH, 2.0);

  //add chi2 of w or charged higgs in hadronic side
  TLorentzVector fit_w_jj = fitting_jet[2] + fitting_jet[3];
  chi2 += TMath::Power(fit_w_jj.M()-par[8], 2.0)/TMath::Power(W_WIDTH, 2.0);

  //add chi2 of t in hadronic side
  TLorentzVector fit_t_jjj = fitting_jet[1] + fitting_jet[2] + fitting_jet[3];
  chi2 += TMath::Power(fit_t_jjj.M()-T_MASS, 2.0)/TMath::Power(T_WIDTH, 2.0);
    
  return chi2;
}//Double_t Kinematic_Fitter::Chi2_Func(const Double_t* par)

//////////
