#include "Kinematic_Fitter_2.h"

//////////

ClassImp(Kinematic_Fitter_2);

//////////

Kinematic_Fitter_2::Kinematic_Fitter_2(const Bool_t& a_chk_debug) : Kinematic_Fitter_Base(a_chk_debug)
{
  cout << "Kinematic_Fitter type : 2" << endl;
}//Kinematic_Fitter_2::Kinematic_Fitter_2()

///////////

Kinematic_Fitter_2::~Kinematic_Fitter_2()
{
}//Kinematic_Fitter_2::~Kinematic_Fitter_2()

//////////

void Kinematic_Fitter_2::Set(const TLorentzVector& a_met, const TLorentzVector& a_lepton, const vector<TLorentzVector>& a_jet_vector, const Bool_t* a_target_jet, const Bool_t* a_b_tag)
{
  measured_met = a_met;
  
  measured_lepton = a_lepton;
  error_lepton_pt = 0.01*measured_lepton.Pt();

  jet_vector = a_jet_vector;
  Int_t n_jet = jet_vector.size();

  n_b_tag = 0;
  Int_t target_index = 0;
  Double_t pt_scalar_sum = 0;
  for(Int_t i=0; i<n_jet; i++)
    {
      
      if(a_target_jet[i]==kTRUE)
	{
	  measured_jet[target_index] = jet_vector.at(i); 
	  measured_b_tag[target_index] = a_b_tag[i];
      
	  if(measured_b_tag[target_index]==kTRUE) n_b_tag++;
	
	  target_index++;
	}//target jet
      else pt_scalar_sum += TMath::Abs(jet_vector.at(i).Pt());
    }
       
  //construct extra jet
  measured_extra_jet.SetPtEtaPhiM(0, 0, 0, 0);
  measured_extra_jet -= measured_met;
  measured_extra_jet -= measured_lepton;
  for(Int_t i=0; i<4; i++){ measured_extra_jet -= measured_jet[i]; }
  
  Double_t value[2];
  ts_correction->Get_Correction(measured_extra_jet, 0, value);
  error_extra_jet = value[1];
  cout << endl;
  cout << "error old = " << error_extra_jet << endl;

  TLorentzVector ts_vector;
  ts_vector.SetPtEtaPhiM(pt_scalar_sum, measured_extra_jet.Eta(), measured_extra_jet.Phi(), measured_extra_jet.M());
  ts_correction->Get_Correction(ts_vector, 0, value);
  error_extra_jet = value[1];
  cout << "pt_scalar_sum = " <<  pt_scalar_sum << ", error new = " << error_extra_jet << endl;

  return;
}//void Kinematic_Fitter_2::S

//////////

void Kinematic_Fitter_2::Fit()
{
  //clear result value storage
  Clear();

  //calculate neutrino P_z solution
  chk_neutrino_pz_real = Sol_Neutrino_Pz();
  
  //take real sol only
  //if(chk_neutrino_pz_real==kFALSE) return;

  //take imaginary sol only
  //if(chk_neutrino_pz_real==kTRUE) return;

  //Apply top specific correction
  Apply_TS_Correction();

  /*Leading four jets & neutrino p_z sol permutation*/
  for(Int_t i=0; i<2; i++)
    {
      //if neutrino_pz is complex, only take real part and skip iteration
      if(chk_neutrino_pz_real==kFALSE && i==1) continue;
      
      for(Int_t j=0; j<24; j++)
	{
	  //permutation on jets
	  Reordering_Jets();

	  //b tag configuration check
	  Bool_t b_tag_config_check = Pass_B_Tag_Configuration();
	  
	  //native top mass check
	  Bool_t native_top_mass_check = Pass_Native_Top_Mass();
	  
	  //if b tag configuration is failed, it is not necessary to proceed minimization step.
	  if(b_tag_config_check==kTRUE && native_top_mass_check==kTRUE)
	    {
	      ROOT::Math::Functor functor(&Chi2_Func, NFIT);	      
	      minimizer->SetFunction(functor);

	      //set parameters for miniminzation
	      Set_Minimizer_Parameters(i);

	      //do minimization
	      minimizer->Minimize();

	      //store results
	      Store_Results(i, j);
	    }//if
	  
	  //jet ordering permutation                                                               
	  std::next_permutation(reordering_index, reordering_index+4);
	}//four jet permutation loop 4! = 24 combination
    }//neutrino p_z loop
  
  return;
}//void Kinematic_Fitter_2::Fit()

//////////

Double_t Kinematic_Fitter_2::Chi2_Func(const Double_t* par)
{
  /*Contructing fit vectors*/

  //contruct vectors for four fit jets
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

  fitting_neutrino.SetPxPyPzE(fitting_neutrino_px, fitting_neutrino_py, fitting_neutrino_pz, fitting_neutrino_e);

  /*Add chi2*/

  //add chi2 of four jets
  for(Int_t i=0; i<4; i++){ f_chi2_piece[i] = TMath::Power(fitting_jet[i].Pt()-reordered_jet[i].Pt(), 2.0)/TMath::Power(error_reordered_jet_pt[i], 2.0); }
    
  //add chi2 of lepton
  f_chi2_piece[4] = TMath::Power(fitting_lepton.Pt()-measured_lepton.Pt(), 2.0)/TMath::Power(error_lepton_pt, 2.0);
  
  //add chi2 of unclustered energy
  f_chi2_piece[5] = TMath::Power(fitting_extra_jet.Px()-measured_extra_jet.Px(), 2.0)/TMath::Power(error_extra_jet, 2.0);
  f_chi2_piece[6] = TMath::Power(fitting_extra_jet.Py()-measured_extra_jet.Py(), 2.0)/TMath::Power(error_extra_jet, 2.0);

  //add chi2 of w in leptonic side
  TLorentzVector fit_w_lnu = fitting_lepton + fitting_neutrino;
  f_chi2_piece[7] = TMath::Power(fit_w_lnu.M()-W_MASS, 2.0)/TMath::Power(W_WIDTH, 2.0);
  
  //add chi2 of t in leptonic side
  TLorentzVector fit_t_lnuj = fitting_lepton + fitting_neutrino + fitting_jet[0];
  f_chi2_piece[8] = TMath::Power(fit_t_lnuj.M()-T_MASS, 2.0)/TMath::Power(T_WIDTH, 2.0);
  
  //add chi2 of t in hadronic side
  TLorentzVector fit_t_jjj = fitting_jet[1] + fitting_jet[2] + fitting_jet[3];
  f_chi2_piece[9] = TMath::Power(fit_t_jjj.M()-T_MASS, 2.0)/TMath::Power(T_WIDTH, 2.0);

  Double_t chi2 = 0;
  for(Int_t i=0; i<N_CHI2_PIECE; i++){ chi2 += f_chi2_piece[i]; }

  return chi2;
}//Double_t Kinematic_Fitter_2::Chi2_Func(const Double_t* par)

//////////
