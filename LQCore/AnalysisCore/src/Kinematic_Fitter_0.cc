#include "Kinematic_Fitter_0.h"

//////////

ClassImp(Kinematic_Fitter_0);

//////////


Kinematic_Fitter_0::Kinematic_Fitter_0(Bool_t a_chk_debug) : Kinematic_Fitter_Base(a_chk_debug)
{
  cout << "Kinematic_Fitter type : 0" << endl;
}//Kinematic_Fitter_0::Kinematic_Fitter_0()

///////////

Kinematic_Fitter_0::~Kinematic_Fitter_0()
{
}//Kinematic_Fitter_0::~Kinematic_Fitter_0()

//////////

void Kinematic_Fitter_0::Fit()
{
  //clear result value storage
  Clear();
  
  //Calculate neutrion p_z solution
  chk_neutrino_pz_real = Sol_Neutrino_Pz();
  
  //Top specific correction
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
}//void Kinematic_Fitter_0::Fit

//////////

void Kinematic_Fitter_0::Set(const TLorentzVector& a_met, const TLorentzVector& a_lepton, const vector<TLorentzVector>& a_jet_vector, const Bool_t* a_target_jet, const Bool_t* a_b_tag)
{
  measured_met = a_met;
  
  measured_lepton = a_lepton;
  error_lepton_pt = 0.01*measured_lepton.Pt();

  Int_t n_jet = a_jet_vector.size();

  n_b_tag = 0;
  sum_extra_jet.SetPtEtaPhiE(0, 0, 0, 0);
  Int_t target_index = 0;
  for(Int_t i=0; i<n_jet; i++)
    {
      if(a_target_jet[i]==kTRUE)
	{
	  measured_jet[target_index] = a_jet_vector.at(i); 
	  measured_b_tag[target_index] = a_b_tag[i];
	  
	  if(measured_b_tag[target_index]==kTRUE) n_b_tag++;
	
	  target_index++;
	}
      
      sum_extra_jet += a_jet_vector.at(i);
    }
  
  //construct unclustered energy
  measured_ue.SetPtEtaPhiE(0, 0, 0, 0);
  measured_ue -= measured_met;
  measured_ue -= measured_lepton;
  for(Int_t i=0; i<4; i++){ measured_ue -= measured_jet[i]; }
  measured_ue -= sum_extra_jet;

  error_ue = 0.5*measured_ue.Et();

  return;
}//void Kinematic_Fitter_0::Set

//////////

Double_t Kinematic_Fitter_0::Chi2_Func(const Double_t* par)
{
  Double_t chi2 = 0;
  
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
  Double_t fitting_ue_px = par[5]*measured_ue.Px();
  Double_t fitting_ue_py = par[6]*measured_ue.Py();
  Double_t fitting_ue_et = TMath::Sqrt(TMath::Power(fitting_ue_px, 2.0)+TMath::Power(fitting_ue_py, 2.0));
  
  TLorentzVector fitting_ue;
  fitting_ue.SetPxPyPzE(fitting_ue_px, fitting_ue_py, 0, fitting_ue_et);
  
  //contruct a vector for fit neutrino
  Double_t fitting_neutrino_px = 0;
  fitting_neutrino_px -= fitting_lepton.Px();
  for(Int_t i=0; i<4; i++){ fitting_neutrino_px -= fitting_jet[i].Px(); }
  fitting_neutrino_px -= fitting_ue.Px();
  fitting_neutrino_px -= sum_extra_jet.Px();

  Double_t fitting_neutrino_py = 0;
  fitting_neutrino_py -= fitting_lepton.Py();
  for(Int_t i=0; i<4; i++){ fitting_neutrino_py -= fitting_jet[i].Py(); }
  fitting_neutrino_py -= fitting_ue.Py();
  fitting_neutrino_py -= sum_extra_jet.Py();
  
  Double_t fitting_neutrino_pz = par[7];

  Double_t fitting_neutrino_e = TMath::Sqrt(TMath::Power(fitting_neutrino_px, 2.0) + TMath::Power(fitting_neutrino_py, 2.0) + TMath::Power(fitting_neutrino_pz, 2.0));

  fitting_neutrino.SetPxPyPzE(fitting_neutrino_px, fitting_neutrino_py, fitting_neutrino_pz, fitting_neutrino_e);

  /*Add chi2*/

  //add chi2 of four jets
  for(Int_t i=0; i<4; i++){ chi2 += TMath::Power(fitting_jet[i].Pt()-reordered_jet[i].Pt(), 2.0)/TMath::Power(error_reordered_jet_pt[i], 2.0); }
  
  //add chi2 of lepton
  chi2 += TMath::Power(fitting_lepton.Pt()-measured_lepton.Pt(), 2.0)/TMath::Power(error_lepton_pt, 2.0);
  
  //add chi2 of unclustered energy
  chi2 += TMath::Power(fitting_ue.Px()-measured_ue.Px(), 2.0)/TMath::Power(error_ue, 2.0);
  chi2 += TMath::Power(fitting_ue.Py()-measured_ue.Py(), 2.0)/TMath::Power(error_ue, 2.0);

  //add chi2 of w in leptonic side
  TLorentzVector fitting_w_lnu = fitting_lepton + fitting_neutrino;
  chi2 += TMath::Power(fitting_w_lnu.M()-W_MASS, 2.0)/TMath::Power(W_WIDTH, 2.0);
  
  //add chi2 of t in leptonic side
  TLorentzVector fitting_t_lnuj = fitting_lepton + fitting_neutrino + fitting_jet[0];
  chi2 += TMath::Power(fitting_t_lnuj.M()-T_MASS, 2.0)/TMath::Power(T_WIDTH, 2.0);

  //add chi2 of w or charged higgs in hadronic side
  TLorentzVector fitting_w_jj = fitting_jet[2] + fitting_jet[3];
  chi2 += TMath::Power(fitting_w_jj.M()-par[8], 2.0)/TMath::Power(W_WIDTH, 2.0);

  //add chi2 of t in hadronic side
  TLorentzVector fitting_t_jjj = fitting_jet[1] + fitting_jet[2] + fitting_jet[3];
  chi2 += TMath::Power(fitting_t_jjj.M()-T_MASS, 2.0)/TMath::Power(T_WIDTH, 2.0);
    
  return chi2;
}//Double_t Kinematic_Fitter_0::Chi2_Func(const Double_t* par)

//////////
