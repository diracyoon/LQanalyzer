#include "Kinematic_Fitter.h"

//////////

ClassImp(Kinematic_Fitter);

//////////

Kinematic_Fitter::Kinematic_Fitter(const Bool_t& a_chk_debug) : Kinematic_Fitter_Base(a_chk_debug)
{
  cout << "Kinematic_Fitter type : New" << endl;
}//Kinematic_Fitter::Kinematic_Fitter()

///////////

Kinematic_Fitter::~Kinematic_Fitter()
{
}//Kinematic_Fitter::~Kinematic_Fitter()

//////////

void Kinematic_Fitter::Set(const TLorentzVector& a_met, const TLorentzVector& a_lepton, const vector<snu::KJet>& a_jet_vector, const Bool_t a_b_tag[4])
{
  measured_met = a_met;
  
  measured_lepton = a_lepton;
  error_lepton_pt = 0.01*measured_lepton.Pt();

  jet_vector = a_jet_vector;
  
  n_b_tag = 0;
  for(Int_t i=0; i<4; i++)
    {
      measured_jet[i] = jet_vector.at(i); 
      measured_b_tag[i] = a_b_tag[i];
      
      if(measured_b_tag[i]==kTRUE) n_b_tag++;
    }
       
  //construct unclustered energy
  measured_extra_jet.SetPtEtaPhiM(0, 0, 0, 0);
  measured_extra_jet -= measured_met;
  measured_extra_jet -= measured_lepton;
  for(Int_t i=0; i<4; i++){ measured_extra_jet -= measured_jet[i]; }

  Double_t value[2];
  ts_correction->Get_Correction(measured_extra_jet, 0, value);
  error_extra_jet = value[1];

  return;
}//void Kinematic_Fitter::S

//////////

void Kinematic_Fitter::Fit()
{
  //clear result value storage
  Clear();

  //calculate neutrino P_z solution
  Sol_Neutrino_Pz();

  //Apply top specific correction
  Apply_TS_Correction();

  /*Leading four jets & neutrino p_z sol permutation*/
  for(Int_t i=0; i<2; i++)
    {
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
}//void Kinematic_Fitter::Fit()

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
