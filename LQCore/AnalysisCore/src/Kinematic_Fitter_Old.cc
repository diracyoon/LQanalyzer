#include "Kinematic_Fitter_Old.h"

//////////

ClassImp(Kinematic_Fitter_Old);

//////////


Kinematic_Fitter_Old::Kinematic_Fitter_Old(Bool_t a_chk_debug) : Kinematic_Fitter_Base(a_chk_debug)
{
  cout << "Kinematic_Fitter type : Old" << endl;
}//Kinematic_Fitter_Old::Kinematic_Fitter_Old()

///////////

Kinematic_Fitter_Old::~Kinematic_Fitter_Old()
{
}//Kinematic_Fitter_Old::~Kinematic_Fitter_Old()

//////////

void Kinematic_Fitter_Old::Fit()
{
  //clear result value storage
  Clear();
  
  //Calculate neutrion p_z solution
  Sol_Neutrino_Pz();
  
  //Top specific correction
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
}//void Kinematic_Fitter_Old::Fit

//////////

void Kinematic_Fitter_Old::Set(const TLorentzVector& a_met, const TLorentzVector& a_lepton, const vector<snu::KJet>& a_jet_vector, const Bool_t a_b_tag[4])
{
  measured_met = a_met;
  
  measured_lepton = a_lepton;
  error_lepton_pt = 0.01*measured_lepton.Pt();

  n_b_tag = 0;
  sum_extra_jet.SetPtEtaPhiE(0, 0, 0, 0);
  Int_t njet = a_jet_vector.size();
  for(Int_t i=0; i<njet; i++)
    {
      if(i<4)
	{
	  measured_jet[i] = a_jet_vector.at(i); 
	  measured_b_tag[i] = a_b_tag[i];
	  
	  if(measured_b_tag[i]==kTRUE) n_b_tag++;
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
}//void Kinematic_Fitter_Old::Set

//////////

Double_t Kinematic_Fitter_Old::Chi2_Func(const Double_t* par)
{
  Double_t chi2 = 0;
  
  /*Contructing fit vectors*/

  //contruct vectors for four fit jets
  TLorentzVector fit_jet[4];
  for(Int_t i=0; i<4; i++)
    {
      Double_t fit_jet_pt = par[i]*reordered_jet[i].Pt();
      Double_t fit_jet_eta = reordered_jet[i].Eta();
      Double_t fit_jet_phi = reordered_jet[i].Phi();
      Double_t fit_jet_mass = reordered_jet[i].M();
     
      fit_jet[i].SetPtEtaPhiM(fit_jet_pt, fit_jet_eta, fit_jet_phi, fit_jet_mass);
    }//loop for four jets

  //contruct a vector for fit lepton
  Double_t fit_lepton_pt = par[4]*measured_lepton.Pt();
  Double_t fit_lepton_eta = measured_lepton.Eta();
  Double_t fit_lepton_phi = measured_lepton.Phi();
  Double_t fit_lepton_mass = measured_lepton.M();
  
  TLorentzVector fit_lepton;
  fit_lepton.SetPtEtaPhiM(fit_lepton_pt, fit_lepton_eta, fit_lepton_phi, fit_lepton_mass);
  
  //contruct a vector for fit ue
  Double_t fit_ue_px = par[5]*measured_ue.Px();
  Double_t fit_ue_py = par[6]*measured_ue.Py();
  Double_t fit_ue_et = TMath::Sqrt(TMath::Power(fit_ue_px, 2.0)+TMath::Power(fit_ue_py, 2.0));
  
  TLorentzVector fit_ue;
  fit_ue.SetPxPyPzE(fit_ue_px, fit_ue_py, 0, fit_ue_et);
  
  //contruct a vector for fit neutrino
  Double_t fit_neutrino_px = 0;
  fit_neutrino_px -= fit_lepton.Px();
  for(Int_t i=0; i<4; i++){ fit_neutrino_px -= fit_jet[i].Px(); }
  fit_neutrino_px -= fit_ue.Px();
  fit_neutrino_px -= sum_extra_jet.Px();

  Double_t fit_neutrino_py = 0;
  fit_neutrino_py -= fit_lepton.Py();
  for(Int_t i=0; i<4; i++){ fit_neutrino_py -= fit_jet[i].Py(); }
  fit_neutrino_py -= fit_ue.Py();
  fit_neutrino_py -= sum_extra_jet.Py();
  
  Double_t fit_neutrino_pz = par[7];

  Double_t fit_neutrino_e = TMath::Sqrt(TMath::Power(fit_neutrino_px, 2.0) + TMath::Power(fit_neutrino_py, 2.0) + TMath::Power(fit_neutrino_pz, 2.0));

  TLorentzVector fit_neutrino;
  fit_neutrino.SetPxPyPzE(fit_neutrino_px, fit_neutrino_py, fit_neutrino_pz, fit_neutrino_e);

  /*Add chi2*/

  //add chi2 of four jets
  for(Int_t i=0; i<4; i++){ chi2 += TMath::Power(fit_jet[i].Pt()-reordered_jet[i].Pt(), 2.0)/TMath::Power(error_reordered_jet_pt[i], 2.0); }
  
  //add chi2 of lepton
  chi2 += TMath::Power(fit_lepton.Pt()-measured_lepton.Pt(), 2.0)/TMath::Power(error_lepton_pt, 2.0);
  
  //add chi2 of unclustered energy
  chi2 += TMath::Power(fit_ue.Px()-measured_ue.Px(), 2.0)/TMath::Power(error_ue, 2.0);
  chi2 += TMath::Power(fit_ue.Py()-measured_ue.Py(), 2.0)/TMath::Power(error_ue, 2.0);

  //add chi2 of w in leptonic side
  TLorentzVector fit_w_lnu = fit_lepton + fit_neutrino;
  chi2 += TMath::Power(fit_w_lnu.M()-W_MASS, 2.0)/TMath::Power(W_WIDTH, 2.0);
  
  //add chi2 of t in leptonic side
  TLorentzVector fit_t_lnuj = fit_lepton + fit_neutrino + fit_jet[0];
  chi2 += TMath::Power(fit_t_lnuj.M()-T_MASS, 2.0)/TMath::Power(T_WIDTH, 2.0);

  //add chi2 of w or charged higgs in hadronic side
  TLorentzVector fit_w_jj = fit_jet[2] + fit_jet[3];
  chi2 += TMath::Power(fit_w_jj.M()-par[8], 2.0)/TMath::Power(W_WIDTH, 2.0);

  //add chi2 of t in hadronic side
  TLorentzVector fit_t_jjj = fit_jet[1] + fit_jet[2] + fit_jet[3];
  chi2 += TMath::Power(fit_t_jjj.M()-T_MASS, 2.0)/TMath::Power(T_WIDTH, 2.0);
    
  return chi2;
}//Double_t Kinematic_Fitter_Old::Chi2_Func(const Double_t* par)

//////////
