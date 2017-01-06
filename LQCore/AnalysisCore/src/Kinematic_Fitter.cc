#include "Kinematic_Fitter.h"

//////////

ClassImp(Kinematic_Fitter);

//////////

TLorentzVector Kinematic_Fitter::measured_met;
TLorentzVector Kinematic_Fitter::measured_lepton;
TLorentzVector Kinematic_Fitter::sum_extra_jet;
TLorentzVector Kinematic_Fitter::measured_ue;
TLorentzVector Kinematic_Fitter::reordered_jet[4];

Double_t Kinematic_Fitter::error_reordered_jet_pt[4];
Double_t Kinematic_Fitter::error_lepton_pt;
Double_t Kinematic_Fitter::error_ue;

//////////

Kinematic_Fitter::Kinematic_Fitter()
{
  const char* min_name = "Minuit";
  const char* algo_name = "";
  
  minimizer = ROOT::Math::Factory::CreateMinimizer(min_name, algo_name);
  
  minimizer->SetMaxFunctionCalls(10000); 
  minimizer->SetMaxIterations(10000);  
  minimizer->SetTolerance(0.001);
  minimizer->SetPrintLevel(1);
}//Kinematic_Fitter::Kinematic_Fitter()

///////////

Kinematic_Fitter::~Kinematic_Fitter()
{
}//Kinematic_Fitter::~Kinematic_Fitter()

//////////

void Kinematic_Fitter::Clear()
{
  return;
}//

//////////

void Kinematic_Fitter::Fit()
{
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
  
  for(Int_t i=0; i<2; i++)
    {
      Int_t reordering_index[4] = {0, 1, 2, 3};

      for(Int_t j=0; j<24; j++)
	{
	  //jet reordering
	  for(Int_t k=0; k<4; k++)
            {
              Double_t pt = measured_jet[reordering_index[k]].Pt();
              Double_t eta = measured_jet[reordering_index[k]].Eta();
              Double_t phi = measured_jet[reordering_index[k]].Phi();

              Double_t q_mass;
              if(k==3) q_mass = C_MASS;
              else q_mass = B_MASS;

              reordered_jet[k].SetPtEtaPhiM(pt, eta, phi, q_mass);
              
	      error_reordered_jet_pt[k] = error_measured_jet_pt[reordering_index[k]];
              
	      reordered_b_tag[k] = chk_b_tag[reordering_index[k]];
            }

	  //b tag configuration check
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
	  
	  //parameters for fit
	  string parameter_name[NFIT] = {"B_Leptonic_Side", "B_Hadronic_Side", "Jet_0_W_CH", "Jet_1_W_CH", "Lepton", "UE_PX", "UE_PY", "Neutrino_Pz", "W_Mass_Or_CH_Mass"};
	  Double_t parameter_start[NFIT] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, neutrino_pz[i], W_MASS};
	  Double_t parameter_step[NFIT] = {0.03, 0.03, 0.03, 0.03, 0.01, 0.1, 0.1, 1.0, 0.1};

	  ROOT::Math::Functor functor(&Chi2_Func, NFIT);
	  minimizer->SetFunction(functor);
	  
	  for(Int_t p=0; p<NFIT; p++){ minimizer->SetVariable(p, parameter_name[p], parameter_start[p], parameter_step[p]); }
	  
	  //dominization
	  minimizer->Minimize();

	  //jet ordering permutation
	  std::next_permutation(reordering_index, reordering_index+4);
	}//four jet permutation loop 4! = 24 combination
    }//neutrino p_z loop
  
  return;
}//void Kinematic_Fitter::Fit

//////////

void Kinematic_Fitter::Set(const TLorentzVector& a_met, const TLorentzVector& a_lepton, const vector<TLorentzVector>& a_jet_vector, const Bool_t a_chk_b_tag[4], const TLorentzVector& a_ue)
{
  measured_met = a_met;
  measured_lepton = a_lepton;
  jet_vector = a_jet_vector;
  
  n_b_tag = 0;
  sum_extra_jet.SetPtEtaPhiE(0, 0, 0, 0);
  Int_t njet = jet_vector.size();
  for(Int_t i=0; i<njet; i++)
    {
      if(i<4)
	{
	  measured_jet[i] = jet_vector.at(i); 
	  chk_b_tag[i] = a_chk_b_tag[i];
	  
	  if(chk_b_tag[i]==kTRUE) n_b_tag++;
	}
      
      sum_extra_jet += jet_vector.at(i);
    }
  
  measured_ue = a_ue;

  return;
}//void Kinematic_Fitter::Set

//////////

Double_t Kinematic_Fitter::Chi2_Func(const Double_t* par)
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
}//Double_t Kinematic_Fitter::Chi2_Func(const Double_t* par)

//////////
