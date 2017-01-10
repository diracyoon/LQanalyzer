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
  
  /*Lepton pt error*/
  error_lepton_pt = 0.01*measured_lepton.Pt();

  /*Unclustered energy error*/
  error_ue = 0.5*measured_ue.Et();
  
  /*Top specific correction*/
  Top_Specific_Correction();
  
  /*Leading four jets & neutrino p_z sol permutation*/
  for(Int_t i=0; i<2; i++)
    {
      Int_t reordering_index[4] = {0, 1, 2, 3};

      for(Int_t j=0; j<24; j++)
	{
	  //jet reordering
	  for(Int_t k=0; k<4; k++)
            {
	      //apply top specific correction
	      Double_t ts_corr = 1;
              if(k==3) ts_corr = ts_corr_value[1][reordering_index[k]];
              else ts_corr = ts_corr_value[0][reordering_index[k]];
	      
              Double_t pt = measured_jet[reordering_index[k]].Pt()*ts_corr;
              Double_t eta = measured_jet[reordering_index[k]].Eta();
              Double_t phi = measured_jet[reordering_index[k]].Phi();

	      //set quark mass
              Double_t q_mass;
              if(k==3) q_mass = C_MASS;
              else q_mass = B_MASS;

	      reordered_jet[k].SetPtEtaPhiM(pt, eta, phi, q_mass);
              
	      //jet pt error 
	      Double_t jet_pt_error = 1;
	      if(k==0) jet_pt_error = ts_corr_error[0][reordering_index[k]];
	      else jet_pt_error = ts_corr_error[1][reordering_index[k]];
	      
	      error_reordered_jet_pt[k] = jet_pt_error;
              
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
	  string parameter_name[NFIT] = {"B_Leptonic_Side", "B_Hadronic_Side", "W_CH_Jet_0", "W_CH_Jet_1", "Lepton", "UE_PX", "UE_PY", "Neutrino_Pz", "W_Or_CH_Mass"};
	  Double_t parameter_start[NFIT] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, neutrino_pz[i], W_MASS};
	  Double_t parameter_step[NFIT] = {0.03, 0.03, 0.03, 0.03, 0.01, 0.1, 0.1, 1.0, 0.1};

	  ROOT::Math::Functor functor(&Chi2_Func, NFIT);
	  minimizer->SetFunction(functor);
	  
	  for(Int_t p=0; p<NFIT; p++)
	    { 
	      minimizer->SetVariable(p, parameter_name[p], parameter_start[p], parameter_step[p]);
	      
	      Double_t lower_limit = -1e4;
	      Double_t upper_limit =  1e4;
	      
	      //neutrino pz
	      if(p==7)
		{
		  lower_limit = 0.5*neutrino_pz[i];
		  upper_limit = 1.5*neutrino_pz[i];
		  
		}
	      //W or CH mass 
	      else if(p==8)
		{
		  lower_limit = 0.8*W_MASS;
		  upper_limit = 1.2*150;//heaviest charged higss mass
		}
	      //rescale part
	      else
		{
		  lower_limit = 0;
		  upper_limit = 1.5; 
		}
		
		minimizer->SetVariableLimits(p, lower_limit, upper_limit);
	    }
	  
	  //do minization
	  minimizer->Minimize();

	  //store results
	  const Double_t* parameter_results = minimizer->X();
	  cout << parameter_results[0] << endl;
	  cout << parameter_results[0] << endl;
	  
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

void Kinematic_Fitter::Top_Specific_Correction()
{
  //correction for 8 TeV analysis
  for(Int_t i=0; i<4; i++)
    {
      Double_t pt = measured_jet[i].Pt();
      Double_t eta = TMath::Abs(measured_jet[i].Eta());
      
      Double_t par0_b=0., par1_b=0., par2_b=0., par3_b=0.;
      Double_t par0_c=0., par1_c=0., par2_c=0., par3_c=0.;
      Double_t par0_l=0., par1_l=0., par2_l=0., par3_l=0.;
      Double_t par0_sigb=0, par1_sigb =0, par2_sigb=0, par3_sigb=0;
      Double_t par0_sigc=0, par1_sigc =0, par2_sigc=0, par3_sigc=0;
      Double_t par0_sigl=0, par1_sigl =0, par2_sigl=0, par3_sigl=0;

      if(eta<0.174)
	{
	  par0_b=-0.5774; par1_b=0.0542; par2_b=17.44; par3_b=-0.001465;
	  par0_c=-0.1995; par1_c=0.01237; par2_c=8.528; par3_c=-0.0005324;
	  par0_l=-0.2979; par1_l=0.02754; par2_l=8.201; par3_l=-0.001326;
	  
	  par0_sigb=-0.3333; par1_sigb=0.05044; par2_sigb=14.5; par3_sigb=-0.001823;
	  par0_sigc=0.4101; par1_sigc=-0.04369; par2_sigc=-0.987; par3_sigc=0.001559;
	  par0_sigl=0.2374; par1_sigl=-0.0199; par2_sigl=1.208; par3_sigl=0.0006941;
	}//if(eta<0.174)
      else if(eta<0.348)
	{
	  par0_b=-0.5821; par1_b=0.0546; par2_b=17.57; par3_b=-0.001468;
	  par0_c=-0.2866; par1_c=0.02688; par2_c=9.53; par3_c=-0.001208;
	  par0_l=-0.2291; par1_l=0.01767; par2_l=7.283; par3_l=-0.0009294;
	  
	  par0_sigb=-0.4744; par1_sigb=0.06758; par2_sigb=16.99; par3_sigb=-0.002373;
	  par0_sigc=0.3315; par1_sigc=-0.03358; par2_sigc=-0.01424; par3_sigc=0.001234;
	  par0_sigl=0.3143; par1_sigl=-0.03392; par2_sigl=0.6005; par3_sigl=0.001332;
	}//else if(eta<0.348)
      else if(eta<0.522)
	{
	  par0_b=-0.6125; par1_b=0.05943; par2_b=17.82; par3_b=-0.001656;
	  par0_c=-0.1168; par1_c=0.001053; par2_c=7.585; par3_c=-0.0001227;
	  par0_l=-0.2434; par1_l=0.02058; par2_l=7.29; par3_l=-0.00107;
	  
	  par0_sigb=-0.5308; par1_sigb=0.07682; par2_sigb=17.08; par3_sigb=-0.002757;
	  par0_sigc=0.1435; par1_sigc=-0.007334; par2_sigc=2.413; par3_sigc=0.0002287;
	  par0_sigl=0.4923; par1_sigl=-0.05751; par2_sigl=-1.827; par3_sigl=0.002178;
	}//else if(eta<0.522)
      else if(eta<0.696)
	{
	  par0_b=-0.5643; par1_b=0.0536; par2_b=16.92; par3_b=-0.001468;
	  par0_c=-0.2361; par1_c=0.02047; par2_c=8.749; par3_c=-0.0009981;
	  par0_l=-0.293; par1_l=0.0274; par2_l=7.966; par3_l=-0.001331;
	  
	  par0_sigb=-0.04004; par1_sigb=0.01117; par2_sigb=9.989; par3_sigb=-0.0003874;
	  par0_sigc=0.5479; par1_sigc=-0.06382; par2_sigc=-2.705; par3_sigc=0.002356;
	  par0_sigl=0.2242; par1_sigl=-0.01989; par2_sigl=1.711; par3_sigl=0.000754;
	}//else if(eta<0.696)
      else if(eta<0.879)
	{
	  par0_b=-0.5693; par1_b=0.05391; par2_b=17.08; par3_b=-0.001469;
	  par0_c=-0.09953; par1_c=-0.0002103; par2_c=7.183; par3_c=-0.0001419;
	  par0_l=-0.2977; par1_l=0.02823; par2_l=7.915; par3_l=-0.001372;
	  
	  par0_sigb=-0.5279; par1_sigb=0.07609; par2_sigb=17.52; par3_sigb=-0.002753;
	  par0_sigc=0.372; par1_sigc=-0.03663; par2_sigc=-0.6704; par3_sigc=0.001218;
	  par0_sigl=0.1221; par1_sigl=-0.004342; par2_sigl=2.967; par3_sigl=0.0001104;
	}//else if(eta<0.879)
      else if(eta<1.131)
	{
	  par0_b=-0.5026; par1_b=0.04535; par2_b=15.99; par3_b=-0.001178;
	  par0_c=-0.1355; par1_c=0.004565; par2_c=7.865; par3_c=-0.0003958;
	  par0_l=-0.3262; par1_l=0.03197; par2_l=8.365; par3_l=-0.001546;
	  
	  par0_sigb=-0.008748; par1_sigb=0.006114; par2_sigb=10.23; par3_sigb=-0.0001456;
	  par0_sigc=0.1049; par1_sigc=-0.000284; par2_sigc=3.248; par3_sigc=-0.0001325;
	  par0_sigl=0.2977; par1_sigl=-0.02945; par2_sigl=1.046; par3_sigl=0.001089;
	}//else if(eta<1.131)
      else if(eta<1.392)
	{
	  par0_b=-0.5189; par1_b=0.04724; par2_b=16.3; par3_b=-0.001244;
	  par0_c=-0.153; par1_c=0.007903; par2_c=7.817; par3_c=-0.0005719;
	  par0_l=-0.3484; par1_l=0.03416; par2_l=8.612; par3_l=-0.001646;
	  
	  par0_sigb=-0.00576; par1_sigb=0.006223; par2_sigb=10.61; par3_sigb=-0.0001532;
	  par0_sigc=0.3699; par1_sigc=-0.03524; par2_sigc=-0.06528; par3_sigc=0.001121;
	  par0_sigl=0.3011; par1_sigl=-0.02928; par2_sigl=1.432; par3_sigl=0.001088;
	}//else if(eta<1.392)
      else if(eta<1.74)
	{
	  par0_b=-0.4733; par1_b=0.04073; par2_b=15.53; par3_b=-0.0009859;
	  par0_c=-0.09843; par1_c=-0.002012; par2_c=6.98; par3_c=-0.00005632;
	  par0_l=-0.3413; par1_l=0.03374; par2_l=8.07; par3_l=-0.001653;
	  
	  par0_sigb=-0.1803; par1_sigb=0.03173; par2_sigb=12.78; par3_sigb=-0.001192;
	  par0_sigc=0.1394; par1_sigc=-0.003827; par2_sigc=3.497; par3_sigc=-0.00005752;
	  par0_sigl=0.5663; par1_sigl=-0.06745; par2_sigl=-1.654; par3_sigl=0.00255;
	}//else if(eta<1.74)
      else
	{
	  par0_b=-0.463; par1_b=0.04362; par2_b=13.5; par3_b=-0.001173;
	  par0_c=0.1074; par1_c=-0.03254; par2_c=3.746; par3_c=0.001216;
	  par0_l=-0.1234; par1_l=0.006638; par2_l=4.315; par3_l=-0.0007936;
	  
	  par0_sigb=-0.1383; par1_sigb=0.01703; par2_sigb=13.63; par3_sigb=-0.0004439;
	  par0_sigc=0.4363; par1_sigc=-0.04966; par2_sigc=-0.5373; par3_sigc=0.001785;
	  par0_sigl=0.3632; par1_sigl=-0.04; par2_sigl=0.8006; par3_sigl=0.001426;
	}//else
      
      Double_t meanb = par0_b+par1_b*TMath::Sqrt(pt)+par2_b/pt+par3_b*pt;
      ts_corr_value[0][i]=(1+meanb);

      Double_t meanc = par0_c+par1_c*TMath::Sqrt(pt)+par2_c/pt+par3_c*pt;
      ts_corr_value[1][i]=(1+meanc);

      Double_t meanl = par0_l+par1_l*TMath::Sqrt(pt)+par2_l/pt+par3_l*pt;
      ts_corr_value[2][i]=(1+meanl);
      
      Double_t sigmab = par0_sigb+par1_sigb*TMath::Sqrt(pt)+par2_sigb/pt+par3_sigb*pt;
      ts_corr_error[0][i] = (1+meanb)*sigmab*pt;
      
      Double_t sigmac = par0_sigc+par1_sigc*TMath::Sqrt(pt)+par2_sigc/pt+par3_sigc*pt;
      ts_corr_error[1][i] = (1+meanc)*sigmac*pt;
      
      Double_t sigmal = par0_sigl+par1_sigl*TMath::Sqrt(pt)+par2_sigl/pt+par3_sigl*pt;
      ts_corr_error[2][i] = (1+meanl)*sigmal*pt;
    }//for loop on leading four jets 
  
  return;
}//

//////////
