#include "DataDrivenBackgrounds.h"

// std includes
#include <iostream>
#include <vector>
#include <string>

// ROOT includes
#include <TFile.h>

/// local includes
#include "Reweight.h"
#include "EventBase.h"


using namespace snu;
using namespace std;


DataDrivenBackgrounds::DataDrivenBackgrounds(){
  
  k_mcperiod=-1;
  corr_isdata=false;
  string lqdir = getenv("LQANALYZER_DIR");
  m_fakeobj = new HNCommonLeptonFakes(lqdir+"/LQAnalysis/src/HNCommonLeptonFakes/share/");
  dd_eventbase= 0;


}


DataDrivenBackgrounds::~DataDrivenBackgrounds(){
  delete m_fakeobj;

}

void DataDrivenBackgrounds::SetEventBase(EventBase* ebase){
  dd_eventbase = ebase;
}

void DataDrivenBackgrounds::CheckEventBase(){
  std::vector<snu::KMuon> muonColl;
  dd_eventbase->GetMuonSel()->SelectMuons(muonColl,"MUON_POG_TIGHT");
  cout << "muon size = " << muonColl.size() << endl;
  if(muonColl.size()  > 0) cout << "new muon size = " << muonColl.size() <<  " " << muonColl.at(0).Pt() << endl;

}



void DataDrivenBackgrounds::SetMCPeriod(int mcperiod){
  k_mcperiod=mcperiod;
}

void DataDrivenBackgrounds::SetIsData(bool isdata){
  corr_isdata=isdata;
}

void DataDrivenBackgrounds::PrintSummary(){
  //// summarize results

  cout << "Functions for ChargeFlip include" << endl;
  cout << "   DataDrivenBackgrounds::SetMCPeriod(int mcperiod) " << endl;
  cout << "   DataDrivenBackgrounds::SetIsData(bool isdata)  " << endl;
  cout << "   DataDrivenBackgrounds::PrintSummary() " << endl;
  cout << "   DataDrivenBackgrounds::WeightCFEvent(std::vector<snu::KElectron> electrons, bool runchargeflip) " << endl;
  cout << "   DataDrivenBackgrounds::ShiftElectronEnergy(std::vector<snu::KElectron> electrons, bool applyshift) " << endl;
  cout << "   DataDrivenBackgrounds::CFRate(snu::KElectron el) " << endl;
  cout << "   DataDrivenBackgrounds::CFRate_Run2(snu::KElectron el, TString el_id) " << endl;
  cout << "   DataDrivenBackgrounds::Get_DataDrivenWeight_EM(vector<snu::KMuon> k_muons, vector<snu::KElectron> k_electrons, vector<bool> istight_e, vector<bool> istight_m, TString cut) " << endl;
  cout << "Functions for NonPrompt include" << endl;
  
  cout << "   DataDrivenBackgrounds::Get_DataDrivenWeight_MM(vector<snu::KMuon> k_muons,vector<bool> istight_m, TString cutID) " << endl;
  cout << "   DataDrivenBackgrounds::Get_DataDrivenWeight_MM(bool geterr, vector<snu::KMuon> k_muons) " << endl;
  cout << "   DataDrivenBackgrounds::Get_DataDrivenWeight_MMM(bool geterr, vector<snu::KMuon> k_muons) " << endl;
  cout << "   DataDrivenBackgrounds::Get_DataDrivenWeight_M(vector<snu::KMuon> k_muons, vector<bool> istight_m, TString cutID) " << endl;
  cout << "   DataDrivenBackgrounds::Get_DataDrivenWeight_E(vector<snu::KElectron> k_electrons,  vector<bool> istight_e) " << endl;
  cout << "   DataDrivenBackgrounds::Get_DataDrivenWeight_EE(vector<snu::KElectron> k_electrons ,  vector<bool> istight_e) " << endl;

  
}


//############################################################################################################
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@  
@  CHARGE FLIP FUNCTIONS
@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

float DataDrivenBackgrounds::WeightCFEvent(std::vector<snu::KElectron> electrons, bool runchargeflip){

  if(electrons.size()!=2) return 0.;
  if(runchargeflip) {
    if(electrons.at(0).Charge() != electrons.at(1).Charge()) {
      float cf1=  CFRate(electrons.at(0));
      float cf2=  CFRate(electrons.at(1));
      return ((cf1/(1.-cf1)) + (cf2/(1.-cf2)));
    }// OS requirement                                                                                                                                                                                                    
    else return 0.;
  }// cf requirement                                                                                                                                                                                                      
  else {
    if(electrons.at(0).Charge() != electrons.at(1).Charge()) return 0.;
  }

  return 1.;

}


std::vector<snu::KElectron> DataDrivenBackgrounds::ShiftElectronEnergy(std::vector<snu::KElectron> electrons, bool applyshift){

  std::vector<snu::KElectron> shiftedel;

  for(unsigned int iel=0; iel < electrons.size(); iel++){
    float scale =0.98;

    if(applyshift)electrons.at(iel).SetPtEtaPhiM(electrons.at(iel).Pt()*scale, electrons.at(iel).Eta(), electrons.at(iel).Phi(), 0.511e-3);
    shiftedel.push_back(electrons.at(iel));
  }
  return shiftedel;
}



float DataDrivenBackgrounds::CFRate(snu::KElectron el){
  if(el.Pt() < 10.) return 0.;
  Double_t frac = 0. ;
  float pt = el.Pt();
  Double_t p0 = 0. ;
  Double_t p1 = 0. ;


  Double_t scale_factor_EE = 1. ;
  Double_t scale_factor_BB = 1. ;

  float eta = el.Eta();

  //--root fitting                                                                                                                                                                                                        
  if( fabs(eta) <= 0.9 ) { // inner BB region                                                                                                                                                                             

    scale_factor_BB = 1.22 ; // BB                                                                                                                                                                                        

    p0 = 3.31e-05 ; p1 = -6.5e-04 ; // root fit                                                                                                                                                                           
    // p0 = 2.8e-05 ; p1 = 0. ;// UK eye fit                                                                                                                                                                              

    frac = p0 + p1*(1./pt) ;
    if( 1./pt < 0.017 ){
      p0 = 1.92e-04 ; p1 = -0.011 ;
      frac = max(p0 + p1*(1./pt), frac);
    }
    frac = max(frac,0.);
    frac *=scale_factor_BB ;

  }else if( fabs(eta) > 0.9 && fabs(eta) <= 1.4442 ){ // outer BB region                                                                                                                                                  
    scale_factor_BB = 1.22 ; // BB                                                                                                                                                                                        
    p0 = 2.21e-04 ; p1 = -5.1e-03 ; // root fit                                                                                                                                                                           
    //    p0 = 1.2e-04 ; p1 = 0. ; // UK eye fit                                                                                                                                                                          
    frac = p0 + p1*(1./pt) ;
    if( 1./pt < 0.02 ){
      p0 = 6.35e-04 ; p1 = -0.027 ;
      frac = max(p0 + p1*(1./pt), frac);
    }
    frac = max(frac,0.);
    frac *=scale_factor_BB ;

  } else {  // fabs(eta) > 1.4                                                                                                                                                                                            


    scale_factor_EE = 1.40 ; //                                                                                                                                                                                           

    //--region:  1/pt > 0.02                                                                                                                                                                                              
    p0 = 4.91e-04 ; p1 = -0.952e-02 ;
    frac = p0 + p1*(1./pt) ;

    if( (1./pt) <= 0.02 ){
      p0 = 2.70e-03 ;  p1 = -1.21e-01 ;
      frac = max(p0 + p1*(1./pt), frac) ;
    }
    frac *= scale_factor_EE ;
  }

  return float(frac) ;

  return 1. ;
}

//CFrate for ELECTRON_HN_TIGHT 13 TeV                                                                                                                                                                                     
float DataDrivenBackgrounds::CFRate_Run2(snu::KElectron el, TString el_id){
  if(el.Pt() < 20.) return 0.;

  Double_t p_B1_1[2] = {0.};
  Double_t p_B2_1[2] = {0.};
  Double_t p_E_1[2] = {0.};

  Double_t p_B1_2[2] = {0.};
  Double_t p_B2_2[2] = {0.};
  Double_t p_E_2[2] = {0.};

  Double_t scale_factor_EE = 1. ;
  Double_t scale_factor_BB = 1. ;

  if(el_id == "ELECTRON_HN_TIGHT"){
    p_B1_1 = {0.0001261, -0.005951};
    p_B2_1 = {0.00063, -0.01965};
    p_E_1 = {0.006827, -0.2198};

    p_B1_2 = {1.584e-05, 5.337e-05};
    p_B2_2 = {3.946e-05, 0.0008439};
    p_E_2 = {0.003153, -0.003891};
  }

  Double_t frac = 0. ;
  float pt_inv = 1. / el.Pt();
  float eta = el.Eta();

  //--root fitting                                                                                                                                                                                                        
  if( fabs(eta) <= 0.9 ) { // inner BB region                                                                                                                                                                             
    //scale_factor_BB = 1.22 ; // BB                                                                                                                                                                                      
    frac = p_B1_2[0] + p_B1_2[1] * pt_inv;
    if( pt_inv < 0.02 ){
      frac = max(p_B1_1[0] + p_B1_1[1] * pt_inv, frac);
    }
    frac = max(frac,0.);
    frac *=scale_factor_BB ;

  }else if( fabs(eta) > 0.9 && fabs(eta) <= 1.4442 ){ // outer BB region                                                                                                                                                  
    //scale_factor_BB = 1.22 ; // BB                                                                                                                                                                                      
    frac = p_B2_2[0] + p_B2_2[1] * pt_inv;
    if( pt_inv < 0.025 ){
      frac = max(p_B2_1[0] + p_B2_1[1] * pt_inv, frac);
    }
    frac = max(frac,0.);
    frac *=scale_factor_BB ;

  } else {  // fabs(eta) > 1.4                                                                                                                                                                                            
    //scale_factor_EE = 1.40;                                                                                                                                                                                             
    frac = p_E_2[0] + p_E_2[1] * pt_inv;
    if( pt_inv <= 0.02 ){
      frac = max(p_E_1[0] + p_E_1[1] * pt_inv, frac);
    }
    frac *= scale_factor_EE ;
  }

  return float(frac) ;
  return 1. ;

}



/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@
@  Non prompt function                                                                                                                                                              
@
@  Get_DataDrivenWeight X will return the weight for channel X using final rates
@
@  Get_TempDataDrivenWeightX will return the weight for functions testing/optimising for dilepton
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */



float DataDrivenBackgrounds::Get_DataDrivenWeight_EM(bool geterr, vector<snu::KMuon> k_muons, vector<snu::KElectron> k_electrons, TString IDmu, TString IDel, TString method){

  // geterr = true : function returns error not event weight  
  // k_electrons are loose electrons defined in analysis code 
  // k_muons are loose muons defined in analysis code 
  // ID defines which id is used for tight leptons              (currently only "POGTight" for electrons) 
  // method defines the techinique used to measure fake rates   (currently only "dijet" for electrons) 
  // @@@@@ CURRENTLY ONLY  POG ID IS USED HERE                                                                                                                                       

  
  if (method != "dijet") return 0;

  if(k_muons.size()!=1) return 0.;
  if(k_electrons.size()!=1) return 0.;

  //// MAKE THIS TAKE IN TSTRING 
  
  /// fix later
  bool is_mu1_tight    =  dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(0),IDmu);
  bool is_el1_tight    =  dd_eventbase->GetElectronSel()->ElectronPass(k_electrons.at(0),IDel);

  vector<TLorentzVector> muons= MakeTLorentz(k_muons);
  vector<TLorentzVector> electrons=MakeTLorentz(k_electrons);
  
  //// cut filled when em fake rates are 

  /// last input is for systematics 
  float em_weight = m_fakeobj->get_dilepton_em_eventweight(geterr,muons, electrons, is_mu1_tight, is_el1_tight);
  
  return em_weight;
}

float DataDrivenBackgrounds::Get_DataDrivenWeight_MM(bool geterr, vector<snu::KMuon> k_muons,   TString IDmu, TString method){
  if(method == "dxy" && IDmu == "MUON_HN_TRI_TIGHT") return Get_DataDrivenWeight_MM(geterr, k_muons);
  
  float mm_weight = 0.;

  if(k_muons.size()!=2) return 0.;
  bool is_mu1_tight=dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(0),IDmu);
  bool is_mu2_tight=dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(1),IDmu);
  vector<TLorentzVector> muons=MakeTLorentz(k_muons);
  mm_weight =m_fakeobj->get_dilepton_mm_eventweight("dijet",geterr, muons, is_mu1_tight,is_mu2_tight,IDmu);

  return 1.;
}

float DataDrivenBackgrounds::Get_DataDrivenWeight_MM(bool geterr, vector<snu::KMuon> k_muons){

  float mm_weight = 0.;

  if(k_muons.size()==2){
    bool is_mu1_tight = dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(0),"MUON_HN_TRI_TIGHT");
    bool is_mu2_tight = dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(1),"MUON_HN_TRI_TIGHT");
    
    vector<TLorentzVector> muons=MakeTLorentz(k_muons);
    
    mm_weight =m_fakeobj->get_dilepton_mm_eventweight(geterr, muons, is_mu1_tight,is_mu2_tight);

  }
  return mm_weight;

}

float DataDrivenBackgrounds::Get_DataDrivenWeight_MMM(bool geterr, vector<snu::KMuon> k_muons){

  float mmm_weight = 0.;

  if(k_muons.size()==3){

    bool is_mu1_tight = dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(0),"MUON_HN_TRI_TIGHT");
    bool is_mu2_tight = dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(1),"MUON_HN_TRI_TIGHT");
    bool is_mu3_tight = dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(2),"MUON_HN_TRI_TIGHT");

    vector<TLorentzVector> muons=MakeTLorentz(k_muons);

    mmm_weight =m_fakeobj->get_trilepton_mmm_eventweight(geterr, muons, is_mu1_tight,is_mu2_tight, is_mu3_tight);

  }
  return mmm_weight;
}

float DataDrivenBackgrounds::Get_DataDrivenWeight(bool geterr, std::vector<snu::KMuon> k_muons, TString muid, int n_muons, std::vector<snu::KElectron> k_electrons, TString elid, int n_electrons, TString elidloose, TString elmethod){

  float this_weight = 0.;

  if( k_muons.size() != n_muons || k_electrons.size() != n_electrons ){
    //Message("[Get_DataDrivenWeight] number of lepton is wrong..", ERROR);
    return 0.;
  }

  std::vector<bool> isT;
  bool AllTight = true;

  for(unsigned int i=0; i<k_muons.size(); i++){
    if( dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(i), muid) ){
      isT.push_back(true);
    }
    else{
      isT.push_back(false);
      AllTight = false;
    }
  }
  for(unsigned int i=0; i<k_electrons.size(); i++){
    if(dd_eventbase->GetElectronSel()->ElectronPass(k_electrons.at(i), elid)){
      isT.push_back(true);
    }
    else{
      isT.push_back(false);
      AllTight = false;
    }
  }

  if(AllTight){
    //Message("[Get_DataDrivenWeight] All leptons pass Tight. Return 0. weight..", DEBUG);
    return 0.;
  }

  std::vector<TLorentzVector> muons=MakeTLorentz(k_muons);
  std::vector<TLorentzVector> electrons=MakeTLorentz(k_electrons);
  TString elkey = GetElFRKey(elidloose, elid, elmethod);
  this_weight =m_fakeobj->get_eventweight(geterr, muons, muid, electrons, elkey , isT);

  return this_weight;
}

float DataDrivenBackgrounds::Get_DataDrivenWeight_M(bool geterr, vector<snu::KMuon> k_muons, TString IDmu, TString method){

  /// two IDs possible  : 

  // HN          : id = "MUON_HN_TIGHT"
  // POG Tight   : id = "MUON_POG_TIGHT

  if(k_muons.size()!=1) return 0.;
  
  bool is_mu1_tight =false;

  TString s_id="fake_Eff_muon_pog"; /// used to select FR Histogram from Map
  if(IDmu == "MUON_HN_TIGHT")    s_id="fake_Eff_muon_hn";

  is_mu1_tight=dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(0),IDmu);


  vector<TLorentzVector> muons=MakeTLorentz(k_muons);

  float f=  m_fakeobj->getFakeRate_muon(geterr,fabs(k_muons.at(0).Eta()), k_muons.at(0).Pt(), s_id);
  float r=  1.;//m_fakeobj->getPromptRate_muon(geterr,fabs(k_muons.at(0).Eta()), k_muons.at(0).Pt(), s_id);

  float w = m_fakeobj->lepton_weight(!is_mu1_tight, r,f);
  return w;

}

float DataDrivenBackgrounds::Get_DataDrivenWeight_E(bool geterr,vector<snu::KElectron> k_electrons,  TString IDe, TString method){

  if(k_electrons.size()!=1) return 0.;

  bool is_el1_tight    = dd_eventbase->GetElectronSel()->ElectronPass(k_electrons.at(0),IDe);

  vector<TLorentzVector> electrons=MakeTLorentz(k_electrons);

  float f=  m_fakeobj->getFakeRate_electronEta(0,k_electrons.at(0).Pt(),fabs(k_electrons.at(0).Eta()),"pt_eta_40_looseregion1");
  float r=  1.;//m_fakeobj->getPromptRate_electron 

  float w = m_fakeobj->lepton_weight(!is_el1_tight, r,f);
  return w;

}

void  DataDrivenBackgrounds::Test(){

  cout << "HNTight_dxy_b017_e011 B = " << m_fakeobj->getFakeRate_electronEta(0,20., 1.2,"dxysig_pt_eta_HNTight_dxy_b017_e011_40") << endl;
  cout << "HNTight_dxy_b015_e011 B= " << m_fakeobj->getFakeRate_electronEta(0,20., 1.2,"dxysig_pt_eta_HNTight_dxy_b015_e011_40") << endl;
  cout << "HNTight_dxy_b014_e011 B= " << m_fakeobj->getFakeRate_electronEta(0,20., 1.2,"dxysig_pt_eta_HNTight_dxy_b014_e011_40") << endl;
  cout << "HNTight_dxy_b013_e011 B= " << m_fakeobj->getFakeRate_electronEta(0,20., 1.2,"dxysig_pt_eta_HNTight_dxy_b013_e011_40") << endl;
  cout << "HNTight_dxy_b017_e017 E = " << m_fakeobj->getFakeRate_electronEta(0,20., 2.2,"dxysig_pt_eta_HNTight_dxy_b017_e017_40") << endl;
  cout << "HNTight_dxy_b017_e015 E= " << m_fakeobj->getFakeRate_electronEta(0,20., 2.2,"dxysig_pt_eta_HNTight_dxy_b017_e015_40") << endl ;
  cout << "HNTight_dxy_b017_e014 E= " << m_fakeobj->getFakeRate_electronEta(0,20., 2.2,"dxysig_pt_eta_HNTight_dxy_b017_e014_40") << endl;
  cout << "HNTight_dxy_b017_e013 E= " << m_fakeobj->getFakeRate_electronEta(0,20., 2.2,"dxysig_pt_eta_HNTight_dxy_b017_e013_40") << endl;
  cout << "HNTight_dxy_b017_e012 E= " << m_fakeobj->getFakeRate_electronEta(0,20., 2.2,"dxysig_pt_eta_HNTight_dxy_b017_e012_40") << endl;
}
float DataDrivenBackgrounds::Get_DataDrivenWeight_EE(bool geterr,vector<snu::KElectron> k_electrons){

  return Get_DataDrivenWeight_EE( geterr,k_electrons, "ELECTRON16_HN_FAKELOOSE_NOD0","ELECTRON16_HN_TIGHT","dijet_ajet40");
}
float DataDrivenBackgrounds::Get_DataDrivenWeight_EE(bool geterr,vector<snu::KElectron> k_electrons,   TString IDloose,TString IDtight, TString method){



  // geterr = true : function returns error not event weight 
  // electrons are loose electrons defined in analysis code 
  // ID defines which id is used for tight leptons              (currently only "ELECTRON_POGTIGHT" for electrons) 
  // method defines the techinique used to measure fake rates   (currently only "dijet" for electrons) 


  if(k_electrons.size()==0) return 0.;

  if(k_electrons.size()==2){

    bool is_el1_tight    =  dd_eventbase->GetElectronSel()->ElectronPass(k_electrons.at(0),IDtight);
    bool is_el2_tight    = dd_eventbase->GetElectronSel()->ElectronPass(k_electrons.at(1),IDtight);
    vector<TLorentzVector> electrons=MakeTLorentz(k_electrons);

    TString elkey = GetElFRKey(IDloose, IDtight,  method);

    return  m_fakeobj->get_dilepton_ee_eventweight(geterr,electrons, is_el1_tight, is_el2_tight, elkey );
      
    /// "" loose ID needs filling here
  }
  return 1.;
}

TString DataDrivenBackgrounds::GetElFRKey( TString IDloose,TString IDtight, TString method){

  TString cut = "pt_eta_";

  if(method.Contains("opt_dijet")){
    string IDtightmod = string(IDtight);
    if(IDtight.Contains("miniiso")){
      std::string sr = "miniiso_";
      std::string::size_type i = IDtightmod.find(sr);
      if (i != std::string::npos)
	IDtightmod.erase(i, sr.length());
    }
    if(IDtight.Contains("dxysig_")){
      std::string sr = "dxysig_";
      std::string::size_type i = IDtightmod.find(sr);
      if (i != std::string::npos)
        IDtightmod.erase(i, sr.length());
    }
    

    TString sjpt="";
    if(method.Contains("ajet20"))sjpt="20";
    if(method.Contains("ajet30"))sjpt="30";
    if(method.Contains("ajet40"))sjpt="40";
    if(method.Contains("ajet60"))sjpt="60";
    else sjpt="40";
    
    if(IDtight.Contains("HNTight_") && (IDtight.Contains("_b") || IDtight.Contains("_e")) ){
      if(IDtight.Contains("dxysig")) {
	if(IDtight.Contains("miniiso"))  cut="miniiso_dxysig_"+cut;
	else cut="dxysig_"+cut;
      }
      else if(IDtight.Contains("miniiso"))  cut ="miniiso_"+cut;
      cut+=IDtightmod+"_"+sjpt;
    }
    else{
      if(method.Contains("dxysig")) cut="dxysig_"+cut;
      cut+=IDtightmod+"_"+sjpt;
    }

    if(method.Contains("_d0")) cut+="_d0";
    return cut;
    
  }

  else if(method.Contains("dijet")){
    TString sjpt="";
      if(method.Contains("ajet20"))sjpt="20_";
      if(method.Contains("ajet30"))sjpt="30_";
      if(method.Contains("ajet40"))sjpt="40_";
      if(method.Contains("ajet60"))sjpt="60_";
      else sjpt="40_";
      
      TString tag="";

      if(IDtight == "ELECTRON16_HN_TIGHT" && IDloose=="ELECTRON16_HN_FAKELOOSE_NOD0")  cut+=sjpt+IDtight + "_dijet_nod0";
      if(IDtight == "ELECTRON16_HN_TIGHT_DXYSIG" && IDloose=="ELECTRON16_HN_FAKELOOSE_NOD0"&& (method.Contains("miniiso")))  cut+=sjpt+IDtight + "_dijet_nod0_dxysig_miniiso";
      else if(IDtight == "ELECTRON16_HN_TIGHT_DXYSIG" && IDloose=="ELECTRON16_HN_FAKELOOSE_NOD0")  cut+=sjpt+IDtight + "_dijet_nod0_dxysig";
      
      if(IDtight == "ELECTRON16_HN_TIGHT" && IDloose=="ELECTRON16_HN_FAKELOOSE")  cut+=sjpt+IDtight + "_dijet_d0";
      if(IDtight == "ELECTRON16_HN_TIGHT_DXYSIG" && IDloose=="ELECTRON16_HN_FAKELOOSE"&& (method.Contains("miniiso")))  cut+=sjpt+IDtight + "_dijet_d0_dxysig_miniiso";
      else if(IDtight == "ELECTRON16_HN_TIGHT_DXYSIG" && IDloose=="ELECTRON16_HN_FAKELOOSE")  cut+=sjpt+IDtight + "_dijet_d0_dxysig";


      if(IDtight == "ELECTRON16_HN_TIGHT" && IDloose=="ELECTRON16_HN_FAKELOOSE_ISO04")  cut+=sjpt+IDtight + "_dijet_iso04";
      if(IDtight == "ELECTRON16_HN_TIGHT" && IDloose=="ELECTRON16_HN_FAKELOOSE_ISO06")  cut+=sjpt+IDtight + "_dijet_iso06";
      if(IDtight == "ELECTRON16_FR_POG_TIGHT" && IDloose=="ELECTRON16_POG_FAKELOOSE")  cut+=sjpt+IDtight + "_dijet_pog";
      if(IDtight == "ELECTRON16_FR_POG_MEDIUM" && IDloose=="ELECTRON16_POG_MEDIUM_FAKELOOSE")  cut+=sjpt+IDtight + "_dijet_pog";
      if(IDtight == "ELECTRON16_FR_POG_TIGHT_CC" && IDloose=="ELECTRON16_POG_FAKELOOSE_CC")  cut+=sjpt+IDtight + "_dijet_pog";
      if(IDtight == "ELECTRON16_FR_POG_MEDIUM_CC" && IDloose=="ELECTRON16_POG_MEDIUM_FAKELOOSE_CC")  cut+=sjpt+IDtight + "_dijet_pog";
      if(IDtight == "ELECTRON16_FR_MVA_TIGHT_CC" && IDloose=="ELECTRON16_MVA_FAKELOOSE_CC")  cut+=sjpt+IDtight + "_dijet_mva";
      if(IDtight == "ELECTRON16_FR_POG_TIGHT_DXYCC" && IDloose=="ELECTRON16_POG_FAKELOOSE_CC")  cut+=sjpt+IDtight + "_dijet_pog";
      if(IDtight == "ELECTRON16_FR_POG_MEDIUM_DXYCC" && IDloose=="ELECTRON16_POG_MEDIUM_FAKELOOSE_CC")  cut+=sjpt+IDtight + "_dijet_pog";
      if(IDtight == "ELECTRON16_FR_MVA_TIGHT_DXYCC" && IDloose=="ELECTRON16_MVA_FAKELOOSE_CC")  cut+=sjpt+IDtight + "_dijet_mva";

      if(IDtight == "ELECTRON16_FR_POG_TIGHT_CC" && IDloose=="ELECTRON16_POG_FAKELOOSE_CC_d0")  cut+=sjpt+IDtight + "_dijet_pog_d0";
      if(IDtight == "ELECTRON16_FR_POG_MEDIUM_CC" && IDloose=="ELECTRON16_POG_MEDIUM_FAKELOOSE_CC_d0")  cut+=sjpt+IDtight + "_dijet_pog_d0";
      if(IDtight == "ELECTRON16_FR_MVA_TIGHT_CC" && IDloose=="ELECTRON16_MVA_FAKELOOSE_CC_d0")  cut+=sjpt+IDtight + "_dijet_mva_d0";
      if(IDtight == "ELECTRON16_FR_POG_TIGHT_DXYCC" && IDloose=="ELECTRON16_POG_FAKELOOSE_CC_d0")  cut+=sjpt+IDtight + "_dijet_pog_d0";
      if(IDtight == "ELECTRON16_FR_POG_MEDIUM_DXYCC" && IDloose=="ELECTRON16_POG_MEDIUM_FAKELOOSE_CC_d0")  cut+=sjpt+IDtight + "_dijet_pog_d0";
      if(IDtight == "ELECTRON16_FR_MVA_TIGHT_DXYCC" && IDloose=="ELECTRON16_MVA_FAKELOOSE_CC_d0")  cut+=sjpt+IDtight + "_dijet_mva_d0";


      return cut;
      
  }

  return "";

}



/// MISC FUNCTIONS

vector<TLorentzVector> DataDrivenBackgrounds::MakeTLorentz(vector<snu::KElectron> el){

  vector<TLorentzVector> tl_el;
  for(vector<KElectron>::iterator itel=el.begin(); itel!=el.end(); ++itel) {
    TLorentzVector tmp_em;
    tmp_em.SetPtEtaPhiM((*itel).Pt(),(*itel).Eta(),(*itel).Phi(),(*itel).M());
    tl_el.push_back(tmp_em);
  }
  return tl_el;
}

vector<TLorentzVector> DataDrivenBackgrounds::MakeTLorentz(vector<snu::KMuon> mu){

  vector<TLorentzVector> tl_mu;
  for(vector<KMuon>::iterator itmu=mu.begin(); itmu!=mu.end(); ++itmu) {
    TLorentzVector tmp_mu;
    tmp_mu.SetPtEtaPhiM((*itmu).Pt(),(*itmu).Eta(),(*itmu).Phi(),(*itmu).M());
    tl_mu.push_back(tmp_mu);
  }
  return tl_mu;
}


vector<TLorentzVector> DataDrivenBackgrounds::MakeTLorentz(vector<snu::KJet> j){

  vector<TLorentzVector> tl_jet;
  for(vector<KJet>::iterator itj=j.begin(); itj!=j.end(); ++itj) {
    TLorentzVector tmp_j;
    tmp_j.SetPtEtaPhiM((*itj).Pt(),(*itj).Eta(),(*itj).Phi(),(*itj).M());
    tl_jet.push_back(tmp_j);
  }
  return tl_jet;
}


vector<TLorentzVector> DataDrivenBackgrounds::MakeTLorentz(vector<snu::KFatJet> j){

  vector<TLorentzVector> tl_jet;
  for(vector<KFatJet>::iterator itj=j.begin(); itj!=j.end(); ++itj) {
    TLorentzVector tmp_j;
    tmp_j.SetPtEtaPhiM((*itj).Pt(),(*itj).Eta(),(*itj).Phi(),(*itj).M());
    tl_jet.push_back(tmp_j);
  }
  return tl_jet;
}
