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


DataDrivenBackgrounds::DataDrivenBackgrounds( ){

  m_fakeobj=0;
  k_mcperiod=-1;
  corr_isdata=false;
  string lqdir = getenv("LQANALYZER_DIR");
  //m_fakeobjtrilep = new HNCommonLeptonFakesTriLep("");
  dd_eventbase= 0;
  UsePtCone=false;

  _muon_isocone=0.4;
  _electron_isocone=0.3;

}

bool DataDrivenBackgrounds::SetupFake(){

  
  if(m_fakeobj) {
    if (m_fakeobj->IsInitialised()) return false;
  }

  
  map<TString, pair<pair<TString,TString>, pair<float,TString> >  > fake_hists;
  fake_hists["fr_muon_central"] = std::make_pair(std::make_pair("Muon_Data_v7_SIP3p5_FR.root","Muon_Data_v7_SIP3p5_FR_Awayjet40"), std::make_pair(70., "TH2D"));
  fake_hists["fr_muon_awayjet20"] = std::make_pair(std::make_pair("Muon_Data_v7_SIP3p5_FR.root","Muon_Data_v7_SIP3p5_FR_Awayjet20"), std::make_pair(70., "TH2D"));
  fake_hists["fr_muon_awayjet30"] = std::make_pair(std::make_pair("Muon_Data_v7_SIP3p5_FR.root","Muon_Data_v7_SIP3p5_FR_Awayjet30"), std::make_pair(70., "TH2D"));
  fake_hists["fr_muon_awayjet60"] = std::make_pair(std::make_pair("Muon_Data_v7_SIP3p5_FR.root","Muon_Data_v7_SIP3p5_FR_Awayjet60"), std::make_pair(70., "TH2D"));
  fake_hists["fr_electron_central"] = std::make_pair(std::make_pair("Electron_Data_v7_FR.root","Electron_Data_v7_FR_Awayjet40") , std::make_pair(70., "TH2D"));

  m_fakeobj = new HNCommonLeptonFakes("/data1/LQAnalyzer_rootfiles_for_analysis/CATAnalysis2016/Fake/DiLep/",fake_hists);
  m_fakeobj->SetMinPtMuon(10.);
  m_fakeobj->SetMinPtElectron(15.);

}

DataDrivenBackgrounds::DataDrivenBackgrounds(bool setupOwnFakes){
  
  if(!setupOwnFakes){
    DataDrivenBackgrounds();
    return;
  }

  k_mcperiod=-1;
  corr_isdata=false;
  string lqdir = getenv("LQANALYZER_DIR");
  
  dd_eventbase= 0;
  _muon_isocone=0.4;
  _electron_isocone=0.3;
  m_fakeobj=0;
  
  
}

bool DataDrivenBackgrounds::SetupFake(TString path, std::map<TString, std::pair<std::pair<TString,TString>  ,std::pair<float,TString> > >fake_hist_config){


  if(m_fakeobj) {
    if (m_fakeobj->IsInitialised()) return false;
  }
  m_fakeobj = new HNCommonLeptonFakes(path, fake_hist_config);
  
  //m_fakeobj->SetupFake(path, fake_hist_config);
  m_fakeobj->SetMinPtMuon(10.);
  m_fakeobj->SetMinPtElectron(15.);
  return true;

}


 DataDrivenBackgrounds::~DataDrivenBackgrounds(){
   delete m_fakeobj;
   // delete m_fakeobjtrilep;
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

   bool debug(false);


   Double_t p_B1_1[2]  = {0.0001261, -0.005951};
   Double_t p_B2_1[2] =  {0.00063, -0.01965};
   Double_t p_E_1[2] = {0.006827, -0.2198};

   Double_t p_B1_2[2] = {1.584e-05, 5.337e-05};
   Double_t p_B2_2[2] = {3.946e-05, 0.0008439};
   Double_t p_E_2[2] =  {0.003153, -0.003891};

   Double_t scale_factor_EE = 1. ;
   Double_t scale_factor_BB = 1. ;


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

   if(debug){cout << "ID = " << el_id << " CF rate = " << endl;}
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



float DataDrivenBackgrounds::Get_DataDrivenWeight_EM(bool geterr, vector<snu::KElectron> k_electrons, vector<snu::KMuon> k_muons, TString IDel, TString IDmu, TString variable,TString tagE, TString tagM, int FakeType){
   if(k_muons.size()!=1) return 0.;
   if(k_electrons.size()!=1) return 0.;

   /// fix later
   bool is_mu1_tight    =  dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(0),IDmu);
   bool is_el1_tight    =  dd_eventbase->GetElectronSel()->ElectronPass(k_electrons.at(0),IDel);

   vector<TLorentzVector> muons= MakeTLorentz(k_muons);
   vector<TLorentzVector> electrons=MakeTLorentz(k_electrons);

   float val1(0.), val2(0.);
   if(variable=="ptcone"){
     val1 = k_electrons[0].PTCone(_electron_isocone,dd_eventbase->GetElectronSel()->IsoCutValue(k_electrons.at(0),IDel));
     val2 = k_muons[0].PTCone(_muon_isocone,dd_eventbase->GetMuonSel()->IsoCutValue(k_muons.at(0),IDmu));     
   }
   else if(variable=="pt"){
     val1=  electrons[0].Pt();
     val2=  muons[0].Pt();
   }

   /// last input is for systematics 
   float em_weight = m_fakeobj->get_dilepton_em_eventweight(geterr, electrons,muons, is_el1_tight, is_mu1_tight, tagM, tagE,val1,val2, FakeType);

   return em_weight;
 }

float DataDrivenBackgrounds::Get_DataDrivenWeight_LL(bool geterr, vector<snu::KElectron> k_electrons, vector<snu::KMuon> k_muons, TString IDel, TString IDmu, float val1, float val2,TString tagE, TString tagM, int FakeType){
						     
  if(k_muons.size()!=1) return 0.;
  if(k_electrons.size()!=1) return 0.;

  /// fix later                                                                                                                                                                     
  bool is_mu1_tight    =  dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(0),IDmu);
  bool is_el1_tight    =  dd_eventbase->GetElectronSel()->ElectronPass(k_electrons.at(0),IDel);

  vector<TLorentzVector> muons= MakeTLorentz(k_muons);
  vector<TLorentzVector> electrons=MakeTLorentz(k_electrons);

  /// last input is for systematics                                                                                                                                                 
  float em_weight = m_fakeobj->get_dilepton_em_eventweight(geterr, electrons,muons, is_el1_tight, is_mu1_tight, tagM, tagE,val1,val2, FakeType);

  return em_weight;
}





float DataDrivenBackgrounds::Get_DataDrivenWeight_MM(bool geterr, vector<snu::KMuon> k_muons,  TString ID, TString variable, TString tag , int FakeType){
  
  float mm_weight = 0.;

  if(k_muons.size()==2){
    bool is_mu1_tight = dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(0),ID);
    bool is_mu2_tight = dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(1),ID);
    
    vector<TLorentzVector> muons=MakeTLorentz(k_muons);
    
    float val1(0.), val2(0.);
    
    if(variable=="ptcone"){
     val1 = k_muons[0].PTCone(_muon_isocone,dd_eventbase->GetMuonSel()->IsoCutValue(k_muons.at(0),ID));
     val2 = k_muons[1].PTCone(_muon_isocone,dd_eventbase->GetMuonSel()->IsoCutValue(k_muons.at(1),ID));
    }
    else if(variable=="pt"){
      val1=  muons[0].Pt();
      val2=  muons[0].Pt();
    }
    mm_weight =m_fakeobj->get_dilepton_mm_eventweight(geterr, muons, is_mu1_tight,is_mu2_tight,tag, val1, val2, FakeType);
    
  }
  return mm_weight;
  
}

float DataDrivenBackgrounds::Get_DataDrivenWeight_LL(bool geterr, vector<snu::KMuon> k_muons,  TString ID, float val1, float val2, TString tag , int FakeType){

  float ll_weight = 0.;

  if(k_muons.size()==2){
    bool is_mu1_tight = dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(0),ID);
    bool is_mu2_tight = dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(1),ID);

    vector<TLorentzVector> muons=MakeTLorentz(k_muons);
    ll_weight =m_fakeobj->get_dilepton_mm_eventweight(geterr, muons, is_mu1_tight,is_mu2_tight,tag, val1, val2, FakeType);
  }
  
  return ll_weight;

}




float DataDrivenBackgrounds::Get_DataDrivenWeight_MMM(bool geterr, vector<snu::KMuon> k_muons){
  
  float mmm_weight = 0.;
  
  if(k_muons.size()==3){
    
    bool is_mu1_tight = dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(0),"MUON_HN_TRI_TIGHT");
    bool is_mu2_tight = dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(1),"MUON_HN_TRI_TIGHT");
    bool is_mu3_tight = dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(2),"MUON_HN_TRI_TIGHT");
    
    vector<TLorentzVector> muons=MakeTLorentz(k_muons);
    
    if(m_fakeobjtrilep)mmm_weight =m_fakeobjtrilep->get_trilepton_mmm_eventweight(geterr, muons, is_mu1_tight,is_mu2_tight, is_mu3_tight);
    
  }
  return mmm_weight;
}

float DataDrivenBackgrounds::Get_DataDrivenWeight(bool geterr, std::vector<snu::KMuon> k_muons, TString muid,  int n_muons, std::vector<snu::KElectron> k_electrons, TString elid, int n_electrons, TString elidloose, TString method, int HalfSampleErrorDir){
  
  float this_weight = 0.;
  
  if( int(k_muons.size()) != n_muons || int(k_electrons.size()) != n_electrons ){
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
  
  if(m_fakeobjtrilep)m_fakeobjtrilep->SetUsePtCone(UsePtCone);
  
  std::vector<TLorentzVector> muons=MakeTLorentz(k_muons);
  std::vector<TLorentzVector> electrons=MakeTLorentz(k_electrons);
  std::vector<TString> vkeys;
  for(unsigned int ikey=0; ikey < electrons.size(); ikey++){
    TString regel1="reg1";
    if(k_electrons.at(ikey).Pt() > 50.) {
      if(k_electrons.at(ikey).IsEB1()) regel1="reg2";
      if(k_electrons.at(ikey).IsEB2()) regel1="reg3";
      if(k_electrons.at(ikey).IsEE()) regel1="reg4";
    }
    vkeys.push_back(regel1);
  }
  
  std::vector<TString> keys = GetElFRKey(elidloose, elid, method, vkeys);
  
  if(m_fakeobjtrilep)this_weight =m_fakeobjtrilep->get_eventweight(geterr, muons, muid, electrons,  keys , isT, HalfSampleErrorDir);
  
  return this_weight;
}

float DataDrivenBackgrounds::Get_DataDrivenWeight_M(bool geterr, vector<snu::KMuon> k_muons, TString IDmu, TString tag){
  
  /// two IDs possible  : 
  if(k_muons.size()!=1) return 0.;
  
  bool is_mu1_tight =is_mu1_tight=dd_eventbase->GetMuonSel()->MuonPass(k_muons.at(0),IDmu);
  
  vector<TLorentzVector> muons=MakeTLorentz(k_muons);
  
  float f=  m_fakeobj->getFakeRate(k_muons[0].PTCone(_muon_isocone,dd_eventbase->GetMuonSel()->IsoCutValue(k_muons.at(0),IDmu)),fabs(k_muons.at(0).Eta()), tag,10.);
  float r=  1.;
  
  bool debug(false);
  
  float w = m_fakeobj->lepton_weight(!is_mu1_tight, r,f);
  return w;
  
}

float DataDrivenBackgrounds::GetFakeRateEl(float pt, float eta, TString tag){
  return m_fakeobj->getFakeRate(pt,eta, tag, 15.);
}

float DataDrivenBackgrounds::Get_DataDrivenWeight_E(bool geterr,vector<snu::KElectron> k_electrons,  TString IDe, TString tag){
  
  if(k_electrons.size()!=1) return 0.;

  bool is_el1_tight    = dd_eventbase->GetElectronSel()->ElectronPass(k_electrons.at(0),IDe);
  
  vector<TLorentzVector> electrons=MakeTLorentz(k_electrons);

  float f=  m_fakeobj->getFakeRate(k_electrons.at(0).Pt(),fabs(k_electrons.at(0).Eta()),tag,15.);
  float r=  1.;
  
  if(geterr) return m_fakeobj->getFakeRate(k_electrons.at(0).Pt(),fabs(k_electrons.at(0).Eta()),tag,15.,true);

  float w = m_fakeobj->lepton_weight(!is_el1_tight, r,f);
  return w;

}

void  DataDrivenBackgrounds::Test(){

}


float DataDrivenBackgrounds::Get_DataDrivenWeight_EE(bool geterr, vector<snu::KElectron> k_electrons,  TString ID, TString variable, TString tag , int FakeType){

  float ee_weight = 0.;

  if(k_electrons.size()==2){
    bool is_el1_tight = dd_eventbase->GetElectronSel()->ElectronPass(k_electrons.at(0),ID);
    bool is_el2_tight = dd_eventbase->GetElectronSel()->ElectronPass(k_electrons.at(1),ID);

    vector<TLorentzVector> electrons=MakeTLorentz(k_electrons);

    float val1(0.), val2(0.);

    if(variable=="ptcone"){
      val1 = k_electrons[0].PTCone(_electron_isocone,dd_eventbase->GetElectronSel()->IsoCutValue(k_electrons.at(0),ID));
      val2 = k_electrons[1].PTCone(_electron_isocone,dd_eventbase->GetElectronSel()->IsoCutValue(k_electrons.at(1),ID));
    }
    else if(variable=="pt"){
      val1=  k_electrons[0].Pt();
      val2=  k_electrons[0].Pt();
    }
    ee_weight =m_fakeobj->get_dilepton_ee_eventweight(geterr, electrons, is_el1_tight,is_el2_tight,tag, val1, val2, FakeType);

  }
  return ee_weight;
}

float DataDrivenBackgrounds::Get_DataDrivenWeight_LL(bool geterr, vector<snu::KElectron> k_electrons,  TString ID, float val1, float val2, TString tag , int FakeType){

  float ee_weight = 0.;

  if(k_electrons.size()==2){
    bool is_el1_tight = dd_eventbase->GetElectronSel()->ElectronPass(k_electrons.at(0),ID);
    bool is_el2_tight = dd_eventbase->GetElectronSel()->ElectronPass(k_electrons.at(1),ID);

    vector<TLorentzVector> electrons=MakeTLorentz(k_electrons);
    ee_weight =m_fakeobj->get_dilepton_ee_eventweight(geterr, electrons, is_el1_tight,is_el2_tight,tag, val1, val2, FakeType);

  }
  return ee_weight;
}



vector<TString> DataDrivenBackgrounds::GetElFRKey( TString IDloose,TString IDtight, TString method, std::vector<TString> regs1){
  
  if(!method.Contains("mva")) {
    vector<TString> tmpkeys;
    tmpkeys.push_back( GetElFRKey(IDloose, IDtight, method));
    return tmpkeys;
  }

  vector<TString> keys_forfakes;
  TString cut = "";

  for(unsigned int i=0; i < regs1.size(); i++){
    int reg=0;
    if(regs1[i].Contains("reg1"))reg=1;
    if(regs1[i].Contains("reg2"))reg=2;
    if(regs1[i].Contains("reg3"))reg=3;
    if(regs1[i].Contains("reg4"))reg=4;

    if(IDtight == "ELECTRON_HN_TIGHT" &&  IDloose=="ELECTRON_HN_FAKELOOSE"&&reg==1) cut="0.95_iso0.05_dxy0.01_dz0.04";
    if(IDtight == "ELECTRON_HN_TIGHT" &&  IDloose=="ELECTRON_HN_FAKELOOSE"&&reg==2) cut="0.76_iso0.05_dxy0.01_dz0.04";
    if(IDtight == "ELECTRON_HN_TIGHT" &&  IDloose=="ELECTRON_HN_FAKELOOSE"&&reg==3) cut="0.72_iso0.05_dxy0.01_dz0.04";
    if(IDtight == "ELECTRON_HN_TIGHT" &&  IDloose=="ELECTRON_HN_FAKELOOSE"&&reg==4) cut="0.7_iso0.05_dxy0.01_dz0.04";
  
    keys_forfakes.push_back(cut);
  }
  return keys_forfakes;
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

      
      /// tmp
      if(IDtight == "ELECTRON_POG_TIGHT" && IDloose=="ELECTRON_POG_FAKELOOSE")  cut+=sjpt+"ELECTRON16_FR_POG_TIGHT_CC_dijet_pog_d0"; /// NEW from .10 tag POG Tight ID == 

      if(IDtight == "ELECTRON_MVA_TIGHT" && IDloose=="ELECTRON_MVA_FAKELOOSE")  cut+=sjpt+"ELECTRON16_FR_MVA_TIGHT_CC_dijet_mva_d0";
      
      return cut;
      
  }

  return "";

}



/// MISC FUNCTIONS

void DataDrivenBackgrounds::SetUsePtCone(bool b){

  UsePtCone = b;

}

double DataDrivenBackgrounds::MuonConePt(snu::KMuon muon, double tightiso){

  double mu_pt_corr = muon.Pt()*(1+max(0.,(muon.RelIso04()-tightiso))) ;

  return mu_pt_corr;

}

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
    double this_pt = (*itmu).Pt();

    if(UsePtCone){
      this_pt = MuonConePt((*itmu), 0.1);
    }

    tmp_mu.SetPtEtaPhiM(this_pt,(*itmu).Eta(),(*itmu).Phi(),(*itmu).M());
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
