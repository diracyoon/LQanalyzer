#include "SignalPlotsMM.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

SignalPlotsMM::SignalPlotsMM(TString name, int nmu): StdPlots(name){

  TH1::SetDefaultSumw2(true);


  bool singlemu = (nmu < 0)  ?   true : (nmu == 1);
  bool doublemu = (nmu < 0)  ?   true : (nmu == 2);
  bool trimu = (nmu < 0)  ?   true : (nmu == 3);
  bool quadmu = (nmu < 0)  ?   true : (nmu == 4);

  /// Mass plots                                                                                                                                                                                                                                                                                                                              
  if (quadmu){
    map_sig["h_zzmass"]               = SetupHist("h_zzmass_"          + name,"Invariant mass of the os pairs" ,50,100,150,"M_{ll} GeV");
    map_sig["h_llllmass"]                 = SetupHist("h_llllmass_"           + name,"Invariant mass of the two leading muectrons",200,0,1000,"M_{llll} GeV");
    map_sig["h_z1mass"]               = SetupHist("h_z1mass_"          + name,"Invariant mass of the os pairs" ,80,0,400,"M_{ll} GeV");
    map_sig["h_z2mass"]               = SetupHist("h_z2mass_"          + name,"Invariant mass of the os pairs" ,80,0,400,"M_{ll} GeV");
    map_sig["h_llllpt"]                 = SetupHist("h_llllpt_"           + name,"Invariant pt of the two leading muectrons",500,0,1000,"P_T (llll) GeV");
    map_sig["h_z1_z2_dphi"]                 = SetupHist("h_z1_z2_dphi_"   + name,"Delta Phi Z1,Z2", 28, 0., 3.5,"#Delta(#phi) (Z_{1},Z_{2})");
  }

  
  if(doublemu){
    map_sig["h_llmass"]                 = SetupHist("h_llmass_"           + name,"Invariant mass of the two leading muectrons",1000,0,1000,"M_{ll} GeV");
    map_sig["h_llpt"]                 = SetupHist("h_llpt_"           + name,"Invariant pt of the two leading muectrons",500,0,1000,"P_T (ll) GeV");
    map_sig["h_lljjmass"]               = SetupHist("h_lljjmass_"         + name,"Invariant mass of the four particles",200,0,2000,"M_{lljj} GeV");
    map_sig["h_llfjmass"]               = SetupHist("h_llfjmass_"         + name,"Invariant mass of the four particles",200,0,2000,"M_{lljj} GeV");
    map_sig["h_lljjjj_ss_mass"]          = SetupHist("h_lljjjj_ss_mass_"       + name,"Invariant mass of the four particles",200,0,2000,"M_{lljjjj} (ss) GeV");
    map_sig["h_lljjjj_mass"]             = SetupHist("h_lljjjj_mass_"       + name,"Invariant mass of the four particles",200,0,2000,"M_{lljjjj} GeV");
    map_sig["h_lljmass"]                = SetupHist("h_lljmass_"         + name,"Invariant mass of the four particles",200,0,2000,"M_{llj} GeV");
    map_sig["h_l1jj_central_mass"]               = SetupHist("h_l1jj_centralmass_"          + name,"Invariant mass of the two leading jets and leading muectron",100,0,1000);
    map_sig["h_l2jj_central_mass"]               = SetupHist("h_l2jj_centralmass_"          + name,"Invariant mass of the two leading jets and leading muectron",100,0,1000);
    map_sig["h_lljj_central_mass"]               = SetupHist("h_lljj_centralmass_"          + name,"Invariant mass of the two leading jets and leading muectron",100,0,1000);


    map_sig2["h_l1l2jjmass"] = SetupHist2("h_l1l2jjmass_" +  name , "", 100,0,1000,100,0,1000);
    map_sig["h_l1_l2jjmass"] = SetupHist("h_l1_l2jjmass_"       + name,"Invariant mass of the two leading jets and leading muectron",200,0,10);
    map_sig["h_forward_jet_pt"]             = SetupHist("h_forward_jet_pt_"                + name,"h_forward_jet_pt",60,0,300);
    map_sig["h_central_jet_pt"]             = SetupHist("h_central_jet_pt_"                + name,"h_central_jet_pt",60,0,300);
    map_sig["h_forward_jet_eta"]             = SetupHist("h_forward_jet_eta_"                + name,"h_forward_jet_eta",50,-2.5,2.5);
    map_sig["h_central_jet_eta"]             = SetupHist("h_central_jet_eta_"                + name,"h_central_jet_eta",50,-2.5,2.5);


  }
  else{

  }
  map_sig["h_sumeta"]     = SetupHist("h_sumeta_"   + name,"", 100, -10., 10.);


 
  map_sig["h_jjmass"]                 = SetupHist("h_dijetsmass_"        + name,"Invariant mass of the two leading jets",500,0.,500.,"M_{jj} GeV");
  map_sig["h_fjmass"]                 = SetupHist("h_fatjetsmass_"        + name,"Invariant mass of the two leading jets",500,0.,500.,"M_{jj} GeV");
  map_sig["h_fjprunedmass"]                 = SetupHist("h_fatjetprunedmass_"        + name,"Invariant mass of the two leading jets",500,0.,500.,"M_{jj} GeV");
  map_sig["h_leadjetmass"]            = SetupHist("h_leadjetsmass_"      + name,"Invariant mass of the two leading jets",150,0.,750.,"M_{j1j2} GeV");
  map_sig["h_osllmass"]               = SetupHist("h_osllmass_"          + name,"Invariant mass of the two leading os muectrons",100,0,500,"M_{ll} GeV");
  map_sig["h_osll_closestZ_mass"]       =SetupHist("h_osll_closestZ_mass_"          + name,"Invariant mass of the two leading os muectrons",100,0,500,"M_{ll} GeV");

  map_sig["h_l1jjmass"]               = SetupHist("h_l1jjmass_"          + name,"Invariant mass of the two leading jets and leading muectron",100,0,1000,"M_{l1jj} GeV");
  map_sig["h_l1fjmass"]               = SetupHist("h_l1fjmass_"          + name,"Invariant mass of the two leading jets and leading muectron",100,0,1000,"M_{l1jj} GeV");
  if(nmu > 1 || nmu < 0) map_sig["h_l2jjmass"]               = SetupHist("h_l2jjmass_"          + name,"Invariant mass of the two leading jets and second muectron",100,0,1000,"M_{l2jj} GeV");
  if(nmu > 1|| nmu < 0) map_sig["h_l2fjmass"]               = SetupHist("h_l2fjmass_"          + name,"Invariant mass of the two leading jets and second muectron",100,0,1000,"M_{l2jj} GeV");

  if(trimu){
    map_sig["h_lllmass"]                 = SetupHist("h_lllmass_"           + name,"Invariant mass of the two leading muectrons",1000,0,1000,"M_{lll} GeV");
    map_sig["h_lllpt"]                 = SetupHist("h_lllpt_"           + name,"Invariant pt of the two leading muectrons",500,0,1000,"P_T (lll) GeV");
  }



  map_sig["h_l1jjjj_mass"]             = SetupHist("h_l1jjjj_mass_"       + name,"Invariant mass of the four particles",200,0,2000,"M_{l1jjjj} GeV");
  map_sig["h_lljjjj_os_mass"]          = SetupHist("h_lljjjj_os_mass_"       + name,"Invariant mass of the four particles",200,0,2000,"M_{lljjjj} (os) GeV");
  map_sig["h_l1jjjj_os_mass"]          = SetupHist("h_l1jjjj_os_mass_"       + name,"Invariant mass of the four particles",200,0,2000,"M_{l1jjjj} (os) GeV");
  map_sig["h_l1jjjj_ss_mass"]          = SetupHist("h_l1jjjj_ss_mass_"       + name,"Invariant mass of the four particles",200,0,2000,"M_{l1jjjj} (ss) GeV");
  if(nmu!= 1){
    map_sig["h_l2jjjj_mass"]             = SetupHist("h_l2jjjj_mass_"       + name,"Invariant mass of the four particles",200,0,2000,"M_{l2jjjj} GeV");
    map_sig["h_l2jjjj_os_mass"]          = SetupHist("h_l2jjjj_os_mass_"       + name,"Invariant mass of the four particles",200,0,2000,"M_{l2jjjj} (os) GeV");
    map_sig["h_l2jjjj_ss_mass"]          = SetupHist("h_l2jjjj_ss_mass_"       + name,"Invariant mass of the four particles",200,0,2000, "M_{l2jjjj} (s) GeV");
  }


  map_sig["h_st_forward"]          = SetupHist("h_st_forward_"       + name,"ST foward jets", 100, 0., 1000.,"ST_{fj} GeV");
  map_sig["h_st_central"]          = SetupHist("h_st_central_"       + name,"ST central jets", 100, 0., 1000.,"ST_{cj} GeV");
  map_sig["h_st_fc_ratio"]          = SetupHist("h_st_fc_ratio_"       + name,"ST ratio", 100, 0., 2.,"ST (ratio)");

  map_sig["h_ht_forward"]          = SetupHist("h_ht_forward_"       + name,"ST foward jets", 100, 0., 1000.,"HT_{fj} GeV");
  map_sig["h_ht_central"]          = SetupHist("h_ht_central_"       + name,"ST central jets", 100, 0., 1000.,"HT_{cj} GeV");
  map_sig["h_ht_fc_ratio"]          = SetupHist("h_ht_fc_ratio_"       + name,"ST ratio", 1000, 0., 10.,"HT (ratio)");


  map_sig["h_eltype" ]          = SetupHist("h_eltype_"       + name,"Type", 45, 0., 45.,"El Type");
  map_sig["h_lt" ]          = SetupHist("h_lt_"       + name,"LT", 100., 0., 500.,"LT");
  map_sig["h_lt_ht" ]          = SetupHist("h_lt_ht_"       + name,"LT", 100., 0., 5.,"LT");
  map_sig["h_lt_forward" ]          = SetupHist("h_lt_forward_"       + name,"LT", 100., 0., 500.,"LT");
  map_sig["h_lt_central" ]          = SetupHist("h_lt_central_"       + name,"LT", 100., 0., 500.,"LT");
  map_sig["h_lt_rforward" ]          = SetupHist("h_lt_rforward_"       + name,"LT", 100., 0., 5.,"LT");

  map_sig["h_mass_forward"]          = SetupHist("h_mass_forward_"       + name,"Mass Forward jets", 100, 0., 1000.,"M_{fj} GeV");
  map_sig["h_mass_central"]          = SetupHist("h_mass_central_"       + name,"Mass Forward jets", 100, 0., 1000.,"M_{cj} GeV");


  map_sig["h_centralNJets"]          = SetupHist("h_centralNJets_"   + name,"NFJ",10, 0, 10.,"N_{CJ}");
  map_sig["h_forwardNJets"]          = SetupHist("h_forwardNJets_"   + name,"NFJ",10, 0, 10.,"N_{FJ}");

  map_sig["h_max_deta"]  = SetupHist("h_max_deta_"  + name,"max_deta", 100, 0., 10.,"#Delta #Eta");
  map_sig["h_llfj_mass"]          = SetupHist("h_llfj_mass_" + name , "h_llfj_mass_", 100, 0., 1000.,"M_{llfj} GeV");
  map_sig["h_llcj_mass"]          = SetupHist("h_llcj_mass_"  + name,  "h_llfj_mass_",100, 0., 1000.,"M_{llfj} GeV");
  map_sig["h_llDfj_mass"]          = SetupHist("h_llDfj_mass_"  + name, "h_llfj_mass_", 100, 0., 1000.,"M_{llDfj} GeV");

  /// ALL Lepton plots                                                                                                                                                                                                                                                                                                                        
  map_sig["h_LeptonEta"]              = SetupHist("h_LeptonEta_"         + name,"leading lepton eta",60,-3.,3., "#eta_{l}");
  map_sig["h_LeptonPhi"]              = SetupHist("h_LeptonPhi_"         + name,"leading lepton phi",64,-3.2,3.2);
  map_sig["h_LeptonPt"]               = SetupHist("h_LeptonPt_"          + name,"lepton pt",250,0,500);
  map_sig["h_LeptonE"]               = SetupHist("h_LeptonE_"          + name,"lepton e",500,0,500);
  map_sig["h_leadingLeptonPt"]        = SetupHist("h_leadingLeptonPt_"   + name,"leading lepton pt",100,0,500);

  map_sig2["h_lt_jjdr"]  = SetupHist2("h_lt_jjdr_"  + name,"eta " , 100., 0.,2000., 30.,0., 4.5);
  map_sig2["h_lt_CM"]   = SetupHist2("h_lt_CM_"  + name,"eta " , 100., 0.,2000., 100, 0., 1000.);
  map_sig2["h_lt_lljjmass"]   = SetupHist2("h_lt_lljjmass_"  + name,"eta " , 100., 0.,2000., 100, 0.,2000.);
  map_sig2["h_lt_lmeanjjmass"]   = SetupHist2("h_lt_lmeanjjmass_" + name,"eta " , 100., 0.,2000., 100, 0.,1000.);
  map_sig2["h_lt_ldiffmeanjjmass"]   = SetupHist2("h_lt_ldiffmeanjjmass_" + name,"eta " , 100., 0.,2000., 100, 0.,1000.);


  map_sig2["EtaRegion"]           = SetupHist2("EtaRegion_el_"  + name,"eta " , 8, 0., 8., 8, 0.,8.);
  map_sig2["JetEtaRegion"]           = SetupHist2("JetEtaRegion_el_"  + name,"eta " , 8, 0., 8., 8, 0.,8.);

  if(nmu > 1 || nmu < 0) {
    map_sig["h_secondLeptonPt"]         = SetupHist("h_secondLeptonPt_"    + name,"secondary lepton pt",60,0,300);
    map_sig["h_secondLeptonEta"]        = SetupHist("h_secondLeptonEta_"   + name,"second lepton eta",60,-3.,3.);
    map_sig["h_secondLeptonRelIso_B"]     = SetupHist("h_secondLeptonRelIso_B"      + name,"leading lepton relIso",100,0,1.);
    map_sig["h_secondLeptonRelIso_EC"]     = SetupHist("h_secondLeptonRelIso_EC"      + name,"leading lepton relIso",100,0,1.);
    map_sig["h_secondLeptonRelIso"]     = SetupHist("h_secondLeptonRelIso_"      + name,"leading lepton relIso",100,0,1.);
    map_sig["h_secondLeptonMiniRelIso"]     = SetupHist("h_secondLeptonMiniRelIso_"      + name,"leading lepton mini relIso",100,0,1.);
    if(nmu > 2|| nmu < 0 ) {
      map_sig["h_thirdLeptonPt"]         = SetupHist("h_thirdLeptonPt_"    + name,"thirdary lepton pt",60,0,300);
      map_sig["h_thirdLeptonEta"]        = SetupHist("h_thirdLeptonEta_"   + name,"third lepton eta",60,-3.,3.);
      map_sig["h_thirdLeptonRelIso"]     = SetupHist("h_thirdLeptonRelIso_"      + name,"leading lepton relIso",100,0,1.);
      map_sig["h_thirdLeptonMiniRelIso"]     = SetupHist("h_thirdLeptonMiniRelIso_"      + name,"leading lepton mini relIso",100,0,1.);
      if(nmu > 3 || nmu < 0) {
        map_sig["h_fourthLeptonPt"]         = SetupHist("h_fourthLeptonPt_"    + name,"fourthary lepton pt",60,0,300);
	map_sig["h_fourthLeptonEta"]        = SetupHist("h_fourthLeptonEta_"   + name,"fourth lepton eta",60,-3.,3.);
        map_sig["h_fourthLeptonRelIso"]     = SetupHist("h_fourthLeptonRelIso_"      + name,"leading lepton relIso",100,0,1.);
        map_sig["h_fourthLeptonMiniRelIso"]     = SetupHist("h_fourthLeptonMiniRelIso_"      + name,"leading lepton mini relIso",100,0,1.);
      }
    }

  }
  map_sig["h_leadingLeptonEta"]       = SetupHist("h_leadingLeptonEta_"  + name,"leading lepton eta",60,-3.,3.);

  map_sig["h_LeptonIso"]              = SetupHist("h_LeptonIso_"         + name,"leading lepton relIso",50,0,10.);
  map_sig["h_LeptonRelIso"]           = SetupHist("h_LeptonRelIso_"      + name,"leading lepton relIso",100,0,1.);
  map_sig["h_leadingLeptonRelIso_B"]    = SetupHist("h_leadingLeptonRelIso_B"      + name,"leading lepton relIso",100,0,1.);
  map_sig["h_leadingLeptonRelIso_EC"]    = SetupHist("h_leadingLeptonRelIso_EC"      + name,"leading lepton relIso",100,0,1.);
  map_sig["h_leadingLeptonMiniRelIso"]    = SetupHist("h_leadingLeptonMiniRelIso_"      + name,"leading lepton mini relIso",100,0,1.);
  map_sig["h_LeptonDZ_EC"]               = SetupHist("h_LeptonDZ_EC_"          + name," leading lepton DZ", 1000, -10. , 10.);
  map_sig["h_LeptonDZ_B"]               = SetupHist("h_LeptonDZ_B_"          + name," leading lepton DZ", 1000, -10. , 10.);
  map_sig["h_LeptonDXY_EC"]              = SetupHist("h_LeptonDXY_EC_"         + name," lepton DXY_{ec}", 1000, -0.25 , 0.25);
  map_sig["h_LeptonDXY_B"]              = SetupHist("h_LeptonDXY_B_"         + name," lepton DXY_{b}", 1000, -0.25 , 0.25);
  map_sig["h_LeptonDXY"]              = SetupHist("h_LeptonDXY_"         + name," leading lepton DXY", 1000, -1. , 1.25);
  map_sig["h_LeptonDXYSig"]              = SetupHist("h_LeptonDXYSig_"         + name," leading lepton DXYSig", 1000, -50. , 50.);


  /// Lepton/Jet/Lepton dR                                                                                                                                                                                                                                                                                                                    
  map_sig["h_LeptonJetdR"]            = SetupHist("h_LeptonJetdR_"       + name,"lepton jet dR",50,0,5);
  map_sig["h_LeptonAwayJetdR"]        = SetupHist("h_LeptonAwayJetdR_"   + name, "away jet dr", 50,0,5);
  map_sig["h_leadLeptonAwayJetRatio"]        = SetupHist("h_leadLeptonAwayJetRatio_"   + name, "away jet / lepton", 100,0,5);
  map_sig["h_secondLeptonAwayJetRatio"]        = SetupHist("h_secondLeptonAwayJetRatio_"   + name, "away jet / lepton", 100,0,5);
  map_sig["h_leadLeptondR"]           = SetupHist("h_leadLeptondR_"      + name,"leading lepton dR",50,0,5);
  map_sig["h_leadLeptonstarjetdR"]           = SetupHist("h_leadLeptonstarjetdR_"      + name,"leading lepton dR",500,0,25);
  map_sig["h_leadLeptondPhi"]         = SetupHist("h_leadLeptondPhi_"    + name,"leading lepton dPhi",50,-5,5);
  map_sig["h_leadJetdR"]                = SetupHist("h_leadJetdR_"         + name,"leading jet dR",50,0,5);
  map_sig["h_leadJetdPhi"]              = SetupHist("h_leadJetdPhi_"       + name,"leading jet dPhi",50,-5,5);
  map_sig["h_leadJetFatJetdR"]                = SetupHist("h_leadJetFatJetdR_"         + name,"leading jet dR",50,0,5);
  map_sig["h_leadJetFatJetdPhi"]              = SetupHist("h_leadJetFatJetdPhi_"       + name,"leading jet dPhi",50,-5,5);


  //// Jet Plots                                                                                                                                                                

  map_sig["h_FatJetPt"]             = SetupHist("h_FatJetPt_"      + name,"leading jet pt",100,0,1000);
  map_sig["h_FatJetEta"]            = SetupHist("h_FatJetEta_"      + name,"leading jet pt",120,-5,5);
  map_sig["h_FatJetPhi"]            = SetupHist("h_FatJetPhi_"      + name,"leading jet pt",64,-3.2,3.2);

  map_sig["h_leadingJetPt"]           = SetupHist("h_leadingJetPt_"      + name,"leading jet pt",60,0,300);
  map_sig["h_secondJetPt"]            = SetupHist("h_secondJetPt_"       + name,"secondary jet pt",60,0,300);
  map_sig["h_HT"]                     = SetupHist("h_HT_"                + name,"sum jet pt",50,0,1000);
  map_sig["h_ST"]                     = SetupHist("h_ST_"                + name,"sum event pt",50,0,1500);
  map_sig["h_STall"]                     = SetupHist("h_STall_"                + name,"sum event pt",50,0,1500);
  map_sig["h_jets_pt"]                = SetupHist("h_jets_pt_"           + name,"jet pt",60,0,300);
  map_sig["h_jets_eta"]               = SetupHist("h_jets_eta_"          + name,"#eta distribution of the two jets",120,-3,3);
  map_sig["h_jets_phi"]               = SetupHist("h_jets_phi_"          + name,"#phi distribution of the two jets",140,-3.5,3.5);
  map_sig["h_PileupJetIDMVA"]         = SetupHist("h_pileupJetIDMVA_"    + name, "" ,100, -1.,1.);
  map_sig["h_bTag"]                   = SetupHist("h_bTag_"              + name,"bTag discrimant",100,-1,3);
  map_sig["h_Njets"]                  = SetupHist("h_Njets_"             + name,"number of jets",10,0,10);
  map_sig["h_Nfatjets"]                  = SetupHist("h_Nfatjets_"             + name,"number of jets",10,0,10);
  map_sig["h_Nbjets_l"]                 = SetupHist("h_Nbjets_l_"            + name,"number of b jets",5,0,5);
  map_sig["h_Nbjets_m"]                 = SetupHist("h_Nbjets_m_"            + name,"number of b jets",5,0,5);
  map_sig["h_Nbjets_t"]                 = SetupHist("h_Nbjets_t_"            + name,"number of b jets",5,0,5);

  map_sig["h_Nbjets_all_l"]                 = SetupHist("h_Nbjets_all_l_"            + name,"number of b jets",5,0,5);
  map_sig["h_Nbjets_all_m"]                 = SetupHist("h_Nbjets_all_m_"            + name,"number of b jets",5,0,5);
  map_sig["h_Nbjets_all_t"]                 = SetupHist("h_Nbjets_all_t_"            + name,"number of b jets",5,0,5);

  map_sig["h_lep_jet_dphi"]           =      SetupHist("h_lep_jet_dphi_"            + name,"dphi lepton jets", 50, 0., 5.);
  map_sig["h_leading_jets_eta"]               = SetupHist("h_leading_jets_eta_"          + name,"#eta distribution of the two leading_jets",120,-5,5);
  map_sig["h_second_jets_eta"]               = SetupHist("h_second_jets_eta_"          + name,"#eta distribution of the two second_jets",120,-5,5);
  map_sig["h_third_jets_eta"]               = SetupHist("h_third_jets_eta_"          + name,"#eta distribution of the two third_jets",120,-5,5);
  map_sig["h_fourth_jets_eta"]               = SetupHist("h_fourth_jets_eta_"          + name,"#eta distribution of the two fourth_jets",120,-5,5);
  map_sig["h_rest_jets_eta"]               = SetupHist("h_rest_jets_eta_"          + name,"#eta distribution of the two rest_jets",120,-5,5);

  /// dPhi/MT                                                                                                                                                                                                                                                                                                                                 
  map_sig["h_MTlepton"]               = SetupHist("h_MTlepton_"        + name,"Mt",100,0.0,500.0);
  map_sig["h_dphi_METlepton"]         = SetupHist("h_dphi_METlepton_"  + name , "METlepdphi", 50, -5., 5.);

  map_sig["h_contraMTlepton"]               = SetupHist("h_contraMTlepton_"        + name,"Mt",100,0.0,500.0);

  /// MET plots                                                                                                                                                                                                                                                                                                                               
  map_sig["h_PFMET"]                  = SetupHist("h_PFMET_"               + name,"Missing Et",100,0.0,500.0);
  map_sig["h_PFMET_t1"]                  = SetupHist("h_PFMET_t1_"               + name,"Missing Et",100,0.0,500.0);
  map_sig["h_PFMET2_ST"]                  = SetupHist("h_PFMET2_ST_"               + name,"Missing Et",200,0.0,100.0, "E_{T}/S^{2}_{T}");
  map_sig["h_PFMET2_STall"]                  = SetupHist("h_PFMET2_STall_"               + name,"Missing Et",200,0.0,100.0, "E_{T}/S^{2}_{T}");
  map_sig["h_PFMET_phi"]              = SetupHist("h_PFMET_phi_"           + name,"Missing Et",100,-3.2,3.2);
  map_sig["h_PFMET_phi_t1"]              = SetupHist("h_PFMET_phi_t1_"           + name,"Missing Et",100,-3.2,3.2);


  map_sig["h_nVertices"]              = SetupHist("h_nVertices_"         + name,"number of even vertices",60,0.0,60.0);
  /// Charge plot                                                                                                                                                                                                                                                                                                                             
  map_sig["h_sumcharge"]              = SetupHist("h_sumcharge_"         + name,"Charge of the lepton pair",6,-3,3);

  /// Number of objects                                                                                                                                                                                                                                                                                                                       
  map_sig["h_Nmuons"]                 = SetupHist("h_Nmuons_"           + name,"number of mu",5,0,5);

  map_sig["h_Nelectrons"]   = SetupHist("h_Nelectrons_"           + name,"number of el",5,0,5, "N_{electrons}");

  map_sig["h_Nvtx"]                   = SetupHist("h_Nvtx_"            + name,"number of vertices",40,0,40.);




}





void SignalPlotsMM::Fill(snu::KEvent ev, std::vector<snu::KMuon>& muons, std::vector<snu::KElectron>& electrons, std::vector<snu::KJet>& jets, Double_t weight) {
  std::vector<snu::KFatJet> fatjets;
  std::vector<snu::KJet> alljets;
  Fill(ev,muons, electrons, jets, alljets, fatjets,  weight);
}

void SignalPlotsMM::Fill(snu::KEvent ev, std::vector<snu::KMuon>& muons, std::vector<snu::KElectron>& electrons, std::vector<snu::KJet>& jets, std::vector<snu::KJet>& alljets,std::vector<snu::KFatJet>& fatjets, Double_t weight) {
  
  bool debug =false;
  if(debug)cout<< "Plotting [1] " << endl;
  Fill("h_Nelectrons", electrons.size(), weight);
  Fill("h_Nmuons", muons.size(), weight);
  

  //// Jet mass variables
  dijetmass_tmp=dijetmass=9999.9;
  UInt_t m(0),n(0);
  float leadawayjetratio = -999.;
  float secondawayjetratio = -999.;
  float min_ejet_Dr=10000.;  
  float min_eleadawayjet_Dr= -10000.;
  for(UInt_t i=0; i<muons.size(); i++){ 
    for(UInt_t emme=0; emme<jets.size(); emme++){
      float dR =muons[i].DeltaR(jets[emme]);
      //float dphi = fabs(TVector2::Phi_mpi_pi(muons[i].Phi()- jets[emme].Phi()));
      if(dR > 0.5) {
	if(dR > min_eleadawayjet_Dr) {
	  min_eleadawayjet_Dr=  dR;
	  if(i==0)leadawayjetratio = jets[emme].Pt() / muons[i].Pt();
	  if(i==1)secondawayjetratio = jets[emme].Pt() / muons[i].Pt();
	}	
      }
      if(dR< min_ejet_Dr) min_ejet_Dr=dR;

    }
  }
  
    
  float min_ll_Dr=10000.;
  float min_ll_DPhi=10000.;
  for(UInt_t i=0; i<muons.size(); i++){
     for(UInt_t j=0; j<muons.size(); j++){
      if(i==j) continue;
      float dPhi = fabs(TVector2::Phi_mpi_pi(muons[i].Phi() - muons[j].Phi() ));
      float dR =muons[i].DeltaR(muons[j]);
      if(dR < min_ll_Dr) min_ll_Dr = dR;
      if(dPhi < min_ll_DPhi) min_ll_DPhi = dPhi;
     }
  }

  
  float min_jj_Dr=10000.;
  float min_jj_DPhi=10000.;
  for(UInt_t i=0; i<jets.size(); i++){
    for(UInt_t j=0; j<jets.size(); j++){
      if(i==j) continue;
      float dR =jets[i].DeltaR(jets[j]);
      float dPhi = fabs(TVector2::Phi_mpi_pi(jets[i].Phi() - jets[j].Phi()));
      if(dR < min_jj_Dr) min_jj_Dr = dR;
      if(dR < min_jj_DPhi) min_jj_DPhi = dPhi;
    }
  }


  if(muons.size()!=0)Fill("h_leadLeptonstarjetdR",min_jj_Dr*min_ll_Dr,weight);
  if(muons.size()!=0)Fill("h_leadLeptondR",min_ll_Dr,weight);
  if(muons.size()!=0)Fill("h_leadLeptondPhi",min_ll_DPhi,weight);


  if(jets.size()!=0){
    if(muons.size()!=0){
      Fill("h_LeptonJetdR",min_ejet_Dr,weight);
      Fill("h_LeptonAwayJetdR",min_eleadawayjet_Dr,weight);
      Fill("h_leadLeptonAwayJetRatio",leadawayjetratio,weight);
      Fill("h_secondLeptonAwayJetRatio",secondawayjetratio,weight);
    }
    if(jets.size() > 1)Fill("h_leadJetdR",min_jj_Dr,weight);
    if(jets.size() > 1)Fill("h_leadJetdPhi",min_jj_DPhi,weight);
  }

  float min_jfj_Dr=10000.;
  float min_jfj_DPhi=10000.;
  for(UInt_t i=0; i<jets.size(); i++){
    for(UInt_t j=0; j<fatjets.size(); j++){
      float dR =jets[i].DeltaR(fatjets[j]);
      float dPhi = fabs(TVector2::Phi_mpi_pi(jets[i].Phi() - fatjets[j].Phi()));
      if(dR < min_jfj_Dr) min_jfj_Dr = dR;
      if(dR < min_jfj_DPhi) min_jfj_DPhi = dPhi;
    }
  }

  if(jets.size()!=0){
    if(fatjets.size()!=0){
      if(jets.size() > 1)Fill("h_leadJetFatJetdR",min_jj_Dr,weight);
      if(jets.size() > 1)Fill("h_leadJetFatJetdPhi",min_jj_DPhi,weight);
    }
  }


  if(debug)cout<< "Plotting [2] " << endl;  
  int nbjetl=0;
  int nbjetm=0;
  int nbjett=0;
  int nbjetl_all=0;
  int nbjetm_all=0;
  int nbjett_all=0;

  /// use CSVM https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
  float leadjetmass=0.;
  float contramass(-999.);
  int jet_eta1=0;
  int jet_eta2=0;

  for(UInt_t emme=0; emme<jets.size(); emme++){
    
    for(UInt_t enne=1; enne<jets.size(); enne++) {
      if(emme == enne) continue;
      dijetmass_tmp = (jets[emme]+jets[enne]).M();
      if(emme==0 && enne == 1) leadjetmass = dijetmass_tmp;
      if ( fabs(dijetmass_tmp-80.4) < fabs(dijetmass-80.4) ) {
	dijetmass = dijetmass_tmp;
	m = emme;
	n = enne;

	if(jets[m].Eta() < -2.1){
          jet_eta1=1;
        }
        else  if(jets[m].Eta() < -1.5){
          jet_eta1=2;
        }
        else  if(jets[m].Eta() < -0.8){
          jet_eta1=3;
	}
        else  if(jets[m].Eta() < 0.){
          jet_eta1=4;
        }
        else  if(jets[m].Eta() < 0.8){
          jet_eta1=5;
        }
	else  if(jets[m].Eta() < 1.5){
          jet_eta1=6;
        }
        else  if(jets[m].Eta() < 2.1){
          jet_eta1=7;
        }
	else  jet_eta1=8;

	if(jets[n].Eta() < -2.1){
          jet_eta2=1;
        }
        else  if(jets[n].Eta() < -1.5){
          jet_eta2=2;
        }
	else  if(jets[n].Eta() < -0.8){
          jet_eta2=3;
        }
        else  if(jets[n].Eta() < 0.){
          jet_eta2=4;
        }
        else  if(jets[n].Eta() < 0.8){
          jet_eta2=5;
        }
	else  if(jets[n].Eta() < 1.5){
          jet_eta2=6;
        }
        else  if(jets[n].Eta() < 2.1){
          jet_eta2=7;
        }
        else  jet_eta2=8;

        float dPhi = fabs(TVector2::Phi_mpi_pi(jets[emme].Phi() - jets[enne].Phi()));
        contramass=2*jets[emme].Pt()*jets[enne].Pt()*(1+cos(dPhi));
        contramass=sqrt(contramass);

      }
    }
  }
  if(jets.size() > 1)  Fill("h_contraMTlepton",contramass, weight);
  float fmass = 9999.;
  int ifn=-999;
  for(UInt_t ifj =0; ifj < fatjets.size(); ifj++){
    float fmass_tmp=fatjets[ifj].PrunedMass();
    if(fabs(fmass_tmp - 80.4) <  fabs(fmass - 80.4)){
      fmass=fmass_tmp;
      ifn=ifj;
    }
  }



  /// lljjjj plots                                                                                                                                                                                                                                                                                                                            
  if((jets.size()>3) && (muons.size() == 2)){
    int p_forward_jet(0), m_forward_jet(0);
    int index_f(-999), index_b(-999);
    vector<int> central_jets;
    float most_forward=0.;
    float most_backward=0.;



    for(unsigned int ij = 0 ; ij < jets.size(); ij++){
      if(jets[ij].Eta() > 1.5) {
        if(jets[ij].Eta() > most_forward){
          index_f= ij;
          most_forward=jets[ij].Eta();
        }
        p_forward_jet++;
        Fill("h_forward_jet_pt", jets[ij].Pt(),weight);
        Fill("h_forward_jet_eta", jets[ij].Eta(),weight);


      }
      else if(jets[ij].Eta() < -1.5) {
        if(jets[ij].Eta() < most_backward){
          index_b=ij;
          most_backward=jets[ij].Eta();
        }
	m_forward_jet++;

        Fill("h_forward_jet_pt", jets[ij].Pt(),weight);
        Fill("h_forward_jet_eta", jets[ij].Eta(),weight);
      }
      else{
        Fill("h_central_jet_pt", jets[ij].Pt(),weight);
        Fill("h_central_jet_eta", jets[ij].Eta(),weight);
      }
    }
    for(unsigned int ij = 0 ; ij < jets.size(); ij++){
      if(ij == unsigned(index_b)) continue;
      if(ij == unsigned(index_f)) continue;
      central_jets.push_back( ij);
    }

    for(unsigned int ij = 0 ; ij < jets.size(); ij++){
      if(ij == unsigned(index_b)) continue;
      if(ij == unsigned(index_f)) continue;
      central_jets.push_back( ij);
    }
    float ce_dijetmass=9999.9;
    float ce_dijetmass_tmp=9999.9;

    int c_m(0), c_n(0);
    for(UInt_t emme=0; emme<central_jets.size(); emme++){
      for(UInt_t enne=1; enne<central_jets.size(); enne++) {
        if(emme == enne) continue;
        ce_dijetmass_tmp = (jets[central_jets[emme]] + jets[central_jets[enne]]).M();
        if ( fabs(ce_dijetmass_tmp-80.4) < fabs(ce_dijetmass-80.4) ) {
          ce_dijetmass = ce_dijetmass_tmp;
          c_m = central_jets[emme];
          c_n = central_jets[enne];
        }
      }
    }

    Fill("JetEtaRegion",jet_eta1, jet_eta2,  weight);

    if( (p_forward_jet >= 1) && (m_forward_jet >= 1)) {
      
      Fill("h_l2jj_central_mass", (muons[1]+jets[c_m]+jets[c_n]).M(),weight);
      Fill("h_lljj_central_mass", (muons[0] + muons[1]+jets[c_n]+jets[c_m]).M(),weight);
      Fill("h_l1jj_central_mass", (muons[0]+jets[c_n]+jets[c_m]).M(),weight);

    }
  }
  
  /// lljjjj plots                                                                                                                                                               
  if((jets.size()>3) && (muons.size() == 2)){

    
    
    for(unsigned int ijj = 0; ijj < jets.size() ; ijj++){
      for(unsigned int ijj2 = 1; ijj2 < jets.size() ; ijj2++){
        for(unsigned int ijj3 = 2; ijj3< jets.size() ; ijj3++){
          for(unsigned int ijj4 = 3; ijj4 < jets.size() ; ijj4++){
            if(ijj == ijj2) continue;
            if(ijj == ijj3) continue;
            if(ijj == ijj4) continue;
            if(ijj2 == ijj3) continue;
            if(ijj2 == ijj4) continue;
            if(ijj3 == ijj4) continue;

            Fill("h_l1jjjj_mass", (muons[0]+ jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
            Fill("h_l2jjjj_mass", (muons[1]+jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
            Fill("h_lljjjj_mass", (muons[0] + muons[1]+jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
            if(muons[0].Charge() == muons[1].Charge()) {
              Fill("h_l1jjjj_ss_mass", (muons[0]+ jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
              Fill("h_l2jjjj_ss_mass", (muons[1]+jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
              Fill("h_lljjjj_ss_mass", (muons[0] + muons[1]+jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
            }
            else{
              Fill("h_l1jjjj_os_mass", (muons[0]+ jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
              Fill("h_l2jjjj_os_mass", (muons[1]+jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
              Fill("h_lljjjj_os_mass", (muons[0] + muons[1]+jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
            }
          }
        }
      }
    }
  }
  
  float lljjmass=0.;
  float lmeanjjmass=0.;
  float ldiffmeanjjmass=0.;

  if(muons.size() ==2 ){
    Fill("h_llmass", (muons[0]+muons[1]).M(),weight);
    Fill("h_llpt", (muons[0]+muons[1]).Pt(),weight);
    if(jets.size() >= 2){
      Fill("h_l1jjmass", (muons[0] + jets[m] + jets[n]).M(),weight);
      Fill("h_l2jjmass", (muons[1] + jets[m] + jets[n]).M(),weight);
      Fill("h_l1l2jjmass", (muons[0] + jets[m] + jets[n]).M(), (muons[1] + jets[m] + jets[n]).M(),weight);
      
      float l1_l2=fabs( (muons[0] + jets[m] + jets[n]).M() -  (muons[1] + jets[m] + jets[n]).M());
      Fill("h_l1_l2jjmass",  l1_l2, weight);
      lljjmass = (muons[0] +muons[1] + jets[m] + jets[n]).M();
      lmeanjjmass = ((muons[0] + jets[m] + jets[n]).M() +  (muons[1] + jets[m] + jets[n]).M()) / 2.;
      ldiffmeanjjmass = ((muons[0] + jets[m] + jets[n]).M() -  (muons[1] + jets[m] + jets[n]).M()) / 2.;
      Fill("h_lljjmass", (muons[0] +muons[1] + jets[m] + jets[n]).M(),weight);
    }
    if(fatjets.size() > 0){
      Fill("h_l1fjmass", (muons[0] + fatjets[ifn]).M(),weight);
      Fill("h_l2fjmass", (muons[1] + fatjets[ifn]).M(),weight);
      Fill("h_llfjmass", (muons[0] +muons[1] + fatjets[ifn]).M(),weight);
    }

  }

  
  if(leadjetmass != 0.)Fill("h_leadjetmass", leadjetmass,weight);
  
  if(debug)cout<< "Plotting [3] " << endl;  

   
  int sum_charge(0);
  float osmass_closeZ(100000.);
  int nos(0.);
  if(muons.size()>=2){
    for(unsigned int imu1 = 0 ; imu1< muons.size()-1; imu1++){
      for(unsigned int imu2 = imu1+1 ; imu2< muons.size(); imu2++){
        if(muons[imu1].Charge() != muons[imu2].Charge())  {
          nos++;
          if( fabs( (muons[imu1]+muons[imu2]).M() - 91.1) < fabs(osmass_closeZ - 91.1) ){
            osmass_closeZ = (muons[imu1]+muons[imu2]).M();
	  }
          Fill("h_osllmass", (muons[imu1]+muons[imu2]).M(),weight);
	}
      }
    }

    Fill("h_osll_closestZ_mass",osmass_closeZ ,weight);

    if(muons.size()==3){
      snu::KParticle eee= muons[0] + muons[1]+ muons[2];

      Fill("h_lllmass",  eee.M(), weight);
      Fill("h_lllpt",  eee.Pt(), weight);
    }

    if(muons.size() >= 4){

      snu::KParticle eeee= muons[0] + muons[1]+ muons[2]+muons[3];

      Fill("h_llllmass",  eeee.M(), weight);
      Fill("h_llllpt",  eeee.Pt(), weight);

      std::vector<int> vplus;
      std::vector<int> vminus;
      for(unsigned int imu1 = 0 ; imu1< muons.size(); imu1++){
	if(muons[imu1].Charge() < 0 ) vplus.push_back(imu1);
	else vminus.push_back(imu1);
      }
      if ( (vplus.size() == 2) || (vminus.size() == 2) ) {

	snu::KParticle Z1cand;
	snu::KParticle Z2cand;
	if(1){
	  snu::KParticle Z1= muons[vplus[0]] + muons[vminus[0]];
	  snu::KParticle Z2= muons[vplus[1]] + muons[vminus[1]];
	  if(Z1.M() > 4. && Z1.M() < 120) {
	    if(Z1.M() > 40. && Z1.M() < 120) {
	      if(fabs(Z1.M() -90.1) < fabs(Z2.M() -90.1)){
		Z1cand=Z1;
		Z2cand=Z2;
	      }
	      else{
		Z1cand=Z2;
		Z2cand=Z1;
	      }
	    }
	  }
	  if(fabs(Z1.M() - 91.19) <10.) {
	    if(fabs(Z2.M() - 91.19) <10.) {
	      Fill("h_zzmass", Z1.M(),weight);
	      Fill("h_zzmass", Z2.M(),weight);
	    }
	  }
	}
	if(1){
	  snu::KParticle Z1= muons[vplus[0]] + muons[vminus[1]];
	  snu::KParticle Z2= muons[vplus[1]] + muons[vminus[0]];

	  if(Z1.M() > 4. && Z1.M() < 120) {
	    if(Z1.M() > 40. && Z1.M() < 120) {
	      if(fabs(Z1.M() -90.1) < fabs(Z2.M() -90.1)){
		if(fabs(Z1.M() -90.1) < fabs(Z1cand.M() -90.1)){
		  Z1cand=Z1;
		  Z2cand=Z2;
		}
	      }
	      else{
		if(fabs(Z2.M() -90.1) < fabs(Z1cand.M() -90.1)){
		  Z1cand=Z2;
		  Z2cand=Z1;
		}
	      }
	    }
	  }

	  if(fabs(Z1.M() - 91.19) <10.) {
	    if(fabs(Z2.M() - 91.19) <10.) {
	      Fill("h_zzmass", Z1.M(),weight);
	      Fill("h_zzmass", Z2.M(),weight);
	    }
	  }
	}

	Fill("h_z1mass", Z1cand.M(),weight);
	Fill("h_z2mass", Z2cand.M(),weight);
	float ZZdphi =fabs(TVector2::Phi_mpi_pi( Z1cand.Phi() - Z2cand.Phi()));
	Fill("h_z1_z2_dphi", ZZdphi,weight);

      }

    }
    Fill("h_sumcharge", (muons[0] +muons[1]).Charge(), weight);
  }

  int imu(0);

   
   

   
  if(muons.size()==2){
    
    if( jets.size()>=1) Fill("h_lljmass", (muons[0] + muons[1]+jets[0]).M(),weight);
    if(jets.size()>1){
      Fill("h_l2jjmass", (muons[1]+jets[m]+jets[n]).M(),weight);
      Fill("h_lljjmass", (muons[0] + muons[1]+jets[m]+jets[n]).M(),weight);
      Fill("h_l1jjmass", (muons[0]+jets[m]+jets[n]).M(),weight);

    }
  }
  if(muons.size()>=2){
      if(muons[0].Charge() != muons[1].Charge())   Fill("h_osllmass", (muons[0]+muons[1]).M(),weight);
  }
  
  
  int iel(0);
  
  for(std::vector<snu::KMuon>::iterator muit = muons.begin(); muit != muons.end(); muit++, iel++){
  
    Fill("h_LeptonPt", muit->Pt(),weight);
    Fill("h_LeptonE", muit->Energy(),weight);
    Fill("h_LeptonPhi",muit->Phi(),weight);
    Fill("h_LeptonEta",muit->Eta(),weight);
    Fill("h_LeptonDXY", muit->dXY(),weight);
    Fill("h_LeptonDXYSig", muit->dXYSig() ,weight);


    if(fabs(muit->Eta()) > 1.5) {
      Fill("h_LeptonDXY_EC", muit->dXY(),weight);
      Fill("h_LeptonDZ_EC", muit->dZ(),weight);
    }
    else{
      Fill("h_LeptonDXY_B", muit->dXY(),weight);
      Fill("h_LeptonDZ_B", muit->dZ(),weight);
    }
   
    float mu_reliso_04 =  muit->RelIso04();
    float mu_minireliso = muit->RelMiniIso();
    
    Fill("h_LeptonIso", muit->RelIso04(),weight);

    Fill("h_LeptonRelIso", mu_reliso_04,weight);
    if(iel==2){
      Fill("h_thirdLeptonPt", muit->Pt(),weight);
      Fill("h_thirdLeptonRelIso", mu_reliso_04,weight);
      Fill("h_thirdLeptonMiniRelIso", mu_minireliso,weight);
      Fill("h_thirdLeptonEta",muit->Eta(),weight);
    }
    if(iel==3){
      Fill("h_fourthLeptonPt", muit->Pt(),weight);
      Fill("h_fourthLeptonRelIso", mu_reliso_04,weight);
      Fill("h_fourthLeptonMiniRelIso", mu_minireliso,weight);
      Fill("h_fourthLeptonEta",muit->Eta(),weight);
    }

    if(iel==1){
      Fill("h_secondLeptonPt", muit->Pt(),weight);

      if(fabs(muit->Eta()) > 1.5) {
	Fill("h_secondLeptonRelIso_EC", mu_reliso_04,weight);
      }
      else       Fill("h_secondLeptonRelIso_B", mu_reliso_04,weight);
      Fill("h_secondLeptonEta",muit->Eta(),weight);
      Fill("h_secondLeptonMiniRelIso", mu_minireliso,weight);
    }
    if(iel==0){
      Fill("h_leadingLeptonEta",muit->Eta(),weight);
      if(fabs(muit->Eta()) > 1.5) {
	Fill("h_leadingLeptonRelIso_EC", mu_reliso_04,weight);
      }
      else         Fill("h_leadingLeptonRelIso_B", mu_reliso_04,weight);
      Fill("h_leadingLeptonPt", muit->Pt(),weight);

    }
    
    sum_charge += muit->Charge();
  }
  
  Fill("h_sumcharge",sum_charge,weight);
  
  //// Fillplots
  for(unsigned int i=0 ; i < muons.size(); i++){
    float dphi = fabs(TVector2::Phi_mpi_pi(muons.at(i).Phi()- ev.METPhi(snu::KEvent::pfmet)));
    float MT = sqrt(2.* muons.at(i).Et()*ev.PFMET() * (1 - cos( dphi)));
    Fill("h_MTlepton",MT, weight);
    Fill("h_dphi_METlepton",dphi, weight);
  }

  if(debug)cout<< "Plotting [5] " << endl;

  Fill("h_PFMET",ev.PFMET(), weight);
  Fill("h_PFMET_phi",ev.METPhi(snu::KEvent::pfmet), weight);
  Fill("h_PFMET_t1",ev.PFMETType1(), weight);
  Fill("h_PFMET_phi_t1",ev.PFMETType1Phi(), weight);

  Fill("h_nVertices", ev.nVertices(), weight); 
  

  
  Fill("h_Nvtx",ev.nVertices(), weight);
  
  //// Mass plots
  if(jets.size()>1){
    Fill("h_jjmass", (jets[m]+jets[n]).M(),weight); 
  }

  if(fatjets.size()> 0){
    Fill("h_fjmass", (fatjets[ifn]).M(),weight);
    Fill("h_fjprunedmass", (fatjets[ifn]).PrunedMass(),weight);

  }
  
  float ht=0.;
  float htall=0.;

  if(fatjets.size()> 0){
    Fill("h_FatJetPt", fatjets[ifn].Pt(),weight);
    Fill("h_FatJetEta", fatjets[ifn].Eta(),weight);
    Fill("h_FatJetPhi", fatjets[ifn].Phi(),weight);
    for(UInt_t j=0; j < fatjets.size(); j++){
      ht+= fatjets[j].Pt();
    }
  }

  vector<snu::KJet> central_jets;
  vector<snu::KJet> forward_jets;

  float max_deta=0.;
  snu::KParticle fjets;
  
  for(UInt_t j=0; j < alljets.size(); j++){ 
    htall+= alljets[j].Pt();
  }

  for(UInt_t j=0; j < jets.size(); j++){ 

    if(fabs(jets[j].Eta()) < 1.5) central_jets.push_back(jets[j]);
    else forward_jets.push_back(jets[j]);


    ht+= jets[j].Pt();
    if(j==0)Fill("h_leadingJetPt", jets[0].Pt(),weight);
    if(j==1)Fill("h_secondJetPt", jets[1].Pt(),weight);
    if(j==0)     Fill("h_leading_jets_eta",jets[j].Eta(),weight);
    if(j==1)     Fill("h_second_jets_eta",jets[j].Eta(),weight);
    if(j==2)     Fill("h_third_jets_eta",jets[j].Eta(),weight);
    if(j==3)     Fill("h_fourth_jets_eta",jets[j].Eta(),weight);
    Fill("h_jets_pt", jets[j].Pt(),weight);
    Fill("h_jets_eta",jets[j].Eta(),weight);
    Fill("h_PileupJetIDMVA", jets[j].PileupJetIDMVA(),weight);
    Fill("h_jets_phi",jets[j].Phi(),weight);
    Fill("h_bTag",jets[j].CSVInclV2(),weight);
    if(jets.at(j).IsBTagged(KJet::CSVv2, KJet::Loose)) nbjetl++;
    if(jets.at(j).IsBTagged(KJet::CSVv2, KJet::Medium)) nbjetm++;
    if(jets.at(j).IsBTagged(KJet::CSVv2, KJet::Tight)) nbjett++;
  }
  
  for(UInt_t j=0; j < alljets.size(); j++){
    if(alljets.at(j).IsBTagged(KJet::CSVv2, KJet::Loose)) nbjetl_all++;
    if(alljets.at(j).IsBTagged(KJet::CSVv2, KJet::Medium)) nbjetm_all++;
    if(alljets.at(j).IsBTagged(KJet::CSVv2, KJet::Tight)) nbjett_all++;

  }

  float st_foward=0.;
  float st_central=0.;
  snu::KParticle fj;
  snu::KParticle cj;
  float ht_foward=0.;
  float ht_central=0.;
  for(UInt_t j=0; j < forward_jets.size(); j++){
    st_foward=st_foward+forward_jets[j].Pt();
    ht_foward=ht_foward+forward_jets[j].Pt();
    fj+=forward_jets[j];
  }
  for(UInt_t j=0; j < central_jets.size(); j++){
    st_central=st_central+central_jets[j].Pt();
    ht_central=ht_central+central_jets[j].Pt();
    cj+=central_jets[j];
  }


  float lt=0.;
  float ltc=-0.;  float ltf=0.;
  for(unsigned int iel2=0 ; iel2 < muons.size(); iel2++){
    lt+=muons[iel2].Pt();
    if(fabs(muons[iel2].Eta()) < 1.5){
      st_central=st_central+ muons[iel2].Pt();
      ltc+= muons[iel2].Pt();
      cj+=  muons[iel2];

    }
    else{
      st_foward=st_foward+ muons[iel2].Pt();
      ltf+=muons[iel2].Pt();
      fj+=muons[iel2];
    }
  }  int eta1=0.;
  int eta2=0.;
  float sumeta=0.;
  for(unsigned int iel2=0 ; iel2 < muons.size(); iel2++){
    sumeta=sumeta+muons[iel2].Eta();
    int etareg=0.;
    if(muons[iel2].Eta() < -2.1){
      etareg=1;
    }
    else  if(muons[iel2].Eta() < -1.5){
      etareg=2;
    }
    else  if(muons[iel2].Eta() < -0.8){
      etareg=3;
    }
    else  if(muons[iel2].Eta() < 0.){
      etareg=4;
    }
    else  if(muons[iel2].Eta() < 0.8){
      etareg=5;
    }
    else  if(muons[iel2].Eta() < 1.5){
      etareg=6;
    }
    else  if(muons[iel2].Eta() < 2.1){
      etareg=7;
    }
    else etareg=8;

    if(iel2==0) eta1=etareg;
    else eta2=etareg;

  }
  Fill("EtaRegion",eta1, eta2,  weight);
  Fill("h_lt",lt,weight);
  Fill("h_lt_ht",lt/ht,weight);
  if(ltf > 0.)Fill("h_lt_forward",ltf,weight);
  if(ltc > 0.) {
    Fill("h_lt_central",ltc,weight);
    Fill("h_lt_rforward",ltf/ltc,weight);
  }

  Fill("h_lt_jjdr",lt,min_jj_Dr, weight);
  Fill("h_lt_CM",lt,contramass, weight);
  Fill("h_lt_lljjmass",lt,lljjmass,weight);
  Fill("h_lt_lmeanjjmass",lt,lmeanjjmass,weight);
  Fill("h_lt_ldiffmeanjjmass",lt,ldiffmeanjjmass,weight);
  Fill("h_sumeta",sumeta,weight);

  if(st_foward>0.)Fill("h_st_forward",st_foward,weight);
  if(st_central>0.)Fill("h_st_central",st_central,weight);
  if(st_central > 0.)Fill("h_st_fc_ratio",st_foward/st_central,weight);
  if(ht_foward>0.)Fill("h_ht_forward",ht_foward,weight);
  if(ht_central>0.)Fill("h_ht_central",ht_central,weight);
  if(ht_central>0.)Fill("h_ht_fc_ratio",ht_foward/ht_central,weight);
  Fill("h_mass_forward",fj.M(),weight);
  Fill("h_mass_central",cj.M(),weight);
  Fill("h_centralNJets",central_jets.size(), weight);
  Fill("h_forwardNJets",forward_jets.size(), weight);
  Fill("h_max_deta",max_deta,weight);


  if(muons.size() ==2)
    {
      fjets=fjets+muons[0]+muons[1];
      Fill("h_llcj_mass",cj.M(), weight);
      Fill("h_llfj_mass",fj.M(), weight);
      Fill("h_llDfj_mass",fjets.M(), weight);
    }
    

  float st = ht + ev.PFMET();
  float stall = htall + ev.PFMET();
  for(unsigned int i=0 ; i < muons.size(); i++){
    st+= muons.at(i).Pt();
    stall+= muons.at(i).Pt();
  }
  
  Fill("h_HT", ht,weight);
  Fill("h_ST", st,weight);
  Fill("h_STall", stall,weight);
  Fill("h_PFMET2_ST",pow(ev.PFMET(),2.)/ st, weight);
  Fill("h_PFMET2_STall",pow(ev.PFMET(),2.)/ stall, weight);
  Fill("h_Nbjets_l",nbjetl, weight);
  Fill("h_Nbjets_m",nbjetm, weight);
  Fill("h_Nbjets_t",nbjett, weight);
  Fill("h_Nbjets_all_l",nbjetl_all, weight);
  Fill("h_Nbjets_all_m",nbjetm_all, weight);
  Fill("h_Nbjets_all_t",nbjett_all, weight);

  Fill("h_Njets",jets.size(), weight);
  Fill("h_Nfatjets",fatjets.size(), weight);

  
  if(debug)cout<< "Plotting [6] " << endl;
  return;
}/// End of Fill



void SignalPlotsMM::Write() {
 
  for(map<TString, TH1*>::iterator it = map_sig.begin(); it != map_sig.end(); it++){
    it->second->Write();
  }

  for(map<TString, TH2*>::iterator it = map_sig2.begin(); it != map_sig2.end(); it++){
    it->second->Write();
  }


  for(map<TString, TH3*>::iterator it = map_sig3.begin(); it != map_sig3.end(); it++){
    it->second->Write();
  }

}


SignalPlotsMM::SignalPlotsMM(): StdPlots(){
}


/**
 * Copy constructor.
 */
SignalPlotsMM::SignalPlotsMM(const SignalPlotsMM& sp): StdPlots(sp)
{
  for(std::map<TString, TH1*>::iterator mit = map_sig.begin(); mit != map_sig.end() ; mit++){
    std::map<TString, TH1*>::iterator mit2 = sp.GetMap().find(mit->first);
    mit->second = mit2->second;
  }

  for(std::map<TString, TH2*>::iterator mit = map_sig2.begin(); mit != map_sig2.end() ; mit++){
    std::map<TString, TH2*>::iterator mit2 = sp.GetMap2().find(mit->first);
    mit->second = mit2->second;
  }

  for(std::map<TString, TH3*>::iterator mit = map_sig3.begin(); mit != map_sig3.end() ; mit++){
    std::map<TString, TH3*>::iterator mit2 = sp.GetMap3().find(mit->first);
    mit->second = mit2->second;
  }
  
}


SignalPlotsMM& SignalPlotsMM::operator= (const SignalPlotsMM& sp)
{
  if (this != &sp) {

    for(std::map<TString, TH1*>::iterator mit = map_sig.begin(); mit != map_sig.end() ; mit++){
      std::map<TString, TH1*>::iterator mit2 = sp.GetMap().find(mit->first);
      mit->second = mit2->second;
    }
    
    for(std::map<TString, TH2*>::iterator mit = map_sig2.begin(); mit != map_sig2.end() ; mit++){
      std::map<TString, TH2*>::iterator mit2 = sp.GetMap2().find(mit->first);
      mit->second = mit2->second;
    }

    for(std::map<TString, TH3*>::iterator mit = map_sig3.begin(); mit != map_sig3.end() ; mit++){
      std::map<TString, TH3*>::iterator mit2 = sp.GetMap3().find(mit->first);
      mit->second = mit2->second;
    }
  }
  return *this;
}

SignalPlotsMM::~SignalPlotsMM() {
   for(std::map<TString, TH1*>::iterator mit = map_sig.begin(); mit != map_sig.end() ; mit++){
     delete mit->second ;
  }

   for(std::map<TString, TH2*>::iterator mit = map_sig2.begin(); mit != map_sig2.end() ; mit++){
     delete mit->second ;
   }

   for(std::map<TString, TH3*>::iterator mit = map_sig3.begin(); mit != map_sig3.end() ; mit++){
     delete mit->second ;
   }
   
}


TH1D* SignalPlotsMM::SetupHist(TString hname, TString alabel, int nbin, double xmin, double xmax, TString xtitle ){

  TH1D* htmp= new TH1D(hname,alabel,nbin,xmin,xmax);
  htmp->GetXaxis()->SetTitle(xtitle);
  return htmp;
}


TH2D* SignalPlotsMM::SetupHist2(TString hname, TString alabel, int nbinx, double xmin, double xmax, int nbiny, double ymin, double ymax, TString xtitle , TString ytitle){

  TH2D* htmp= new TH2D(hname,alabel,nbinx,xmin,xmax, nbiny, ymin, ymax);
  htmp->GetXaxis()->SetTitle(xtitle);
  htmp->GetYaxis()->SetTitle(ytitle);
  return htmp;
}



void SignalPlotsMM::Fill(TString name, double value, double w){
  std::map<TString, TH1*>::iterator it = map_sig.find(name);
  if(it!= map_sig.end())   it->second->Fill(value, w);

  else cout << name << " not found in map_sig" << endl;
  return;
}
 
void SignalPlotsMM::Fill(TString name, double value1, double value2, double w){
   std::map<TString, TH2*>::iterator it = map_sig2.find(name);
   if(it!= map_sig2.end()) it->second->Fill(value1, value2, w);
   else cout << name << " not found in map_sig" << endl;
   return;
 }



void SignalPlotsMM::Fill(TString name, double value1, double value2, double value3, double w){
  std::map<TString, TH3*>::iterator it = map_sig3.find(name);
  if(it!= map_sig3.end()) it->second->Fill(value1, value2, value3, w);
  else cout << name << " not found in map_sig" << endl;
  return;
}




