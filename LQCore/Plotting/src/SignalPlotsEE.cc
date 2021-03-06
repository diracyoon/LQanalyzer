#include "SignalPlotsEE.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

SignalPlotsEE::SignalPlotsEE(TString name): StdPlots(name){

  TH1::SetDefaultSumw2(true);

 
  /// Mass plots
  map_sig["h_jjmass"]                 =     new TH1D("h_dijetsmass_"        + name,"Invariant mass of the two leading jets",500,0.,500.);
  map_sig["h_leadjetmass"]            =     new TH1D("h_leadjetsmass_"      + name,"Invariant mass of the two leading jets",150,0.,750.);
  map_sig["h_osllmass"]               =     new TH1D("h_osllmass_"          + name,"Invariant mass of the two leading os electrons",100,0,500);
  map_sig["h_l1jjmass"]               =     new TH1D("h_l1jjmass_"          + name,"Invariant mass of the two leading jets and leading electron",100,0,1000);
  map_sig["h_l2jjmass"]               =     new TH1D("h_l2jjmass_"          + name,"Invariant mass of the two leading jets and second electron",100,0,1000);
  map_sig["h_llmass"]                 =     new TH1D("h_llmass_"           + name,"Invariant mass of the two leading electrons",1000,0,1000);
  map_sig["h_llpt"]                 =     new TH1D("h_llpt_"           + name,"Invariant pt of the two leading electrons",500,0,1000);
  map_sig["h_lljjmass"]               =     new TH1D("h_lljjmass_"         + name,"Invariant mass of the four particles",200,0,2000);
  map_sig["h_lljjjj_mass"]             =     new TH1D("h_lljjjj_mass_"       + name,"Invariant mass of the four particles",200,0,2000);
  map_sig["h_l1jjjj_mass"]             =     new TH1D("h_l1jjjj_mass_"       + name,"Invariant mass of the four particles",200,0,2000);
  map_sig["h_l2jjjj_mass"]             =     new TH1D("h_l2jjjj_mass_"       + name,"Invariant mass of the four particles",200,0,2000);
  map_sig["h_lljjjj_os_mass"]          =     new TH1D("h_lljjjj_os_mass_"       + name,"Invariant mass of the four particles",200,0,2000);
  map_sig["h_l1jjjj_os_mass"]          =     new TH1D("h_l1jjjj_os_mass_"       + name,"Invariant mass of the four particles",200,0,2000);
  map_sig["h_l2jjjj_os_mass"]          =     new TH1D("h_l2jjjj_os_mass_"       + name,"Invariant mass of the four particles",200,0,2000);
  map_sig["h_lljjjj_ss_mass"]          =     new TH1D("h_lljjjj_ss_mass_"       + name,"Invariant mass of the four particles",200,0,2000);
  map_sig["h_l1jjjj_ss_mass"]          =     new TH1D("h_l1jjjj_ss_mass_"       + name,"Invariant mass of the four particles",200,0,2000);
  map_sig["h_l2jjjj_ss_mass"]          =     new TH1D("h_l2jjjj_ss_mass_"       + name,"Invariant mass of the four particles",200,0,2000);


  map_sig["h_lljmass"]                =     new TH1D("h_lljmass_"         + name,"Invariant mass of the four particles",200,0,2000);
  
  /// DiLepton plots
  map_sig["h_ll_phi"]                 =     new TH1D("h_ll_phi_"          + name,"dilepton pt",100,-4,4);
  map_sig["h_ll_eta"]                 =     new TH1D("h_ll_eta_"          + name,"dilepton pt",100,-5,5);
  map_sig["h_ll_pt"]                  =     new TH1D("h_ll_pt_"          + name,"dilepton pt",100,0,1000);

  /// 4J plots
  map_sig["h_l1jj_central_mass"]               =     new TH1D("h_l1jj_centralmass_"          + name,"Invariant mass of the two leading jets and leading electron",100,0,1000);
  map_sig["h_l2jj_central_mass"]               =     new TH1D("h_l2jj_centralmass_"          + name,"Invariant mass of the two leading jets and leading electron",100,0,1000);
  map_sig["h_lljj_central_mass"]               =     new TH1D("h_lljj_centralmass_"          + name,"Invariant mass of the two leading jets and leading electron",100,0,1000);

  map_sig["h_forward_jet_pt"]             =     new TH1D("h_forward_jet_pt_"                + name,"h_forward_jet_pt",60,0,300);
  map_sig["h_central_jet_pt"]             =     new TH1D("h_central_jet_pt_"                + name,"h_central_jet_pt",60,0,300);
  map_sig["h_forward_jet_eta"]             =     new TH1D("h_forward_jet_eta_"                + name,"h_forward_jet_eta",50,-2.5,2.5);
  map_sig["h_central_jet_eta"]             =     new TH1D("h_central_jet_eta_"                + name,"h_central_jet_eta",50,-2.5,2.5);


  
  /// Lepton plots  
  map_sig["h_LeptonEta"]              =     new TH1D("h_LeptonEta_"         + name,"leading lepton eta",60,-3.,3.);
  map_sig["h_LeptonPhi"]              =     new TH1D("h_LeptonPhi_"         + name,"leading lepton phi",64,-3.2,3.2);
  map_sig["h_LeptonPt"]               =     new TH1D("h_LeptonPt_"          + name,"lepton pt",100,0,500);
  map_sig["h_leadingLeptonPt"]        =     new TH1D("h_leadingLeptonPt_"   + name,"leading lepton pt",100,0,500);
  map_sig["h_secondLeptonPt"]         =     new TH1D("h_secondLeptonPt_"    + name,"secondary lepton pt",60,0,300);
  map_sig["h_leadingLeptonEta"]       =     new TH1D("h_leadingLeptonEta_"  + name,"leading lepton eta",60,-3.,3.);
  map_sig["h_secondLeptonEta"]        =     new TH1D("h_secondLeptonEta_"   + name,"second lepton eta",60,-3.,3.);
  map_sig["h_LeptonIso"]              =     new TH1D("h_LeptonIso_"         + name,"leading lepton relIso",50,0,10.);
  map_sig["h_LeptonRelIso"]           =     new TH1D("h_LeptonRelIso_"      + name,"leading lepton relIso",100,0,1.);
  map_sig["h_leadingLeptonRelIso"]    =     new TH1D("h_leadingLeptonRelIso_"      + name,"leading lepton relIso",100,0,1.);
  map_sig["h_secondLeptonRelIso"]     =     new TH1D("h_secondLeptonRelIso_"      + name,"leading lepton relIso",100,0,1.);
  map_sig["h_LeptonDZ"]               =     new TH1D("h_LeptonDZ_"          + name," leading lepton DZ", 400, -0.5 , 0.5);
  map_sig["h_LeptonDXY"]              =     new TH1D("h_LeptonDXY_"         + name," leading lepton DXY", 1000, -0.25 , 0.25);
  
  
  /// Lepton/Jet/Lepton dR
  map_sig["h_LeptonJetdR"]            =     new TH1D("h_LeptonJetdR_"       + name,"lepton jet dR",50,0,5);
  map_sig["h_LeptonAwayJetdR"]        =     new TH1D("h_LeptonAwayJetdR_"   + name, "away jet dr", 50,0,5);
  map_sig["h_leadLeptondR"]           =     new TH1D("h_leadLeptondR_"      + name,"leading lepton dR",50,0,5);
  map_sig["h_leadLeptondPhi"]         =     new TH1D("h_leadLeptondPhi_"    + name,"leading lepton dPhi",50,-5,5);
  map_sig["h_leadJetdR"]                =     new TH1D("h_leadJetdR_"         + name,"leading jet dR",50,0,5);
  map_sig["h_leadJetdPhi"]              =     new TH1D("h_leadJetdPhi_"       + name,"leading jet dPhi",50,-5,5);

  //// Jet Plots
  map_sig["h_leadingJetPt"]           =     new TH1D("h_leadingJetPt_"      + name,"leading jet pt",60,0,300);
  map_sig["h_secondJetPt"]            =     new TH1D("h_secondJetPt_"       + name,"secondary jet pt",60,0,300);
  map_sig["h_HT"]                     =     new TH1D("h_HT_"                + name,"sum jet pt",50,0,1000);
  map_sig["h_ST"]                     =     new TH1D("h_ST_"                + name,"sum event pt",50,0,1000);
  map_sig["h_jets_pt"]                =     new TH1D("h_jets_pt_"           + name,"jet pt",60,0,300);
  map_sig["h_jets_eta"]               =     new TH1D("h_jets_eta_"          + name,"#eta distribution of the two jets",120,-3,3);
  map_sig["h_jets_phi"]               =     new TH1D("h_jets_phi_"          + name,"#phi distribution of the two jets",140,-3.5,3.5);
  map_sig["h_PileupJetIDMVA"]         =     new TH1D("h_pileupJetIDMVA_"    + name, "" ,100, -1.,1.);
  map_sig["h_bTag"]                   =     new TH1D("h_bTag_"              + name,"bTag discrimant",100,-1,3);
  map_sig["h_Njets"]                  =     new TH1D("h_Njets_"             + name,"number of jets",10,0,10);
  map_sig["h_Nbjets"]                 =     new TH1D("h_Nbjets_"            + name,"number of b jets",5,0,5);
  map_sig["h_lep_jet_dphi"]           =      new TH1D("h_lep_jet_dphi_"            + name,"dphi lepton jets", 50, 0., 5.);
  map_sig["h_awayjet_chargedem_frac"]         =      new TH1D("h_awayjet_chargedem_frac_"   + name,"charged em frac", 50, 0., 1.);

  map_sig["h_leading_jets_eta"]               =     new TH1D("h_leading_jets_eta_"          + name,"#eta distribution of the two leading_jets",120,-5,5);
  map_sig["h_second_jets_eta"]               =     new TH1D("h_second_jets_eta_"          + name,"#eta distribution of the two second_jets",120,-5,5);
  map_sig["h_third_jets_eta"]               =     new TH1D("h_third_jets_eta_"          + name,"#eta distribution of the two third_jets",120,-5,5);
  map_sig["h_fourth_jets_eta"]               =     new TH1D("h_fourth_jets_eta_"          + name,"#eta distribution of the two fourth_jets",120,-5,5);
  map_sig["h_rest_jets_eta"]               =     new TH1D("h_rest_jets_eta_"          + name,"#eta distribution of the two rest_jets",120,-5,5);

  /// dPhi/MT
  map_sig["h_MTlepton"]               =     new TH1D("h_MTlepton_"        + name,"Mt",100,0.0,500.0);
  map_sig["h_dphi_METlepton"]         =     new TH1D("h_dphi_METlepton_"  + name , "METlepdphi", 50, -5., 5.);
 
  
  /// MET plots
  map_sig["h_PFMET"]                  =     new TH1D("h_PFMET_"               + name,"Missing Et",100,0.0,500.0);
  map_sig["h_PFMET_phi"]              =     new TH1D("h_PFMET_phi_"           + name,"Missing Et",100,-3.2,3.2);


  map_sig["h_nVertices"]              =     new TH1D("h_nVertices_"         + name,"number of even vertices",60,0.0,60.0);
  /// Charge plot
  map_sig["h_sumcharge"]              =     new TH1D("h_sumcharge_"         + name,"Charge of the lepton pair",6,-3,3);

  /// Number of objects
  map_sig["h_Nmuons"]                 =     new TH1D("h_Nmuons_"           + name,"number of mu",5,0,5);
  map_sig["h_Nelectrons"]             =     new TH1D("h_Nelectrons_"           + name,"number of el",5,0,5);


  map_sig["h_Nvtx"]                   =     new TH1D("h_Nvtx_"            + name,"number of vertices",40,0,40.);


}




void SignalPlotsEE::Fill(snu::KEvent ev, std::vector<snu::KMuon>& muons, std::vector<snu::KElectron>& electrons, std::vector<snu::KJet>& jets, Double_t weight) {
  
  bool debug =false;
  if(debug)cout<< "Plotting [1] " << endl;
  Fill("h_Nelectrons", electrons.size(), weight);
  Fill("h_Nmuons", muons.size(), weight);
  
  

  //// Jet mass variables
  dijetmass_tmp=dijetmass=9999.9;
  UInt_t m(0),n(0);
  
  float min_ejet_Dr=10000.;  
  float min_eleadawayjet_Dr= -10000.;
  for(UInt_t i=0; i<electrons.size(); i++){ 
    for(UInt_t emme=0; emme<jets.size(); emme++){
      float dR =electrons[i].DeltaR(jets[emme]);
      //float dphi = fabs(TVector2::Phi_mpi_pi(electrons[i].Phi()- jets[emme].Phi()));
      if(dR > 0.5) {
	if(dR > min_eleadawayjet_Dr) min_eleadawayjet_Dr=  dR;
      }
      if(dR< min_ejet_Dr) min_ejet_Dr=dR;

    }
  }
  
  

    
  float min_ee_Dr=10000.;
  float min_ee_DPhi=10000.;
  for(UInt_t i=0; i<electrons.size(); i++){
     for(UInt_t j=0; j<electrons.size(); j++){
      if(i==j) continue;
      float dPhi = fabs(TVector2::Phi_mpi_pi(electrons[i].Phi() - electrons[j].Phi() ));
      float dR =electrons[i].DeltaR(electrons[j]);
      if(dR < min_ee_Dr) min_ee_Dr = dR;
      if(dPhi < min_ee_DPhi) min_ee_DPhi = dPhi;
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

  
  if(electrons.size()!=0)Fill("h_leadLeptondR",min_ee_Dr,weight);
  if(electrons.size()!=0)Fill("h_leadLeptondPhi",min_ee_DPhi,weight);


  if(jets.size()!=0){
    if(electrons.size()!=0){
      Fill("h_LeptonJetdR",min_ejet_Dr,weight);
      Fill("h_LeptonAwayJetdR",min_eleadawayjet_Dr,weight);
    }
    if(jets.size() > 1)Fill("h_leadJetdR",min_jj_Dr,weight);
    if(jets.size() > 1)Fill("h_leadJetdPhi",min_jj_DPhi,weight);
  }

  if(debug)cout<< "Plotting [2] " << endl;  
  int nbjet=0;
  /// use CSVM https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
  float leadjetmass=0.;
  for(UInt_t emme=0; emme<jets.size(); emme++){
    
    for(UInt_t enne=1; enne<jets.size(); enne++) {
      if(emme == enne) continue;
      dijetmass_tmp = (jets[emme]+jets[enne]).M();
      if(emme==0 && enne == 1) leadjetmass = dijetmass_tmp;
      if ( fabs(dijetmass_tmp-80.4) < fabs(dijetmass-80.4) ) {
	dijetmass = dijetmass_tmp;
	m = emme;
	n = enne;
      }
    }
  }

  if((jets.size()>3) && (electrons.size() == 2)){
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
    
    if( (p_forward_jet >= 1) && (m_forward_jet >= 1)) {
      
      
      Fill("h_l2jj_central_mass", (electrons[1]+jets[c_m]+jets[c_n]).M(),weight);
      Fill("h_lljj_central_mass", (electrons[0] + electrons[1]+jets[c_n]+jets[c_m]).M(),weight);
      Fill("h_l1jj_central_mass", (electrons[0]+jets[c_n]+jets[c_m]).M(),weight);
      
    }
  }

  /// lljjjj plots
  if((jets.size()>3) && (electrons.size() == 2)){
    
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

	    Fill("h_l1jjjj_mass", (electrons[0]+ jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
	    Fill("h_l2jjjj_mass", (electrons[1]+jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
	    Fill("h_lljjjj_mass", (electrons[0] + electrons[1]+jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
	    if(electrons[0].Charge() == electrons[1].Charge()) {
	      Fill("h_l1jjjj_ss_mass", (electrons[0]+ jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
	      Fill("h_l2jjjj_ss_mass", (electrons[1]+jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
	      Fill("h_lljjjj_ss_mass", (electrons[0] + electrons[1]+jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
	    }
	    else{
	      Fill("h_l1jjjj_os_mass", (electrons[0]+ jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
	      Fill("h_l2jjjj_os_mass", (electrons[1]+jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
	      Fill("h_lljjjj_os_mass", (electrons[0] + electrons[1]+jets[ijj]+jets[ijj2]+jets[ijj3]+jets[ijj4]).M(),weight);
	    }
	  }
	}
      }
    }
  }

  if(electrons.size() ==2 ){
    Fill("h_llmass", (electrons[0]+electrons[1]).M(),weight);
    Fill("h_llpt", (electrons[0]+electrons[1]).Pt(),weight);
    if(jets.size() >= 2){
      Fill("h_l1jjmass", (electrons[0] + jets[m] + jets[n]).M(),weight);
      Fill("h_l2jjmass", (electrons[1] + jets[m] + jets[n]).M(),weight);
      Fill("h_lljjmass", (electrons[0] +electrons[1] + jets[m] + jets[n]).M(),weight);
    }
  }

  
  if(leadjetmass != 0.)Fill("h_leadjetmass", leadjetmass,weight);
  
  if(debug)cout<< "Plotting [3] " << endl;  

   int sum_charge(0);
   


  if(electrons.size()>=2){
      if(electrons[0].Charge() != electrons[1].Charge())   Fill("h_osllmass", (electrons[0]+electrons[1]).M(),weight);
      Fill("h_sumcharge", (electrons[0] +electrons[1]).Charge(), weight);
  }
  
  int iel(0);
 
  for(std::vector<snu::KElectron>::iterator elit = electrons.begin(); elit != electrons.end(); elit++, iel++){
  
    Fill("h_LeptonPt", elit->Pt(),weight);
    Fill("h_LeptonPhi",elit->Phi(),weight);
    Fill("h_LeptonEta",elit->Eta(),weight);
    Fill("h_LeptonDXY", elit->dxy(),weight);
    Fill("h_LeptonDZ", elit->dz(),weight);
     
    float el_reliso_03 =  elit->PFRelIso(0.3);
    float el_iso_03 = elit->PFAbsIso(0.3);
    
    
    Fill("h_LeptonIso", el_iso_03,weight);
    Fill("h_LeptonRelIso", el_reliso_03,weight);
    if(iel==1){
      Fill("h_secondLeptonPt", elit->Pt(),weight);
      Fill("h_secondLeptonRelIso", el_reliso_03,weight);
      Fill("h_secondLeptonEta",elit->Eta(),weight);
    }
    if(iel==0){
      Fill("h_leadingLeptonEta",elit->Eta(),weight);
      Fill("h_leadingLeptonRelIso", el_reliso_03,weight);
      Fill("h_leadingLeptonPt", elit->Pt(),weight);
    }
    
    sum_charge += elit->Charge();
  }
  
  Fill("h_sumcharge",sum_charge,weight);
  
  //// Fillplots
  for(unsigned int i=0 ; i < electrons.size(); i++){
    float dphi = fabs(TVector2::Phi_mpi_pi(electrons.at(i).Phi()- ev.METPhi(snu::KEvent::pfmet)));
    float MT = sqrt(2.* electrons.at(i).Et()*ev.PFMET() * (1 - cos( dphi)));
    Fill("h_MTlepton",MT, weight);
    Fill("h_dphi_METlepton",dphi, weight);
  }

  if(debug)cout<< "Plotting [5] " << endl;

  Fill("h_PFMET",ev.PFMET(), weight);
  Fill("h_PFMET_phi",ev.METPhi(snu::KEvent::pfmet), weight);
  Fill("h_nVertices", ev.nVertices(), weight);
  
  
  Fill("h_Nvtx",ev.PileUpInteractionsTrue(), weight);
  
  //// Mass plots
  if(jets.size()>1){
    Fill("h_jjmass", (jets[m]+jets[n]).M(),weight); 
  }
  float ht=0.;
  for(UInt_t j=0; j < jets.size(); j++){ 
    ht+= jets[j].Pt();
    if(j==0)Fill("h_leadingJetPt", jets[0].Pt(),weight);
    if(j==1)Fill("h_secondJetPt", jets[1].Pt(),weight);
    Fill("h_jets_pt", jets[j].Pt(),weight);
    Fill("h_jets_eta",jets[j].Eta(),weight);
    if(j==0)     Fill("h_leading_jets_eta",jets[j].Eta(),weight);
    if(j==1)     Fill("h_second_jets_eta",jets[j].Eta(),weight);
    if(j==2)     Fill("h_third_jets_eta",jets[j].Eta(),weight);
    if(j==3)     Fill("h_fourth_jets_eta",jets[j].Eta(),weight);

    Fill("h_PileupJetIDMVA", jets[j].PileupJetIDMVA(),weight);
    Fill("h_jets_phi",jets[j].Phi(),weight);
    Fill("h_bTag",jets[j].CSVInclV2(),weight);
    if(jets.at(j).IsBTagged(KJet::CSVv2, KJet::Medium)) nbjet++; 
    for(unsigned int iel2=0 ; iel2 < electrons.size(); iel2++){
      float dphi =fabs(TVector2::Phi_mpi_pi(electrons.at(iel2).Phi()- jets.at(j).Phi()));
      Fill("h_lep_jet_dphi", dphi,weight);
      if(dphi > 2.5)  Fill("h_awayjet_chargedem_frac", (jets.at(j).ChargedEMEnergyFraction()), weight);
    }
  }
  
  float st = ht + ev.PFMET();
  for(unsigned int i=0 ; i < electrons.size(); i++){
    st+= electrons.at(i).Pt();
  }
  
  Fill("h_HT", ht,weight);
  Fill("h_ST", st,weight);
  Fill("h_Nbjets",nbjet, weight);
  Fill("h_Njets",jets.size(), weight);
  
  if(debug)cout<< "Plotting [6] " << endl;
  return;
}/// End of Fill



void SignalPlotsEE::Write() {
 
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


SignalPlotsEE::SignalPlotsEE(): StdPlots(){
}


/**
 * Copy constructor.
 */
SignalPlotsEE::SignalPlotsEE(const SignalPlotsEE& sp): StdPlots(sp)
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


SignalPlotsEE& SignalPlotsEE::operator= (const SignalPlotsEE& sp)
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

SignalPlotsEE::~SignalPlotsEE() {
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

void SignalPlotsEE::Fill(TString name, double value, double w){
  std::map<TString, TH1*>::iterator it = map_sig.find(name);
  if(it!= map_sig.end())   it->second->Fill(value, w);

  else cout << name << " not found in map_sig" << endl;
  return;
}
 
void SignalPlotsEE::Fill(TString name, double value1, double value2, double w){
   std::map<TString, TH2*>::iterator it = map_sig2.find(name);
   if(it!= map_sig2.end()) it->second->Fill(value1, value2, w);
   else cout << name << " not found in map_sig" << endl;
   return;
 }



void SignalPlotsEE::Fill(TString name, double value1, double value2, double value3, double w){
  std::map<TString, TH3*>::iterator it = map_sig3.find(name);
  if(it!= map_sig3.end()) it->second->Fill(value1, value2, value3, w);
  else cout << name << " not found in map_sig" << endl;
  return;
}




