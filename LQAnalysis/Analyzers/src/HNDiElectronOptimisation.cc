// $Id: ExampleAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNDiElectronOptimisation Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *d
 ***************************************************************************/

/// Local includes
#include "HNDiElectronOptimisation.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HNDiElectronOptimisation);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
HNDiElectronOptimisation::HNDiElectronOptimisation() :  AnalyzerCore(),  out_electrons(0) {


  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNDiElectronOptimisation");

  Message("In HNDiElectronOptimisation constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void HNDiElectronOptimisation::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //

   Message("Making clever hists for Z ->ll test code", INFO);

   //// Initialise Plotting class functions
   /// MakeCleverHistograms ( type, "label")  type can be muhist/elhist/jethist/sighist



   return;
}


void HNDiElectronOptimisation::ExecuteEvents()throw( LQError ){
  
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;

  if(!isData){
    weight*= MCweight;
    weight*=WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi);
  }

  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; 


  /// Get Efficiency of signal + Bkg using multiple triggers                                                                                                                                                                                                                  


  /// Dimuon for mmu channel
  std::vector<TString> triggerlist;
  triggerlist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  float trigger_ps_weight= WeightByTrigger(triggerlist, TargetLumi);
  
  bool _dimu =   isData ?  (k_channel.Contains("DoubleMuon")) : true ;

  //if(PassTriggerOR(triggerlist)&& _dimu)   OptimiseID(true, false, true, weight);
  
  TString analysis_trigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
  
  //if(!PassTrigger(analysis_trigger))   return;
  
  
  OptimiseID(true, true,  true, weight);
  //if(k_running_nonprompt)OptimiseID(true, true,false,  weight);
  return;

  std::vector<snu::KElectron> electronVetoColl=GetElectrons("ELECTRON_HN_VETO");
  std::vector<snu::KMuon> muonVetoColl=GetMuons("MUON_HN_VETO");

  if ((electronVetoColl.size() + muonVetoColl.size()) >2) return;
  
  CheckJetsCloseToLeptons(GetElectrons("ELECTRON_HN_VETO"), GetJets("JET_NOLEPTONVETO"));
  std::vector<snu::KFatJet> fatjetcoll = GetFatJets("FATJET_HN");
  std::vector<snu::KJet> jets = GetJets("JET_HN");
  FillHist("NFatJets", fatjetcoll.size(), weight, 0., 10., 10);
  FillHist("2DJets", fatjetcoll.size(), jets.size(),  weight, 0., 10., 10,  0., 10., 10);
  for(unsigned int ifjet=0; ifjet < fatjetcoll.size(); ifjet++){
    FillHist(("tau21"), fatjetcoll[ifjet].Tau2()/fatjetcoll[ifjet].Tau1(), weight, 0., 1., 100);
    FillHist(("PrunedMass"), fatjetcoll[ifjet].PrunedMass(),  weight, 0., 200., 100);
    FillHist(("SoftDropMass"), fatjetcoll[ifjet].SoftDropMass(),  weight, 0., 200., 100);
  }
  

}


void HNDiElectronOptimisation::CheckJetsCloseToLeptons(std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets){

  if(electrons.size() == 2) {
    if(SameCharge(electrons)) {
      snu::KParticle looseee = electrons.at(0) + electrons.at(1);
      if(fabs(looseee.M() - 90.) > 20. ){
	for(unsigned int ijet=0; ijet < jets.size(); ijet++){
	  for(unsigned int iel=0; iel < electrons.size(); iel++){
	    if( electrons[iel].DeltaR(jets.at(ijet)) < 0.4){
	      /// Jet usually vetoe
	      FillHist(("ptratio_jetel"), electrons[iel].Pt()/jets[ijet].Pt(), weight, 0., 1., 100);
	    }
	  }
	}
      }
    }
  }
  return;
}



void HNDiElectronOptimisation::OptimiseID(bool isss, bool dilep, bool removed0, float w){
  
  bool getpogid(true);
  bool getpogidcc(true);
  
  bool optmva(true);
  /// Correct MC for pileup                                                                                                                                                                                                                                                     
  if (!isData) {
    w*= eventbase->GetEvent().PileUpWeight();
  }
  

  std::vector<snu::KJet> jets = GetJets("LPJET_HN");
  /*  std::vector<snu::KJet> jets_pu = GetJets("LPJET_HN_PU");
  std::vector<snu::KJet> jets_pu1 = GetJets("LPJET_HN_PU1");
  std::vector<snu::KJet> jets_pu2 = GetJets("LPJET_HN_PU2");
  std::vector<snu::KJet> jets_pu3 = GetJets("LPJET_HN_PU3");
  */
  
  std::vector<snu::KMuon> trilepmuons;
  double this_RelIso = 0.4;

  if(k_running_nonprompt){
    trilepmuons= GetHNTriMuonsByLooseRelIso(this_RelIso, false);
  }
  else trilepmuons = GetMuons("MUON_HN_TRI_TIGHT");
  
  
  if(optmva){

    /// ELECTRON_HN_FAKELOOSEST has loosest dZ/dXY that is optimised
    //| ID                                   |  ptmin  | |etamax| |  isomax03_b   | isomax03_e | IsTight(POG)  | IsMedium(POG)  | IsLoose(POG) | IsVeto(POG)   | |dxymax_b|  | |dxymax_e| |  |dzmax_b| |  |dzmax_e| | GsfCtfScPix | convveto | |dxysigmin|  |  |dxysigmax| | UseMiniIso |  IsTight(MVA)  |  IsMedium(MVA)  | IsZZ(MVA) |

    // CATConfig/SelectionConfig/user_electrons.sel:| ELECTRON_HN_FAKELOOSEST              | 5.     | 2.5     |   0.5     | 0.50      | false         | false          |  false       | false         | 0.05    | 0.05  | 0.10   | 0.20   | true        | true     | 999.         | 4.            | false      |      false     |  false           | false |

    
    std::vector<snu::KElectron> tmploose_el = GetElectrons("ELECTRON_HN_FAKELOOSEST");
    
    /// apply trigger emulation to electrons
    /// also apply loose mva cuts
    std::vector<snu::KElectron> looseveto;

    for(unsigned int iel=0; iel<tmploose_el.size(); iel++){
      bool pass_trigger_emulation=true;
      if(tmploose_el[iel].Pt() < 15.){
        if(!tmploose_el[iel].PassHLTID()) pass_trigger_emulation=false;
      }
      else{
        if(!tmploose_el[iel].IsTrigMVAValid()) pass_trigger_emulation=false;
      }
      if(!pass_trigger_emulation) continue;
      if(fabs(tmploose_el[iel].dz()) > 0.2) continue;
      if(fabs(tmploose_el[iel].dxy()) > 0.1) continue;

      if(fabs(tmploose_el[iel].SCEta())<0.8 ){
        if(tmploose_el[iel].MVA() < -0.02) continue;
      }
      else  if(fabs(tmploose_el[iel].SCEta())<1.479 ){
        if(tmploose_el[iel].MVA() < -0.52) continue;
      }
      else {
        if(tmploose_el[iel].MVA() < -0.52) continue;
      }
      looseveto.push_back(tmploose_el[iel]);
    }


    /// require 2 leptons (veto 3 leoton events)
    if(looseveto.size() != 2) return;;

    ///// Setup cuts to optimise                                                                                                                                    
    vector<float> vcut_mva;
    vector<TString> vcut_mva_s;

    int nmva1=5;
    int nmva2=5;
    int nmva3=10;
    int nmva4=6;
    int ndxy=1;
    int ndz=1;
    int niso=5;

    Float_t  mvabins [29] = {0., 0.1,0.2,0.3,0.4, 0.5, 0.54,0.58,0.62,0.66,0.7,0.72,0.74,0.76,0.78,0.80,0.82,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98};
    
    for(unsigned int imva=0; imva < nmva1; imva++){
      float cut_dmva = float(imva)*0.10+ 0.;
      vcut_mva.push_back(cut_dmva);
      stringstream ss;
      ss << cut_dmva;
      vcut_mva_s.push_back(TString(ss.str()));
    }
    for(unsigned int imva=0; imva < nmva2; imva++){
      float cut_dmva = float(imva)*0.04+ 0.5;
      vcut_mva.push_back(cut_dmva);
      stringstream ss;
      ss << cut_dmva;
      vcut_mva_s.push_back(TString(ss.str()));
    }
    for(unsigned int imva=0; imva < nmva3; imva++){
      float cut_dmva = float(imva)*0.02 + 0.7;
      vcut_mva.push_back(cut_dmva);
      stringstream ss;
      ss << cut_dmva;
      vcut_mva_s.push_back(TString(ss.str()));
    }
    for(unsigned int imva=0; imva < nmva4; imva++){
      float cut_dmva = float(imva)*0.01 + 0.9;
      vcut_mva.push_back(cut_dmva);
      stringstream ss;
      ss << cut_dmva;
      vcut_mva_s.push_back(TString(ss.str()));
    }


    vector<float> vcut_dxy_b;
    vector<TString> vcut_dxy_b_s;
    vector<TString> vcut_dxy_b_s_fake;

    for(unsigned int dxy_b=0;dxy_b < ndxy; dxy_b++){
      float cut_dxy_b =  float(dxy_b)*0.005 + 0.01;
      float cut_dxy_b2 =  cut_dxy_b ;
      if(dxy_b % 2) cut_dxy_b2 =  cut_dxy_b + 0.005;
      vcut_dxy_b.push_back(cut_dxy_b);
      stringstream ss;
      stringstream ss_fake;
      ss <<cut_dxy_b;
      ss_fake <<cut_dxy_b2;
      vcut_dxy_b_s.push_back(TString(ss.str()));
      vcut_dxy_b_s_fake.push_back(TString(ss_fake.str()));
    }
    

    vector<float> vcut_dz_b;
    vector<TString> vcut_dz_b_s;

    for(unsigned int dz_b=0;dz_b < ndz; dz_b++){
      float cut_dz_b =  float(dz_b)*0.02 + 0.04;
      vcut_dz_b.push_back(cut_dz_b);
      stringstream ss;
      ss <<cut_dz_b;
      vcut_dz_b_s.push_back(TString(ss.str()));
    }

    vector<float> vcut_iso_b;
    vector<TString> vcut_iso_b_s;
    for(unsigned int iso_b=0;iso_b < niso; iso_b++){
      float cut_iso_b = float(iso_b)*0.01 + 0.05;
      vcut_iso_b.push_back(cut_iso_b);
      stringstream ss;
      ss <<cut_iso_b;
      vcut_iso_b_s.push_back(TString(ss.str()));
    }
    
    if(SameCharge(looseveto))FillHist("IDREF",0.  , w, 0.,2.,2);
    
    if(!SameCharge(looseveto)) return;
    
    if(looseveto[1].Pt() > 10. && looseveto[0].Pt() > 20.){
      if(PassTrigger("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v")) FillHist("LowMassTrigger",0., weight, 0., 3., 3.);
    }
    if(looseveto[1].Pt() > 10. && looseveto[0].Pt() > 18.){
      if(PassTrigger("HLT_Ele15_IsoVVVL_PFHT600_v")) FillHist("LowMassTrigger",1., weight, 0., 3., 3.);
    }
    
    if(looseveto[1].Pt() > 10. && looseveto[0].Pt() > 25.){
      if(PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) FillHist("LowMassTrigger",2., weight, 0., 3., 3.);
      snu::KParticle Z = looseveto[0] + looseveto[1];
      /*
      if(Z.M() < 80 ||  Z.M() > 100){
	if(jets.size() > 1) {
	  if(PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) FillCLHist(sighist_ee, "SIGNALOPT", eventbase->GetEvent(),  GetMuons("MUON_NOCUT"), looseveto,jets, weight);
	  else if(PassTrigger("HLT_Ele25_eta2p1_WPTight_Gsf_v"))  FillCLHist(sighist_ee, "SIGNALOPTEXTRA", eventbase->GetEvent(),  GetMuons("MUON_NOCUT"), looseveto,jets, weight);
	  else if(PassTrigger("HLT_Ele27_WPTight_Gsf_v"))  FillCLHist(sighist_ee, "SIGNALOPTEXTRA2", eventbase->GetEvent(),  GetMuons("MUON_NOCUT"), looseveto,jets, weight);
	}
									   
      */
	//if(jets.size() == 1) FillCLHist(sighist_ee, "1JETSIGNALOPT", eventbase->GetEvent(),  GetMuons("MUON_NOCUT"), looseveto,jets, weight);
      
    }
    /// only use signal like FOFO events
    bool use=false;
    if(CheckSignalRegion(isss,looseveto, jets,"Low", w)) {
      use=true;
      FillCLHist(sighist_ee, "PASSSIGNALOPTLow", eventbase->GetEvent(),  GetMuons("MUON_NOCUT"), looseveto,jets, weight);
    }
    
    if(CheckSignalRegion(isss,looseveto, jets,"",w)){
      use=true;
      FillCLHist(sighist_ee, "PASSSIGNALOPTMid", eventbase->GetEvent(),  GetMuons("MUON_NOCUT"), looseveto,jets, weight);
    }

    if(CheckSignalRegion(isss,looseveto, jets,"High", w)) {
      FillCLHist(sighist_ee, "PASSSIGNALOPTHigh", eventbase->GetEvent(),  GetMuons("MUON_NOCUT"), looseveto,jets, weight);
      use=true;
    }
    if(!use)  return;
    
    for(unsigned int dxy_b=0; dxy_b < vcut_dxy_b.size(); dxy_b++){
      for(unsigned int dxy_e=dxy_b; dxy_e < vcut_dxy_b.size(); dxy_e++){
	for(unsigned int dz_b=0; dz_b < vcut_dz_b.size(); dz_b++){

	  
	  /// loose_el here is same as loose_veto collection but IP are tightened when we change IP cuts in tight

	    std::vector<snu::KElectron> loose_el;
	    /// only neede loose_el for fake bkg
	    for(unsigned int loose_iel=0; loose_iel<looseveto.size(); loose_iel++){
	      float reliso = looseveto[loose_iel].PFRelIso(0.3);
	      if(fabs(looseveto[loose_iel].dz()) > vcut_dz_b[dz_b]) continue;
	      
	      if(fabs(looseveto[loose_iel].SCEta())<0.8 ){
		if(fabs(looseveto[loose_iel].dxy()) > vcut_dxy_b[dxy_b]) continue;
	      }
	      else  if(fabs(looseveto[loose_iel].SCEta())<1.479 ){
		if(fabs(looseveto[loose_iel].dxy()) > vcut_dxy_b[dxy_b]) continue;
	      }
	      else {
		if(fabs(looseveto[loose_iel].dxy()) > vcut_dxy_b[dxy_e]) continue;
	      }
	      loose_el.push_back(looseveto[loose_iel]);
	    }
	    if(loose_el.size() < 2) continue;
	    for(unsigned int ipt=0; ipt < 4; ipt++){
	      float ptmin(0.);
	      float ptmax(0.);
	      TString ptlabel="";
	      
	      if(ipt==0){
		ptmin=10.;
		ptmax=15.;
		ptlabel="pt1";
	      }
	      if(ipt==1){
		ptmin=15.;
		ptmax=25.;
		ptlabel="pt2";
	      }
	      if(ipt==2){
		ptmin=25.;
		ptmax=50.;
		ptlabel="pt3";
	      }
	      if(ipt==3){
		ptmin=50.;
		ptmax=1000.;
		ptlabel="pt4";
	      }
	    
	    
	      
	      for(unsigned int iso_b=0; iso_b < vcut_iso_b.size(); iso_b++){
		for(unsigned int iso_e=0; iso_e < vcut_iso_b.size(); iso_e++){
		  for(unsigned int imva3=0; imva3 < vcut_mva.size(); imva3++){
		    for(unsigned int imva2=imva3; imva2 < vcut_mva.size(); imva2++){
		      for(unsigned int imva1=imva2; imva1 < vcut_mva.size(); imva1++){
			
			int etareg_1(0);
			int etareg_2(0);
			float mva_1(0);
			float mva_2(0);
			float mva_3(0);
			float el_mva1(0);
			float el_mva2(0);
			std::vector<snu::KElectron> tight_el;
			bool tight1=false;
			bool tight2=false; 

			
			
			for(unsigned int Liel=0; Liel<loose_el.size(); Liel++){
			  if(loose_el[Liel].Pt() < ptmax && loose_el[Liel].Pt()> ptmin){
			    mva_1=vcut_mva[imva1];
			    mva_2=vcut_mva[imva2];
			    mva_3=vcut_mva[imva3];
			  }
			  else{
			    mva_1=0.941;
			    mva_2=0.899;
			    mva_3=0.758;
			  }

			  float reliso = loose_el[Liel].PFRelIso(0.3);

			  if(fabs(loose_el[Liel].SCEta())<0.8 ){
                            if(Liel==0)el_mva1=mva_1;
                            if(Liel==1)el_mva2=mva_1;
                            if(Liel==0)etareg_1=1;
                            if(Liel==1)etareg_2=1;
			  }
			  else  if(fabs(loose_el[Liel].SCEta())<1.479 ){
			    if(Liel==0)el_mva1=mva_2;
                            if(Liel==1)el_mva2=mva_2;
                            if(Liel==0)etareg_1=2;
                            if(Liel==1)etareg_2=2;
			  }
			  else {
                            if(Liel==0)el_mva1=mva_3;
                            if(Liel==1)el_mva2=mva_3;
                            if(Liel==0)etareg_1=3;
                            if(Liel==1)etareg_2=3;

			  }

			  
			  bool mvatight=false;
			  if(fabs(loose_el[Liel].SCEta())<0.8 ){
			    if(loose_el[Liel].MVA() > mva_1) mvatight=true;
			  }
			  else  if(fabs(loose_el[Liel].SCEta())<1.479 ){
			    if(loose_el[Liel].MVA() > mva_2) mvatight=true;
			  }
			  else {
			    if(loose_el[Liel].MVA() > mva_3)mvatight=true;
			  }
			  if(!mvatight) continue;
			  if(fabs(loose_el[Liel].SCEta())<1.479 ){
			    if(reliso > vcut_iso_b[iso_b]) continue;
			  }
			  else{
			    if(reliso > vcut_iso_b[iso_e]) continue;
			  }

			  if(Liel==0)tight1=true;
			  if(Liel==1)tight2=true;

			  tight_el.push_back(loose_el[Liel]);
			  
			}
			if(!k_running_nonprompt){
			  if(CheckSignalRegion(isss,tight_el, jets,"Low", w))FillHist(("LowIDREF_"+ptlabel+"_"+vcut_mva_s[imva1]+"iso"+vcut_iso_b_s[iso_b] +"_"+vcut_iso_b_s[iso_e]+"_dxy"+vcut_dxy_b_s[dxy_b]+ vcut_dxy_b_s[dxy_e]+"_dz"+vcut_dz_b_s[dz_b]),vcut_mva[imva2],vcut_mva[imva3] ,w, mvabins,28,mvabins,28);
			  if(CheckSignalRegion(isss,tight_el, jets,"", w))FillHist(("MidIDREF_"+ptlabel+"_"+vcut_mva_s[imva1]+"iso"+vcut_iso_b_s[iso_b] +"_"+vcut_iso_b_s[iso_e]+"_dxy"+vcut_dxy_b_s[dxy_b]+ vcut_dxy_b_s[dxy_e]+"_dz"+vcut_dz_b_s[dz_b]),vcut_mva[imva2],vcut_mva[imva3] ,w, mvabins,28,mvabins,28);
			  if(CheckSignalRegion(isss,tight_el, jets,"High", w))FillHist(("HighIDREF_"+ptlabel+"_"+vcut_mva_s[imva1]+"iso"+vcut_iso_b_s[iso_b] +"_"+vcut_iso_b_s[iso_e]+"_dxy"+vcut_dxy_b_s[dxy_b]+ vcut_dxy_b_s[dxy_e]+"_dz"+vcut_dz_b_s[dz_b]),vcut_mva[imva2],vcut_mva[imva3] ,w, mvabins,28,mvabins,28);

			}
			else{
				      
			  float mva_el(0.);
			  TString mvakey(""), isokey(""), dxykey("");
			  TString mvakey2(""), isokey2(""), dxykey2("");
			  if(etareg_1==1){
			    mvakey=vcut_mva_s[imva1];
			    isokey=vcut_iso_b_s[iso_b];
			    dxykey=vcut_dxy_b_s_fake[dxy_b];
			  }
			  if(etareg_2==1){
			    mvakey2=vcut_mva_s[imva1];
			    isokey2=vcut_iso_b_s[iso_b];
			    dxykey2=vcut_dxy_b_s_fake[dxy_b];
			  }
			  if(etareg_1==2){
			    mvakey=vcut_mva_s[imva2];
			    isokey=vcut_iso_b_s[iso_b];
			    dxykey=vcut_dxy_b_s_fake[dxy_b];
			  }
			  if(etareg_2==2){
			    mvakey2=vcut_mva_s[imva2];
			    isokey2=vcut_iso_b_s[iso_b];
			    dxykey2=vcut_dxy_b_s_fake[dxy_b];
			  }
			  if(etareg_1==3){
			    mvakey=vcut_mva_s[imva3];
			    isokey=vcut_iso_b_s[iso_e];
			    dxykey=vcut_dxy_b_s_fake[dxy_e];
			  }
			  if(etareg_2==3){
			    mvakey2=vcut_mva_s[imva3];
			    isokey2=vcut_iso_b_s[iso_e];
			    dxykey2=vcut_dxy_b_s_fake[dxy_e];
			  }

			  if(CheckSignalRegion(isss,loose_el, jets,"Low", w)) {

			    float ee_weight_mva =  1.;//m_datadriven_bkg->Get_DataDrivenWeight_EEmva(false, loose_el, tight1,tight2,mvakey+"_iso" +isokey+"_dxy"+dxykey+"_dz"+vcut_dz_b_s[dz_b], mvakey2+"_iso" +isokey2+"_dxy"+dxykey2+"_dz"+vcut_dz_b_s[dz_b]); 															
			    FillHist(("LowIDREF_"+ptlabel+"_"+vcut_mva_s[imva1]+"iso"+vcut_iso_b_s[iso_b]+"_"+vcut_iso_b_s[iso_e]+"_dxy"+vcut_dxy_b_s[dxy_b]+ vcut_dxy_b_s[dxy_e]+"_dz"+vcut_dz_b_s[dz_b]),vcut_mva[imva2],vcut_mva[imva3] , ee_weight_mva, mvabins,28,mvabins,28);
			  }
			  
			  if(CheckSignalRegion(isss,loose_el, jets,"", w)) {
			    
                            float ee_weight_mva = 1.;// m_datadriven_bkg->Get_DataDrivenWeight_EEmva(false, loose_el, tight1,tight2,mvakey+"_iso" +isokey+"_dxy"+dxykey+"_dz"+vcut_dz_b_s[dz_b], mvakey2+"_iso" +isokey2+"_dxy"+dxykey2+"_dz"+vcut_dz_b_s[dz_b]);
												
			    
                            FillHist(("MidIDREF_"+ptlabel+"_"+vcut_mva_s[imva1]+"iso"+vcut_iso_b_s[iso_b]+"_"+vcut_iso_b_s[iso_e]+"_dxy"+vcut_dxy_b_s[dxy_b]+ vcut_dxy_b_s[dxy_e]+"_dz"+vcut_dz_b_s[dz_b]),vcut_mva[imva2],vcut_mva[imva3] , ee_weight_mva, mvabins,28,mvabins,28);
				      
                          }
			  if(CheckSignalRegion(isss,loose_el, jets,"High", w)) {

                            float ee_weight_mva = 1.;// m_datadriven_bkg->Get_DataDrivenWeight_EEmva(false, loose_el, tight1,tight2,mvakey+"_iso" +isokey+"_dxy"+dxykey+"_dz"+vcut_dz_b_s[dz_b], mvakey2+"_iso" +isokey2+"_dxy"+dxykey2+"_dz"+vcut_dz_b_s[dz_b]); 
												
			    
                            FillHist(("HighIDREF_"+ptlabel+"_"+vcut_mva_s[imva1]+"iso"+vcut_iso_b_s[iso_b]+"_"+vcut_iso_b_s[iso_e]+"_dxy"+vcut_dxy_b_s[dxy_b]+ vcut_dxy_b_s[dxy_e]+"_dz"+vcut_dz_b_s[dz_b]),vcut_mva[imva2],vcut_mva[imva3] , ee_weight_mva, mvabins,28,mvabins,28);
				      
			    }
			  
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  return;;

  if(getpogidcc){
    std::vector<snu::KElectron> electronLooseColl_tight = GetElectrons("ELECTRON16_POG_FAKELOOSE_CC");
    std::vector<snu::KElectron> electronLooseColl_medium = GetElectrons("ELECTRON16_POG_MEDIUM_FAKELOOSE_CC");
    std::vector<snu::KElectron> electronLooseColl_mva = GetElectrons("ELECTRON16_MVA_FAKELOOSE_CC");
    
    std::vector<snu::KElectron> electronMedium_chargeconst  = GetElectrons("ELECTRON16_FR_POG_MEDIUM_CC");
    std::vector<snu::KElectron> electronTight_chargeconst  = GetElectrons("ELECTRON16_FR_POG_TIGHT_CC");
    std::vector<snu::KElectron>  electronMVAColl_mva = GetElectrons("ELECTRON16_FR_MVA_TIGHT_CC");

    std::vector<snu::KElectron> electronMedium_chargeconst_dxy  = GetElectrons("ELECTRON16_FR_POG_MEDIUM_DXYCC");
    std::vector<snu::KElectron> electronTight_chargeconst_dxy  = GetElectrons("ELECTRON16_FR_POG_TIGHT_DXYCC");
    std::vector<snu::KElectron>  electronMVAColl_mva_dxy = GetElectrons("ELECTRON16_FR_MVA_TIGHT_DXYCC");

    if(!k_running_nonprompt){
      
      if(dilep){
	FillHist("IDREF",0.  , w, 0.,2.,2);
	if(CheckSignalRegion(isss,electronLooseColl_tight, jets,"Signal_Mediumlooseiso_nod0", w)) FillHist("IDREF",1.  , w, 0.,2.,2);
	if(CheckSignalRegion(isss,electronMedium_chargeconst, jets,"Signal_Mediumlooseiso_nod0", w)) FillHist("IDcutflow",0.  , w, 0.,7.,7);
	if(CheckSignalRegion(isss,electronTight_chargeconst, jets,"Signal_Tightlooseiso_nod0", w))   FillHist("IDcutflow",1.  , w, 0.,7.,7);
	if(CheckSignalRegion(isss,electronMVAColl_mva, jets,"Signal_MVAlooseiso_nod0", w))   FillHist("IDcutflow",2.  , w, 0.,7.,7);
	
	if(CheckSignalRegion(isss,electronMedium_chargeconst_dxy, jets,"Signal_Mediumlooseiso_nod0", w)) FillHist("IDcutflow",3.  , w, 0.,7.,7);
	if(CheckSignalRegion(isss,electronTight_chargeconst_dxy, jets,"Signal_Tightlooseiso_nod0", w))   FillHist("IDcutflow",4.  , w, 0.,7.,7);
	if(CheckSignalRegion(isss,electronMVAColl_mva_dxy, jets,"Signal_MVAlooseiso_nod0", w))   FillHist("IDcutflow",5.  , w, 0.,7.,7);

	FillHist("IDREF_NN",0.  , w, 0.,2.,2);

        if(CheckSignalRegionNN(isss,electronLooseColl_tight, jets,"Signal_Mediumlooseiso_nod0", w)) FillHist("NNIDREF",0.  , w, 0.,2.,2);
        if(CheckSignalRegionNN(isss,electronMedium_chargeconst, jets,"Signal_Mediumlooseiso_nod0", w)) FillHist("NNIDcutflow",0.  , w, 0.,7.,7);
        if(CheckSignalRegionNN(isss,electronTight_chargeconst, jets,"Signal_Tightlooseiso_nod0", w))   FillHist("NNIDcutflow",1.  , w, 0.,7.,7);
        if(CheckSignalRegionNN(isss,electronMVAColl_mva, jets,"Signal_MVAlooseiso_nod0", w))   FillHist("NNIDcutflow",2.  , w, 0.,7.,7);

        if(CheckSignalRegionNN(isss,electronMedium_chargeconst_dxy, jets,"Signal_Mediumlooseiso_nod0", w)) FillHist("NNIDcutflow",3.  , w, 0.,7.,7);
        if(CheckSignalRegionNN(isss,electronTight_chargeconst_dxy, jets,"Signal_Tightlooseiso_nod0", w))   FillHist("NNIDcutflow",4.  , w, 0.,7.,7);
        if(CheckSignalRegionNN(isss,electronMVAColl_mva_dxy, jets,"Signal_MVAlooseiso_nod0", w))   FillHist("NNIDcutflow",5.  , w, 0.,7.,7);


      }
      else {
	
	
	/// mm+ew
	FillHist("IDREF_trilep",0.  , w, 0.,2.,2);
	
	if(CheckSignalRegionTriLep(trilepmuons, electronLooseColl_tight, jets,"Signal_Mediumlooseiso_nod0", w)) FillHist("IDREF_trilep",0.  , w, 0.,2.,2);
	if(CheckSignalRegionTriLep(trilepmuons, electronMedium_chargeconst, jets,"Signal_Mediumlooseiso_nod0", w)) FillHist("IDcutflow_trilep",0.  , w, 0.,7.,7);
	if(CheckSignalRegionTriLep(trilepmuons, electronTight_chargeconst, jets,"Signal_Tightlooseiso_nod0", w))   FillHist("IDcutflow_trilep",1.  , w, 0.,7.,7);
	if(CheckSignalRegionTriLep(trilepmuons, electronMVAColl_mva, jets,"Signal_MVAlooseiso_nod0", w))   FillHist("IDcutflow_trilep",2.  , w, 0.,7.,7);
	
	if(CheckSignalRegionTriLep(trilepmuons, electronMedium_chargeconst_dxy, jets,"Signal_Mediumlooseiso_nod0", w)) FillHist("IDcutflow_trilep",3.  , w, 0.,7.,7);
	if(CheckSignalRegionTriLep(trilepmuons, electronTight_chargeconst_dxy, jets,"Signal_Tightlooseiso_nod0", w))   FillHist("IDcutflow_trilep",4.  , w, 0.,7.,7);
	if(CheckSignalRegionTriLep(trilepmuons, electronMVAColl_mva_dxy, jets,"Signal_MVAlooseiso_nod0", w))   FillHist("IDcutflow_trilep",5.  , w, 0.,7.,7);

      }
      
    }
    else{
      
      TString idtag="";
      

      if(dilep){
	TString medium_looselabel="ELECTRON16_POG_MEDIUM_FAKELOOSE_CC";
	TString tight_looselabel="ELECTRON16_POG_FAKELOOSE_CC";
	TString mva_looselabel="ELECTRON16_MVA_FAKELOOSE_CC";
	if(!removed0){
	  idtag="_d0";
	  medium_looselabel="ELECTRON16_POG_MEDIUM_FAKELOOSE_CC_d0";
	  tight_looselabel="ELECTRON16_POG_FAKELOOSE_CC_d0";
	  mva_looselabel="ELECTRON16_MVA_FAKELOOSE_CC_d0";
	}


	float ee_weight_medium = m_datadriven_bkg->Get_DataDrivenWeight_EE(false, electronLooseColl_medium,medium_looselabel,"ELECTRON16_FR_POG_MEDIUM_CC","dijet_ajet40"+idtag);
	float ee_weight_tight =   m_datadriven_bkg->Get_DataDrivenWeight_EE(false, electronLooseColl_tight,tight_looselabel,"ELECTRON16_FR_POG_TIGHT_CC","dijet_ajet40"+idtag);
	float ee_weight_mva =  m_datadriven_bkg->Get_DataDrivenWeight_EE(false, electronLooseColl_mva, mva_looselabel, "ELECTRON16_FR_MVA_TIGHT_CC","dijet_ajet40"+idtag);
	
	float ee_weight_medium_dxy = m_datadriven_bkg->Get_DataDrivenWeight_EE(false, electronLooseColl_medium,medium_looselabel,"ELECTRON16_FR_POG_MEDIUM_DXYCC","dijet_ajet40"+idtag);
	float ee_weight_tight_dxy =   m_datadriven_bkg->Get_DataDrivenWeight_EE(false, electronLooseColl_tight,tight_looselabel,"ELECTRON16_FR_POG_TIGHT_DXYCC","dijet_ajet40"+idtag);
	float ee_weight_mva_dxy =  m_datadriven_bkg->Get_DataDrivenWeight_EE(false, electronLooseColl_mva, mva_looselabel, "ELECTRON16_FR_MVA_TIGHT_DXYCC","dijet_ajet40"+idtag);
	
	if(CheckSignalRegion(isss, electronLooseColl_medium, jets,"Signal_Mediumlooseiso_nod0", w)) FillHist("IDcutflow"+idtag,0.  , ee_weight_medium, 0.,7.,7);
	if(CheckSignalRegion(isss,electronLooseColl_tight, jets,"Signal_Tightlooseiso_nod0", w))   FillHist("IDcutflow"+idtag,1.  , ee_weight_tight, 0.,7.,7);
	if(CheckSignalRegion(isss,electronLooseColl_mva, jets,"Signal_MVAlooseiso_nod0", w))   FillHist("IDcutflow"+idtag,2.  , ee_weight_mva, 0.,7.,7);
	
	if(CheckSignalRegion(isss, electronLooseColl_medium, jets,"Signal_Mediumlooseiso_nod0", w)) FillHist("IDcutflow"+idtag,3.  , ee_weight_medium_dxy, 0.,7.,7);
	if(CheckSignalRegion(isss,electronLooseColl_tight, jets,"Signal_Tightlooseiso_nod0", w))   FillHist("IDcutflow"+idtag,4.  , ee_weight_tight_dxy, 0.,7.,7);
	if(CheckSignalRegion(isss,electronLooseColl_mva, jets,"Signal_MVAlooseiso_nod0", w))   FillHist("IDcutflow"+idtag,5.  , ee_weight_mva_dxy, 0.,7.,7);

	
	if(CheckSignalRegionNN(isss, electronLooseColl_medium, jets,"Signal_Mediumlooseiso_nod0", w)) FillHist("NNIDcutflow"+idtag,0.  , ee_weight_medium, 0.,7.,7);
        if(CheckSignalRegionNN(isss,electronLooseColl_tight, jets,"Signal_Tightlooseiso_nod0", w))   FillHist("NNIDcutflow"+idtag,1.  , ee_weight_tight, 0.,7.,7);
        if(CheckSignalRegionNN(isss,electronLooseColl_mva, jets,"Signal_MVAlooseiso_nod0", w))   FillHist("NNIDcutflow"+idtag,2.  , ee_weight_mva, 0.,7.,7);

        if(CheckSignalRegionNN(isss, electronLooseColl_medium, jets,"Signal_Mediumlooseiso_nod0", w)) FillHist("NNIDcutflow"+idtag,3.  , ee_weight_medium_dxy, 0.,7.,7);
        if(CheckSignalRegionNN(isss,electronLooseColl_tight, jets,"Signal_Tightlooseiso_nod0", w))   FillHist("NNIDcutflow"+idtag,4.  , ee_weight_tight_dxy, 0.,7.,7);
        if(CheckSignalRegionNN(isss,electronLooseColl_mva, jets,"Signal_MVAlooseiso_nod0", w))   FillHist("NNIDcutflow"+idtag,5.  , ee_weight_mva_dxy, 0.,7.,7);


	
      }
      else{
	
	TString medium_looselabel="ELECTRON16_POG_MEDIUM_FAKELOOSE_CC";
	TString tight_looselabel="ELECTRON16_POG_FAKELOOSE_CC";
	TString mva_looselabel="ELECTRON16_MVA_FAKELOOSE_CC";
	if(!removed0){
	  idtag="_d0";
	  medium_looselabel="ELECTRON16_POG_MEDIUM_FAKELOOSE_CC_d0";
	  tight_looselabel="ELECTRON16_POG_FAKELOOSE_CC_d0";
	  mva_looselabel="ELECTRON16_MVA_FAKELOOSE_CC_d0";
	}
	
	float mme_weight_medium = m_datadriven_bkg->Get_DataDrivenWeight(false, trilepmuons, "MUON_HN_TRI_TIGHT", 2,electronLooseColl_medium, "ELECTRON16_FR_POG_MEDIUM_CC",1,medium_looselabel, "dijet_ajet40"+idtag);
	float mme_weight_tight =   m_datadriven_bkg->Get_DataDrivenWeight(false,trilepmuons, "MUON_HN_TRI_TIGHT", 2,electronLooseColl_tight, "ELECTRON16_FR_POG_TIGHT_CC",1,tight_looselabel,"dijet_ajet40"+idtag);
	float mme_weight_mva =  m_datadriven_bkg->Get_DataDrivenWeight(false, trilepmuons, "MUON_HN_TRI_TIGHT", 2,electronLooseColl_mva, "ELECTRON16_FR_MVA_TIGHT_CC", 1,mva_looselabel,"dijet_ajet40"+idtag);
	
	float mme_weight_medium_dxy = m_datadriven_bkg->Get_DataDrivenWeight(false, trilepmuons, "MUON_HN_TRI_TIGHT", 2,electronLooseColl_medium, "ELECTRON16_FR_POG_MEDIUM_DXYCC",1,medium_looselabel,"dijet_ajet40"+idtag);
	float mme_weight_tight_dxy =   m_datadriven_bkg->Get_DataDrivenWeight(false, trilepmuons, "MUON_HN_TRI_TIGHT", 2,electronLooseColl_tight, "ELECTRON16_FR_POG_TIGHT_DXYCC",1,tight_looselabel,"dijet_ajet40"+idtag);
	float mme_weight_mva_dxy =  m_datadriven_bkg->Get_DataDrivenWeight(false, trilepmuons, "MUON_HN_TRI_TIGHT", 2,electronLooseColl_mva,"ELECTRON16_FR_MVA_TIGHT_DXYCC", 1,mva_looselabel,"dijet_ajet40"+idtag);
	
	if(CheckSignalRegionTriLep( trilepmuons, electronLooseColl_medium, jets,"Signal_Mediumlooseiso_nod0", w)) FillHist("IDcutflow_trilep"+idtag,0.  , mme_weight_medium, 0.,7.,7);
	if(CheckSignalRegionTriLep( trilepmuons, electronLooseColl_tight, jets,"Signal_Tightlooseiso_nod0", w))   FillHist("IDcutflow_trilep"+idtag,1.  , mme_weight_tight, 0.,7.,7);
	if(CheckSignalRegionTriLep( trilepmuons, electronLooseColl_mva, jets,"Signal_MVAlooseiso_nod0", w))   FillHist("IDcutflow_trilep"+idtag,2.  , mme_weight_mva, 0.,7.,7);
	
	if(CheckSignalRegionTriLep( trilepmuons, electronLooseColl_medium, jets,"Signal_Mediumlooseiso_nod0", w)) FillHist("IDcutflow_trilep"+idtag,3.  , mme_weight_medium_dxy, 0.,7.,7);
	if(CheckSignalRegionTriLep( trilepmuons, electronLooseColl_tight, jets,"Signal_Tightlooseiso_nod0", w))   FillHist("IDcutflow_trilep"+idtag,4.  , mme_weight_tight_dxy, 0.,7.,7);
	if(CheckSignalRegionTriLep( trilepmuons, electronLooseColl_mva, jets,"Signal_MVAlooseiso_nod0", w))   FillHist("IDcutflow_trilep"+idtag,5.  , mme_weight_mva_dxy, 0.,7.,7);
	
      }
    }
  }
  
  if(getpogid){
    std::vector<snu::KElectron> electronLooseColl_tight = GetElectrons("ELECTRON16_POG_FAKELOOSE");
    std::vector<snu::KElectron> electronLooseColl_medium = GetElectrons("ELECTRON16_POG_MEDIUM_FAKELOOSE");
    
    std::vector<snu::KElectron> electronMedium_chargeconst  = GetElectrons("ELECTRON16_FR_POG_MEDIUM");
    std::vector<snu::KElectron> electronTight_chargeconst  = GetElectrons("ELECTRON16_FR_POG_TIGHT");
    
    std::vector<snu::KJet> jets = GetJets("JET_HN");
    
    if(!k_running_nonprompt){
      
      if(dilep){
	if(CheckSignalRegion(isss,electronLooseColl_tight, jets,"Signal_ooseiso_d0", w)) FillHist("IDREF_nocc",0.  , w, 0.,1.,1);
	if(CheckSignalRegion(isss,electronMedium_chargeconst, jets,"Signal_Mediumlooseiso_d0", w)) FillHist("IDcutflow_nocc",0.  , w, 0.,5.,5);
	if(CheckSignalRegion(isss,electronTight_chargeconst, jets,"Signal_Tightlooseiso_d0", w))   FillHist("IDcutflow_nocc",1.  , w, 0.,5.,5);

	if(CheckSignalRegionNN(isss,electronLooseColl_tight, jets,"Signal_ooseiso_d0", w)) FillHist("NNIDREF_nocc",0.  , w, 0.,1.,1);
        if(CheckSignalRegionNN(isss,electronMedium_chargeconst, jets,"Signal_Mediumlooseiso_d0", w)) FillHist("NNIDcutflow_nocc",0.  , w, 0.,5.,5);
        if(CheckSignalRegionNN(isss,electronTight_chargeconst, jets,"Signal_Tightlooseiso_d0", w))   FillHist("NNIDcutflow_nocc",1.  , w, 0.,5.,5);

      }
      else{
	if(CheckSignalRegionTriLep(trilepmuons,electronLooseColl_tight, jets,"Signal_Mediumlooseiso_d0", w)) FillHist("IDREF_nocc_trilep",0.  , w, 0.,1.,1);
        if(CheckSignalRegionTriLep(trilepmuons,electronMedium_chargeconst, jets,"Signal_Mediumlooseiso_d0", w)) FillHist("IDcutflow_nocc_trilep",0.  , w, 0.,5.,5);
        if(CheckSignalRegionTriLep(trilepmuons,electronTight_chargeconst, jets,"Signal_Tightlooseiso_d0", w))   FillHist("IDcutflow_nocc_trilep",1.  , w, 0.,5.,5);

      }
      

    }
    else{

      if(dilep){

	float ee_weight_medium = m_datadriven_bkg->Get_DataDrivenWeight_EE(false, electronLooseColl_medium,"ELECTRON16_POG_MEDIUM_FAKELOOSE","ELECTRON16_FR_POG_MEDIUM","dijet_ajet40");
	float ee_weight_tight =   m_datadriven_bkg->Get_DataDrivenWeight_EE(false, electronLooseColl_tight,"ELECTRON16_POG_FAKELOOSE","ELECTRON16_FR_POG_TIGHT","dijet_ajet40");

	if(CheckSignalRegion(isss, electronLooseColl_medium, jets,"Signal_Mediumlooseiso_d0", w)) FillHist("IDcutflow_nocc",0.  , ee_weight_medium, 0.,2.,2);
	if(CheckSignalRegion(isss,electronLooseColl_tight, jets,"Signal_Tightlooseiso_d0", w))   FillHist("IDcutflow_nocc",1.  , ee_weight_tight, 0.,2.,2);
	
	if(CheckSignalRegionNN(isss, electronLooseColl_medium, jets,"Signal_Mediumlooseiso_d0", w)) FillHist("NNIDcutflow_nocc",0.  , ee_weight_medium, 0.,2.,2);
        if(CheckSignalRegionNN(isss,electronLooseColl_tight, jets,"Signal_Tightlooseiso_d0", w))   FillHist("NNIDcutflow_nocc",1.  , ee_weight_tight, 0.,2.,2);

      }
      else{

	float mme_weight_medium = m_datadriven_bkg->Get_DataDrivenWeight(false, trilepmuons, "MUON_HN_TRI_TIGHT", 2, electronLooseColl_medium,"ELECTRON16_FR_POG_MEDIUM", 1,"ELECTRON16_POG_MEDIUM_FAKELOOSE","dijet_ajet40");
	float mme_weight_tight =   m_datadriven_bkg->Get_DataDrivenWeight(false, trilepmuons, "MUON_HN_TRI_TIGHT", 2,electronLooseColl_tight,"ELECTRON16_FR_POG_TIGHT", 1, "ELECTRON16_POG_FAKELOOSE","dijet_ajet40");

	if(CheckSignalRegionTriLep( trilepmuons, electronLooseColl_medium, jets,"Signal_Mediumlooseiso_d0", w)) FillHist("IDcutflow_nocc_trilep",0.  , mme_weight_medium, 0.,2.,2);
        if(CheckSignalRegionTriLep( trilepmuons, electronLooseColl_tight, jets,"Signal_Tightlooseiso_d0", w))   FillHist("IDcutflow_nocc_trilep",1.  , mme_weight_tight, 0.,2.,2);

      }
    }
  }
  
  std::vector<TString> isotag;
  isotag.push_back("");
  isotag.push_back("dxysig_");
  isotag.push_back("miniiso_dxysig_");
  
  if(dilep){
    if(CheckSignalRegion(isss, GetElectrons("HNTight_isoref"),jets,"", w)) FillHist("ISOREF", 0., w, 0.,1.,1);
    if(CheckSignalRegion(isss, GetElectrons("HNTight_isoref_dxysig"),jets,"", w)) FillHist("ISOREFdxysig", 0., w, 0.,1.,1);
    if(CheckSignalRegionNN(isss, GetElectrons("HNTight_isoref"),jets,"", w)) FillHist("NNISOREF", 0., w, 0.,1.,1);
  }
  else {
    if(CheckSignalRegionTriLep( trilepmuons, GetElectrons("HNTight_isoref"),jets,"", w)) FillHist("ISOREF_trilep", 0., w, 0.,1.,1);
  }
  for(unsigned int s = 0; s < isotag.size() ; s++){
    TString itag=isotag[s];
    
    
    if(!k_running_nonprompt){
      std::vector <TString> isocuts;
      isocuts.push_back("050");
      isocuts.push_back("0525");
      isocuts.push_back("055");
      isocuts.push_back("060");
      isocuts.push_back("065");
      isocuts.push_back("075");
      isocuts.push_back("100");
      isocuts.push_back("125");
      int iisocut(0);
      for(unsigned int iiso = 0 ; iiso < isocuts.size(); iiso++){
	for(unsigned int iiso2 = 0 ; iiso2 < isocuts.size(); iiso2++, iisocut++){
	  if(dilep)	{
	    if(CheckSignalRegion(isss, GetElectrons(("HNTight_"+itag+"b"+isocuts[iiso]+"_e"+isocuts[iiso2]).Data()),jets,"", w)) FillHist("ISOcutflow"+itag, iisocut , w,0., 64., 64);
	    if(CheckSignalRegionNN(isss, GetElectrons(("HNTight_"+itag+"b"+isocuts[iiso]+"_e"+isocuts[iiso2]).Data()),jets,"", w)) FillHist("NNISOcutflow"+itag, iisocut , w,0., 64., 64);
	  }
	  else if(CheckSignalRegionTriLep(trilepmuons, GetElectrons(("HNTight_"+itag+"b"+isocuts[iiso]+"_e"+isocuts[iiso2]).Data()),jets,"", w)) FillHist("ISOcutflow"+itag+"_tchan", iisocut , w,0., 64., 64);
	  
	}
      }
    }
    else{
      
      std::vector <TString> isocuts;
      isocuts.push_back("050");
      isocuts.push_back("0525");
      isocuts.push_back("055");
      isocuts.push_back("060");
      isocuts.push_back("065");
      isocuts.push_back("075");
      isocuts.push_back("100");
      isocuts.push_back("125");
      
      int iisocut(0);
      for(unsigned int iiso = 0 ; iiso < isocuts.size(); iiso++){
	for(unsigned int iiso2 = 0 ; iiso2 < isocuts.size(); iiso2++, iisocut++){
	  std::vector<snu::KElectron>  eltight = GetElectrons(("HNTight_"+itag+"b"+isocuts[iiso]+"_e"+isocuts[iiso2]).Data());

	  TString idtag="";
	  TString looseid="ELECTRON16_HN_FAKELOOSE_NOD0";
	  if(!removed0){
	    looseid="ELECTRON16_HN_FAKELOOSE_D0";
	    idtag="_d0";
	  }

	  


	  std::vector<snu::KElectron>  elloosetmp=  GetElectrons(looseid);
	  std::vector<snu::KElectron>  elloose;
	  if(!removed0){
	    for(unsigned int iel=0; iel < elloosetmp.size(); iel++){
	      if(fabs(elloosetmp[iel].dxy()) < GetDXYCut(elloosetmp[iel],"b020_e020")){
		elloose.push_back(elloosetmp[iel]);
	      }
	    }
	  }
	  else elloose = elloosetmp;

	  if(dilep)  {

	    float evw = m_datadriven_bkg->Get_DataDrivenWeight_EE(false,  elloose, looseid,("HNTight_"+itag+"b"+isocuts[iiso]+"_e"+isocuts[iiso2]).Data(), "opt_dijet_ajet40"+idtag);
	    
	    if(CheckSignalRegion(isss, elloose,jets,"", w)) FillHist("ISOcutflow"+itag+idtag, iisocut , evw,0., 64., 64);
	    if(CheckSignalRegionNN(isss, elloose,jets,"", w)) FillHist("NNISOcutflow"+itag+idtag, iisocut , evw,0., 64., 64);
	  }
	  else {

	    float mme_evw =   m_datadriven_bkg->Get_DataDrivenWeight(false, trilepmuons, "MUON_HN_TRI_TIGHT", 2, elloose,("HNTight_"+itag+"b"+isocuts[iiso]+"_e"+isocuts[iiso2]).Data(),1, looseid, "opt_dijet_ajet40"+idtag);
	    if(CheckSignalRegionTriLep(trilepmuons, elloose,jets,"", w)) FillHist("ISOcutflow"+itag+"_trilep"+idtag, iisocut , mme_evw,0., 64., 64);
	  }
	}
      }
    }
  } 
  
  if(dilep){
    if(CheckSignalRegion(isss,GetElectrons("HNTight_dxy_ref"),jets,"", w)) FillHist("d0cutflow", 0., w, 0.,30.,30);
    if(CheckSignalRegion(isss,GetElectrons("HNTight_dxy_miniiso_ref"),jets,"", w)) FillHist("d0cutflow_miniiso", 0., w, 0.,30.,30);
    if(CheckSignalRegionNN(isss,GetElectrons("HNTight_dxy_ref"),jets,"", w)) FillHist("NNd0cutflow", 0., w, 0.,30.,30);
  }
  else     if(CheckSignalRegionTriLep( trilepmuons, GetElectrons("HNTight_dxy_ref"),jets,"", w)) FillHist("d0cutflow_trilep", 0., w, 0.,30.,30);

  if(!k_running_nonprompt){

    std::vector<TString> dxyb;
    dxyb.push_back("050");
    dxyb.push_back("025");
    dxyb.push_back("015");
    dxyb.push_back("010");
    std::vector<TString> dxye;
    dxye.push_back("100");
    dxye.push_back("050");
    dxye.push_back("040");
    dxye.push_back("025");
    dxye.push_back("020");
    dxye.push_back("010");

    std::vector<TString> dxyb2;
    dxyb2.push_back("017");
    dxyb2.push_back("015");
    dxyb2.push_back("014");
    dxyb2.push_back("013");
    dxyb2.push_back("012");
    dxyb2.push_back("011");
    dxyb2.push_back("010");
    std::vector<TString> dxye2;
    dxye2.push_back("017");
    dxye2.push_back("015");
    dxye2.push_back("014");
    dxye2.push_back("013");
    dxye2.push_back("012");
    dxye2.push_back("011");
    dxye2.push_back("010");

    if(1){
      int idxy2(1);
      for(unsigned int ib=0; ib <dxyb2.size(); ib++){
	for(unsigned int ie=0; ie <dxye2.size(); ie++, idxy2++){
	  if(dilep){
	    if(CheckSignalRegion(isss, GetElectrons(("HNTight_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data()) ,jets,"", w)) FillHist("d0cutflowzoomed", idxy2, w,  0.,55.,55);
	    if(CheckSignalRegion(isss, GetElectrons(("HNTight_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data()) ,jets,"", w)) FillHist("d0cutflowzoomed_now", idxy2, 1.,  0.,55.,55);
	    if(CheckSignalRegionNN(isss, GetElectrons(("HNTight_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data()) ,jets,"", w)) FillHist("NNd0cutflowzoomed", idxy2, w,  0.,55.,55);
            if(CheckSignalRegionNN(isss, GetElectrons(("HNTight_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data()) ,jets,"", w)) FillHist("NNd0cutflowzoomed_now", idxy2, 1.,  0.,55.,55);
	  }
	  else{
	    if(CheckSignalRegionTriLep(trilepmuons,  GetElectrons(("HNTight_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data()) ,jets,"", w)) FillHist("d0cutflowzoomed_trilep", idxy2, w,  0.,55.,55);
            if(CheckSignalRegionTriLep(trilepmuons,  GetElectrons(("HNTight_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data()) ,jets,"", w)) FillHist("d0cutflowzoomed_now_trilep", idxy2, 1.,  0.,55.,55);
	  }
	}
      }
    }
    
    int idxy(1);
    for(unsigned int ib=0; ib <dxyb.size(); ib++){
      bool passlast=true;
      for(unsigned int ie=0; ie <dxye.size(); ie++, idxy++){
	if(dilep){
	  if(CheckSignalRegion(isss, GetElectrons(("HNTight_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data()) ,jets,"", w)) FillHist("d0cutflow", idxy, w,  0.,30.,30);
	  if(CheckSignalRegion(isss, GetElectrons(("HNTight_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data()) ,jets,"", w)) FillHist("d0cutflow_now", idxy, 1.,  0.,30.,30);
	  if(CheckSignalRegionNN(isss, GetElectrons(("HNTight_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data()) ,jets,"", w)) FillHist("NNd0cutflow", idxy, w,  0.,30.,30);
          if(CheckSignalRegionNN(isss, GetElectrons(("HNTight_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data()) ,jets,"", w)) FillHist("NNd0cutflow_now", idxy, 1.,  0.,30.,30);
	}
	else{
	  if(CheckSignalRegionTriLep(trilepmuons,  GetElectrons(("HNTight_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data()) ,jets,"", w)) FillHist("d0cutflow_trilep", idxy, w,  0.,30.,30);
          if(CheckSignalRegionTriLep(trilepmuons, GetElectrons(("HNTight_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data()) ,jets,"", w)) FillHist("d0cutflow_now_trilep", idxy, 1.,  0.,30.,30);
	}
      }
    }
  }
  else{
    std::vector<TString> dxyb;
    dxyb.push_back("050");
    dxyb.push_back("025");
    dxyb.push_back("015");
    dxyb.push_back("010");
    std::vector<TString> dxye;
    dxye.push_back("100");
    dxye.push_back("050");
    dxye.push_back("040");
    dxye.push_back("025");
    dxye.push_back("020");
    dxye.push_back("010");

    std::vector<TString> dxyb2;
    dxyb2.push_back("017");
    dxyb2.push_back("015");
    dxyb2.push_back("014");
    dxyb2.push_back("013");
    dxyb2.push_back("012");
    dxyb2.push_back("011");
    dxyb2.push_back("010");
    std::vector<TString> dxye2;
    dxye2.push_back("017");
    dxye2.push_back("015");
    dxye2.push_back("014");
    dxye2.push_back("013");
    dxye2.push_back("012");
    dxye2.push_back("011");
    dxye2.push_back("010");

    if(1){
      int idxy2(1);
      for(unsigned int ib=0; ib <dxyb2.size(); ib++){
        for(unsigned int ie=0; ie <dxye2.size(); ie++, idxy2++){
	  
	  std::vector<snu::KElectron>  eltight = GetElectrons(("HNTight_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data());

	  TString idtag="";
	  TString looseid="ELECTRON16_HN_FAKELOOSE_NOD0";
	  if(!removed0){
	    looseid="ELECTRON16_HN_FAKELOOSE_D0";
	    idtag="_d0";
	  }

	  std::vector<snu::KElectron>  elloosetmp=  GetElectrons(looseid);
	  std::vector<snu::KElectron>  elloose;
          if(!removed0){
            for(unsigned int iel=0; iel < elloosetmp.size(); iel++){
              if(fabs(elloosetmp[iel].dxy()) < GetDXYCut(elloosetmp[iel],"b"+dxyb2[ib]+"_e"+dxye2[ie])){
               elloose.push_back(elloosetmp[iel]);
	      } 
            }
          }
          else elloose = elloosetmp;


	  if(dilep){
	    
	    float evw = m_datadriven_bkg->Get_DataDrivenWeight_EE(false,  elloose, looseid,("HNTight_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data(), "opt_dijet_ajet40"+idtag);
	    
	    
	    if(CheckSignalRegion(isss, elloose ,jets,"", evw)) FillHist("d0cutflowzoomed"+idtag, idxy2, evw,  0.,55.,55);
	    if(CheckSignalRegion(isss, elloose  ,jets,"", evw)){
	      FillHist("d0cutflowzoomed_now"+idtag, idxy2, 1.,  0.,55.,55);
	      if(PassID(elloose[0], "HNTight_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie])  &&PassID(elloose[1], "HNTight_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]) )               FillHist("TTd0cutflowzoomed_now"+idtag, idxy2, 1.,  0.,55.,55);
	      if(PassID(elloose[0], "HNTight_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie])  &&!PassID(elloose[1], "HNTight_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]) )               FillHist("TLd0cutflowzoomed_now"+idtag, idxy2, 1.,  0.,55.,55);
	      if(!PassID(elloose[0], "HNTight_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie])  &&PassID(elloose[1], "HNTight_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]) )               FillHist("LTd0cutflowzoomed_now"+idtag, idxy2, 1.,  0.,55.,55);
	      if(!PassID(elloose[0], "HNTight_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie])  &&!PassID(elloose[1], "HNTight_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]) )               FillHist("LLd0cutflowzoomed_now"+idtag, idxy2, 1.,  0.,55.,55);
	      
	    }
	    if(CheckSignalRegionNN(isss, elloose ,jets,"", evw)) FillHist("NNd0cutflowzoomed"+idtag, idxy2, evw,  0.,55.,55);
            if(CheckSignalRegionNN(isss, elloose  ,jets,"", evw)) FillHist("NNd0cutflowzoomed_now"+idtag, idxy2, 1.,  0.,55.,55);

	  }
	  else{

	    float mme_evw =   m_datadriven_bkg->Get_DataDrivenWeight(false, trilepmuons, "MUON_HN_TRI_TIGHT", 2, elloose, ("HNTight_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data(),1, looseid, "opt_dijet_ajet40"+idtag);
	    if(CheckSignalRegionTriLep(trilepmuons, elloose ,jets,"", mme_evw)) FillHist("d0cutflowzoomed_trilep"+idtag, idxy2, mme_evw,  0.,55.,55);
            if(CheckSignalRegionTriLep(trilepmuons, elloose  ,jets,"", mme_evw)) FillHist("d0cutflowzoomed_now_trilep"+idtag, idxy2, 1.,  0.,55.,55);

	  }
	}
      }
    }


    int idxy(1);
    for(unsigned int ib=0; ib <dxyb.size(); ib++){
      for(unsigned int ie=0; ie <dxye.size(); ie++, idxy++){
	std::vector<snu::KElectron>  eltight = GetElectrons(("HNTight_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data());

	TString idtag="";
	TString looseid="ELECTRON16_HN_FAKELOOSE_NOD0";
	if(!removed0){
	  looseid="ELECTRON16_HN_FAKELOOSE_D0";
	  idtag="_d0";
	}

	std::vector<snu::KElectron>  elloosetmp=  GetElectrons(looseid);
	std::vector<snu::KElectron>  elloose;
	if(!removed0){
	  for(unsigned int iel=0; iel < elloosetmp.size(); iel++){
	    if(fabs(elloosetmp[iel].dxy()) < GetDXYCut(elloosetmp[iel],"b"+dxyb[ib]+"_e"+dxye[ie])){
	     elloose.push_back(elloosetmp[iel]);
	    }
	  }
	}
	else elloose = elloosetmp;

	
	if(dilep){
	  float evw = m_datadriven_bkg->Get_DataDrivenWeight_EE(false,  elloose, looseid,("HNTight_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data(), "opt_dijet_ajet40"+idtag);
	  
	  if(CheckSignalRegion(isss, elloose ,jets,"", evw)) FillHist("d0cutflow"+idtag, idxy, evw,  0.,30.,30);
	  if(CheckSignalRegion(isss, elloose  ,jets,"", evw)) FillHist("d0cutflow_now"+idtag, idxy, 1.,  0.,30.,30);
	  if(CheckSignalRegionNN(isss, elloose ,jets,"", evw)) FillHist("NNd0cutflow"+idtag, idxy, evw,  0.,30.,30);
          if(CheckSignalRegionNN(isss, elloose  ,jets,"", evw)) FillHist("NNd0cutflow_now"+idtag, idxy, 1.,  0.,30.,30);
	}
	else{
	  float mme_evw =   m_datadriven_bkg->Get_DataDrivenWeight(false, trilepmuons, "MUON_HN_TRI_TIGHT", 2, elloose,("HNTight_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data(),  1, looseid,"opt_dijet_ajet40"+idtag);
	  if(CheckSignalRegionTriLep(trilepmuons, elloose ,jets,"", mme_evw)) FillHist("d0cutflow_trilep"+idtag, idxy, mme_evw,  0.,30.,30);
          if(CheckSignalRegionTriLep(trilepmuons, elloose  ,jets,"", mme_evw)) FillHist("d0cutflow_now_trilep"+idtag, idxy, 1.,  0.,30.,30);
	}
      }
    }
  }


  if(dilep){
    if(CheckSignalRegion(isss,GetElectrons("HNTight_dxy_ref"),jets,"", w)) FillHist("d0sigcutflow", 0., w, 0.,30.,30);
    if(CheckSignalRegionNN(isss,GetElectrons("HNTight_dxy_ref"),jets,"", w)) FillHist("NNd0sigcutflow", 0., w, 0.,30.,30);
  }
  else{
    if(CheckSignalRegionTriLep(trilepmuons, GetElectrons("HNTight_dxy_ref"),jets,"", w)) FillHist("d0sigcutflow_trilep", 0., w, 0.,30.,30);
  }
  
  if(!k_running_nonprompt){
    
    std::vector<TString> dxyb;
    dxyb.push_back("050");
    dxyb.push_back("025");
    dxyb.push_back("015");
    dxyb.push_back("010");
    std::vector<TString> dxye;
    dxye.push_back("100");
    dxye.push_back("050");
    dxye.push_back("040");
    dxye.push_back("025");
    dxye.push_back("020");
    dxye.push_back("010");

    std::vector<TString> dxyb2;
    dxyb2.push_back("017");
    dxyb2.push_back("015");
    dxyb2.push_back("014");
    dxyb2.push_back("013");
    dxyb2.push_back("012");
    dxyb2.push_back("011");
    dxyb2.push_back("010");
    std::vector<TString> dxye2;
    dxye2.push_back("017");
    dxye2.push_back("015");
    dxye2.push_back("014");
    dxye2.push_back("013");
    dxye2.push_back("012");
    dxye2.push_back("011");
    dxye2.push_back("010");
    
    if(1){
      int idxy2(1);
      for(unsigned int ib=0; ib <dxyb2.size(); ib++){
        for(unsigned int ie=0; ie <dxye2.size(); ie++, idxy2++){

	  if(dilep){
	    if(CheckSignalRegion(isss, GetElectrons(("HNTight_dxysig_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data()) ,jets,"", w)) FillHist("d0sigcutflowzoomed", idxy2, w,  0.,55.,55);
	    if(CheckSignalRegion(isss, GetElectrons(("HNTight_dxysig_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data()) ,jets,"", w)) FillHist("d0sigcutflowzoomed_now", idxy2, 1.,  0.,55.,55);

	    if(CheckSignalRegionNN(isss, GetElectrons(("HNTight_dxysig_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data()) ,jets,"", w)) FillHist("NNd0sigcutflowzoomed", idxy2, w,  0.,55.,55);
																     
            if(CheckSignalRegionNN(isss, GetElectrons(("HNTight_dxysig_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data()) ,jets,"", w)) FillHist("NNd0sigcutflowzoomed_now", idxy2, 1.,  0.,55.,55);																     
																     

	  }
	  else{
	    if(CheckSignalRegionTriLep(trilepmuons, GetElectrons(("HNTight_dxysig_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data()) ,jets,"", w)) FillHist("d0sigcutflowzoomed_trilep", idxy2, w,  0.,55.,55);
							       								     
            if(CheckSignalRegionTriLep(trilepmuons, GetElectrons(("HNTight_dxysig_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data()) ,jets,"", w)) FillHist("d0sigcutflowzoomed_trilep_now", idxy2, 1.,  0.,55.,55);
	    

	  }
        }
      }
    }



    int idxy(1);
    for(unsigned int ib=0; ib <dxyb.size(); ib++){
      for(unsigned int ie=0; ie <dxye.size(); ie++, idxy++){
	if(dilep){
	  if(CheckSignalRegion(isss, GetElectrons(("HNTight_dxysig_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data()) ,jets,"", w)) FillHist("d0sigcutflow", idxy, w,  0.,30.,30);
	  if(CheckSignalRegion(isss, GetElectrons(("HNTight_dxysig_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data()) ,jets,"", w)) FillHist("d0sigcutflow_now", idxy, 1.,  0.,30.,30);
	  if(CheckSignalRegionNN(isss, GetElectrons(("HNTight_dxysig_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data()) ,jets,"", w)) FillHist("NNd0sigcutflow", idxy, w,  0.,30.,30);
          if(CheckSignalRegionNN(isss, GetElectrons(("HNTight_dxysig_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data()) ,jets,"", w)) FillHist("NNd0sigcutflow_now", idxy, 1.,  0.,30.,30);

	}
	else{
	  if(CheckSignalRegionTriLep(trilepmuons, GetElectrons(("HNTight_dxysig_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data()) ,jets,"", w)) FillHist("d0sigcutflow_trilep", idxy, w,  0.,30.,30);
          if(CheckSignalRegionTriLep(trilepmuons, GetElectrons(("HNTight_dxysig_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data()) ,jets,"", w)) FillHist("d0sigcutflow_now_trilep", idxy, 1.,  0.,30.,30);

	}
      }
    }
  }
  else{
    std::vector<TString> dxyb;
    dxyb.push_back("050");
    dxyb.push_back("025");
    dxyb.push_back("015");
    dxyb.push_back("010");
    std::vector<TString> dxye;
    dxye.push_back("100");
    dxye.push_back("050");
    dxye.push_back("040");
    dxye.push_back("025");
    dxye.push_back("020");
    dxye.push_back("010");

    std::vector<TString> dxyb2;
    dxyb2.push_back("017");
    dxyb2.push_back("015");
    dxyb2.push_back("014");
    dxyb2.push_back("013");
    dxyb2.push_back("012");
    dxyb2.push_back("011");
    dxyb2.push_back("010");
    std::vector<TString> dxye2;
    dxye2.push_back("017");
    dxye2.push_back("015");
    dxye2.push_back("014");
    dxye2.push_back("013");
    dxye2.push_back("012");
    dxye2.push_back("011");
    dxye2.push_back("010");

    if(1){
      int idxy2(1);
      for(unsigned int ib=0; ib <dxyb2.size(); ib++){
        for(unsigned int ie=0; ie <dxye2.size(); ie++, idxy2++){

	  std::vector<snu::KElectron>  eltight = GetElectrons(("HNTight_dxysig_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data());

	  TString idtag="";
	  TString looseid="ELECTRON16_HN_FAKELOOSE_NOD0";
	  if(!removed0){
	    looseid="ELECTRON16_HN_FAKELOOSE_D0";
	    idtag="_d0";
	  }

	  std::vector<snu::KElectron>  elloosetmp=  GetElectrons(looseid);
	  std::vector<snu::KElectron>  elloose;
          if(!removed0){
            for(unsigned int iel=0; iel < elloosetmp.size(); iel++){
              if(fabs(elloosetmp[iel].dxy()) < GetDXYCut(elloosetmp[iel],"b"+dxyb2[ib]+"_e"+dxye2[ie])){
               elloose.push_back(elloosetmp[iel]);
              }
            }
          }
	  else elloose = elloosetmp;

	  if(dilep){
	    float evw = m_datadriven_bkg->Get_DataDrivenWeight_EE(false,  elloose, looseid,("HNTight_dxysig_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data(), "opt_dijet_ajet40"+idtag);
	    
	    if(CheckSignalRegion(isss, elloose ,jets,"", evw)) FillHist("d0sigcutflowzoomed"+idtag, idxy2, evw,  0.,55.,55);
	    if(CheckSignalRegion(isss, elloose  ,jets,"", evw)) FillHist("d0sigcutflowzoomed_now"+idtag, idxy2, 1.,  0.,55.,55);
	    if(CheckSignalRegionNN(isss, elloose ,jets,"", evw)) FillHist("NNd0sigcutflowzoomed"+idtag, idxy2, evw,  0.,55.,55);
            if(CheckSignalRegionNN(isss, elloose  ,jets,"", evw)) FillHist("NNd0sigcutflowzoomed_now"+idtag, idxy2, 1.,  0.,55.,55);
	  }
	  else{
	    float mme_evw =   m_datadriven_bkg->Get_DataDrivenWeight(false, trilepmuons, "MUON_HN_TRI_TIGHT", 2, elloose, ("HNTight_dxysig_dxy_b"+dxyb2[ib]+"_e"+dxye2[ie]).Data(), 1, looseid,"opt_dijet_ajet40"+idtag);
	    if(CheckSignalRegionTriLep(trilepmuons, elloose ,jets,"", mme_evw)) FillHist("d0sigcutflowzoomed_trilep"+idtag, idxy2, mme_evw,  0.,55.,55);
            if(CheckSignalRegionTriLep(trilepmuons, elloose ,jets,"", mme_evw)) FillHist("d0sigcutflowzoomed_now_trilep"+idtag, idxy2, 1.,  0.,55.,55);

	  }
        }
      }
    }


    int idxy(1);
    for(unsigned int ib=0; ib <dxyb.size(); ib++){
      for(unsigned int ie=0; ie <dxye.size(); ie++, idxy++){
	std::vector<snu::KElectron>  eltight = GetElectrons(("HNTight_dxysig_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data());

	TString idtag="";
	TString looseid="ELECTRON16_HN_FAKELOOSE_NOD0";
	if(!removed0){
	  looseid="ELECTRON16_HN_FAKELOOSE_D0";
	  idtag="_d0";
	}

	std::vector<snu::KElectron>  elloosetmp=  GetElectrons(looseid);
	std::vector<snu::KElectron>  elloose;
	if(!removed0){
	  for(unsigned int iel=0; iel < elloosetmp.size(); iel++){
	    if(fabs(elloosetmp[iel].dxy()) < GetDXYCut(elloosetmp[iel],"b"+dxyb[ib]+"_e"+dxye[ie])){
	     elloose.push_back(elloosetmp[iel]);
	    }
	  }
	}
	else elloose = elloosetmp;

								 
	if(dilep){
	  float evw = m_datadriven_bkg->Get_DataDrivenWeight_EE(false,  elloose, looseid,("HNTight_dxysig_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data(), "opt_dijet_ajet40"+idtag);

	  if(CheckSignalRegion(isss, elloose ,jets,"", w))FillHist("d0sigcutflow"+idtag, idxy, evw,  0.,30.,30);
	  if(CheckSignalRegion(isss, elloose ,jets,"", w))FillHist("d0sigcutflow_now"+idtag, idxy, 1.,  0.,30.,30);

          if(CheckSignalRegionNN(isss, elloose ,jets,"", w))FillHist("NNd0sigcutflow"+idtag, idxy, evw,  0.,30.,30);
          if(CheckSignalRegionNN(isss, elloose ,jets,"", w))FillHist("NNd0sigcutflow_now"+idtag, idxy, 1.,  0.,30.,30);

	}
	else{
	  float mme_evw =   m_datadriven_bkg->Get_DataDrivenWeight(false, trilepmuons, "MUON_HN_TRI_TIGHT", 2, elloose,("HNTight_dxysig_dxy_b"+dxyb[ib]+"_e"+dxye[ie]).Data(),1, looseid, "opt_dijet_ajet40"+idtag);
	  if(CheckSignalRegionTriLep(trilepmuons, elloose ,jets,"", w))FillHist("d0sigcutflow_trilep"+idtag, idxy, mme_evw,  0.,30.,30);
          if(CheckSignalRegionTriLep(trilepmuons, elloose ,jets,"", w))FillHist("d0sigcutflow_now_trilep"+idtag, idxy, 1.,  0.,30.,30);

	}
      }
    }
  }
  
  
}

float  HNDiElectronOptimisation::GetTightWeight(){
  return 1.;
}

float HNDiElectronOptimisation::GetMediumWeight(){
  
  return 1.;
}

void HNDiElectronOptimisation::GetSSSignalEfficiency(float w){
  
  // w = lumi weight * MC gen weight


  /// check truth matching eff.


  //// Check efficiency of Selecting two same sign electrons
  if(SameCharge(GetElectrons("ELECTRON_NOCUT")))   FillCutFlow("SS_NoCut",w);
  if(SameCharge(GetElectrons("ELECTRON_PTETA")))   FillCutFlow("SS_PtEta",w);
  // TIGHT POG
  if(SameCharge(GetElectrons("ELECTRON16_POG_TIGHT_SC")))    FillCutFlow("SS_Tight",w*GetTightWeight());
  if(SameCharge(GetElectrons("ELECTRON16_POG_TIGHT_D0"))) FillCutFlow("SS_Tight_d0veto",w*GetTightWeight());
  if(SameCharge(GetElectrons("ELECTRON16_POG_TIGHT_RELISO")))     FillCutFlow("SS_Tight_reliso", w*GetTightWeight());  
  if(SameCharge(GetElectrons("ELECTRON16_POG_TIGHT_CHARGECONS"))) FillCutFlow("SS_Tight_chargeconst", w*GetTightWeight());
  if(SameCharge(GetElectrons("ELECTRON16_POG_TIGHT_DXYSIG")))     FillCutFlow("SS_Tight_dxysig",w*GetTightWeight());

  TString analysis_trigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
  /// Trigger List (unprescaled)                                                                                                                                                                                                                                                
  std::vector<TString> triggerslist;
  triggerslist.push_back(analysis_trigger);

  if(PassTrigger(analysis_trigger) && SameCharge(GetElectrons("ELECTRON16_POG_TIGHT_DXYSIG")))  FillCutFlow("SS_Tight_trigger",w*GetTightWeight());



  // MEDIUM POG
  if(SameCharge(GetElectrons("ELECTRON16_POG_MEDIUM")))    FillCutFlow("SS_Medium",w*GetMediumWeight());
  if(SameCharge(GetElectrons("ELECTRON16_POG_MEDIUM_D0"))) FillCutFlow("SS_Medium_d0veto",w*GetMediumWeight());
  if(SameCharge(GetElectrons("ELECTRON16_POG_MEDIUM_RELISO")))     FillCutFlow("SS_Medium_reliso", w*GetMediumWeight());
  if(SameCharge(GetElectrons("ELECTRON16_POG_MEDIUM_CHARGECONS"))) FillCutFlow("SS_Medium_chargeconst", w*GetMediumWeight());
  if(SameCharge(GetElectrons("ELECTRON16_POG_MEDIUM_DXYSIG")))     FillCutFlow("SS_Medium_dxysig",w*GetMediumWeight());

  if(PassTrigger(analysis_trigger) && SameCharge(GetElectrons("ELECTRON16_POG_MEDIUM_DXYSIG")))  FillCutFlow("SS_Medium_trigger",w*GetMediumWeight());

  /// HN ID
  if(SameCharge(GetElectrons("ELECTRON_HN_TIGHT"))) FillCutFlow("SS_anal_el", w);
  
  std::vector<snu::KJet> jets = GetJets("JET_HN");
  std::vector<snu::KJet> jets_pu = GetJets("JET_HN_PU");
  if(SameCharge(GetElectrons("ELECTRON_HN_TIGHT"))){
    if(jets.size() > 0)  FillCutFlow("ss_anal_el_jet", w*GetTightWeight());
    if(jets.size() > 1)  FillCutFlow("ss_anal_el_dijet", w*GetTightWeight());
    if(PassTrigger(analysis_trigger)){
      if(jets.size() > 1)  FillCutFlow("ss_trig_anal_el_dijet", w);
    }
    if(jets_pu.size() > 0)  FillCutFlow("ss_anal_el_pujet", w*GetTightWeight());
    if(jets_pu.size() > 1)  FillCutFlow("ss_anal_el_pudijet", w*GetTightWeight());
    if(PassTrigger(analysis_trigger)){
      if(jets_pu.size() > 1) FillCutFlow("ss_trig_anal_el_pudijet", w*GetTightWeight());
    }
  }
  
}

void HNDiElectronOptimisation::GetOSSignalEfficiency(float w){
  
  // w = lumi weight * MC gen weight                                                                                                                                                                                                                                            
  //// Check efficiency of Selecting two same sign electrons                                                                                                                                                                                                                    
  
  if(OppositeCharge(GetElectrons("ELECTRON_NOCUT")))   FillCutFlow("OS_NoCut",w);
  if(OppositeCharge(GetElectrons("ELECTRON_PTETA")))   FillCutFlow("OS_PtEta",w);
  // TIGHT POG                                                                                                                                                                                                                                                                  
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_TIGHT_SC")))    FillCutFlow("OS_Tight",w*GetTightWeight());
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_TIGHT_D0"))) FillCutFlow("OS_Tight_d0veto",w*GetTightWeight());
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_TIGHT_RELISO")))     FillCutFlow("OS_Tight_reliso", w*GetTightWeight());
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_TIGHT_CHARGECONS"))) FillCutFlow("OS_Tight_chargeconst", w*GetTightWeight());
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_TIGHT_DXYSIG")))     FillCutFlow("OS_Tight_dxysig",w*GetTightWeight());

  TString analysis_trigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
  /// Trigger List (unprescaled)                                                                                                                                                                                                                                             
  std::vector<TString> triggerslist;
  triggerslist.push_back(analysis_trigger);
  
  if(PassTrigger(analysis_trigger) && OppositeCharge(GetElectrons("ELECTRON16_POG_TIGHT_DXYSIG")))  FillCutFlow("OS_Tight_trigger",w*GetTightWeight());


  // MEDIUM POG                                                                                                                                                                                                                                                                 
  
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_MEDIUM")))    FillCutFlow("OS_Medium",w*GetMediumWeight());
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_MEDIUM_D0"))) FillCutFlow("OS_Medium_d0veto",w*GetMediumWeight());
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_MEDIUM_RELISO")))     FillCutFlow("OS_Medium_reliso", w*GetMediumWeight());
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_MEDIUM_CHARGECONS"))) FillCutFlow("OS_Medium_chargeconst", w*GetMediumWeight());
  if(OppositeCharge(GetElectrons("ELECTRON16_POG_MEDIUM_DXYSIG")))     FillCutFlow("OS_Medium_dxysig",w*GetMediumWeight());

  if(PassTrigger(analysis_trigger) && OppositeCharge(GetElectrons("ELECTRON16_POG_MEDIUM_DXYSIG")))  FillCutFlow("OS_Medium_trigger",w*GetMediumWeight());
  /// HN ID                                                                                                                                                                                                                                                                     
  if(OppositeCharge(GetElectrons("ELECTRON_HN_TIGHT"))) FillCutFlow("OS_anal_el", w);
  
  std::vector<snu::KJet> jets = GetJets("JET_HN");
  std::vector<snu::KJet> jets_pu = GetJets("JET_HN_PU");
  if(OppositeCharge(GetElectrons("ELECTRON_HN_TIGHT"))){
    if(jets.size() > 0)  FillCutFlow("os_anal_el_jet", w);
    if(jets.size() > 1)  FillCutFlow("os_anal_el_dijet", w);
    if(PassTrigger(analysis_trigger)){
	if(jets.size() > 1)  FillCutFlow("os_trig_anal_el_dijet", w);
    }
    if(jets_pu.size() > 0)  FillCutFlow("os_anal_el_pujet", w);
    if(jets_pu.size() > 1)  FillCutFlow("os_anal_el_pudijet", w);
    if(PassTrigger(analysis_trigger)){
      if(jets_pu.size() > 1)  FillCutFlow("os_trig_anal_el_pudijet", w);
    }
  }
}
										   
								   

  

void HNDiElectronOptimisation::SignalValidation(){
  
  // Check jet properties                                                                                                                                                                                                                                                                                                                                                         
  FillCLHist(sighist_ee, "SIGNALVALIDATION_EE", eventbase->GetEvent(),  GetMuons("MUON_NOCUT"), GetElectrons("ELECTRON_NOCUT"),GetJets("JET_HN"), weight);
  FillCLHist(sighist_ee, "SIGNALVALIDATION_TIGHT_EE", eventbase->GetEvent(),  GetMuons("MUON_NOCUT"), GetElectrons("ELECTRON_POG_TIGHT"),GetJets("JET_HN"), weight);
  FillCLHist(sighist_mm, "SIGNALVALIDATION_MM", eventbase->GetEvent(),  GetMuons("MUON_NOCUT"), GetElectrons("ELECTRON_NOCUT"),GetJets("JET_HN"), weight);
  if(GetJets("JET_HN").size() > 1){
    if(GetElectrons("ELECTRON_POG_TIGHT").size() == 2){
      FillCLHist(sighist_ee, "SIGNALVALIDATION_EE_PRESEL", eventbase->GetEvent(),  GetMuons("MUON_NOCUT"), GetElectrons("ELECTRON_POG_TIGHT"),GetJets("JET_HN"), weight);
    }
    if(GetMuons("MUON_POG_TIGHT").size() == 2){
      FillCLHist(sighist_mm, "SIGNALVALIDATION_MM_PRESEL", eventbase->GetEvent(),  GetMuons("MUON_POG_TIGHT"), GetElectrons("ELECTRON_NOCUT"),GetJets("JET_HN"), weight);
    }
  }
  

}
 

void HNDiElectronOptimisation::GetTriggEfficiency(){
  //ListTriggersAvailable();                                                                                                                                                                                                                                                                                                                                                      
  vector<int> pt1;
  pt1.push_back(18);
  pt1.push_back(35);
  pt1.push_back(25);
  pt1.push_back(25);
  pt1.push_back(25);
  pt1.push_back(30);
  pt1.push_back(30);
  pt1.push_back(20);
  pt1.push_back(25);
  pt1.push_back(30);
  pt1.push_back(120);
  pt1.push_back(20);
  pt1.push_back(25);
  pt1.push_back(25);
  pt1.push_back(30);
  pt1.push_back(35);
  pt1.push_back(28);
  vector<int>pt2;
  pt2.push_back(10);
  pt2.push_back(35);
  pt2.push_back(25);
  pt2.push_back(10);
  pt2.push_back(10);
  pt2.push_back(10);
  pt2.push_back(10);
  pt2.push_back(15);
  pt2.push_back(15);
  pt2.push_back(10);
  pt2.push_back(10);
  pt2.push_back(15);
  pt2.push_back(15);
  pt2.push_back(15);
  pt2.push_back(10);
  pt2.push_back(10);
  pt2.push_back(10);


  std::vector<TString> lists_triggers;
  lists_triggers.push_back("HLT_Ele15_IsoVVVL_PFHT600_v");
  lists_triggers.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v");
  lists_triggers.push_back("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v");
  lists_triggers.push_back("HLT_Ele22_eta2p1_WPLoose_Gsf_v");
  lists_triggers.push_back("HLT_Ele23_WPLoose_Gsf_v");
  lists_triggers.push_back("HLT_Ele27_WPLoose_Gsf_v");
  lists_triggers.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_v");
  lists_triggers.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v");
  lists_triggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");
  lists_triggers.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_HT200");
  lists_triggers.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");
  lists_triggers.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  lists_triggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  lists_triggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded_v1");
  
  lists_triggers.push_back("HLT_Ele27_WPTight_Gsf_v");
  lists_triggers.push_back("HLT_Ele32_eta2p1_WPTight_Gsf_v");
  lists_triggers.push_back("HLT_Ele25_eta2p1_WPTight_Gsf_v");

  std::vector<TString> lists_triggerORs;
  lists_triggerORs.push_back("HLT_Ele15_IsoVVVL_PFHT600_v");
  lists_triggerORs.push_back("HLT_Ele25_eta2p1_WPTight_Gsf_v");
  lists_triggerORs.push_back("HLT_Ele27_WPTight_Gsf_v");
  lists_triggerORs.push_back("HLT_Ele32_eta2p1_WPTight_Gsf");

  for(unsigned int i=0; i < lists_triggers.size(); i++){
    FillTriggerEfficiency(lists_triggers.at(i), weight, "denominator_nojet", lists_triggers );
  }

  std::vector<snu::KJet> jets=GetJets("JET_HN");
  std::vector<snu::KElectron> electrons= GetElectrons("ELECTRON_POG_TIGHT");
  if(jets.size() > 1){
    for(unsigned int i=0; i < lists_triggers.size(); i++){
      TString trig=lists_triggers.at(i);
      FillTriggerEfficiency(lists_triggers.at(i), weight, "denominator", lists_triggers );
      if(PassTrigger(trig))  {
	
	FillTriggerEfficiency(lists_triggers.at(i), weight, "numerator",lists_triggers );
	
	if(electrons.size() ==2) {
	  
	  FillTriggerEfficiency(lists_triggers.at(i), weight, "numerator_dimuon",lists_triggers );
	  if(electrons.at(0).Pt() > pt1.at(i) && electrons.at(1).Pt() > pt2.at(i))  FillTriggerEfficiency(lists_triggers.at(i), weight, "numerator_dimuon_pt",lists_triggers );
	}
      }
      for(int itor = 0; itor < lists_triggerORs.size(); itor++){

	if(PassTrigger(trig) || PassTrigger(lists_triggerORs.at(itor))){
	  FillTriggerEfficiency(lists_triggers.at(i), weight, "numerator_"+lists_triggerORs.at(itor),lists_triggers );
	  if(electrons.size() ==2) {
	    
	    FillTriggerEfficiency(lists_triggers.at(i), weight, "numerator_dimuon_"+lists_triggerORs.at(itor),lists_triggers );
	    if(electrons.at(0).Pt() > pt1.at(i) && electrons.at(1).Pt() > 5.)  FillTriggerEfficiency(lists_triggers.at(i), weight, "numerator_dimuon_pt5_"+lists_triggerORs.at(itor),lists_triggers );
	    if(electrons.at(0).Pt() > pt1.at(i) && electrons.at(1).Pt() > 10.)  FillTriggerEfficiency(lists_triggers.at(i), weight, "numerator_dimuon_pt10_"+lists_triggerORs.at(itor),lists_triggers );
	  }
	}
      }
    }
  }
}

void HNDiElectronOptimisation::FillTriggerEfficiency(TString cut, float weight, TString label, std::vector<TString> list){

  if(GetHist("TriggerEfficiency_" + label)) {
    GetHist("TriggerEfficiency_"+label)->Fill(cut,weight);

  }
  else{
    int ntrig = list.size();
    AnalyzerCore::MakeHistograms("TriggerEfficiency_"+label,ntrig,0.,float(ntrig));

    for(unsigned int it=0; it < list.size(); it++){
      GetHist("TriggerEfficiency_"+label)->GetXaxis()->SetBinLabel(it+1,list.at(it));
    }
  }
  
}

bool HNDiElectronOptimisation::CheckSignalRegionTriLep(  std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, TString name, float w){
  
  
  if(electrons.size() != 1 ) return false ;
  if(muons.size() != 2 ) return false ;
  // Set by trigger                                                                                                                                                                                                                        
  if(electrons.at(0).Pt() < 20.)  {cout << "Faile " << endl; return false;}
  if(muons.at(0).Pt() < 20.) return false;
  if(muons.at(1).Pt() < 15.) return false;


  
  snu::KParticle emu = electrons.at(0) + muons.at(0);
  snu::KParticle emu2 = electrons.at(0) + muons.at(1);
  snu::KParticle mm = muons.at(0) + muons.at(1);
  if(emu.Charge() ==0){
    if(emu.M() < 4.) return false;
    if(fabs(emu.M() - 91.1) < 15.) return false;
  }
  if(emu2.Charge() ==0){
    if(emu2.M() < 4.) return false;
    if(fabs(emu2.M() - 91.1) < 15.) return false;
  }
  if(mm.Charge() ==0){

    if(mm.M() < 4.) return false;
    if(fabs(mm.M() - 91.1) < 15.) return false;
  }
  
  snu::KParticle emm = muons.at(0) + muons.at(1) + electrons.at(0);
  
  if(fabs(emm.M() - 91.1) < 15.) return false;

  if(eventbase->GetEvent().PFMET() >  40.) return false;

  int nbjet=0;
  for(unsigned int ij=0; ij <jets.size(); ij++){
    
    if( jets[ij].IsBTagged(  snu::KJet::CSVv2, snu::KJet::Medium))  nbjet++;
  }
  if(nbjet > 0) return false;
  
  return true;

}


bool HNDiElectronOptimisation::CheckSignalRegion( bool isss,  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, TString name, float w){

  bool debug=false;
  if(electrons.size() != 2 ) {if(debug)cout << "Fail el size" << endl; return false ;}
  // Set by trigger
  if(electrons.at(0).Pt() < 25.) {if(debug)cout << "Fail pt1 " << endl; return false;}
  if(electrons.at(1).Pt() < 10.)  {if(debug)cout << "Fail pt2  " << endl; return false;}

  if(isss&&!SameCharge(electrons)) {if(debug)cout << "Fail ss " << endl; return false;}

  if(!isss&&SameCharge(electrons)) {if(debug)cout << "Fail os  " << endl; return false;}
  if(jets.size() < 2) {if(debug)cout << "Fail jets " << endl; return false;}

  snu::KParticle ee = electrons.at(0) + electrons.at(1);
  if(ee.M()  < 10.) {if(debug)cout << "Fail mee  " << endl; return false;}
  
  if(name.Contains("Low")){
    if(ee.M()  > 70.) {if(debug)cout << "Fail mee  " << endl; return false;}
  }
  
  float dijetmass_tmp=999.;
  float dijetmass=9990000.;
  int m=-999;
  int n=-999;
  for(UInt_t emme=0; emme<jets.size(); emme++){
    for(UInt_t enne=1; enne<jets.size(); enne++) {
      if(emme == enne) continue;
      dijetmass_tmp = (jets[emme]+jets[enne]).M();
      if ( fabs(dijetmass_tmp-80.4) < fabs(dijetmass-80.4) ) {
        dijetmass = dijetmass_tmp;
        m = emme;
        n = enne;
      }
    }
  }
  snu::KParticle jj = jets.at(m) + jets.at(n) ;

  if(name.Contains("Low")){
    if(jj.M() > 150.) {if(debug)cout << "Fail mjj  " << endl; return false;}
  }
  
  else{
    if(jj.M() > 120.) {if(debug)cout << "Fail mjj  " << endl; return false;}
    if(jj.M() < 50.) {if(debug)cout << "Fail mjj  " << endl; return false;}

  }
  
  if(name.Contains("Low")){
    if(electrons.at(0).DeltaR(electrons.at(1)) > 3.5) {if(debug)cout << "Fail dr  " << endl; return false;}
  }
  if(name.Contains("Z")){
    if(!((ee.M() > 80.) && (ee.M() < 100.))) {if(debug)cout << "Fail mZ  " << endl; return false;}
  }
  else{
    if((ee.M() > 80.) && (ee.M() < 100.) ) {if(debug)cout << "Fail mZ  " << endl; return false;}
  }
  float ST = electrons[0].Pt() + electrons[1].Pt();
  float looseST = electrons[0].Pt() + electrons[1].Pt();

  std::vector<snu::KJet> loosejets=GetJets("JET_NOCUT");
  

  float HT=0.;
  for(unsigned int ij=0; ij <jets.size(); ij++){
    ST+= jets[ij].Pt();
    HT+= jets[ij].Pt();

  }
  for(unsigned int ij=0; ij <loosejets.size(); ij++){
    looseST+= loosejets[ij].Pt();
  }

  if(!name.Contains("Low")){
    
    if((pow(eventbase->GetEvent().PFMET(),2.)/ looseST) > 10.) {if(debug)cout << "Fail met  " << endl; return false;}
  }
  else{
    if(eventbase->GetEvent().PFMET() > 60) return false;
  }

  if(name.Contains("Low")){
    if(HT >  200.) return false;
  }
  float dphi_1 = fabs(TVector2::Phi_mpi_pi(electrons.at(0).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet)));
  float MT_1 = sqrt(2.* electrons.at(0).Et()*eventbase->GetEvent().PFMET() * (1 - cos( dphi_1)));

  float dphi_2 = fabs(TVector2::Phi_mpi_pi(electrons.at(0).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet)));
  float MT_2 = sqrt(2.* electrons.at(1).Et()*eventbase->GetEvent().PFMET() * (1 - cos( dphi_2)));

  if(name.Contains("Low")){
    if(MT_1 > 100.) return false;
    if(MT_2 > 80.) return false;
    if(ST > 350.)  return false;
    snu::KParticle lljj = jets.at(m) + jets.at(n) + electrons.at(0) + electrons.at(1);
    if(lljj.M() > 300.) return false;

  }
  else{
    snu::KParticle lljj = jets.at(m) + jets.at(n) + electrons.at(0) + electrons.at(1);
    if(jets.size() > 3){
      snu::KParticle lljjjj = jets.at(0) + jets.at(1) + jets.at(2)+jets.at(3) + electrons.at(0) + electrons.at(1);
      if (lljjjj.M() > 1000.) return false;
    }
    if(ST < 200.)  return false;
    if(lljj.M() < 200.) return false;
    if(electrons.at(0).Pt() < 50.) return false;
    if(electrons.at(1).Pt() < 25.) return false;
    if((electrons.at(0).Pt() + electrons.at(1).Pt() ) < 100.) return false;

  }

  if(name.Contains("High")){
    if((electrons.at(0).Pt() + electrons.at(1).Pt() ) < 150.) return false;
    if (ST < 400. ) return false;

    float mindRjj=999.;
    for(unsigned int ij=0; ij <jets.size()-1; ij++){
      for(unsigned int ij2=ij+1; ij2 <jets.size(); ij2++){
	if(jets.at(ij).DeltaR(jets.at(ij2)) < mindRjj) mindRjj=jets.at(ij).DeltaR(jets.at(ij2)) ;
      }
    }
    //if(mindRjj > 1.5) {if(debug)cout << "Fail dr  " << endl; return false;}
  }

  int nbjet=0;
  for(unsigned int ij=0; ij <jets.size(); ij++){
    if( jets[ij].IsBTagged(  snu::KJet::CSVv2, snu::KJet::Medium))  nbjet++;
  }
  if(nbjet > 0) {if(debug)cout << "Fail nbjet  " << endl; return false;}
  
  if(debug)cout << "PASSES ID" << endl;
  return true;

}

bool HNDiElectronOptimisation::CheckSignalRegionNN( bool isss,  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, TString name, float w){

  if(electrons.size() != 2 ) return false ;
  // Set by trigger                                                                                                                                                             
  if(electrons.at(0).Pt() < 25.) return false;
  if(electrons.at(1).Pt() < 15.) return false;
  if(isss&&!SameCharge(electrons)) return false;

  if(!isss&&SameCharge(electrons)) return false;
  if(jets.size() < 4) return false;
  snu::KParticle ee = electrons.at(0) + electrons.at(1);
  if(ee.M()  < 10.) return false;

  if((ee.M() > 80.) && (ee.M() < 100.) ) return false;

  int nbjet=0;
  for(unsigned int ij=0; ij <jets.size(); ij++){
    if( jets[ij].IsBTagged(  snu::KJet::CSVv2, snu::KJet::Medium))  nbjet++;

  }
  if(nbjet > 0) return false;

  return true;

}



bool HNDiElectronOptimisation::OptMassCheckSignalRegion(std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, float evmet, bool runchargeflip , std::vector<float> cuts , TString additional_option){
  
  return OptMassCheckSignalRegion(electrons, jets, 1, evmet, runchargeflip, cuts, additional_option);
}

bool HNDiElectronOptimisation::OptMassCheckSignalRegion(std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, int btagwp,  float evmet, bool runchargeflip , std::vector<float> cuts, TString additional_option){
  
  if(additional_option.Contains("btag_t")) btagwp = 2; 
  if(additional_option.Contains("btag_l")) btagwp = 0; 

  bool removeclosejet = false;
  if(additional_option.Contains("closejet")) removeclosejet = true;
  
  bool cut_on_me2jj=true;
  if(additional_option.Contains("noMe2jj")) cut_on_me2jj= false;
  
  bool cut_on_ee=true;
  if(additional_option.Contains("noeeupper")) cut_on_ee= false;

  bool cut_on_eejj=true;
  if(additional_option.Contains("noeejjupper")) cut_on_eejj= false;

  bool cut_on_pt=true;
  if(additional_option.Contains("pt10")) cut_on_pt= false;

  
  int nbjet=0;
  snu::KJet::WORKING_POINT wp;
  if(btagwp==0) wp=snu::KJet::Loose;
  if(btagwp==1) wp=snu::KJet::Medium;
  if(btagwp==2) wp=snu::KJet::Tight;

  for(unsigned int ij=0; ij <jets.size(); ij++){
    if(IsBTagged(jets.at(ij), snu::KJet::CSVv2, snu::KJet::Medium, GetPeriod())) nbjet++;
  }

  if(additional_option.Contains("lowmass")) return (nbjet == 0 && (evmet < 30.) && LowMassCheckSignalRegion( electrons,  jets, runchargeflip));
  if(additional_option.Contains("midmass"))  return (nbjet == 0 && (evmet < 35.) && MidMassCheckSignalRegion( electrons,  jets, runchargeflip));
  if(additional_option.Contains("highmass")) return (nbjet == 0 && (evmet < 35.) && HighMassCheckSignalRegion( electrons,  jets, runchargeflip));
  

  if(additional_option.Contains("basic")){
    if(electrons.size() != 2) return false;
    if(electrons.at(0).Pt() < 20.) return false;
    if(electrons.at(1).Pt() < 15.) return false;
    if(jets.size() < 2) return false;
    if(jets.at(0).Pt() < cuts.at(9)) return false;
    if(!runchargeflip){
      if(!SameCharge(electrons)) return false;
    }
    float wmassjj= 100000.;
    int indexj1=0;
    int indexj2=0;

    for(unsigned int ij=0; ij < jets.size()-1; ij++){
      for(unsigned int ij2=ij+1; ij2 < jets.size(); ij2++){
	snu::KParticle jjtmp = jets.at(ij) + jets.at(ij2) ;
	if(fabs(jjtmp.M() - 80.4) < wmassjj) {
	  wmassjj = fabs(jjtmp.M() - 80.4);
	  indexj1=ij;
	  indexj2=ij2;
	}
      }
    }

    snu::KParticle jj = jets.at(indexj1) + jets.at(indexj2);
    if(jj.M() > 110.) return false;
    if(jj.M() < 50.) return false;
    if(evmet > 35.) return false;
    
    return true;
  }
  
  if(electrons.size() != 2 ) return false;
  
  if(electrons.at(0).Pt() < cuts.at(0)) return false;
  if(!cut_on_pt){
    if(electrons.at(1).Pt() < 10.) return false;
  }
  else if(electrons.at(1).Pt() < cuts.at(1)) return false;
  
  if(!runchargeflip){
    if(!SameCharge(electrons)) return false;
  }
  if(jets.size() < 2) return false;
  
  
  snu::KParticle ee = electrons.at(0) + electrons.at(1);
  float mass=ee.M();
  
  if(mass < cuts.at(2)) return false;

  if(!cut_on_ee) {
    if(mass > 100000000.) return false;
  }
  else if(mass > cuts.at(3)) return false;
  if(fabs(mass-90.) < 10.) return false;

  if(jets.at(0).Pt() < cuts.at(9)) return false;
   
  float wmassjj= 100000.;
  int indexj1=0;
  int indexj2=0;

  for(unsigned int ij=0; ij < jets.size()-1; ij++){
    for(unsigned int ij2=ij+1; ij2 < jets.size(); ij2++){
      snu::KParticle jjtmp = jets.at(ij) + jets.at(ij2) ;
      if(fabs(jjtmp.M() - 80.4) < wmassjj) {
        wmassjj = fabs(jjtmp.M() - 80.4);
        indexj1=ij;
        indexj2=ij2;
      }
    }
  }
  snu::KParticle jj = jets.at(indexj1) + jets.at(indexj2);
  if(jj.M() > cuts.at(5)) return false;
  if(jj.M() < cuts.at(4)) return false;

  snu::KParticle eejj = electrons.at(0) + electrons.at(1)+ jets.at(indexj1) + jets.at(indexj2) ;

  if(!cut_on_eejj) {
    if(eejj.M()  > 200.) return false;
  }
  else  if(eejj.M()  > cuts.at(7)) return false;
  if(eejj.M()  < cuts.at(6)) return false;

  snu::KParticle e1jj = electrons.at(0)+ jets.at(indexj1) + jets.at(indexj2) ;
  //if(e1jj.M() > cuts.at(13)) return false;
  //if(e1jj.M() < cuts.at(12)) return false;

  snu::KParticle e2jj = electrons.at(1)+ jets.at(indexj1) + jets.at(indexj2) ;
  
  if(cut_on_me2jj){
    if(e2jj.M() > cuts.at(11)) return false;
    if(e2jj.M() < cuts.at(10)) return false;
  }

  if(removeclosejet){
    if(HasCloseBJet( electrons.at(0)) ) return false;
    if(HasCloseBJet( electrons.at(1)) ) return false;
  }

  if(evmet > cuts.at(8)) return false;

  if(nbjet > 0) return false;

  return true;

}


bool HNDiElectronOptimisation::LowMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runchargeflip){
  if(electrons.size() != 2 ) return false;
  if(electrons.at(0).Pt() < 20.) return false;
  if(electrons.at(1).Pt() < 15.) return false;
  if(!runchargeflip){
    if(!SameCharge(electrons)) return false;
  }
  if(jets.size() < 2) return false;
  snu::KParticle ee = electrons.at(0) + electrons.at(1);
  if(ee.M()  < 10.) return false;
  if(fabs(ee.M()-90.)  < 10.) return false;
  
  float wmassjj= 100000.;
  int indexj1=0;
  int indexj2=0;
  for(unsigned int ij=0; ij < jets.size()-1; ij++){
    for(unsigned int ij2=ij+1; ij2 < jets.size(); ij2++){
      snu::KParticle jjtmp = jets.at(ij) + jets.at(ij2);
      if(fabs(jjtmp.M() - 80.4) < wmassjj) {
        wmassjj = fabs(jjtmp.M() - 80.4);
        indexj1=ij;
        indexj2=ij2;
      }
    }
  }
  snu::KParticle jj = jets.at(indexj1) + jets.at(indexj2);
  if(jj.M() > 120.) return false;

  snu::KParticle eejj = electrons.at(0) + electrons.at(1)+ jets.at(indexj1) + jets.at(indexj2) ;
  if(eejj.M() > 200.)return false;

  //if(eventbase->GetEvent().PFMET() > 30.) return false;
  
  int nbjet=0;
  for(unsigned int ij=0; ij <jets.size(); ij++){
    if(IsBTagged(jets.at(ij), snu::KJet::CSVv2, snu::KJet::Medium, GetPeriod())) nbjet++;
  }
  //if(nbjet > 0) return false;

  return true;
  
}

bool HNDiElectronOptimisation::MidMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runchargeflip){
  
  if(electrons.size() != 2 ) return false;
  if(electrons.at(0).Pt() < 20.) return false;
  if(electrons.at(1).Pt() < 15.) return false;
  if(!runchargeflip){
    if(!SameCharge(electrons)) return false;
  }
  if(jets.size() < 2) return false;
  snu::KParticle ee = electrons.at(0) + electrons.at(1);
  if(ee.M()  < 15.) return false;
  if(fabs(ee.M()-90.)  < 10.) return false;

  if(jets.at(0).Pt() < 30.) return false;

  float wmassjj= 100000.;
  int indexj1=0;
  int indexj2=0;
  for(unsigned int ij=0; ij < jets.size()-1; ij++){
    for(unsigned int ij2=ij+1; ij2 < jets.size(); ij2++){
      snu::KParticle jjtmp = jets.at(ij) + jets.at(ij2) ;
      if(fabs(jjtmp.M() - 80.4) < wmassjj) {
        wmassjj = fabs(jjtmp.M() - 80.4);
	indexj1=ij;
        indexj2=ij2;
      }
    }
  }
  snu::KParticle jj = jets.at(indexj1) + jets.at(indexj2);
  if(jj.M() > 110.) return false;
  if(jj.M() < 50.) return false;

  
  //if(eventbase->GetEvent().PFMET() > 35.) return false;

  int nbjet=0;
  for(unsigned int ij=0; ij <jets.size(); ij++){
    if(IsBTagged(jets.at(ij), snu::KJet::CSVv2, snu::KJet::Medium, GetPeriod())) nbjet++;
  }
  //if(nbjet > 0) return false;

  return true;
}



	
	

bool HNDiElectronOptimisation::HighMassCheckSignalRegion(  std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets, bool runchargeflip){
  if(electrons.size() != 2 ) return false;
  if(electrons.at(0).Pt() < 20.) return false;
  if(electrons.at(1).Pt() < 15.) return false;
  if(!runchargeflip){
    if(!SameCharge(electrons)) return false;
  }
  if(jets.size() < 2) return false;
  snu::KParticle ee = electrons.at(0) + electrons.at(1);
  if(ee.M()  < 15.) return false;
  if(fabs(ee.M()-90.)  < 10.) return false;

  float wmassjj= 100000.;
  int indexj1=0;
  int indexj2=0;
  for(unsigned int ij=0; ij < jets.size()-1; ij++){
    for(unsigned int ij2=ij+1; ij2 < jets.size(); ij2++){
      snu::KParticle jjtmp = jets.at(ij) + jets.at(ij2) ;
      if(fabs(jjtmp.M() - 80.4) < wmassjj) {
        wmassjj = fabs(jjtmp.M() - 80.4);
        indexj1=ij;
        indexj2=ij2;
      }
    }
  }
  snu::KParticle jj = jets.at(indexj1) + jets.at(indexj2);
  if(jj.M() > 110.) return false;
  if(jj.M() < 50.) return false;

  //  if(eventbase->GetEvent().PFMET() > 35.) return false;
  if(jets.at(0).Pt() < 40.) return false;
  
  int nbjet=0;
  for(unsigned int ij=0; ij <jets.size(); ij++){
    if(IsBTagged(jets.at(ij), snu::KJet::CSVv2, snu::KJet::Medium, GetPeriod())) nbjet++;
  }
  //  if(nbjet > 0) return false;

  return true;
}





void HNDiElectronOptimisation::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger << DEBUG << "END OF CYCLE: isdata=" << isData <<  LQLogger::endmsg;

}


void HNDiElectronOptimisation::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  

  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //  DeclareVariable(out_muons, "Signal_Muons");

  
  return;
  
}

HNDiElectronOptimisation::~HNDiElectronOptimisation() {
  
  Message("In HNDiElectronOptimisation Destructor" , INFO);
  
}



void HNDiElectronOptimisation::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void HNDiElectronOptimisation::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNDiElectronOptimisationCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HNDiElectronOptimisation::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



