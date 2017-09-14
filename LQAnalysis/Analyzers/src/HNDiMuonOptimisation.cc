// $Id: ExampleAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNDiMuonOptimisation Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *d
 ***************************************************************************/

/// Local includes
#include "HNDiMuonOptimisation.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HNDiMuonOptimisation);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
HNDiMuonOptimisation::HNDiMuonOptimisation() :  AnalyzerCore(),  out_electrons(0) {


  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNDiMuonOptimisation");

  Message("In HNDiMuonOptimisation constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void HNDiMuonOptimisation::InitialiseAnalysis() throw( LQError ) {
  
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


void HNDiMuonOptimisation::ExecuteEvents()throw( LQError ){
  
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;


  if(!PassMETFilter()) return;

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

  std::vector<TString> triggerlist_single;
  triggerlist_single.push_back("HLT_IsoMu24_v");
  triggerlist_single.push_back("HLT_IsoTkMu24_v");
  
  float trigger_ps_weight= WeightByTrigger(triggerlist, TargetLumi);
  
  bool _dimu =   isData ?  (k_channel.Contains("DoubleMuon")) : true ;
  bool _singlemu =   isData ?  (k_channel.Contains("SingleMuon")) : true ;
  
  if(_dimu){
    if(PassTrigger(triggerlist[0]) || PassTrigger(triggerlist[1])) OptimiseID(true, true,  true, weight);
  }
  if(_singlemu){
    if (!(PassTrigger(triggerlist[0]) || PassTrigger(triggerlist[1]))){
      if(PassTrigger(triggerlist_single[0]) || PassTrigger(triggerlist_single[1])) OptimiseID(true, true,  false, weight);
    }
  }

  std::vector<snu::KElectron> electronVetoColl=GetElectrons("ELECTRON_HN_VETO");
  std::vector<snu::KMuon> muonVetoColl=GetMuons("MUON_HN_VETO");

  if ((electronVetoColl.size() + muonVetoColl.size()) >2) return;
  
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

void HNDiMuonOptimisation::OptimiseID(bool isss, bool dilep, bool isdileptrig, float w){
  
  
  if(isdileptrig)FillHist("NCIDREFSNU_single", 1,w, 0.,2., 2);
  

  std::vector<snu::KElectron> elColl = GetElectrons("ELECTRON_HN_VETO");
  std::vector<snu::KMuon> muonColl = GetMuons("MUON_HN_VETO");
  std::vector<snu::KJet> jets           = GetJets("JET_HN");
  std::vector<snu::KJet> alljets           = GetJets("JET_NOLEPTONVETO");  

  if(muonColl.size() > 2) return;

  if (!isData) {
    w*= eventbase->GetEvent().PileUpWeight();
  }
  
  for (int nLoops = 0 ; nLoops < 2; nLoops++){
    std::vector<snu::KMuon> snu_loose ;
    TString muonTag;
    if(nLoops ==0) {
      snu_loose  = GetMuons("MUON_HN_Tight_FAKELOOSEST",false);
      muonTag="Tight";
    }
    if(nLoops ==1) {
      snu_loose = GetMuons("MUON_HN_Medium_FAKELOOSEST",false);
      muonTag="Medium";
    }


    std::vector<snu::KElectron> el = GetElectrons("ELECTRON_HN_VETO");
    
    vector<float> vcut_dxy_b;
    vector<TString> vcut_dxy_b_s;

    for(unsigned int dxy_b=0;dxy_b < 1; dxy_b++){
      float cut_dxy_b =  float(dxy_b)*0.005 + 0.005;
      vcut_dxy_b.push_back(cut_dxy_b);
      stringstream ss;
      ss <<cut_dxy_b;
      vcut_dxy_b_s.push_back(TString(ss.str()));
    }
    vector<float> vcut_dxysig_b;
    vcut_dxysig_b.push_back(3.);
    vcut_dxysig_b.push_back(4.);
    vector<TString> vcut_dxysig_b_s;
    vcut_dxysig_b_s.push_back("3");
    vcut_dxysig_b_s.push_back("4");

    vector<float> vcut_dz_b;
    vector<TString> vcut_dz_b_s;

    for(unsigned int dz_b=0;dz_b < 4; dz_b++){
      float cut_dz_b =  float(dz_b)*0.02 + 0.02;
      vcut_dz_b.push_back(cut_dz_b);
      stringstream ss;
      ss <<cut_dz_b;
      vcut_dz_b_s.push_back(TString(ss.str()));
    }

    vector<float> vcut_iso_b;
    vector<TString> vcut_iso_b_s;
    for(unsigned int iso_b=0;iso_b < 5; iso_b++){
      float cut_iso_b = float(iso_b)*0.01 + 0.05;
      vcut_iso_b.push_back(cut_iso_b);
      stringstream ss;
      ss <<cut_iso_b;
      vcut_iso_b_s.push_back(TString(ss.str()));
    }

    int ncut=0;
    //// Loop over cuts and fill loose and tight el and get fake rates for ID                                                                    
    for(unsigned int dxy_b=0; dxy_b < vcut_dxy_b.size(); dxy_b++){
      for(unsigned int dxysig_b=0; dxysig_b < vcut_dxysig_b.size(); dxysig_b++){
        for(unsigned int dz_b=0; dz_b < vcut_dz_b.size(); dz_b++){
          for(unsigned int iso_b=0; iso_b < vcut_iso_b.size(); iso_b++){
	    for(unsigned int dz_e=0; dz_e < vcut_dz_b.size(); dz_e++){
	      for(unsigned int iso_e=0; iso_e < vcut_iso_b.size(); iso_e++, ncut++){
		
		
		std::vector<snu::KMuon> loose_mu;
		std::vector<snu::KMuon> tight_mu;
		bool tight1=false;
		bool tight2=false;
		
		TString fake_tag1("");
		TString fake_tag2("");
		int imupass=0;
		for(unsigned int imu=0; imu<snu_loose.size(); imu++){
		  float reliso = snu_loose[imu].RelIso04();
		  
		  if(fabs(snu_loose[imu].Eta()) < 1.5) {
		    if(fabs(snu_loose[imu].dZ()) > vcut_dz_b[dz_b]) continue;
		    
		    if(imupass==0) fake_tag1=vcut_iso_b_s[iso_b]+"_"+vcut_dxy_b_s[dxy_b]+"_"+vcut_dxysig_b_s[dxysig_b]+"_"+vcut_dz_b_s[dz_b];
		    if(imupass==1) fake_tag2=vcut_iso_b_s[iso_b]+"_"+vcut_dxy_b_s[dxy_b]+"_"+vcut_dxysig_b_s[dxysig_b]+"_"+vcut_dz_b_s[dz_b];
		    //cout << imu << " B " << vcut_iso_b_s[iso_b]+"_"+vcut_dxy_b_s[dxy_b]+"_"+vcut_dxysig_b_s[dxysig_b]+"_"+vcut_dz_b_s[dz_b] << endl;
		    imupass++;
		    loose_mu.push_back(snu_loose[imu]);
		    if(fabs(snu_loose[imu].dXYSig()) > vcut_dxysig_b[dxysig_b]) continue;
		    if(snu_loose[imu].GlobalChi2()  > 10.) continue;
		    if(fabs(snu_loose[imu].dXY()) > vcut_dxy_b[dxy_b]) continue;
		    if(reliso > vcut_iso_b[iso_b]) continue;
		    if(imupass==0) tight1=true;
		    if(imupass==1) tight2=true;
		    tight_mu.push_back(snu_loose[imu]);
		  }
		  else{
		    if(fabs(snu_loose[imu].dZ()) > vcut_dz_b[dz_e]) continue;
		    
		    if(imupass==0) fake_tag1=vcut_iso_b_s[iso_e]+"_"+vcut_dxy_b_s[dxy_b]+"_"+vcut_dxysig_b_s[dxysig_b]+"_"+vcut_dz_b_s[dz_e];
		    if(imupass==1) fake_tag2=vcut_iso_b_s[iso_e]+"_"+vcut_dxy_b_s[dxy_b]+"_"+vcut_dxysig_b_s[dxysig_b]+"_"+vcut_dz_b_s[dz_e];
		    
		    loose_mu.push_back(snu_loose[imu]);
		    imupass++;
		    if(fabs(snu_loose[imu].dXYSig()) > vcut_dxysig_b[dxysig_b]) continue;
		    if(snu_loose[imu].GlobalChi2()  > 10.) continue;
		    if(fabs(snu_loose[imu].dXY()) > vcut_dxy_b[dxy_b]) continue;
		    if(reliso > vcut_iso_b[iso_e]) continue;
		    if(imupass==0) tight1=true;
		    if(imupass==1) tight2=true;
		    
		    tight_mu.push_back(snu_loose[imu]);
		    
		  }
		}
		
		bool cb_1(false);
		bool cb_2(false);
		
		
		if(loose_mu.size()==2){
		  for(unsigned int ij =0 ; ij < alljets.size() ; ij++){
		    if(loose_mu.at(0).DeltaR(alljets.at(ij)) < 0.5) {
		      if(alljets.at(ij).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) cb_1=true;
		    }
		    if(loose_mu.at(1).DeltaR(alljets.at(ij)) < 0.5) {
		      if(alljets.at(ij).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) cb_1=true;
		    }
		  }
		}
		    
		
		bool passtrigcuts_double=false;
		bool passtrigcuts_single=false;
		
		if(loose_mu.size() ==2){
		  if(isdileptrig){
		    if(loose_mu[0].Pt() > 20){
		      if(loose_mu[1].Pt() > 10){
			    passtrigcuts_double=true;
		      }
		    }
		  }
		  else{
		    if(loose_mu[0].Pt() > 25){
		      if(loose_mu[1].Pt() > 5){
			passtrigcuts_single=true;
		      }
		    }
		  }
		}
		
		if(passtrigcuts_single){
		  if(!k_running_nonprompt){
		    
		    if(ncut==0)FillHist("IDREFSNU_single" ,1,w, 0.,2., 2);
		    if(CheckSignalRegion(true,tight_mu,el, jets,alljets,"Low", w))FillHist(("LowIDREFSNU"+muonTag+"_single") ,ncut,w, 0.,1000., 1000);
		    if(CheckSignalRegion(true,tight_mu,el, jets,alljets,"", w))FillHist(("MediumIDREFSNU"+muonTag+"_single") ,ncut,w, 0.,1000., 1000);
		    if(CheckSignalRegion(true,tight_mu,el, jets,alljets,"High", w))FillHist(("HighIDREFSNU"+muonTag+"_single"),ncut,w, 0.,1000., 1000);
		    
		  }
		  else{
		    if(loose_mu.size() ==2){
		      float ev_weight=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, loose_mu, tight1, tight2,muonTag+fake_tag1, muonTag+fake_tag2,cb_1,cb_2 ,"ptcorr_eta", vcut_iso_b[iso_b],vcut_iso_b[iso_e],false, true);
		      if(CheckSignalRegion(true,loose_mu,el, jets,alljets,"Low", ev_weight))FillHist(("LowIDREFSNU"+muonTag+"_single"),ncut,ev_weight, 0.,1000., 1000);
		      if(CheckSignalRegion(true,loose_mu,el, jets,alljets,"", ev_weight))FillHist(("MediumIDREFSNU"+muonTag+"_single") ,ncut,ev_weight, 0.,1000., 1000);
		      if(CheckSignalRegion(true,loose_mu,el, jets,alljets,"High", ev_weight))FillHist(("HighIDREFSNU"+muonTag+"_single"),ncut,ev_weight, 0.,1000., 1000);

		    }
		  }
		}
		if(passtrigcuts_double){
		  if(!k_running_nonprompt){
		    
		    if(ncut==0)FillHist("IDREFSNU_double" ,1,w, 0.,2., 2);
		    if(CheckSignalRegion(true,tight_mu,el, jets,alljets,"Low", w))FillHist(("LowIDREFSNU"+muonTag+"_double") ,ncut,w, 0.,1000., 1000);
		    if(CheckSignalRegion(true,tight_mu,el, jets,alljets,"", w))FillHist(("MediumIDREFSNU"+muonTag+"_double") ,ncut,w, 0.,1000., 1000);
		    if(CheckSignalRegion(true,tight_mu,el, jets,alljets,"High", w))FillHist(("HighIDREFSNU"+muonTag+"_double"),ncut,w, 0.,1000., 1000);
		  }
		  else{
		    if(loose_mu.size()==2){
		      float ev_weight=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, loose_mu, tight1, tight2,muonTag+fake_tag1, muonTag+fake_tag2 , cb_1, cb_2,"ptcorr_eta", vcut_iso_b[iso_b],vcut_iso_b[iso_e],false, false);
		      
		      if(CheckSignalRegion(true,loose_mu,el, jets,alljets,"Low", ev_weight))FillHist(("LowIDREFSNU"+muonTag+"_double"), ncut,ev_weight, 0.,1000., 1000);
		      if(CheckSignalRegion(true,loose_mu,el, jets,alljets,"", ev_weight))FillHist(("MediumIDREFSNU"+muonTag+"_double") ,ncut,ev_weight, 0.,1000., 1000);
		      if(CheckSignalRegion(true,loose_mu,el, jets,alljets,"High", ev_weight))FillHist(("HighIDREFSNU"+muonTag+"_double"),ncut,ev_weight, 0.,1000., 1000);
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
  if(1){
    std::vector<snu::KJet> jets           = GetJets("JET_HN");
    std::vector<snu::KJet> alljets           = GetJets("JET_NOLEPTONVETO");
    
    
    std::vector<snu::KMuon> pogmedium_loose  = GetMuons("MUON_POG_FAKEMEDIUM",false);
    std::vector<snu::KMuon> pogtight_loose  = GetMuons("MUON_POG_FAKETIGHT",false);
    std::vector<snu::KElectron> el = GetElectrons("ELECTRON_HN_VETO");
    std::vector<snu::KMuon> pogmedium_tight  = GetMuons("MUON_POG_MEDIUM",false);
    std::vector<snu::KMuon> pogtight_tight  = GetMuons("MUON_POG_TIGHT",false);
    
    std::vector<snu::KMuon> gent_loose =  GetMuons("MUON_HNGENT_LOOSE",false);
    std::vector<snu::KMuon> gent_tight= GetMuons("MUON_HNGENT_TIGHT",false);
    
    std::vector<snu::KMuon>  hn_loose =  GetMuons("MUON_HN_LOOSE",false);
    std::vector<snu::KMuon>  hn_tight =  GetMuons("MUON_HN_TIGHT",false);
    
    



    bool passtrigcuts_double_pogm=false;
    bool passtrigcuts_single_pogm=false;
    
    if(pogmedium_loose.size() ==2){
      if(isdileptrig){
	if(pogmedium_loose[0].Pt() > 20){
	  if(pogmedium_loose[1].Pt() > 10){
	    passtrigcuts_double_pogm=true;
	  }
	}
      }
      else{
	if(pogmedium_loose[0].Pt() > 25){
	  if(pogmedium_loose[1].Pt() > 5){
	    passtrigcuts_single_pogm=true;
	  }
	}
      }
    }
    

    bool passtrigcuts_double_pogt=false;
    bool passtrigcuts_single_pogt=false;

    if(pogtight_loose.size() ==2){
      if(isdileptrig){
        if(pogtight_loose[0].Pt() > 20){
          if(pogtight_loose[1].Pt() > 10){
            passtrigcuts_double_pogt=true;
          }
        }
      }
      else{
        if(pogtight_loose[0].Pt() > 25){
          if(pogtight_loose[1].Pt() > 5){
            passtrigcuts_single_pogt=true;
          }
        }
      }
    }
    

    bool passtrigcuts_double_gent=false;
    bool passtrigcuts_single_gent=false;

    if(gent_loose.size() ==2){
      if(isdileptrig){
        if(gent_loose[0].Pt() > 20){
          if(gent_loose[1].Pt() > 10){
            passtrigcuts_double_gent=true;
          }
        }
      }
      else{
        if(gent_loose[0].Pt() > 25){
          if(gent_loose[1].Pt() > 5){
            passtrigcuts_single_gent=true;
          }
        }
      }
    }

    bool passtrigcuts_double_hn=false;
    bool passtrigcuts_single_hn=false;

    if(hn_loose.size() ==2){
      if(isdileptrig){
	if(hn_loose[0].Pt() > 20){
          if(hn_loose[1].Pt() > 10){
            passtrigcuts_double_hn=true;
          }
        }
      }
      else{
	if(hn_loose[0].Pt() > 25){
          if(hn_loose[1].Pt() > 5){
            passtrigcuts_single_hn=true;
          }
        }
      }
    }




    if(!isdileptrig)    FillHist("POGIDREFSNU_single", 1,w, 0.,2., 2);
    else FillHist("POGIDREFSNU_double", 1,w, 0.,2., 2);
    

    if(!k_running_nonprompt){
      
      if(passtrigcuts_single_pogm){
	if(CheckSignalRegion(true,pogmedium_tight,el, jets,alljets,"Low", w))FillHist(("LowIDREFsingle") ,0,w, 0.,4., 4);
	if(CheckSignalRegion(true,pogmedium_tight,el, jets,alljets,"", w))FillHist(("MediumIDREFsingle") ,0,w, 0.,4., 4);
	if(CheckSignalRegion(true,pogmedium_tight,el, jets,alljets,"High", w))FillHist(("HighIDREFsingle") ,0,w, 0.,4., 4);
      }
      if(passtrigcuts_single_pogt){
        if(SameCharge(pogtight_tight)&&elColl.size()==0 &&  muonColl.size() == 2 && jets.size() > 1) FillCLHist(sighist_mm, "passtrigcuts_single_pogt", eventbase->GetEvent(),  pogtight_tight, elColl,jets, alljets, w);

	if(CheckSignalRegion(true,pogtight_tight,el, jets,alljets,"Low", w))FillHist(("LowIDREFsingle") ,1,w, 0.,4., 4);
	if(CheckSignalRegion(true,pogtight_tight,el, jets,alljets,"", w))FillHist(("MediumIDREFsingle") ,1,w, 0.,4., 4);
	if(CheckSignalRegion(true,pogtight_tight,el, jets,alljets,"High", w))FillHist(("HighIDREFsingle") ,1,w, 0.,4., 4);
      }
      if(passtrigcuts_single_gent){
	
	if(CheckSignalRegion(true,gent_tight,el, jets,alljets,"Low", w))FillHist(("LowIDREFsingle") ,2,w, 0.,4., 4);
	if(CheckSignalRegion(true,gent_tight,el, jets,alljets,"", w))FillHist(("MediumIDREFsingle") ,2,w, 0.,4., 4);
	if(CheckSignalRegion(true,gent_tight,el, jets,alljets,"High", w))FillHist(("HighIDREFsingle") ,2,w, 0.,4., 4);
      }
      if(passtrigcuts_single_hn){

        if(CheckSignalRegion(true,hn_tight,el, jets,alljets,"Low", w))FillHist(("LowIDREFsingle") ,3,w, 0.,4., 4);
	if(CheckSignalRegion(true,hn_tight,el, jets,alljets,"", w))FillHist(("MediumIDREFsingle") ,3,w, 0.,4., 4);
        if(CheckSignalRegion(true,hn_tight,el, jets,alljets,"High", w))FillHist(("HighIDREFsingle") ,3,w, 0.,4., 4);
      }

      if(passtrigcuts_double_pogm){ 
	if(CheckSignalRegion(true,pogmedium_tight,el, jets,alljets,"Low", w))FillHist(("LowIDREFdouble") ,0,w, 0.,4., 4);
	if(CheckSignalRegion(true,pogmedium_tight,el, jets,alljets,"", w))FillHist(("MediumIDREFdouble") ,0,w, 0.,4., 4);
	if(CheckSignalRegion(true,pogmedium_tight,el, jets,alljets,"High", w))FillHist(("HighIDREFdouble") ,0,w, 0.,4., 4);
      }
      if(passtrigcuts_double_pogt){ 
	if(SameCharge(pogtight_tight)&&elColl.size()==0 && muonColl.size() == 2 && jets.size() > 1) FillCLHist(sighist_mm, "passtrigcuts_double_pogt", eventbase->GetEvent(),  pogtight_tight, elColl,jets, alljets,w);
	
	if(CheckSignalRegion(true,pogtight_tight,el, jets,alljets,"Low", w)){
	  FillCLHist(sighist_mm, "passtrigcuts_double_pogt_lowmass", eventbase->GetEvent(), pogtight_tight, elColl,jets, alljets,w);
	  snu::KParticle Z = pogtight_tight[0]+pogtight_tight[1];
          bool pass1=false;
          bool pass2=false;
          bool pass3=false;
          float ct_lt(0.);
          if(fabs(pogtight_tight[0].Eta()) < 1.5) ct_lt+= pogtight_tight[0].Pt();
          if(fabs(pogtight_tight[1].Eta()) < 1.5) ct_lt+= pogtight_tight[1].Pt();
	  
          if((Z.Pt() > 5. && Z.Pt() < 25.)) pass1=true;
          if(pogtight_tight[0].DeltaR(pogtight_tight[1]) > 2.) pass2=true;
          if(ct_lt < 25) pass3=true;
          if(pass1&&pass2&&pass3)           FillCLHist(sighist_mm, "passtrigcuts_double_pogt_pa3", eventbase->GetEvent(),  pogtight_tight, elColl,jets, alljets, w);
          if(pass1||pass2||pass3)           FillCLHist(sighist_mm, "passtrigcuts_double_pogt_po3", eventbase->GetEvent(),  pogtight_tight, elColl,jets, alljets, w);
	}
										       

	if(CheckSignalRegion(true,pogtight_tight,el, jets,alljets,"Low", w))FillHist(("LowIDREFdouble") ,1,w, 0.,4., 4);
	if(CheckSignalRegion(true,pogtight_tight,el, jets,alljets,"", w))FillHist(("MediumIDREFdouble") ,1,w, 0.,4., 4);
	if(CheckSignalRegion(true,pogtight_tight,el, jets,alljets,"High", w))FillHist(("HighIDREFdouble") ,1,w, 0.,4., 4);
      }
      if(passtrigcuts_double_gent){ 
	if(CheckSignalRegion(true,gent_tight,el, jets,alljets,"Low", w))FillHist(("LowIDREFdouble") ,2,w, 0.,4., 4);
	if(CheckSignalRegion(true,gent_tight,el, jets,alljets,"", w))FillHist(("MediumIDREFdouble") ,2,w, 0.,4., 4);
	if(CheckSignalRegion(true,gent_tight,el, jets,alljets,"High", w))FillHist(("HighIDREFdouble") ,2,w, 0.,4., 4);
      }
      if(passtrigcuts_double_hn){
	if(CheckSignalRegion(true,hn_tight,el, jets,alljets,"Low", w))FillHist(("LowIDREFdouble") ,3,w, 0.,4., 4);
        if(CheckSignalRegion(true,hn_tight,el, jets,alljets,"", w))FillHist(("MediumIDREFdouble") ,3,w, 0.,4., 4);
	if(CheckSignalRegion(true,hn_tight,el, jets,alljets,"High", w))FillHist(("HighIDREFdouble") ,3,w, 0.,4., 4);
	if(CheckSignalRegion(true,hn_tight,el, jets,alljets,"Low", w))FillCLHist(sighist_mm, "passtrigcuts_double_hn_lowmass", eventbase->GetEvent(), hn_tight, elColl,jets, alljets,w);
	
      }
    }
    else{
 
      bool cb_1(false);
      bool cb_2(false);
      
      if(passtrigcuts_single_pogm){
	float ev_weight_med=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, pogmedium_loose, PassID(pogmedium_loose[0],"MUON_POG_MEDIUM"),  PassID(pogmedium_loose[1],"MUON_POG_MEDIUM"), "pogmedium","pogmedium", cb_1, cb_2,"ptcorr_eta",0.25,0.25,false,true);
	
	if(CheckSignalRegion(true,pogmedium_loose,el, jets,alljets,"Low", w))FillHist(("LowIDREFsingle"),0 ,ev_weight_med, 0.,4., 4);
	if(CheckSignalRegion(true,pogmedium_loose,el, jets,alljets,"", w))FillHist(("MediumIDREFsingle"),0 ,ev_weight_med, 0.,4., 4);
	if(CheckSignalRegion(true,pogmedium_loose,el, jets,alljets,"High", w))FillHist(("HighIDREFsingle") ,0,ev_weight_med, 0.,4., 4);
	
      }
      if(passtrigcuts_single_pogt){
	
	float ev_weight_tight=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, pogtight_loose, PassID(pogtight_loose[0],"MUON_POG_TIGHT"),  PassID(pogtight_loose[1],"MUON_POG_TIGHT"), "pogtight","pogtight",cb_1, cb_2, "ptcorr_eta",0.15,0.15,false,true);
	if(CheckSignalRegion(true,pogtight_loose,el, jets,alljets,"Low", w))FillHist(("LowIDREFsingle"),1 ,ev_weight_tight, 0.,4., 4);
	if(CheckSignalRegion(true,pogtight_loose,el, jets,alljets,"", w))FillHist(("MediumIDREFsingle"),1 ,ev_weight_tight, 0.,4., 4);
	if(CheckSignalRegion(true,pogtight_loose,el, jets,alljets,"High", w))FillHist(("HighIDREFsingle") ,1,ev_weight_tight, 0.,4., 4);
	if(SameCharge(pogtight_loose)&&elColl.size() == 0 && muonColl.size() == 2&& jets.size() > 1){
	  FillCLHist(sighist_mm, "pastrigcuts_single_pogt", eventbase->GetEvent(),  pogtight_loose, elColl,jets, alljets, ev_weight_tight);
	}
      }
      if(passtrigcuts_single_gent){
	
	float ev_weight_gent=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, gent_loose, PassID(gent_loose[0],"MUON_HNGENT_TIGHT"),  PassID(gent_loose[1],"MUON_HNGENT_TIGHT"), "gent","gent",cb_1, cb_2, "ptcorr_eta",0.1,0.1,false,true);  
	
	if(CheckSignalRegion(true,gent_loose,el, jets,alljets,"Low", w))FillHist(("LowIDREFsingle") ,2,ev_weight_gent, 0.,4., 4);
	if(CheckSignalRegion(true,gent_loose,el, jets,alljets,"", w))FillHist(("MediumIDREFsingle") ,2,ev_weight_gent, 0.,4., 4);
	if(CheckSignalRegion(true,gent_loose,el, jets,alljets,"High", w))FillHist(("HighIDREFsingle") ,2,ev_weight_gent, 0.,4., 4);
      }
      if(passtrigcuts_single_hn){

        float ev_weight_hn=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, hn_loose, PassID(hn_loose[0],"MUON_HN_TIGHT"),  PassID(hn_loose[1],"MUON_HN_TIGHT"), "Tight0.09_0.005_3_0.04","Tight0.09_0.005_3_0.04",cb_1, cb_2, "ptcorr_eta",0.09,0.09,false,true);

        if(CheckSignalRegion(true,hn_loose,el, jets,alljets,"Low", w))FillHist(("LowIDREFsingle") ,3,ev_weight_hn, 0.,4., 4);
        if(CheckSignalRegion(true,hn_loose,el, jets,alljets,"", w))FillHist(("MediumIDREFsingle") ,3,ev_weight_hn, 0.,4., 4);
        if(CheckSignalRegion(true,hn_loose,el, jets,alljets,"High", w))FillHist(("HighIDREFsingle") ,3,ev_weight_hn, 0.,4., 4);
      }
      
      
      if(passtrigcuts_double_pogm){
	float ev_weight_med=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, pogmedium_loose, PassID(pogmedium_loose[0],"MUON_POG_MEDIUM"),  PassID(pogmedium_loose[1],"MUON_POG_MEDIUM"), "pogmedium","pogmedium", cb_1, cb_2,"ptcorr_eta",0.25,0.25,false,false);
	if(CheckSignalRegion(true,pogmedium_loose,el, jets,alljets,"Low", w))FillHist(("LowIDREFdouble"),0 ,ev_weight_med, 0.,4., 4);
	if(CheckSignalRegion(true,pogmedium_loose,el, jets,alljets,"", w))FillHist(("MediumIDREFdouble"),0 ,ev_weight_med, 0.,4., 4);
	if(CheckSignalRegion(true,pogmedium_loose,el, jets,alljets,"High", w))FillHist(("HighIDREFdouble") ,0,ev_weight_med, 0.,4., 4);

        if(CheckSignalRegion(true,pogmedium_loose,el, jets,alljets,"Low", w))FillHist(("LowIDREFdoublenoW"),0 ,1., 0.,4., 4);
	
      }

      if(passtrigcuts_double_pogt){
	

	
	float ev_weight_tight=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, pogtight_loose, PassID(pogtight_loose[0],"MUON_POG_TIGHT"),  PassID(pogtight_loose[1],"MUON_POG_TIGHT"), "pogtight","pogtight",cb_1, cb_2, "ptcorr_eta",0.15,0.15,false,false);

	if(CheckSignalRegion(true,pogtight_loose,el, jets,alljets,"Low", w))FillHist(("LowIDREFdouble"),1 ,ev_weight_tight, 0.,4., 4);
	if(CheckSignalRegion(true,pogtight_loose,el, jets,alljets,"", w))FillHist(("MediumIDREFdouble"),1 ,ev_weight_tight, 0.,4., 4);
	if(CheckSignalRegion(true,pogtight_loose,el, jets,alljets,"High", w))FillHist(("HighIDREFdouble") ,1,ev_weight_tight, 0.,4., 4);
	
	if(CheckSignalRegion(true,pogtight_loose,el, jets,alljets,"Low", w))FillHist(("LowIDREFdoublenoW"),1 ,1., 0.,4., 4);


	if(SameCharge(pogtight_loose) && elColl.size() == 0 &&muonColl.size() == 2&& jets.size() > 1)FillCLHist(sighist_mm, "passtrigcuts_double_pogt", eventbase->GetEvent(),  pogtight_loose, elColl,jets, alljets, ev_weight_tight);
	if(CheckSignalRegion(true,pogtight_loose,el, jets,alljets,"Low", w)){
	  FillCLHist(sighist_mm, "passtrigcuts_double_pogt_lowmass", eventbase->GetEvent(), pogtight_loose, elColl,jets, alljets, ev_weight_tight);
	  snu::KParticle Z = pogtight_loose[0]+pogtight_loose[1];
          bool pass1=false;
          bool pass2=false;
          bool pass3=false;
          float ct_lt(0.);
          if(fabs(pogtight_loose[0].Eta()) < 1.5) ct_lt+= pogtight_loose[0].Pt();
          if(fabs(pogtight_loose[1].Eta()) < 1.5) ct_lt+= pogtight_loose[1].Pt();

          if((Z.Pt() > 5. && Z.Pt() < 25.)) pass1=true;
          if(pogtight_loose[0].DeltaR(pogtight_loose[1]) > 2.) pass2=true;
          if(ct_lt < 25) pass3=true;
          if(pass1&&pass2&&pass3)           FillCLHist(sighist_mm, "passtrigcuts_double_pogt_pa3", eventbase->GetEvent(),  pogtight_loose, elColl,jets, alljets, ev_weight_tight);
          if(pass1||pass2||pass3)           FillCLHist(sighist_mm, "passtrigcuts_double_pogt_po3", eventbase->GetEvent(),  pogtight_loose, elColl,jets, alljets, ev_weight_tight);


	}
										       

      }

      if(passtrigcuts_double_gent){
	
	float ev_weight_gent=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, gent_loose, PassID(gent_loose[0],"MUON_HNGENT_TIGHT"),  PassID(gent_loose[1],"MUON_HNGENT_TIGHT"), "gent","gent",cb_1, cb_2, "ptcorr_eta",0.1,0.1,false,false);
	
	if(CheckSignalRegion(true,gent_loose,el, jets,alljets,"Low", w))FillHist(("LowIDREFdouble") ,2,ev_weight_gent, 0.,4., 4);
	if(CheckSignalRegion(true,gent_loose,el, jets,alljets,"", w))FillHist(("MediumIDREFdouble") ,2,ev_weight_gent, 0.,4., 4);
	if(CheckSignalRegion(true,gent_loose,el, jets,alljets,"High", w))FillHist(("HighIDREFdouble") ,2,ev_weight_gent, 0.,4., 4);
	
	if(CheckSignalRegion(true,gent_loose,el, jets,alljets,"Low", w))FillHist(("LowIDREFdoublenoW") ,2,1., 0.,4., 4);

      }
      if(passtrigcuts_double_hn){
	
	if(hn_loose.size() ==2){
	  float ev_weight_hn=m_datadriven_bkg->Get_DataDrivenWeight_MM(false, hn_loose, PassID(hn_loose[0],"MUON_HN_TIGHT"),  PassID(hn_loose[1],"MUON_HN_TIGHT"), "Tight0.07_0.005_3_0.04","Tight0.07_0.005_3_0.04",cb_1, cb_2, "ptcorr_eta",0.07,0.07,false,false);
	  
	  if(CheckSignalRegion(true,hn_loose,el, jets,alljets,"Low", w))FillHist(("LowIDREFdouble") ,3,ev_weight_hn, 0.,4., 4);
	  if(CheckSignalRegion(true,hn_loose,el, jets,alljets,"", w))FillHist(("MediumIDREFdouble") ,3,ev_weight_hn, 0.,4., 4);
	  if(CheckSignalRegion(true,hn_loose,el, jets,alljets,"High", w))FillHist(("HighIDREFdouble") ,3,ev_weight_hn, 0.,4., 4);

	  if(CheckSignalRegion(true,hn_loose,el, jets,alljets,"Low", w))FillHist(("LowIDREFdoublenoW") ,3,1., 0.,4., 4);
	  if(CheckSignalRegion(true,hn_loose,el, jets,alljets,"Low", w))FillCLHist(sighist_mm, "passtrigcuts_double_hn_lowmass", eventbase->GetEvent(), hn_loose, elColl,jets, alljets,ev_weight_hn);
	}
      }
    }
  }
      

  
}
	

void HNDiMuonOptimisation::GetTriggEfficiency(){
  //ListTriggersAvailable();                                                                                                                                                                                                                                                                                                                                                      
  vector<int> pt1;

  vector<int>pt2;


  std::vector<TString> lists_triggers;
  //lists_triggers.push_back("HLT_Ele25_eta2p1_WPTight_Gsf_v");

  std::vector<TString> lists_triggerORs;
  //lists_triggerORs.push_back("HLT_Ele15_IsoVVVL_PFHT600_v");

  for(unsigned int i=0; i < lists_triggers.size(); i++){
    FillTriggerEfficiency(lists_triggers.at(i), weight, "denominator_nojet", lists_triggers );
  }

}

void HNDiMuonOptimisation::FillTriggerEfficiency(TString cut, float weight, TString label, std::vector<TString> list){

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


bool HNDiMuonOptimisation::CheckSignalRegion( bool isss,  std::vector<snu::KMuon> muons, std::vector<snu::KElectron> el,  std::vector<snu::KJet> jets,  std::vector<snu::KJet> alljets, TString name, float w){

  bool debug=false;
  
  
  if(el.size() > 0) return false ;
  if(muons.size() != 2 ) {if(debug)cout << "Fail el size" << endl; return false ;}
  // Set by trigger
  if(muons.at(0).Pt() < 20.) {if(debug)cout << "Fail pt1 " << endl; return false;}
  if(muons.at(1).Pt() < 10.)  {if(debug)cout << "Fail pt2  " << endl; return false;}
  
  if(isss&&!SameCharge(muons)) {if(debug)cout << "Fail ss " << endl; return false;}

  if(!isss&&SameCharge(muons)) {if(debug)cout << "Fail os  " << endl; return false;}
  if(jets.size() < 2) {if(debug)cout << "Fail jets " << endl; return false;}

  snu::KParticle mm = muons.at(0) + muons.at(1);
  if(mm.M()  < 10.) {if(debug)cout << "Fail mee  " << endl; return false;}
  
  if(name.Contains("Low")){
    if(mm.M()  > 70.) {if(debug)cout << "Fail mee  " << endl; return false;}
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
  
  float dPhi = fabs(TVector2::Phi_mpi_pi(jets[m].Phi() - jets[n].Phi()));
  float contramass=2*jets[m].Pt()*jets[n].Pt()*(1+cos(dPhi));
  contramass=sqrt(contramass);
  
  if(name.Contains("Low")){
    if(contramass > 100) return false;
    if((jets[0]+jets[1]).M() > 200.) return false;
    if((jets[0] + muons[0] + muons[1]).M() > 250.)  return false;
  }
  else{
    if(jj.M() > 120.) {if(debug)cout << "Fail mjj  " << endl; return false;}
    if(jj.M() < 50.) {if(debug)cout << "Fail mjj  " << endl; return false;}
    
  }
  
  if(name.Contains("Low")){
    if(muons.at(0).DeltaR(muons.at(1)) > 3.5) {if(debug)cout << "Fail dr  " << endl; return false;}
  }
  
  float ST = muons[0].Pt() + muons[1].Pt();
  float looseST = muons[0].Pt() + muons[1].Pt();

  std::vector<snu::KJet> loosejets=GetJets("JET_NOCUT");
  
  
  float HT=0.;
  for(unsigned int ij=0; ij <alljets.size(); ij++){
    ST+= alljets[ij].Pt();

  }
  for(unsigned int ij=0; ij <jets.size(); ij++){

    HT+= jets[ij].Pt();
    
  }

  float LT = muons[0].Pt() + muons[1].Pt();
  if(name.Contains("Low")){
    
    if((pow(eventbase->GetEvent().PFMET(),2.)/ ST) > 15.) {if(debug)cout << "Fail met  " << endl; return false;}
    if(eventbase->GetEvent().PFMET() > 80) return false;
  }
  else if(name.Contains("High")){
    if((pow(eventbase->GetEvent().PFMET(),2.)/ ST) > 8.) {if(debug)cout << "Fail met  " << endl; return false;}
  }
  else{
    if((pow(eventbase->GetEvent().PFMET(),2.)/ ST) > 12.5) {if(debug)cout << "Fail met  " << endl; return false;}
    if(eventbase->GetEvent().PFMET() > 80) return false;

  }
  

  float dphi_1 = fabs(TVector2::Phi_mpi_pi(muons.at(0).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet)));
  float MT_1 = sqrt(2.* muons.at(0).Et()*eventbase->GetEvent().PFMET() * (1 - cos( dphi_1)));
  
  float dphi_2 = fabs(TVector2::Phi_mpi_pi(muons.at(0).Phi()- eventbase->GetEvent().METPhi(snu::KEvent::pfmet)));
  float MT_2 = sqrt(2.* muons.at(1).Et()*eventbase->GetEvent().PFMET() * (1 - cos( dphi_2)));
  
  if(name.Contains("Low")){
    if(MT_1 > 100.) return false;
    if(MT_2 > 100.) return false;
    if(ST > 450.)  return false;
    snu::KParticle lljj = jets.at(m) + jets.at(n) + muons.at(0) + muons.at(1);
    //if(lljj.M() > 300.) return false;
    
  }
  else   if(name.Contains("High")){
    snu::KParticle lljj = jets.at(m) + jets.at(n) + muons.at(0) + muons.at(1);

    if(ST < 400.)  return false;
    if(lljj.M() < 200.) return false;
    if(muons.at(0).Pt() < 50.) return false;
    if(muons.at(1).Pt() < 25.) return false;
  }
  else{
    
    if(LT/HT < 0.2)  return false;
    if(LT < 35)  return false;
    if(HT < 35)  return false;

    float dRmj=-999.;
    float ratiomj=0.;
    for(unsigned int im=0; im <muons.size(); im++){
      for(unsigned int ij2=0; ij2 <jets.size(); ij2++){
        if(muons.at(im).DeltaR(jets.at(ij2)) >  0.5){
	  if(muons.at(im).DeltaR(jets.at(ij2)) > dRmj){
	    dRmj=muons.at(im).DeltaR(jets.at(ij2)) ;
	    if(im==0) ratiomj = jets.at(ij2).Pt() / muons[im].Pt();
	  }
	}
      }
    }
    if(ratiomj > 3.) return false;
    if(dRmj > 4.5) return false;
    
    if(muons[0].DeltaR(muons[1]) > 4.) return false;
    
  }
  
  if(name.Contains("High")){
    if((muons.at(0).Pt() + muons.at(1).Pt() ) < 50.) return false;
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
  for(unsigned int ij=0; ij <alljets.size(); ij++){
    if( alljets[ij].IsBTagged(  snu::KJet::CSVv2, snu::KJet::Medium))  nbjet++;
  }
  if(nbjet > 0) {if(debug)cout << "Fail nbjet  " << endl; return false;}
  
  if(debug)cout << "PASSES ID" << endl;
  return true;
}

bool HNDiMuonOptimisation::CheckSignalRegionNN( bool isss,  std::vector<snu::KMuon> muons, std::vector<snu::KJet> jets, TString name, float w){

  if(muons.size() != 2 ) return false ;
  // Set by trigger                                                                                                                                                             
  if(muons.at(0).Pt() < 25.) return false;
  if(muons.at(1).Pt() < 15.) return false;
  if(isss&&!SameCharge(muons)) return false;

  if(!isss&&SameCharge(muons)) return false;
  if(jets.size() < 4) return false;
  snu::KParticle ee = muons.at(0) + muons.at(1);
  if(ee.M()  < 10.) return false;

  if((ee.M() > 80.) && (ee.M() < 100.) ) return false;

  int nbjet=0;
  for(unsigned int ij=0; ij <jets.size(); ij++){
    if( jets[ij].IsBTagged(  snu::KJet::CSVv2, snu::KJet::Medium))  nbjet++;

  }
  if(nbjet > 0) return false;

  return true;

}




void HNDiMuonOptimisation::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger << DEBUG << "END OF CYCLE: isdata=" << isData <<  LQLogger::endmsg;

}


void HNDiMuonOptimisation::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  

  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Muons", "LQTree");
  //  DeclareVariable(out_muons, "Signal_Muons");

  
  return;
  
}

HNDiMuonOptimisation::~HNDiMuonOptimisation() {
  
  Message("In HNDiMuonOptimisation Destructor" , INFO);
  
}



void HNDiMuonOptimisation::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void HNDiMuonOptimisation::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNDiMuonOptimisationCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HNDiMuonOptimisation::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



