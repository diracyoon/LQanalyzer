// $Id: ExampleAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNSignalEfficiencies Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *d
 ***************************************************************************/

/// Local includes
#include "HNSignalEfficiencies.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HNSignalEfficiencies);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
HNSignalEfficiencies::HNSignalEfficiencies() :  AnalyzerCore(),  out_electrons(0) {


  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNSignalEfficiencies");

  Message("In HNSignalEfficiencies constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void HNSignalEfficiencies::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //

   Message("Making clever hists for Z ->ll test code", INFO);

   //// Initialise Plotting class functions
   /// MakeCleverHistograms ( type, "label")  type can be muhist/elhist/jethist/sighist


   vector<TString> labels;
   labels.push_back("hn"); 
   //labels.push_back("medium");   // pog numbers
   labels.push_back("tight");    // pog numbers 
   
   for(unsigned int il = 0; il < labels.size(); il++){
     TString label = labels.at(il);
     continue;
     MakeCleverHistograms(sighist_ee,label + "_SSee_1jet");
     MakeCleverHistograms(sighist_ee,label + "_SSee_DiJet");
     MakeCleverHistograms(sighist_ee,label + "_SSPreselection");
     MakeCleverHistograms(sighist_ee,label + "_OSee_1jet");
     MakeCleverHistograms(sighist_ee,label + "_OSee_DiJet");
     MakeCleverHistograms(sighist_ee,label + "_OSPreselection");
     MakeCleverHistograms(sighist_ee,label + "_LowMass");
     MakeCleverHistograms(sighist_ee,label + "_MediumMass");
     MakeCleverHistograms(sighist_ee,label + "_HighMass");
     MakeCleverHistograms(sighist_ee,label + "_LowMassCR");
     MakeCleverHistograms(sighist_ee,label + "_MediumMassCR");
     MakeCleverHistograms(sighist_ee,label + "_HighMassCR");
   }
   
   

   /// Validation signalplots
   /*
   MakeCleverHistograms(sighist_ee, "SIGNALVALIDATION_EE" );
   MakeCleverHistograms(sighist_mm, "SIGNALVALIDATION_MM" );
   */

   return;
}


void HNSignalEfficiencies::ExecuteEvents()throw( LQError ){
  
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;


  counter("All",1.);
  //if(FailHNDataSetCheck()) return;

  if(!isData){
    weight*= MCweight;
    weight*=WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi);
  }
  
  vector<TString> ElTriggers = GetHNDiLepElTriggers();

  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; 

  counter("Vertex",1.);


  TString eltag = "ELECTRON_NOCUT"; 
  if(isData)eltag = "ELECTRON_HN_HIGHDXY_FAKELOOSE";

  TString mutag = "MUON_NOCUT";
  if(isData)mutag="MUON_HN_TRI_HIGHDXY_LOOSE";

  std::vector<snu::KMuon> muons_nc=GetMuons(mutag,false);
  //if(muons_nc.size() > 0) cout << muons_nc[0].Pt() << " " << muons_nc[0].Eta() << endl;                                                                            
  /*
  if(muons_nc.size()  == 3){
    cout << "TRI " << endl;
    TruthPrintOut();
    for(unsigned int iel =0 ; iel < muons_nc.size(); iel++){
      if(muons_nc[iel].Pt() < 5.) continue;
      cout << "El pt = " << muons_nc[iel].Pt() << " eta ="  << muons_nc[iel].Eta() << " " << muons_nc[iel].Phi() << endl;
      cout << "type      :  Mother pdgid    :   Isprompt(SK)    :  prompt flag  "<<  endl;
      cout << muons_nc[iel].GetType() << " " << muons_nc[iel].MotherPdgId() << " " << muons_nc[iel].MCIsPrompt() << " " << muons_nc[iel].IsPromptFlag() << endl;
      cout << "Ext conv     : sk conv" << endl;
    }
  }
  else  if(SameCharge(muons_nc)){
    TruthPrintOut();
    for(unsigned int iel =0 ; iel < muons_nc.size(); iel++){
      if(muons_nc[iel].Pt() < 5.) continue;
      cout << "El pt = " << muons_nc[iel].Pt() << " eta ="  << muons_nc[iel].Eta() << " " << muons_nc[iel].Phi() << endl;
      cout << "type      :  Mother pdgid    :   Isprompt(SK)    :  prompt flag  "<<  endl;
      cout << muons_nc[iel].GetType() << " " << muons_nc[iel].MotherPdgId() << " " << muons_nc[iel].MCIsPrompt() << " " << muons_nc[iel].IsPromptFlag() << endl;
      cout << "Ext conv     : sk conv" << endl;
    }
  }
  else return;

  */
  std::vector<snu::KElectron> electrons_nc=GetElectrons(false,false,eltag);
  //if(electrons_nc.size() > 0) cout << electrons_nc[0].Pt() << " " << electrons_nc[0].Eta() << endl;

  if(electrons_nc.size() ==3){
    cout << "TRIEL" << endl;
    TruthPrintOut();
    
    for(unsigned int iel =0 ; iel < electrons_nc.size(); iel++){
      if(electrons_nc[iel].Pt() < 10.) continue;
      cout << "El pt = " << electrons_nc[iel].Pt() << " eta ="  << electrons_nc[iel].Eta() << " " << electrons_nc[iel].Phi() << endl;
      cout << "type      :  Mother pdgid    :   Isprompt(SK)    :  prompt flag  "<<  endl;
      cout << electrons_nc[iel].GetType() << " " << electrons_nc[iel].MotherPdgId() << " " << electrons_nc[iel].MCIsPrompt() << " " << electrons_nc[iel].IsPromptFlag() << endl;
      cout << "Ext conv     : sk conv" << endl;
      cout << electrons_nc[iel].MCIsExternalConversion() << "  " << electrons_nc[iel].MCIsFromConversion() << endl;
    }

  }
  if(SameCharge(electrons_nc)){
    TruthPrintOut();
    for(unsigned int iel =0 ; iel < electrons_nc.size(); iel++){
      if(electrons_nc[iel].Pt() < 10.) continue;
      cout << "El pt = " << electrons_nc[iel].Pt() << " eta ="  << electrons_nc[iel].Eta() << " " << electrons_nc[iel].Phi() << endl;
      cout << "type      :  Mother pdgid    :   Isprompt(SK)    :  prompt flag  "<<  endl;
      cout << electrons_nc[iel].GetType() << " " << electrons_nc[iel].MotherPdgId() << " " << electrons_nc[iel].MCIsPrompt() << " " << electrons_nc[iel].IsPromptFlag() << endl;
      cout << "Ext conv     : sk conv" << endl;
      cout << electrons_nc[iel].MCIsExternalConversion() << "  " << electrons_nc[iel].MCIsFromConversion() << endl;
    }
  }
  else return;

  

  std::vector<snu::KElectron> electrons_tm_nc;

  std::vector<snu::KElectron> looseveto = GetElectrons(false,false,"ELECTRON_HN_FAKELOOSEST");
  
  if(looseveto.size() < 2) return;
  if(!PassTrigger("HLT_Ele27_WPTight_Gsf_v")) return;
  if(!SameCharge(looseveto)) return;
  for(unsigned int iel=0; iel<looseveto.size(); iel++){
    
    if(looseveto.at(iel).TriggerMatched("HLT_Ele27_WPTight_Gsf_v")) {
      if(fabs(looseveto[iel].dxy()) < 0.01){
	if(looseveto[iel].Pt() > 30){
	  if(looseveto[iel].PFRelIso(0.3) > 0.05){
	    cout << "None tight el " << endl;
	  }
	  else             cout << "Ttight el " << endl;

	}
      }
    }
    float el_pt_corr = looseveto.at(iel).Pt()*(1+max(0.,(looseveto.at(iel).PFRelIso(0.3)-0.05))) ; /// will need changing for systematics                                         

    float reliso = looseveto[iel].PFRelIso(0.3);
    bool pass_trigger_emulation=true;
    if(looseveto[iel].Pt() < 15.){
      if(!looseveto[iel].PassHLTID()) pass_trigger_emulation=false;
    }
    else{
      if(!looseveto[iel].IsTrigMVAValid()) pass_trigger_emulation=false;
    }
    if(!pass_trigger_emulation) continue;
    if(fabs(looseveto[iel].dz()) > 0.2) continue;
    if(fabs(looseveto[iel].dxy()) > 0.1) continue;

    if(fabs(looseveto[iel].SCEta())<0.8 ){
      if(looseveto[iel].MVA() < -0.02) continue;
    }
    else  if(fabs(looseveto[iel].SCEta())<1.479 ){
      if(looseveto[iel].MVA() < -0.52) continue;
    }
    else {
      if(looseveto[iel].MVA() < -0.52) continue;
    }
    electrons_tm_nc.push_back(looseveto[iel]);
  }



  std::vector<snu::KMuon> muons_tm_nc=GetMuons("MUON_NOCUT",false);

  if((muons_tm_nc.size() + electrons_tm_nc.size()) > 1)  counter("SkimReq",1.);
  if((muons_tm_nc.size() + electrons_tm_nc.size())== 1){
    if(muons_tm_nc.size()  > 0 ){
      if(muons_tm_nc[0].Pt() > 5.)  counter("SkimReq",1.);
    }
  }

  for(unsigned int i=0; i < electrons_nc.size(); i++){
    //  cout << "el " << i+1 << electrons_nc[i].Pt() << " " << electrons_nc[i].Eta() << " " << electrons_nc[i].Phi() << " " << electrons_nc[i].MCMatched() << endl;
  }
  
  if(electrons_tm_nc.size() >=2)counter("NCDiEl",1.);

  
  std::vector<snu::KJet> jets_nc = GetJets("JET_HN");  

  if(!isData)GetTriggEfficiency();

  std::vector<snu::KElectron>  hn_electrons=GetElectrons(false,false,"ELECTRON_HN_TIGHT");
  std::vector<snu::KElectron>  pog_electrons=GetElectrons(false,false,"ELECTRON_POG_TIGHT");

  std::vector<snu::KMuon>  hn_muons=GetMuons("MUON_HN_TIGHT",false);
  std::vector<snu::KMuon>  pog_muons=GetMuons("MUON_POG_TIGHT",false);


  std::vector<snu::KFatJet> fatjetcoll = GetFatJets("FATJET_HN");
  std::vector<snu::KJet> jets = GetJets("JET_NOLEPTONVETO");
  std::vector<snu::KJet> hn04jets = GetJets("JET_HN");
  std::vector<snu::KJet> hnjets =  GetJetsWFT("JET_HN","FATJET_HN");
  


  FillCLHist(sighist_ee, "NoCut", eventbase->GetEvent(), muons_tm_nc,electrons_tm_nc,jets_nc, jets,fatjetcoll,weight);
  if(SameCharge(electrons_tm_nc))   FillCLHist(sighist_ee, "SSNoCut", eventbase->GetEvent(), muons_tm_nc,electrons_tm_nc,jets_nc, jets,fatjetcoll,weight);
  if(SameCharge(electrons_nc)) {
    if(electrons_nc[1].Pt() > 10){
      if(electrons_nc[0].MCMatched()&&!electrons_nc[1].MCMatched()){
	FillHist("El_typePF", electrons_nc[1].GetType(),weight, 0., 40., 40);
	FillHist("El_motherpdgidPF", fabs(electrons_nc[1].MotherPdgId()),weight, 0., 450., 450);
	FillHist("El_motherpdgidPF_typePF",electrons_nc[1].GetType(), fabs(electrons_nc[1].MotherPdgId()),weight, 0., 40., 40, 0., 450., 450);

      }
      else if(!electrons_nc[0].MCMatched()&&electrons_nc[1].MCMatched()){
	FillHist("El_typeFP", electrons_nc[0].GetType(),weight, 0., 40., 40);
	FillHist("El_motherpdgidFP", fabs(electrons_nc[0].MotherPdgId()),weight, 0., 450., 450);
        FillHist("El_motherpdgidFP_typeFP",electrons_nc[0].GetType(), fabs(electrons_nc[0].MotherPdgId()),weight, 0., 40., 40, 0., 450., 450);

      }
      else if(!electrons_nc[0].MCMatched()&&!electrons_nc[1].MCMatched()){
	FillHist("El_typeFF", electrons_nc[0].GetType(),weight, 0., 40., 40);
	FillHist("El_typeFF", electrons_nc[1].GetType(),weight, 0., 40., 40);
	FillHist("El_motherpdgidFF", fabs(electrons_nc[0].MotherPdgId()),weight, 0., 450., 450);
        FillHist("El_motherpdgidFF_typeFF",electrons_nc[0].GetType(), fabs(electrons_nc[0].MotherPdgId()),weight, 0., 40., 40, 0., 450., 450);
	FillHist("El_motherpdgidFF", fabs(electrons_nc[1].MotherPdgId()),weight, 0., 450., 450);
        FillHist("El_motherpdgidFF_typeFF",electrons_nc[1].GetType(), fabs(electrons_nc[1].MotherPdgId()),weight, 0., 40., 40, 0., 450., 450);

      }
    }
  }

  if(hn_electrons.size() == 2){
    counter("DiEl",1.);
    if(SameCharge(hn_electrons))     counter("SSDiEl",1.);
    else counter("OSDiEl",1.);
  }
  else if(hn_electrons.size() == 3)  counter("TriEl",1.);
  else if(hn_electrons.size() == 1)  counter("SingleEl",1.);
  
  if(hn_muons.size() == 2){
    counter("DiMu",1.);
    if(SameCharge(hn_muons))     counter("SSDiMu",1.);
    else counter("OSDiMu",1.);
  }
  else if(hn_muons.size() == 3)  counter("TriMu",1.);
  else if(hn_muons.size() == 1)  counter("SingleMu",1.);

  
  int nmatched=0;
    
  for(unsigned int iel =0 ; iel < electrons_nc.size(); iel++){
    if(electrons_nc[iel].Pt() < 20.) continue;
    
      
    if((electrons_nc[iel].Eta()) < 1.5)FillHist("El_energy_b_reco", electrons_nc[iel].Energy(), weight,   0., 1000, 1000); 
    else  FillHist("El_energy_ec_reco", electrons_nc[iel].Energy(), weight,   0., 1000, 1000); 
    if(!isData)FillHist("El_energy_truth", GetTruthMatchedParticle(electrons_nc[iel]).E(), weight,   0., 1000, 1000);

    bool prompt = electrons_nc[iel].MCMatched();
    
    if(isData){
      prompt = PassID(electrons_nc[iel], "ELECTRON_HN_HIGHDXY_TIGHT");
    }
    
    if(prompt){
      nmatched++;
      FillHist("PromptElectron_iso_pt", electrons_nc[iel].Pt(),electrons_nc[iel].PFRelIso(0.3) ,weight,   0., 1000, 100, 0.,2., 200);
      FillHist("PromptElectron_pt", electrons_nc[iel].Pt() ,weight,  0., 1000, 1000);

      float massdrop=MassDrop(electrons_nc[iel], jets);
      float closejetpt= electrons_nc[iel].Pt()/ GetJetsCloseToLeptonPt(electrons_nc[iel], jets);
      float closejetrelpt= GetPtRelLepTJet(electrons_nc[iel], jets);
      
      if(electrons_nc[iel].Pt()> 50){
	FillHist("Prompt_massdrop", massdrop,  weight, 0., 1000, 1000);
	FillHist("Prompt_massdrop_ratio",closejetpt, massdrop,  weight, 0., 1000, 1000,  0., 2., 100);
	FillHist("Prompt_massdrop_rel",closejetrelpt, massdrop,  weight, 0., 1000, 1000,  0., 200, 200);
      }
      FillHist("Prompt_ratio" ,electrons_nc[iel].Pt(), closejetpt,  weight, 0., 1000, 100,  0., 2., 100);
      FillHist("Prompt_relpt" ,electrons_nc[iel].Pt(), closejetrelpt,  weight, 0., 1000, 100,  0., 200, 200); 

      if(electrons_nc[iel].PFRelIso(0.3) > 0.05){
	FillHist("NonIsoPromptElectron_pt", electrons_nc[iel].Pt() ,weight,  0., 1000, 1000);
	FillHist("NonIsoPrompt_massdrop", massdrop,  weight, 0., 1000, 1000);
	FillHist("NonIsoPrompt_massdrop_ratio",closejetpt, massdrop,  weight, 0., 1000, 1000,  0., 2., 100);
	FillHist("NonIsoPrompt_massdrop_rel",closejetrelpt, massdrop,  weight, 0., 1000, 1000,  0., 200, 200);
	FillHist("NonIsoPrompt_ratio_rel",closejetpt, closejetrelpt,  weight, 0., 2., 100,  0., 100, 100);
	
	if(closejetpt > 0.8 && closejetrelpt > 8){
	  FillHist("NonIsoFix1PromptElectron_pt", electrons_nc[iel].Pt() ,weight,  0., 1000, 1000);
	  FillHist("NonIsoFix1Prompt_massdrop", massdrop,  weight, 0., 1000, 1000);
	}
	if(massdrop > 10.){
	  FillHist("NonIsoFix2PromptElectron_pt", electrons_nc[iel].Pt() ,weight,  0., 1000, 1000);
	  FillHist("NonIsoFix2Prompt_ratio" ,electrons_nc[iel].Pt(), closejetpt,  weight, 0., 1000, 100,  0., 2., 100);
	  FillHist("NonIsoFix2Prompt_relpt" ,electrons_nc[iel].Pt(), closejetrelpt,  weight, 0., 1000, 100,  0., 200, 200);
	  FillHist("NonIsoFix2Prompt_ratio_rel",closejetpt, closejetrelpt,  weight, 0., 2., 100,  0., 100, 100);

	}
	FillHist("NonIsoPrompt_ratio" ,electrons_nc[iel].Pt(), closejetpt,  weight, 0., 1000, 100,  0., 2., 100);
	FillHist("NonIsoPrompt_relpt" ,electrons_nc[iel].Pt(), closejetrelpt,  weight, 0., 1000, 100,  0., 200, 200);
      }    
    }
    else{
      FillHist("NPromptElectron_iso_pt", electrons_nc[iel].Pt(),electrons_nc[iel].PFRelIso(0.3) ,weight,  0., 1000, 100, 0.,2., 200);
      FillHist("NPromptElectron_pt", electrons_nc[iel].Pt() ,weight,  0., 1000, 1000);

      float massdrop=MassDrop(electrons_nc[iel], jets);
      float closejetpt= electrons_nc[iel].Pt()/GetJetsCloseToLeptonPt(electrons_nc[iel], jets);
      float closejetrelpt= GetPtRelLepTJet(electrons_nc[iel], jets);
      FillHist("NPrompt_massdrop", massdrop,  weight, 0., 1000, 1000);
      FillHist("NPrompt_massdrop_ratio",closejetpt, massdrop,  weight, 0., 1000, 1000,  0., 2., 100);
      FillHist("NPrompt_massdrop_rel",closejetrelpt, massdrop,  weight, 0., 1000, 1000,  0., 200, 200);
      FillHist("NPrompt_ratio" ,electrons_nc[iel].Pt(), closejetpt,  weight, 0., 1000, 100,  0., 2., 100);
      FillHist("NPrompt_relpt" ,electrons_nc[iel].Pt(), closejetrelpt,  weight, 0., 1000, 100,  0., 100, 100);
    }
  }
  
    return;
  for(unsigned int iel =0 ; iel < muons_nc.size(); iel++){

    bool prompt = muons_nc[iel].MCMatched();

    if(isData){
      prompt = PassID(muons_nc[iel], "MUON_HN_TRI_HIGHDXY_TIGHT");
    }
    
    if(muons_nc[iel].MCMatched()){
      nmatched++;

      FillHist("PromptMuon_iso_pt", muons_nc[iel].Pt(),muons_nc[iel].RelIso04() ,weight,  0., 1000, 100, 0.,2., 200);
      float massdrop=MassDrop(muons_nc[iel], jets);
      float closejetpt=muons_nc[iel].Pt()/GetJetsCloseToLeptonPt(muons_nc[iel], jets);
      float closejetrelpt= GetPtRelLepTJet(muons_nc[iel], jets);
      FillHist("PromptMuon_pt", muons_nc[iel].Pt() ,weight,  0., 1000, 1000);
      FillHist("Prompt_massdrop", massdrop,  weight, 0., 1000, 1000);
      FillHist("Prompt_ratio" ,muons_nc[iel].Pt(), closejetpt,  weight, 0., 1000, 100,  0., 100, 100);
      FillHist("Prompt_relpt" ,muons_nc[iel].Pt(), closejetrelpt,  weight, 0., 1000, 100,  0., 100, 100);
      
      }
    else{
      FillHist("NPromptMuon_pt", muons_nc[iel].Pt() ,weight,  0., 1000, 1000);
      FillHist("NPromptMuon_iso_pt", muons_nc[iel].Pt(),muons_nc[iel].RelIso04() ,weight,  0., 1000, 100, 0.,2., 200);
      float massdrop=MassDrop(muons_nc[iel], jets);
      float closejetpt=muons_nc[iel].Pt()/GetJetsCloseToLeptonPt(muons_nc[iel], jets);
      float closejetrelpt= GetPtRelLepTJet(muons_nc[iel], jets);
      FillHist("NPrompt_massdrop", massdrop,  weight, 0., 1000, 1000);
      FillHist("NPrompt_ratio" ,muons_nc[iel].Pt(), closejetpt,  weight, 0., 1000, 100,  0., 100, 100);
      FillHist("NPrompt_relpt" ,muons_nc[iel].Pt(), closejetrelpt,  weight, 0., 1000, 100,  0., 100, 100);
    }
  }

  if(SameCharge(hn_electrons)) { 
    counter("SSHNee",1);
    FillCLHist(sighist_ee, "SSHNee", eventbase->GetEvent(), hn_muons,hn_electrons,hn04jets, jets,fatjetcoll,weight);
  }
  if(SameCharge(pog_electrons)) {
    counter("SSPOGee",1);
    FillCLHist(sighist_ee, "SSPOGee", eventbase->GetEvent(), pog_muons,pog_electrons,hn04jets, jets,fatjetcoll,weight);
  }
 
  if(SameCharge(hn_muons)) {
    counter("SSHNmm",1);
    FillCLHist(sighist_mm, "SSHNmm", eventbase->GetEvent(), hn_muons,hn_electrons,hn04jets, jets,fatjetcoll,weight);
  }
  if(SameCharge(pog_muons)) {
    counter("SSPOGmm",1);
    FillCLHist(sighist_mm, "SSPOGmm", eventbase->GetEvent(), pog_muons,pog_electrons,hn04jets, jets,fatjetcoll,weight);
  }
 
}


void HNSignalEfficiencies::counter(TString cut, float w){

  map<TString,float>::iterator itmapcounter = mapcounter.find(cut) ;
  if (itmapcounter == mapcounter.end()){
    mapcounter[cut] = w;
  }
  else{
    float sum = itmapcounter->second;
    mapcounter[cut] = sum+w;
  }

  FillCutFlow(cut,w);
}

void HNSignalEfficiencies::GetTriggEfficiency(){
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

void HNSignalEfficiencies::FillTriggerEfficiency(TString cut, float weight, TString label, std::vector<TString> list){

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


void HNSignalEfficiencies::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger << DEBUG << "END OF CYCLE: isdata=" << isData <<  LQLogger::endmsg;


  for ( map<TString,float>::iterator itmapcounter = mapcounter.begin(); itmapcounter != mapcounter.end(); itmapcounter++){
    cout << itmapcounter->first << " has count "<< itmapcounter->second << endl;
  }

}



void HNSignalEfficiencies::BeginCycle() throw( LQError ){
  
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

HNSignalEfficiencies::~HNSignalEfficiencies() {
  
  Message("In HNSignalEfficiencies Destructor" , INFO);
  
}



void HNSignalEfficiencies::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void HNSignalEfficiencies::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNSignalEfficienciesCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HNSignalEfficiencies::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



