// $id: ExampleAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNEMu Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "HNEMu.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HNEMu);

/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
HNEMu::HNEMu() :  AnalyzerCore(),  out_electrons(0) {

  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNEMu");

  Message("In HNEMu constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  
  k_met=0., k_emumass=0., k_emujjmass=0., k_l1jjmass=0., k_l2jjmass=0., k_njet=0;
  k_nbjet_m=-1;
  k_emujjmass_lowmass=0., k_l1jjmass_lowmass=0., k_l2jjmass_lowmass=0.;
  k_weight=0.;
  k_l1pt=0., k_l2pt=0., k_j1pt=0., k_jjmass=0., k_jjmass_lowmass=0.;
  k_l1eta=0., k_l2eta=0.;
  
  //MakeCleverHistograms(sighist,"SSemu_0BJet");

}


void HNEMu::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //

   return;
}


void HNEMu::ExecuteEvents()throw( LQError ){
    


  vector<int> run_num;
  vector<int> ev_num;

  run_num.push_back(191226);
  run_num.push_back(193336);


  ev_num.push_back(1261781261);
  ev_num.push_back(80512298);

  bool check_event=false;
  for(int iev =0; iev < ev_num.size() ; iev++){
    if(eventbase->GetEvent().RunNumber() == run_num.at(iev)){
      if(eventbase->GetEvent().EventNumber() == ev_num.at(iev))check_event = true;
    }
  }
  //if(!check_event) return;

  // m_logger << INFO << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;



  k_met=0., k_emumass=0., k_emujjmass=0., k_l1jjmass=0., k_l2jjmass=0., k_njet=0;
  k_emujjmass_lowmass=0., k_l1jjmass_lowmass=0., k_l2jjmass_lowmass=0.;

  k_nbjet_m=-1;
  k_weight=0.;
  k_l1pt=0., k_l2pt=0., k_j1pt=0., k_jjmass=0.,k_jjmass_lowmass=0.;
  k_l1eta=0., k_l2eta=0.;
  
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
  
  std::vector<TString> triggerslist2;
  triggerslist2.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  triggerslist2.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  
  //if(PassTrigger(triggerslist2, prescale)) tmp_weight *= TriggerScaleFactorEMu();
  
  FillEventCutFlow("NoCut","", 1.);

  if(!PassMETFilter())  throw LQError( "Fails basic cuts",  LQError::SkipEvent );
  
  FillEventCutFlow("EventCut", "",weight);

  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) throw LQError( "Fails basic cuts",  LQError::SkipEvent );
  FillEventCutFlow("VertexCut","", weight);
    
  
  
}// End of exeucte event loop


float HNEMu::WeightCFEvent(std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons, bool runchargeflip){

  if(electrons.size()!=1) return 0.;
  if(muons.size()!=1) return 0.;
  
  if(runchargeflip) {
    if(electrons.at(0).Charge() != muons.at(0).Charge()) {
      float cf1=  m_datadriven_bkg->CFRate(electrons.at(0));  /// FIX AFTER SETTLED

      return  cf1;
    }// OS requirement
    else return 0.;
  }// cf requirement
  else {
    if(electrons.at(0).Charge() != muons.at(0).Charge()) return 0.;
  }
  
  return 1.;
  
}


void HNEMu::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger << INFO << "Number of os mc events = " << m_os_Z_nw  << LQLogger::endmsg; 
  m_logger << INFO << "Number of os mc events (weighted) = " << m_os_Z  << LQLogger::endmsg; 
  m_logger << INFO << "Number of ss mc events = " << m_ss_Z_nw  << LQLogger::endmsg; 
  m_logger << INFO << "Number of ss mc events (weighted)= " << m_ss_Z  << LQLogger::endmsg; 
}


void HNEMu::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  

  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  
  DeclareVariable(k_met, "met", "MyTree");
  DeclareVariable(k_emumass, "emu_mass", "MyTree");
  DeclareVariable(k_emujjmass, "emujj_mass", "MyTree");
  DeclareVariable(k_l1jjmass, "l1jj_mass", "MyTree");
  DeclareVariable(k_l2jjmass, "l2jj_mass", "MyTree");
  DeclareVariable(k_emujjmass_lowmass, "emujj_mass_lowmass", "MyTree");
  DeclareVariable(k_l1jjmass_lowmass, "l1jj_mass_lowmass", "MyTree");
  DeclareVariable(k_l2jjmass_lowmass, "l2jj_mass_lowmass", "MyTree");
  DeclareVariable(k_njet, "njet", "MyTree");
  DeclareVariable(k_nbjet_m , "nbjet_m",  "MyTree");
  DeclareVariable(k_jjmass, "jj_mass", "MyTree");
  DeclareVariable(k_jjmass_lowmass, "jj_mass_lowmass", "MyTree");
  DeclareVariable(k_l1pt, "l1_pt", "MyTree");
  DeclareVariable(k_l1eta, "l1_eta", "MyTree");
  DeclareVariable(k_l2eta, "l2_eta", "MyTree");
  DeclareVariable(k_l2pt, "l2_pt", "MyTree");

  DeclareVariable(k_j1pt, "jet1_pt", "MyTree");
  DeclareVariable(k_weight, "weight", "MyTree");

  return;
  
}

HNEMu::~HNEMu() {
  
  Message("In HNEMu Destructor" , INFO);
 }
     

void HNEMu::FillEventCutFlow(TString cut, TString label , float weight){

  if(GetHist(label + "_eventcutflow")) {
    GetHist(label + "_eventcutflow")->Fill(cut,weight);

  }
  else{
    AnalyzerCore::MakeHistograms(label + "_eventcutflow",19,0.,19.);

    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(2,"NoCut_w");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(3,"eventcut");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(4,"TriggerCut");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(5,"VertexCut");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(6,"DiEl");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(7,"eedR");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(8,"SSDiEl");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(9,"SS_lepveto");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(10,"DiJet");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(11,"Presel");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(12,"Presel_noZ");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(13,"Presel_nobjet");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(14,"lowmass");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(15,"lowmassCR");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(16,"mediummass");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(17,"mediummassCR");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(18,"highmass");
    GetHist(label + "_eventcutflow")->GetXaxis()->SetBinLabel(19,"highmassCR");
  }
  
}

     
void HNEMu::FillCutFlow(TString cut, float weight){
  
  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow",16,0.,16.);
    
    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"SS_NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"SS_Tight");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"SS_Tight_convveto");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"SS_Tight_d0veto");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"SS_Tight_reliso");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"SS_Medium_chargeconst");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(7,"SS_Tight_chargeconst");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(8,"SS_Tight_noclosejet");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(9,"SS_anal_el");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(10,"Signal_anal");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(11,"Signal_Tightlooseiso_d0");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(12,"Signal_Mediumlooseiso_d0");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(13,"Signal_drcut1");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(14,"Signal_drcut2");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(15,"Signal_anal_dr1");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(16,"Signal_anal_dr2");
  }
}
     
     
void HNEMu::FillIsoCutFlow(TString cut, float weight){
       
  
  if(GetHist("isocutflow")) {
    GetHist("isocutflow")->Fill(cut,weight);
    
  }
  else{
    AnalyzerCore::MakeHistograms("isocutflow",36,0.,36.);
    
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(1,"iso_d0_03_iso3_60");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(2,"iso_d0_03_iso3_50");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(3,"iso_d0_03_iso3_40");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(4,"iso_d0_03_iso3_30");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(5,"iso_d0_03_iso3_20");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(6,"iso_d0_03_iso3_10");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(7,"iso_d0_03_iso3_09");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(8,"iso_d0_03_iso3_075");
    
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(10,"iso_d0_02_iso3_60");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(11,"iso_d0_02_iso3_50");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(12,"iso_d0_02_iso3_40");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(13,"iso_d0_02_iso3_30");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(14,"iso_d0_02_iso3_20");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(15,"iso_d0_02_iso3_10");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(16,"iso_d0_02_iso3_09");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(17,"iso_d0_02_iso3_075");

    GetHist("isocutflow")->GetXaxis()->SetBinLabel(19,"iso_d0_01_iso3_60");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(20,"iso_d0_01_iso3_50");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(21,"iso_d0_01_iso3_40");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(22,"iso_d0_01_iso3_30");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(23,"iso_d0_01_iso3_20");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(24,"iso_d0_01_iso3_10");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(25,"iso_d0_01_iso3_09");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(26,"iso_d0_01_iso3_075");

    GetHist("isocutflow")->GetXaxis()->SetBinLabel(28,"iso_d0_005_iso3_60");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(29,"iso_d0_005_iso3_50");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(30,"iso_d0_005_iso3_40");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(31,"iso_d0_005_iso3_30");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(32,"iso_d0_005_iso3_20");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(33,"iso_d0_005_iso3_10");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(34,"iso_d0_005_iso3_09");
    GetHist("isocutflow")->GetXaxis()->SetBinLabel(35,"iso_d0_005_iso3_075");

    
    
  }
}


void HNEMu::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void HNEMu::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNEMuCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HNEMu::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}




//  LocalWords:  masscuts jetResdown
