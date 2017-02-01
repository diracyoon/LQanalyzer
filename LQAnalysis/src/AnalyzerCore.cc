/***************************************************************************
 * @Project: LQAnalyzer Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes 
#include "AnalyzerCore.h"
#include "EventBase.h"
#include "Reweight.h"

//Plotting                                                      
#include "MuonPlots.h"
#include "ElectronPlots.h"
#include "JetPlots.h"
#include "SignalPlotsEE.h"
#include "SignalPlotsMM.h"
#include "SignalPlotsEM.h"
#include "TriLeptonPlots.h"

//ROOT includes
#include <TFile.h>


AnalyzerCore::AnalyzerCore() : LQCycleBase(), n_cutflowcuts(0), MCweight(-999.),reset_lumi_mask(false),changed_target_lumi(false) {

  bool debug(false);
  
  TH1::SetDefaultSumw2(true);  
  /// clear list of triggers stored in KTrigger
  triggerlist.clear();
  CorrectionMap.clear();
  CorrectionMapGraph.clear();
  cutflow_list.clear();
  // If running on LQNtuples this is not important.
  // If creating an SKTree ntuple this controls what triggers are accessible
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("AnalyzerCore");

  Message("In AnalyzerCore constructor", INFO);
  
  /////////////////////////////////////////////////////////////////////// 
  //////// For HN analysis  /////////////////////////////////////////////  
  //////////////////////////////////////////////////////////////////////  
  //// FakeRate Input file           
  //////////////////////////////////////////////////////////////////////                                                                                                   
  //string analysisdir = getenv("IDFILEDIR");
  string pileupdir = getenv("PILEUPFILEDIR");
  //string triggerdir = getenv("TRIGGERFILEDIR");

  FillCorrectionHists();

  string lqdir = getenv("LQANALYZER_DIR");
  m_fakeobj = new HNCommonLeptonFakes(lqdir+"/LQAnalysis/src/HNCommonLeptonFakes/share/");
  
  if (TString(getenv("CATVERSION")).Contains("v7-6-6")){
    rmcor = new rochcor2015();
  }
  /// Currently only have csvv2 or cMVAv2 btaggers: In HN we use csvv2 
  /// List of taggers
  std::vector<TString> vtaggers;
  vtaggers.push_back("CSVv2");
  //  vtaggers.push_back("cMVAv2");

  if (TString(getenv("CATVERSION")).Contains("v7-6-6")){
    cout << "Setting up 2015 selection " << endl;
    SetupSelectionMuon(lqdir + "/CATConfig/SelectionConfig/muons_2015.sel");
    SetupSelectionMuon(lqdir + "/CATConfig/SelectionConfig/user_muons_2015.sel");
    
    SetupSelectionElectron(lqdir + "/CATConfig/SelectionConfig/electrons_2015.sel");
    SetupSelectionElectron(lqdir + "/CATConfig/SelectionConfig/user_electrons_2015.sel");
    
    SetupSelectionJet(lqdir + "/CATConfig/SelectionConfig/jets_2015.sel");
    SetupSelectionJet(lqdir + "/CATConfig/SelectionConfig/user_jets_2015.sel");
  }
  else{
    cout << "Setting up 2016 selection " << endl;
    SetupSelectionMuon(lqdir + "/CATConfig/SelectionConfig/muons.sel");
    SetupSelectionMuon(lqdir + "/CATConfig/SelectionConfig/user_muons.sel");

    SetupSelectionElectron(lqdir + "/CATConfig/SelectionConfig/electrons.sel");
    SetupSelectionElectron(lqdir + "/CATConfig/SelectionConfig/user_electrons.sel");

    SetupSelectionJet(lqdir + "/CATConfig/SelectionConfig/jets.sel");
    SetupSelectionJet(lqdir + "/CATConfig/SelectionConfig/user_jets.sel");
  }
  
  if(debug){
    for( map<TString,vector<pair<TString,float> > >::iterator it=  selectionIDMapfMuon.begin() ; it !=  selectionIDMapfMuon.end(); it++){
      cout << it->first << endl;
      for (unsigned int i=0 ; i < it->second.size(); i++){
	cout << it->second.at(i).first << " " << it->second.at(i).second << endl;
      }
    }
  }
  Message("SetupSelection DONE", DEBUG); 

  
  // List of working points
  std::vector<TString> v_wps;
  v_wps.push_back("Loose");
  v_wps.push_back("Medium");
  v_wps.push_back("Tight");
  MapBTagSF = SetupBTagger(vtaggers,v_wps);
  cout << "                                                  " << endl;
  cout << "   ########    ###       ###   ###  ###           " << endl;
  cout << "   ########    ####      ###   ###  ###           " << endl;
  cout << "   ###         #####     ###   ###  ###           " << endl;
  cout << "   ###         ### ##    ###   ###  ###           " << endl;
  cout << "   ########    ###  ##   ###   ###  ###           " << endl;
  cout << "   ########    ###   ##  ###   ###  ###           " << endl;
  cout << "        ###    ###    ## ###   ###  ###           " << endl;
  cout << "        ###    ###     #####   ###  ###           " << endl;
  cout << "   ########    ###      ####   ########           " << endl;
  cout << "   ########    ###       ###   ########           " << endl;
  cout << "                                                  " << endl;

  cout << "##################################################" << endl;
  if(1){
    ifstream runlumi((lqdir + "/data/Luminosity/"+getenv("yeartag")+"/lumi_catversion_" + getenv("CATVERSION")+".txt").c_str());
    if(!runlumi) {
      cerr << "Did not find "+lqdir + "/data/Luminosity/"+getenv("yeartag")+"/lumi_catversion_" + getenv("CATVERSION")+".txt'), exiting ..." << endl;
      exit(EXIT_FAILURE);
    }
    string lline;
    int x=1;
    while(getline(runlumi,lline) ){
      std::istringstream is( lline );
      
      string trigname;
      float trig_lumi;
      int run;
      is >> trigname;
      if(trigname=="###" ) continue;
      if(trigname=="END") break;
      if(trigname=="run" ){
	is >> run;
	is >> trig_lumi;
	cout << "Run number "<< run <<" ; Muon trigger (unprescaled) luminosity  " << trig_lumi << ";" << endl;
	
	mapLumi2016[run] = trig_lumi;
	continue;
      }
      if(trigname=="block" ){
	is >> run;
	is >> trig_lumi;
	cout << "mapLumi[" << run <<" ] = " << trig_lumi << ";" << endl;
	
	mapLumiPerBlock2016[run] = trig_lumi;
	ostringstream ss;
	ss << x;
	mapLumiNamePerBlock2016[run]="Lumi"+ss.str();
	x++;
	continue;
      }

    }
    runlumi.close();
  }
  
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;

  cout << "Reading Luminosity File" << endl;

  if(1){
    ifstream triglumi2016((lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+".txt").c_str());
    if(!triglumi2016) {
      cerr << "Did not find "+lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+".txt'), exiting ..." << endl;
      exit(EXIT_FAILURE);
    }
    string sline2016;
    
    cout << "Trigname : Lumi pb-1" << endl;
    while(getline(triglumi2016,sline2016) ){
      std::istringstream is( sline2016 );
      
      string trigname;
      float trig_lumi;
      is >> trigname;
      if(trigname=="###" ) continue;

      
      if(trigname=="END") break;
      if(!TString(trigname).Contains("Lumi")){
	is >> trig_lumi;
	cout << trigname << " " << trig_lumi << endl;
      }
      else{
	string tmp;
	is >> tmp;
	is >> trig_lumi;
	cout << trigname << " " << trig_lumi << endl;
      }
      trigger_lumi_map_cat2016[TString(trigname)] = trig_lumi;
      continue;
    }
    triglumi2016.close();
  }
  
  cout <<  "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  
  if(TString(getenv("CATVERSION")).Contains("v8")){
    if(TString(getenv("PUFILE")).Contains(".root")){
      if(!k_isdata) reweightPU = new Reweight((pileupdir + "/" + getenv("PUFILE")).c_str());       
    }
  }


}



void AnalyzerCore::SetupSelectionJet(std::string path_sel){
  Message("SetupSelectionJet", DEBUG);
  ifstream jetselconfig(path_sel.c_str());
  if(!jetselconfig) {
    cerr << "Did not find "+ path_sel+", exiting ..." << endl;
    exit(EXIT_FAILURE);
  }
  string jetline;
  int ncuts=0;
  vector<TString> cutnames;
  while(getline(jetselconfig,jetline) ){
    vector<pair<TString,TString> > string_jetsel;
    vector<pair<TString,float> > float_jetsel;
    TString idlabel;
    if (TString(jetline).Contains("webpage")) continue;
    if (TString(jetline).Contains("###")) continue;
    if (TString(jetline).Contains("ncut")) {
      std::istringstream is( jetline );
      string tmp;
      int itmp;
      is >> tmp;
      is >> itmp;
      ncuts = 2*(itmp +1);
      continue;
    }
    if (TString(jetline).Contains("ptmin")) {
      std::istringstream is( jetline );
      string tmp;
      for (int x =0; x < ncuts; x++){
        is >> tmp;
        cutnames.push_back(TString(tmp));
      }
    }
    else{
      std::istringstream is( jetline );
      string tmp;
      float tmpf;
      for (int x =0; x < ncuts; x++){
        if ( x%2 ==0) {
          is >> tmp;
          continue;
        }

        if (x > 6 && x < 18){
          is >> tmp;
          string_jetsel.push_back(make_pair(cutnames.at(x),tmp) );

        }
        else if ( x ==1) {is >> idlabel;}
        else {
          is >> tmpf;
          float_jetsel.push_back(make_pair(cutnames.at(x),tmpf));
        }
      }
    }
    selectionIDMapsJet[idlabel] = string_jetsel;
    selectionIDMapfJet[idlabel] = float_jetsel;
  }
}

void AnalyzerCore::SetupSelectionMuon(std::string path_sel){
  Message("SetupSelectionMuon", DEBUG);
  ifstream muonselconfig(path_sel.c_str());
  if(!muonselconfig) {
    cerr << "Did not find "+ path_sel+", exiting ..." << endl;
    exit(EXIT_FAILURE);
  }
  string muonline;
  int ncuts=0;
  vector<TString> cutnames;
  while(getline(muonselconfig,muonline) ){
    vector<pair<TString,TString> > string_muonsel;
    vector<pair<TString,float> > float_muonsel;
    TString idlabel;
    if (TString(muonline).Contains("webpage")) continue;
    if (TString(muonline).Contains("###")) continue;
    if (TString(muonline).Contains("ncut")) {
      std::istringstream is( muonline );
      string tmp;
      int itmp;
      is >> tmp;
      is >> itmp;
      ncuts = 2*(itmp +1);
      continue;
    }
    if (TString(muonline).Contains("ptmin")) {
      std::istringstream is( muonline );
      string tmp;
      for (int x =0; x < ncuts; x++){
	is >> tmp;
	cutnames.push_back(TString(tmp));
      }
    }
    else{
      std::istringstream is( muonline );
      string tmp;
      float tmpf;
      for (int x =0; x < ncuts; x++){
	if ( x%2 ==0) {
	  is >> tmp;
	  continue;
	}

	if (x > 10 && x < 16){
	  is >> tmp;
	  string_muonsel.push_back(make_pair(cutnames.at(x),tmp) );

	}
	else if ( x ==1) {is >> idlabel;}
	else {
	  is >> tmpf;
	  float_muonsel.push_back(make_pair(cutnames.at(x),tmpf));
	}
      }
    }
    selectionIDMapsMuon[idlabel] = string_muonsel;
    selectionIDMapfMuon[idlabel] = float_muonsel;
  }
}


void AnalyzerCore::SetupSelectionElectron(std::string path_sel){

  Message("SetupSelectionElectron", DEBUG);

  //// currently hard coded set of cuts

  ifstream elselconfig(path_sel.c_str());
  if(!elselconfig) {
    cerr << "Did not find "+ path_sel+", exiting ..." << endl;
    exit(EXIT_FAILURE);
  }
  string elline;
  vector<TString> cutnames;
  int ncuts=0;
  while(getline(elselconfig,elline) ){
    vector<pair<TString,TString> > string_elsel;
    vector<pair<TString,float> > float_elsel;
    TString idlabel;
    if (TString(elline).Contains("webpage")) continue;
    if (TString(elline).Contains("###")) continue;
    if (TString(elline).Contains("ncut")) {
      std::istringstream is( elline );
      string tmp;
      int itmp;
      is >> tmp;
      is >> itmp;
      ncuts = 2*(itmp +1);
      continue;
    }

    if (TString(elline).Contains("ptmin")) {
      std::istringstream is( elline );
      string tmp;
      for (int x =0; x < ncuts; x++){
        is >> tmp;
        cutnames.push_back(TString(tmp));
      }
    }
    else{
      std::istringstream is( elline );
      string tmp;
      float tmpf;

      for (int x =0; x < ncuts; x++){
        if ( x%2 ==0) {
          is >> tmp;
          continue;
        }

        if (x > 10 && x < 18){
          is >> tmp;
          string_elsel.push_back(make_pair(cutnames.at(x),tmp) );
        }
	else  if (x > 26 && x < 30){
          is >> tmp;
          string_elsel.push_back(make_pair(cutnames.at(x),tmp) );
        }
        else if ( x ==1) {is >> idlabel;}
        else {
          is >> tmpf;
          float_elsel.push_back(make_pair(cutnames.at(x),tmpf));
        }
      }
    }
    selectionIDMapsElectron[idlabel] = string_elsel;
    selectionIDMapfElectron[idlabel] = float_elsel;
  }
}


std::map<TString,BTagSFUtil*> AnalyzerCore::SetupBTagger(std::vector<TString> taggers, std::vector<TString> wps){
  
  std::map<TString,BTagSFUtil*>  tmpmap;
  for(std::vector<TString>::const_iterator it = taggers.begin(); it != taggers.end(); it++){
    for(std::vector<TString>::const_iterator it2 = wps.begin(); it2 != wps.end(); it2++){
      tmpmap[*it + "_" + *it2 + "_lf"]= new BTagSFUtil("incl", it->Data(), it2->Data());
      tmpmap[*it +  "_" + *it2 + "_hf"]= new BTagSFUtil("mujets", it->Data(), it2->Data());
    }
  }

  return tmpmap;
}


float AnalyzerCore::GetDiLepMass(std::vector<snu::KElectron> electrons){

  if(electrons.size() != 2) return 0.;
  snu::KParticle p = electrons.at(0) + electrons.at(1);
  return p.M();
}

float AnalyzerCore::GetDiLepMass(std::vector<snu::KMuon> muons){

  if(muons.size() != 2) return 0.;
  snu::KParticle p = muons.at(0) + muons.at(1);
  return p.M();
}


bool AnalyzerCore::EtaRegion(TString reg,  std::vector<snu::KElectron> electrons){
  if(electrons.size() != 2) return false;

  if(reg.Contains("EE")){
    if(fabs(electrons.at(0).Eta()) < 1.5) return false;
    if(fabs(electrons.at(1).Eta()) < 1.5) return false;
    return true;
  }
  if(reg.Contains("BB")){
    if(fabs(electrons.at(0).Eta()) > 1.5) return false;
    if(fabs(electrons.at(1).Eta()) > 1.5) return false;
    return true;
  }
  if(reg.Contains("EB")){
    if(fabs(electrons.at(0).Eta()) > 1.5){
      if(fabs(electrons.at(1).Eta()) > 1.5) return false;
    }
    if(fabs(electrons.at(0).Eta()) < 1.5){
      if(fabs(electrons.at(1).Eta()) < 1.5) return false;
    }
    return true;
  }
  return false;
  
}


bool AnalyzerCore::EtaRegion(TString reg,  std::vector<snu::KMuon> muons){
  
  if(muons.size() != 2) return false;
  if(reg.Contains("EE")){
    if(muons.at(0).Eta() < 1.5) return false;
    if(muons.at(1).Eta() < 1.5) return false;
    return true;
  }
  if(reg.Contains("BB")){
    if(muons.at(0).Eta() > 1.5) return false;
    if(muons.at(1).Eta() > 1.5) return false;
    return true;
  }
  if(reg.Contains("EB")){
    if(muons.at(0).Eta() > 1.5){
      if(muons.at(1).Eta() > 1.5) return false;
    }
    if(muons.at(0).Eta() < 1.5){
      if(muons.at(1).Eta() < 1.5) return false;
    }
    return true;
  }
  return false;
}




std::vector<snu::KJet> AnalyzerCore::GetJets(BaseSelection::ID jetid, float ptcut, float etacut){
  return GetJets(GetStringID(jetid), ptcut, etacut);
}

std::vector<snu::KMuon> AnalyzerCore::GetMuons(BaseSelection::ID muonid, float ptcut, float etacut){
  return GetMuons(GetStringID(muonid), ptcut, etacut);
}

std::vector<snu::KMuon> AnalyzerCore::GetMuons(BaseSelection::ID muonid, bool keepfake, float ptcut, float etacut){
  return GetMuons(GetStringID(muonid), keepfake, ptcut, etacut);
}

std::vector<snu::KElectron> AnalyzerCore::GetElectrons(BaseSelection::ID electronid, float ptcut, float etacut){
  return GetElectrons(true, true,GetStringID(electronid), ptcut, etacut);
}

std::vector<snu::KElectron> AnalyzerCore::GetElectrons(bool keepcf, bool keepfake, BaseSelection::ID electronid, float ptcut, float etacut){
  return GetElectrons(keepcf, keepfake,GetStringID(electronid), ptcut, etacut);
}



TString AnalyzerCore::GetStringID(BaseSelection::ID id){
  if(id == BaseSelection::ELECTRON_POG_VETO ) return "ELECTRON_POG_VETO";
  if(id == BaseSelection::ELECTRON_POG_LOOSE) return "ELECTRON_POG_LOOSE";
  if(id == BaseSelection::ELECTRON_POG_MEDIUM) return "ELECTRON_POG_MEDIUM";
  if(id == BaseSelection::ELECTRON_POG_TIGHT           ) return "ELECTRON_POG_TIGHT";
  if(id == BaseSelection::ELECTRON_POG_MVATrig) return "ELECTRON_POG_MVATrig";
  if(id == BaseSelection::ELECTRON_POG_MVANonTrig    ) return "ELECTRON_POG_MVANonTrig";
  if(id == BaseSelection::ELECTRON_ECAL_FIDUCIAL     ) return "ELECTRON_ECAL_FIDUCIAL";
  if(id == BaseSelection::ELECTRON_HN_VETO           ) return "ELECTRON_HN_VETO";
  if(id == BaseSelection::ELECTRON_HN_TIGHT          ) return "ELECTRON_HN_TIGHT";
  if(id == BaseSelection::ELECTRON_HN_FAKELOOSE      ) return "ELECTRON_HN_FAKELOOSE";
  if(id == BaseSelection::ELECTRON_HN_FAKELOOSE_NOD0 ) return "ELECTRON_HN_FAKELOOSE_NOD0";
  if(id == BaseSelection::ELECTRON_TOP_VETO          ) return "ELECTRON_TOP_VETO";
  if(id == BaseSelection::ELECTRON_TOP_LOOSE         ) return "ELECTRON_TOP_LOOSE";
  if(id == BaseSelection::ELECTRON_TOP_TIGHT         ) return "ELECTRON_TOP_TIGHT";
  if(id == BaseSelection::ELECTRON_PTETA             ) return "ELECTRON_PTETA";
  if(id == BaseSelection::ELECTRON_NOCUT             ) return "ELECTRON_NOCUT";
  if(id == BaseSelection::MUON_POG_LOOSE             ) return "MUON_POG_LOOSE";
  if(id == BaseSelection::MUON_POG_MEDIUM            ) return "MUON_POG_MEDIUM";
  if(id == BaseSelection::MUON_POG_TIGHT             ) return "MUON_POG_TIGHT";
  if(id == BaseSelection::MUON_HN_VETO               ) return "MUON_HN_VETO";
 if(id == BaseSelection::MUON_HN_FAKELOOSE          ) return "MUON_HN_FAKELOOSE";
 if(id == BaseSelection::MUON_HN_TIGHT              ) return "MUON_HN_TIGHT";
 if(id == BaseSelection::MUON_FAKELOOSE             ) return "MUON_FAKELOOSE";
 if(id == BaseSelection::MUON_TOP_VETO              ) return "MUON_TOP_VETO";
 if(id == BaseSelection::MUON_TOP_LOOSE             ) return "MUON_TOP_LOOSE";
 if(id == BaseSelection::MUON_TOP_TIGHT             ) return "MUON_TOP_TIGHT";
 if(id == BaseSelection::MUON_PTETA                 ) return "MUON_PTETA";
 if(id == BaseSelection::MUON_NOCUT                 ) return "MUON_NOCUT";
 if(id == BaseSelection::PFJET_LOOSE                ) return "PFJET_LOOSE";
 if(id == BaseSelection::PFJET_MEDIUM               ) return "PFJET_MEDIUM";
 if(id == BaseSelection::PFJET_TIGHT                ) return "PFJET_TIGHT";
 if(id == BaseSelection::JET_HN                     ) return "JET_HN";
 if(id == BaseSelection::JET_HN_TChannel            ) return "JET_HN_TChannel";
 if(id == BaseSelection::JET_NOLEPTONVETO           ) return "JET_NOLEPTONVETO";
 if(id == BaseSelection::JET_LOOSE                  ) return "JET_LOOSE";
 if(id == BaseSelection::JET_TIGHT                  ) return "JET_TIGHT";
 if(id == BaseSelection::PHOTON_POG_LOOSE           ) return "PHOTON_POG_LOOSE";
 if(id == BaseSelection::PHOTON_POG_MEDIUM          ) return "PHOTON_POG_MEDIUM";
 if(id == BaseSelection::PHOTON_POG_TIGHT           ) return "PHOTON_POG_TIGHT";
 cerr << " ID [--] not found" << endl; exit(EXIT_FAILURE);

}


void AnalyzerCore::FillCorrectionHists(){

  string file = string(getenv("LQANALYZER_DIR")) + "/CATConfig/CattupleConfig/corrections_"+string(getenv("CATVERSION"))+".txt" ;
  cout << "Correction file = " << file << endl;
  ifstream corr_file(file.c_str());
  string sline;
  while(getline(corr_file,sline) ){
    std::istringstream is( sline );
    if(sline.find("###")!=string::npos) {  
      string tmp;
      is >> tmp;
      string label;
      is >> label;
      string dirname;
      is >> dirname;
      string filename;
      is >> filename;
      string histsname;
      is >> histsname;
      string histtype;
      is >> histtype;
      FillCorrectionHist(label,dirname,filename,histsname,histtype);
    }
    if(sline=="#END") break;
  }
}



void AnalyzerCore::FillCorrectionHist(string label, string dirname, string filename, string histsname, string histtype){

  TDirectory* origDir = gDirectory;

  if(TString(histtype).Contains("TGraphAsymmErrors")) {
    TFile *infile_sf = TFile::Open((string(getenv(dirname.c_str()))+ "/" + filename).c_str());
    CheckFile(infile_sf);
    TDirectory* tempDir = getTemporaryDirectory();
    tempDir->cd();
    TGraphAsymmErrors* tmp = (TGraphAsymmErrors*)infile_sf->Get(histsname.c_str())->Clone();
    CorrectionMapGraph[label] = tmp;
    infile_sf->Close();
    delete infile_sf;
    origDir->cd();
  }

  if(!TString(histtype).Contains("TH2")) return;
  TFile *infile_sf = TFile::Open((string(getenv(dirname.c_str()))+ "/" + filename).c_str());
  CheckFile(infile_sf);
  TDirectory* tempDir = getTemporaryDirectory();
  tempDir->cd();
  TH2F* tmp =  dynamic_cast<TH2F*> (( infile_sf->Get(histsname.c_str()))->Clone());
  CorrectionMap[label] = tmp;
  cout << "CorrectionMap["<<label <<"] = " << histsname << " (from " << getenv(dirname.c_str())<< "/" << filename<<")" << endl;
  infile_sf->Close();
  delete infile_sf;
  origDir->cd();
}




std::vector<snu::KJet> AnalyzerCore::GetJets(TString jetid,  bool smearjets,float ptcut, float etacut){
  
  std::vector<snu::KJet> jetColl;
  
  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsJet.find(jetid);
  if(it== selectionIDMapsJet.end()){
    cerr << "Jet ID ["+jetid+"] not found" << endl; exit(EXIT_FAILURE);
  }
  else {
    TString muontag="";
    TString eltag="";
    for (unsigned int i=0; i  < it->second.size(); i++){
      if ( it->second.at(i).first == "remove_near_muonID") muontag =  it->second.at(i).second;
      if ( it->second.at(i).first == "remove_near_electronID") eltag =  it->second.at(i).second;
    }
    if(smearjets){
     if (muontag.Contains("NONE") && eltag.Contains("NONE"))  eventbase->GetJetSel()->SelectJets(isData,jetColl,jetid, ptcut,etacut);
      else if (muontag.Contains("NONE") || eltag.Contains("NONE")) {    cerr << "cannot choose to remove jets near only one lepton" << endl; exit(EXIT_FAILURE);}
      else eventbase->GetJetSel()->SelectJets(isData,jetColl, GetMuons(muontag), GetElectrons(eltag) ,jetid, ptcut,etacut);
    
    }
    else{
      if (muontag.Contains("NONE") && eltag.Contains("NONE"))  eventbase->GetJetSel()->SelectJets(isData,jetColl,jetid, ptcut,etacut,false);
      else if (muontag.Contains("NONE") || eltag.Contains("NONE")) {    cerr << "cannot choose to remove jets near only one lepton" << endl; exit(EXIT_FAILURE);}
      else eventbase->GetJetSel()->SelectJets(isData,jetColl, GetMuons(muontag), GetElectrons(eltag) ,jetid, ptcut,etacut,false);

    }
  }

  return jetColl;
  
}

std::vector<snu::KMuon> AnalyzerCore::GetMuons(TString muid, float ptcut, float etacut){
  return GetMuons(muid, true, ptcut, etacut);
}

std::vector<snu::KMuon> AnalyzerCore::GetMuons(TString muid, bool keepfakes, float ptcut, float etacut){

  std::vector<snu::KMuon> muonColl;
  
  if(muid.Contains("NONE")) return muonColl;
  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsMuon.find(muid);
  if(it== selectionIDMapsMuon.end()){
    cerr << "Muon ID ["+muid+"] not found" << endl; exit(EXIT_FAILURE);
  }
  else {
    if (ptcut == -999.)  eventbase->GetMuonSel()->SelectMuons(muonColl,muid);
    else eventbase->GetMuonSel()->SelectMuons(muonColl,muid, ptcut, etacut);
  }
  //if(k_running_nonprompt) eventbase->GetMuonSel()->SelectMuons(muonColl, BaseSelection::MUON_HN_FAKELOOSE, 15., 2.4);
  //else eventbase->GetMuonSel()->SelectMuons(muonColl, BaseSelection::MUON_HN_TIGHT, 15., 2.4);
  
  return  GetTruePrompt(muonColl, keepfakes);
  
}


std::vector<snu::KElectron> AnalyzerCore::GetElectrons(TString elid,float ptcut, float etacut){
  return GetElectrons( true,  true, elid, ptcut, etacut);
}

std::vector<snu::KElectron> AnalyzerCore::GetElectrons(bool keepcf, bool keepfake, TString elid,float ptcut, float etacut){
  
  std::vector<snu::KElectron> electronColl;
  
  if(elid.Contains("NONE")) return electronColl;


  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsElectron.find(elid);
  if(it== selectionIDMapsElectron.end()){
    cerr << "Electron ID ["+elid+"] not found" << endl; exit(EXIT_FAILURE);
  }
  else {
    if (ptcut == -999.)  eventbase->GetElectronSel()->SelectElectrons(electronColl,elid);
    else eventbase->GetElectronSel()->SelectElectrons(electronColl,elid, ptcut, etacut);
  }
  
  //if(elid == "ELECTRON_HN_TIGHT"){
    /// This is the vector of electrons with optimie cuts
    //std::vector<snu::KElectron> _electronColl;
    //if(k_running_nonprompt) eventbase->GetElectronSel()->SelectElectrons(_electronColl, BaseSelection::ELECTRON_HN_FAKELOOSE_NOD0, 15., 2.5);
    //else eventbase->GetElectronSel()->SelectElectrons(_electronColl,BaseSelection::ELECTRON_HN_TIGHT, 15., 2.5);
    //electronColl =ShiftElectronEnergy(_electronColl, k_running_chargeflip);
  // }
   

  return  GetTruePrompt(electronColl, keepcf, keepfake); 

}




bool AnalyzerCore::HasCloseBJet(snu::KElectron el, KJet::Tagger tag, KJet::WORKING_POINT wp){

  std::vector<snu::KJet> alljets = GetJets("JET_NOLEPTONVETO");

  bool cl = false;
  for(unsigned int ij =0; ij < alljets.size(); ij++){
    if(el.DeltaR(alljets.at(ij)) < 0.5){
      if(alljets.at(ij).IsBTagged(tag,wp )) cl = true;
    }
  }

  return cl;
}


float AnalyzerCore::WeightCFEvent(std::vector<snu::KElectron> electrons, bool runchargeflip){

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



TDirectory* AnalyzerCore::getTemporaryDirectory(void) const
{

  
  gROOT->cd();
  TDirectory* tempDir = 0;
  int counter = 0;
  while (not tempDir) {
    // First, let's find a directory name that doesn't exist yet:               
    std::stringstream dirname;
    dirname << "AnalyzerCore_%i" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    // Let's try to make this directory:                                        
    tempDir = gROOT->mkdir((dirname.str()).c_str());

  }

  return tempDir;
}


double AnalyzerCore::MuonTrackingEffScaleFactor(vector<snu::KMuon> mu){
  float sf= 1.;
  if(isData) return 1.;
  for(vector<KMuon>::iterator itmu=mu.begin(); itmu!=mu.end(); ++itmu) {
    float mueta=itmu->Eta();
    if(CheckCorrectionGraph("TRACKING_EFF")){
      sf*= GetCorrectionGraph("TRACKING_EFF")->Eval(mueta);
    }
  }
  return sf;
}

double AnalyzerCore::MuonISOScaleFactor(TString muid, vector<snu::KMuon> mu,int sys){
  float sf= 1.;
  float sferr=1.;
  if(isData) return 1.;
  for(vector<KMuon>::iterator itmu=mu.begin(); itmu!=mu.end(); ++itmu) {
    float mupt=itmu->Pt();
    if(itmu->Pt() >120. )mupt=119.;
    if(itmu->Pt() < 20.) mupt=21.;
    if(CheckCorrectionHist("ISO_" + muid)){
	sferr = double(sys)*GetCorrectionHist("ISO_" + muid)->GetBinError( GetCorrectionHist("ISO_" + muid)->FindBin( fabs(itmu->Eta()), mupt) );
	
	sf*= (1. + sferr)*GetCorrectionHist("ISO_" + muid)->GetBinContent( GetCorrectionHist("ISO_" + muid)->FindBin( fabs(itmu->Eta()), mupt) );
    }
    //else{ m_logger << ERROR << "No ISO scalefactors are available for "<< muid << LQLogger::endmsg; exit(0);}  
  }
  return sf;
}

double AnalyzerCore::MuonScaleFactor(TString muid, vector<snu::KMuon> mu,int sys){
  float sf= 1.;
  float sferr=1.;
  if(isData) return 1.;
  if(mu.size() == 0) return 1.;
  for(vector<KMuon>::iterator itmu=mu.begin(); itmu!=mu.end(); ++itmu) {
    float mupt=itmu->Pt();
    if(itmu->Pt() <20.) mupt= 21.;
    if(itmu->Pt() >120.) mupt= 119.;
    if(CheckCorrectionHist("ID_" + muid)){
      sferr = double(sys)*GetCorrectionHist("ID_" + muid)->GetBinError( GetCorrectionHist("ID_" + muid)->FindBin( fabs(itmu->Eta()), mupt) );
      
      sf*=  (1. + sferr)* GetCorrectionHist("ID_" + muid)->GetBinContent( GetCorrectionHist("ID_" + muid)->FindBin( fabs(itmu->Eta()), mupt) );
    }
  }

  return sf;
}

bool AnalyzerCore::TriggerMatch(TString trigname, vector<snu::KMuon> mu){
  
  if(mu.size() == 2){
    if(!mu.at(0).TriggerMatched(trigname)) return false;
    if(!mu.at(1).TriggerMatched(trigname)) return false;
  }
  return true;
}

bool AnalyzerCore::Is2015Analysis(){
  if(TString(getenv("running2015")).Contains("True")) return true;
  else return false;
}

double AnalyzerCore::TriggerScaleFactor( vector<snu::KElectron> el, vector<snu::KMuon> mu,  TString trigname, int direction){
  
  if(isData) return 1.;
  
  //// ONLY SINGLE MUON TRIGGER SCAKE FACTOR ARE UPDATED, OTHERS ARE 2015 SFS

  ///https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgHLTScaleFactorMeasurements
  
  if(Is2015Analysis()){
    if(el.size() == 2){
      //if (trigname.Contains("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) return 0.997*0.997*0.998;
      if (trigname.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v")) return 0.995*0.998;
      if (trigname.Contains("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v")) return 0.997;
      if (trigname.Contains("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) {
	float aeta1 = fabs(el.at(0).SCEta());
	float aeta2 = fabs(el.at(1).SCEta());
	if ( aeta1 < 1.2 ) {
	  if      ( aeta2 < 1.2 ) return 0.953 + direction*0.023;
	  else if ( aeta2 < 2.4 ) return 0.957 + 0.5*direction*0.027;
	}
	else if ( aeta1 < 2.4 ) {
	  if      ( aeta2 < 1.2 ) return 0.967 + 0.5*direction*0.029;
	  else if ( aeta2 < 2.4 ) return 0.989 - direction*0.031;
	}
      }
    }    
  
  if(mu.size() == 1){
    if (trigname.Contains("HLT_IsoMu20") || trigname.Contains("HLT_IsoTkMu20"))  {
      ///https://indico.cern.ch/event/462268/contributions/1979019/attachments/1188638/1724574/2015.11.17_MuonPOG_SingleMuTrigEff_SF_KPLee_v2.pdf
      
      float mupt=mu.at(0).Pt();
      if(mupt > 120.) mupt = 119.;
 if(mupt < 20.) mupt = 21.;
      if(eventbase->GetEvent().RunNumber()  < 255032){
	if(CheckCorrectionHist("MUON_TRIGGER_C")){
	  float sferr = double(direction)*GetCorrectionHist("MUON_TRIGGER_C")->GetBinError(GetCorrectionHist("MUON_TRIGGER_C")->FindBin( fabs(mu.at(0).Eta()), mupt) );
	  return  (1. + sferr)*GetCorrectionHist("MUON_TRIGGER_C")->GetBinContent( GetCorrectionHist("MUON_TRIGGER_C")->FindBin(  fabs(mu.at(0).Eta()), mupt) );
	}
      }
      else if(eventbase->GetEvent().RunNumber()  < 257820){
	if(CheckCorrectionHist("MUON_TRIGGER_D1")){
	  float sferr = double(direction)*GetCorrectionHist("MUON_TRIGGER_D1")->GetBinError( GetCorrectionHist("MUON_TRIGGER_D1")->FindBin( fabs(mu.at(0).Eta()), mupt) );
	  return  (1. + sferr)*GetCorrectionHist("MUON_TRIGGER_D1")->GetBinContent( GetCorrectionHist("MUON_TRIGGER_D1")->FindBin(  fabs(mu.at(0).Eta()), mupt) );
	}
      }
      else {
	if(CheckCorrectionHist("MUON_TRIGGER_D2")){
	  float sferr = double(direction)*GetCorrectionHist("MUON_TRIGGER_D2")->GetBinError( GetCorrectionHist("MUON_TRIGGER_D2")->FindBin( fabs(mu.at(0).Eta()), mupt) );
	  return  (1. + sferr)*GetCorrectionHist("MUON_TRIGGER_D2")->GetBinContent( GetCorrectionHist("MUON_TRIGGER_D2")->FindBin(  fabs(mu.at(0).Eta()), mupt) );
	}
      }
    }
  }
  if(el.size() == 1){

    // http://fcouderc.web.cern.ch/fcouderc/EGamma/scaleFactors/moriond2016_76X/eleMVATrig/ScaleFactor_GsfElectronToRECO_passingTrigWP80.txt.egamma_egammaPlots.pdf
    if (trigname.Contains("HLT_Ele23_WPLoose_Gsf_v")) {
      float el_eta = el.at(0).SCEta();
      float el_pt = el.at(0).Pt();
      if(el_pt < 25.){
	if(el_eta < -1.56) return 0.928;
	else if(el_eta < -1.4442) return 0.949;
	else if(el_eta < -0.8) return 0.964;
	else if(el_eta < 0.) return 0.963;
	else if(el_eta < 0.8) return 0.985;
	else if(el_eta < 1.4442) return 0.966;
	else if(el_eta < 1.56) return 0.903;
	else return 0.906;
      }
      else  if(el_pt < 35.){
	if(el_eta < -1.56) return 0.937;
	else if(el_eta < -1.4442) return 0.973;
	else if(el_eta < -0.8) return 0.983;
	else if(el_eta < 0.) return 0.985;
	else if(el_eta < 0.8) return 1.005;
        else if(el_eta < 1.4442) return 0.984;
	else if(el_eta < 1.56) return 0.931;
	else return 0.920;
      } 
      else  if(el_pt < 45.){
	if(el_eta < -1.56) return 0.953;
	else if(el_eta < -1.4442) return 0.966;
	else if(el_eta < -0.8) return 0.985;
	else if(el_eta < 0.) return 0.986;
	else if(el_eta < 0.8) return 1.001;
        else if(el_eta < 1.4442) return 0.986;
	else if(el_eta < 1.56) return 0.955;
	else return 0.941;
      } 
      else  if(el_pt < 55.){
	if(el_eta < -1.56) return 0.963;
	else if(el_eta < -1.4442) return 0.972;
	else if(el_eta < -0.8) return 0.987;
	else if(el_eta < 0.) return 0.990;
	else if(el_eta < 0.8) return 1.003;
        else if(el_eta < 1.4442) return 0.986;
	else if(el_eta < 1.56) return 0.950;
	else return 0.950;
      } 
      else{
	if(el_eta < -1.56) return 0.973;
	else if(el_eta < -1.4442) return 0.975;
	else if(el_eta < -0.8) return 0.983;
	else if(el_eta < 0.) return 0.992;
	else if(el_eta < 0.8) return 1.;
        else if(el_eta < 1.4442) return 0.985 ;
	else if(el_eta < 1.56) return 0.991;
	else return 0.958;
      }
      
    }
  } 

  if(mu.size() == 2){
    //if (trigname.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")) return (0.993*0.984*0.980*0.968);
    //if (trigname.Contains("HLT_Mu17_TrkIsoVVL_Mu8_OR_TkMu8_TrkIsoVVL_DZ_v")) return 0.982*0.985*.973;
    //if (trigname.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")) return (0.993*0.984*0.980*0.968);
    //  new SF frmo AN2015_309_v7.pdf
    if( (trigname.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")) || (trigname.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"))){
      float aeta1 = fabs(mu.at(0).Eta());
      float aeta2 = fabs(mu.at(1).Eta());
      
      if ( aeta1 < 1.2 ) {
	if      ( aeta2 < 1.2 ) return 0.926 + direction*0.022;
	else if ( aeta2 < 2.4 ) return 0.943 + 0.5*direction*0.026;
      }
      else if ( aeta1 < 2.4 ) {
	if      ( aeta2 < 1.2 ) return 0.958 + 0.5*direction*0.027;
	else if ( aeta2 < 2.4 ) return 0.926 - direction*0.027;
      }
    }
  }
       
  if( (el.size() == 1) && (mu.size() == 1)){
    if ((trigname.Contains("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")) ||  (trigname.Contains("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v"))){
      float aeta1 = fabs(mu.at(0).Eta());
      float aeta2 = fabs(el.at(0).SCEta());
      
      if ( aeta1 < 1.2 ) {
	if      ( aeta2 < 1.2 ) return 0.966 + direction*0.022;
	else if ( aeta2 < 2.4 ) return 0.973 + 0.5*direction*0.025;
      }
      else if ( aeta1 < 2.4 ) {
	if      ( aeta2 < 1.2 ) return 0.981 + 0.5*direction*0.025;
	else if ( aeta2 < 2.4 ) return 0.984 - direction*0.030;
      }
    }
  }
  
  
  /// Trilepton
  if (trigname.Contains("HLT_TripleMu_12_10_5")) return  0.992*0.986*0.981* 0.982;
  if (trigname.Contains("HLT_DiMu9_Ele9_CaloIdL_TrackIdL")) return  0.981*0.95;
  if (trigname.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL")) return 0.966*1.005;
  if (trigname.Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL")) return 0.990*1.003*0.991*0.998;
  
  }
  
  if(el.size() == 2){
    return 1.;
  }

  if(mu.size() == 1){
    if (trigname.Contains("HLT_IsoMu22") || trigname.Contains("HLT_IsoTkMu22"))  {
      /// https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults
      
     return 1.;

    }
  }
  if(el.size() == 1){
    return 1.;
  } 

  if(mu.size() == 2){
    return 1.;
  }
  
  
  

  return 1.;
  
}



double AnalyzerCore::ElectronScaleFactor( TString elid, vector<snu::KElectron> el, int sys){
  float sf= 1.;
  if(isData) return 1.;
  
  std::string sid= "";
  
  for(vector<KElectron>::iterator itel=el.begin(); itel!=el.end(); ++itel) {
    float elpt=itel->Pt();
    if(elpt > 200.) elpt= 199.;
    if(elpt < 20.) elpt= 21.;
    
    if(CheckCorrectionHist("EL_GSF")){
      sf *= GetCorrectionHist("EL_GSF")->GetBinContent(GetCorrectionHist("EL_GSF")->FindBin(fabs(itel->SCEta()), elpt));
    }
    
    if(CheckCorrectionHist("ID_" + elid)){
	int bin =  GetCorrectionHist("ID_" + elid)->FindBin(fabs(itel->SCEta()), elpt);
	sf *= GetCorrectionHist("ID_" + elid)->GetBinContent(bin);
    }
  }

 
  return sf;
}

double AnalyzerCore::ElectronRecoScaleFactor(vector<snu::KElectron> el){
  
  if(isData) return 1.;
  float sf= 1.;
  for(vector<KElectron>::iterator itel=el.begin(); itel!=el.end(); ++itel) {
    float elpt= itel->Pt() ;
    if(itel->Pt() > 200.) elpt=119.;
    if(itel->Pt() < 20.) elpt=21.;
    if(CheckCorrectionHist("EL_RECO")){
      int bin =  GetCorrectionHist("EL_RECO")->FindBin(fabs(itel->SCEta()), elpt);
      sf *= GetCorrectionHist("EL_RECO")->GetBinContent(bin);
    }
  }
  
  return sf;
}

float AnalyzerCore::TempPileupWeight(){
  
  if(!TString(getenv("CATVERSION")).Contains("v7-6-6")){
    if(isData) return 1.;
    else {
      if(TString(getenv("PUFILE")).Contains(".root")){
	return reweightPU->GetWeight(eventbase->GetEvent().nVertices(),TString(getenv("CATVERSION")));
      }
    }
  }
  else{
    if(isData) return 1.; 
    return eventbase->GetEvent().AltPileUpWeight();
  }
  return 1.;
}

float AnalyzerCore::WeightByTrigger(vector<TString> triggernames, float tlumi){

  if(isData){
    for(unsigned int i=0; i < triggernames.size() ; i++){
      //// code here sets weight to -99999. if user tries to use incorrect dataset in data 
      /// datasets for each trigger in 2015 can be seen in following googledoc:
      //// https://docs.google.com/spreadsheets/d/1BkgAHCC4UtP5sddTZ5G5iWY16BxleuK7rqT-Iz2LHiM/pubhtml?gid=0&single=true
       if(!k_channel.Contains("DoubleMuon")){
	if(triggernames.at(i).Contains("HLT_Mu17_Mu8_DZ_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_Mu8_SameSign_DZ_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu20_Mu10_SameSign_DZ_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu8_TrkIsoVVL_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_TripleMu_12_10_5")) return -99999.;
      }
      if(!k_channel.Contains("SingleMuon")){
	if(triggernames.at(i).Contains("HLT_Mu8_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu20_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu50_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoMu22_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoTkMu22_")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoMu20_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoTkMu20_")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoMu24_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoTkMu24_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu24_eta2p1_v" )) return -99999.;
      }
      if(!k_channel.Contains("DoubleEG")){
	if(triggernames.at(i).Contains("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v")) return -99999.;

      }
      if(!k_channel.Contains("SingleElectron")){
	// Single Electon                                                                                                                                                       
	if(triggernames.at(i).Contains("HLT_Ele23_WPLoose_Gsf_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele27_WPTight_Gsf_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele27_eta2p1_WPLoose_Gsf_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele25_eta2p1_WPTight_Gsf_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele45_WPLoose_Gsf_v")) return -99999.;
      }
      if(!k_channel.Contains("MuonEG")){
	if(triggernames.at(i).Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v")) return -99999.;
      }
    }
    return 1.;
  }


  float trigps= -1.;
  for(unsigned int i=0; i < triggernames.size() ; i++){
    if(WeightByTrigger(triggernames.at(i), tlumi) > trigps) trigps = WeightByTrigger(triggernames.at(i),tlumi) ;
  }

  //if(trigps  > 1.) {m_logger << ERROR << "Error in getting weight for trigger prescale. It cannot be > 1, this means trigger lumi >> total lumi"  << LQLogger::endmsg; exit(0);}
  if(trigps  < 0.) {m_logger << ERROR << "Error in getting weight for trigger " << triggernames.at(0) << " prescale. It cannot be < 0, this means trigger lumi >> total lumi"  << LQLogger::endmsg; exit(0);}
  
  return trigps;
 
}

float AnalyzerCore::WeightByTrigger(TString triggername, float tlumi){

  /// Function applies weight to MC 
  /// Depends on trigger 
  /// 

  if(isData) return 1.;
  //  brilcalc lumi -u /pb 
  // --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/moriond16_normtag.json 
  // -i jsonfiles/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt 
  //--hltpath ""
  
  /// NUMBERS FROM GetLumi_Triggers.py SCRIPT
  
  /// In v766 path lumi is corrected for removal of bad beamspot LS
  // https://github.com/vallot/CATTools/commit/aae3e60b194b1bacf2595a33c8fa27f411dac16b
  for(map<TString, float>::iterator mit = trigger_lumi_map_cat2016.begin(); mit != trigger_lumi_map_cat2016.end(); mit++){
    if(triggername.Contains(mit->first)) return (mit->second / tlumi);
  }
  m_logger << ERROR << "Error in getting weight for trigger  " << triggername << "  prescale. Trigname is not correct or not in map"  << LQLogger::endmsg; exit(0);

  return 1.;
}


void AnalyzerCore::AddTriggerToList(TString triggername){
  
  triggerlist.push_back(triggername);
}

AnalyzerCore::~AnalyzerCore(){
  

  Message("In AnalyzerCore Destructor" , INFO);

  if(TString(getenv("PUFILE")).Contains(".root")){
    if(!k_isdata)  delete reweightPU;
  }
  if(FRHist) delete FRHist;

  for(map<TString, TH1*>::iterator it = maphist.begin(); it!= maphist.end(); it++){
    delete it->second;
  }
  maphist.clear();

  for(map<TString, TH2*>::iterator it = maphist2D.begin(); it!= maphist2D.end(); it++){
    delete it->second;
  }
  maphist2D.clear();



  for(map<TString, MuonPlots*>::iterator it = mapCLhistMu.begin(); it != mapCLhistMu.end(); it++){
    delete it->second;
  }
  mapCLhistMu.clear();
  

  for(map<TString, JetPlots*>::iterator it = mapCLhistJet.begin(); it != mapCLhistJet.end(); it++){
    delete it->second;
  }
  mapCLhistJet.clear();

  for(map<TString, ElectronPlots*>::iterator it = mapCLhistEl.begin(); it != mapCLhistEl.end(); it++){
    delete it->second;
  }
  mapCLhistEl.clear();

  for(map<TString, SignalPlotsEE*>::iterator it = mapCLhistSigEE.begin(); it != mapCLhistSigEE.end(); it++){
    delete it->second;
  }
  mapCLhistSigEE.clear();

  for(map<TString, SignalPlotsMM*>::iterator it = mapCLhistSigMM.begin(); it != mapCLhistSigMM.end(); it++){
    delete it->second;
  }
  mapCLhistSigMM.clear();


  for(map<TString, SignalPlotsEM*>::iterator it = mapCLhistSigEM.begin(); it != mapCLhistSigEM.end(); it++){
    delete it->second;
  }
  mapCLhistSigEM.clear();

  
  for(map<TString, TriLeptonPlots*>::iterator it = mapCLhistTriLep.begin(); it != mapCLhistTriLep.end(); it++){
    delete it->second;
  }
  mapCLhistTriLep.clear();


  for(map<TString,TNtupleD*>::iterator it = mapntp.begin(); it!= mapntp.end(); it++){ 
    delete it->second;
  }
  mapntp.clear();
  CorrectionMap.clear();
  CorrectionMapGraph.clear();

  delete m_fakeobj;
  if (TString(getenv("CATVERSION")).Contains("v7-6-6")){
    delete rmcor;
  }
  
  for(map<TString,TH2F*>::iterator it = CorrectionMap.begin(); it != CorrectionMap.end(); it++){
    delete it->second;
  }
  for(map<TString,TGraphAsymmErrors*>::iterator it = CorrectionMapGraph.begin(); it != CorrectionMapGraph.end(); it++){
    delete it->second;
  }

  for(std::map<TString,BTagSFUtil*>::iterator it = MapBTagSF.begin(); it!= MapBTagSF.end(); it++){
    delete it->second;
  }
  MapBTagSF.clear();

}

//###
//###   IMPORTANT BASE FUNCTION: SETS UP EVENT FOR ALL CYCLES
//###

void AnalyzerCore::SetUpEvent(Long64_t entry, float ev_weight) throw( LQError ) {
  
  Message("In SetUpEvent(Long64_t entry) " , DEBUG);
  m_logger << DEBUG << "This is entry " << entry << LQLogger::endmsg;
  if (!fChain) throw LQError( "Chain is not initialized",  LQError::SkipCycle );     
  
  if(LQinput){

    m_logger << DEBUG << "k_isdata = " << k_isdata << " and isData = " << isData << LQLogger::endmsg;
    if(k_isdata != isData) throw LQError( "!!! Event is confused. It does not know if it is data or MC", LQError::SkipCycle );
  }
  else isData = k_isdata;
  
  if (!(entry % output_interval)) {
    m_logger << INFO <<  "Processing entry " << entry <<  "/" << nentries << LQLogger::endmsg;

  }
  
  snu::KEvent eventinfo = GetEventInfo();
  
  if(k_isdata){
    if(ev_weight!=1.) Message("ERROR in setting weights. This is Data...", INFO);
    MCweight=1.;
    weight = 1.;
  }
  else {
    MCweight = eventinfo.MCWeight(); //Get MC weight here FIX ME                                                              
    weight= ev_weight; 
  }

 //
  // creates object that stores all SKTree classes	
  //                                                                                                        
  
  snu::KTrigger triggerinfo = GetTriggerInfo(triggerlist);
    
  std::vector<snu::KJet> skjets= GetAllJets();
  std::vector<snu::KGenJet> skgenjets=GetAllGenJets();
  
   
  /// np == numberof particles you want to store at truth info. 30 is default unless running nocut sktree OR signal
  int np =  AssignnNumberOfTruth();
  
  LQEvent lqevent(GetAllMuons(), GetAllElectrons(), GetAllPhotons(), skjets, skgenjets,GetTruthParticles(np), triggerinfo,eventinfo);
  
  //  eventbase is master class to use in analysis 
  //
  
  eventbase = new EventBase(lqevent);
    
  eventbase->GetElectronSel()->SetIDSMap(selectionIDMapsElectron);
  eventbase->GetElectronSel()->SetIDFMap(selectionIDMapfElectron);
  eventbase->GetMuonSel()->SetIDSMap(selectionIDMapsMuon);
  eventbase->GetMuonSel()->SetIDFMap(selectionIDMapfMuon);
  eventbase->GetJetSel()->SetIDSMap(selectionIDMapsJet);
  eventbase->GetJetSel()->SetIDFMap(selectionIDMapfJet);

  if(!k_isdata){
    if(!changed_target_lumi){
      changed_target_lumi=true;
    }
  }

}


int AnalyzerCore::VersionStamp(TString cversion){
  
  if(cversion.Contains("v7-4-4")) return 1;
  else if(cversion.Contains("v7-4-5")) return 2;
  else if(cversion.Contains("v7-6-2") || cversion.Contains("v7-6-3") || cversion.Contains("v7-6-4")   ) return 3;
  else if((cversion.Contains("v7-6-5") || cversion.Contains("v7-6-6"))) return 4;
  else if((cversion.Contains("v8-0-1"))) return 5;
  else if((cversion.Contains("v8-0-2"))) return 6;
  
  return 5;
 
}

int AnalyzerCore::AssignnNumberOfTruth(){
  int np = 1000;
  if(k_classname.Contains("SKTreeMaker")) np = 1000;
  if(k_classname.Contains("SKTreeMakerDiLep")) np = 0;
  if(k_classname.Contains("SKTreeMakerTriLep")) np = 0;

  if(k_classname.Contains("SKTreeMaker")){
    if(k_sample_name.Contains("QCD") && !k_sample_name.Contains("mad")) np = 0;
  }

  /// List of signal samples
  /// G.Yu needs to add signal here
  
  if(IsSignal()) np = 1000;
  
  return np;
}

bool AnalyzerCore::IsSignal(){
  
  if(k_sample_name.Contains("Majornana")) return true;
  if(k_sample_name.Contains("HN")) return true;
  return false;
}


void AnalyzerCore::ClassInfo(){
  
  /*if(eventinfo.CatVersion().empty()){ 
    m_logger << INFO << "Catuple version is v7-4-X. Only basic infomation is available." << LQLogger::endmsg;
    
  }
    
  else if(TString(eventinfo.CatVersion()).Contains("v7-6-2")){
    m_logger << INFO <<  "Running on catuples version " << eventinfo.CatVersion() << LQLogger::endmsg;
    
    
    }*/

}

float AnalyzerCore::SumPt( std::vector<snu::KJet> particles){

  float sumpt=0.;
  for(std::vector<snu::KJet>::iterator it = particles.begin(); it != particles.end(); it++){
    sumpt += it->Pt();
  }
  return sumpt;
}
  


std::vector<snu::KElectron> AnalyzerCore::ShiftElectronEnergy(std::vector<snu::KElectron> electrons, bool applyshift){
  
  std::vector<snu::KElectron> shiftedel;

  for(unsigned int iel=0; iel < electrons.size(); iel++){
    float scale =0.98;

    if(applyshift)electrons.at(iel).SetPtEtaPhiM(electrons.at(iel).Pt()*scale, electrons.at(iel).Eta(), electrons.at(iel).Phi(), 0.511e-3);
    shiftedel.push_back(electrons.at(iel));
  }    
  return shiftedel;
}



bool AnalyzerCore::IsDiEl(){
  if(isData) return false;
  int iel(0);
  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){
    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;
    if(fabs(eventbase->GetTruth().at(ig).PdgId()) == 11) iel++;
  }
  if(iel >1) return true;
  else return false;
}
void AnalyzerCore::TruthPrintOut(){
  if(isData) return;
  m_logger << INFO<< "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  cout << "Particle Index |  PdgId  | GenStatus   | Mother PdgId |  Part_Eta | Part_Pt | Part_Phi | Mother Index |   " << endl;
  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){
    
    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;
    if (eventbase->GetTruth().at(ig).PdgId() == 2212)  cout << ig << " | " << eventbase->GetTruth().at(ig).PdgId() << "  |               |         |        |         |        |         |" << endl;

    cout << ig << " |  " <<  eventbase->GetTruth().at(ig).PdgId() << " |  " << eventbase->GetTruth().at(ig).GenStatus() << " |  " << eventbase->GetTruth().at(eventbase->GetTruth().at(ig).IndexMother()).PdgId()<< " |   " << eventbase->GetTruth().at(ig).Eta() << " | " << eventbase->GetTruth().at(ig).Pt() << " | " << eventbase->GetTruth().at(ig).Phi() << " |   " << eventbase->GetTruth().at(ig).IndexMother()  << endl;
  }
}


bool AnalyzerCore::isPrompt(long pdgid) {
  /// mother pdgid
  pdgid = abs(pdgid);
  if (pdgid == 24) return true; // Z
  else if (pdgid == 23) return true; // W
  else if (pdgid == 15) return true; // taus
  else if (pdgid == 90) return true; // N
  else return false;
}

void AnalyzerCore::EndEvent()throw( LQError ){
  
  delete eventbase;                                                                                                            

}
  
void AnalyzerCore::CheckFile(TFile* file)throw( LQError ){

  if(file) m_logger << INFO << "Analyzer: File " << file->GetName() << " was found." << LQLogger::endmsg;
  else m_logger  << INFO <<"Analyzer  " << file->GetName()  << "  : ERROR Rootfile failed to open." << LQLogger::endmsg;

  if(!file) exit(EXIT_FAILURE);
  return;
}


void AnalyzerCore::ListTriggersAvailable(){
  cout << "Set of triggers you can use are: " << endl;
  for(unsigned int i=0; i < eventbase->GetTrigger().GetHLTInsideDatasetTriggerNames().size(); i++){
    cout << eventbase->GetTrigger().GetHLTInsideDatasetTriggerNames().at(i)<< " has prescale " << eventbase->GetTrigger().GetHLTInsideDatasetTriggerPrescales().at(i)<< endl;
  }
  return;
}





float AnalyzerCore::TriggerEff(TString trigname,  std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons){

  Message("In TriggerEff", DEBUG);

  if(isData){
    //if(PassTrigger(trigname)) return 1.;
    return 1.;
  }
  if(electrons.size() >= 1 && muons.size() >= 1){
    float trig_eff(1.);
    if(trigname.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")){
      trig_eff*= GetEff(electrons.at(0), "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v");
      trig_eff*= GetEff(muons.at(0), "HLT_Mu8");
      return trig_eff;
    }
    else  return 0.;
  }
  else return 0.;
  Message("END CorrectMuonMomentum", DEBUG);

  return 1.;
}



float AnalyzerCore::TriggerEff(TString trigname, std::vector<snu::KElectron> electrons){
  
  if(isData){
    //if(PassTrigger(trigname)) return 1.;
    return 1.;
  }
  
  if(electrons.size() >=2){
    float trig_eff(1.);
    if(trigname.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")){
      
      trig_eff*= GetEff(electrons.at(0), "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v");
      trig_eff*= GetEff(electrons.at(1), "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v");
      trig_eff *= 0.995; // DZ efficiency AN2016_228
      return trig_eff;
    }
    else  return 0.;
  }
  else return 0.;  

  return 1.;
}



float AnalyzerCore::GetEff(snu::KElectron el, TString trigname){
  float sceta = fabs(el.SCEta());
  float pt = el.Pt();
  if (trigname == "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v"){
    if(sceta < 0.8){
      if( pt < 20)return 0.005;
      else if( pt < 25)return 0.418;
      else if( pt < 30)return 0.964;
      else if( pt < 35)return 0.975;
      else if( pt < 40)return 0.981;
      else if( pt < 50)return 0.984;
      else if( pt < 75)return 0.986;
      else if( pt < 100)return 0.988;
      else if( pt < 200)return 0.985;
      else if( pt < 300)return 0.978;
      else  return 0.966;
    }
    else   if(sceta < 1.477){
      if( pt < 20)return 0.005;
      else if( pt < 25)return 0.357;
      else if( pt < 30)return 0.966;
      else if( pt < 35)return 0.981;
      else if( pt < 40)return 0.985;
      else if( pt < 50)return 0.987;
      else if( pt < 75)return 0.988;
      else if( pt < 100)return 0.990;
      else if( pt < 200)return 0.988;
      else if( pt < 300)return 0.991;
      else return 0.993;
    }
    else   if(sceta < 2.0){
      if( pt < 20)return 0.12;
      else if( pt < 25)return 0.372;
      else if( pt < 30)return 0.980;
      else if( pt < 35)return 0.991;
      else if( pt < 40)return 0.993;
      else if( pt < 50)return 0.995;
      else if( pt < 75)return 0.995;
      else if( pt < 100)return 0.996;
      else if( pt < 200)return 0.997;
      else if( pt < 300)return 0.994;
      else if( pt < 2000)return 1.00;
    }
    else {
      if( pt < 20)return 0.003;
      else if( pt < 25)return 0.198;
      else if( pt < 30)return 0.939;
      else if( pt < 35)return 0.984;
      else if( pt < 40)return 0.989;
      else if( pt < 50)return 0.990;
      else if( pt < 75)return 0.993;
      else if( pt < 100)return 0.997;
      else if( pt < 200)return 0.996;
      else if( pt < 300)return 1.00;
      else if( pt < 2000)return 1.00;
    }
  }
  
  if (trigname == "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v"){
    if(sceta < 0.8){
      if( pt < 10)return 0.630;
      else if( pt < 20)return 0.946;
      else if( pt < 25)return 0.959;
      else if( pt < 30)return 0.966;
      else if( pt < 35)return 0.975;
      else if( pt < 40)return 0.981;
      else if( pt < 50)return 0.984;
      else if( pt < 75)return 0.986;
      else if( pt < 100)return 0.988;
      else if( pt < 200)return 0.985;
      else if( pt < 300)return 0.978;
      else return 0.966;
    }
    else   if(sceta < 1.477){
      if( pt < 10)return 0.578;
      else if( pt < 20)return 0.953;
      else if( pt < 25)return 0.971;
      else if( pt < 30)return 0.976;
      else if( pt < 35)return 0.981;
      else if( pt < 40)return 0.985;
      else if( pt < 50)return 0.987;
      else if( pt < 75)return 0.988;
      else if( pt < 100)return 0.990;
      else if( pt < 200)return 0.988;
      else if( pt < 300)return 0.991;
      else  return 0.993;
    }
    else   if(sceta < 2.0){
      if( pt < 10)return 0.517;
      else if( pt < 20)return 0.934;
      else if( pt < 25)return 0.984;
      else if( pt < 30)return 0.991;
      else if( pt < 35)return 0.991;
      else if( pt < 40)return 0.993;
      else if( pt < 50)return 0.995;
      else if( pt < 75)return 0.995;
      else if( pt < 100)return 0.996;
      else if( pt < 200)return 0.997;
      else if( pt < 300)return 0.994;
      else  return 1.0;
    }
    else {
      if( pt < 10)return 0.516;
      else if( pt < 20)return 0.934;
      else if( pt < 25)return 0.966;
      else if( pt < 30)return 0.977;
      else if( pt < 35)return 0.984;
      else if( pt < 40)return 0.989;
      else if( pt < 50)return 0.990;
      else if( pt < 75)return 0.993;
      else if( pt < 100)return 0.997;
      else if( pt < 200)return 0.996;
      else if( pt < 300)return 1.00;
      else return 1.00;
    }
  }
  
  return 1.;
}



float AnalyzerCore::TriggerEff(TString trigname, std::vector<snu::KMuon> muons){
  
  if(isData){
    ///    if(PassTrigger(trigname)) return 1.;
    return 1.;
  }
  
  if(muons.size() >= 2){
    float trig_eff(1.);
    if(trigname.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")){
      trig_eff*= GetEff(muons.at(0), "HLT_Mu17");
      trig_eff*= GetEff(muons.at(1), "HLT_Mu8");
      return trig_eff;
    }
    else return 0.;
  }
  
  if(muons.size() >= 1){
    if (trigname.Contains("HLT_IsoMu22") || trigname.Contains("HLT_IsoTkMu22"))  {
      /// https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults                                                                                                                                                            
      float mupt=muons.at(0).Pt();
      if(mupt > 500.) mupt = 499.;
      if(mupt < 25.) mupt = 25.;
      if(eventbase->GetEvent().RunNumber()  < 274093){
	return  GetCorrectionHist("SingleMuon_274093")->GetBinContent( GetCorrectionHist("SingleMuon_274093")->FindBin(  fabs(muons.at(0).Eta()), mupt) );
      }
      else {
	return GetCorrectionHist("SingleMuon_276097")->GetBinContent( GetCorrectionHist("SingleMuon_276097")->FindBin(  fabs(muons.at(0).Eta()), mupt) );
      }
    }
    else return 0.;
  }
  else return 0.;
  return 1.;
}

bool AnalyzerCore::CheckCorrectionHist(TString label){
  map<TString, TH2F*>::iterator  mapit = CorrectionMap.find(label);
  if (mapit!= CorrectionMap.end()){
    return true;
  }
  else return false;
}
bool AnalyzerCore::CheckCorrectionGraph(TString label){
  map<TString, TGraphAsymmErrors*>::iterator  mapit = CorrectionMapGraph.find(label);
  if (mapit!= CorrectionMapGraph.end()){
    return true;
  }
  else return false;
}

TGraphAsymmErrors* AnalyzerCore::GetCorrectionGraph(TString label){
  map<TString, TGraphAsymmErrors*>::iterator mapit = CorrectionMapGraph.find(label);
  if (mapit!= CorrectionMapGraph.end()){
    return mapit->second;
  }
  else{
    m_logger << ERROR << "Could not find correction file with label label" <<  LQLogger::endmsg;
    exit(0);
  }
}

TH2F* AnalyzerCore::GetCorrectionHist(TString label){
  map<TString, TH2F*>::iterator mapit = CorrectionMap.find(label);
  if (mapit!= CorrectionMap.end()){
    return mapit->second;
  }
  else{
    m_logger << ERROR << "Could not find correction file with label label" <<  LQLogger::endmsg;
    exit(0);
  }
}


float AnalyzerCore::GetEff(snu::KMuon mu, TString trigname){
  float eta = fabs(mu.Eta());
  float pt = mu.Pt();
  if (trigname == "HLT_Mu17"){
    if(eta < 0.8){
      if( pt < 20)return 0.612;
      else if( pt < 25)return 0.930;
      else if( pt < 30)return 0.931;
      else if( pt < 35)return 0.932;
      else if( pt < 40)return 0.932;
      else if( pt < 50)return 0.932;
      else if( pt < 75)return 0.929;
      else if( pt < 100)return 0.923;
      else if( pt < 200)return 0.920;
      else if( pt < 300)return 0.901;
      else  return 0.899;
    }
    else   if(eta < 1.25){
      if( pt < 20)return 0.591;
      else if( pt < 25)return 0.914;
      else if( pt < 30)return 0.917;
      else if( pt < 35)return 0.916;
      else if( pt < 40)return 0.916;
      else if( pt < 50)return 0.917;
      else if( pt < 75)return 0.914;
      else if( pt < 100)return 0.911;
      else if( pt < 200)return 0.903;
      else if( pt < 300)return 0.871;
      else return 0.901;
    }
    else   if(eta < 1.6){
      if( pt < 20)return 0.603;
      else if( pt < 25)return 0.947;
      else if( pt < 30)return 0.956;
      else if( pt < 35)return 0.960;
      else if( pt < 40)return 0.962;
      else if( pt < 50)return 0.961;
      else if( pt < 75)return 0.956;
      else if( pt < 100)return 0.936;
      else if( pt < 200)return 0.898;
      else if( pt < 300)return 0.809;
      else return 0.718;
    }
    else   if(eta < 2.2){
      if( pt < 20)return 0.565;
      else if( pt < 25)return 0.902;
      else if( pt < 30)return 0.920;
      else if( pt < 35)return 0.928;
      else if( pt < 40)return 0.934;
      else if( pt < 50)return 0.936;
      else if( pt < 75)return 0.936;
      else if( pt < 100)return 0.928;
      else if( pt < 200)return 0.896;
      else if( pt < 300)return 0.786;
      else return 0.703;
    }
    else {
      if( pt < 20)return 0.470;
      else if( pt < 25)return 0.778;
      else if( pt < 30)return 0.817;
      else if( pt < 35)return 0.842;
      else if( pt < 40)return 0.858;
      else if( pt < 50)return 0.876;
      else if( pt < 75)return 0.886;
      else if( pt < 100)return 0.886;
      else if( pt < 200)return 0.872;
      else if( pt < 300)return 0.837;
      else return 0.700;
    }
  }
  
  if (trigname == "HLT_Mu8"){
    if(eta < 0.8){
      if( pt < 10)return 0.917;
      else if( pt < 20)return 0.927;
      else if( pt < 25)return 0.931;
      else if( pt < 30)return 0.932;
      else if( pt < 35)return 0.932;
      else if( pt < 40)return 0.933;
      else if( pt < 50)return 0.933;
      else if( pt < 75)return 0.929;
      else if( pt < 100)return 0.924;
      else if( pt < 200)return 0.920;
      else if( pt < 300)return 0.901;
      else  return 0.899;
    }
    else   if(eta < 1.25){
      if( pt < 10)return 0.914;
      else if( pt < 20)return 0.918;
      else if( pt < 25)return 0.921;
      else if( pt < 30)return 0.922;
      else if( pt < 35)return 0.920;
      else if( pt < 40)return 0.920;
      else if( pt < 50)return 0.920;
      else if( pt < 75)return 0.917;
      else if( pt < 100)return 0.913;
      else if( pt < 200)return 0.905;
      else if( pt < 300)return 0.876;
      else return 0.901;
    }
    else   if(eta < 1.6){
      if( pt < 10)return 0.908;
      else if( pt < 20)return 0.946;
      else if( pt < 25)return 0.957;
      else if( pt < 30)return 0.962;
      else if( pt < 35)return 0.964;
      else if( pt < 40)return 0.965;
      else if( pt < 50)return 0.964;
      else if( pt < 75)return 0.958;
      else if( pt < 100)return 0.938;
      else if( pt < 200)return 0.899;
      else if( pt < 300)return 0.812;
      else return 0.718;
    }
    else   if(eta < 2.2){
      if( pt < 10)return 0.869;
      else if( pt < 20)return 0.915;
      else if( pt < 25)return 0.932;
      else if( pt < 30)return 0.937;
      else if( pt < 35)return 0.941;
      else if( pt < 40)return 0.943;
      else if( pt < 50)return 0.944;
      else if( pt < 75)return 0.941;
      else if( pt < 100)return 0.932;
      else if( pt < 200)return 0.899;
      else if( pt < 300)return 0.791;
      else return 0.711;
    }
  }
  return 1.;
}

     

bool AnalyzerCore::PassTrigger(vector<pair<TString,TString> > list){
  
  vector<TString> triglist;
  for(vector<pair<TString,TString>  >::iterator it = list.begin(); it != list.end(); it++){
    if(k_channel.Contains(it->second)){
      for(unsigned int i=0; i <  triglist.size(); i++){
	if(PassTrigger(it->first)) return false;
      }
      if(PassTrigger(it->first)) return true;
      triglist.push_back(it->first);
    }
  }
  return false;
}



bool AnalyzerCore::PassTrigger(TString trig){
  vector<TString> list;
  list.push_back(trig);
  int  prescaler=1.;
  return TriggerSelector(list, eventbase->GetTrigger().GetHLTInsideDatasetTriggerNames(), eventbase->GetTrigger().GetHLTInsideDatasetTriggerDecisions(), eventbase->GetTrigger().GetHLTInsideDatasetTriggerPrescales(), prescaler);
 
 
}




TDirectory* AnalyzerCore::GetTemporaryDirectory(void) const
{
  gROOT->cd();
  TDirectory* tempDir = 0;
  int counter = 0;
  while (not tempDir) {
    // First, let's find a directory name that doesn't exist yet:                                              
    std::stringstream dirname;
    dirname << "HNCommonLeptonFakes_%i" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    // Let's try to make this directory:                                                                       
    tempDir = gROOT->mkdir((dirname.str()).c_str());

  }

  return tempDir;

}


void AnalyzerCore::Message(TString message, LQMsgType type){
  m_logger <<  type << message << LQLogger::endmsg;
}


void AnalyzerCore::MakeCleverHistograms(histtype type, TString clhistname ){
  
  //// ELECTRON PLOTs                                                                                          
  if(type==elhist) mapCLhistEl[clhistname] = new ElectronPlots(clhistname);
  //// BaseSelection::MUON PLOTs                                                                                              
  if(type==muhist) mapCLhistMu[clhistname] = new MuonPlots(clhistname);
  /// JET PLOTs                                                                                                
  if(type==jethist) mapCLhistJet[clhistname] = new JetPlots(clhistname);
  /// Signal plots                                                                                             
  if(type==sighist_ee)  mapCLhistSigEE[clhistname] = new SignalPlotsEE(clhistname);
  if(type==sighist_mm)  mapCLhistSigMM[clhistname] = new SignalPlotsMM(clhistname);
  if(type==sighist_em)  mapCLhistSigEM[clhistname] = new SignalPlotsEM(clhistname);

  if(type==trilephist)  mapCLhistTriLep[clhistname] = new TriLeptonPlots(clhistname);
      
  return;
}

void AnalyzerCore::MakeHistograms(){
  //// Additional plots to make                                                                                
  maphist.clear();
  maphist2D.clear();

    
}

void AnalyzerCore::MakeHistograms(TString hname, int nbins, float xbins[], TString label){
  maphist[hname] =  new TH1D(hname.Data(),hname.Data(),nbins,xbins);
  maphist[hname]->GetXaxis()->SetTitle(label);
}

void AnalyzerCore::MakeHistograms(TString hname, int nbins, float xmin, float xmax, TString label){

  maphist[hname] =  new TH1D(hname.Data(),hname.Data(),nbins,xmin,xmax);
  //maphist[hname]->GetXaxis()->SetTitle("TEST");
  maphist[hname]->GetXaxis()->SetTitle(label);
}

void AnalyzerCore::MakeHistograms2D(TString hname, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax, TString label) {

  maphist2D[hname] =  new TH2D(hname.Data(),hname.Data(),nbinsx,xmin,xmax, nbinsy,ymin,ymax);
  maphist2D[hname]->GetXaxis()->SetTitle(label);
}

void AnalyzerCore::MakeHistograms2D(TString hname, int nbinsx,  float xbins[], int nbinsy,  float ybins[], TString label) {

  maphist2D[hname] =  new TH2D(hname.Data(),hname.Data(),nbinsx , xbins, nbinsy,ybins);
  maphist2D[hname]->GetXaxis()->SetTitle(label);
}

bool AnalyzerCore::PassMETFilter(){
  
  bool pass (true);
  
  ///https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
    
  if(Is2015Analysis()){
    if (!eventbase->GetEvent().PassCSCHaloFilterTight()) {
      pass = false;
      m_logger << DEBUG << "Event Fails PassCSCHaloFilterTight " << LQLogger::endmsg;
    }

    if (!eventbase->GetEvent().PassHBHENoiseFilter()) {
      pass = false; 
      m_logger << DEBUG << "Event Fails PassHBHENoiseFilter " << LQLogger::endmsg;
    }

    if(!eventbase->GetEvent().PassEcalDeadCellTriggerPrimitiveFilter()) {
      pass = false;
      m_logger << DEBUG << "Event Fails PassEcalDeadCellTriggerPrimitiveFilter" << LQLogger::endmsg;
    }

    //Bad EE Supercrystal filter (post-ICHEP: extend to include an additional problematic SC --only for 2012)
    if (!eventbase->GetEvent().PassBadEESupercrystalFilter()) {
      pass = false;
      m_logger << DEBUG << "Event Fails PassBadEESupercrystalFilter" << LQLogger::endmsg;
    }
  }
  else{
    
    if (!eventbase->GetEvent().PassCSCHaloFilterTight()) {
      pass = false;
      m_logger << DEBUG << "Event Fails PassCSCHaloFilterTight " << LQLogger::endmsg;
    }
    
    if (!eventbase->GetEvent().PassTightHalo2016Filter()) {
      pass = false;
      m_logger << DEBUG << "Event Fails PassTightHalo2016Filter " << LQLogger::endmsg;
    }
    
    //if(isData){
    if (!eventbase->GetEvent().PassHBHENoiseFilter()) {
      pass = false; 
      m_logger << DEBUG << "Event Fails PassHBHENoiseFilter " << LQLogger::endmsg;
    }
    if (!eventbase->GetEvent().PassHBHENoiseIsoFilter()) {
      pass = false;
      m_logger << DEBUG << "Event Fails PassHBHENoiseIsoFilter " << LQLogger::endmsg;
    }
    if(!eventbase->GetEvent().PassEcalDeadCellTriggerPrimitiveFilter()) {
      pass = false;
      m_logger << DEBUG << "Event Fails PassEcalDeadCellTriggerPrimitiveFilter" << LQLogger::endmsg;
    }
  }
    //Bad EE Supercrystal filter (post-ICHEP: extend to include an additional problematic SC --only for 2012)
    //    if (!eventbase->GetEvent().PassBadEESupercrystalFilter()) {
    //      pass = false;
    //      m_logger << DEBUG << "Event Fails PassBadEESupercrystalFilter" << LQLogger::endmsg;
    //    }
    //}
  return pass;
}



void AnalyzerCore::FillHist(TString histname, float value, float w, float xbins[], int nbins , TString label){
  m_logger << DEBUG << "FillHist : " << histname << LQLogger::endmsg;

  if(GetHist(histname)) GetHist(histname)->Fill(value, w);
  
  else{
    if (nbins < 0) {
      m_logger << ERROR << histname << " was NOT found. Nbins was not set also... please configure histogram maker correctly" << LQLogger::endmsg;
      exit(0);
    }
    m_logger << DEBUG << "Making the histogram" << LQLogger::endmsg;
    MakeHistograms(histname, nbins, xbins, label);
    if(GetHist(histname)) GetHist(histname)->GetXaxis()->SetTitle(label);
    if(GetHist(histname)) GetHist(histname)->Fill(value, w);
  }

}


void AnalyzerCore::FillHistPerLumi(TString histname, float value, float w, float xmin, float xmax,int nbins, int nlumibins){
  
  if(!Is2015Analysis()){
    if(nlumibins==10){
      
      if(!GetHist(histname+"_perlumi")) {
	MakeHistograms(histname+"_perlumi", 9, 0., 9.);
	int nbin=0;

	for(std::map<int,TString>::iterator it = mapLumiNamePerBlock.begin(); it != mapLumiNamePerBlock.end(); it++){
	  nbin++;
	  GetHist(histname+"_perlumi")->GetXaxis()->SetBinLabel(nbin,it->second);
	}
      }
      
      for(map<int,TString>::iterator it = mapLumiNamePerBlock.begin(); it != mapLumiNamePerBlock.end(); it++){
	
        if(!GetHist(histname+"_"+it->second)) {
          MakeHistograms(histname+"_"+it->second, nbins, xmin, xmax);
	}
      } 


      for(map<int,TString>::iterator it = mapLumiNamePerBlock.begin(); it != mapLumiNamePerBlock.end(); it++){
	if(eventbase->GetEvent().RunNumber()  < it->first) {
	  map<int,float>::iterator it2 = mapLumiPerBlock.find(it->first);
	  
	  if(isData){
	    float neww= w /it2->second;
	    if(GetHist(histname+"_perlumi")) GetHist(histname+"_perlumi")->Fill(it->second, neww);
	    if(GetHist(histname+"_"+it->second)) GetHist(histname+"_"+it->second)->Fill(value,neww);
	  }
	  else{
	    //	    float neww = w * (it2->second/TargetLumi);
	    float neww = w /TargetLumi;
	    if(GetHist(histname+"_perlumi")) GetHist(histname+"_perlumi")->Fill(it->second, neww);
	    if(GetHist(histname+"_"+it->second)) GetHist(histname+"_"+it->second)->Fill(value,neww);
	  }
	} 
      }
    }// nbins
  }

}

void AnalyzerCore::FillHist(TString histname, float value, float w, float xmin, float xmax, int nbins , TString label){
  
  m_logger << DEBUG << "FillHist : " << histname << LQLogger::endmsg;
  if(GetHist(histname)) GetHist(histname)->Fill(value, w);  
  else{
    if (nbins < 0) {
      m_logger << ERROR << histname << " was NOT found. Nbins was not set also... please configure histogram maker correctly" << LQLogger::endmsg;
      exit(0);
    }
    m_logger << DEBUG << "Making the histogram" << LQLogger::endmsg;
    MakeHistograms(histname, nbins, xmin, xmax, label);
    if(GetHist(histname)) GetHist(histname)->GetXaxis()->SetTitle(label);
    if(GetHist(histname)) GetHist(histname)->Fill(value, w);
  }
  
}

void AnalyzerCore::FillHist(TString histname, float value1, float value2, float w, float xmin, float xmax, int nbinsx, float ymin, float ymax, int nbinsy , TString label){

  m_logger << DEBUG << "FillHist : " << histname << LQLogger::endmsg;
  if(GetHist2D(histname)) GetHist2D(histname)->Fill(value1,value2, w);
  else{
    if (nbinsx < 0) {
      m_logger << ERROR << histname << " was NOT found. Nbins was not set also... please configure histogram maker correctly" << LQLogger::endmsg;
      exit(0);
    }
    m_logger << DEBUG << "Making the histogram" << LQLogger::endmsg;
    MakeHistograms2D(histname, nbinsx, xmin, xmax,nbinsy, ymin, ymax , label);
    if(GetHist2D(histname)) GetHist2D(histname)->GetXaxis()->SetTitle(label);
    if(GetHist2D(histname)) GetHist2D(histname)->Fill(value1,value2, w);
  }

}

void AnalyzerCore::FillHist(TString histname, float valuex, float valuey, float w, float xbins[], int nxbins, float ybins[], int nybins , TString label){
  m_logger << DEBUG << "FillHist : " << histname << LQLogger::endmsg;
  if(GetHist2D(histname)) GetHist2D(histname)->Fill(valuex,valuey, w);

  else{
    if (nxbins < 0) {
      m_logger << ERROR << histname << " was NOT found. Nbins was not set also... please configure histogram maker correctly" << LQLogger::endmsg;
      exit(0);
    }
    m_logger << DEBUG << "Making the histogram" << LQLogger::endmsg;
    MakeHistograms2D(histname, nxbins, xbins, nybins, ybins , label);
    if(GetHist2D(histname)) GetHist2D(histname)->GetXaxis()->SetTitle(label);
    
    if(GetHist2D(histname)) GetHist2D(histname)->Fill(valuex, valuey, w);
  }

}



void AnalyzerCore::FillHist(TString histname, float value, float w , TString label){

  if(GetHist(histname)){
    GetHist(histname)->Fill(value, w);  /// Plots Z peak        
    GetHist(histname)->GetXaxis()->SetTitle(label);
  }
  
  else m_logger << INFO << histname << " was NOT found. Will add the histogram to the hist map on first event." << LQLogger::endmsg;
  
  
  return;
}

void AnalyzerCore::FillCLHist(histtype type, TString hist, vector<snu::KMuon> muons, double w){

  if(type==muhist){
    map<TString, MuonPlots*>::iterator mupit = mapCLhistMu.find(hist);
    if(mupit != mapCLhistMu.end()) mupit->second->Fill(w,muons);
    else m_logger << INFO  << hist << " not found in mapCLhistMu" << LQLogger::endmsg;
  }
  else  m_logger << INFO  << "Type not set to muhist, is this a mistake?" << LQLogger::endmsg;

}


void AnalyzerCore::FillCLHist(histtype type, TString hist, vector<snu::KElectron> electrons,double w){

  if(type==elhist){
    map<TString, ElectronPlots*>::iterator elpit = mapCLhistEl.find(hist);
    if(elpit !=mapCLhistEl.end()) elpit->second->Fill(w,electrons);
    else m_logger << INFO  << hist << " not found in mapCLhistEl" <<LQLogger::endmsg;
  }
  else  m_logger << INFO  << "Type not set to elhist, is this a mistake?" << LQLogger::endmsg;
}

void AnalyzerCore::FillCLHist(histtype type, TString hist, vector<snu::KJet> jets, double w){
  
  if(type==jethist){
    map<TString, JetPlots*>::iterator jetpit = mapCLhistJet.find(hist);
    if(jetpit !=mapCLhistJet.end()) jetpit->second->Fill(w,jets);
    else m_logger << INFO  << hist << " not found in mapCLhistJet" <<LQLogger::endmsg;
  }
  else  m_logger << INFO  <<"Type not set to jethist, is this a mistake?" << LQLogger::endmsg;

}

void AnalyzerCore::FillCLHist(histtype type, TString hist, snu::KEvent ev,vector<snu::KMuon> muons, vector<snu::KElectron> electrons, vector<snu::KJet> jets,double w){

  if(type==trilephist){

    map<TString, TriLeptonPlots*>::iterator trilepit = mapCLhistTriLep.find(hist);
    if(trilepit !=mapCLhistTriLep.end()) trilepit->second->Fill(ev, muons, electrons, jets,w);
    else {
      mapCLhistTriLep[hist] = new TriLeptonPlots(hist);
      trilepit = mapCLhistTriLep.find(hist);
      trilepit->second->Fill(ev, muons, electrons, jets,w);
    }
  }
  else if(type==sighist_ee){

    map<TString, SignalPlotsEE*>::iterator sigpit_ee = mapCLhistSigEE.find(hist);
    if(sigpit_ee !=mapCLhistSigEE.end()) sigpit_ee->second->Fill(ev, muons, electrons, jets,w);
    else {
      mapCLhistSigEE[hist] = new SignalPlotsEE(hist);
      sigpit_ee = mapCLhistSigEE.find(hist);
      sigpit_ee->second->Fill(ev, muons, electrons, jets,w);
    }
  }
  else if(type==sighist_mm){

    map<TString, SignalPlotsMM*>::iterator sigpit_mm = mapCLhistSigMM.find(hist);
    if(sigpit_mm !=mapCLhistSigMM.end()) sigpit_mm->second->Fill(ev, muons, electrons, jets,w);
    else {
      mapCLhistSigMM[hist] = new SignalPlotsMM(hist);
      sigpit_mm = mapCLhistSigMM.find(hist);
      sigpit_mm->second->Fill(ev, muons, electrons, jets,w);
    }
  }
  else if(type==sighist_em){

    map<TString, SignalPlotsEM*>::iterator sigpit_em = mapCLhistSigEM.find(hist);
    if(sigpit_em !=mapCLhistSigEM.end()) sigpit_em->second->Fill(ev, muons, electrons, jets,w);
    else {
      mapCLhistSigEM[hist] = new SignalPlotsEM(hist);
      sigpit_em = mapCLhistSigEM.find(hist);
      sigpit_em->second->Fill(ev, muons, electrons, jets,w);
    }
  }
 else  m_logger << INFO  <<"Type not set to sighist, is this a mistake?" << LQLogger::endmsg;


}


void AnalyzerCore::WriteHistograms() throw (LQError){
  // This function is called after the cycle is ran. It wrues all histograms to the output file. This function is not used by user. But by the contrioller code.
  WriteHists();
  WriteCLHists();
}

  
void AnalyzerCore::WriteCLHists(){

  for(map<TString, MuonPlots*>::iterator mupit = mapCLhistMu.begin(); mupit != mapCLhistMu.end(); mupit++){

    Dir = m_outputFile->mkdir(mupit->first);
    m_outputFile->cd( Dir->GetName() );
    mupit->second->Write();
    m_outputFile->cd();
  }

  for(map<TString, ElectronPlots*>::iterator elpit = mapCLhistEl.begin(); elpit != mapCLhistEl.end(); elpit++)\
    {

      Dir = m_outputFile->mkdir(elpit->first);
      m_outputFile->cd( Dir->GetName() );
      elpit->second->Write();
      m_outputFile->cd();
    }

  for(map<TString, JetPlots*>::iterator jetpit = mapCLhistJet.begin(); jetpit != mapCLhistJet.end(); jetpit++)\
    {
      
      Dir = m_outputFile->mkdir(jetpit->first);
      m_outputFile->cd( Dir->GetName() );
      jetpit->second->Write();
      m_outputFile->cd();
    }
  for(map<TString, SignalPlotsEE*>::iterator sigpit_ee = mapCLhistSigEE.begin(); sigpit_ee != mapCLhistSigEE.end(); sigpit_ee++){
    
    Dir = m_outputFile->mkdir(sigpit_ee->first);
    m_outputFile->cd( Dir->GetName() );
    sigpit_ee->second->Write();
    m_outputFile->cd();
  }
  for(map<TString, SignalPlotsMM*>::iterator sigpit_mm = mapCLhistSigMM.begin(); sigpit_mm != mapCLhistSigMM.end(); sigpit_mm++){

    Dir = m_outputFile->mkdir(sigpit_mm->first);
    m_outputFile->cd( Dir->GetName() );
    sigpit_mm->second->Write();
    m_outputFile->cd();
  }
  for(map<TString, SignalPlotsEM*>::iterator sigpit_em = mapCLhistSigEM.begin(); sigpit_em != mapCLhistSigEM.end(); sigpit_em++){

    Dir = m_outputFile->mkdir(sigpit_em->first);
    m_outputFile->cd( Dir->GetName() );
    sigpit_em->second->Write();
    m_outputFile->cd();
  }
  
  for(map<TString, TriLeptonPlots*>::iterator trilepit = mapCLhistTriLep.begin(); trilepit != mapCLhistTriLep.end(); trilepit++){

    Dir = m_outputFile->mkdir(trilepit->first);
    m_outputFile->cd( Dir->GetName() );
    trilepit->second->Write();
    m_outputFile->cd();
  }



  return;
}

void AnalyzerCore::WriteHists(){

  /// Open Output rootfile
  m_outputFile->cd();

  for(map<TString, TH1*>::iterator mapit = maphist.begin(); mapit != maphist.end(); mapit++){
    
    
    
    if(mapit->first.Contains("closejet")){
      if(!m_outputFile->GetDirectory( "closejet" )){
	Dir = m_outputFile->mkdir("closejet");
	m_outputFile->cd( Dir->GetName() );
      }
      else  m_outputFile->cd("closejet");
      mapit->second->Write();
      m_outputFile->cd();
    }

    
    
    else {
      mapit->second->Write();
    }
  }
  
  for(map<TString, TH2*>::iterator mapit = maphist2D.begin(); mapit != maphist2D.end(); mapit++){
    mapit->second->Write();
  }


  return;
}

TH1* AnalyzerCore::GetHist(TString hname){

  TH1* h = NULL;
  std::map<TString, TH1*>::iterator mapit = maphist.find(hname);
  if(mapit != maphist.end()) return mapit->second;
  else m_logger << DEBUG  << hname << " was not found in map" << LQLogger::endmsg;

  return h;
}



TH2* AnalyzerCore::GetHist2D(TString hname){

  TH2* h = NULL;
  std::map<TString, TH2*>::iterator mapit = maphist2D.find(hname);
  if(mapit != maphist2D.end()) return mapit->second;
  else m_logger << DEBUG  << hname << " was not found in map" << LQLogger::endmsg;

  return h;
}




void AnalyzerCore::FillCutFlow(TString cut, float weight){

  bool cut_exists=false;
  m_logger << DEBUG  << "FillCutFlow" <<  LQLogger::endmsg;
  for(unsigned int i=0; i < cutflow_list.size();  i++){
    if (cut == cutflow_list.at(i)) cut_exists=true;
  }

  if (!cut_exists){
    m_logger << DEBUG  << "FillCutFlow: " << cut <<  LQLogger::endmsg;

    n_cutflowcuts=n_cutflowcuts+1;
    if(GetHist("cutflow")) {
      m_logger << DEBUG  << "FillCutFlow: mod " << cut <<  LQLogger::endmsg;

      map<TString, TH1*>::iterator it = maphist.find("cutflow");
      vector<float> counters;
      for(int j = 0; j  < n_cutflowcuts ; j++){
	counters.push_back(it->second->GetBinContent(1+j));
      }
      delete it->second;
      
      AnalyzerCore::MakeHistograms("cutflow", n_cutflowcuts,0.,float(n_cutflowcuts));
      for(unsigned int i=0;i < cutflow_list.size();  i++){
	GetHist("cutflow")->GetXaxis()->SetBinLabel(i+1,cutflow_list.at(i));
	GetHist("cutflow")->Fill(cutflow_list.at(i), weight);
      }
      GetHist("cutflow")->GetXaxis()->SetBinLabel(n_cutflowcuts,cut);
      cutflow_list.push_back(cut);
      
    }
    else{
      m_logger << DEBUG  << "FillCutFlow:Fill " << cut <<  LQLogger::endmsg;
      AnalyzerCore::MakeHistograms("cutflow", n_cutflowcuts,0.,float(n_cutflowcuts));
      GetHist("cutflow")->GetXaxis()->SetBinLabel(1,cut);
      GetHist("cutflow")->Fill(cut,weight);
      cutflow_list.push_back(cut);
    }
  }
  else {
    if(GetHist("cutflow")) {
      GetHist("cutflow")->Fill(cut,weight);
    }
  }
}


bool AnalyzerCore::Zcandidate(std::vector<snu::KMuon> muons, float interval, bool require_os){
  
  if(muons.size()!=2) return false;
  if(require_os&&SameCharge(muons)) return false;
  
  KParticle Z = muons.at(0) + muons.at(1);
  if(fabs(Z.M() - 90.) <  interval) return true;
  else return false;
  
}
  
bool AnalyzerCore::SameCharge(std::vector<snu::KMuon> muons){
  
  if(muons.size()!=2) return false;
  if(muons.at(0).Charge() == muons.at(1).Charge()) return true;
  return false;
}


bool AnalyzerCore::Zcandidate(std::vector<snu::KElectron> electrons, float interval, bool require_os){

  if(electrons.size()!=2) return false;
  if(require_os&&SameCharge(electrons)) return false;

  KParticle Z = electrons.at(0) + electrons.at(1);
  if(fabs(Z.M() - 90.) <  interval) return true;
  else return false;

}

bool AnalyzerCore::SameCharge(std::vector<snu::KElectron> electrons, bool runningcf){
  
  if(electrons.size() > 2){
    int p_charge=0;
    int m_charge=0;
    for(unsigned int iel = 0 ; iel < electrons.size() ; iel++){
      if(electrons.at(iel).Charge() < 0 ) m_charge++;
      if(electrons.at(iel).Charge() > 0 ) p_charge++;
    }
    if(p_charge > 1) return true;
    if(m_charge > 1) return true;
  }
  if(electrons.size()!=2) return false;


  if(!runningcf){
    if(electrons.at(0).Charge() == electrons.at(1).Charge()) return true;
  }
  else     if(electrons.at(0).Charge() != electrons.at(1).Charge()) return true;

  return false;
}


int AnalyzerCore::NBJet(std::vector<snu::KJet> jets,  KJet::Tagger tag, KJet::WORKING_POINT wp){

  int nbjet=0;

  TString btag_key_lf("") , btag_key_hf("");
  TString wp_string="";
  if(wp == snu::KJet::Loose)wp_string = "Loose";
  if(wp == snu::KJet::Medium)wp_string = "Medium";
  if(wp == snu::KJet::Tight)wp_string = "Tight";

  TString tag_string="";
  if(tag== snu::KJet::CSVv2) tag_string ="CSVv2";
  if(tag== snu::KJet::cMVAv2) tag_string = "cMVAv2"; 
   
  btag_key_lf = tag_string+"_"+wp_string+"_lf";
  btag_key_hf = tag_string+"_"+wp_string+"_hf";
  std::map<TString,BTagSFUtil*>::iterator it_lf = MapBTagSF.find(btag_key_lf);
  std::map<TString,BTagSFUtil*>::iterator it_hf = MapBTagSF.find(btag_key_hf);
  
  if(it_lf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}
  if(it_hf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}

  /// systematics allowed are +-1 and +-3 for HN analysis 
  if ( tag == snu::KJet::JETPROB) return -999;
  for(unsigned int ij=0; ij <jets.size(); ij++){

    if(IsBTagged(jets.at(ij), tag, wp)) nbjet++;
    continue;

    bool isBtag=false;
    if (isData) {

      if (it_lf->second->IsTagged(jets.at(ij).BJetTaggerValue(tag),  -999999, jets.at(ij).Pt(), jets.at(ij).Eta()))
	isBtag=true;
    }
    else if (jets.at(ij).HadronFlavour() > 1){
      if (it_hf->second->IsTagged(jets.at(ij).BJetTaggerValue(tag),  jets.at(ij).HadronFlavour(),jets.at(ij).Pt(), jets.at(ij).Eta()))
        isBtag=true;
    }
    else{
      if (it_lf->second->IsTagged(jets.at(ij).BJetTaggerValue(tag),  jets.at(ij).HadronFlavour(),jets.at(ij).Pt(), jets.at(ij).Eta()))
	isBtag=true;
    }
    
    if(isBtag )nbjet++;
  }
  return nbjet;
}


int AnalyzerCore::IsBTagged(snu::KJet jet,  KJet::Tagger tag, KJet::WORKING_POINT wp){

  int nbjet=0;

  TString btag_key_lf("") , btag_key_hf("");
  TString wp_string="";
  if(wp == snu::KJet::Loose)wp_string = "Loose";
  if(wp == snu::KJet::Medium)wp_string = "Medium";
  if(wp == snu::KJet::Tight)wp_string = "Tight";

  TString tag_string="";
  if(tag== snu::KJet::CSVv2) tag_string ="CSVv2";
  if(tag== snu::KJet::cMVAv2) tag_string = "cMVAv2";

  btag_key_lf = tag_string+"_"+wp_string+"_lf";
  btag_key_hf = tag_string+"_"+wp_string+"_hf";
  std::map<TString,BTagSFUtil*>::iterator it_lf = MapBTagSF.find(btag_key_lf);
  std::map<TString,BTagSFUtil*>::iterator it_hf = MapBTagSF.find(btag_key_hf);

  if(it_lf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}
  if(it_hf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}

  /// systematics allowed are +-1 and +-3 for HN analysis
  if ( tag == snu::KJet::JETPROB) return -999;
  
  bool isBtag=false;
  if (isData) {
    
    if (it_lf->second->IsTagged(jet.BJetTaggerValue(tag),  -999999, jet.Pt(), jet.Eta()))
      isBtag=true;
  }
    else if (jet.HadronFlavour() > 1){
      if (it_hf->second->IsTagged(jet.BJetTaggerValue(tag),  jet.HadronFlavour(),jet.Pt(), jet.Eta()))
        isBtag=true;
    }
    else{
      if (it_lf->second->IsTagged(jet.BJetTaggerValue(tag),  jet.HadronFlavour(),jet.Pt(), jet.Eta()))
        isBtag=true;
    }
  
  if(isBtag )nbjet++;
  
  return nbjet;
}

double AnalyzerCore::MuonDYMassCorrection(std::vector<snu::KMuon> mu, double w){
  
  if(mu.size()< 2) return 0.;
  snu::KParticle Z = mu.at(0) + mu.at(1);
  
  double factor (1.);
  if(Z.M() > 90.){
    factor = 8.37401e-01 + 1.61277e-03*Z.M();
  }
  return w*factor;
}

float AnalyzerCore::CFRate(snu::KElectron el){
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

bool AnalyzerCore::IsTight(snu::KMuon muon){
  /// ADD TIGHT BaseSelection::MUON REQUIREMENT
  float reliso= muon.RelIso04();

  if(( reliso >= 0.1)) return false;
  if(fabs(muon.dXY()) >= 0.05) return false; 
  return true;
}


bool AnalyzerCore::IsTight(snu::KElectron el){
  
  return eventbase->GetElectronSel()->PassUserID("ELECTRON_HN_TIGHT",el);

}
  
bool AnalyzerCore::IsCF(snu::KElectron el){
  vector<snu::KTruth> truth =  eventbase->GetTruth();
  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){
    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;
    if(fabs(eventbase->GetTruth().at(ig).PdgId()) == 11){
      if(fabs(eventbase->GetTruth().at(eventbase->GetTruth().at(ig).IndexMother()).PdgId()) == 23 ||
	 fabs(eventbase->GetTruth().at(eventbase->GetTruth().at(ig).IndexMother()).PdgId()) == 24){
	if(eventbase->GetTruth().at(ig).PdgId() * el.Charge() > 0 ) return true;
	else return false;
      }
    }
  }
  return false;
}

vector<snu::KElectron> AnalyzerCore::GetTruePrompt(vector<snu::KElectron> electrons, bool keep_chargeflip, bool keepfake){
  if(electrons.size() == 0)
    return electrons;
  
  vector<snu::KElectron> prompt_electrons;
  for(unsigned int i = 0; i < electrons.size(); i++){

    if(!k_isdata){
      if(keepfake&&keep_chargeflip) prompt_electrons.push_back(electrons.at(i));
      else if(keep_chargeflip&&electrons.at(i).MCMatched()) prompt_electrons.push_back(electrons.at(i));
      else if(keepfake&&! electrons.at(i).MCIsCF()) prompt_electrons.push_back(electrons.at(i)); 
      else if(electrons.at(i).MCMatched() && !electrons.at(i).MCIsCF()) prompt_electrons.push_back(electrons.at(i));
    }// Data
    else prompt_electrons.push_back(electrons.at(i));
  }/// loop

  return prompt_electrons;


}

vector<snu::KMuon> AnalyzerCore::GetTruePrompt(vector<snu::KMuon> muons, bool keepfake){
  if(muons.size()==0)return muons;

  vector<snu::KMuon> prompt_muons;

  for(unsigned int i = 0; i < muons.size(); i++){
    if(!k_isdata){

      if(keepfake) prompt_muons.push_back(muons.at(i));
      else if(muons.at(i).MCMatched()) prompt_muons.push_back(muons.at(i));
    }// Data
    else prompt_muons.push_back(muons.at(i));
  }/// loop
  return prompt_muons;

}


void AnalyzerCore::CorrectMuonMomentum(vector<snu::KMuon>& k_muons){
  
  if(Is2015Analysis()){
    vector<TLorentzVector> tlv_muons = MakeTLorentz(k_muons);
    int imu(0);
    for(std::vector<snu::KMuon>::iterator it = k_muons.begin(); it != k_muons.end(); it++, imu++){
      float qter =1.; /// uncertainty

      if(k_isdata)rmcor->momcor_data(tlv_muons[imu], float(it->Charge()), eventbase->GetEvent().RunNumber(), qter);
      else rmcor->momcor_mc(tlv_muons[imu], float(it->Charge()), it->ActiveLayer(), qter);

      it->SetPtEtaPhiM(tlv_muons[imu].Pt(),tlv_muons[imu].Eta(), tlv_muons[imu].Phi(), tlv_muons[imu].M());
      //it->SetPtEtaPhiM(tlv_muons[imu].Pt(), it->Eta(), it->Phi() , it->M());
      //it->scale(tlv_muons[imu].E()/it->E());
    }
    
  }
  else{
    Message("In CorrectMuonMomentum", DEBUG);
    
    vector<TLorentzVector> tlv_muons = MakeTLorentz(k_muons);
    int imu(0);
    for(std::vector<snu::KMuon>::iterator it = k_muons.begin(); it != k_muons.end(); it++, imu++){
      if(Is2015Analysis()) continue;
      it->SetPtEtaPhiE(it->RochPt(), it->RochEta(), it->RochPhi(), it->RochE()); 
    }
  }
  Message("END CorrectMuonMomentum",DEBUG);

}



float AnalyzerCore::Get_DataDrivenWeight_EM(vector<snu::KMuon> k_muons, vector<snu::KElectron> k_electrons, TString cut){

  Message("In Get_DataDrivenWeight_EM", DEBUG);
  
  if(k_muons.size()!=1) return 0.;
  if(k_electrons.size()!=1) return 0.;

  bool is_mu1_tight    = IsTight(k_muons.at(0));
  bool is_el1_tight    = IsTight(k_electrons.at(0));

  vector<TLorentzVector> muons=MakeTLorentz(k_muons);
  vector<TLorentzVector> electrons=MakeTLorentz(k_electrons);

  float em_weight = m_fakeobj->get_dilepton_em_eventweight(muons, electrons, is_mu1_tight, is_el1_tight, cut, 0);

  return em_weight;
}

float AnalyzerCore::Get_DataDrivenWeight_MM(vector<snu::KMuon> k_muons, TString cutID){

  Message("In Get_DataDrivenWeight_MM", DEBUG);
  if(k_muons.size()==0) return 0.;
  float mm_weight = 0.;

  if(k_muons.size()==2){

    bool is_mu1_tight    = IsTight(k_muons.at(0));
    bool is_mu2_tight    = IsTight(k_muons.at(1));
    
    vector<TLorentzVector> muons=MakeTLorentz(k_muons);

    mm_weight =m_fakeobj->get_dilepton_mm_eventweight(muons, is_mu1_tight,is_mu2_tight, cutID);
    
  }
  return mm_weight;
}


float AnalyzerCore::Get_DataDrivenWeight_M(vector<snu::KMuon> k_muons, TString cutID){

  Message("In Get_DataDrivenWeight_M", DEBUG);

  if(k_muons.size()!=1) return 0.;

  bool is_mu1_tight    = IsTight(k_muons.at(0));
  vector<TLorentzVector> muons=MakeTLorentz(k_muons);

  float r = 1.;
  if( k_muons.at(0).Pt() < 20.) r = 0.82;
  else if( k_muons.at(0).Pt() < 30.) r = 0.85;
  else if( k_muons.at(0).Pt() < 40.) r = 0.96;
  else if( k_muons.at(0).Pt() < 50.) r = 0.92;
  else r = 0.95;

  TString id="fake_Eff_muon_pog";
  if(cutID.Contains("HN")) id= "fake_Eff_muon_hn";
  float f=  m_fakeobj->getFakeRate_muon(0,fabs(k_muons.at(0).Eta()), k_muons.at(0).Pt(), id);

  float w = m_fakeobj->lepton_weight(!is_mu1_tight, r,f);
  return w;

}

float AnalyzerCore::Get_DataDrivenWeight_E(vector<snu::KElectron> k_electrons){
  Message("In Get_DataDrivenWeight_EE", DEBUG);
  return 0.;
  if(k_electrons.size()!=1) return 0.;

  bool is_el1_tight    = IsTight(k_electrons.at(0));
  vector<TLorentzVector> electrons=MakeTLorentz(k_electrons);

  float r = 1.;
  if( k_electrons.at(0).Pt() < 20.) r = 0.82;
  else if( k_electrons.at(0).Pt() < 30.) r = 0.85;
  else if( k_electrons.at(0).Pt() < 40.) r = 0.96;
  else if( k_electrons.at(0).Pt() < 50.) r = 0.92;
  else r = 0.95;

  TString cut1  = "pt_eta_40_looseregion1";
  
  float f=  m_fakeobj->getFakeRate_electronEta(0,k_electrons.at(0).Pt(),fabs(k_electrons.at(0).Eta()),cut1);
  
  float w = m_fakeobj->lepton_weight(!is_el1_tight, r,f);
  return w;

}



float AnalyzerCore::Get_DataDrivenWeight_EE(vector<snu::KElectron> k_electrons){

  Message("In Get_DataDrivenWeight_EE", DEBUG);

  if(k_electrons.size()==0) return 0.;

  float ee_weight = 0.;
  if(k_electrons.size()==2){
    
    

    bool is_el1_tight    = IsTight(k_electrons.at(0));
    bool is_el2_tight    = IsTight(k_electrons.at(1));

    vector<TLorentzVector> electrons=MakeTLorentz(k_electrons);
    ee_weight =m_fakeobj->get_dilepton_ee_eventweight(electrons, is_el1_tight,is_el2_tight); 
    
  }
  return ee_weight;

}



void AnalyzerCore::MakeNtp(TString hname, TString myvar){

  mapntp[hname] =  new TNtupleD(hname.Data(),hname.Data(),myvar.Data());
}


void AnalyzerCore::FillNtp(TString hname, Double_t myinput[]){

  if (GetNtp(hname)) GetNtp(hname)->Fill(myinput);
  else m_logger << INFO << hname << " was NOT found. Check you ntp. " << LQLogger::endmsg;

  return;
}

// //
void AnalyzerCore::WriteNtp(){

  /// Open Output rootfile
  m_outputFile->cd();

  for(map<TString, TNtupleD*>::iterator mapit = mapntp.begin(); mapit != mapntp.end(); mapit++){
    mapit->second->Write();
  }

  return;
}
// //



// //
TNtupleD* AnalyzerCore::GetNtp(TString hname){

  TNtupleD* n = NULL;
  std::map<TString, TNtupleD*>::iterator mapit = mapntp.find(hname);
  if (mapit != mapntp.end()) return mapit->second;
  else m_logger << INFO << hname << " was not found in map" << LQLogger::endmsg;

  return n;
}


vector<TLorentzVector> AnalyzerCore::MakeTLorentz(vector<snu::KElectron> el){

  vector<TLorentzVector> tl_el;
  for(vector<KElectron>::iterator itel=el.begin(); itel!=el.end(); ++itel) {
    TLorentzVector tmp_em;
    tmp_em.SetPtEtaPhiM((*itel).Pt(),(*itel).Eta(),(*itel).Phi(),(*itel).M());
    tl_el.push_back(tmp_em);
  }
  return tl_el;
}

vector<TLorentzVector> AnalyzerCore::MakeTLorentz(vector<snu::KMuon> mu){
  
  vector<TLorentzVector> tl_mu;
  for(vector<KMuon>::iterator itmu=mu.begin(); itmu!=mu.end(); ++itmu) {
    TLorentzVector tmp_mu;
    tmp_mu.SetPtEtaPhiM((*itmu).Pt(),(*itmu).Eta(),(*itmu).Phi(),(*itmu).M());
    tl_mu.push_back(tmp_mu);
  }
  return tl_mu;
}


vector<TLorentzVector> AnalyzerCore::MakeTLorentz(vector<snu::KJet> j){

  vector<TLorentzVector> tl_jet;
  for(vector<KJet>::iterator itj=j.begin(); itj!=j.end(); ++itj) {
    TLorentzVector tmp_j;
    tmp_j.SetPtEtaPhiM((*itj).Pt(),(*itj).Eta(),(*itj).Phi(),(*itj).M());
    tl_jet.push_back(tmp_j);
  }
  return tl_jet;
}

bool AnalyzerCore::Compare_Jet_Pt(const snu::KJet& jet0, const snu::KJet& jet1)
{
  return jet0.Pt() > jet1.Pt();
}



