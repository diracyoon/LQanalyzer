#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include <iomanip>
#include <sstream>
#include "TSystem.h"

#include <map>
#include "SampleMap.C"
float GetEventsProcessed(std::string file);
float GetEventsPassed(std::string file);

// defined in boost/utility.hpp, by the way
template <typename Iter>
Iter next(Iter iter)
{
  return ++iter;
}

void MakeInputListForSubmitScript(){
    
  map<TString, TString> lqmap_tmp = GetLQMap();
  map<TString, TString> lqmap;
  for(std::map<TString, TString>::iterator mit =lqmap_tmp.begin(); mit != lqmap_tmp.end();++mit){
    TString def_version = TString(getenv("CATVERSION"));
    
    TString dir = "ls /data2/DATA/cattoflat/MC/" + def_version + "/"+ mit->first + "/*.root > inputlist_lists.txt";
    system(dir.Data());
    std::ifstream fin("inputlist_lists.txt");
    std::string word;
    float number_events_processed(0.);
    float number_events_passed(0.);
    while ( fin >> word ) {
      number_events_processed+= GetEventsProcessed(word);
      number_events_passed+= GetEventsPassed(word);
    }
    if(number_events_processed !=0 ) lqmap[mit->first] = mit->second;
    system("rm inputlist_lists.txt");
  }

    
    
  ofstream lumi_file;
  string lfile =  "list_all_mc.sh";
  
  lumi_file.open(lfile.c_str());
  lumi_file.setf(ios::fixed,ios::floatfield); 
  lumi_file.precision(1);
  lumi_file << "#!/bin/sh" << endl;
  lumi_file << "" << endl;
  lumi_file << "########################" << endl;
  lumi_file << "### SAMPLE LIST ########## "  << endl;
  lumi_file << "#######################" << endl;
  lumi_file << "" << endl;
  lumi_file << "" << endl;

  lumi_file << "source ${LQANALYZER_DIR}/LQRun/txt/list_user_mc.sh" << endl;
  lumi_file << "" << endl;
  lumi_file << "" << endl;
  
  lumi_file << " declare -a input_samples=('WZ_pythia8')" << endl;;
  lumi_file << "" << endl;
  lumi_file << " declare -a all_mc=('" ;
  

  for(std::map<TString, TString>::iterator mit =lqmap.begin(); mit != lqmap.end();++mit){
    
    cout << mit->first << " " << mit->second << endl;
    if(next(mit)!= lqmap.end()) lumi_file << mit->second << "' '"   ;
    else  lumi_file << mit->second ;
  }
  lumi_file << "') " << endl; 
  
  
  lumi_file << "" << endl;
  lumi_file << " declare -a mc_noqcd=('" ;
  for(std::map<TString, TString>::iterator mit =lqmap.begin(); mit != lqmap.end();++mit){
    if(mit->second.Contains("QCD")) continue;
    cout << mit->first << " " << mit->second << endl;
    if(next(mit)!= lqmap.end()) lumi_file << mit->second << "' '"   ;
    else  lumi_file << mit->second ;
    
  }
    lumi_file << "') " << endl;


    
  lumi_file << "" << endl;
  lumi_file << "" << endl;
  lumi_file << "declare -a diboson_pythia=('WZ_pythia8' 'ZZ_pythia8' 'WW_pythia8')" << endl;
  lumi_file << "" << endl;
  lumi_file << "declare -a dy_mcatnlo=('DY10to50_MCatNLO' 'DY50plus_MCatNLO') " << endl;
  lumi_file << "" << endl;
  lumi_file << "declare -a dilepton_list=('DY10to50_MCatNLO' 'DY50plus_MCatNLO' 'WJets_MCatNLO' 'WZ_pythia8' 'ZZ_pythia8' 'WW_pythia8'  'TT_MG5' 'singletop_s_MCatNLO' 'singletop_tbar_Powheg' 'singletop_t_Powheg' 'singletop_tbarW_Powheg' 'singletop_tW_Powheg' 'ttWJetsToLNu_MCatNLO' 'ttWJetsToQQ_MCatNLO' 'ttZToQQ_MCatNLO')" << endl;
  
  lumi_file << "declare -a trilepton_list=('DY10to50_MCatNLO' 'DY50plus_MCatNLO' 'WJets_MCatNLO' 'WZ_pythia8' 'ZZ_pythia8' 'WW_pythia8'  'TT_MG5' 'singletop_s_MCatNLO' 'singletop_tbar_Powheg' 'singletop_t_Powheg' 'singletop_tbarW_Powheg' 'singletop_tW_Powheg' 'ZZ_llnunu_powheg' 'ZZ_llqq_MCatNLO' 'ZZ_llll_MCatNLO' 'ZZ_llll_powheg' 'WJets_MCatNLO' 'WW_llnn_powheg' 'WZ_lllnu_powheg' 'WZ_lllnu_MCatNLO' 'WZ_llqq_MCatNLO' 'WN_lllnu_powheg' 'WZZ_MCatNLO' )" << endl;

  lumi_file << "" << endl;


  lumi_file << "declare -a qcd=('" ;
  for(std::map<TString, TString>::iterator mit =lqmap.begin(); mit != lqmap.end();++mit){
    if(!mit->second.Contains("QCD")) continue;
    bool last=true;
    for(std::map<TString, TString>::iterator mit2 =next(mit); mit2 != lqmap.end();++mit2){
      if(mit2->second.Contains("QCD")) last=false;
    }
    if(!last)lumi_file << mit->second << "' '"   ;
    else  lumi_file << mit->second ;
  } lumi_file << "') " << endl;
  lumi_file << "" << endl;
  
  lumi_file << "declare -a qcd_mu=('" ;
  for(std::map<TString, TString>::iterator mit =lqmap.begin(); mit != lqmap.end();++mit){
    if(!mit->second.Contains("QCD_mu")) continue;
    bool last=true;
    for(std::map<TString, TString>::iterator mit2 =next(mit); mit2 != lqmap.end();++mit2){
      if(mit2->second.Contains("QCD_mu")) last=false;
    }   
    if(!last)lumi_file << mit->second << "' '"   ;
    else  lumi_file << mit->second ;
  } lumi_file << "') " << endl;
  lumi_file << "" << endl;
  
  lumi_file << "declare -a qcd_eg=('" ;
  for(std::map<TString, TString>::iterator mit =lqmap.begin(); mit != lqmap.end();++mit){
    if(!mit->second.Contains("QCD")) continue;
    if(mit->second.Contains("QCD_mu")) continue;
    bool last=true;
    for(std::map<TString, TString>::iterator mit2 =next(mit); mit2 != lqmap.end();++mit2){
      if(mit2->second.Contains("QCD")) {
	if(!mit->second.Contains("QCD_mu")) {
	  last=false;
	}
      }
    }
    if(!last)lumi_file << mit->second << "' '"   ;
    else  lumi_file << mit->second ;
  } lumi_file << "') " << endl;
  lumi_file << "" << endl;


  lumi_file << "" << endl;
  lumi_file << "declare -a hn_mm=('WZ_pythia8' 'ZZ_pythia8' 'WpWp_madgraph' 'WpWp_qcd_madgraph'  'ttWJetsToLNu_MCatNLO' 'ttWJetsToQQ_MCatNLO' 'ttZToLLNuNu_MCatNLO' 'ttZToQQ_MCatNLO') " << endl;
  lumi_file << "" << endl;
  lumi_file << "declare -a hn_ee=('WZ_pythia8' 'ZZ_pythia8' 'WpWp_madgraph' 'WpWp_qcd_madgraph'  'ttWJetsToLNu_MCatNLO' 'ttWJetsToQQ_MCatNLO' 'ttZToLLNuNu_MCatNLO' 'ttZToQQ_MCatNLO' 'DY10to50_MCatNLO' 'DY50plus_MCatNLO') " << endl;
  lumi_file << "" << endl;
  lumi_file << "declare -a hn_fakeee=('DY10to50_MCatNLO' 'DY50plus_MCatNLO' 'WJets_MCatNLO'  'TT_MG5')" << endl;
  lumi_file << "" << endl;
  lumi_file << "" << endl;
  lumi_file << "declare -a singletop=('singletop_s_MCatNLO' 'singletop_tbar_Powheg' 'singletop_t_Powheg' 'singletop_tbarW_Powheg' 'singletop_tW_Powheg') " << endl;


  string lqdir = getenv("LQANALYZER_DIR");
  TString user = TString(getenv("USER"));


  string lfile2 =  lqdir + "/LQRun/txt/list_all_mc.sh";

  if(user.Contains("jalmond"))
    gSystem->Exec(("cp " + lfile + "  /data1/LQAnalyzer_rootfiles_for_analysis/CATAnalysis/").c_str());
  gSystem->Exec(("mv " + lfile +" " + lfile2).c_str());


  return;
  
}

float GetEventsProcessed(std::string filename){
  TFile* file = TFile::Open(filename.c_str());
  //  cout << file << endl;
  TH1F*  EventCounter = (TH1F*) (file ->Get("ntuple/hNEvent"));

  float value = EventCounter->GetBinContent(1);
  file->Close();
  return value;
}

float GetEventsPassed(std::string filename){
  TFile* file = TFile::Open(filename.c_str());

  TH1F*  EventCounter = (TH1F*) (file ->Get("ntuple/hNEvent"));

  float value = EventCounter->GetBinContent(1);
  file->Close();
  return value;
}