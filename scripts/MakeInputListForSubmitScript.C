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
#include "Luminosity/SampleMap.C"
float GetEventsProcessed(std::string file);
float GetEventsPassed(std::string file);

// defined in boost/utility.hpp, by the way
template <typename Iter>
Iter next(Iter iter)
{
  return ++iter;
}

void MakeInputListForSubmitScript(){
    
  TString path_of_list=TString(getenv("LQANALYZER_DATASET_DIR"))+ "/cattuplist_"+TString(getenv("CATVERSION"))+".txt";  

  map<TString, TString> lqmap_tmp = GetLQMap2016(path_of_list);
  map<TString, TString> lqmap;
  for(std::map<TString, TString>::iterator mit =lqmap_tmp.begin(); mit != lqmap_tmp.end();++mit){
    TString def_version = TString(getenv("CATVERSION"));
    
    TString analysisdir = TString(getenv("HOSTNAME"));
    bool cluster = false;
    if(analysisdir.Contains("cmscluster.snu.ac.kr")) cluster=true;

    TString dir = "ls /data2/DATA/cattoflat/MC/" + def_version + "/"+ mit->first + "/*.root > inputlist_lists.txt";
    if(cluster) dir = "ls  /data4/DATA/FlatCatuples/MC/" + def_version + "/"+ mit->first + "/*.root > inputlist_lists.txt";
    
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

  
  TString def_version = TString(getenv("CATVERSION"));
  
  ofstream lumi_file;
  TString lfile =  "list_all_mc_" + def_version+".sh";
  
  lumi_file.open(lfile.Data());
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
  
  lumi_file << "declare -a input_samples=('WZ')" << endl;;
  lumi_file << "" << endl;
  lumi_file << "declare -a all_mc=('" ;
  

  for(std::map<TString, TString>::iterator mit =lqmap.begin(); mit != lqmap.end();++mit){
    
    cout << mit->first << " " << mit->second << endl;
    if(next(mit)!= lqmap.end()) lumi_file << mit->second << "' '"   ;
    else  lumi_file << mit->second ;
  }
  lumi_file << "') " << endl; 
  
  
  lumi_file << "" << endl;
  lumi_file << "declare -a mc_noqcd=('" ;
  for(std::map<TString, TString>::iterator mit =lqmap.begin(); mit != lqmap.end();++mit){
    if(mit->second.Contains("QCD")) continue;
    cout << mit->first << " " << mit->second << endl;
    if(next(mit)!= lqmap.end()) lumi_file << mit->second << "' '"   ;
    else  lumi_file << mit->second ;
    
  }
  lumi_file << "') " << endl;
  
  lumi_file << "" << endl;
 lumi_file << "declare -a ch_wa=('" ;
 for(std::map<TString, TString>::iterator mit =lqmap.begin(); mit != lqmap.end();++mit){
   if(!mit->second.Contains("tthwA")) continue;

   cout << mit->first << " " << mit->second << endl;
   if(next(mit)!= lqmap.end()) lumi_file << mit->second << "' '"   ;
   else  lumi_file << mit->second ;

 }

 lumi_file << "') " << endl;
 lumi_file << "" << endl;
 
 lumi_file << "" << endl;
 lumi_file << "declare -a ch_wz=('" ;
 for(std::map<TString, TString>::iterator mit =lqmap.begin(); mit != lqmap.end();++mit){
   if(!mit->second.Contains("tthwZ")) continue;
   
   cout << mit->first << " " << mit->second << endl;
   if(next(mit)!= lqmap.end()) lumi_file << mit->second << "' '"   ;
   else  lumi_file << mit->second ;

 }




  lumi_file << "') " << endl;
  lumi_file << "" << endl;

 
  lumi_file << "declare -a hn_ll_ee=('" ;
  for(std::map<TString, TString>::iterator mit =lqmap.begin(); mit != lqmap.end();++mit){
    if(!mit->second.Contains("Schan")) continue;
    if(mit->second.Contains("_Mu")) continue;
    if(mit->second.Contains("_EpM")) continue;
    if(mit->second.Contains("_EmM")) continue;

    cout << mit->first << " " << mit->second << endl;
    if(next(mit)!= lqmap.end()) lumi_file << mit->second << "' '"   ;
    else  lumi_file << mit->second ;

  }
  lumi_file << "') " << endl;

  lumi_file << "declare -a hn_ll_mm=('" ;
  for(std::map<TString, TString>::iterator mit =lqmap.begin(); mit != lqmap.end();++mit){
    if(!mit->second.Contains("Schan")) continue;
    if(mit->second.Contains("_E")) continue;
    if(mit->second.Contains("_EpM")) continue;
    if(mit->second.Contains("_EmM")) continue;
    
    cout << mit->first << " " << mit->second << endl;
    if(next(mit)!= lqmap.end()) lumi_file << mit->second << "' '"   ;
    else  lumi_file << mit->second ;

  }

  lumi_file << "') " << endl;

  lumi_file << "" << endl;
  lumi_file << "declare -a hn_ll_em=('" ;
  for(std::map<TString, TString>::iterator mit =lqmap.begin(); mit != lqmap.end();++mit){
    if(!mit->second.Contains("Schan")) continue;
    if(mit->second.Contains("_EpE")) continue;
    if(mit->second.Contains("_EmE")) continue;
    if(mit->second.Contains("_Mu")) continue;

    cout << mit->first << " " << mit->second << endl;
    if(next(mit)!= lqmap.end()) lumi_file << mit->second << "' '"   ;
    else  lumi_file << mit->second ;

  }

  

  lumi_file << "') " << endl;





  lumi_file << "" << endl;
  lumi_file << "declare -a hn_mmm=('" ;
  for(std::map<TString, TString>::iterator mit =lqmap.begin(); mit != lqmap.end();++mit){
    if(!mit->second.Contains("_mumumu")) continue;
    if(!mit->second.Contains("HN")) continue;

    cout << mit->first << " " << mit->second << endl;
    if(next(mit)!= lqmap.end()) lumi_file << mit->second << "' '"   ;
    else  lumi_file << mit->second ;

  }



  lumi_file << "') " << endl;


  lumi_file << "" << endl;
  lumi_file << "declare -a hn_ll_tchann_ee=('" ;
  for(std::map<TString, TString>::iterator mit =lqmap.begin(); mit != lqmap.end();++mit){
    if(!mit->second.Contains("HN_Tchan")) continue;
    if(mit->second.Contains("mum")) continue;
    if(mit->second.Contains("mup")) continue;
    bool last=true;
    for(std::map<TString, TString>::iterator mit2 =next(mit); mit2 != lqmap.end();++mit2){
      if(mit2->second.Contains("HN_Tchan") && mit2->second.Contains("mu")) last=false;
    }
    if (last) lumi_file << mit->second << "' " ;
    else{
      if(next(mit)!= lqmap.end()) lumi_file << mit->second << "' '"   ;
      else  lumi_file << mit->second ;
    }
  }

  
  lumi_file << "') " << endl;



  lumi_file << "" << endl;
  lumi_file << "declare -a hn_ll_tchann_mm=('" ;
  for(std::map<TString, TString>::iterator mit =lqmap.begin(); mit != lqmap.end();++mit){
    if(!mit->second.Contains("HN_Tchan")) continue;
    if(mit->second.Contains("emem")) continue;
    if(mit->second.Contains("epep")) continue;
    bool last=true;
    for(std::map<TString, TString>::iterator mit2 =next(mit); mit2 != lqmap.end();++mit2){
      if(mit2->second.Contains("HN_Tchan") &&(mit2->second.Contains("em") ||  mit2->second.Contains("ep")) ) last=false;
    }
    if (last) lumi_file << mit->second << "' " ;
    else{
      if(next(mit)!= lqmap.end()) lumi_file << mit->second << "' '"   ;
      else  lumi_file << mit->second ;
    }


  }


  lumi_file << "') " << endl;


  
  lumi_file << "" << endl;
  lumi_file << "" << endl;
  lumi_file << "declare -a diboson_pythia=('WZ' 'ZZ' 'WW' 'WZTo3LNu_powheg' 'ZZTo4L_powheg' )" << endl;
  lumi_file << "" << endl;
  lumi_file << "declare -a dy_mcatnlo=('DYJets_10to50' 'DYJets') " << endl;
  lumi_file << "" << endl;
  lumi_file << "declare -a dilepton_list=('DYJets_10to50' 'DYJets' 'WJets' 'WpWpEWK' 'WpWpQCD'  'WZ' 'ZZ' 'WW'  'TTJets_MG'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t' 'SingleTbar_tW_noHadron'  'SingleTbar_tW' 'SingleTop_tW_noHadron' 'SingleTop_tW' 'WWW' 'ttW' 'ttZ' 'ttH_nonbb' 'ttH_bb' 'ttbb' 'ZZZ' 'WZZ'  'VBF_HToMuMu' 'WGtoLNuG' 'WGtoLNuEE' 'WGtoLNuMM' 'ZGto2LG' )" << endl;
  
  lumi_file << "declare -a trilepton_list=('DYJets_10to50' 'DYJets' 'WJets' 'WZ' 'ZZ' 'WW'  'TTJets_MG' 'SingleTop_s' 'SingleTbar_t' 'SingleTop_t' 'SingleTbar_tW_noHadron'  'SingleTbar_tW' 'SingleTop_tW_noHadron' 'SingleTop_tW' 'WZZ' 'ZZZ' 'WWW'  'WZTo3LNu_powheg' 'ZZTo4L_powheg' )" << endl;

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
    if(mit->second.Contains("WpWp")) continue;    
    if(!mit->second.Contains("MuEnriched")) continue;
    bool last=true;
    for(std::map<TString, TString>::iterator mit2 =next(mit); mit2 != lqmap.end();++mit2){
      if(mit2->second.Contains("MuEnriched") && mit2->second.Contains("QCD") ) last=false;
    }   
    if(!last)lumi_file << mit->second << "' '"   ;
    else  lumi_file << mit->second ;
  } lumi_file << "') " << endl;
  lumi_file << "" << endl;
  
  lumi_file << "declare -a qcd_eg=('" ;
  for(std::map<TString, TString>::iterator mit =lqmap.begin(); mit != lqmap.end();++mit){
    if(!mit->second.Contains("QCD")) continue;
    if(mit->second.Contains("WpWp")) continue;
    if(mit->second.Contains("MuEnriched")) continue;
    bool last=true;
    for(std::map<TString, TString>::iterator mit2 =next(mit); mit2 != lqmap.end();++mit2){
      if(mit2->second.Contains("QCD") && mit2->second.Contains("EMEnriched")) last=false;
    }
    if(!last)lumi_file << mit->second << "' '"   ;
    else  lumi_file << mit->second ;
  } lumi_file << "') " << endl;
  lumi_file << "" << endl;


  lumi_file << "" << endl;
  lumi_file << "declare -a hn_mm=('WJets' 'WZ' 'ZZ' 'WW'  'TTJets_MG' 'WpWpEWK' 'WpWpQCD'  'WZZ' 'ZZZ' 'WWW'  'WZTo3LNu_powheg' 'ZZTo4L_powheg' 'ttZ' 'ttW'  'ttH_nonbb' 'ttH_bb' 'ttbb' 'ZZZ' 'WZZ'  'VBF_HToMuMu' 'WGtoLNuG' 'WGtoLNuEE' 'WGtoLNuMM' 'ZGto2LG' )  " << endl;    
  lumi_file << "" << endl;
  lumi_file << "declare -a hn_ee=('DYJets_10to50' 'DYJets' 'WJets' 'WZ' 'ZZ' 'WW'  'TTJets_MG' 'WpWpEWK' 'WpWpQCD'  'WZZ' 'ZZZ' 'WWW'  'WZTo3LNu_powheg' 'ZZTo4L_powheg' 'ttZ' 'ttW'  'ttH_nonbb' 'ttH_bb' 'ttbb' 'ZZZ' 'WZZ'  'VBF_HToMuMu' 'WGtoLNuG' 'WGtoLNuEE' 'WGtoLNuMM' 'ZGto2LG' )" << endl;    
  lumi_file << "" << endl;
  lumi_file << "declare -a hn_fakeee=('DYJets_10to50' 'DYJets' 'WJets'  'TTJets_MG')" << endl;
  lumi_file << "" << endl;
  lumi_file << "" << endl;

//lumi_file << "declare -a singletop=('SingleTop_s' 'SingleTop_tW_noHadron' 'SingleTop_tW' 'SingleTbar_s' 'SingleTbar_tW_noHadron' 'SingleTbar_tW') " << endl;

  lumi_file << "declare -a singletop=('" ;
  for(std::map<TString, TString>::iterator mit =lqmap.begin(); mit != lqmap.end();++mit){
    if(!mit->second.Contains("SingleTop")) continue;
    bool last=true;
    for(std::map<TString, TString>::iterator mit2 =next(mit); mit2 != lqmap.end();++mit2){
      if(mit2->second.Contains("SingleTop")) last=false;
    }
    if(!last)lumi_file << mit->second << "' '"   ;
    else  lumi_file << mit->second ;
  } lumi_file << "') " << endl;
  lumi_file << "" << endl;



  string lqdir = getenv("LQANALYZER_DIR");
  TString user = TString(getenv("USER"));


  string lfile2 =  lqdir + "/LQRun/txt/list_all_mc_" + string(def_version) +".sh";


  gSystem->Exec(("cp " + lfile + " "+ string(getenv("LQANALYZER_DATASET_DIR"))).Data());
  if (user.Contains("jalmond"))  gSystem->Exec(("chmod 777 " + string(getenv("LQANALYZER_DATASET_DIR")) + lfile).Data());
  gSystem->Exec(("mv " + lfile +" " + lfile2).Data());


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
