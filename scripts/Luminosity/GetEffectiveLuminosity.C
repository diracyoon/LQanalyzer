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

#include <unistd.h>

using namespace std;

#include <map>

#include "SampleMap.C"

float GetEventsProcessed(std::string file);
float GetEventsPassed(std::string file);
float GetSumWeights(std::string filename);

map<TString, Double_t> map_lumi;
map<TString, Double_t> neventmap;
map<TString, Double_t> n_w_eventmap;
void GetEffectiveLuminosity(TString path_of_list, int sig,  TString tag,TString version="") {
  
  bool debug = true;
  map_lumi.clear();
  neventmap.clear();
  n_w_eventmap.clear();
  std::vector<TString> missing_samples;
  std::vector<TString> missing_samples2;

  cout << "path_of_list = " << path_of_list << endl;
  cout << "tag =  " << tag << endl;
  cout << "version = " << version << endl;
  
  //TString path_of_list="/data1/LQAnalyzer_rootfiles_for_analysis/DataSetLists/cattuplist_"+TString(getenv("CATVERSION"))+".txt";
  
  bool NewList(true);
  if (path_of_list.Contains(string(getenv("LQANALYZER_DATASET_DIR")) +"/ca")) NewList=false;

  if (NewList) cout << "New list " << endl;
  else cout <<"Not new list" << endl;

  if(CheckMaps(path_of_list, sig)) return;
  TString def_version = TString(getenv("CATVERSION"));
  if(!version.Contains("v8") ) version = def_version;

  map<TString, TString> missing_map;
  vector<TString> vec_available ;
  if(!NewList){
    missing_map= GetMissingMap(version,path_of_list, sig);
    vec_available = GetAvailableMap(version,path_of_list, sig);
  }

  cout << "GetEffectiveLuminosity: Using file: " << path_of_list << " @ " << endl;
  map<TString, TString> datasets =  GetDatasetNames(version,path_of_list, sig);
  
  
  map<TString, Double_t> dirmap = GetXSecMap2016(path_of_list, sig); 
  map<TString, TString> lqmap = GetLQMap2016(path_of_list, sig);
  map<TString, TString> trilepmap = GetTriLepMap2016(path_of_list, sig);
  
  TString analysisdir = TString(getenv("HOSTNAME"));
  bool cluster = false;
  
  if(analysisdir.Contains("cmscluster.snu.ac.kr")) cluster=true;
  
  for(std::map<TString, Double_t>::iterator mit =dirmap.begin(); mit != dirmap.end();++mit){

    TString path_file=" /data2/DATA/cattoflat/MC/";

    ///   IF NEED LATER TO REDICE DISK SPACE CAN MAKE GENSKIM TO REPLACE CATTUPLE
    //    if(mit->first.Contains("amcatnlo")) path_file="/data2/DATA/cattoflat/genskim/";
      
    TString dir = "ls "+path_file + version + "/"+ mit->first + "/*.root > inputlist_efflumi.txt";
    if(cluster) dir = "ls /data4/DATA/FlatCatuples/MC/" + version + "/"+ mit->first + "/*.root > inputlist_efflumi.txt";
    
    bool use_sum_genweight(false);
    cout << "Checking " << mit->first << endl;
    if(mit->first.Contains("amcatnlo")) use_sum_genweight=true;
    else if(mit->first.Contains("_Schannel_")){
      if(mit->first.Contains("HN") || mit->first.Contains("Majorana") )  use_sum_genweight=true;
    }
    else if(mit->first.Contains("eavyNeutrino_trilepton")){
      use_sum_genweight=true;
    }
      else use_sum_genweight=false;
    
    
    cout << "Running " << dir << endl;
    system(dir.Data());
    
    cout << "Finished Running " << dir << endl;

    std::ifstream fin("inputlist_efflumi.txt");
    std::string word;
    
    float number_events_processed(0.);
    float number_events_passed(0.);
    float sum_of_weights(0.);
    
    if(!use_sum_genweight) {
      while ( fin >> word ) {
	number_events_processed+= GetEventsProcessed(word);
	number_events_passed+= GetEventsPassed(word);
      }
      sum_of_weights = number_events_processed;
      if(number_events_processed ==0 ) {
	missing_samples.push_back(mit->first);
	missing_samples2.push_back(mit->second);
	continue;
      }
    }
    else{
      std::vector<std::string> filelist;
      while ( fin >> word ) {
	number_events_processed+= GetEventsProcessed(word);;
        number_events_passed+= GetEventsPassed(word);
	filelist.push_back(word);
      }
      if(number_events_processed ==0 ) {
	missing_samples.push_back(mit->first);
	missing_samples2.push_back(mit->second);
	continue;
      }
      
      
      TString command1 = "rm -r /data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/" + mit->first;
      TString command2 = "mkdir /data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"+ mit->first;
      TString command2b = "mkdir /data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"+ mit->first + "/output/";
      
      cout << "Running " << command1.Data() << endl;
      system(command1.Data());
      cout << "Running " << command2.Data() << endl;

      system(command2.Data());
      cout << "Running " << command2b.Data() << endl;
      system(command2b.Data());
      
      for(unsigned int i=0; i < filelist.size(); i++){
	std::string istr;
	std::stringstream out;
	out << i;
	istr = out.str();
	
	string lqdir = getenv("LQANALYZER_DIR");
	ofstream outputFile;
	outputFile.open(("Run_" + istr + ".sh").c_str());
	outputFile << "#!/bin/sh" << endl;
	outputFile <<"\n"<< endl;
	outputFile <<"#$ -S /bin/bash \n"<< endl;
	outputFile <<"#$ -N Run_" + istr  + "  \n"<< endl;
	outputFile <<"#$ -wd  /data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"+mit->first + " \n"<< endl;
	outputFile <<"#$ -o  /data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"+mit->first + " \n"<< endl;
	outputFile <<"#$ -e  /data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"+mit->first + " \n"<< endl;
	outputFile <<"echo 'Job started at ' `date` \n"<< endl;
	outputFile <<"cd /share/apps/root_v5-34-32/root/ \n"<< endl;
	outputFile <<". bin/thisroot.sh \n"<< endl;
	outputFile <<"cd " + lqdir + "\n"<< endl;
	
	outputFile <<"source setup.sh \n"<< endl;
	outputFile <<"echo 'PWD= '$PWD \n"<< endl;
	outputFile <<"cd /data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"+mit->first + " \n"<< endl;
	outputFile <<"root -l -b -q \'CountGenWeights.C(\"/data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"+mit->first+ "\",\""+filelist.at(i)+"\",\""+ "hist" + TString(istr) +".root\")\' "<< endl;
	outputFile <<"echo 'Ran macro 2' \n"<< endl;
	outputFile <<""<< endl;
	TString command4a="cp Run_" + istr + ".sh  /data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"+mit->first ;
        system(command4a.Data());
	TString command4b="cp CountGenWeights.C /data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"+mit->first ;
	system(command4b.Data());
	
	TString command3 = "qsub -V  /data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"+mit->first+ "/Run_" + istr + ".sh";
	if(debug)cout << command3 << endl;
	system(command3.Data());
	stringstream ss;
	ss << istr;
	string strI = ss.str();

	system(("rm  Run_" + strI + ".sh").c_str());
      }
      
      system("qstat");
      system("qstat > logqstat");

      cout << "Checking if job is complete...." << endl;
      bool jobComplete=false;
      int nsubmits=0;
      while (!jobComplete){
	jobComplete=true;
	for(unsigned int i=0; i < filelist.size(); i++){
	  std::string istr;
	  std::stringstream out;
	  out << i;

	  istr = out.str();
	  std::ifstream infile("/data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"+mit->first +"/output/hist" + TString(istr) +".root");
	  if(!infile.good()) {
	      jobComplete=false;
	      if(nsubmits>200){
		cout << "File /data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"+mit->first +"/output/hist" + TString(istr) +".root does not exist" << endl;  
		
		system("root -l -b -q \'CountGenWeights.C(\"/data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"+mit->first+ "\",\""+filelist.at(i)+"\",\""+ "hist" + TString(istr) +".root\")\'");
	      }
	      sleep(5);
	      nsubmits++;
	      break;
	  }
	}
      }
      cout << "Job complete" << endl;

      
      TString command4 = "rm log/checkoutput.txt";
      TString command5= "ls   "+path_file + version + "/"+  mit->first + "/  > log/checkoutput.txt";
      if(cluster) command5= "ls   /data4/DATA/FlatCatuples/MC/" + version + "/"+ mit->first + "/  > log/checkoutput.txt";
      system(command4.Data());
      system(command5.Data());
      TString filename_fs = "log/checkoutput.txt";
      ifstream name_file(filename_fs.Data());
      int counter=0;
      while(!name_file.eof()) {
	string filen;
	name_file >> filen;
	if(TString(filen).Contains(".root"))counter++;
      }
      if(debug){
	if(!cluster)cout << "Number of files in "+path_file + version + "/"+  mit->first + "/ = " << counter << endl;
	else cout << "Number of files in //data4/DATA/FlatCatuples/MC/" + version + "/"+  mit->first + "/ = " << counter << endl;
      }
      bool JobDone=false;
      while (JobDone==false){
	TString command6= "ls   /data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"  +  mit->first + "/output/ > log/checkcounted.txt";
	system(command6.Data());
	
	TString filename_counted = "log/checkcounted.txt";
	ifstream name_countedfile(filename_counted.Data());
	int counter_counted=0;
	while(!name_countedfile.eof()) {
	  string filen;
	  name_countedfile >> filen;
	  if(TString(filen).Contains(".root"))counter_counted++;
	}
	if(counter_counted == counter) JobDone=true;
      }
      TString haddcommand = "hadd   /data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"+mit->first + "/output/Output.root  /data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"+mit->first + "/output/*.root ";
      system(haddcommand.Data());
      TFile* file = TFile::Open(( "/data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"+mit->first + "/output/Output.root").Data());
      TH1F*  SumWCounter = (TH1F*) (file ->Get("sumweight"));
      sum_of_weights = SumWCounter->Integral();
    }

    if(number_events_processed ==0 ) continue;
    fin.close();
       
    float lumi = number_events_processed / mit->second;
    if(use_sum_genweight) lumi = sum_of_weights / mit->second;
    

    std::cout.precision(10);
    
    //if(debug){
      std::cout << "/data2/LQ_SKTreeOutput/Lumi/" + TString(getenv("USER")) + "/"+mit->first << "    nevents =  " << number_events_processed << " sum of weights =  " << sum_of_weights << " eff lumi = " << lumi <<std::endl;
      //}

    system("rm inputlist_efflumi.txt");
    
    map_lumi[ mit->first] = lumi;
    neventmap[ mit->first]=number_events_processed;
    n_w_eventmap[mit->first]=sum_of_weights;
  } 

  
  ofstream lumi_file;
  string lfile =  "datasets_snu_nonsig_CAT_mc_" + string(version.Data()) + string(tag)+".txt";
  if(sig==1) lfile =  "datasets_snu_sig_CAT_mc_" + string(version.Data()) + string(tag)+".txt";
  if(sig==0) lfile =  "datasets_snu_CAT_mc_" + string(version.Data()) + string(tag)+".txt";

  if(cluster) lfile =  "datasets_snu_cluster_CAT_mc_" + string(version.Data()) + ".txt";

  lumi_file.open(lfile.c_str());
  lumi_file.setf(ios::fixed,ios::floatfield); 
  lumi_file.precision(5);


  if(!NewList){
    lumi_file << "########################" << endl;
    lumi_file << "### CATTUPLES ########## "  << endl;
    lumi_file << "#######################" << endl;
    lumi_file << "" << endl;
    lumi_file << "" << endl;
    lumi_file << "#######################################################################################################################################" << endl;
    lumi_file << "# sample     Nevents  Sum(GenWeights)  xsec*Eff   Effective lumi    LocalPath #####################################################################################" << endl;
    lumi_file << "#######################################################################################################################################" << endl;
    lumi_file << "" << endl;
    lumi_file << "" << endl;
    lumi_file << "#### CATTuples" << endl;
    for(std::map<TString, Double_t>::iterator mit =map_lumi.begin(); mit != map_lumi.end();++mit){
      std::map<TString, TString>::iterator mit2 = lqmap.find(mit->first);    
      std::map<TString, Double_t>::iterator mit3 = dirmap.find(mit->first);    
      std::map<TString, Double_t>::iterator mit4 = neventmap.find(mit->first);    
      std::map<TString, Double_t>::iterator mit5 = n_w_eventmap.find(mit->first);    
      if(!cluster)lumi_file <<  mit2->second << "  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << "  /data2/DATA/cattoflat/MC/" << version <<"/"  << mit->first << "/" <<endl;
      else lumi_file <<  mit2->second << "  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << "  /data4/DATA/FlatCatuples/MC/" << version <<"/"  << mit->first << "/" <<endl;
      
    }
    
    lumi_file << "" << endl;
    lumi_file << "" << endl;
    lumi_file << "#### Single_lepton_skims:_SKTrees" << endl;
    for(std::map<TString, Double_t>::iterator mit =map_lumi.begin(); mit != map_lumi.end();++mit){
      std::map<TString, TString>::iterator mit2 = lqmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit3 = dirmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit4 = neventmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit5 = n_w_eventmap.find(mit->first);
      if(cluster)lumi_file <<  "SK" << mit2->second << "  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << " /data4/LocalNtuples/SKTrees13TeV/" + string(version.Data()) +"/SKTrees/MC/" << mit2->second << "/" <<endl;
      else lumi_file <<  "SK" << mit2->second << "  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << " /data2/CatNtuples/" + string(version.Data()) +"/SKTrees/MC/" << mit2->second << "/" <<endl;
      
    }
    
    
    lumi_file << "" << endl;
    lumi_file << "" << endl;
    lumi_file << "#### Dilepton_skims:_SKTrees" << endl;
    for(std::map<TString, Double_t>::iterator mit =map_lumi.begin(); mit != map_lumi.end();++mit){
      std::map<TString, TString>::iterator mit2 = lqmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit3 = dirmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit4 = neventmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit5 = n_w_eventmap.find(mit->first);
      
      lumi_file <<  "SK" << mit2->second << "_dilep  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << " /data2/CatNtuples/" + string(version.Data()) +"/SKTrees/MCDiLep/" <<  mit2->second << "/" <<endl;
    }
    lumi_file << "" << endl;
    lumi_file << "" << endl;
    lumi_file << "#### HNDilepton_skims:_SKTrees" << endl;
    for(std::map<TString, Double_t>::iterator mit =map_lumi.begin(); mit != map_lumi.end();++mit){
      std::map<TString, TString>::iterator mit2 = lqmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit3 = dirmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit4 = neventmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit5 = n_w_eventmap.find(mit->first);

      lumi_file <<  "SK" << mit2->second << "_hndilep  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << " /data2/CatNtuples/" + string(version.Data()) +"/SKTrees/MCDiLep/" <<  mit2->second << "/" <<endl;
    }
    
    lumi_file << "" << endl;
    lumi_file << "#### HNFake_skims:_SKTrees" << endl;

    for(std::map<TString, Double_t>::iterator mit =map_lumi.begin(); mit != map_lumi.end();++mit){
      std::map<TString, TString>::iterator mit2 = lqmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit3 = dirmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit4 = neventmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit5 = n_w_eventmap.find(mit->first);
      
      lumi_file <<  "SK" << mit2->second << "_hnfake  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << " /data2/CatNtuples/" + string(version.Data()) +"/SKTrees/MCHNFake/" <<  mit2->second << "/" <<endl;
    }


    lumi_file << "" << endl;
    lumi_file << "" << endl;
    lumi_file << "#### Trilepton_skims:_SKTrees" << endl;
    for(std::map<TString, Double_t>::iterator mit =map_lumi.begin(); mit != map_lumi.end();++mit){
      std::map<TString, TString>::iterator mit2 = lqmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit3 = dirmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit4 = neventmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit5 = n_w_eventmap.find(mit->first);
      std::map<TString, TString>::iterator check_trilep= trilepmap.find(mit->first);
      //if(check_trilep != trilepmap.end())
     lumi_file <<  "SK" << mit2->second << "_trilep  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << " /data2/CatNtuples/" + string(version.Data()) +"/SKTrees/MCTriLep/" <<  mit2->second << "/" <<endl;
    }
    
    
    
    lumi_file << "" << endl;
    lumi_file << "#### NoCut skims: SKTrees" << endl;
    for(std::map<TString, Double_t>::iterator mit =map_lumi.begin(); mit != map_lumi.end();++mit){
      std::map<TString, TString>::iterator mit2 = lqmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit3 = dirmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit4 = neventmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit5 = n_w_eventmap.find(mit->first);
      lumi_file <<  "SK" << mit2->second << "_nocut  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << " /data2/CatNtuples/" + string(version.Data()) +"/SKTrees/MCNoCut/" <<  mit2->second << "/" <<endl;
    }
    
    lumi_file << "" << endl;
    lumi_file << "##################################################################" << endl;
    lumi_file << "#### Missing/Not produced samples in this version are listed below " << endl;
    lumi_file << "#### Missing : means miniAOD not available " << endl;
    lumi_file << "#### Available: means miniAOD/catuples available but no flatcatuples/sktree made. " << endl;
    lumi_file << "##################################################################" << endl;
    
    for(map<TString, TString>::iterator it = missing_map.begin(); it!= missing_map.end(); it++){
      lumi_file << "Missing: " << it->first << "   " << it->second << endl;
    }
    for(vector<TString>::iterator it= vec_available.begin(); it!= vec_available.end(); it++){
      lumi_file <<"Available[not produced]: " << *it << endl;
    }
    
    
    lumi_file << "" << endl;
    lumi_file << "##################################################################" << endl;
    lumi_file << "#### FULL Name of MINIAODS " << endl;
    lumi_file << "##################################################################" << endl;
    
    for(map<TString, TString>::iterator it = datasets.begin(); it!= datasets.end(); it++){
      lumi_file << "DATASET: " << it->first << "   " << it->second << endl;
    }
    
    lumi_file << "################################################################ " << endl;
    lumi_file << "#### Private produced samples : Not made in batch at kisti" << endl;
    lumi_file << "##################################################################" << endl;
  }
  else{
    for(std::map<TString, Double_t>::iterator mit =map_lumi.begin(); mit != map_lumi.end();++mit){
      std::map<TString, TString>::iterator mit2 = lqmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit3 = dirmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit4 = neventmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit5 = n_w_eventmap.find(mit->first);
      if(!cluster)lumi_file <<  mit2->second << "  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << "  /data2/DATA/cattoflat/MC/" << version <<"/"  << mit->first << "/" <<endl;
      else lumi_file <<  mit2->second << "  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << "  /data4/DATA/FlatCatuples/MC/" << version <<"/"  << mit->first << "/" <<endl;

    }

    for(std::map<TString, Double_t>::iterator mit =map_lumi.begin(); mit != map_lumi.end();++mit){
      std::map<TString, TString>::iterator mit2 = lqmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit3 = dirmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit4 = neventmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit5 = n_w_eventmap.find(mit->first);
      if(cluster)lumi_file <<  "SK" << mit2->second << "  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << " /data4/LocalNtuples/SKTrees13TeV/" + string(version.Data()) +"/SKTrees/MC/" << mit2->second << "/" <<endl;
      else lumi_file <<  "SK" << mit2->second << "  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << " /data2/CatNtuples/" + string(version.Data()) +"/SKTrees/MC/" << mit2->second << "/" <<endl;

    }
    for(std::map<TString, Double_t>::iterator mit =map_lumi.begin(); mit != map_lumi.end();++mit){
      std::map<TString, TString>::iterator mit2 = lqmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit3 = dirmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit4 = neventmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit5 = n_w_eventmap.find(mit->first);

      lumi_file <<  "SK" << mit2->second << "_dilep  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << " /data2/CatNtuples/" + string(version.Data()) +"/SKTrees/MCDiLep/" <<  mit2->second << "/" <<endl;
    }
    for(std::map<TString, Double_t>::iterator mit =map_lumi.begin(); mit != map_lumi.end();++mit){
      std::map<TString, TString>::iterator mit2 = lqmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit3 = dirmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit4 = neventmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit5 = n_w_eventmap.find(mit->first);

      lumi_file <<  "SK" << mit2->second << "_trilep  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << " /data2/CatNtuples/" + string(version.Data()) +"/SKTrees/MCTriLep/" <<  mit2->second << "/" <<endl;
    }
    for(std::map<TString, Double_t>::iterator mit =map_lumi.begin(); mit != map_lumi.end();++mit){
      std::map<TString, TString>::iterator mit2 = lqmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit3 = dirmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit4 = neventmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit5 = n_w_eventmap.find(mit->first);
      
      lumi_file <<  "SK" << mit2->second << "_hndilep  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << " /data2/CatNtuples/" + string(version.Data()) +"/SKTrees/MCHNDiLep/" <<  mit2->second << "/" <<endl;
																					       
    }
    
    for(std::map<TString, Double_t>::iterator mit =map_lumi.begin(); mit != map_lumi.end();++mit){
      std::map<TString, TString>::iterator mit2 = lqmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit3 = dirmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit4 = neventmap.find(mit->first);
      std::map<TString, Double_t>::iterator mit5 = n_w_eventmap.find(mit->first);

      lumi_file <<  "SK" << mit2->second << "_hnfake  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << " /data2/CatNtuples/" + string(version.Data()) +"/SKTrees/MCHNFake/" <<  mit2->second << "/" <<endl;
    }

  }


  string lqdir = getenv("LQANALYZER_DIR");
  string lfile2 =  lqdir + "/LQRun/txt/datasets_snu_nonsig_CAT_mc_" + string(version.Data()) + ".txt";
  if(sig==1) lfile2 =  lqdir + "/LQRun/txt/datasets_snu_sig_CAT_mc_" + string(version.Data()) + ".txt";
  if(sig==0) lfile2 =  lqdir + "/LQRun/txt/datasets_snu_CAT_mc_" + string(version.Data()) + ".txt";

  if(cluster) lfile2 =  lqdir + "/LQRun/txt/Cluster/datasets_snu_cluster_CAT_mc_" + string(version.Data()) + ".txt";

  TString user = TString(getenv("USER"));
  if(!cluster)gSystem->Exec(("cp " + lfile + "  " + getenv("LQANALYZER_DATASETFILE_DIR")+ "/").c_str()); 
  else gSystem->Exec(("cp " + lfile + "  /data4/LocalNtuples/LQAnalyzer_rootfiles_for_analysis/CATAnalysis2016/").c_str());

  //gSystem->Exec(("mv " + lfile +" " + lfile2).c_str());

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

float GetSumWeights(std::string filename){
  TFile* file = TFile::Open(filename.c_str());

  //cout << "GetSumWeights :: " << endl;
  TDirectory * dir = (TDirectory*)file->Get(TString(filename) + ":/ntuple");
  TTree * tree;
  dir->GetObject("event",tree);
  
  float genWeight=0;
  TBranch        *b_genWeight; 
  
  tree->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
  
  float sum_weight=0.;
  Long64_t nentries = tree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    tree->LoadTree(jentry);
    nb = tree->GetEntry(jentry);
    sum_weight += genWeight;
  }
  //cout << "sum_weight = " << sum_weight << endl;
  return sum_weight;
}
