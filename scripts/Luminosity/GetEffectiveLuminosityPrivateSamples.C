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

float GetEventsProcessed(std::string file);
float GetEventsPassed(std::string file);
float GetSumWeights(std::string filename);

map<TString, TString>  map_cv;
map<TString, Double_t> map_lumi;
map<TString, Double_t> neventmap;
map<TString, Double_t> n_w_eventmap;
void GetEffectiveLuminosityPrivateSamples(){
  

  TString version="";;
   
  
  map_cv.clear();
  map_lumi.clear();
  neventmap.clear();
  n_w_eventmap.clear();
  
  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////
  /// This is an example... xsmap/lqmap and catversion_map need to be filled  with private sample information
  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////  ////
  
  map<TString, Double_t> xsmap;
  /// key=  "Name of FlatCatuple dir"  = xsec;
  //xsmap["VBF_HToMuMu_M125_13TeV_powheg_pythia8"] =0.000828308;
  //xsmap["13TeV_HN40_mumumu_width_AUTO_VmuN_0p1"] = 8.605;
  //xsmap["13TeV_HN60_mumumu_width_AUTO_VmuN_0p1"] = 3.416;
  //xsmap["13TeV_HN150_mumumu_width_AUTO_VmuN_0p1"] = 0.0004561;
  //xsmap["13TeV_HN700_mumumu_width_AUTO_VmuN_0p1"] = 0.000002644;
  xsmap["CHToCB_M125_madgraph_13TeV_2016"] = 1;

  map<TString, TString> lqmap;
  //  key=  "Name of FlatCatuple dir"  = name to use in sktree -i <SAMPLENAME> 
  //lqmap["VBF_HToMuMu_M125_13TeV_powheg_pythia8"] ="vhf_Htomm_Powheg";
  // lqmap["13TeV_HN40_mumumu_width_AUTO_VmuN_0p1"] ="HN40_mumumu_VmuN_0p1";
  // lqmap["13TeV_HN60_mumumu_width_AUTO_VmuN_0p1"] ="HN60_mumumu_VmuN_0p1";
  // lqmap["13TeV_HN150_mumumu_width_AUTO_VmuN_0p1"] ="HN150_mumumu_VmuN_0p1";
  // lqmap["13TeV_HN700_mumumu_width_AUTO_VmuN_0p1"] ="HN700_mumumu_VmuN_0p1";
  lqmap["CHToCB_M125_madgraph_13TeV_2016"] = "CHToCB_M125_madgraph_13TeV_2016";

  //// Version must be the version in the directoy name of your flat catuples 
  map<TString, TString> catversion_map;
  
  //catversion_map["VBF_HToMuMu_M125_13TeV_powheg_pythia8"] = "v7-6-4"; /// This specifies what directory to look in
  // catversion_map["13TeV_HN40_mumumu_width_AUTO_VmuN_0p1"] = "v7-6-4";
  // catversion_map["13TeV_HN60_mumumu_width_AUTO_VmuN_0p1"] = "v7-6-4";
  // catversion_map["13TeV_HN150_mumumu_width_AUTO_VmuN_0p1"] = "v7-6-4";
  // catversion_map["13TeV_HN700_mumumu_width_AUTO_VmuN_0p1"] = "v7-6-4";
  catversion_map["CHToCB_M125_madgraph_13TeV_2016"] = "v8-0-6";
  
  for(std::map<TString, Double_t>::iterator mit =xsmap.begin(); mit != xsmap.end();++mit){
    std::map<TString, TString>::iterator mitv = catversion_map.find(mit->first);
    if(mitv== catversion_map.end()) {cout << "Error in naming datasets in map" << endl; return;}
    version= mitv->second;
    
    TString dir = "ls /data2/DATA/cattoflat/MC/" + version + "/"+ mit->first + "/*.root > inputlist_private.txt";
    
    bool use_sum_genweight(false);
    if(mit->first.Contains("amcatnlo")) use_sum_genweight=true;
    else use_sum_genweight=false;
    
    system(dir.Data());
    
    
    std::ifstream fin("inputlist_private.txt");
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
	cout << "Problem with input files. Are they in the correct local path?" << endl;
        cout <<"FlatCatuples: /data2/DATA/cattoflat/MC/CATVERSION/" << endl;
        cout <<"SKTrees: Lepton Skim : /data2/CatNtuples/CATVERSION/SKTrees/MC/" << endl;
        cout <<"SKTrees:DiLepton Skim : /data2/CatNtuples/CATVERSION/SKTrees/MCDiLep/" << endl;
        return;
      }
    }
    else{
      std::vector<std::string> filelist;
      while ( fin >> word ) {
	number_events_processed+= GetEventsProcessed(word);
        number_events_passed+= GetEventsPassed(word);
	filelist.push_back(word);
      }
      if(number_events_processed ==0 ) {
	cout << "Problem with input files. Are they in the correct local path?" << endl;
	cout <<"FlatCatuples: /data2/DATA/cattoflat/MC/CATVERSION/" << endl;
	cout <<"SKTrees: Lepton Skim : /data2/CatNtuples/CATVERSION/SKTrees/MC/" << endl;
	cout <<"SKTrees:DiLepton Skim : /data2/CatNtuples/CATVERSION/SKTrees/MCDiLep/" << endl;
	return;
      }
      TString command1 = "rm -r "+ mit->first;
      TString command2 = "mkdir "+ mit->first;
      TString command2b = "mkdir "+ mit->first + "/output/";

      system(command1.Data());
      system(command2.Data());
      system(command2b.Data());
      
      for(unsigned int i=0; i < filelist.size(); i++){
      std::string istr;
      std::stringstream out;
      out << i;
      istr = out.str();
      TString command3 = "root -l -b -q \'CountGenWeights.C(\""+mit->first+ "\",\""+filelist.at(i)+"\",\""+ "hist" + TString(istr) +".root\")\' &";
      TString command3b = "root -l -b -q \'CountGenWeights.C(\""+mit->first+ "\",\""+filelist.at(i)+"\",\""+ "hist" + TString(istr) +".root\")\' ";
	cout << command3 << endl;
	if(i==0) system(command3.Data());
	else if(i%10) system(command3.Data());
	else system(command3b.Data());
      }
      
      
      TString command4 = "rm log/checkoutput_private.txt";
      TString command5= "ls   /data2/DATA/cattoflat/MC/" + version + "/"+ mit->first + "/  > log/checkoutput_private.txt";
      system(command4.Data());
      system(command5.Data());
      TString filename_fs = "log/checkoutput_private.txt";
      ifstream name_file(filename_fs.Data());
      int counter=0;
      while(!name_file.eof()) {
	string filen;
	name_file >> filen;
	if(TString(filen).Contains(".root"))counter++;
      }
      cout << "Number of files in /data2/DATA/cattoflat/MC/" + version + "/"+ mit->first + "/ = " << counter << endl;
     
      bool JobDone=false;
      while (JobDone==false){
	TString command6= "ls  "  + mit->first + "/output/ > log/checkcounted_private.txt";
	system(command6.Data());
	
	TString filename_counted = "log/checkcounted_private.txt";
	ifstream name_countedfile(filename_counted.Data());
	int counter_counted=0;
	while(!name_countedfile.eof()) {
	  string filen;
	  name_countedfile >> filen;
	  if(TString(filen).Contains(".root"))counter_counted++;
	}
	if(counter_counted == counter) JobDone=true;
      }
      TString haddcommand = "hadd  " + mit->first + "/output/Output.root " +  mit->first + "/output/*.root ";
      system(haddcommand.Data());
      TFile* file = TFile::Open((mit->first + "/output/Output.root").Data());
      TH1F*  SumWCounter = (TH1F*) (file ->Get("sumweight"));
      sum_of_weights = SumWCounter->Integral();
    }

    if(number_events_processed ==0 ) continue;
    fin.close();
       
    float lumi = number_events_processed / mit->second;
    if(use_sum_genweight) lumi = sum_of_weights / mit->second;
    
    std::cout.precision(10);
    std::cout <<mit->first << "    nevents =  " << number_events_processed << " sum of weights =  " << sum_of_weights << " eff lumi = " << lumi <<std::endl;
    
    system("rm inputlist_private.txt");
    system(("rm -r " +  mit->first).Data());
    
    map_lumi[mit->first] = lumi;
    neventmap[mit->first]=number_events_processed;
    n_w_eventmap[mit->first]=sum_of_weights;
    map_cv[mit->first] = version;
  } 

  
  ofstream lumi_file;
  lumi_file.precision(10);

  string lfile =  "datasets_snu_CAT_mc_" + string(version.Data()) + ".txt";
  lumi_file.open(lfile.c_str());
  lumi_file.setf(ios::fixed,ios::floatfield); 
  lumi_file.precision(1);

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
  for(std::map<TString, Double_t>::iterator mit =map_lumi.begin(); mit != map_lumi.end();++mit){
    std::map<TString, TString>::iterator mit2 = lqmap.find(mit->first);    
    std::map<TString, Double_t>::iterator mit3 = xsmap.find(mit->first);    
    std::map<TString, Double_t>::iterator mit4 = neventmap.find(mit->first);    
    std::map<TString, Double_t>::iterator mit5 = n_w_eventmap.find(mit->first);    
    lumi_file <<  mit2->second << "  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << "  /data2/DATA/cattoflat/MC/" << version <<"/"  << mit->first << "/" <<endl;
  }

  lumi_file << "" << endl;
  lumi_file << "" << endl;
  lumi_file << "#### Single lepton skims: SKTrees" << endl;
  for(std::map<TString, Double_t>::iterator mit =map_lumi.begin(); mit != map_lumi.end();++mit){
    std::map<TString, TString>::iterator mit2 = lqmap.find(mit->first);
    std::map<TString, Double_t>::iterator mit3 = xsmap.find(mit->first);
    std::map<TString, Double_t>::iterator mit4 = neventmap.find(mit->first);
    std::map<TString, Double_t>::iterator mit5 = n_w_eventmap.find(mit->first);
    lumi_file <<  "SK" << mit2->second << "  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << " /data2/CatNtuples/" + string(version.Data()) +"/SKTrees/MC/" << mit2->second << "/" <<endl;
  }

  
  lumi_file << "" << endl;
  lumi_file << "" << endl;
  lumi_file << "#### Single lepton skims: SKTrees" << endl;
  for(std::map<TString, Double_t>::iterator mit =map_lumi.begin(); mit != map_lumi.end();++mit){
    std::map<TString, TString>::iterator mit2 = lqmap.find(mit->first);
    std::map<TString, Double_t>::iterator mit3 = xsmap.find(mit->first);
    std::map<TString, Double_t>::iterator mit4 = neventmap.find(mit->first);
    std::map<TString, Double_t>::iterator mit5 = n_w_eventmap.find(mit->first);
    lumi_file <<  "SK" << mit2->second << "_dilep  " << mit4->second << " " << mit5->second << " " <<  mit3->second <<" "  << mit->second << " /data2/CatNtuples/" + string(version.Data()) +"/SKTrees/MCDiLep/" <<  mit2->second << "/" <<endl;
  }
  

  lumi_file << "" << endl;
  
  string lqdir = getenv("LQANALYZER_DIR");
  string lfile2 =  lqdir + "/LQRun/txt/datasets_snu_CAT_mc_private.txt";

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

float GetSumWeights(std::string filename){
  TFile* file = TFile::Open(filename.c_str());

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
  cout << "sum_weight = " << sum_weight << endl;
  return sum_weight;
}
