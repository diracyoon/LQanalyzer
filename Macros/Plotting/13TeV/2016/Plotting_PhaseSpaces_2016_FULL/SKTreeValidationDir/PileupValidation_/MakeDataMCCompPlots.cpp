//// Code makes directory of histograms and cutflow. 
#include "MakeDataMCCompPlots.h"
#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "CMS_lumi.h"


int main(int argc, char *argv[]) {
  
  /////////////////////////////////////////////////////
  //
  //  3 stages to setup : 
  //  1: McatNloWZ specifies with WZ MC we plot
  //  2: reg: specifies which fake region is plotted
  //  3: BM specifies if we plot 1.2 or 4.7 fb-1
  //
  //////////////////////////////////////////////////////
  
  TH1::SetDefaultSumw2(true);
  
  if(argc == 1) {
    cout << "No config file set" << endl;
    return 0;  
  }
  
  /// Set Plotting style
  setTDRStyle();
  gStyle->SetPalette(1);
  
  //read in config
  
  for (int i = 1; i < argc; ++i) {
    string configfile = argv[i];   
    SetUpMasterConfig(configfile);
    int a =MakeCutFlow_Plots(configfile);
  }
  
  system(("scp -r " + output_path + " jalmond@lxplus.cern.ch:~/www/SNU/WebPlots/").c_str());

  cout << "Open plots in " << output_index_path << endl;

  return 0;
}

int MakeCutFlow_Plots(string configfile){
  
  std::string pname = "/home/jalmond/WebPlots/"+ path + "/indexCMS.html";
  std::string phistname = "/home/jalmond/WebPlots/"+ path + "/histograms/" + histdir  + "/indexCMS.html";
  output_path = "/home/jalmond/WebPlots/"+ path ;
  output_index_path = string("https://jalmond.web.cern.ch/jalmond/SNU/WebPlots/")+ path + "/histograms/" + histdir  + "/indexCMS.html";
  system(("mkdir /home/jalmond/WebPlots/" + path).c_str());
  system(("mkdir /home/jalmond/WebPlots/" + path+ "/histograms/").c_str());
  system(("mkdir /home/jalmond/WebPlots/" + path+"/histograms/" + histdir + "/").c_str());

  histpage.open(phistname.c_str());
  page.open(pname.c_str());
  
  page << "<html><font face=\"Helvetica\"><head><title> HvyN Analysis </title></head>" << endl;
  page << "<body>" << endl;
  page << "<h1> HvyN Analysis Plots </h1>" << endl;
  page << "<br> <font size=\"4\"><b> " << message <<  " </b></font> <br><br>" << endl;
  page << "<a href=\"histograms/" +histdir + "/indexCMS.html\">"+ histdir + "</a><br>"; 
  
  //MakeCutFlow(histdir);  
  int M=MakePlots(histdir);  

  return 1;

}


int MakePlots(string hist) {


  
  ////////////////////// ////////////////
  ////  MAIN PART OF CODE for user/
  ///////////////////////////////////////
  //// What samples to use in histogram
  vector<pair<pair<vector<pair<TString,float> >, int >, TString > > samples;  
  vector<pair<pair<vector<pair<TString,float> >, int >, TString > > samples_ss;  
  vector<string> cut_label;
  //// Sets flags for using CF/NP/logY axis/plot data/ and which mc samples to use
  
  SetUpConfig( samples, samples_ss, cut_label);  
  cuts.clear();

  // ----------Get list of cuts to plot  ----------------------
  ifstream cut_name_file(cutfile.c_str());
  if(!cut_name_file) {
    cerr << "Did not find " + cutfile + ", exiting ..." << endl;
    return 1;
  }
  while(!cut_name_file.eof()) {
    string cutname;
    cut_name_file >> cutname;
    if(cutname=="END") break;
    //cutname = cutname + "POGTight";
    allcuts.push_back(cutname);
  }
  

  ifstream histo_name_file(histfile.c_str());
  if(!histo_name_file) {
    cerr << "Did not find " << histfile << ", exiting ..." << endl;
    return 1;
  }
  
  histpage << "<table border = 1><tr>"
	   << "<th> <a name=\"PlotName (variable_Cut)\"> PlotName (variable_Cut) </a> </th>"
	   << "<th> Data/MC Plot </th>"
	   << "<th> Data/MC LogPlots </th>"
	   << "</tr>" << endl;
  
  while(!histo_name_file.eof()) {
    string h_name;
    int rebin;
    double xmin,xmax;
    histo_name_file >> h_name;
    if(repeat(h_name))continue;
    if(h_name=="END") break;
    histo_name_file >> rebin;
    histo_name_file >> xmin;
    histo_name_file >> xmax;
    
    if(h_name.find("#")!=string::npos) continue;
    
    for(unsigned int ncut=0; ncut<allcuts.size();  ncut++){
      string name = allcuts.at(ncut) + "_" + h_name;
      if(TString(allcuts.at(ncut)).Contains("purw")) name =  h_name+ "_" + allcuts.at(ncut);
      bool isSS(false);
      //if(TString(allcuts.at(ncut)).Contains("SS"))isSS = true;

      cout << "###################### " << name << " ###########################"<< endl;
      
	
	/// Make nominal histogram stack
	map<TString, TH1*> legmap;
	THStack* mstack;
	if(!isSS) mstack= MakeStack(samples , "Nominal",name, xmin, xmax, legmap, rebin , true);
	else  mstack=MakeStack(samples_ss, "Nominal",name, xmin, xmax, legmap, rebin , true);

	THStack* mstack_nostat;
	if(!isSS)mstack_nostat = MakeStack(samples , "Nominal",name, xmin, xmax, legmap, rebin , false);
	else mstack_nostat = MakeStack(samples_ss , "Nominal",name, xmin, xmax, legmap, rebin , false);

	//// mhist sets error config
	map<TString,TH1*> mhist;
	mhist["Nominal"] = MakeSumHist(mstack);
	
	map<TString,TH1*> mhist_nostat;
	mhist_nostat["Nominal"] = MakeSumHist(mstack_nostat);

	TH1* hup = MakeStackUp(mhist, name+"UP");
	TH1* hup_nostat = MakeStackUp(mhist_nostat, name+"UPnostat");
	TH1* hdown = MakeStackDown(mhist, name+"DOWN");
	
	cout << "Final Background Integral = " <<  MakeSumHist(mstack)->Integral() << " : Up = " << hup->Integral() << " : Down= " << hdown->Integral() << endl;
	
	/// Make data histogram
	ylog=false;
	if(TString(name).Contains("llmass")){ylog=true;}
	if(TString(name).Contains("LeptonPt")){ylog=true;}
	
	if(TString(name).Contains("Tri")) {ylog=false;}
        //if(TString(name).Contains("SSE")) {ylog=false;}

	TH1* hdata = MakeDataHist(name, xmin, xmax, hup, ylog, rebin);
	CheckHist(hdata);	
	float ymin (0.001), ymax( 0.);
	
	ymax = GetMaximum(hdata, hup, ylog, name, xmax);
	
	if(showdata)cout << "Total data = " <<  hdata->Integral() << endl;
	scale = 1.;

	/// Make legend
	TLegend* legend = MakeLegend(legmap, hdata, showdata, ylog, ymax, xmax);       		

        vector<THStack*> vstack;		
	vstack.push_back(mstack);   	
	vstack.push_back(mstack_nostat);   	

	
	TCanvas* c = CompDataMC(hdata,vstack,hup,hdown, hup_nostat, legend,name,rebin,xmin,xmax, ymin,ymax, path, histdir,ylog, showdata, channel);      	

	string canvasname = c->GetName();
	canvasname.erase(0,4);

	
	PrintCanvas(c, histdir, canvasname, c->GetName());
    }
  }            
  page.close();
  
  return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///   CODE FOR CUTFLOW
///%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void MakeCutFlow(string type){

  
  vector<string> cut_label;  
  vector<pair<pair<vector<pair<TString,float> >, int >, TString > > cfsamples;  
  vector<pair<pair<vector<pair<TString,float> >, int >, TString > > cfsamples_ss;  
  SetUpConfig( cfsamples, cfsamples_ss, cut_label);

  cuts.clear();    
  // ----------Get list of cuts to plot  ----------------------
  cout << cutfile << endl;
  ifstream cut_name_file(cutfile.c_str());
  while(!cut_name_file.eof()) {
    string cutname;
    cut_name_file >> cutname;
    if(cutname=="END") break;
    //cutname = cutname+"POGTight";
    cut_label.push_back(cutname);
    cout << "$$$$$$$$$ " << hist << " $$$$$$$$$$$$$$$$$$$$" << endl;
    cout << cutname << endl;
    if(TString(cutname).Contains("purw")) return;
    else    cuts.push_back((cutname+ hist +"_" + cutname).c_str());    
  }

 
  vector<float> totalnumbers;
  vector<float> totalnumbersup;
  vector<float> totalnumbersdown;  
  
  int i_cut(0);
  for(vector<string>::iterator it = cuts.begin(); it!=cuts.end(); it++, i_cut++){
    bool isSS=false;
    //if(TString(*it).Contains("SS"))isSS = true;

    vector<pair<vector<pair<TString,float> >, TString> > samples;   
    if(!isSS){
      for(vector<pair<pair<vector<pair<TString,float> >, int >, TString > >::iterator it2 = cfsamples.begin(); it2!=cfsamples.end(); it2++){
	samples.push_back(make_pair(it2->first.first,it2->second));      
      }
    }
    else{
      for(vector<pair<pair<vector<pair<TString,float> >, int >, TString > >::iterator it2 = cfsamples_ss.begin(); it2!=cfsamples_ss.end(); it2++){
        samples.push_back(make_pair(it2->first.first,it2->second));
      }
    }
    
      
    /// Vectors for cutflow	table
    map<TString,float> samples_numbers;
    map<TString,float> samples_numbers_staterr;
    map<TString,float> samples_numbers_up;
    map<TString,float> samples_numbers_down;
    
    /// Vector for systematic table
    map<TString,float> syst_stat;
    map<TString,float> syst_norm;
    map<TString,float> syst_total;

    
    float totalbkg(0.),totalbkgdown(0.), totalbkgup(0.); 
    for(vector<pair<vector<pair<TString,float> >, TString> >::iterator it2 = samples.begin() ; it2!= samples.end(); it2++){
      TString cutname = *it;
      totalbkg+= Calculate(*it,"Normal",*it2);
      totalbkgup+= Calculate(*it,"Up",*it2);
      totalbkgdown+= Calculate(*it,"Down",*it2);
      
      samples_numbers[it2->second] = Calculate(cutname,"Normal",*it2);
      samples_numbers_up[it2->second] = Calculate(cutname,"Down",*it2);
      samples_numbers_down[it2->second] = Calculate(cutname,"Up",*it2);
      samples_numbers_staterr[it2->second] = Calculate(cutname,"StatErr",*it2);
    }	
    
    
    float totaldata(0.);
    if(showdata) totaldata = GetIntegral(*it,"data","data");
    float errdata(0.);
    if(showdata) errdata= GetError(*it,"data","data");
    
    
    float totalerr_up(0.),totalerr_down(0.),totalerrup(0.),totalerrdown(0.),total_staterr(0.);
    for(map<TString,float>::iterator mapit = samples_numbers.begin(); mapit!= samples_numbers.end(); mapit++){
      
      map<TString,float>::iterator mapit_up;
      mapit_up = samples_numbers_up.find(mapit->first);
      map<TString,float>::iterator mapit_down;
      mapit_down = samples_numbers_down.find(mapit->first);
      map<TString,float>::iterator mapit_stat;
      mapit_stat = samples_numbers_staterr.find(mapit->first);
      totalerr_up += ( mapit_up->second*mapit_up->second+ mapit_stat->second*mapit_stat->second); 
      totalerr_down += (mapit_down->second*mapit_down->second + mapit_stat->second*mapit_stat->second); 
      totalerrup += (mapit_up->second*mapit_up->second); 
      totalerrdown += (mapit_down->second*mapit_down->second); 
      total_staterr += mapit_stat->second*mapit_stat->second;
      TString sample = mapit->first;
      if(sample.Contains("t#bar{t},t/#bar{t},t/#bar{t}W,t#bar{t}V")) sample = "t$\\bar{t}$,t/$\\bar{t}$,t/$\\bar{t}$W,t$\\bar{t}$V";
      else if(sample.Contains("t#bar{t}")) sample = "t$\bar{t}$";
      else if(sample.Contains("t#bar{t}V")) sample = "t$\bar{t}$+V";
      else if(sample.Contains("t/#bar{t}")) sample = "t/$\bar{t}$";
      
      if(sample.Contains("DY#rightarrow ll; 10 < m(ll) < 50")) sample = "DY$\\rightarrow$ ll; 10 < m(ll) < 50";
      if(sample.Contains("DY#rightarrow ll; m(ll) > 50")) sample = "DY$\\rightarrow$ ll; m(ll) > 50";
      if(sample.Contains("DY#rightarrow ll")) sample = "DY$\\rightarrow$ ll";
      cout << sample << " background = " << mapit->second<< " +- " << mapit_stat->second << " + " << mapit_up->second << " - " << mapit_down->second <<  endl;      
   
      
    }
  
    
    totalerr_up = sqrt(totalerr_up);
    totalerr_down = sqrt(totalerr_down);
    total_staterr = sqrt(total_staterr);
    totalerrup = sqrt(totalerrup);
    totalerrdown = sqrt(totalerrdown);
    
    cout << "Total Bkg   = " << totalbkg << "+- " << total_staterr << " + " << totalerrup << " - " << totalerrdown << endl;
    if(showdata)cout << "Total Data  = " << totaldata << endl;
    cout << "-------------" << endl;
    if(totaldata > totalbkg)cout <<"Significance = " << (totaldata - totalbkg) / (sqrt( (errdata*errdata) + (totalerr_up*totalerr_up))) << endl;
    else cout <<"Significance = " << (totaldata - totalbkg) / (sqrt( (errdata*errdata) + (totalerr_down*totalerr_down))) << endl;
    float significance = (totaldata - totalbkg) / (sqrt( (errdata*errdata) + (totalerr_up*totalerr_up))) ;
    
    if(significance < 0.) significance = (totaldata - totalbkg) / (sqrt( (errdata*errdata) + (totalerr_down*totalerr_down))) ;
    

    //// Make TEX file
    ofstream ofile_tex;
    string latex_file =  "Tables/" + cut_label.at(i_cut) + ".tex";
    ofile_tex.open(latex_file.c_str());
    ofile_tex.setf(ios::fixed,ios::floatfield);
    ofile_tex << "\\documentclass[10pt]{article}" << endl;
    ofile_tex << "\\usepackage{epsfig,subfigure,setspace,xtab,xcolor,array,colortbl}" << endl;

    ofile_tex << "\\begin{document}" << endl;
    ofile_tex << "\\input{Tables/" + cut_label.at(i_cut)  + "Table.txt}" << endl;
    ofile_tex << "\\end{document}" << endl;
    
    /// Make text file
    ofstream ofile;
    string latex =  "Tables/" + cut_label.at(i_cut) + "Table.txt";
    
    ofile.open(latex.c_str());
    ofile.setf(ios::fixed,ios::floatfield); 
    
    ofile.precision(1);
    ofile << "\\begin{table}[h]" << endl;
    ofile << "\\begin{center}" << endl;
    ofile << "\\begin{tabular}{lr@{\\hspace{0.5mm}}c@{\\hspace{0.5mm}}c@{\\hspace{0.5mm}}l}" << endl;
    ofile << "\\hline" << endl;
    ofile << "\\hline" << endl;   
    //ofile << "Source & \\multicolumn{4}{c}{$\\mu\\mu\\mu$} \\"<<"\\" << endl;
    ofile << "Source & \\multicolumn{4}{c}{" << columnname << "} \\"<<"\\" << endl;
    ofile << "\\hline" << endl;   
	 
    
    for(map<TString,float>::iterator mapit = samples_numbers.begin(); mapit!= samples_numbers.end(); mapit++){
      
      map<TString,float>::iterator mapit_up;
      mapit_up = samples_numbers_up.find(mapit->first);
      map<TString,float>::iterator mapit_down;
      mapit_down = samples_numbers_down.find(mapit->first);
      map<TString,float>::iterator mapit_stat;
      mapit_stat = samples_numbers_staterr.find(mapit->first);
      
      TString sample = mapit->first;

      if(sample.Contains("t#bar{t},t/#bar{t},t/#bar{t}W,t#bar{t}V")) sample = "t$\\bar{t}$,t/$\\bar{t}$,t/$\\bar{t}$W,t$\\bar{t}$V";
      else if(sample.Contains("t/#bar{t}")) sample = "t/$\\bar{t}$";
      else if(sample.Contains("t#bar{t}V")) sample = "t$\\bar{t}$+V";
      else if(sample.Contains("t#bar{t}")) sample = "t$\\bar{t}$";

      if(sample.Contains("DY#rightarrow ll; 10 < m(ll) < 50")) sample = "DY$\\rightarrow$ ll; 10 < m(ll) < 50";
      if(sample.Contains("DY#rightarrow ll; m(ll) > 50")) sample = "DY$\\rightarrow$ ll; m(ll) > 50";
      if(sample.Contains("DY#rightarrow ll")) sample = "DY$\\rightarrow$ ll";

      if(mapit->second!=0.0){
	ofile << sample + "&" <<  mapit->second << "& $\\pm$& "  << mapit_stat->second <<  "&$^{+" <<  mapit_up->second << "}_{-" <<  mapit_down->second  << "}$" ; 
	ofile  <<  "\\"  << "\\" << endl;	   
      }
    }
    ofile << "\\hline" << endl;
    ofile << "Total&" << totalbkg << "& $\\pm$&"  << total_staterr << "&$^{+" << totalerrup  << "}_{-" << totalerrdown << "}$" ; 
    ofile  <<  "\\"  << "\\" << endl;
    ofile << "\\hline" << endl;
    
    ofile << "Data&  \\multicolumn{4}{c}{$" << totaldata << "$}\\" << "\\" <<endl;
    ofile << "\\hline" << endl;
    if(significance < 0) ofile << "Signficance&  \\multicolumn{4}{c}{$" << significance << "\\sigma$}\\" << "\\" <<endl;
    if(significance > 0) ofile << "Signficance&  \\multicolumn{4}{c}{$+" << significance << "\\sigma$}\\" << "\\" <<endl;
    ofile << "\\hline" << endl;   
    ofile << "\\hline" << endl;   
    ofile << "\\end{tabular}" << endl;
    ofile << "\\caption{" << caption << "}" << endl;
    ofile << "\\end{center}" << endl;
    ofile << "\\end{table}" << endl;    
     
  
    string latex_command = "latex Tables/" + cut_label.at(i_cut) +".tex";
    string dvi_command = "dvipdf " + cut_label.at(i_cut) +".dvi";
    string mv_command = "mv " + cut_label.at(i_cut) +".pdf /home/jalmond/WebPlots/" + path +"/histograms/"+ histdir ;
    
    system((latex_command.c_str()));
    system((dvi_command.c_str()));
    system((mv_command.c_str()));
    system(("rm *aux"));
    system(("rm *log"));
    system(("rm *dvi"));
    
    string cftitle = cut_label.at(i_cut);
    
    histpage << "<tr><td>"<< "cutflow " + cut_label.at(i_cut)  <<"</td>"<<endl;
    histpage <<"<td>"<<endl;
    histpage << "<a href=\"" << cut_label.at(i_cut)  << ".pdf\">";
    histpage << "Cutflow: " + cut_label.at(i_cut)  + ".pdf</p>";
    histpage << "</td>" << endl;
  }

  return;
}

bool repeat (string hname){
  map<string,int>::iterator mit = norepeatplot.find(hname);
  if(mit!=norepeatplot.end())return true;
  else{
    norepeatplot[hname]=1;
    return false;
  }
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PrintCanvas(TCanvas* c1, string folder, string plot_description, string title){

  std::string tpdf = "/home/jalmond/WebPlots/"+ path +  "/histograms/"+folder+"/"+title;
  
  if(plot_description.empty())plot_description=title;
  histpage << "<tr><td>"<< plot_description <<"</td>"<<endl;
  histpage <<"<td>"<<endl;
  histpage << "<a href=\"" << title.c_str() << ".png\">";
  histpage << "<img src=\"" << title.c_str() << ".png\" width=\"100%\"/>";
  histpage << "</td>" << endl;
  histpage <<"<td>"<<endl;
  histpage << "<a href=\"" << title.c_str() << "_log.png\">";
  histpage << "<img src=\"" << title.c_str() << "_log.png\" width=\"100%\"/>";
  histpage << "</td>" << endl;

  
  return;
  
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


TLegend* MakeLegend(map<TString, TH1*> map_legend,TH1* hlegdata,  bool rundata , bool logy, float ymax, float xmax){
  
  double x1 = 0.5;
  double y1 = 0.5;
  double x2 = 0.6;
  double y2 = 0.9;

  int nbinsX=hlegdata->GetNbinsX();
  
  int max_bin = 0;
  int ymax_bin = 0;
  float r_max_tmp = 0.;
  for(unsigned int ibin=1; ibin < hlegdata->GetNbinsX()+1; ibin++){
    float r_max = hlegdata->GetBinContent(ibin) / hlegdata->GetMaximum();
    if(r_max > r_max_tmp){
      r_max_tmp= r_max;
      ymax_bin= ibin;
    }
    if( xmax >= hlegdata->GetBinLowEdge(ibin) && xmax <  hlegdata->GetBinLowEdge(ibin+1)) max_bin = ibin;
  }
  if(max_bin == 0) max_bin = hlegdata->GetNbinsX()+1;

  if(ymax_bin/max_bin > 0.5){
    x1 = 0.2;
    y1 = 0.6;
    x2 = 0.5;
    y2 = 0.9;
  }
  else{
    x1 = 0.6;
    y1 = 0.7;
    x2 = 0.95;
    y2 = 0.9;
  }
  
  cout << "Test" << endl;
  TLegend* legendH = new TLegend(x1,y1,x2,y2);
  legendH->SetFillColor(kWhite);
  legendH->SetTextFont(42);
  
  legendH->SetBorderSize(0);
  legendH->SetTextSize(0.025);
  //legendH->SetNColumns(2);

  if(rundata) 	legendH->AddEntry(hlegdata,"Data","pE");
  
  //  for(map<TString, TH1*>::iterator it = map_legend.begin(); it!= map_legend.end(); it++){
  
  vector<TString> legorder;
  legorder.push_back("DY#rightarrow ll");
  legorder.push_back("Wjets");
  //legorder.push_back("t#bar{t},t/#bar{t},t/#bar{t}W,t#bar{t}V");
  legorder.push_back("Top");
  legorder.push_back("VV");
  //  legorder.push_back("Diboson");
  //legorder.push_back("QCD");
  //  legorder.push_back("Other");
  //legorder.push_back("Mismeas. Charge Background");
  //legorder.push_back("Misid. Lepton Background");
  
  map<double, TString> order_hists;
  for(map<TString, TH1*>::iterator it = map_legend.begin(); it!= map_legend.end(); it++){
    order_hists[it->second->Integral()] = it->first;
  }
  
  
  if(map_legend.size()  < 3){
    for(map<TString, TH1*>::iterator it = map_legend.begin(); it!= map_legend.end(); it++) {
      legendH->AddEntry(it->second, it->first,"f");
    }
  }else{
    //for(unsigned int ileg = 0; ileg < legorder.size() ; ileg++){
    //map<TString, TH1*>::iterator it = map_legend.find(legorder.at(ileg));
    //if(it->second)legendH->AddEntry(it->second,it->first.Data(),"f");    
      
    // }
    for(map<double, TString>::iterator it =order_hists.begin(); it!= order_hists.end(); it++) {
      map<TString, TH1*>::iterator it2 = map_legend.find(it->second);
      if(it->first > 0)legendH->AddEntry(it2->second,it->second.Data(),"f");
    }
    
  }
  legendH->SetFillColor(kWhite);
  legendH->SetTextFont(42);

  return legendH;
  
}




TH1* MakeDataHist(string name, double xmin, double xmax, TH1* hup, bool ylog, int rebin){

  /// Make data histogram
  TFile* file_data =  TFile::Open((dataloc).c_str());
  TH1* hdata = dynamic_cast<TH1*> ((file_data->Get(name.c_str()))->Clone());
  
  hdata->Rebin(rebin);

  float ymin (0.2), ymax( 10000.);
  ymax = GetMaximum(hdata, hup, ylog, name, xmax);
  

  /// Set Ranges / overflows

  FixOverUnderFlows(hdata, xmax);  
   
  hdata->GetXaxis()->SetRangeUser(xmin,xmax);
  hdata->GetYaxis()->SetRangeUser(ymin, ymax);

  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(1.6);
  
  //// X title  
  SetTitles(hdata, name);
  
  return hdata;

}

void CheckHist(TH1* h){
  if(!h) {
    cout << "Not able to find data histogram" << endl;
  }
}


vector<pair<TString,float> >  InitSample (TString sample){
  
  vector<pair<TString,float> > list;  
  
  if(sample.Contains("dylow_")){
    list.push_back(make_pair("DYJets_10to50",0.2));
  }

  if(sample.Contains("dyhigh_")){
    list.push_back(make_pair("DYJets",0.2));    
  }
  if(sample.Contains("dy_")){
    //list.push_back(make_pair("DY10to50_MCatNLO",0.2));
    //list.push_back(make_pair("DY50plus_MCatNLO",0.2));
    list.push_back(make_pair("DYJets_MG_10to50",0.2));
    list.push_back(make_pair("DYJets",0.2));
  }



  if(sample.Contains("wjet_")){
    list.push_back(make_pair("WJets",0.2));
  }      
  ///// Top samples //////////////////    
  if(sample.Contains("top")){
    // https://github.com/ekenn003/PlotMaker/blob/master/python/xsections.py
    //list.push_back(make_pair("SingleTop_tW",0.25));
    //list.push_back(make_pair("SingleTop_s",0.25));
    //list.push_back(make_pair("SingleTop_t",0.25));
    //list.push_back(make_pair("SingleTbar_t",0.25));
    //list.push_back(make_pair("SingleTbar_tW",0.25));
    list.push_back(make_pair("TT_powheg",0.25));
  }

  
  if(sample.Contains("qcd"))
    {
      list.push_back(make_pair("QCD_Pt-20to30_EMEnriched",0.30));
      list.push_back(make_pair("QCD_Pt-300toInf_EMEnriched",0.30));
      list.push_back(make_pair("QCD_Pt-30to50_EMEnriched",0.30));
      list.push_back(make_pair("QCD_Pt-170to300_EMEnriched",0.30));
      list.push_back(make_pair("QCD_Pt-50to80_EMEnriched",0.30));
      list.push_back(make_pair("QCD_Pt-120to170_EMEnriched",0.30));
      list.push_back(make_pair("QCD_Pt-80to120_EMEnriched",0.30));
    }


  //////// Diboson ////////
  if(sample.Contains("vv")){    
    list.push_back(make_pair("WZTo3LNu_powheg",0.15));
    list.push_back(make_pair("WJets",0.2));
    
    //list.push_back(make_pair("WWTo2L2Nu",0.15));
    //list.push_back(make_pair("WWToLNuQQ",0.15));
    //list.push_back(make_pair("ZZTo4L_powheg",0.15));
    list.push_back(make_pair("ZZ",0.15));
    list.push_back(make_pair("WW",0.15));
    
  }
  if(sample.Contains("other_")){
    list.push_back(make_pair("WpWpQCD",0.15));
    list.push_back(make_pair("WpWpEWK",0.15));
    list.push_back(make_pair("WWW",0.15));
    list.push_back(make_pair("WGtoLNuG",0.15));
    list.push_back(make_pair("ZGto2LG",0.15));
    list.push_back(make_pair("ZZZ",0.15));
    list.push_back(make_pair("ttH_nonbb",0.15));
    list.push_back(make_pair("ttH_bb",0.15));
  }
  if(sample.Contains("prompt")){
    list.push_back(make_pair("WWW",0.25));
    list.push_back(make_pair("TTWW",0.25));
    list.push_back(make_pair("TTG",0.25));
    list.push_back(make_pair("ZZZ",0.25));
    list.push_back(make_pair("WWG",0.25));
    list.push_back(make_pair("HtoWW",0.25));
    list.push_back(make_pair("HtoTauTau",0.22));
    list.push_back(make_pair("ggHtoZZ",0.22));
    list.push_back(make_pair("WZ_py",0.12));
    list.push_back(make_pair("ZZ_py",0.09));
    list.push_back(make_pair("WW_dp",0.5));
    list.push_back(make_pair("ttW",0.25));
    list.push_back(make_pair("ttZ",0.25));
  }

  if(sample.Contains("nonprompt_DoubleEG")){
    list.push_back(make_pair("nonprompt_DoubleEG",0.34));
  }

  if(sample.Contains("chargeflip")){
    list.push_back(make_pair("chargeflip",0.12));
  }

  if(list.size()==0) cout << "Error in making lists " << sample<< endl;
 
  return list;
}

void CheckSamples(int nsamples){
  if(nsamples==0) {
    cout << "No sample in this vector" << endl;
    exit(1);
  }
  return;
}




THStack* MakeStack(vector<pair<pair<vector<pair<TString,float> >, int >, TString > > sample, TString type, string name, float xmin, float xmax,map<TString, TH1*>& legmap , int rebin, bool include_syst_err){
  

  string clonename = name;	
  if(include_syst_err) clonename += "_andsysyerr";
  
  THStack* stack = new THStack(clonename.c_str(), clonename.c_str());
  
  bool debug(false);
  if(type.Contains("Nominal")) debug=true;

  TString fileloc = "";
 
  TDirectory* origDir = gDirectory;

  float sum_integral=0.;

  for(vector<pair<pair<vector<pair<TString,float> >, int >, TString > >::iterator it = sample.begin() ; it!= sample.end(); it++){
    if(type.Contains("Nominal")) fileloc = mcloc;
        
    if(!type.Contains("Nominal")) {
      if(it->first.first.at(0).first.Contains("nonprompt_DoubleEG"))fileloc=mcloc;
    }    
    
    CheckSamples( it->first.first.size() );
    
    int isample=0;
    TFile* file =  TFile::Open((fileloc+ fileprefix + it->first.first.at(isample).first + filepostfix).Data());
    cout << fileloc+ fileprefix + it->first.first.at(isample).first + filepostfix << endl;
    if(!file) cout << "Could not open " << fileloc+ fileprefix + it->first.first.at(isample).first + filepostfix << endl;
    cout << "file = " << file << endl;
    
    gROOT->cd();
    TDirectory* tempDir = 0;
    int counter = 0;
    while (not tempDir) {
      std::stringstream dirname;
      dirname << "WRHNCommonLeptonFakes_%i" << counter;
      if (gROOT->GetDirectory((dirname.str()).c_str())) {
	++counter;
	continue;
      }      
      tempDir = gROOT->mkdir((dirname.str()).c_str());      
    }
            
    tempDir->cd();
    
    cout <<  clonename.c_str() << endl;
    TH1* h_tmp = dynamic_cast<TH1*> ((file->Get(name.c_str()))->Clone(clonename.c_str()));
  cout << h_tmp << endl;
  while(!h_tmp) {
      isample++;
      gROOT->cd();
      file =  TFile::Open((fileloc+ fileprefix + it->first.first.at(isample).first + filepostfix).Data());
      tempDir->cd();
      h_tmp = dynamic_cast<TH1*> ((file->Get(name.c_str()))->Clone(clonename.c_str()));
    }
    
    CheckHist(h_tmp);

    if(debug)cout <<  it->second <<  "  contribution " << 1 << "/" << it->first.first.size()  << " is from " << fileprefix + it->first.first.at(isample).first + filepostfix <<" : Integral = " <<h_tmp->Integral() << " " << fileloc << endl;
    
    
    for(unsigned int i=isample+1; i < it->first.first.size(); i++){	    
      clonename+="A";
             
      origDir->cd();
      
      TFile* file_loop =  TFile::Open((fileloc+ fileprefix + it->first.first.at(i).first + filepostfix).Data());	    

      tempDir->cd();
      TH1* h_loop = dynamic_cast<TH1*> ((file_loop->Get(name.c_str()))->Clone(clonename.c_str()));	    	    	    
      if(!h_loop) continue;
      CheckHist(h_loop);

      h_tmp->Add(h_loop);	  	    	    
            
      if(debug)cout <<  it->second <<  "  contribution " <<i+1 <<"/" << it->first.first.size()  << " is from ExampleAnalyzer_SK" << it->first.first.at(i).first << ".NTUP_SMWZ.Reco.root : Integral = " <<h_loop->Integral() << " sum integral = " << h_tmp->Integral()    << endl;
      file_loop->Close();
    }	  	  
    
    /// TH1* is now made. Now make pretty
    FixOverUnderFlows(h_tmp, xmax);	  
    ///Set colors
    h_tmp->SetFillColor(it->first.second);
    //h_tmp->SetLineColor(it->first.second);	  
    h_tmp->SetLineColor(kBlack);
    
    if(!stack->GetHists()) {
      stack->SetName( (string("s_") + name).c_str() );
      stack->SetTitle( (string("s_") + name).c_str() );
      SetTitles(h_tmp, name);
  
    }//stack empt   
    
    h_tmp->Rebin(rebin);
    h_tmp->SetLineWidth(3.0);
    
    SetErrors(h_tmp, it->first.first.at(0).second, include_syst_err );


    stack->Add(h_tmp);
    sum_integral+=h_tmp->Integral();
    
    if(type.Contains("Nominal")) {
      legmap[it->second] = h_tmp;
    }
    
    file->Close();	  
    origDir->cd();
  }
  
  
  return stack;
}


TH1* MakeStackUp(map<TString, TH1*> map_of_stacks, TString clonename){
  
  map<TString, TH1*>::iterator it =  map_of_stacks.find("Nominal");
  
  float norm=1;
  TH1* h_up = dynamic_cast<TH1*>(it->second->Clone(clonename.Data())); // copy of nominal
  
  for(int binx=1 ; binx < h_up->GetNbinsX()+1; binx++){
    float nom_content = h_up->GetBinContent(binx);
    float nom_error = h_up->GetBinError(binx);
    
    /// Now use 5%+ norm + stat error
    float errup2 =  nom_error*nom_error + 0.05*0.05*nom_content*nom_content;
    if(clonename.Contains("stat")) errup2 =  nom_error*nom_error;

    /// add rest of systs

    float new_bin = nom_content + sqrt(errup2);
    h_up->SetBinContent(binx,new_bin);
   
  }
  
  return  h_up;
  
}


TH1* MakeStackDown(map<TString, TH1*> map_of_stacks, TString clonename){
  
  map<TString, TH1*>::iterator it =  map_of_stacks.find("Nominal");
  
  float norm=1;
  TH1* h_down = dynamic_cast<TH1*>(it->second->Clone(clonename.Data())); // copy of nominal
  
  for(int binx=1 ; binx < h_down->GetNbinsX()+1; binx++){
    float nom_content = h_down->GetBinContent(binx);
    float nom_error = h_down->GetBinError(binx);
    
    
    //// nom_error = stat err + normalisation error, set previously on nom hist
    float errdown2 =  nom_error*nom_error;
    
    /// Now use 5%+ norm + stat error
    float errup2 =  nom_error*nom_error + 0.05*0.05*nom_content*nom_content;
    
    float new_bin = nom_content  - sqrt(errdown2);
    h_down->SetBinContent(binx,new_bin);
   
  }
  
  return  h_down;


}


TH1* MakeSumHist(THStack* thestack){
  
  TH1* hsum=0;  
  TList* list = thestack->GetHists();
  TIter it(list, true);
  TObject* obj=0;
  while( (obj = it.Next()) ) {
    TH1* h = dynamic_cast<TH1*>(obj);
    
    if(!hsum) hsum = (TH1*)h->Clone( (string(h->GetName()) + "_sum").c_str() );
    else {
      hsum->Add(h, 1.0);
    }
  }//hist loop
  
  return hsum;
}


void SetErrors(TH1* hist, float normerr, bool includestaterr ){

  for(int binx =1; binx < hist->GetNbinsX()+1; binx++){
    float newbinerr = hist->GetBinError(binx)*hist->GetBinError(binx) + hist->GetBinContent(binx)*hist->GetBinContent(binx)*normerr*normerr;
    if(!includestaterr)  newbinerr =hist->GetBinError(binx)*hist->GetBinError(binx) ;
  
    hist->SetBinError(binx, sqrt(newbinerr));
  }
  
  return;

}



void SetTitles(TH1* hist, string name){
  
  string xtitle ="";
  string ytitle ="Counts";

  float binedge_up = hist->GetBinLowEdge(2);
  float binedge_down = hist->GetBinLowEdge(1);
  
  float width = binedge_up - binedge_down;
  
  std::ostringstream str_width;
  str_width<< int(width);
  
  cout << "Hist " << name << " Add GeV to y title = " << HistInGev(name) << endl;
  if(HistInGev(name)) ytitle = "Counts / " +str_width.str() + " GeV";
  
  if(name.find("h_MET")!=string::npos){
    xtitle="E^{miss}_{T} (GeV)"; 
    if(name.find("phi")!=string::npos){
      xtitle="#phi_{E^{miss}_{T}} "; 
    }
  }
  if(name.find("PFMet")!=string::npos){
    xtitle="E^{miss}_{T} (GeV)";
  }
  if(name.find("h_MT")!=string::npos) xtitle="M_{T} (GeV)";
  if(name.find("h_dphi_METe")!=string::npos) xtitle="#Delta (#phi_{E^{miss}_{T}} - #phi_{el})";
  if(name.find("h_dphi_METm")!=string::npos) xtitle="#Delta (#phi_{E^{miss}_{T}} - #phi_{mu})";
  if(name.find("h_NoHFMET")!=string::npos) xtitle="NoHF E^{miss}_{T} (GeV)";
  if(name.find("h_PFMET")!=string::npos) xtitle="E^{miss}_{T} (GeV)";

  if(name.find("h_jet_emfra")!=string::npos) xtitle="Jet EMFrac";
  if(name.find("jet_el_ptratio")!=string::npos) xtitle="El p_{T}/ Jet p_{T}";

    
  if(name.find("muons_eta")!=string::npos)xtitle="Muon #eta";
  if(name.find("muons_phi")!=string::npos)xtitle="Muon #phi";
  if(name.find("MuonPt")!=string::npos)xtitle="Muon p_{T} (GeV)";
  if(name.find("MuonD0")!=string::npos)xtitle="d0";
  if(name.find("MuonD0Sig")!=string::npos)xtitle="d0/#Sigma_{d0}";
  if(name.find("LeptonPt")!=string::npos)xtitle="Lepton p_{T} (GeV/c)";
  if(name.find("leadingLeptonPt")!=string::npos)xtitle="Leading lepton p_{T} (GeV/c)";
  if(name.find("secondLeptonPt")!=string::npos)xtitle="Second lepton p_{T} (GeV/c)";
  if(name.find("leadingMuonPt")!=string::npos)xtitle="Lead p_{T} (GeV)";
  if(name.find("secondMuonPt")!=string::npos)xtitle="Second p_{T} (GeV)";
  if(name.find("thirdMuonPt")!=string::npos)xtitle="Third p_{T} (GeV)";
  
  if(name.find("jets_pt")!=string::npos)xtitle="Jet p_{T} (GeV)";
  if(name.find("LeptonEta")!=string::npos)xtitle="Lepton  #eta";
  
  if(name.find("electrons_eta")!=string::npos)xtitle="Electron #eta";
  if(name.find("electrons_phi")!=string::npos)xtitle="Electron #phi";
  if(name.find("el_pt")!=string::npos)xtitle="Electron p_{T} (GeV)";
  if(name.find("leadingElectronPt")!=string::npos)xtitle="Leading electron p_{T} (GeV/c)";
  if(name.find("secondElectronPt")!=string::npos)xtitle="Trailing electron p_{T} (GeV/c)";
  if(name.find("thirdELectronPt")!=string::npos)xtitle="Third electron p_{T} (GeV)";
  
  if(name.find("lljjmass")!=string::npos)xtitle="e^{#pm}#mu^{#pm}jj invariant mass (GeV/c^{2})";
  if(name.find("llmass")!=string::npos)xtitle="ll invariant mass (GeV/c^{2})";
  if(name.find("l1jjmass")!=string::npos)xtitle="l_{1}jj invariant mass (GeV/c^{2})";
  if(name.find("l2jjmass")!=string::npos)xtitle="l_{2}jj invariant mass (GeV/c^{2})";

  if(name.find("charge")!=string::npos)xtitle="sum of lepton charge";

  if(name.find("mumumass")!=string::npos)xtitle="m(#mu#mu) (GeV)";
  if(name.find("eemass")!=string::npos)xtitle="e^{#pm}e^{#pm} invariant mass (GeV)";
  if(name.find("emumass")!=string::npos)xtitle="e^{#pm}#mu^{#pm} invariant mass (GeV)";
  
  if(name.find("jets_eta")!=string::npos)xtitle="jet #eta";
  if(name.find("jets_phi")!=string::npos)xtitle="jet #phi";
  if(name.find("Njets")!=string::npos)xtitle="Number of jets";
  if(name.find("Nbjet")!=string::npos)xtitle="Number of bjets";
  if(name.find("bTag")!=string::npos)xtitle="CSV";

  if(name.find("el1jet_mindr")!=string::npos)xtitle="min#Delta R(e_{1}j)";
  if(name.find("el2jet_mindr")!=string::npos)xtitle="min#Delta R(e_{2}j)";

  if(name.find("leadMuonJetdR")!=string::npos)xtitle="min#Delta R(#mu j)";
  if(name.find("leadJetdR")!=string::npos)xtitle="min#Delta R(jj)";
  if(name.find("mu1jjmass")!=string::npos)xtitle="m(#mu_{1}jj) (GeV)";
  if(name.find("mu2jjmass")!=string::npos)xtitle="m(#mu_{2}jj) (GeV)";
  if(name.find("mumujjmass")!=string::npos)xtitle="m(#mu#mujj) (GeV)";

  if(name.find("leadElectronJetdR")!=string::npos)xtitle="min#Delta R(e_j)";
  if(name.find("e1jjmass")!=string::npos)xtitle="e_{1}jj invariant mass (GeV/c^{2})";
  if(name.find("e2jjmass")!=string::npos)xtitle="e_{2}jj invariant mass (GeV/c^{2})";
  if(name.find("eejjmass")!=string::npos)xtitle="e^{#pm}e^{#pm}jj invariant mass (GeV/c^{2})";

  if(name.find("leadingMuonIso")!=string::npos)xtitle="PF Iso #mu_{1} (GeV)";
  if(name.find("secondMuonIso")!=string::npos)xtitle="PF Iso #mu_{2} (GeV)";

  if(name.find("leadingElectronIso")!=string::npos)xtitle="PF Iso e_{1} (GeV)";
  if(name.find("secondELectronIso")!=string::npos)xtitle="PF Iso e_{2} (GeV)";

  if(name.find("MuonD0_")!=string::npos)xtitle="d0";
  if(name.find("MuonD0Sig")!=string::npos)xtitle="d0sig";
  
  if(name.find("D0_")!=string::npos)xtitle="d0";
  if(name.find("D0Sig")!=string::npos)xtitle="d0sig";

  if(name.find("nVertice")!=string::npos)xtitle="Number of vertices";

  if(name.find("MuonJetdR")!=string::npos)xtitle="mindR(#mu,jet)";
  if(name.find("ElectronJetdR")!=string::npos)xtitle="mindR(e,jet)";
  if(name.find("LeadJetdR")!=string::npos)xtitle="mindR(jet,jet)";
  if(name.find("LeadMuondR")!=string::npos)xtitle="mindR(#mu,#mu)";
  if(name.find("LeadElectrondR")!=string::npos)xtitle="mindR(e,e)";
  

  if(name.find("muon_deta_")!=string::npos)xtitle="#Delta #eta (#mu,#mu)";
  if(name.find("el_deta_")!=string::npos)xtitle="#Delta #eta (e,e)";
  if(name.find("leaddimudeltaR_")!=string::npos)xtitle="#Delta R (#mu,#mu)";
  if(name.find("leaddieldeltaR_")!=string::npos)xtitle="#Delta R (e,e)";

  if(name.find("dijetsmass")!=string::npos)xtitle="m(j_{1}j_{2}) (GeV/c^{2})";
  if(name.find("leaddijetdr")!=string::npos)xtitle="#Delta R(j_{1}j_{2})";
  if(name.find("leadingJetPt")!=string::npos)xtitle="jet1 p_{T} (GeV)";
  if(name.find("secondJetPt")!=string::npos)xtitle="jet2 p_{T} (GeV)";
  if(name.find("_Nvtx")!=string::npos)xtitle="# Reco Vertices";


  hist->GetXaxis()->SetTitle(xtitle.c_str());
  hist->GetYaxis()->SetTitle(ytitle.c_str());

  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  return;
}


bool HistInGev(string name){
  
  bool ingev=false;
  cout << name << " : ingev =" << ingev << endl;
  if(name.find("ElectronPt")!=string::npos){ ingev=true;cout << "Found ElectronPt" << endl;}
  if(name.find("_pt_")!=string::npos){ ingev=true;cout << "Found _pt_" << endl;}
  if(name.find("pt")!=string::npos){
    if(name.find("LeptonPt")!=string::npos) {ingev=true;cout << "Found pt" << endl;}
    else if(name.find("!Lepton")!=string::npos){ ingev=true;cout << "Found pt" << endl;}
  }
  if(name.find("mass_")!=string::npos){ ingev=true;cout << "Found mass" << endl;}
  if(name.find("MET")!=string::npos){ ingev=true;cout << "Found MET" << endl;}
  if(name.find("Nvtx")!=string::npos) { ingev=false;}
  cout << name << " : ingev =" << ingev<< endl;

  return ingev;

}


float  GetMaximum(TH1* h_data, TH1* h_up, bool ylog, string name, float xmax){
  
  float yscale= 1.2;
  
  int max_bin = 0;
  int ymax_bin = 0;
  float r_max_tmp = 0.;
  
  for(unsigned int ibin=1; ibin < h_data->GetNbinsX()+1; ibin++){
    float r_max = h_data->GetBinContent(ibin) / h_data->GetMaximum();
    if(r_max > r_max_tmp){
      r_max_tmp= r_max;
      ymax_bin= ibin;
    }
    if( xmax >= h_data->GetBinLowEdge(ibin) && xmax <  h_data->GetBinLowEdge(ibin+1)) max_bin = ibin;
  }
  if(max_bin == 0) max_bin = h_data->GetNbinsX()+1;
  
  bool scale_up=false;
  if(ymax_bin/max_bin > 0.5){
  }
  else{
    /*x1 = 0.6;
    y1 = 0.7;
    x2 = 0.95;
    y2 = 0.9;*/
    int bin_leg_min = int(0.5*max_bin);
    for(unsigned int bin=bin_leg_min; bin < h_data->GetNbinsX()+1; bin++){
      if(h_data->GetBinContent(bin) / h_data->GetMaximum() > 0.01 ) scale_up = true; 
    }
  }
  
  if(name.find("llmass")!=string::npos) yscale*=1.3;
  
  if(ylog){
    float scale_for_log=1.;
    if(h_data->GetMaximum() > 100000) scale_for_log = 10;
    else if(h_data->GetMaximum() > 10000) scale_for_log = 10;
    else if(h_data->GetMaximum() > 1000) scale_for_log = 10;
    yscale*=scale_for_log;
    if(scale_up) yscale*= 10000;
  }
  else{
    yscale=1.2;
  }

  if(name.find("Tri")!=string::npos) {
    yscale*=1.5;
    if(ylog) yscale*=100.;
  }

  //  if(name.find("Eta")!=string::npos) yscale/=2.;
  float max_data = h_data->GetMaximum()*yscale;
  float max_bkg = h_up->GetMaximum()*yscale;
  
  if(max_data > max_bkg) return max_data;
  else return max_bkg;

  
  if(name.find("eemass")!=string::npos) yscale*=1.3;
  if(name.find("eta")!=string::npos) yscale*=2.5;
  if(name.find("MET")!=string::npos) yscale*=1.2;
  if(name.find("e1jj")!=string::npos) yscale*=1.2;
  if(name.find("Nje")!=string::npos) yscale*=1.2;
  if(name.find("l2jj")!=string::npos) yscale*=1.3;
  if(name.find("charge")!=string::npos) yscale*=1.5;
  if(name.find("deltaR")!=string::npos) yscale*=1.5;
  if(name.find("bTag")!=string::npos) yscale*=2.5;
  if(name.find("emujj")!=string::npos) yscale*=1.3;
  if(name.find("dijetmass")!=string::npos) yscale*=1.5;
  if(name.find("LeptonPt")!=string::npos) yscale*=0.7;
  if(name.find("secondElectronPt")!=string::npos) yscale*=1.2;
  
  
  return -1000.;
}



float GetTotal(TString cut, vector<pair<TString,float> > samples){
  
  float total(0.);
  for(vector<pair<TString,float> >::iterator it = samples.begin(); it!=samples.end(); it++){
    total += GetIntegral(cut,(*it).first,"MC");
  }
  
  return total;    
}


float GetStatError2(TString cut, vector<pair<TString,float> > samples){  
  float err = GetStatError(cut,samples);
  err = err*err;
  return err;
}

float GetStatError(TString cut, vector<pair<TString,float> > samples){  

  TString path  = mcloc;
  
  int isample=0;
  TFile* f0 =  TFile::Open((path+ fileprefix + samples.at(isample).first +  filepostfix).Data());  

  TH1* h_tmp=NULL;
  if((f0->Get(cut.Data()))){
    h_tmp = dynamic_cast<TH1*> ((f0->Get(cut.Data()))->Clone());
  }
    float stat_error(-99999.);

  while(!h_tmp){
    isample++;
    f0 =  TFile::Open((path+ fileprefix + samples.at(isample).first +  filepostfix).Data());
    if((f0->Get(cut.Data()))) h_tmp = dynamic_cast<TH1*> ((f0->Get(cut.Data()))->Clone());
  }
    
  for(unsigned int i=isample+1; i < samples.size(); i++){
    TFile* f =  TFile::Open((path+ fileprefix + samples.at(i).first +  filepostfix).Data());
    if(!f->Get(cut.Data())) continue;
    TH1* h = dynamic_cast<TH1*> ((f->Get(cut.Data()))->Clone());
    h_tmp->Add(h);
    f->Close();
    if(i == (samples.size()- 1)){
      stat_error = Error(h_tmp);
    }
  }

  if(samples.size()==1) stat_error = Error(h_tmp);
  
  f0->Close();
  return stat_error;
  
}




float GetIntegral(TString cut, TString isample, TString type){
  
  TString filepath = mcloc  + fileprefix + isample + filepostfix;
  if(type.Contains("data")) filepath=dataloc;
  TFile* f =  TFile::Open(( filepath.Data()));
    
  if(!((f->Get(cut.Data())))){
    cout << "Histogram " << cut << " in "  << mcloc+  fileprefix + isample + filepostfix << " not found" << endl;
    return 0.;
  }
  

  TH1* h = dynamic_cast<TH1*> ((f->Get(cut.Data())->Clone()));
  
  
  if(!h) {
    cout << "Histogram " << cut << " in "  <<  (mcloc + fileprefix + isample + filepostfix) << " not found" << endl;
    return 0.;
  }
  float integral = h->Integral();
 
  f->Close();

  return integral;
  
}

float GetError(TString cut, TString isample, TString type){
  
  TString filepath = mcloc  + fileprefix + isample + filepostfix;
  if(type.Contains("data")) filepath=dataloc;
  TFile* f =  TFile::Open(( filepath.Data()));
  
  TH1* h = dynamic_cast<TH1*> ((f->Get(cut.Data())->Clone()));
  
  if(!h) {
    cout << "Histogram " << cut << " in "  << (path+ fileprefix + isample + filepostfix) << " not found" << endl;
    return 0.;
  }
  
      
  float err = Error(h);
  
  f->Close();

  return err;
  
}


float GetNormErr(TString cut, vector<pair<TString,float> > samples){
  
  float norm_err(0.);
  int i=0;
  for( vector<pair<TString,float> >::iterator it = samples.begin(); it!=samples.end(); it++, i++){
    norm_err += GetIntegral(cut,it->first,"MC")* it->second*GetIntegral(cut,it->first,"MC")* it->second;
  }
  
  return sqrt(norm_err);
}

float GetNormErr2(TString cut,  vector<pair<TString,float> > samples){

  float err = GetNormErr(cut,samples);
  err = err*err;
  return err;
  
}



float GetErr2(TString cut, vector<pair<TString,float> > samples, TString err_type,TString var){
  
  float err = GetErr(cut,samples,err_type,var);
  err = err*err;
  return err;

}

float GetErr(TString cut, vector<pair<TString,float> > samples, TString err_type,TString var){

  float total(0.);
  float total_witherr(0.);
  
  int i=0;
  for(vector<pair<TString,float> >::iterator it = samples.begin(); it!=samples.end(); it++, i++){
    total += GetIntegral(cut,it->first,("MC/"));
    total_witherr += GetIntegral(cut,it->first,("MC_" + err_type));
  }
  
  
  float err=  total - total_witherr;      
  

  if(err_type.Contains("MM")){    
    if((var.Contains("UP")||var.Contains("Up"))) {
      if(( ( total - total_witherr)< 0.)) return fabs(err);
      else return total*0.1;
    }
    if((var.Contains("DOWN")||var.Contains("Down")) && ( ( total - total_witherr > 0.))) return err;
    else return total*0.1;  
  }
  
  if((var.Contains("UP")||var.Contains("Up"))) {
    if(( ( total - total_witherr)< 0.)) return fabs(err);
    else return 0.;
  }
  if((var.Contains("DOWN")||var.Contains("Down"))){
    if( ( total - total_witherr > 0.)) return err;
    else return 0.;  
  }

    return err;

}



float GetSystPercent(TString cut, TString syst, pair<vector<pair<TString,float> >,TString > samples ){

  return ( 100.*(GetSyst(cut, syst,samples )/ Calculate(cut,"Normal",samples)));
}



float GetSyst(TString cut, TString syst, pair<vector<pair<TString,float> >,TString > samples ){
  
  if(syst.Contains("Stat"))return GetStatError(cut,samples.first) ;    
  if(syst.Contains("Normalisation")) return GetNormErr(cut,samples.first);

  if(syst.Contains("TOTAL"))return sqrt(  GetStatError2(cut,samples.first)+
					  GetNormErr2(cut,samples.first));  
  return -9999.;
}


float Calculate(TString cut, TString variance, pair<vector<pair<TString,float> >,TString > samples ){
  
  if(samples.second.Contains("Nonprompt_DoubleEG")){
    if(variance.Contains("Normal"))  return GetTotal(cut,samples.first) ;  
    if(variance.Contains("StatErr")) return GetStatError(cut,samples.first) ;  
  }
  
  
  if(variance.Contains("Up") || variance.Contains("UP")||variance.Contains("Down") || variance.Contains("DOWN")  ){
    
    return fabs(sqrt( (GetNormErr2(cut,samples.first))));
    
  }
  
  if(variance.Contains("Normal")) return GetTotal(cut,samples.first) ;    
  
  if(variance.Contains("StatErr")) return GetStatError(cut,samples.first) ;    
    
  return -999999.;

}




float Error(TH1* h){
  double err ;
  double integral = h->IntegralAndError(0,h->GetNbinsX(),err,"");
  
  return err;
}



void SetUpMasterConfig(string name){
  
  // Get list of cuts to plot
  ifstream master_config_name_file(name.c_str());
  if(!master_config_name_file) {
    cerr << "Did not find " + name + ", exiting ..." << endl;
    return;
  }
  while(!master_config_name_file.eof()) {
    string tmp;
    string tmppath;
    master_config_name_file >> tmp;
    master_config_name_file >> tmppath;
    
    if(tmp=="END") break;
    if(tmp.find("#")!=string::npos) continue;
    cout << "tmp = " << tmp << endl;
    
    if(tmp=="mcpath") mcloc = tmppath;
    if(tmp=="datapath") dataloc = tmppath;
    if(tmp=="nonpromptpath") nonpromptloc = tmppath;
    if(tmp=="prefix") fileprefix = tmppath;
    if(tmp=="postfix") filepostfix = tmppath;
    if(tmp=="plottingpath") plotloc = tmppath;
    if(tmp=="cutpath")  cutloc = tmppath;
    if(tmp=="outputdir")    path = tmppath;

    if(tmp=="showdata")    {
      if (tmppath == "true") showdata=true;
      else showdata=false;
    }
    if(tmp=="ylog") {
      if (tmppath == "true")ylog=true;
      else ylog = false;
    }

    if(tmp=="usenp"){
      if (tmppath == "true")usenp = true;
      else usenp=false;
    }
    
    if(tmp=="samples_ss"){
      listofsamples_ss.push_back(tmppath);
    }
    if(tmp=="samples"){
      listofsamples.push_back(tmppath);
    }
    
    
    if(tmp=="histdir") histdir = tmppath;
    
    cutfile = cutloc;
    histfile =  plotloc;
  }
}

void  SetUpConfig(vector<pair<pair<vector<pair<TString,float> >, int >, TString > >& samples, vector<pair<pair<vector<pair<TString,float> >, int >, TString > >& samples_ss, vector<string>& cut_label){
    
  /// colours of histograms
  int tcol(0), zzcol(0), fcol(0), ttcol(0), zcol(0), wzcol(0), sscol(0),  wwcol(0), wcol(0),  ttvcol(0), higgscol(0), vvvcol(0), vvcol(0), vgammacol(0);
  
  // Get list of cuts to plot  
  ifstream colour_name_file("Config/colour.txt");
  if(!colour_name_file) {
    cerr << "Did not find Config/colour.txt, exiting ..." << endl;
    return;
  }
  while(!colour_name_file.eof()) {        
    string histname;
    int col;
    colour_name_file >> histname;    
    if(histname=="END") break;
    colour_name_file >> col;


    if(histname=="tcol")  tcol =col;
    if(histname=="ttcol") ttcol =col;
    if(histname=="zzcol") zzcol =col;
    if(histname=="fcol")  fcol =col;
    if(histname=="zcol")  zcol =col;
    if(histname=="wzcol") wzcol =col;
    if(histname=="wwcol") wwcol =col;
    if(histname=="wcol")  wcol =col;
    if(histname=="higgscol") higgscol =col;
    if(histname=="ttvcol") ttvcol = col;
    
  }
  
  /// Setup list of samples: grouped into different processes 
  vector<pair<TString,float> > top = InitSample("top");
  vector<pair<TString,float> > vv = InitSample("vv");
  
  // Zjet
  vector<pair<TString,float> > zlow = InitSample("dylow_");
  vector<pair<TString,float> > zhigh = InitSample("dyhigh_");

  vector<pair<TString,float> > w = InitSample("wjet_");
  vector<pair<TString,float> > other = InitSample("other_");

  /// QCD samples
  vector<pair<TString,float> > QCD = InitSample("qcd");
  /// ALL same sign processes

  vector<pair<TString,float> > prompt   = InitSample("prompt");


  /// NP is nonprompt
  vector<pair<TString,float> > np;
  np.push_back(make_pair("nonprompt_DoubleEG",0.34));
  
  vector<pair<TString,float> > cf;
  cf.push_back(make_pair("chargeflip",0.12));
  for( unsigned int i = 0; i < listofsamples.size(); i++){
    cout << listofsamples.at(i) << endl;
    if(listofsamples.at(i) =="other")samples.push_back(make_pair(make_pair(other,zzcol),"Other"));
    if(listofsamples.at(i) =="vv")samples.push_back(make_pair(make_pair(vv,zzcol),"VV"));
    if(listofsamples.at(i) =="dy")samples.push_back(make_pair(make_pair(zhigh,zcol),"DY#rightarrow ll"));
    if(listofsamples.at(i) =="dylow")samples.push_back(make_pair(make_pair(zlow,tcol),"DY#rightarrow ll; 10 < m(ll) < 50"));
    if(listofsamples.at(i) =="dyhigh")samples.push_back(make_pair(make_pair(zhigh,zcol),"DY#rightarrow ll; m(ll) > 50"));
    //if(listofsamples.at(i) =="top")samples.push_back(make_pair(make_pair(top,ttcol),"t#bar{t},t/#bar{t},t/#bar{t}W,t#bar{t}V"));
    if(listofsamples.at(i) =="top")samples.push_back(make_pair(make_pair(top,ttcol),"Top"));
    if(listofsamples.at(i) =="wjets")samples.push_back(make_pair(make_pair(w,wcol),"Wjets"));
    
    if(listofsamples.at(i) =="prompt")samples.push_back(make_pair(make_pair(prompt,vvcol),"Prompt Background"));
    if(listofsamples.at(i) =="qcd")samples.push_back(make_pair(make_pair(QCD,fcol),"QCD"));
    if(listofsamples.at(i) =="nonprompt_DoubleEG")samples.push_back(make_pair(make_pair(np,fcol),"Misid. Lepton Background"));   
    if(listofsamples.at(i) =="chargeflip")samples.push_back(make_pair(make_pair(cf,zcol),"Mismeas. Charge Background"));   
  }
  for( unsigned int i = 0; i < listofsamples_ss.size(); i++){
    if(listofsamples_ss.at(i) =="vv")samples_ss.push_back(make_pair(make_pair(vv,zzcol),"WZ,ZZ,WW"));
    if(listofsamples_ss.at(i) =="dylow")samples_ss.push_back(make_pair(make_pair(zlow,tcol),"DY#rightarrow ll; 10 < m(ll) < 50"));
    if(listofsamples_ss.at(i) =="dyhigh")samples_ss.push_back(make_pair(make_pair(zhigh,zcol),"DY#rightarrow ll; m(ll) > 50"));
    if(listofsamples_ss.at(i) =="top")samples_ss.push_back(make_pair(make_pair(top,ttcol),"t#bar{t},t/#bar{t},t/#bar{t}W,t#bar{t}V"));
    if(listofsamples_ss.at(i) =="wjets")samples_ss.push_back(make_pair(make_pair(w,wcol),"Wjets"));
    
    if(listofsamples_ss.at(i) =="prompt")samples_ss.push_back(make_pair(make_pair(prompt,vvcol),"Prompt Background"));
    if(listofsamples_ss.at(i) =="qcd")samples_ss.push_back(make_pair(make_pair(QCD,fcol),"QCD"));
    if(listofsamples_ss.at(i) =="nonprompt_DoubleEG")samples_ss.push_back(make_pair(make_pair(np,fcol),"Misid. Lepton Background"));
    if(listofsamples_ss.at(i) =="chargeflip")samples_ss.push_back(make_pair(make_pair(cf,zcol),"Mismeas. Charge Background"));
  }


  ///// Fix cut flow code
  caption="Number of events containing two prompt electrons and two jets, with Z peak removed, in 19 fb$^{-1}$ of CMS data at 8~TeV";
  hist = "/h_Nelectrons";
  columnname="";

  return;

}


TCanvas* CompDataMC(TH1* hdata, vector<THStack*> mcstack,TH1* hup, TH1* hdown,TH1* hup_nostat,TLegend* legend, const string hname, const  int rebin, double xmin, double xmax,double ymin, double ymax,string path , string folder, bool logy, bool usedata, TString channel) {
  
  ymax = GetMaximum(hdata, hup, ylog, hname, xmax);
  
  string cname;
  if(hdata) cname= string("c_") + hdata->GetName();
  else cname = string("c_") + ((TNamed*)mcstack.at(0)->GetHists()->First())->GetName();
  
  string label_plot_type = "";
  //Create Canvases

  unsigned int outputWidth = 1200;
  unsigned int outputHeight = 1400;

  TCanvas* canvas = new TCanvas((cname+ label_plot_type).c_str(), (cname+label_plot_type).c_str(), outputWidth,outputHeight);
  TCanvas* canvas_log = new TCanvas((cname+ label_plot_type+"log").c_str(), (cname+label_plot_type+"log").c_str(), outputWidth,outputHeight);
  
  

  // references for T, B, L, R                                                                                                                                         
  float T = 0.07*outputHeight;
  float B = 0.15*outputHeight;
  float L = 0.17*outputWidth;
  float R = 0.08*outputWidth;
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetLeftMargin( L/outputWidth );
  canvas->SetRightMargin( R/outputWidth );
  canvas->SetTopMargin( T/outputHeight );
  canvas->SetBottomMargin( B/outputHeight );
  canvas->SetTickx(0);
  canvas->SetTicky(0);

  
  std::string title=canvas->GetName();
  std::string tpdf = "/home/jalmond/WebPlots/"+ path + "/histograms/"+folder+"/"+title+".png";
  std::string tlogpdf = "/home/jalmond/WebPlots/"+ path + "/histograms/"+folder+"/"+title+"_log.png";
  
  ///####################   Standard plot

  if(!TString(hname).Contains("Tri")) {
    //if(!TString(hname).Contains("SSE")) {
      
      if(TString(hname).Contains("llmass")){canvas_log->SetLogy();canvas->SetLogy();}
      if(TString(hname).Contains("LeptonPt")){canvas_log->SetLogy();canvas->SetLogy();}
      //}
  }
  
  
  canvas->cd();
  
  //// %%%%%%%%%% TOP HALF OF PLOT %%%%%%%%%%%%%%%%%%
  TH1* h_nominal = MakeSumHist2(mcstack.at(0));
  SetNomBinError(h_nominal, hup, hdown);

  
  hdata->SetLineColor(kBlack);
  
  // draw data hist to get axis settings
  //hdata->GetYaxis()->SetTitleOffset(1.5);
  hdata->Draw("p9hist");
  TLatex label;
  label.SetTextSize(0.04);
  label.SetTextColor(2);
  label.SetTextFont(42);
  label.SetNDC();
  label.SetTextColor(1);
  //  //#  label.DrawLatex(0.6 ,0.34,"High Mass Region");
  //  label.DrawLatex(0.6 ,0.4,"Nvtx reweighted");
  

  //return canvas;  
  mcstack.at(0)->Draw("HIST9same");
  
  // draw axis on same canvas
  hdata->Draw("axis same");

  vector<float> err_up_tmp;
  vector<float> err_down_tmp;

  // get graph of errors for data., Set error 1.8 for 0 entry bins
  const double alpha = 1 - 0.6827;
  TGraphAsymmErrors * g = new TGraphAsymmErrors(hdata);
  for (int i = 0; i < g->GetN(); ++i) {
    int N = g->GetY()[i];
    double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
    double U =  (N==0) ?  ( ROOT::Math::gamma_quantile_c(alpha,N+1,1) ) :
      ( ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) );
      if ( N!=0 ) {
	g->SetPointEYlow(i, N-L );
	g->SetPointEXlow(i, 0);
	g->SetPointEYhigh(i, U-N );
	g->SetPointEXhigh(i, 0);
	err_down_tmp.push_back(N-L);
        err_up_tmp.push_back(U-N);
	
      }
      else {
	g->SetPointEYlow(i, 0.1);
	g->SetPointEXlow(i, 0.);
	g->SetPointEYhigh(i, 1.8);
	g->SetPointEXhigh(i, 0.);
	err_down_tmp.push_back(0.);
        err_up_tmp.push_back(1.8);
      }
  }
  

  gPad->Update();

  g->SetLineWidth(2.0);
  g->SetMarkerSize(0.);
  g->Draw(" p0" );

  //return canvas;  
  //hdata->GetYaxis()->SetRangeUser(0.01, ymax*10.);
  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(1.6);
  hdata->SetLineWidth(2.);
  hdata->Draw( ("same9p9hist"));
  

  legend->Draw();
  //return canvas;  
  CMS_lumi( canvas, 4, 11 );
  //return canvas;  
  canvas->Update();
  canvas->RedrawAxis();

  canvas->Print(tpdf.c_str(), ".png");


  //// %%%%%%%%%% PRINT ON LOG
  canvas_log->cd();
  canvas_log->SetLogy();     
  
  gPad->SetLogz(1);
  //// %%%%%%%%%% TOP HALF OF PLOT %%%%%%%%%%%%%%%%%%
  
  float scale_for_log = 1.;
  if(!ylog){
    ymax = GetMaximum(hdata, hup, !ylog, hname, xmax);
    hdata->GetYaxis()->SetRangeUser(0.2, ymax*scale_for_log);
  }
  //hdata->GetYaxis()->SetTitleOffset(1.6);
  hdata->Draw("p9hist");
  
  mcstack.at(0)->Draw("9HIST same");
  hdata->Draw("9samep9hist");
  hdata->Draw("axis same");
  
  gPad->Update();


  g->Draw(" p0" );

  
  legend->Draw();
  
  /// Make significance hist
  
  
  TH1* h_significance=(TH1F*)hdata->Clone("sig");
  TH1* h_divup=(TH1F*)hup->Clone("divup");
  TH1* h_divdown=(TH1F*)hdown->Clone("divdown");
  
  TH1F* hdev = (TH1F*)hdata->Clone("hdev");
  TH1F* hdev_err = (TH1F*)hdata->Clone("hdev_err");
  TH1F* hdev_err_stat = (TH1F*)hdata->Clone("hdev_err_stat");
    
  TH1* errorbandratio = (TH1*)h_nominal->Clone("AAA");
    
  hdata->GetXaxis()->SetLabelSize(0.); ///
  hdata->GetXaxis()->SetTitle("");
  
  h_divup->Divide(h_nominal);
  h_divdown->Divide(h_nominal);
  
  
  
  // How large fraction that will be taken up by the data/MC ratio part
  double FIGURE2_RATIO = 0.3;
  double SUBFIGURE_MARGIN = 0.0;
  canvas_log->SetBottomMargin(FIGURE2_RATIO);
  TPad *p = new TPad( "p_test", "", 0, 0, 1, 1.0 - SUBFIGURE_MARGIN, 0, 0, 0);  // create new pad, fullsize to have equal font-sizes in both plots
  p->SetTopMargin(1-FIGURE2_RATIO);   // top-boundary (should be 1 - thePad->GetBottomMargin() )
  p->SetFillStyle(0);     // needs to be transparent
  p->Draw();
  p->cd();
  p->SetTicks(0,1);

  
  
  Double_t *staterror;
  
  for (Int_t i=1;i<=hdev->GetNbinsX()+1;i++) {
    hdev_err->SetBinContent(i, 1.0);
    if(h_nominal->GetBinContent(i) > 0 &&  hdev->GetBinContent(i) > 0){
      hdev_err->SetBinError(i, (hup_nostat->GetBinContent(i)- h_nominal->GetBinContent(i))/h_nominal->GetBinContent(i) );
    }
    else{
      hdev_err->SetBinError(i, 0.0);
    }
  }
  for (Int_t i=1;i<=hdev->GetNbinsX()+1;i++) {
    hdev_err_stat->SetBinContent(i, 1.0);
    if(h_nominal->GetBinContent(i) > 0 &&  hdev->GetBinContent(i) > 0){
      hdev_err_stat->SetBinError(i, (hup->GetBinContent(i)-h_nominal->GetBinContent(i))/h_nominal->GetBinContent(i) );
    }
    else{
      hdev_err_stat->SetBinError(i, 0.0);
    }
  } 
  
  
  for (Int_t i=1;i<=hdev->GetNbinsX()+1;i++) {
    if(h_nominal->GetBinContent(i) > 0 &&  hdev->GetBinContent(i) > 0){
      hdev->SetBinContent(i, hdev->GetBinContent(i)/ h_nominal->GetBinContent(i));
      //hdev->SetBinContent(i, h_nominal->GetBinContent(i)/ h_nominal->GetBinContent(i));
      hdev->SetBinError(i, 0.01);
    }
    else {
      hdev->SetBinContent(i, -99);
      hdev->SetBinError(i, 0.);
    }
  }
  
  /// set errors for datamc plot
  TGraphAsymmErrors * gratio = new TGraphAsymmErrors(hdev);

  for (int i = 0; i < gratio->GetN(); ++i) {
    
    if(err_down_tmp.at(i)  !=0.) {
      gratio->SetPointEYlow(i, err_down_tmp.at(i) / h_nominal->GetBinContent(i+1) );
      gratio->SetPointEXlow(i, 0);
      gratio->SetPointEYhigh(i, err_up_tmp.at(i) /h_nominal->GetBinContent(i+1));
      gratio->SetPointEXhigh(i, 0);
    }
    else{
      gratio->SetPointEYlow(i, 0);
      gratio->SetPointEXlow(i, 0);
      gratio->SetPointEYhigh(i, 1.8 / h_nominal->GetBinContent(i+1));
      gratio->SetPointEXhigh(i, 0);
    }
  }
  
  
  //////////// Plot all
  
  
  hdev->GetYaxis()->SetTitle( "Data / #Sigma MC" );
  hdev->GetYaxis()->SetRangeUser(0.25,+1.75);
  //hdev->GetYaxis()->SetLabelOffset(0.005, "Y");
  hdev->GetYaxis()->SetLabelSize(0.02);
  hdev->GetYaxis()->SetTitleSize(0.02);
  hdev->GetYaxis()->SetTitleOffset(1.);

  

  hdev->GetYaxis()->SetNdivisions(3);
  hdev->SetMarkerStyle(20);
  //hdev->SetMarkerSize(2.3);
  hdev_err_stat->SetMarkerSize(0.);
  hdev_err->SetMarkerSize(0.);
  hdev->SetLineColor(kBlack);
  hdev_err->SetFillColor(kRed);
  hdev_err->SetLineColor(kRed);
  hdev_err->SetFillStyle(3444);
  hdev_err_stat->SetFillColor(kOrange-9);
  hdev_err_stat->SetLineColor(kOrange-9);
  hdev->Draw("phist");
  
  hdev_err_stat->Draw("sameE4");
  hdev_err->Draw("sameE4");
  gratio->SetLineWidth(2.0);
  gratio->SetMarkerSize(0.);
  gratio->Draw(" p0" );
  hdev->Draw("same p hist");
      
    
  TLine *devz = new TLine(hdev->GetBinLowEdge(hdev->GetXaxis()->GetFirst()),1.0,hdev->GetBinLowEdge(hdev->GetXaxis()->GetLast()+1),1.0  );
  devz->SetLineWidth(1);
  devz->SetLineStyle(1);
  devz->Draw("SAME");
  
  
  
  

  CMS_lumi( canvas_log, 4, 11 );
  canvas_log->Update();
  canvas_log->RedrawAxis();
  canvas_log->Print(tlogpdf.c_str(), ".png");
  gPad->RedrawAxis();
  
  return canvas;

}



TH1* MakeSumHist2(THStack* thestack){

  TH1* hsum=0;
  TList* list = thestack->GetHists();
  TIter it(list, true);
  TObject* obj=0;
  while( (obj = it.Next()) ) {
    TH1* h = dynamic_cast<TH1*>(obj);
    if(!hsum) hsum = (TH1*)h->Clone( (string(h->GetName()) + "_sum").c_str() );
    else {
      hsum->Add(h, 1.0);
    }
  }//hist loop

  return hsum;
}



void SetNomBinError(TH1* hnom, TH1* hup, TH1* hdown){

  for(int i=1; i < hnom->GetNbinsX()+1; i++){

    float err1 = fabs(hnom->GetBinContent(i)- hup->GetBinContent(i));
    float err2 = fabs(hnom->GetBinContent(i)- hdown->GetBinContent(i));

    if(err1 > err2 ) hnom->SetBinError(i, err1);
    if(err2 > err1 ) hnom->SetBinError(i, err2);
  }
  return;
}



TH1* MakeErrorBand(TH1* hnom, TH1* hup, TH1* hdown){

  TH1* errorband = (TH1*)hnom->Clone("aa");
  for(int i=1; i < errorband->GetNbinsX()+1; i++){

    float bin_content = (hup->GetBinContent(i)+ hdown->GetBinContent(i))/2.;
    float bin_error = (hup->GetBinContent(i)- hdown->GetBinContent(i))/2.;

    errorband->SetBinContent(i,bin_content);
    errorband->SetBinError(i,bin_error);
  }

  errorband->SetFillStyle(3444);
  errorband->SetFillColor(kBlue-8);
  errorband->SetMarkerSize(0);
  errorband->SetMarkerStyle(0);
  errorband->SetLineColor(kWhite);
  errorband->Draw("E2Same");

  return errorband;

}


void MakeLabel(float rhcol_x, float rhcol_y){
  TLatex label;
  label.SetTextSize(0.04);
  label.SetTextColor(2);
  label.SetTextFont(42);
  label.SetNDC();
  label.SetTextColor(1);
  label.DrawLatex(rhcol_x,rhcol_y,"#int L dt = 20.4 fb^{-1}");
  label.DrawLatex(rhcol_x + 0.2,rhcol_y ,"#sqrt{s}= 8 TeV");
  label.SetTextSize(0.045);

  label.DrawLatex(rhcol_x+0.115, rhcol_y + 0.09,"Work In Progress");
  label.SetTextFont(72);
  label.DrawLatex(rhcol_x, rhcol_y + 0.09,"CMS");

  return;
}



void
CMS_lumi( TPad* pad, int iPeriod, int iPosX )
{
  bool outOfFrame    = false;
  if( iPosX/10==0 )
    {
      outOfFrame = true;
    }
  outOfFrame=true;
  int alignY_=3;
  int alignX_=2;
  if( iPosX/10==0 ) alignX_=1;
  if( iPosX==0    ) alignY_=1;
  if( iPosX/10==1 ) alignX_=1;
  if( iPosX/10==2 ) alignX_=2;
  if( iPosX/10==3 ) alignX_=3;
  int align_ = 10*alignX_ + alignY_;

  float H = pad->GetWh();
  float W = pad->GetWw();
  float l = pad->GetLeftMargin();
  float t = pad->GetTopMargin();
  float r = pad->GetRightMargin();
  float b = pad->GetBottomMargin();
  float e = 0.025;

  pad->cd();

  TString lumiText;
  if( iPeriod==1 )
    {
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
    }
  else if ( iPeriod==2 )
    {
      lumiText += lumi_8TeV;
      lumiText += " (8 TeV)";
    }
  else if( iPeriod==3 )
    {
      lumiText = lumi_8TeV;
      lumiText += " (8 TeV)";
      lumiText += " + ";
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
    }
  else if ( iPeriod==4 )
    {
      lumiText += lumi_13TeV_2016_eg_BtoH;
      lumiText += " (13 TeV)";
    }
  else if ( iPeriod==7 )
    {
      if( outOfFrame ) lumiText += "#scale[0.85]{";
      lumiText += lumi_13TeV_2016_eg_BtoH;
      lumiText += " (13 TeV)";
      lumiText += " + ";
      lumiText += lumi_8TeV;
      lumiText += " (8 TeV)";
      lumiText += " + ";
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
      if( outOfFrame) lumiText += "}";
    }
  else if ( iPeriod==12 )
    {
      lumiText += "8 TeV";
    }


  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);

  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  latex.SetTextFont(42);
  latex.SetTextAlign(31);
  latex.SetTextSize(lumiTextSize*t);
  latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  if( outOfFrame )
    {
      latex.SetTextFont(cmsTextFont);
      latex.SetTextAlign(11);
      latex.SetTextSize(cmsTextSize*t);
      latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);

      latex.SetTextFont(extraTextFont);
      latex.SetTextAlign(align_);
      latex.SetTextSize(extraTextSize*t);
      latex.DrawLatex(l + relPosX*(1-l-r), 1-t+lumiTextOffset*t , extraText);

    }
  return;
  pad->cd();

  //  writeExtraText=true;
  float posX_;
  if( iPosX%10<=1 )
    {
      posX_ =   l + relPosX*(1-l-r);
    }
  else if( iPosX%10==2 )
    {
      posX_ =  l + 0.5*(1-l-r);
    }
  else if( iPosX%10==3 )
    {
      posX_ =  1-r - relPosX*(1-l-r);
    }
  float posY_ = 1-t - relPosY*(1-t-b);
  if( !outOfFrame )
    {
      if( drawLogo )
        {
          posX_ =   l + 0.045*(1-l-r)*W/H;
          //posY_ = 1-t - 0.045*(1-t-b);
	  posY_ =  1-t+lumiTextOffset*t;

          float xl_0 = posX_;
          float yl_0 = posY_ - 0.15;
          float xl_1 = posX_ + 0.15*H/W;
          float yl_1 = posY_;
          TASImage* CMS_logo = new TASImage("CMS-BW-label.png");
          TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
          pad_logo->Draw();
          pad_logo->cd();
          CMS_logo->Draw("X");
          pad_logo->Modified();
          pad->cd();
        }
      else
        {
          latex.SetTextFont(cmsTextFont);
          latex.SetTextSize(cmsTextSize*t);
          latex.SetTextAlign(align_);
          //latex.DrawLatex(posX_, posY_, cmsText);
          if( writeExtraText )
            {
              latex.SetTextFont(extraTextFont);
              latex.SetTextAlign(align_);
              latex.SetTextSize(extraTextSize*t);
	      //latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText);
	      //latex.DrawLatex(posX_, 1-t+lumiTextOffset*t, extraText); 
            }
        }
    }
  else if( writeExtraText )
    {
      if( iPosX==0)
        {
          posX_ =   l +  relPosX*(1-l-r);
          posY_ =   1-t+lumiTextOffset*t;
        }
      latex.SetTextFont(extraTextFont);
      latex.SetTextSize(extraTextSize*t);
      latex.SetTextAlign(align_);
      posX_ =   l +  relPosX*(1-l-r);
      posY_ = 1-t+lumiTextOffset*t;

      
      latex.DrawLatex(posX_, 1-t+lumiTextOffset*t , extraText);
      //latex.DrawLatex(posX_, 1-t+lumiTextOffset*t, extraText);
      
    }
  return;
}





