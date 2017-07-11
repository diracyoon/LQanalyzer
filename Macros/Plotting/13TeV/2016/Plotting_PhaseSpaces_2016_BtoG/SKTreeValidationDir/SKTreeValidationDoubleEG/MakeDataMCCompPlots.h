#ifndef MAKEDATAMCCOMPLOTS_H__
#define MAKEDATAMCCOMPLOTS_H__

////////////////////////////////////////////////////////////
// STD includes
////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <map>
#include <sstream>

////////////////////////////////////////////////////////////
/// Local includes
////////////////////////////////////////////////////////////
#include "HistUtils.hpp"

////////////////////////////////////////////////////////////
///ROOT includes
////////////////////////////////////////////////////////////
#include "TH1.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include <TStyle.h>
#include "TLatex.h"

using namespace std;

int MakePlots(std::string hist);

void MakeCutFlow(std::string hist);
int MakeCutFlow_Plots(string configfile);
void PrintCanvas(TCanvas* c1, std::string folder, std::string plotdesciption,  std::string title);
bool repeat(string hname);
TLegend* MakeLegend(map<TString, TH1*> legmap,TH1* h_legdata, bool rundata, bool log, float ymax, float xmax);
TH1* MakeDataHist(string name, double xmin, double xmax, TH1* h_up,bool ylog , int rebin);
void CheckHist(TH1* h);
void CheckSamples(int nsamples);
vector<pair<TString,float> >  InitSample (TString sample);
THStack* MakeStack(vector<pair<pair<vector<pair<TString,float> >, int >, TString > > sample, TString type, string name, float xmin, float xmax,map<TString, TH1*>& legmap, int rebin, bool include_syst);
void SetErrors(TH1* hist, float normerr, bool usestat);
TH1* MakeStackUp(map<TString, TH1*> map_of_stacks, TString clonename);
TH1* MakeStackDown(map<TString, TH1*> map_of_stacks, TString clonename);
TH1* MakeSumHist(THStack* thestack);
float  GetMaximum(TH1* h_data, TH1* h_up, bool ylog, string name, float xmax);
void SetTitles(TH1* hist, string name);
bool HistInGev(string name);
void fixOverlay();
void setTDRStyle();
TCanvas* CompDataMC(TH1* hdata,  vector<THStack*> mcstack,TH1* hup, TH1* hdown, TH1* hupnostat,TLegend* legend, const string hname, const int rebin, double xmin, double xmax,double ymin, double ymax,string path , string folder, bool logy, bool usedata, TString channel);
TH1* MakeSumHist2(THStack* thestack);
TH1* MakeErrorBand(TH1* hnom, TH1* hup, TH1* hdown);
void SetNomBinError(TH1* hnom, TH1* hup, TH1* hdown);
void MakeLabel( float rhcol_x, float rhcol_y);
std::map<std::string,std::string> _htmls;

float GetSyst(TString cut, TString syst, pair<vector<pair<TString,float> >,TString > samples );
float GetSystPercent(TString cut, TString syst, pair<vector<pair<TString,float> >,TString > samples );
void setZ(bool useAlpgen);
float Calculate(TString cut, TString variance,  pair<vector<pair<TString,float> >, TString >  samples);
void  SetUpConfig(vector<pair<pair<vector<pair<TString,float> >, int >, TString > >& samples , vector<pair<pair<vector<pair<TString,float> >, int >, TString > >& samples_ss, vector<string>& cut_label);
void  SetUpMasterConfig(std::string filename);

////// For cutflow
float Error(TH1* h);
float GetTotal(TString cut, vector<pair<TString,float> > sample);
float GetStatError(TString cut, vector<pair<TString,float> > sample);
float GetStatError2(TString cut, vector<pair<TString,float> > sample);
float GetIntegral(TString cut, TString isample, TString type);
float GetNormErr(TString cut,  vector<pair<TString,float> > samples);
float GetNormErr2(TString cut,  vector<pair<TString,float> > samples);
float GetErr(TString cut,  vector<pair<TString,float> > samples, TString err_type,TString var);
float GetErr2(TString cut,  vector<pair<TString,float> > samples, TString err_type,TString var);
float GetError(TString cut, TString isample, TString type);

//// GLOBAL VARIABLES
int isig=0;
map<string,int> norepeatplot;
TString columnname="";
TString caption="";
float scale=1.;

std::string hist;
bool showdata=true;
std::string cutfile;
std::string histfile;
bool ylog;
bool usenp(false);

TString channel;

std::string path;
std::string message;
std::string fileprefix="";
std::string filepostfix = "";

std::ofstream page;
std::ofstream histpage;

vector<string> cuts;
vector<string> allcuts;
vector<string> listofsamples;
vector<string> listofsamples_ss;

//// Standard bkg folders
string  mcloc="";
/// Data folder
string dataloc = "";
/// data driven
string nonpromptloc= "";
string plotloc ="";
string cutloc ="";
string histdir="";
string output_index_path="";
string output_path = "";


void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);

  tdrStyle->SetEndErrorSize(2);
  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.1);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

  // For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);

  // For the axis titles:
  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.07, "XYZ");
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);

  // For the axis labels:
  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:
  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  
  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->cd();

  

}

// fixOverlay: Redraws the axis
void fixOverlay() {
  gPad->RedrawAxis();
}




#endif // MakeDataMCCompPlots.h
