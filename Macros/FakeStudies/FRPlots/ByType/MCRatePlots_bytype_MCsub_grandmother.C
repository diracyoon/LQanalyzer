#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TKey.h"
#include <iostream>
#include <fstream>
#include <TStyle.h>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include <sstream>      // std::stringstream


//#include "Macro.h"
void setTDRStyle();
TH1D* MakeFRPlot(TString sample, TString id, TString hist);
void MCRatePlots_bytype_MCsub(TString type);

void MCRatePlots_bytype_MCsub(){
  MCRatePlots_bytype_MCsub("10");
  MCRatePlots_bytype_MCsub("11");
  MCRatePlots_bytype_MCsub("12");
  MCRatePlots_bytype_MCsub("16");
  MCRatePlots_bytype_MCsub("17");
  MCRatePlots_bytype_MCsub("18");
  MCRatePlots_bytype_MCsub("19");
  MCRatePlots_bytype_MCsub("1");
  MCRatePlots_bytype_MCsub("21");
  MCRatePlots_bytype_MCsub("22");
  MCRatePlots_bytype_MCsub("23");
  MCRatePlots_bytype_MCsub("24");
  MCRatePlots_bytype_MCsub("25");
  MCRatePlots_bytype_MCsub("28");
  MCRatePlots_bytype_MCsub("2");
  MCRatePlots_bytype_MCsub("6");
  MCRatePlots_bytype_MCsub("7");
  MCRatePlots_bytype_MCsub("8");
  MCRatePlots_bytype_MCsub("9");
}


void MCRatePlots_bytype_MCsub(TString type){
 
  setTDRStyle();

  TString sample="WJets";
  TString id="MMHNTIGHT";
  TString hist="ByType_grandmother_Full_"+type;
    
  TH1D* h = MakeFRPlot(sample, id, hist);
  TH1D* ht = MakeFRPlot(sample, id, hist+"_tight");
  
  ht->Divide(h);
  
  ht->SetLineColor(kRed-2);

  TCanvas* c1 = new TCanvas("TL", "TL", 800, 600);
  
  
  TLegend *legend = new TLegend(.65, 0.5, .9, 0.7);
  legend->SetFillColor(10);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.03);

  ht->GetYaxis()->SetRangeUser(0.,2.);
  ht->Draw("histEtext");
  
  //legend->Draw();
  
  ht->GetXaxis()->SetRangeUser(0.,100.);
  c1->SaveAs((sample + "_"+id + "_"+hist+"_pdgid_0_100"+type+".pdf").Data());
  ht->GetXaxis()->SetRangeUser(100.,400.);
  c1->SaveAs((sample + "_"+id + "_"+hist+"_pdgid_100_400"+type+".pdf").Data());
  ht->GetXaxis()->SetRangeUser(400.,500.);
  c1->SaveAs((sample + "_"+id + "_"+hist+"_pdgid_400_500"+type+".pdf").Data());
  ht->GetXaxis()->SetRangeUser(500.,600.);
  c1->SaveAs((sample + "_"+id + "_"+hist+"_pdgid_500_600"+type+".pdf").Data());
  ht->GetXaxis()->SetRangeUser(600.,1000.);
  c1->SaveAs((sample + "_"+id + "_"+hist+"_pdgid_600_1000"+type+".pdf").Data());
  ht->GetXaxis()->SetRangeUser(1000.,2000.);
  c1->SaveAs((sample + "_"+id + "_"+hist+"_pdgid_1000_5000"+type+".pdf").Data());




}

TH1D* MakeFRPlot(TString sample, TString id, TString hist){
  
  TH1::SetDefaultSumw2(); 
  TString datapath="/data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer/data/output/CAT/HNDiElectron/periodBtoH/HNDiElectron_"+sample+"_cat_v8-0-7.root";
  TFile * fdata  = new TFile(datapath);
  
  TH1D* h= (TH1D*)fdata->Get(id+hist);
  cout << id+hist << endl;
  if(!h) exit(1);

  return h;
}




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
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
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
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);
  tdrStyle->SetEndErrorSize(2);
  //  tdrStyle->SetErrorMarker(20);
  //  tdrStyle->SetErrorX(0.);

  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);

  // tdrStyle->SetDateY(Float_t y = 0.01);

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
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

  // For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);

  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.4);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset



  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.4);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

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
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

}
