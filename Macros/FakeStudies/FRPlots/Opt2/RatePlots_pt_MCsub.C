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
TH1D* MakePRPlot(TString sample, TString hist);
TH1D* MakeFRPlot(TString sample, TString hist);


void FR_pt(){
 
  setTDRStyle();
  TString cut="ElELECTRON16_FR_MVA_TIGHT_CC_dijet_mva_40_pt_eta";
  TString cut1="ElMVA_TIGHT_40_pt_eta";//Eldijet_mva0.940963_iso0.05_dxy0.01_dz0.1_40_pt_eta";
  TString cut2="Elloose2_dijet_mva0.940963_iso0.05_dxy0.05_dz0.1_40_pt_eta";
  
  TH1D* h = MakeFRPlot("",cut1);
  TH1D* h2 = MakeFRPlot("",cut2);
    
  h->SetLineColor(kRed-2);
  h->GetXaxis()->SetRangeUser(15.,90.);
  
  TCanvas* c1 = new TCanvas("TL", "TL", 800, 600);
  
  
  TLegend *legend = new TLegend(.65, 0.5, .9, 0.7);
  legend->SetFillColor(10);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.03);

  h->GetYaxis()->SetRangeUser(0.,1.);
  h->Draw("hist");
  h2->Draw("histsame");
  legend->AddEntry(h,"MVA(tight)_0.01","l");
  legend->AddEntry(h2,"MVA(tight)_0.05","l");
  legend->Draw();
  
  c1->SaveAs("TL_mcsub_pt_pog_mva.pdf");
}


TH1D* MakePRPlot(TString sample, TString hist){

  TH1::SetDefaultSumw2();
  TString datapath="/data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer/data/output/CAT/FakeRateCalculator_El/periodBtoH/FakeRateCalculator_El_data_DoubleEG_cat_v8-0-7.root";
  TString path="/data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer/data/output/CAT/FakeRateCalculator_El/periodBtoH/";

  TFile * fdata  = new TFile(datapath);

  TH2D* h2loose = (TH2D*)fdata->Get(("Prompt_LooseEl_"+hist).Data());
  TH2D* h2tight = (TH2D*)fdata->Get("Prompt_TightEl_"+hist);
  
  
  if(!h2loose) exit(1);
  if(!h2tight) exit(1);

  TH1D *hloose = h2loose->ProjectionX();
  TH1D *htight = h2tight->ProjectionX();
  

  if(!hloose) exit(1);
  if(!htight) exit(1);
  htight->Divide(hloose);

  htight->SetLineWidth(3.5);

  return htight;
}

TH1D* MakeFRPlot(TString sample, TString hist){

  TH1::SetDefaultSumw2();
  TString datapath="/data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer/data/output/CAT/FakeRateCalculator_El/periodBtoH/FakeRateCalculator_El_data_DoubleEG_cat_v8-0-7.root";
  TString path="/data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer/data/output/CAT/FakeRateCalculator_El/periodBtoH";

  TFile * fmc= new TFile((path+"/FakeRateCalculator_El_mc_v8-0-7.root"));

  TFile * fdata  = new TFile(datapath);


  TH2D* h2loose = (TH2D*)fdata->Get(("Loose"+hist).Data());
  TH2D* h2tight = (TH2D*)fdata->Get("Tight"+hist);


  TH2D* h2mcloose = dynamic_cast<TH2D*>(fmc->Get("Loose"+hist)->Clone("Loose"+hist));
  TH2D* h2mctight = dynamic_cast<TH2D*>(fmc->Get("Tight"+hist)->Clone("Tight"+hist));


  if(!h2loose) exit(1);
  if(!h2tight) exit(1);

  //h2loose->Add(h2mcloose,-1.2);
  //h2tight->Add(h2mctight,-1.2);

  TH1D *hloose = h2loose->ProjectionX();
  TH1D *htight = h2tight->ProjectionX();


  if(!hloose) exit(1);
  if(!htight) exit(1);
  htight->Divide(hloose);

  htight->SetLineWidth(3.5);

  return htight;
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
