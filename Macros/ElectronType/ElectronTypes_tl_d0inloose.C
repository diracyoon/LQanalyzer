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
TH1F* MakePlot(int bin1, int bin2,TString sample, int shift);
void MakeRangePlot(int bin1, int bin2, TString label);


void ElectronType(){
 
  MakeRangePlot(0,33,"range0-33");
  MakeRangePlot(0,8,"range1");
  MakeRangePlot(8,16,"range2");
  MakeRangePlot(16,25,"range3");
  MakeRangePlot(25,33,"range4");

}

void MakeRangePlot(int bin1, int bin2 , TString label){  
  setTDRStyle();

  TH1F* h = MakePlot(bin1,bin2,"DYJets",2);
  h->SetMarkerColor(kRed-2);
  TH1F* h2 = MakePlot(bin1,bin2,"TT_powheg",3);
  TH1F* h3 = MakePlot(bin1,bin2,"WJets",4);
  TH1F* h4 = MakePlot(bin1,bin2,"ZGto2LG",5);
  TH1F* h5 = MakePlot(bin1,bin2,"QCD_Pt-30to50_EMEnriched",6);
  TH1F* h6 = MakePlot(bin1,bin2,"qcd_30to80_bctoe",7);
  TH1F* h7 = MakePlot(bin1,bin2,"QCD_DoubleEMEnriched_30-inf_mgg40to80",8);
  
  TH1F* hall = MakePlot(bin1,bin2,"mc",0);

  
  h2->SetMarkerColor(kCyan);
  h3->SetMarkerColor(kBlue-3);
  h4->SetMarkerColor(kGreen+4);
  h5->SetMarkerColor(kGreen-2);
  h6->SetMarkerColor(kRed-2);
  h7->SetMarkerColor(kRed+2);
  
  h2->SetLineColor(kCyan);
  h3->SetLineColor(kBlue-3);
  h4->SetLineColor(kGreen+4);
  h5->SetLineColor(kGreen-2);
  h6->SetLineColor(kRed-2);
  h7->SetLineColor(kRed-2);
  
  hall->SetLineWidth(3.);
  hall->SetLineColor(kRed);
  hall->SetLineStyle(5);

  TCanvas* c1 = new TCanvas("TL", "TL", 800, 600);
  
  
  TLegend *legend = new TLegend(.3, 0.6, .7, 0.9);
  legend->SetFillColor(10);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.04);
  
  h2->SetMarkerStyle(2.);
  h3->SetMarkerStyle(3.);
  h4->SetMarkerStyle(4.);  h5->SetMarkerStyle(5.);
  h6->SetMarkerStyle(2.);
  h7->SetMarkerStyle(2.);
  
  h->GetYaxis()->SetRangeUser(0.,2.);
  //h->GetXaxis()->SetRangeUser( h->GetBinLowEdge(bin1), h->GetBinLowEdge(bin2+1));
  
  h->Draw("pE");
  h2->Draw("psameE");
  h3->Draw("psameE");
  h4->Draw("psameE");
  h5->Draw("psameE");
  h6->Draw("psameE");
  h7->Draw("psameE");
  hall->Draw("histsame");

  legend->AddEntry(h,"DY","p");
  legend->AddEntry(h2,"Top","p");
  legend->AddEntry(h3,"W","p");
  legend->AddEntry(h4,"ZG","p");
  legend->AddEntry(h5,"QCD/EM","p");
  legend->AddEntry(h6,"QCD/bc","p");
  legend->AddEntry(h7,"DiElQCD","p");
  legend->AddEntry(hall,"All","l");
  legend->Draw();
  
  c1->SaveAs("TL(dxy)_"+label+".pdf");
}


TH1F* MakePlot(int bin1, int bin2,TString sample, int shift){
  
  TH1::SetDefaultSumw2(); 
  int nbins=(bin2-bin1)*10;
  if(shift==0) nbins=(bin2-bin1);
  TH1F* h =  new TH1F ("h2","h2", nbins, double(bin1), double(bin2));
  h->GetXaxis()->SetTitle("GetType()");
  h->GetYaxis()->SetTitle("T/L");
  int xbin(0);
  for(int j=bin1-1; j < bin2; j++, xbin++){
    
    
    stringstream ss;
    ss << j;
    string str = "ELType"+ss.str();
    
    TString path= "/data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer//data/output/CAT/ElectronTypes/periodBtoH/";
    TFile * fmc = new TFile(path+"ElectronTypes_"+sample+"_cat_v8-0-6.root");
    
    TH1* hpt= (TH1*)fmc->Get(( str + "_loose_dxy_pt").c_str());
    TH1* h_tightpt= (TH1*)fmc->Get(( str + "_tight_pt").c_str());
    
    if(!h_tightpt) {
      h->Fill(j, 0.);
      
    }
    else if(!hpt){
      h->Fill(j, 0.);
      
    }
    else{
      
      float tl = h_tightpt->Integral()/hpt->Integral();
      double err;
      //h_tightpt->Divide(hpt);
      //h_tightpt->Scale(2.);
      
      h_tightpt->IntegralAndError(0, 34,err);
      err=err/hpt->Integral();
      int bin=( (xbin) * 10) + shift;
      if(shift==0) bin = xbin+1;
      h->SetBinContent(bin, tl);
      h->SetBinError(bin, err);
    }
  }
  
  return h;
}


TH1* makehist(TString histname, TString type,TString path, TString cut, TString label ){

  TLegend *legend = new TLegend(.7, 0.7, .9, 0.9);
  legend->SetFillColor(10);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.04);


  TFile * fmc = new TFile(path+"ElectronTypes_mc_cat_v8-0-6.root");

  TH1* h= (TH1*)fmc->Get(( type + cut));

  if(!h) return h;
  h->Scale(1./ h->Integral());
  h->GetXaxis()->SetTitle(label);
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
