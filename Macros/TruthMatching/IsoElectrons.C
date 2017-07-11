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
#include <TStyle.h>
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"

void setTDRStyle();

TH1F* MakeHist(TString path, TString name);
float CoutEff(TH1F* h);

float CoutEff(TH1F* h){


  int bin1 = h->FindBin(0.0);
  int bin2 = h->FindBin(0.05);
  int bin3 = h->FindBin(2.);
  return (h->Integral(bin1,bin2)/ h->Integral() );


}

TH1F* MakeHist(TString path, TString name){

  TFile * fmc = new TFile(path);

  if(!fmc){
    cout << "No MC file" << endl;
  }
  
  
  TH1F* hist = (TH1F*)fmc->Get(name);
  if(!hist) return hist;
  cout << name << " eff . = " << CoutEff(hist) << endl;


  hist->Scale(1/hist->Integral());
  hist->SetLineWidth(2.);
  return hist;
  
}

void IsoElectrons(TString path){
  
  path="/data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer//data/output/CAT/HNDiElectron/periodBtoH/HNDiElectron_DYJets_cat_v8-0-6.root";
  /// Set Plotting style
  setTDRStyle();
  gStyle->SetPalette(1);
  

  
  if(1){
    TH1F* hz_ss_conv = MakeHist(path,"SSMatched_conv_iso");
    TH1F* hz_ss_tau = MakeHist(path,"SSMatched_tau_iso");
    TH1F* hz_ss_cftau = MakeHist(path,"SSMatched_cftau_iso");
    TH1F* hz_ss_x  =MakeHist(path,"SSMatched_x_iso");
    TH1F* h =MakeHist(path,"ZMatched_iso");
    TH1F * h_m_wz = MakeHist(path,"Matched_W_Z_iso");
    TH1F * h_m_cf = MakeHist(path,"Matched_CF_iso");
    TH1F * h_m_cfconv = MakeHist(path,"Matched_CF_conv_iso");
    TH1F * h_m_cfnotconv = MakeHist(path,"Matched_CF_notconv_iso");
    TH1F * h_m_notcf= MakeHist(path,"Matched_noCF_iso");
    TH1F * h_m_conv=  MakeHist(path,"Matched_CONV_iso");
    TH1F * h_m_convnocf =  MakeHist(path,"Matched_CONV_notCF_iso");
    TH1F * h_m_noconv=  MakeHist(path,"Matched_noCONV_iso");
    TH1F * h_m_tau=   MakeHist(path,"Matched_TAUD_iso");
    TH1F * h_m_notau=   MakeHist(path,"Matched_noTAUD_iso");
    TH1F * h_m_ph=   MakeHist(path,"Matched_PH_iso");

    TH1F * h_nm_cf = MakeHist(path,"NonMatched_CF_iso");
    TH1F * h_nm_cfconv = MakeHist(path,"NonMatched_CF_conv_iso");
    TH1F * h_nm_cfnotconv = MakeHist(path,"NonMatched_CF_notconv_iso");
    TH1F * h_nm_notcf= MakeHist(path,"NonMatched_noCF_iso");
    TH1F * h_nm_conv=  MakeHist(path,"NonMatched_CONV_iso");
    TH1F * h_nm_convnocf =  MakeHist(path,"NonMatched_CONV_notCF_iso");
    TH1F * h_nm_noconv=  MakeHist(path,"NonMatched_noCONV_iso");
    TH1F * h_nm_tau=   MakeHist(path,"NonMatched_TAUD_iso");
    TH1F * h_nm_notau=   MakeHist(path,"NonMatched_noTAUD_iso");
    TH1F * h_nm_ph=   MakeHist(path,"NonMatched_PH_iso");
    
    
    return;
    TCanvas* c1 = new TCanvas(("Plotisoy"), "Plot", 800, 600);

    c1->SetLogy();
    
    h->Draw();
    hz_ss_conv->Draw("histsame");
    hz_ss_tau->Draw("histsame");
    hz_ss_cftau->Draw("histsame");
    hz_ss_x->Draw("histsame");

    c1->SaveAs(("dxy_ss.png"));


  }

  if(0){

    /*    
    int bin1 = hziso->FindBin(-0.05);
    int bin2 = hziso->FindBin(0.05);
    int bin3 = hziso->FindBin(0.4);
    int bin4 = hziso->FindBin(100.);
    cout << "hziso eff . = " << hziso->Integral(bin1,bin2)/ hziso->Integral(bin1,bin4) << endl;
    cout << "hcfconvxy eff . = " << hcfconvxy->Integral(bin1,bin2)/ hcfconvxy->Integral(bin1,bin4) << endl;
    cout << "hcfnonconviso eff . = " << hcfnonconviso->Integral(bin1,bin2)/ hcfnonconviso->Integral(bin1,bin4) << endl;
    cout << "h_conv_noCF eff . = " << h_conv_noCF->Integral(bin1,bin2)/ h_conv_noCF->Integral(bin1,bin4) << endl;
    cout << "h_Tau eff . = " << h_Tau->Integral(bin1,bin2)/ h_Tau->Integral(bin1,bin4) << endl;
    cout << "nonmatched" << endl;
    cout << "h_nonm_CONV eff . = " <<h_nonm_CONV->Integral(bin1,bin2)/h_nonm_CONV->Integral(bin1,bin4) << endl;
    cout << "h_nonm_CF eff . = " <<h_nonm_CF->Integral(bin1,bin2)/h_nonm_CF->Integral(bin1,bin4) << endl;
    cout << "h_nonm_TAUD eff . = " <<h_nonm_TAUD->Integral(bin1,bin2)/h_nonm_TAUD->Integral(bin1,bin4) << endl;
    
    
    cout << "h_t1 = " << h_t1->Integral(bin1,bin2)/h_t1->Integral(bin1,bin4) << endl;
    cout << "h_t2 = " << h_t2->Integral(bin1,bin2)/h_t2->Integral(bin1,bin4) << endl;
    cout << "h_t3 = " << h_t3->Integral(bin1,bin2)/h_t3->Integral(bin1,bin4) << endl;
    
    cout << "h_nt1 = " << h_nt1->Integral(bin1,bin2)/h_nt1->Integral(bin1,bin4) << endl;
    cout << "h_nt2 = " << h_nt2->Integral(bin1,bin2)/h_nt2->Integral(bin1,bin4) << endl;
    
    */
  }
    
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
  tdrStyle->SetPadBottomMargin(0.12);
  tdrStyle->SetPadLeftMargin(0.12);
  tdrStyle->SetPadRightMargin(0.1);

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

  // Change for log plots:loat_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

}
