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
#include "TLegend.h"
#include "TLatex.h"

//#include "Macro.h"
void setTDRStyle();
TH1* makehist(TString histname, TString type,TString path, TString sample, TString cut, TString label );


void ElectronTypeStack(){
 
  
  setTDRStyle();
  
  std::vector<TString> vlabel;
  vlabel.push_back("Z");
  vlabel.push_back("W");
  vlabel.push_back("ZorW");
  std::vector<TString> Types;
  std::vector<TString> etypes;
  Types.push_back("notprompt_ph_mother_l50");
  Types.push_back("notprompt_ph_mother_g50");
  Types.push_back("_notphotonmatched_notprompt");
  Types.push_back("_cf");
  Types.push_back("_cf_conv");
  Types.push_back("_notcf_conv");
  Types.push_back("_prompt");
  Types.push_back("_notcf_notprompt_conv");
  Types.push_back("_notcf_notprompt_notconv");
  
  
  for(int j=0; j < 0; j++){
    TString sel="";
    if (j==1) sel="_mme";
    
    for(int i=0; i < 3; i++){
      for(unsigned int it =0 ; it < Types.size() ; it++){
	etypes.push_back(vlabel[i] + Types[it]+ sel);
      }
    }
    
    etypes.push_back("matched_nottau_notphoton_notfromZ"+sel);
    etypes.push_back("matched_nottau_notfromZ"+sel);
    etypes.push_back("matched_tau"+sel);
    etypes.push_back("nonmatched_tau"+sel);
    etypes.push_back("nonmatched_photon_pdgidl50"+sel);
    etypes.push_back("nonmatched_photon_pdgidg50"+sel);
    etypes.push_back("nonmatched_nottau_notphoton"+sel);  
    etypes.push_back("nm_cf"+sel);  
    etypes.push_back("nm_cf_conv"+sel);  
    etypes.push_back("nm_cf_notconv"+sel);  
    etypes.push_back("nm_photon_notcf_conv"+sel);  
    etypes.push_back("nm_notphoton_notcf_conv"+sel);  
    etypes.push_back("nm_notprompt_photon"+sel);  
    etypes.push_back("nm_notprompt"+sel);  
  }
  
  
  etypes.push_back("PROMPT");
  etypes.push_back("NP_mpdgid1");
  etypes.push_back("NP_mpdgid2");
  etypes.push_back("NP_mpdgid3");
  etypes.push_back("NP_mpdgid4");
  etypes.push_back("NP_mpdgid5");
  etypes.push_back("NP_mpdgid6");
  etypes.push_back("NP_mpdgid15");
  etypes.push_back("NP_mpdgid22");


  for (int pl =0 ; pl < etypes.size(); pl++){
    
    
    TLegend *legend = new TLegend(.7, 0.7, .9, 0.9);
    legend->SetFillColor(10);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    
    TString path= "/data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer//data/output/CAT/ElectronTypes/periodBtoH/";
    
    //makehist(etypes[pl] + sample +"_eta", etypes[pl], path,  sample, "_eta" , "#eta");
    
    std::vector<TString> samples;
    samples.push_back("DYJets");
    samples.push_back("QCD_DoubleEMEnriched_30-inf_mgg40to80");
    samples.push_back("qcd_30to80_bctoe");
    samples.push_back("WJets");
    samples.push_back("ZGto2LG");
    samples.push_back("TT_powheg");

    TCanvas* c1 = new TCanvas("Plot", "Plot", 800, 600);
    
    TH1* hdy = makehist(etypes[pl] + samples[0] +"_dxy", etypes[pl], path,  samples[0], "_dxy" , "dxy");
    if(!hdy) continue;
    hdy->Draw();
    c1->SetLogy();
    for(int s=1; s< samples.size() ; s++){
      TH1* h = makehist(etypes[pl] + samples[s] +"_dxy", etypes[pl], path,  samples[s], "_dxy" , "dxy");
      if(!h) continue;
      if(s==1) h->SetLineColor(kRed);
      if(s==2) h->SetLineColor(kCyan);
      if(s==3) h->SetLineColor(kGreen);
      if(s==4) h->SetLineColor(kBlue);
      if(s==5) h->SetLineColor(kOrange);
      h->SetLineWidth(2.);
      h->Draw("histsame");
      legend->AddEntry(h, samples[s], "l");
    }
    legend->Draw();
    c1->SaveAs("~/ElectronTypes/" + etypes[pl] + "_dxy.pdf");
    
    continue;
          
  }
}


TH1* makehist(TString histname, TString type,TString path, TString sample, TString cut, TString label ){

  TFile * fmc = new TFile(path+"ElectronTypes_"+sample+"_cat_v8-0-6.root");

  TH1F* h= (TH1F*)fmc->Get(( type + cut));

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
