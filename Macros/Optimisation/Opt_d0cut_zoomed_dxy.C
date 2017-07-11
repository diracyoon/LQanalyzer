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


void Opt_d0cut(){
 
  
  setTDRStyle();
  //gStyle->SetPalette(1);
  
  std::vector<TString> masses;
  masses.push_back("40");
  masses.push_back("50");
  masses.push_back("60");
  masses.push_back("100");
  masses.push_back("200");
  masses.push_back("500");
  masses.push_back("1100");
  masses.push_back("1500");

  TCanvas* c1 = new TCanvas("Plot", "Plot", 800, 600);
      
  TLegend *legend = new TLegend(.7, 0.7, .9, 0.9);
  legend->SetFillColor(10);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.04);

  TString path= "/data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer//data/output/CAT/HNDiElectronOptimisation/periodBtoH/";
  TFile * fnp = new TFile(path+"HNDiElectronOptimisation_DoubleEG_SKnonprompt_dilep_cat_v8-0-6.root");
  TFile * fcf = new TFile(path+"HNDiElectronOptimisation_SKchargeflip_dilep_cat_v8-0-6.root");
  TFile * fmc = new TFile(path+"HNDiElectronOptimisation_mc_dilep_cat_v8-0-6.root");


  TH1F* h_npcutflow= (TH1F*)fnp->Get(("d0sigcutflowzoomed"));
  TH1F* h_cfcutflow= (TH1F*)fcf->Get(("d0sigcutflowzoomed"));
  TH1F* h_mccutflow= (TH1F*)fmc->Get(("d0sigcutflowzoomed"));

  for(unsigned int im=0; im < masses.size(); im++){
    TFile * fsig = new TFile((path+"HNDiElectronOptimisation_HNEmEm_"+ masses.at(im)+"_cat_v8-0-6.root"));
    TFile * fsig2 = new TFile((path+"HNDiElectronOptimisation_HNEpEp_"+ masses.at(im)+"_cat_v8-0-6.root"));

    
    cout << "\n ----------------" << endl;
    cout <<  "New mass : " << masses.at(im) << endl;
    TH1F* h_cutflowref= (TH1F*)fsig->Get(("d0sigcutflow"));
    TH1F* h_cutflowref2= (TH1F*)fsig2->Get(("d0sigcutflow"));
    
    h_cutflowref->Add(h_cutflowref2);
    
    TH1F* h_cutflow= (TH1F*)fsig->Get(("d0cutflowzoomed"));
    TH1F* h_cutflow2= (TH1F*)fsig2->Get(("d0cutflowzoomed"));

    h_cutflow->Add(h_cutflow2);


    vector<TString> binlabels;


    std::vector<TString> dxyb2;
    dxyb2.push_back("017");
    dxyb2.push_back("015");
    dxyb2.push_back("014");
    dxyb2.push_back("013");
    dxyb2.push_back("012");
    dxyb2.push_back("011");
    dxyb2.push_back("010");
    std::vector<TString> dxye2;
    dxye2.push_back("017");
    dxye2.push_back("015");
    dxye2.push_back("014");
    dxye2.push_back("013");
    dxye2.push_back("012");
    dxye2.push_back("011");
    dxye2.push_back("010");

    


    TH1F* h_mass_d0opt = new TH1F(("d0cut_opt_HN"+ masses.at(im)).Data(),("d0cut_opt_HN"+ masses.at(im)).Data(), 55,0.,55.);
    if(1){
      int idxy2(1);
      for(unsigned int ib=0; ib <dxyb2.size(); ib++){
        for(unsigned int ie=0; ie <dxye2.size(); ie++,idxy2++){
	  h_mass_d0opt->GetXaxis()->SetBinLabel(idxy2,dxyb2[ib]+"_"+dxye2[ie]);
	  binlabels.push_back(dxyb2[ib]+"_"+dxye2[ie]);
        }
      } 
    }

    
    int ibin(0);
    float maxcont(0.);
    for(int i= 0; i < 55;i++){
      float sig_eff = h_cutflow->GetBinContent(i+2)/ h_cutflowref->GetBinContent(1);
      //float tot_bkg =  h_npcutflow->GetBinContent(i+2) + h_mccutflow->GetBinContent(i+2);           
      cout << i+1 << " " << h_cutflow->GetBinContent(i+2) << " " << h_cutflowref->GetBinContent(1)<< " "<<  sig_eff << " " << h_npcutflow->GetBinContent(i+2) << " " << h_cfcutflow->GetBinContent(i+2) << " " << h_mccutflow->GetBinContent(i+2) << endl;
      float tot_bkg = h_npcutflow->GetBinContent(i+2) + h_cfcutflow->GetBinContent(i+2)+ h_mccutflow->GetBinContent(i+2);
      float bkgtmp = tot_bkg + (0.3*h_npcutflow->GetBinContent(i+2))*(0.3*h_npcutflow->GetBinContent(i+2));
      float denom= 1. + sqrt(bkgtmp);
      
      float punzi = sig_eff/denom;
      if(punzi  > maxcont){
        maxcont=punzi;
	ibin=i;
      }
      cout << punzi << endl;
      h_mass_d0opt->SetBinContent(i+1,punzi);
    } 
    cout << "Max punzi = " << maxcont << " bin = " << binlabels[ibin] << endl;
    h_mass_d0opt->GetYaxis()->SetTitle("Punzi"); 
    h_mass_d0opt->GetXaxis()->SetTitle("d_{xy} cut");
    
    h_mass_d0opt->SetLineWidth(3.);
    if(im==0)  h_mass_d0opt->SetLineColor(kRed);
    if(im==1)  h_mass_d0opt->SetLineColor(kBlue);
    if(im==2)  h_mass_d0opt->SetLineColor(kGreen+4);
    if(im==3)  h_mass_d0opt->SetLineColor(kCyan);
    if(im==4)  h_mass_d0opt->SetLineColor(kSpring-2);
    if(im==5)  h_mass_d0opt->SetLineColor(kOrange);

    if(im==4)  h_mass_d0opt->SetLineStyle(2.);
    if(im==5)  h_mass_d0opt->SetLineStyle(4.);    

    h_mass_d0opt->GetYaxis()->SetRangeUser(0.001, 0.05);
    h_mass_d0opt->GetYaxis()->SetTitleSize(0.05);
    h_mass_d0opt->GetYaxis()->SetTitleOffset(1.5);
    h_mass_d0opt->GetYaxis()->SetLabelSize(0.04);
    
    h_mass_d0opt->GetXaxis()->SetTitleSize(0.05);
    h_mass_d0opt->GetXaxis()->SetLabelSize(0.04);

    if(im==0)h_mass_d0opt->Draw("hist");
    else h_mass_d0opt->Draw("histsame");
    
    legend->AddEntry(h_mass_d0opt,("m_{N}= "+ masses.at(im)).Data(),"l");
  }
  
  legend->Draw("same");
  c1->SaveAs(("d0sigopt.pdf"));
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
