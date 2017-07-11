#include <string.h>
#include "TGaxis.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TKey.h"
#include <iostream>
#include <sstream>
#include <TStyle.h>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

//#include "Macro.h"
void setTDRStyle();


void Opt_id_select(){
 
  
  setTDRStyle();
  gStyle->SetPalette(1);
  
  std::vector<TString> masses;
  masses.push_back("40");
  masses.push_back("60");
  masses.push_back("100");


  TCanvas* c1 = new TCanvas("Plot", "Plot", 800, 600);
      
  TLegend *legend = new TLegend(.4, 0.4, .6, 0.7);
  legend->SetFillColor(10);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.04);
  //legend->SetTitle("Signal Eff.");

  TLegend *legend2 = new TLegend(.6, 0.7, .8, 0.9);
  legend2->SetFillColor(10);
  legend2->SetBorderSize(0);
  legend2->SetTextSize(0.04);


  vector<float> vcut_dxy_b;
  vector<TString> vcut_dxy_b_s;

  for(unsigned int dxy_b=0;dxy_b < 1; dxy_b++){
    float cut_dxy_b =  float(dxy_b)*0.005 + 0.005;
    vcut_dxy_b.push_back(cut_dxy_b);
    stringstream ss;
    ss <<cut_dxy_b;
    vcut_dxy_b_s.push_back(TString(ss.str()));
  }
  vector<float> vcut_dxysig_b;
  vcut_dxysig_b.push_back(3.);
  vcut_dxysig_b.push_back(4.);
  vector<TString> vcut_dxysig_b_s;
  vcut_dxysig_b_s.push_back("3");
  vcut_dxysig_b_s.push_back("4");

  vector<float> vcut_dz_b;
  vector<TString> vcut_dz_b_s;

  for(unsigned int dz_b=0;dz_b < 4; dz_b++){
    float cut_dz_b =  float(dz_b)*0.02 + 0.02;
    vcut_dz_b.push_back(cut_dz_b);
    stringstream ss;
    ss <<cut_dz_b;
    vcut_dz_b_s.push_back(TString(ss.str()));
  }

  vector<float> vcut_iso_b;
  vector<TString> vcut_iso_b_s;
  for(unsigned int iso_b=0;iso_b < 5; iso_b++){
    float cut_iso_b = float(iso_b)*0.01 + 0.05;
    vcut_iso_b.push_back(cut_iso_b);
    stringstream ss;
    ss <<cut_iso_b;
    vcut_iso_b_s.push_back(TString(ss.str()));
  }
  
  vector<TString> cuts;
  int ncut=0;
  for(unsigned int dxy_b=0; dxy_b < vcut_dxy_b.size(); dxy_b++){
    for(unsigned int dxysig_b=0; dxysig_b < vcut_dxysig_b.size(); dxysig_b++){
      for(unsigned int dz_b=0; dz_b < vcut_dz_b.size(); dz_b++){
	for(unsigned int iso_b=0; iso_b < vcut_iso_b.size(); iso_b++){
	  for(unsigned int dz_e=0; dz_e < vcut_dz_b.size(); dz_e++){
	    for(unsigned int iso_e=0; iso_e < vcut_iso_b.size(); iso_e++, ncut++){

	      cuts.push_back(vcut_iso_b_s[iso_b]+"_"+vcut_dxy_b_s[dxy_b]+"_"+vcut_dxysig_b_s[dxysig_b]+"_"+vcut_dz_b_s[dz_b]+ "_"+vcut_dz_b_s[dz_e] + vcut_iso_b_s[iso_e]);
	    }
	  }
	}
      }
    }
  }
  

  TString path= "/afs/cern.ch/work/j/jalmond/CAT/HNDiMuonOptimisation/periodBtoH/";
  TFile * fnp = new TFile(path+"HNDiMuonOptimisation_DoubleMuon_SKnonprompt_hndilep_cat_v8-0-7.root");
  TFile * fmc = new TFile(path+"HNDiMuonOptimisation_MC_trilep_cat_v8-0-7.root");

  TFile * fdata = new TFile(path+"HNDiMuonOptimisation_data_DoubleMuon_cat_v8-0-7.root");

  TH1F* h_npcutflow= (TH1F*)fnp->Get(("LowIDREFSNUMedium_double"));
  TH1F* h_mccutflow= (TH1F*)fmc->Get(("LowIDREFSNUMedium_double"));

  TH1F* h_datacutflow= (TH1F*)fdata->Get(("LowIDREFSNUMedium_double"));


  double FIGURE2_RATIO = 0.35;
  double SUBFIGURE_MARGIN = 0.15;
  c1->SetBottomMargin(FIGURE2_RATIO);
  TPad *p = new TPad( "p_test", "", 0, 0, 1, 1.0 - SUBFIGURE_MARGIN, 0, 0, 0);  // create new pad, fullsize to have equal font-sizes in both plots                        
 

  for(unsigned int im=0; im < masses.size(); im++){
    TFile * fsig = new TFile((path+"HNDiMuonOptimisation_HNMumMum_"+masses[im]+"_cat_v8-0-7.root").Data());
    TFile * fsig2 = new TFile((path+"HNDiMuonOptimisation_HNMupMup_"+masses[im]+"_cat_v8-0-7.root").Data());
    
    
    cout << "\n ----------------" << endl;
    cout <<  "New mass : " << masses.at(im) << endl;

    TH1F* h_cutflowref= (TH1F*)fsig->Get(("POGIDREFSNU_double")); 
    TH1F* h_cutflowref2= (TH1F*)fsig->Get(("POGIDREFSNU_double")); 
    h_cutflowref->Add(h_cutflowref2);
    TH1F* h_cutflow= (TH1F*)fsig->Get(("LowIDREFSNUMedium_double"));
    TH1F* h_cutflow2= (TH1F*)fsig2->Get(("LowIDREFSNUMedium_double"));

    h_cutflow->Add(h_cutflow2);

    for(int i=1; i < h_cutflow->GetNbinsX()+1; i++){
      cout <<i << " : " <<  h_cutflow->GetBinContent(i) << endl;
    }

    TH1F* h_mass_d0opt = new TH1F(("IDcut_opt_HN"+ masses.at(im)).Data(),("IDcut_opt_HN"+ masses.at(im)).Data(),900,0.,900.);
    
    TH1F* h_mass_data = new TH1F(("DataIDcut_opt_HN"+ masses.at(im)).Data(),("DataIDcut_opt_HN"+ masses.at(im)).Data(),900,0.,900.);
    
    TH1F* h_mass_eff = new TH1F(("EffIDcut_opt_HN"+ masses.at(im)).Data(),("EffIDcut_opt_HN"+ masses.at(im)).Data(),900,0.,900.);
    
    TH1F* h_mass_bkg = new TH1F(("BkgIDcut_opt_HN"+ masses.at(im)).Data(),("BkgIDcut_opt_HN"+ masses.at(im)).Data(),900,0.,900.);



    float punzi_max(0.);
    float s_rb_max(0.);
    int i_punzi(-99);
    int i_sb(-999);
    for(int i= 1; i < 810; i++){
      float sig_eff = h_cutflow->GetBinContent(i)/h_cutflowref->GetBinContent(2);
      cout << "Sig eff. = " << sig_eff << endl;
      float tot_bkg = h_npcutflow->GetBinContent(i) +  h_mccutflow->GetBinContent(i);
      cout << "Fakes " << h_npcutflow->GetBinContent(i) <<" MC: " << h_mccutflow->GetBinContent(i) << endl;
      float bkgtmp = tot_bkg + (0.3*h_npcutflow->GetBinContent(i))*(0.3*h_npcutflow->GetBinContent(i));
      float denom= 1. + sqrt(bkgtmp);
      cout << "Data = " << h_datacutflow->GetBinContent(i) << endl;

      float punzi = sig_eff/denom;
      float s_rb = sig_eff/sqrt(tot_bkg);
      cout << "punzi = " << punzi << endl;
      if(punzi_max < punzi) {
	punzi_max=punzi;
	i_punzi = i;
      }
      if(s_rb_max < s_rb){
	s_rb_max=s_rb;
	i_sb=i;
      }
      h_mass_d0opt->SetBinContent(i,punzi);
      h_mass_eff->SetBinContent(i,sig_eff);
      h_mass_bkg->SetBinContent(i,tot_bkg);
      h_mass_data->SetBinContent(i, h_datacutflow->GetBinContent(i));
    } 
    
    cout << "punzi_max = " << punzi_max << " id = " << cuts[i_punzi] << endl;
    cout << "s/sqrtb_max = " << s_rb_max << " id = " << cuts[i_sb] << endl;
    
    //h_mass_eff->GetXaxis()->SetTitle("ID");

    h_mass_d0opt->GetXaxis()->SetTitle("ID"); 
    h_mass_d0opt->GetYaxis()->SetTitle("Punzi"); 
    h_mass_d0opt->SetLineWidth(3.);
    //h_mass_d0opt->GetYaxis()->SetRangeUser(0.,0.004);
    if(im==0)  h_mass_d0opt->SetLineColor(kRed);
    if(im==1)  h_mass_d0opt->SetLineColor(kBlue);
    if(im==2)  h_mass_d0opt->SetLineColor(kGreen+4);

    h_mass_eff->SetLineWidth(3.);
    if(im==0)  h_mass_eff->SetLineColor(kRed);
    if(im==1)  h_mass_eff->SetLineColor(kBlue);
    if(im==2)  h_mass_eff->SetLineColor(kGreen+4);

    
    h_mass_eff->GetYaxis()->SetRangeUser(0.001, 2.);
    h_mass_eff->GetYaxis()->SetTitleSize(0.05);
    h_mass_eff->GetYaxis()->SetTitleOffset(0.7);
    h_mass_eff->GetYaxis()->SetLabelSize(0.05);
    
    h_mass_eff->GetXaxis()->SetTitleSize(0.05);
    h_mass_eff->GetXaxis()->SetLabelSize(0.03);
    h_mass_eff->GetYaxis()->SetTitle("Signal Eff.");

    h_mass_eff->GetXaxis()->SetLabelSize(0.); 
    h_mass_eff->GetYaxis()->SetLabelSize(0.04); 
    h_mass_d0opt->GetYaxis()->SetLabelSize(0.03); 
    h_mass_d0opt->GetYaxis()->SetTitle("Punzi");
    h_mass_d0opt->GetYaxis()->SetTitleSize(0.05);
    h_mass_d0opt->GetYaxis()->SetTitleOffset(1.);
    if(im==0){

      h_mass_eff->Draw("hist");
      
      p->SetTopMargin(1-FIGURE2_RATIO);   // top-boundary (should be 1 - thePad->GetBottomMargin() )
      p->SetFillStyle(0);     // needs to be transparent
      p->Draw();
      p->cd();
      
      //if(1.1*h_mass_d0opt->GetMaximum() > gPad->GetUymax()) h_mass_d0opt->GetYaxis()->SetRangeUser(0., h_mass_d0opt->GetMaximum()*1.1);
      h_mass_d0opt->GetYaxis()->SetRangeUser(0., 0.01);

      h_mass_d0opt->Draw("hist");
      c1->cd();
      c1->Update();
    }
    else{
      h_mass_eff->Draw("histsame");
      p->cd();
      h_mass_d0opt->Draw("histsame");
      c1->cd();
    }
    c1->Update();

    if(im==2){
      //scale hint1 to the pad coordinates
      Float_t rightmax = 1.1*h_mass_bkg->GetMaximum();
      Float_t scale = gPad->GetUymax()/rightmax;
      h_mass_bkg->SetLineWidth(4.);
      h_mass_bkg->SetLineStyle(3.);
      h_mass_bkg->SetLineColor(kRed);
      h_mass_bkg->Scale(scale);
      
      h_mass_bkg->Draw("histsame");
      h_mass_data->SetMarkerColor(kRed);

      h_mass_data->Scale(scale);
      //h_mass_data->Draw("psame");
      cout << "rightmax = " << rightmax << endl;
      
      TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
				gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
      
      axis->SetLineColor(kRed);
      axis->SetLabelColor(kRed);
      axis->SetTitleOffset(1.);
      axis->SetLabelSize(0.03);
      axis->SetTitle("#Bkg");
      axis->SetTitleColor(kRed);
      axis->Draw();
      
      
      
    }
    
    legend->AddEntry(h_mass_d0opt,("Eff. m_{N}= "+ masses.at(im)+ " GeV").Data(),"l");
    if(im==2)legend->AddEntry(h_mass_bkg,"N_{Bkg}","l");
  }
  
  legend->Draw("same");
  //legend2->Draw("same");
  
  c1->SaveAs(("IDmmoptimise_snu_numbers.pdf"));

  
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
  tdrStyle->SetPadLeftMargin(0.1);
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
  tdrStyle->SetPadTickX(0);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(0);

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



