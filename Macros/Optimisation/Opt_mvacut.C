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
#include <sstream>      // std::stringstream


//#include "Macro.h"
void setTDRStyle();


void Opt_mvacut(){
 
  
  setTDRStyle();
  //gStyle->SetPalette(1);
  
  std::vector<TString> masses;
  masses.push_back("40");
  /*masses.push_back("50");
  masses.push_back("60");
  masses.push_back("100");
  masses.push_back("200");
  masses.push_back("500");
  masses.push_back("1100");
  masses.push_back("1500");*/

  TCanvas* c1 = new TCanvas("Plot", "Plot", 800, 600);
      
  TLegend *legend = new TLegend(.7, 0.7, .9, 0.9);
  legend->SetFillColor(10);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.04);

  TString path= "/data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer//data/output/CAT/HNDiElectronOptimisation/periodBtoH/";
  TFile * fnp = new TFile(path+"HNDiElectronOptimisation_DoubleEG_SKnonprompt_hndilep_cat_v8-0-7.root");
  TFile * fcf = new TFile(path+"HNDiElectronOptimisation_SKchargeflip_dilep_cat_v8-0-7.root");
  TFile * fmc = new TFile(path+"HNDiElectronOptimisation_mc_dilep_cat_v8-0-7.root");

  TFile * fdata = new TFile(path+"HNDiElectronOptimisation_data_DoubleEG_cat_v8-0-7.root");


  vector<TString> ptbins;
  ptbins.push_back("pt1");
  ptbins.push_back("pt2");
  ptbins.push_back("pt3");
  ptbins.push_back("pt4");
  
  vector<TString> isobins;
  isobins.push_back("0.05");
  isobins.push_back("0.06");
  isobins.push_back("0.07");
  isobins.push_back("0.08");
  isobins.push_back("0.09");
  
  vector<TString> isobins_ec;
  isobins_ec.push_back("0.05");
  isobins_ec.push_back("0.06");
  isobins_ec.push_back("0.07");
  isobins_ec.push_back("0.08");
  isobins_ec.push_back("0.09");
  
  vector<TString> dxybins;
  dxybins.push_back("0.01");
  //dxybins.push_back("0.015");
  //  dxybins.push_back("0.02");
  //dxybins.push_back("0.025");
  //  dxybins.push_back("0.03");

  vector<TString> dxybins_ec;
  dxybins_ec.push_back("0.01");
  //dxybins_ec.push_back("0.015");
  //dxybins_ec.push_back("0.02");
  //dxybins_ec.push_back("0.025");
  //dxybins_ec.push_back("0.03");
  
  vector<TString> dzbins;
  dzbins.push_back("0.04");
  
  Float_t  mvabins [29] = {0., 0.1,0.2,0.3,0.4, 0.5, 0.54,0.58,0.62,0.66,0.7,0.72,0.74,0.76,0.78,0.80,0.82,0.84,0.86,0.88,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98};
  
  for(unsigned int im=0; im < masses.size(); im++){
    TFile * fsig = new TFile((path+"HNDiElectronOptimisation_HNEmEm_"+ masses.at(im)+"_cat_v8-0-7.root"));
    TFile * fsig2 = new TFile((path+"HNDiElectronOptimisation_HNEpEp_"+ masses.at(im)+"_cat_v8-0-7.root"));

    vector<TString> vmax_names;
    vector<float> vmax_bkg;
    vector<float> vmax_punzi;
    vector<float> vmax_sig_eff;
    vector<float> vmax_binx;
    vector<float> vmax_biny;


    for(unsigned int ipt=0; ipt < ptbins.size(); ipt++){
      TString bin_pt=ptbins[ipt];
      TString max_name="";
      float max_sig_eff=0.;
      float max_bkg=0.;
      float max_punzi=0.;
      float max_binx=0.;
      float max_biny=0.;
      
      
      for(unsigned int imva=1; imva < 26; imva++){
	stringstream ss;
	ss << mvabins[imva];
	
	TString  bin_mva=TString(ss.str());
	
	for(unsigned int iiso_b=0; iiso_b < isobins.size(); iiso_b++){
	  TString bin_iso_B=isobins[iiso_b];
	  for(unsigned int iiso_e=0; iiso_e < isobins_ec.size(); iiso_e++){
	    TString bin_iso_E=isobins_ec[iiso_e];
	    
	    for(unsigned int idxy_b=0; idxy_b < dxybins.size(); idxy_b++){
	      TString bin_dxy_B=dxybins[idxy_b];
	      for(unsigned int idxy_e=0; idxy_e < dxybins_ec.size(); idxy_e++){
		TString bin_dxy_E=dxybins_ec[idxy_e];
		if(idxy_e < idxy_b) continue;
		for(unsigned int idz=0; idz < dzbins.size(); idz++){
		  TString bin_dz=dzbins[idz];
		  
		  
		  TString histname="LowIDREF_"+bin_pt+"_"+bin_mva+"iso"+bin_iso_B+"_"+bin_iso_E+"_dxy"+bin_dxy_B+bin_dxy_E+"_dz"+bin_dz;
		  cout << histname << endl;
		  TH2F* h_npcutflow= (TH2F*)fnp->Get((histname));
		  TH2F* h_cfcutflow= (TH2F*)fcf->Get((histname));
		  TH2F* h_mccutflow= (TH2F*)fmc->Get((histname));
		  //TH1F* h_datacutflow= (TH1F*)fdata->Get((histname));
		  
		  TH2F* h_cutflow= (TH2F*)fsig->Get((histname));
		  TH2F* h_cutflow2= (TH2F*)fsig2->Get((histname));
		  
		  //h_cutflow->Add(h_cutflow2);
		  
		  
		  TH1F* h_cutflow_ref= (TH1F*)fsig->Get(("IDREF"));
		  
		  float incet=0.3;
		  float binnumberx(-1.);
		  float binnumbery(-1.);
		  int ibin(0);
		  float maxcont(0.);
		  float maxtot_bkg(0.);
		  float maxtot_sig(0.);
		  for(int i= 1; i < h_cutflow->GetNbinsX() +1;i++){
		    for(int j= 1; j < h_cutflow->GetNbinsY() +1;j++){
		      float sig_eff = (h_cutflow->GetBinContent(i,j) )/ h_cutflow_ref->GetBinContent(1)  ;
		      float sig_err = (h_cutflow->GetBinError(i,j) )/ h_cutflow_ref->GetBinContent(1)  ;
		      cout << h_cutflow->GetBinContent(i,j) << " " << h_cutflow_ref->GetBinContent(1) << endl;
		      if(sig_eff==0.) break;
		      float tot_bkg = h_npcutflow->GetBinContent(i,j) + h_cfcutflow->GetBinContent(i,j)+ h_mccutflow->GetBinContent(i,j);
		      //float data=h_datacutflow->GetBinContent(i,j);
		      float bkgtmp = tot_bkg + (incet*h_npcutflow->GetBinContent(i,j))*(incet*h_npcutflow->GetBinContent(i,j));
		      float denom= 1. + sqrt(bkgtmp);
		      float punzi = sig_eff/denom;
		      
		      if(punzi  > maxcont){
			maxcont=punzi;
			ibin=i;
			maxtot_bkg=tot_bkg;
			maxtot_sig=sig_eff;
			binnumberx=h_npcutflow->GetXaxis()->GetBinLowEdge(i);
			binnumbery=h_npcutflow->GetYaxis()->GetBinLowEdge(j);
			cout << "Max[init] punzi = " << maxcont << " bin = " << histname << " bkg = " << tot_bkg << " ["<< h_npcutflow->GetBinContent(i,j) << " , " << h_cfcutflow->GetBinContent(i,j) << " , " << h_mccutflow->GetBinContent(i,j)<<"] "<< " sig = " << sig_eff << " +- " << sig_err <<" mva2 = " << binnumberx << " mva3= " << binnumbery << " data=" << 0 << endl;
			  
			
		      }
		      //h_mass_isoopt->SetBinContent(i+1,punzi);
		    } 
		  }
		  if(maxcont > max_punzi){
		    max_punzi=maxcont;
		    max_name=histname;
		    max_bkg=maxtot_bkg;
		    max_sig_eff=maxtot_sig;
		    max_binx=binnumberx;
		    max_biny=binnumbery;
		    cout << "Max punzi = " << maxcont << " bin = " << max_name << " bkg = " << max_bkg << " sig = " << max_sig_eff << " mva2 = " << max_binx << " mva3= " << max_biny << endl;
		  }
		  //h_mass_isoopt->GetYaxis()->SetTitle("Punzi"); 
		  //h_mass_isoopt->GetXaxis()->SetTitle("d_{xy} cut");
		  /*
		  h_mass_isoopt->SetLineWidth(3.);
		  if(im==0)  h_mass_isoopt->SetLineColor(kRed);
		  if(im==1)  h_mass_isoopt->SetLineColor(kBlue);
		  if(im==2)  h_mass_isoopt->SetLineColor(kGreen+4);
		  if(im==3)  h_mass_isoopt->SetLineColor(kCyan);
		  if(im==4)  h_mass_isoopt->SetLineColor(kSpring-2);
		  if(im==5)  h_mass_isoopt->SetLineColor(kOrange);
		  
		  if(im==4)  h_mass_isoopt->SetLineStyle(2.);
		  if(im==5)  h_mass_isoopt->SetLineStyle(4.);    
		  
		  //h_mass_isoopt->GetYaxis()->SetRangeUser(0.001, 0.005);
		  h_mass_isoopt->GetYaxis()->SetTitleSize(0.05);
		  h_mass_isoopt->GetYaxis()->SetTitleOffset(1.5);
		  h_mass_isoopt->GetYaxis()->SetLabelSize(0.04);
		  
		  h_mass_isoopt->GetXaxis()->SetTitleSize(0.05);
		  h_mass_isoopt->GetXaxis()->SetLabelSize(0.04);
		  
		  if(im==0)h_mass_isoopt->Draw("hist");
		  else h_mass_isoopt->Draw("histsame");
    
		  //legend->AddEntry(h_mass_isoopt,("m_{N}= "+ masses.at(im)).Data(),"l");
		  */
		}
	      }
	    }
	  }
	}
      }
    
      vmax_names.push_back(max_name);
      vmax_punzi.push_back(max_punzi);
      vmax_bkg.push_back(max_bkg);
      vmax_sig_eff.push_back(max_sig_eff);
      vmax_binx.push_back(max_binx);
      vmax_biny.push_back(max_biny);
    }
    cout << "Max" << endl;
    for(int i=0; i < vmax_names.size(); i++){
      cout << vmax_names[i] <<  " " << vmax_bkg[i] << " " <<  vmax_sig_eff[i]<<  " max punzi=" << vmax_punzi[i] << " mva2= " << vmax_binx[i] <<  " mva3= " << vmax_biny[i] << endl;
    }
    
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
