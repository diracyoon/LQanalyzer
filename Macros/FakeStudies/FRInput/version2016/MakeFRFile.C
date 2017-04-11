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

#include "TString.h"

void setTDRStyle();
bool CheckFile(TFile* f);
bool CheckHist(TH2* h);


void MakeFRRootFile(){
  
  TString path= "/data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer//data/output/CAT/FakeRateCalculator_El/periodBtoH/";

  TFile * fdata = new TFile(path + "FakeRateCalculator_El_data_DoubleEG_cat_v8-0-6.root");
  TFile * fmc = new TFile(path + "FakeRateCalculator_El_mc_v8-0-6.root");
  cout << path + "FakeRateCalculator_El_data_DoubleEG_cat_v8-0-6.root" << endl;
  if(!fdata)cout << "No Data" << endl;
  if (!fmc) cout << "No MC" << endl;
  /// Set Plotting style
  setTDRStyle();
  gStyle->SetPalette(1);
    
  TString outfile = "FakeRate13TeV_2016_opt.root";
  TFile* fout = new TFile(outfile.Data(),"RECREATE");
  fout->cd();

  std::vector<TString> fakes40;
  fakes40.push_back("20_pt_eta");
  fakes40.push_back("60_pt_eta");
  fakes40.push_back("40_pt_eta");
  fakes40.push_back("30_pt_eta");

  std::vector<TString> isocut;
  isocut.push_back("dxy_b050_e100");
  isocut.push_back("dxy_b050_e050");
  isocut.push_back("dxy_b050_e040");
  isocut.push_back("dxy_b050_e025");
  isocut.push_back("dxy_b050_e020");
  isocut.push_back("dxy_b050_e010");

  isocut.push_back("dxy_b025_e100");
  isocut.push_back("dxy_b025_e050");
  isocut.push_back("dxy_b025_e040");
  isocut.push_back("dxy_b025_e025");
  isocut.push_back("dxy_b025_e020");
  isocut.push_back("dxy_b025_e010");

  isocut.push_back("dxy_b015_e100");
  isocut.push_back("dxy_b015_e050");
  isocut.push_back("dxy_b015_e040");
  isocut.push_back("dxy_b015_e025");
  isocut.push_back("dxy_b015_e020");
  isocut.push_back("dxy_b015_e010");

  isocut.push_back("dxy_b010_e100");
  isocut.push_back("dxy_b010_e050");
  isocut.push_back("dxy_b010_e040");
  isocut.push_back("dxy_b010_e025");
  isocut.push_back("dxy_b010_e020");
  isocut.push_back("dxy_b010_e010");
  
  isocut.push_back("dxy_b017_e017");
  isocut.push_back("dxy_b017_e015");
  isocut.push_back("dxy_b017_e014");
  isocut.push_back("dxy_b017_e013");
  isocut.push_back("dxy_b017_e012");
  isocut.push_back("dxy_b017_e011");
  isocut.push_back("dxy_b017_e010");

  isocut.push_back("dxy_b015_e017");
  isocut.push_back("dxy_b015_e015");
  isocut.push_back("dxy_b015_e014");
  isocut.push_back("dxy_b015_e013");
  isocut.push_back("dxy_b015_e012");
  isocut.push_back("dxy_b015_e011");
  isocut.push_back("dxy_b015_e010");

  isocut.push_back("dxy_b014_e017");
  isocut.push_back("dxy_b014_e015");
  isocut.push_back("dxy_b014_e014");
  isocut.push_back("dxy_b014_e013");
  isocut.push_back("dxy_b014_e012");
  isocut.push_back("dxy_b014_e011");
  isocut.push_back("dxy_b014_e010");

  isocut.push_back("dxy_b013_e017");
  isocut.push_back("dxy_b013_e015");
  isocut.push_back("dxy_b013_e014");
  isocut.push_back("dxy_b013_e013");
  isocut.push_back("dxy_b013_e012");
  isocut.push_back("dxy_b013_e011");
  isocut.push_back("dxy_b013_e010");

  isocut.push_back("dxy_b012_e017");
  isocut.push_back("dxy_b012_e015");
  isocut.push_back("dxy_b012_e014");
  isocut.push_back("dxy_b012_e013");
  isocut.push_back("dxy_b012_e012");
  isocut.push_back("dxy_b012_e011");
  isocut.push_back("dxy_b012_e010");

  isocut.push_back("dxy_b011_e017");
  isocut.push_back("dxy_b011_e015");
  isocut.push_back("dxy_b011_e014");
  isocut.push_back("dxy_b011_e013");
  isocut.push_back("dxy_b011_e012");
  isocut.push_back("dxy_b011_e011");
  isocut.push_back("dxy_b011_e010");

  isocut.push_back("dxy_b010_e017");
  isocut.push_back("dxy_b010_e015");
  isocut.push_back("dxy_b010_e014");
  isocut.push_back("dxy_b010_e013");
  isocut.push_back("dxy_b010_e012");
  isocut.push_back("dxy_b010_e011");


  isocut.push_back("b050_e050");
  isocut.push_back("b050_e0525");
  isocut.push_back("b050_e055");
  isocut.push_back("b050_e060");
  isocut.push_back("b050_e065");
  isocut.push_back("b050_e075");
  isocut.push_back("b050_e100");
  isocut.push_back("b050_e125");

  isocut.push_back("b0525_e050");
  isocut.push_back("b0525_e0525");
  isocut.push_back("b0525_e055");
  isocut.push_back("b0525_e060");
  isocut.push_back("b0525_e065");
  isocut.push_back("b0525_e075");
  isocut.push_back("b0525_e100");
  isocut.push_back("b0525_e125");

  isocut.push_back("b055_e050");
  isocut.push_back("b055_e0525");
  isocut.push_back("b055_e055");
  isocut.push_back("b055_e060");
  isocut.push_back("b055_e065");
  isocut.push_back("b055_e075");
  isocut.push_back("b055_e100");
  isocut.push_back("b055_e125");

  isocut.push_back("b065_e050");
  isocut.push_back("b065_e0525");
  isocut.push_back("b065_e055");
  isocut.push_back("b065_e060");
  isocut.push_back("b065_e065");
  isocut.push_back("b065_e075");
  isocut.push_back("b065_e100");
  isocut.push_back("b065_e125");

  isocut.push_back("b060_e050");
  isocut.push_back("b060_e0525");
  isocut.push_back("b060_e055");
  isocut.push_back("b060_e060");
  isocut.push_back("b060_e065");
  isocut.push_back("b060_e075");
  isocut.push_back("b060_e100");
  isocut.push_back("b060_e125");

  isocut.push_back("b075_e050");
  isocut.push_back("b075_e0525");
  isocut.push_back("b075_e055");
  isocut.push_back("b075_e060");
  isocut.push_back("b075_e065");
  isocut.push_back("b075_e075");
  isocut.push_back("b075_e100");
  isocut.push_back("b075_e125");

  isocut.push_back("b100_e050");
  isocut.push_back("b100_e0525");
  isocut.push_back("b100_e055");
  isocut.push_back("b100_e060");
  isocut.push_back("b100_e065");
  isocut.push_back("b100_e075");
  isocut.push_back("b100_e100");
  isocut.push_back("b100_e125");

  isocut.push_back("b125_e050");
  isocut.push_back("b125_e0525");
  isocut.push_back("b125_e055");
  isocut.push_back("b125_e060");
  isocut.push_back("b125_e065");
  isocut.push_back("b125_e075");
  isocut.push_back("b125_e100");
  isocut.push_back("b125_e125");


  std::vector<TString> fakes;
  fakes.push_back("ELECTRON16_HN_TIGHT_dijet_d0");
  fakes.push_back("ELECTRON16_HN_TIGHT_DXYSIG_dijet_d0_dxysig");
  fakes.push_back("ELECTRON16_HN_TIGHT_dijet_nod0");
  fakes.push_back("ELECTRON16_HN_TIGHT_DXYSIG_dijet_nod0_dxysig");
  fakes.push_back("ELECTRON16_HN_TIGHT_dijet_nod0_B");
  fakes.push_back("ELECTRON16_HN_TIGHT_dijet_nod0_C");
  fakes.push_back("ELECTRON16_HN_TIGHT_dijet_nod0_D");
  fakes.push_back("ELECTRON16_HN_TIGHT_dijet_nod0_E");
  fakes.push_back("ELECTRON16_HN_TIGHT_dijet_nod0_F");
  fakes.push_back("ELECTRON16_HN_TIGHT_dijet_nod0_G");
  fakes.push_back("ELECTRON16_HN_TIGHT_dijet_nod0_H");
  fakes.push_back("ELECTRON16_HN_TIGHT_dijet_iso04");
  fakes.push_back("ELECTRON16_HN_TIGHT_dijet_iso06");
  fakes.push_back("ELECTRON16_FR_POG_TIGHT_dijet_pog");
  fakes.push_back("ELECTRON16_FR_POG_MEDIUM_dijet_pog");
  fakes.push_back("ELECTRON16_FR_POG_TIGHT_CC_dijet_pog");
  fakes.push_back("ELECTRON16_FR_POG_MEDIUM_CC_dijet_pog");
  fakes.push_back("ELECTRON16_FR_MVA_TIGHT_CC_dijet_mva");
  fakes.push_back("ELECTRON16_FR_POG_TIGHT_DXYCC_dijet_pog");
  fakes.push_back("ELECTRON16_FR_POG_MEDIUM_DXYCC_dijet_pog");
  fakes.push_back("ELECTRON16_FR_MVA_TIGHT_DXYCC_dijet_mva");

  fakes.push_back("ELECTRON16_FR_POG_TIGHT_CC_dijet_pog_d0");
  fakes.push_back("ELECTRON16_FR_POG_MEDIUM_CC_dijet_pog_d0");
  fakes.push_back("ELECTRON16_FR_MVA_TIGHT_CC_dijet_mva_d0");
  fakes.push_back("ELECTRON16_FR_POG_TIGHT_DXYCC_dijet_pog_d0");
  fakes.push_back("ELECTRON16_FR_POG_MEDIUM_DXYCC_dijet_pog_d0");
  fakes.push_back("ELECTRON16_FR_MVA_TIGHT_DXYCC_dijet_mva_d0");


  std::vector<TString> fakes_opt;
  fakes_opt.push_back("HNTight_");
  //fakes_opt.push_back("HNTight_miniiso_");
  fakes_opt.push_back("HNTight_dxysig_");
  fakes_opt.push_back("HNTight_miniiso_dxysig_");
  
  for(vector<TString>::iterator it2 = fakes40.begin(); it2!=fakes40.end(); ++it2){
    for(vector<TString>::iterator it = fakes.begin(); it!=fakes.end(); ++it){
      cout << *it2 << endl;
      if(!CheckFile(fdata))return;
     if(!CheckFile(fmc))return;

      TString denom ="LooseEl" + *it + "_"+ *it2;
      TString num ="TightEl" + *it +  "_"+ *it2;
      //if (!denom.Contains("0")){
      //denom ="LooseEl" + *it+ *it2;
      //num ="TightEl" + *it +  *it2;
      
      //}
      TH2D* h_pt_num= (TH2D*)fdata->Get(num.Data());
      TH2D* h_pt_denom= (TH2D*)fdata->Get(denom.Data());
      
      
      cout << h_pt_num << " " << h_pt_denom << endl;
      cout << num << " " << denom << endl;
      CheckHist(h_pt_denom);
      CheckHist(h_pt_num);
      TH2D* h_mcpt_num= (TH2D*)fmc->Get(num.Data());
      TH2D* h_mcpt_denom= (TH2D*)fmc->Get(denom.Data());
      CheckHist(h_mcpt_denom);
      CheckHist(h_mcpt_num);
      cout << "tets" << endl;
      
      TString name = *it+ *it2 ;
      
      TH2D* eff_rate = (TH2D*)h_pt_num->Clone(("FakeRate_" + name).Data());
      cout << "tets" << endl;
      TH2D* hratedenom = (TH2D*)h_pt_denom->Clone((name +"_denom").Data());
      eff_rate->Add(h_mcpt_num,-1.);
      hratedenom->Add(h_mcpt_denom, -1.);
      eff_rate->Divide(eff_rate,hratedenom,1.,1.,"cl=0.683 b(1,1) mode");
      eff_rate->Write();
      
      //TCanvas* c1 = new TCanvas(("Plot"), "Plot", 1600, 1200);
      //eff_rate->Draw("colz textE");
      //c1->SaveAs(("/home/jalmond/WebPlots/13TeV/Fakes/2D_ps_hntight_aj40.pdf"));
    }
  }



  std::vector<TString> fakesdxy;
  fakesdxy.push_back("ELECTRON_HN_HIGHDXY_TIGHT_eldxy");

  
  for(vector<TString>::iterator it = fakesdxy.begin(); it!=fakesdxy.end(); ++it){
    if(!CheckFile(fdata))return;
    if(!CheckFile(fmc))return;
    TString denom ="LooseEl" + *it  +"_pt_eta";
    TString num ="TightEl" + *it  +"_pt_eta";
    TH2D* h_pt_num= (TH2D*)fdata->Get(num.Data());
    TH2D* h_pt_denom= (TH2D*)fdata->Get(denom.Data());
    
    cout << h_pt_num << " " << h_pt_denom << endl;
    cout << num << " " << denom << endl;
    CheckHist(h_pt_denom);
    CheckHist(h_pt_num);
    TH2D* h_mcpt_num= (TH2D*)fmc->Get(num.Data());
    TH2D* h_mcpt_denom= (TH2D*)fmc->Get(denom.Data());
    CheckHist(h_mcpt_denom);
    CheckHist(h_mcpt_num);
    cout << "tets" << endl;
    
    TString name = *it ;
    
    TH2D* eff_rate = (TH2D*)h_pt_num->Clone(("FakeRate_" + name).Data());
    TH2D* hratedenom = (TH2D*)h_pt_denom->Clone((name +"_denom").Data());
    eff_rate->Add(h_mcpt_num,-1.);
    hratedenom->Add(h_mcpt_denom, -1.);
    eff_rate->Divide(eff_rate,hratedenom,1.,1.,"cl=0.683 b(1,1) mode");
    eff_rate->Write();
  }
  


  

  for(vector<TString>::iterator it2 = fakes40.begin(); it2!=fakes40.end(); ++it2){
    for(vector<TString>::iterator it3 = isocut.begin(); it3!=isocut.end(); ++it3){
      for(vector<TString>::iterator it = fakes_opt.begin(); it!=fakes_opt.end(); ++it){

	if(it3->Contains("dxy") && it->Contains("mini")) continue;

        if(!CheckFile(fdata))return;
        if(!CheckFile(fmc))return;
        TString denom ="LooseEl" + *it + *it3 +"_dijet_nod0_"+ *it2;
        TString num ="TightEl" + *it +  *it3 +"_dijet_nod0_"+ *it2;
	if(it->Contains("dxy")){
	  denom ="LooseEl" + *it + *it3 +"_dijet_nod0_dxysig_"+ *it2;
	  num ="TightEl" + *it +  *it3 +"_dijet_nod0_dxysig_"+ *it2;
	}
	if(it->Contains("miniiso")){
          denom ="LooseEl" + *it + *it3+ "_dijet_nod0_dxysig_miniiso_"+ *it2;
          num ="TightEl" + *it + *it3+  "_dijet_nod0_dxysig_miniiso_"+ *it2;
        }

        TH2D* h_pt_num= (TH2D*)fdata->Get(num.Data());
        TH2D* h_pt_denom= (TH2D*)fdata->Get(denom.Data());

        cout << h_pt_num << " " << h_pt_denom << endl;
        cout << num << " " << denom << endl;
        CheckHist(h_pt_denom);
        CheckHist(h_pt_num);
        TH2D* h_mcpt_num= (TH2D*)fmc->Get(num.Data());
        TH2D* h_mcpt_denom= (TH2D*)fmc->Get(denom.Data());
        CheckHist(h_mcpt_denom);
        CheckHist(h_mcpt_num);

        TString name = *it +*it2 + *it3;

        TH2D* eff_rate = (TH2D*)h_pt_num->Clone(("FakeRate_" + name).Data());
        TH2D* hratedenom = (TH2D*)h_pt_denom->Clone((name +"_denom").Data());
        eff_rate->Add(h_mcpt_num,-1.);
        hratedenom->Add(h_mcpt_denom, -1.);
        eff_rate->Divide(eff_rate,hratedenom,1.,1.,"cl=0.683 b(1,1) mode");
        eff_rate->Write();
      }
    }
  }


  for(vector<TString>::iterator it2 = fakes40.begin(); it2!=fakes40.end(); ++it2){
    for(vector<TString>::iterator it3 = isocut.begin(); it3!=isocut.end(); ++it3){
      for(vector<TString>::iterator it = fakes_opt.begin(); it!=fakes_opt.end(); ++it){

        if(it3->Contains("dxy") && it->Contains("mini")) continue;

        if(!CheckFile(fdata))return;
        if(!CheckFile(fmc))return;
        TString denom ="LooseEl" + *it + *it3 +"_dijet_d0_"+ *it2;
        TString num ="TightEl" + *it +  *it3 +"_dijet_d0_"+ *it2;
        if(it->Contains("dxy")){
          denom ="LooseEl" + *it + *it3 +"_dijet_d0_dxysig_"+ *it2;
          num ="TightEl" + *it +  *it3 +"_dijet_d0_dxysig_"+ *it2;
        }
        if(it->Contains("miniiso")){
          denom ="LooseEl" + *it + *it3+ "_dijet_d0_dxysig_miniiso_"+ *it2;
          num ="TightEl" + *it + *it3+  "_dijet_d0_dxysig_miniiso_"+ *it2;
        }

        TH2D* h_pt_num= (TH2D*)fdata->Get(num.Data());
        TH2D* h_pt_denom= (TH2D*)fdata->Get(denom.Data());

        cout << h_pt_num << " " << h_pt_denom << endl;
        cout << num << " " << denom << endl;
        CheckHist(h_pt_denom);
        CheckHist(h_pt_num);
        TH2D* h_mcpt_num= (TH2D*)fmc->Get(num.Data());
        TH2D* h_mcpt_denom= (TH2D*)fmc->Get(denom.Data());
        CheckHist(h_mcpt_denom);
        CheckHist(h_mcpt_num);

        TString name = *it +*it2 + *it3+"_d0";

        TH2D* eff_rate = (TH2D*)h_pt_num->Clone(("FakeRate_" + name).Data());
        TH2D* hratedenom = (TH2D*)h_pt_denom->Clone((name +"_denom").Data());
        eff_rate->Add(h_mcpt_num,-1.);
        hratedenom->Add(h_mcpt_denom, -1.);
        eff_rate->Divide(eff_rate,hratedenom,1.,1.,"cl=0.683 b(1,1) mode");
        eff_rate->Write();
      }
    }
  }



  return;

}

  
bool CheckFile(TFile* f ){
    bool file_exist = true;
    if(!f){
      cout << "File " << f->GetName() << " does not exist. Exiting " << endl;
      file_exist = false;
    }
    
    return file_exist;
}

bool CheckHist(TH2* h ){
  bool hist_exist = true;
  if(!h){
    cout << "No histogram with name " << h->GetName() << endl;
    hist_exist= false;
  }
  return hist_exist;
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
