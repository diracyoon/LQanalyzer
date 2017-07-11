#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TTree.h"
#include "TMath.h"
#include "TKey.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <TStyle.h>
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TLatex.h"

//#include "Macro.h"
void setTDRStyle();
void FixOverFlows(TH1* hist, float max_x) ;
void FixUnderFlows(TH1* hist) ;
void FixOverUnderFlows(TH1* hist, float max_x) ;
TH1* MakeSumHist2(THStack* thestack);

void PlotBkgMM(){
 
  
  setTDRStyle();
  gStyle->SetPalette(1);
  
  std::vector<TString> masses;
  masses.push_back("40");
  masses.push_back("100");
  //masses.push_back("200");
  //masses.push_back("500");
  //masses.push_back("1100");


  TCanvas* c1 = new TCanvas("Plot", "Plot", 800, 600);
      
  TLegend *legend = new TLegend(.7, 0.7, .9, 0.9);
  legend->SetFillColor(10);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.04);

  TString path= "/afs/cern.ch/work/j/jalmond/CAT/HNDiMuonOptimisation/periodBtoH/";
  TFile * fnp = new TFile(path+"HNDiMuonOptimisation_DoubleMuon_SKnonprompt_hndilep_cat_v8-0-7.root");
  TFile * fwz = new TFile(path+"HNDiMuonOptimisation_MC_trilep_cat_v8-0-7.root");
  
  gSystem->Exec("rm list.txt");
  gSystem->Exec("python ~/scripts/listkeys.py -f " + path + "HNDiMuonOptimisation_DoubleMuon_SKnonprompt_hndilep_cat_v8-0-7.root >> list.txt");

  ifstream rlist("list.txt");

  if (!rlist){
    exit(EXIT_FAILURE);
  }

  std::vector<TString>  cuts;
  TString cut = "/passtrigcuts_double_pogt";
  
  std::vector<TString> hists;  
  std::vector<TString> hists2D;  

  string rline;
  while(getline(rlist, rline)){
    TString type;
    TString name;
    std::istringstream is(rline);
    is >> type;
    is >> name;
    if(!name.Contains(cut)) continue;
    
    if(type == "TH1D") hists.push_back(name);
    if(type == "TH2D") hists.push_back(name);
  }
  

  cout << "List of plots = " << endl;
  
  std::vector<TString> fillhists;
  std::vector<TString> emptyhists;
  for(int ih=0; ih < hists.size() ; ih++){

    TH1D* h_np= (TH1D*)(fnp->Get(hists[ih]))->Clone(hists[ih]+"tmp");
    TH1D* h_mc= (TH1D*)(fwz->Get(hists[ih]))->Clone(hists[ih]+"tmp");
    if((h_np->Integral() + h_mc->Integral()) > 0) fillhists.push_back(hists[ih]);
    else {
      emptyhists.push_back(hists[ih]);
      continue;
    }
    cout << hists[ih] + "...." << endl;
  }

  cout << " " << endl;
  cout << "Following hitsograms for cut " << cut << "  are empty " << endl;
  
  for(int ih=0; ih < emptyhists.size() ; ih++){
    cout << emptyhists[ih] << " ####### " << endl;
  }

  for(int ih=0; ih < fillhists.size() ; ih++){
    
    TH1D* h_np= (TH1D*)(fnp->Get(fillhists[ih]))->Clone(hists[ih]+"_f");
    TH1D* h_mc= (TH1D*)(fwz->Get(fillhists[ih]))->Clone(hists[ih]+"_f");

    if(!h_mc) return ;
    if(!h_np) return ;

    if(ih==0){
      legend->AddEntry(h_np,"NonPrompt","f");
      legend->AddEntry(h_mc,"Prompt","f");
    }
    
    
    h_np->SetLineColor(870);
    h_np->SetFillColor(870);

    h_mc->SetLineColor(74);
    h_mc->SetFillColor(74);

    string clonename = fillhists[ih].Data();



    /// Fix ranges
    int nbin = h_np->GetNbinsX();
    
    float xmax = h_np->GetBinLowEdge(nbin+1);

    float xmax_init=xmax;
    float xmin_init=h_np->GetBinLowEdge(1);
    FixOverUnderFlows(h_np,xmax);
    FixOverUnderFlows(h_mc,xmax);
    
    int ibinx(1);
    float xmin_con =h_np->GetBinContent(ibinx) + h_mc->GetBinContent(ibinx);
    float xmin  = h_np->GetBinLowEdge(ibinx);
    while(xmin_con <= 0){
      ibinx++;
      xmin_con =h_np->GetBinContent(ibinx) + h_mc->GetBinContent(ibinx);
      xmin  = h_np->GetBinLowEdge(ibinx);
      if(ibinx == nbin) break;
    }
    
    if(xmin_con==0) continue;
    
    h_np->GetXaxis()->SetRangeUser(xmin, xmax);
    h_mc->GetXaxis()->SetRangeUser(xmin, xmax);

    float range_x=xmax_init-xmin_init;
    
    cout << "range = " << xmax_init-xmin << " nbin = " << nbin << endl;
    float bins_size = (xmax_init-xmin_init) / float(nbin);
        
    int ndiv=1;

    float tmp_bins_size= bins_size;
    cout << "Bin size = " << tmp_bins_size <<  " rangex = " <<  range_x<<endl;
    if(range_x > 400){
      while(tmp_bins_size < 5) {
	if (ndiv > 10) break;
	cout << "In while "<<  endl;
	tmp_bins_size= bins_size;
	ndiv++;    
	cout << nbin << " " << ndiv << " " <<  nbin%ndiv << endl;
	while(nbin%ndiv){
	  cout << "change bin size " << ndiv << " " << tmp_bins_size << endl;
	  tmp_bins_size = tmp_bins_size*ndiv;
	  if (tmp_bins_size> 5) break;
	  ndiv++;    
	}
	if (tmp_bins_size> 5)break;
      }
    }
    cout << "pass" << endl;
    
    if(xmin!= 0  || xmax_init!= xmax)cout << "%%Format ::: Changing xaxis range from 0 - " << xmax_init  << " to " << xmin << " _ " << xmax << endl;
    if(ndiv!=1)cout << "%%Format ::: Changing Number of bins for " << fillhists[ih].Data() << " from " << h_np->GetNbinsX() << " to : "  << h_np->GetNbinsX() / ndiv << endl;
    
    h_np->GetXaxis()->SetRangeUser(xmin, xmax);
    h_mc->GetXaxis()->SetRangeUser(xmin, xmax);
    
    if(xmin < 0) xmax = fabs(xmin);

    h_np->Rebin(ndiv);
    h_mc->Rebin(ndiv);
    
    THStack* stack = new THStack(clonename.c_str(), clonename.c_str());

    stack->Add(h_np);
    stack->Add(h_mc);

    TH1* h_nominal = MakeSumHist2(stack);
    h_nominal->SetFillColor(kWhite);
    h_nominal->SetLineColor(kBlack);
    
    TString xtit = TString(clonename).ReplaceAll(cut,"");
    h_nominal->GetXaxis()->SetTitle(xtit);
    h_nominal->Draw("hist");
    h_nominal->GetXaxis()->SetRangeUser(xmin, xmax);
    
    stack->Draw("HIST9same");
    

    for(unsigned int im=0; im < masses.size(); im++){
      TFile * fsig = new TFile((path+"HNDiMuonOptimisation_HNMumMum_"+masses[im]+"_cat_v8-0-7.root").Data());
      
      TH1F* h_cutflow= (TH1F*)fsig->Get(fillhists[ih]);


      if(!h_cutflow) return;
      h_cutflow->Rebin(ndiv);
      h_cutflow->GetXaxis()->SetRangeUser(xmin, xmax);
      
      h_cutflow->Scale(0.5*h_np->Integral() / h_cutflow->Integral());
      h_cutflow->SetLineWidth(4.);
      h_cutflow->Draw("histsame");
      
      h_cutflow->SetLineStyle(2.);
      if(im==0)h_cutflow->SetLineColor(kCyan);
      if(im==1)h_cutflow->SetLineColor(kGreen-2);
      if(im==2)h_cutflow->SetLineColor(kBlue);
      if(im==3)h_cutflow->SetLineColor(kBlue+4);
      if(im==4)h_cutflow->SetLineColor(kViolet);
      if(ih==0){
	legend->AddEntry(h_cutflow,"HN(--) " + masses.at(im),"l");
      }
    }
    
    legend->Draw();
    TString savename="plots/"+fillhists[ih].ReplaceAll("/","_");
    c1->SaveAs(savename+".pdf");
    
  }
  
  
}


TH1* MakeSumHist2(THStack* thestack){

  TH1* hsum=0;
  TList* list = thestack->GetHists();
  TIter it(list, true);
  TObject* obj=0;
  while( (obj = it.Next()) ) {
    TH1* h = dynamic_cast<TH1*>(obj);
    if(!hsum) hsum = (TH1*)h->Clone( (string(h->GetName()) + "_sum").c_str() );
    else {
      hsum->Add(h, 1.0);
    }
  }//hist loop

  return hsum;
}



void FixOverUnderFlows(TH1* hist, float max_x) {
  FixUnderFlows(hist);
  FixOverFlows(hist, max_x);
}


void FixOverFlows(TH1* hist, float max_x) {

  max_x = max_x - 0.01;
  const int nbins = hist->GetNbinsX();
  const int bin_max_x = hist->FindBin(max_x);
  
  if(hist->GetBinContent(nbins+1) > 0) { //of content 
    
    const double of = hist->GetBinContent(nbins+1);
    const double of_err = hist->GetBinError(nbins+1);
    const double last = hist->GetBinContent(nbins);
    const double last_err = hist->GetBinError(nbins);

    const double err = TMath::Sqrt(of_err*of_err + last_err*last_err);

    hist->SetBinContent(nbins+1,0.0);
    hist->SetBinContent(nbins, last + of);
    hist->SetBinError(nbins+1,0.0);
    hist->SetBinError(nbins,err);

  }

  for(int i = (nbins+1); i > (bin_max_x); i--){
    const double of = hist->GetBinContent(i);
    const double of_err = hist->GetBinError(i);
    const double last = hist->GetBinContent(i-1);
    const double last_err = hist->GetBinError(i-1);

    const double err = TMath::Sqrt(of_err*of_err + last_err*last_err);

    hist->SetBinContent(i,0.0);
    hist->SetBinContent(i-1, last + of);
    hist->SetBinError(i,0.0);
    hist->SetBinError(i-1,err);
  }


}//FixOverFlows

void FixUnderFlows(TH1* hist) {

  if(hist->GetBinContent(0) > 0) { //uf content 
    const double uf = hist->GetBinContent(0);
    const double uf_err = hist->GetBinError(0);
    const double first = hist->GetBinContent(1);
    const double first_err = hist->GetBinError(1);

    const double err = TMath::Sqrt(uf_err*uf_err + first_err*first_err);

    hist->SetBinContent(0, 0.0);
    hist->SetBinContent(1, first + uf);
    hist->SetBinError(0, 0.0);
    hist->SetBinError(1, err);

  }

}//FixUnderFlows



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



