#include "Macro.h"
#include "TGraphAsymmErrors.h"
#include <vector>

#include "CMS_lumi.h"

void mva(){

  TString path = "/data2/CAT_SKTreeOutput/JobOutPut/jalmond/LQanalyzer//data/output/CAT/HNSignalEfficiencies/periodBtoH/";
  std::vector<TString> masses;
  masses.push_back("40");
  masses.push_back("100");
  masses.push_back("500");
  masses.push_back("1500");


  setTDRStyle();

  TLegend* legendH = new TLegend(0.2, 0.7, 0.4, 0.9);
  std::vector<TH1F*> hists;
  for(int i=0; i < masses.size(); i++){
    TFile * file = new TFile(path + masses[i] + "_cat_v7-6-3.root");

    TH1F* h_den= (TH1F*)file->Get(("TriggerEfficiency_denominator_nojet"));
    TH1F* h_num= (TH1F*)file->Get(("TriggerEfficiency_denominator"));
    
    h_num->Divide(h_den);
    h_num->SetMarkerSize(2.);
    hists.push_back(h_num);
    legendH->AddEntry(h_num, "m_{N}="+masses[i]+ " GeV" , "p");
  }
  

  legendH->SetFillColor(kWhite);
  legendH->SetTextSize(0.03);
  
  hists.at(0)->SetLineColor(kRed);
  hists.at(0)->SetMarkerColor(kRed);
  hists.at(0)->SetMarkerStyle(20.);
  
  hists.at(1)->SetLineColor(kBlue);
  hists.at(1)->SetMarkerColor(kBlue);
  hists.at(1)->SetMarkerStyle(21.);

  hists.at(2)->SetLineColor(kGreen);
  hists.at(2)->SetMarkerColor(kGreen);
  hists.at(2)->SetMarkerStyle(22.);

  hists.at(3)->SetLineColor(kCyan);
  hists.at(3)->SetMarkerColor(kCyan);
  hists.at(3)->SetMarkerStyle(22.);
  TCanvas* c1 = new TCanvas(("Plot"), "Plot", 1600, 1200);
  
  hists.at(0)->GetYaxis()->SetRangeUser(0., 1.5);
  hists.at(0)->Draw("lp");
  hists.at(1)->Draw("lpsame");
  hists.at(2)->Draw("lpsame");
  hists.at(3)->Draw("lpsame");

  
  legendH->Draw();



  
  CMS_lumi( c1, 2, 11 );
  c1->Update();
  c1->RedrawAxis();


  
  c1->SaveAs(("/home/jalmond/WebPlots/13TeV/SignalEff/JetEff.pdf" ));
}


void CMS_lumi( TPad* pad, int iPeriod, int iPosX )
{
  bool outOfFrame    = false;
  if( iPosX/10==0 )
    {
      outOfFrame = true;
    }
  int alignY_=3;
  int alignX_=2;
  if( iPosX/10==0 ) alignX_=1;
  if( iPosX==0    ) alignY_=1;
  if( iPosX/10==1 ) alignX_=1;
  if( iPosX/10==2 ) alignX_=2;
  if( iPosX/10==3 ) alignX_=3;
  int align_ = 10*alignX_ + alignY_;

  float H = pad->GetWh();
  float W = pad->GetWw();
  float l = pad->GetLeftMargin();
  float t = pad->GetTopMargin();
  float r = pad->GetRightMargin();
  float b = pad->GetBottomMargin();
  float e = 0.025;

  pad->cd();

  TString lumiText;
  if( iPeriod==1 )
    {
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
    }
  else if ( iPeriod==2 )
    {
      lumiText += " (8 TeV)";
    }
  else if( iPeriod==3 )
    {
      lumiText = lumi_8TeV;
      lumiText += " (8 TeV)";
      lumiText += " + ";
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
    }
  else if ( iPeriod==4 )
    {
      lumiText += lumi_13TeV;
      lumiText += " (13 TeV)";
    }
  else if ( iPeriod==7 )
    {
      if( outOfFrame ) lumiText += "#scale[0.85]{";
      lumiText += lumi_13TeV;
      lumiText += " (13 TeV)";
      lumiText += " + ";
      lumiText += lumi_8TeV;
      lumiText += " (8 TeV)";
      lumiText += " + ";
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
      if( outOfFrame) lumiText += "}";
    }
  else if ( iPeriod==12 )
    {
      lumiText += "8 TeV";
    }

  cout << lumiText << endl;

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);

  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  latex.SetTextFont(42);
  latex.SetTextAlign(31);
  latex.SetTextSize(lumiTextSize*t);
  latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  if( outOfFrame )
    {
      latex.SetTextFont(cmsTextFont);
      latex.SetTextAlign(11);
      latex.SetTextSize(cmsTextSize*t);
      latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);
    }

  pad->cd();

  float posX_;
  if( iPosX%10<=1 )
    {
      posX_ =   l + relPosX*(1-l-r);
    }
  else if( iPosX%10==2 )
    {
      posX_ =  l + 0.5*(1-l-r);
    }
  else if( iPosX%10==3 )
    {
      posX_ =  1-r - relPosX*(1-l-r);
    }
  float posY_ = 1-t - relPosY*(1-t-b);
  if( !outOfFrame )
    {
      if( drawLogo )
        {
          posX_ =   l + 0.045*(1-l-r)*W/H;
          posY_ = 1-t - 0.045*(1-t-b);
          float xl_0 = posX_;
          float yl_0 = posY_ - 0.15;
          float xl_1 = posX_ + 0.15*H/W;
          float yl_1 = posY_;
          TASImage* CMS_logo = new TASImage("CMS-BW-label.png");
          TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
          pad_logo->Draw();
          pad_logo->cd();
          CMS_logo->Draw("X");
          pad_logo->Modified();
          pad->cd();
        }
      else
        {
          latex.SetTextFont(cmsTextFont);
          latex.SetTextSize(cmsTextSize*t);
          latex.SetTextAlign(align_);
          latex.DrawLatex(posX_, posY_, cmsText);
          if( writeExtraText )
            {
              latex.SetTextFont(extraTextFont);
              latex.SetTextAlign(align_);
              latex.SetTextSize(extraTextSize*t);
              latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText);
            }
        }
    }
  else if( writeExtraText )
    {
      if( iPosX==0)
        {
          posX_ =   l +  relPosX*(1-l-r);
          posY_ =   1-t+lumiTextOffset*t;
        }
      latex.SetTextFont(extraTextFont);
      latex.SetTextSize(extraTextSize*t);
      latex.SetTextAlign(align_);
      latex.DrawLatex(posX_, posY_, extraText);
    }
  return;
}


