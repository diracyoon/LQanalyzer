#include "TH1.h"
#include "THStack.h"
#include "TMath.h"
#include "HistUtils.hpp"
#include "TList.h"
#include "TCollection.h"
#include "TObject.h"
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <map>
#include <sstream>


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



double IntegralStack(THStack* stack) {

  double val=0;
  TList* list = stack->GetHists();
  TIter it(list, true);
  TObject* obj=0;
  while( (obj = it.Next()) ) {
    TH1* h = dynamic_cast<TH1*>(obj);
    val += h->Integral();
  }
  return val;
}
