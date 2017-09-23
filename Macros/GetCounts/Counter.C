#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TROOT.h"
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


void Counter(TString path){
  TFile * fnp = new TFile(path);
  
  TIter next(fnp->GetListOfKeys());
  TKey *key;
  
  while ((key = (TKey*)next())) {
    
    if(key->IsFolder()){
      TDirectory * d = (TDirectory*)(fnp->Get(key->GetName()));
      TIter nextDIR(d->GetListOfKeys());
      TKey *keyDIR;
      
      while ((keyDIR = (TKey*)nextDIR())) {
	TClass *cl = gROOT->GetClass(keyDIR->GetClassName());
	if (!cl->InheritsFrom("TH1D")) continue;
	TH1 *h = (TH1*)keyDIR->ReadObj();
	if (TString(keyDIR->GetName()).Contains("h_Nelectrons")){
	  if(!h->Integral()) continue;
	  string tmpname  = keyDIR->GetName();
	  TString ts_tmpname = TString(tmpname);
	  ts_tmpname=ts_tmpname.ReplaceAll("h_Nelectrons","");
	  double err(0.);
	  double integral = h->IntegralAndError(0, h->GetNbinsX()+1,err);
	  cout << ts_tmpname << "         =  " << integral << " " << err  << endl;;
	  //if(h->Integral()> 0)cout <<  tmpname << " : " << h->Integral() << endl;
	}
      }
    }
  }
  

}
