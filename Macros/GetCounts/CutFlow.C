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


void CutFlow(TString path){
  TFile * fnp = new TFile(path);
  
  TIter next(fnp->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    if(!key->IsFolder()){
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1D")) continue;
      TH1 *h = (TH1*)key->ReadObj();
      if (TString(key->GetName()).Contains("cutflow")){
	if(!h->Integral()) continue;
	string tmpname  = key->GetName();
	TString ts_tmpname = TString(tmpname);
	cout << "Cutflow name = " << ts_tmpname << endl;
	for(unsigned int i=1; i < h->GetNbinsX()+1; i++){
	  cout << "Bin label = " << h->GetXaxis()->GetBinLabel(i) << " integral= " <<  h->GetBinContent(i) << endl;
	}
      }
    }
  }
  
}
