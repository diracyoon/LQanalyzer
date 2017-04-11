// LOCAL includes
#include "Reweight.h"

/// c++ std libraries
#include <sstream>
#include <iostream>
#include <stdio.h>

// ROOT include(s): 
#include "TFile.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH1F.h"

using namespace std;

Reweight::Reweight(TString filename){
  
  std::cout << "Initialising Reweight class " << std::endl;

  TDirectory* origDir = gDirectory;
  
  fileData_ = TFile::Open(filename, "READ");
  
  fileMC_   = new TFile("/data1/LQAnalyzer_rootfiles_for_analysis/CATAnalysis2016/Pileup/MCPileUp2016Moriond.root", "READ");
  
  
  TDirectory* tempDir = getTemporaryDirectory();
  tempDir->cd();
  
  h_Data_ = 0;
  h_Data_0j = 0;
  h_Data_1j = 0;
  
  h_Data_ = dynamic_cast<TH1D*>((fileData_->Get("Nvtx_nocut_data"))->Clone());;
  h_Data_0j = dynamic_cast<TH1D*>((fileData_->Get("Nvtx_0j_nocut_data"))->Clone());;
  h_Data_1j = dynamic_cast<TH1D*>((fileData_->Get("Nvtx_1j_nocut_data"))->Clone());;
  


  h_MCmod_ = (TH1D*)fileMC_->Get("Nvtx_nocut_mc");
  h_MCmod_0j = (TH1D*)fileMC_->Get("Nvtx_0j_nocut_mc");
  h_MCmod_1j = (TH1D*)fileMC_->Get("Nvtx_1j_nocut_mc");

  double int_MC_ = h_MCmod_->Integral();
  double int_Data_ = h_Data_->Integral();

  h_Data_->Divide(h_MCmod_);
  h_Data_->Scale(int_MC_ / int_Data_);

  double int_MC_0j = h_MCmod_0j->Integral();
  double int_Data_0j = h_Data_0j->Integral();

  h_Data_0j->Divide(h_MCmod_0j);
  h_Data_0j->Scale(int_MC_0j / int_Data_0j);

  double int_MC_1j = h_MCmod_1j->Integral();
  double int_Data_1j = h_Data_1j->Integral();

  h_Data_1j->Divide(h_MCmod_1j);
  h_Data_1j->Scale(int_MC_1j / int_Data_1j);

  
  fileData_->Close();
  delete fileData_;
  fileMC_->Close();
  delete fileMC_;


  origDir->cd();
  


}

Reweight::~Reweight(){
  delete h_MCmod_;
  delete h_Data_;

}


double Reweight::GetWeight(Int_t nvtx, TString version, int njet){

  if (njet < 0)   return h_Data_->GetBinContent( h_Data_->FindBin(nvtx)  );
  else if(njet ==0)   return h_Data_0j->GetBinContent( h_Data_0j->FindBin(nvtx)  );
  else return h_Data_1j->GetBinContent( h_Data_1j->FindBin(nvtx)  );

  return h_Data_->GetBinContent( h_Data_->FindBin(nvtx)  );
  
}


TDirectory* Reweight::getTemporaryDirectory(void) const
{

  gROOT->cd();
  TDirectory* tempDir = 0;
  int counter = 0;
  while (not tempDir) {
    // First, let's find a directory name that doesn't exist yet:                                                                                                                 
    std::stringstream dirname;
    dirname << "Reweight_%i" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    // Let's try to make this directory:                                                                                                                                          
    tempDir = gROOT->mkdir((dirname.str()).c_str());
    
  }
  
  return tempDir;
}
