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
  
  fileData_ = TFile::Open( filename, "READ");
  if (!fileData_) cout << "\n\nAt least one of the Nvtx reweighting files could not be opened!\n\n";

  TString filenamedown =  "/data1/LQAnalyzer_rootfiles_for_analysis/CatOct15/SNUCAT_Pileup.root";
  TString filenameup = "/data1/LQAnalyzer_rootfiles_for_analysis/CatOct15/SNUCAT_Pileup.root";
  //fileDataUP_ = TFile::Open(filenameup , "READ");
  //if (!fileData_) cout << "\n\nAt least one of the Nvtx reweighting files could not be opened!\n\n";
  //fileDataDOWN_ = TFile::Open( filenamedown, "READ");
  //if (!fileData_) cout << "\n\nAt least one of the Nvtx reweighting files could not be opened!\n\n";

  TDirectory* tempDir = getTemporaryDirectory();
  tempDir->cd();
  
  h_Data_ = 0;
  h_Data_ = dynamic_cast<TH1F*>((fileData_->Get("reweight_nvtx_periodC"))->Clone());;  
  h_Data_up_ = dynamic_cast<TH1F*>((fileData_->Get("reweight_nvtx_periodD"))->Clone());;  
  //h_Data_down_ = dynamic_cast<TH1F*>((fileDataDOWN_->Get("reweight_nvtx_down"))->Clone());;  
  
  
  //h_MCmod_ = dynamic_cast<TH1D*>((h_Data_->Clone("h_MCmod_")));
  
  // Now we can close the file:
  fileData_->Close();
  delete fileData_;
  //fileDataUP_->Close();
  //delete fileDataUP_;
  //fileDataDOWN_->Close();
  //delete fileDataDOWN_;
  // Return to the directory we were in before the function call: 
  origDir->cd();
  
  //  Double_t Summer2012_S10[60] = {

  //};

  //  for (Int_t i = 1; i < 61; i++) {
    //    h_MCmod_->SetBinContent(i, Summer2012_S10[i-1]);
  //}
  
  //  double int_MC_ = h_MCmod_->Integral();
  //double int_Data_ = h_Data_->Integral();
  //double int_Dataup_ = h_Data_up_->Integral();
  //double int_Datadown_ = h_Data_down_->Integral();
  
  //h_Data_->Divide(h_MCmod_);
  //h_Data_->Scale(int_MC_ / int_Data_);
  
  //h_Data_up_->Divide(h_MCmod_);
  //h_Data_up_->Scale(int_MC_ / int_Dataup_);
  
  //h_Data_down_->Divide(h_MCmod_);
  //h_Data_down_->Scale(int_MC_ / int_Datadown_);

}

Reweight::~Reweight(){
}

double Reweight::GetWeight(int nvtx, int period){

  if(period == 1) return h_Data_->GetBinContent( nvtx  ); 
  else return h_Data_up_->GetBinContent( nvtx  );
  
  return h_Data_->GetBinContent( nvtx  );
  
}

double Reweight::GetWeight(int nvtx){
  return GetWeight(nvtx,0);
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
