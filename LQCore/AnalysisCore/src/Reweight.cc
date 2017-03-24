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
  if (!fileData_){
    cout << "\n\nAt least one of the Nvtx reweighting files could not be opened!\n\n";
    fileData_ = TFile::Open( "/data1/LQAnalyzer_rootfiles_for_analysis/CATAnalysis2016/Pileup/DataPileUp_BtoG_2016.root", "READ");

  }
  if(filename.Contains("BtoE")){
    fileMC_   = new TFile("/data1/LQAnalyzer_rootfiles_for_analysis/CATAnalysis2016/Pileup/MCPileUp2016.root", "READ");
  }
  if(filename.Contains("BtoG")){
    fileMC_   = new TFile("/data1/LQAnalyzer_rootfiles_for_analysis/CATAnalysis2016/Pileup/MCPileUp2016_802.root", "READ");
  }
  if(filename.Contains("Moriond")){
    fileMC_   = new TFile("/data1/LQAnalyzer_rootfiles_for_analysis/CATAnalysis2016/Pileup/MCPileUp2016Moriond.root", "READ");
  }


  TDirectory* tempDir = getTemporaryDirectory();
  tempDir->cd();
  //cout << " filename = " << filename << endl;
  h_Data_ = 0;
  if(filename.Contains("BtoE")){
    h_DataF_ = dynamic_cast<TH1F*>((fileData_->Get("Nvtx_nocut_data"))->Clone());;  
  }
  if(filename.Contains("BtoG")){
    h_Data_ = dynamic_cast<TH1D*>((fileData_->Get("Nvtx_nocut_data"))->Clone());;
  }

  if(filename.Contains("Moriond")){
    h_Data_ = dynamic_cast<TH1D*>((fileData_->Get("Nvtx_nocut_data"))->Clone());;
  }


  //
  if(filename.Contains("BtoE")){
    h_MCmodF_ = (TH1F*)fileMC_->Get("h_VertexNoReweight");
  }
  else{
    h_MCmod_ = (TH1D*)fileMC_->Get("Nvtx_nocut_mc");
  }
  double int_MC_ = h_MCmod_->Integral();
  double int_Data_ = h_Data_->Integral();

  h_Data_->Divide(h_MCmod_);
  h_Data_->Scale(int_MC_ / int_Data_);

  
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


double Reweight::GetWeight(Int_t nvtx, TString version){

  if(version.Contains("v8-0-1")){
    return h_DataF_->GetBinContent( h_DataF_->FindBin(nvtx)  );
  }
  if(version.Contains("v8-0-2")){
    return h_Data_->GetBinContent( h_Data_->FindBin(nvtx)  );
  }  
  return 1.;
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
