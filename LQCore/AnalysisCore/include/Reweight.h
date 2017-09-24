#ifndef Reweight_h
#define Reweight_h

class TFile;
class TH1D;
class TH1F;
class TString;
class TDirectory;

class Reweight {

 public:
  
  //// constructors
  Reweight(TString filename);
  ~Reweight();
  
  double GetWeight(int nvtx, TString version, int njet = -999);
  TDirectory* getTemporaryDirectory(void) const;

 private:

  TFile* fileData_;
  TFile* fileMC_;

  TH1D* h_MCmod_;
  TH1D* h_Data_;
  TH1D* h_MCmod_0j;
  TH1D* h_Data_0j;
  TH1D* h_MCmod_1j;
  TH1D* h_Data_1j;
  TH1F* h_MCmodF_;
  TH1F* h_DataF_;


};

#endif

