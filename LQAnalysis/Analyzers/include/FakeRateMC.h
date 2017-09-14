#ifndef FakeRateMC_h
#define FakeRateMC_h

#include "AnalyzerCore.h"


class FakeRateMC : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  FakeRateMC();
  ~FakeRateMC();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  

  void ExecuteEventsMuon(TString idloose, TString idtight, TString tag, float iso, double w);
  void ExecuteEventsElectron(TString idloose, TString idtight, TString tag, float iso, double w);

  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();

  TDirectory* getTemporaryDirectory(void) const;


  float GetPrescaleMu( std::vector<snu::KMuon> muons,bool pass3, bool pass2, bool pass1, float fake_total_lum );
  float GetPrescaleEl( std::vector<snu::KElectron> electrons,bool pass5,  bool pass4, bool pass3, bool pass2, bool pass1, float fake_total_lum );

  void GetFakeRates(std::vector<snu::KMuon> loose_el, std::vector<snu::KMuon> tight_el,TString tightlabel,  std::vector<snu::KJet> jets,std::vector<snu::KJet> alljets, TString tag, double w, float isocut,bool makebasicplots);
  

  int CloseJetType(snu::KElectron el, std::vector<snu::KJet> jets);
  int CloseJetType(snu::KMuon mu, std::vector<snu::KJet> jets);
  int CloseJetTypeHAD(snu::KMuon mu, std::vector<snu::KJet> jets);
  int CloseJetTypeHAD(snu::KElectron el, std::vector<snu::KJet> jets);
  int AwayJetType(snu::KElectron el, std::vector<snu::KJet> jets);


  void MakeMCFakes(std::vector<snu::KElectron> fake_electrons, TString tag, std::vector<snu::KJet> jets,  std::vector<snu::KJet> alljets, TString tightid, float w);
  

  void MakeMCFakeratePlots(TString label, bool pass_single_trigger, std::vector<snu::KMuon> muons, std::vector<snu::KJet> jets, std::vector<snu::KJet> alljets, double w);
  void MakePlotsMCAwaJetPt(TString label,float awayjetptcut, std::vector<snu::KMuon>, std::vector<snu::KJet> jets , std::vector<snu::KJet> alljets, double w);

  
  void MakeMCPlots(TString label, snu::KMuon muon, std::vector<snu::KJet> jets, std::vector<snu::KJet> alljets, double w);
  float GetPrescale( std::vector<snu::KMuon> muons,bool passlowest, bool passlow,bool passhigh, float fake_total_lum);
  bool  UseEvent(std::vector<snu::KMuon> muons,  vector<snu::KJet> jets, float pcut,float precale_weight, float w);
  bool  UseEventAllMT(std::vector<snu::KMuon> muons,  vector<snu::KJet> jets, float pcut,float precale_weight, float w);



  void MakeFakeRatePlots(TString label, TString eltag,   std::vector<snu::KMuon> tightmuons,std::vector<snu::KMuon> muons,   std::vector <snu::KJet> jets, std::vector<snu::KJet> alljets,  float precale_w, float w, float isocut, bool makebasicplots);



  void MakeSingleMuonCRPlots(TString looseid, TString eltag, TString tightid, float w, bool usepujetid);

  void GetFakeRateAndPromptRates(std::vector<snu::KMuon> muons, TString eltag, std::vector<snu::KMuon> tightmuons, float w, float isocut,bool usepujetid, bool runall);

  void GetFakeRateAndPromptRatesPerPeriod(TString looseid, TString eltag, TString tightid, float w, bool usepujetid, bool runall);
private:
  
  
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_electrons;
  std::vector<snu::KMuon> out_muons;

  int n_17_pass;
  int n_17_17_jet_pass;
  int n_17_jet_pass;

  TH2D* MuonFR;
  TH2D* MuonFRcorr;
  TH2D*   MuonFRcbj;
  TH2D* MuonFRcbjcorr;
  TH2D*   MuonFRncbj;
  TH2D* MuonFRncbjcorr;

  TH2D* ElFR;
  TH2D* ElFRcorr;
  TH2D* ElFRcbj;
  TH2D* ElFRcbjcorr;
  TH2D* ElFRncbj;
  TH2D* ElFRncbjcorr;

  
  ClassDef ( FakeRateMC, 1);
};
#endif
