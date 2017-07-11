#ifndef FakeRateCalculator_Mu_h
#define FakeRateCalculator_Mu_h

#include "AnalyzerCore.h"


class FakeRateCalculator_Mu : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  FakeRateCalculator_Mu();
  ~FakeRateCalculator_Mu();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();


  void RunFakes(TString tag, TString ID);

  void GetFakeRates(std::vector<snu::KMuon> loose_el, std::vector<snu::KMuon> tight_el,TString tightlabel,  std::vector<snu::KJet> jets,std::vector<snu::KJet> alljets, TString tag, double w, float isocut,bool makebasicplots);
  
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
  
  ClassDef ( FakeRateCalculator_Mu, 1);
};
#endif
