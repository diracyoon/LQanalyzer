#ifndef __AnalyzerCore_H__
#define __AnalyzerCore_H__

//forward declarations                                                                                
class Reweight;
class EventBase;
class MuonPlots;
class ElectronPlots;
class JetPlots;
class SignalPlotsEE;
class SignalPlotsMM;
class SignalPlotsEM;
class TriLeptonPlots;
class HNpairPlotsMM;
class EventBase;

#include "BaseSelection.h"
#include "LQCycleBase.h"
#include "HNCommonLeptonFakes/HNCommonLeptonFakes/HNCommonLeptonFakes.h"
#include "rochcor2015/rochcor2015.h"
#include "rochcor2015/RoccoR.h"
#include "BTag/BTagSFUtil.h"
#include "TNtupleD.h"
#include "TNtuple.h"
#include "TGraphAsymmErrors.h"

class AnalyzerCore : public LQCycleBase 
{
 public:
 
  enum period  {C=0, D=1, CtoD=2, BtoE_2016=3};
  
  AnalyzerCore();
  virtual ~AnalyzerCore();
  
  //SetUpEvent CORE function: accesses event in ntuple
  virtual void SetUpEvent(Long64_t entry, float ev_weight) throw(LQError);
  virtual void EndEvent() throw(LQError);
  virtual void WriteHistograms() throw(LQError);

  TDirectory*   getTemporaryDirectory(void) const;
 
  TString GetStringID(BaseSelection::ID id);
  std::vector<snu::KJet>  GetJets(BaseSelection::ID jetid, float ptcut=-999., float etacut = -999.);
  std::vector<snu::KFatJet>  GetFatJets(BaseSelection::ID jetid, float ptcut=-999., float etacut = -999.);
  std::vector<snu::KMuon> GetMuons(BaseSelection::ID muid, float ptcut=-999., float etacut = -999.);
  std::vector<snu::KElectron> GetElectrons( BaseSelection::ID elid , float ptcut=-999., float etacut = -999.);


  bool IsDiEl();

  std::vector<snu::KMuon> GetMuons(BaseSelection::ID muid,bool keepfakes, float ptcut=-999., float etacut = -999.);
  std::vector<snu::KElectron> GetElectrons(bool keepcf, bool keepfake, BaseSelection::ID elid , float ptcut=-999., float etacut = -999.);

  std::vector<snu::KJet>  GetJets(TString jetid, bool smearjets=true, float ptcut=-999., float etacut = -999.);
  std::vector<snu::KFatJet>  GetFatJets(TString jetid, bool smearjets=true, float ptcut=-999., float etacut = -999.);
  std::vector<snu::KMuon> GetMuons(TString muid, float ptcut=-999., float etacut = -999.);
  std::vector<snu::KMuon> GetMuons(TString muid, bool keepfakes, float ptcut=-999., float etacut = -999.);
  std::vector<snu::KElectron> GetElectrons(bool keepcf, bool keepfake, TString elid, float ptcut=-999., float etacut = -999.);
  std::vector<snu::KElectron> GetElectrons( TString elid , float ptcut=-999., float etacut = -999.);

  bool Is2015Analysis();
  void SetupSelectionMuon(std::string path_sel);
  void SetupSelectionJet(std::string path_sel);
  void SetupSelectionFatJet(std::string path_sel);
  void SetupSelectionElectron(std::string path_sel);

  void FillCorrectionHist(string label, string dirname, string filename, string histsname, string histtype);
  void FillCorrectionHists();
  TH2F* GetCorrectionHist(TString label);
  bool CheckCorrectionHist(TString label);

  TGraphAsymmErrors* GetCorrectionGraph(TString label);
  bool CheckCorrectionGraph(TString label);

  void FillCutFlow(TString cut, float weight);

  bool TriggerMatch(TString trigname, vector<snu::KMuon> mu);

  bool EtaRegion(TString reg,  std::vector<snu::KMuon> muons);
  bool EtaRegion(TString reg,  std::vector<snu::KElectron> electrons);

  void FillHistPerLumi(TString histname, float value, float w, float xmin, float xmax,int nbins, int nlumibins);

  bool HasCloseBJet(snu::KElectron el, snu::KJet::Tagger tag=snu::KJet::CSVv2 , snu::KJet::WORKING_POINT wp= snu::KJet::Medium);
  int NBJet(std::vector<snu::KJet> jets,  snu::KJet::Tagger tag=snu::KJet::CSVv2, snu::KJet::WORKING_POINT wp = snu::KJet::Medium);

  int IsBTagged(snu::KJet jet,  snu::KJet::Tagger tag, snu::KJet::WORKING_POINT wp);

  int AssignnNumberOfTruth();
  bool IsSignal();

  void ClassInfo();

  int VersionStamp(TString cversion);


  float WeightCFEvent(std::vector<snu::KElectron> electrons, bool runchargeflip);
  bool IsCF(snu::KElectron el);


  double TriggerScaleFactor( vector<snu::KElectron> el, vector<snu::KMuon> mu, TString trigname, int direction=0);;

  float GetDiLepMass(std::vector<snu::KMuon> muons);
  float GetDiLepMass(std::vector<snu::KElectron> electrons);

  double ElectronScaleFactor( TString  elid, vector<snu::KElectron> el, int sys=0);
  double ElectronRecoScaleFactor(vector<snu::KElectron> el);

  double MuonScaleFactor(TString  muid, vector<snu::KMuon> mu, int sys=0);
  double MuonISOScaleFactor(TString muid, vector<snu::KMuon> mu,int sys=0);
  double MuonTrackingEffScaleFactor(vector<snu::KMuon> mu);

  float  JetResCorr(snu::KJet jet, std::vector<snu::KGenJet> genjets);
  float SumPt( std::vector<snu::KJet> particles);
  float SumPt( std::vector<snu::KFatJet> particles);
  bool isPrompt(long pdgid);
  void TruthPrintOut();
  bool IsTight(snu::KElectron electron);
  bool IsTight(snu::KMuon muon);
  std::vector<snu::KElectron> GetTruePrompt(vector<snu::KElectron> electrons,  bool keep_chargeflip, bool keepfake);
  std::vector<snu::KMuon> GetTruePrompt(vector<snu::KMuon> muons,   bool keepfake);

  bool Zcandidate(vector<snu::KMuon> muons, float interval, bool require_os=true);
  bool Zcandidate(vector<snu::KElectron> electrons, float interval, bool require_os=true);
  bool SameCharge(std::vector<snu::KMuon> muons);
  bool SameCharge(std::vector<snu::KElectron> electrons, bool runcf=false);
  
  float CFRate(snu::KElectron el);
  float CFRate_Run2(snu::KElectron el, TString el_id);
  std::vector<snu::KElectron>  ShiftElectronEnergy(std::vector<snu::KElectron> el, bool applyshift);

  float Get_DataDrivenWeight_EE(vector<snu::KElectron> k_electrons);
  float Get_DataDrivenWeight_MM(vector<snu::KMuon> k_muons, TString cutID="HN");
  float Get_DataDrivenWeight_E(vector<snu::KElectron> k_electrons);
  float Get_DataDrivenWeight_M(vector<snu::KMuon> k_muons, TString cutID="HN");
  float Get_DataDrivenWeight_EM(vector<snu::KMuon> k_muons, vector<snu::KElectron> k_electrons, TString cut="HN");

  void CorrectMuonMomentum(vector<snu::KMuon>& k_muons);

  double MuonDYMassCorrection(std::vector<snu::KMuon> mu, double w);
  
  vector<TLorentzVector> MakeTLorentz( vector<snu::KElectron> el);
  vector<TLorentzVector> MakeTLorentz( vector<snu::KMuon> mu);
  vector<TLorentzVector> MakeTLorentz( vector<snu::KJet> jet);
  vector<TLorentzVector> MakeTLorentz( vector<snu::KFatJet> jet);
  // enum for plotting functions/classes
  enum histtype {muhist, elhist, jethist, sighist_ee, sighist_mm, sighist_em, trilephist, hnpairmm};
  
  //Useful message function 
  void Message(TString message, LQMsgType type=INFO);

  //Specify which triggers will be avaiable in KTrigger
  void AddTriggerToList(TString triggername);
  
  //Pileup Reweighting class
  static const Bool_t MC_pu = true;
  Reweight *reweightPU;

  //Event base pointer. Used to get all objects for analysis
  EventBase* eventbase;
  
  UInt_t numberVertices;
  Bool_t *goodVerticiesB;

  TDirectory *Dir;
  map<TString, TH1*> maphist;
  map<TString, TH2*> maphist2D;
  map<TString, TNtupleD*> mapntp;

  map<TString, TH2F*>  CorrectionMap;
  map<TString, TGraphAsymmErrors*>  CorrectionMapGraph;

  map<int, float> mapLumi; 
  map<int, float> mapBadLumi; 
  map<int, float> mapLumiPerBlock;
  map<int, TString> mapLumiNamePerBlock;
  map<TString,float> trigger_lumi_map_cat2015;

  map<TString,vector<pair<TString,TString> > > selectionIDMapsMuon;
  map<TString,vector<pair<TString,float> > > selectionIDMapfMuon;
  map<TString,vector<pair<TString,TString> > > selectionIDMapsElectron;
  map<TString,vector<pair<TString,float> > > selectionIDMapfElectron;
  map<TString,vector<pair<TString,TString> > > selectionIDMapsJet;
  map<TString,vector<pair<TString,float> > > selectionIDMapfJet;
  map<TString,vector<pair<TString,TString> > > selectionIDMapsFatJet;
  map<TString,vector<pair<TString,float> > > selectionIDMapfFatJet;

  map<int, float> mapLumi2016;
  map<int, float> mapBadLumi2016;
  map<int, float> mapLumiPerBlock2016;
  map<int, TString> mapLumiNamePerBlock2016;
  map<TString,float> trigger_lumi_map_cat2016;

  std::vector<TString> cutflow_list;

  TH2F* FRHist;
  TH2F* MuonSF;
  TH2F* SingleMuon_274093;
  TH2F* SingleMuon_276097;
  TH2F* ElectronSF_Tight;
  TH2F* ElectronSF_GSF;
  TH2F* ElectronSF_Medium;
  TH2F* ElectronSF_Loose;
  TH2F* ElectronSF_Veto;
  TH2F* ElectronRECO;
  TH2F* MuonID_tight;
  TH2F* MuonID_medium;
  TH2F* MuonID_loose;
  TH2F* MuonISO_tight_tightID;
  TH2F* MuonISO_tight_mediumID;
  TH2F* MuonISO_loose_tightID;
  TH2F* MuonISO_loose_mediumID;
  TH2F* MuonISO_loose_looseID;
  HNCommonLeptonFakes* m_fakeobj;

  int n_cutflowcuts;

  //Event weights
  Double_t MCweight, weight;

  bool reset_lumi_mask;
  bool changed_target_lumi;

  //used to get trigger prescale
  Int_t prescale;
  
  std::vector<TString> triggerlist;

  //// Making cleaver hist maps
  map<TString, TriLeptonPlots*> mapCLhistTriLep;
  map<TString, HNpairPlotsMM*> mapCLhistHNpairMM;
  map<TString, SignalPlotsEE*> mapCLhistSigEE;
  map<TString, SignalPlotsMM*> mapCLhistSigMM;
  map<TString, SignalPlotsEM*> mapCLhistSigEM;
  map<TString, ElectronPlots*> mapCLhistEl;
  map<TString, MuonPlots*> mapCLhistMu;
  map<TString, JetPlots*> mapCLhistJet;
  
  float WeightByTrigger(TString triggername, float tlumi);
  float WeightByTrigger(vector<TString> triggername, float tlumi);  

  // Function that closes rootfile
  void CloseFiles();
  
  //Make Histograms and fill maphist
  void MakeHistograms();
  void MakeHistograms(TString hname, int nbins, float xmin, float xmax, TString label="");
  void MakeHistograms(TString hname, int nbins, float xbins[], TString label="");
  void MakeHistograms2D(TString hname, int nbinsx, float xbins[], int nbinsy, float ybins[], TString label="");
  void MakeHistograms2D(TString hname, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax, TString label="");
    
  //Makes temporary dir
  TDirectory* GetTemporaryDirectory(void) const;                                                                                                                                 
  
  //Checks if a file exists
  void CheckFile(TFile* file) throw( LQError );
  
  //Plotting 
  TH1* GetHist(TString hname);
  TH2* GetHist2D(TString hname);

  /// Changed  Default json file

  /// Fills hist in maphist
  void FillHist(TString histname, float value, float w , TString label="");
  void FillHist(TString histname, float value, float w , float xmin, float xmax, int nbins=0 , TString label="");
  void FillHist(TString histname, float value, float w , float xmin[], int nbins=0 , TString label="");
  void FillHist(TString histname, float value1, float value2, float w , float x[], int nbinsx, float y[], int nbinsy , TString label="");
  void FillHist(TString histname, float value1,  float value2, float w , float xmin, float xmax, int nbinsx,  float ymin, float ymax, int nbinsy , TString label="");

  /// Fills clever hists
  void FillCLHist(histtype type, TString hist, snu::KEvent ev,vector<snu::KMuon> muons, vector<snu::KElectron> electrons, vector<snu::KJet> jets,double weight);
  void FillCLHist(histtype type, TString hist, snu::KEvent ev,vector<snu::KMuon> muons, vector<snu::KElectron> electrons, vector<snu::KJet> jets,double weight,int nbjet);

  void FillCLHist(histtype type, TString hist, vector<snu::KMuon> muons , double weight);
  void FillCLHist(histtype type, TString hist, vector<snu::KElectron> electrons , double weight);
  void FillCLHist(histtype type, TString hist, vector<snu::KJet> jets , double weight);
  std::map<TString,BTagSFUtil*> SetupBTagger(std::vector<TString> taggers, std::vector<TString> wps);

  // Makes clever histograms
  void MakeCleverHistograms(histtype type, TString clhistname );

  /// File related                                                                                                                                                
  void OpenPutputFile();
  void WriteHists();
  void WriteCLHists();

  void MakeNtp(TString hname, TString myvar);
  TNtupleD* GetNtp(TString hname);
  void FillNtp(TString hname, Double_t myinput[]);  
  void WriteNtp();


  //// Event related                                                                                                                                              
  float TempPileupWeight();

  bool  PassTrigger(TString list);
  bool  PassTrigger(std::vector<std::pair<TString,TString> > list);
  float TriggerEff(TString trigname, std::vector<snu::KElectron> electrons);
  float TriggerEff(TString trigname, std::vector<snu::KMuon> muons);
  float TriggerEff(TString trigname,  std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons);

  float GetEff(snu::KMuon mu, TString trigname);
  float GetEff(snu::KElectron el, TString trigname);


  void ListTriggersAvailable();
  bool PassMETFilter();

  std::map<TString,BTagSFUtil*> MapBTagSF;

  //BTagSFUtil *lBTagSF, *hBTagSF;
  rochcor2015 *rmcor;

  //private methods                                                        
  static bool Compare_Jet_Pt(const snu::KJet& jet0, const snu::KJet& jet1);  
};

#endif
