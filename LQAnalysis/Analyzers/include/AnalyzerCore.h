#ifndef AnalyzerCore_H
#define AnalyzerCore_H

//forward declarations                                                                                                                                            
class EventBase;
class MuonPlots;
class ElectronPlots;
class JetPlots;
class SignalPlotsEE;
class SignalPlotsMM;
class SignalPlotsEM;
class TriLeptonPlots;
class HNpairPlotsMM;
class HNTriLeptonPlots;
class EventBase;

#include "BaseSelection.h"
#include "LQCycleBase.h"
#include "BTag/BTagSFUtil.h"
#include "TNtupleD.h"
#include "TNtuple.h"
#include "TGraphAsymmErrors.h"
#include "MCDataCorrections.h"
#include "DataDrivenBackgrounds.h"
#include "HNGenMatching.h"

#define CSV_THRESHOLD_LOOSE 0.460
#define CSV_THRESHOLD_MEDIUM 0.800
#define CSV_THRESHOLD_TIGHT 0.935

#define DISTANCE_MATCH 0.2

#define BLANK -999

class AnalyzerCore : public LQCycleBase {
  
 public:
 
  enum period  {C=0,
		D=1,
		CtoD=2,
		BtoH_2016=3};
  
  //Default constructor
   
  AnalyzerCore();

  //destructor
  virtual ~AnalyzerCore();

  // SetUpEvent CORE function: accesses event in ntuple
  virtual void SetUpEvent(Long64_t entry, float ev_weight) throw( LQError );
  virtual void EndEvent()throw( LQError );
  virtual void WriteHistograms()throw( LQError );


  void setTDRStyle();
  TString GetStringID(BaseSelection::ID id);
  double GetIsoCut(snu::KElectron el, TString curlable);
  double GetDXYCut(snu::KElectron el, TString curlable);

  std::vector<snu::KJet>  GetJets(BaseSelection::ID jetid, float ptcut=-999., float etacut = -999.);
  std::vector<snu::KJet> GetJetsWFT(TString jetid,TString fatjetid,float ptcut=-999., float etacut = -999.);
  std::vector<snu::KFatJet>  GetFatJets(BaseSelection::ID jetid, float ptcut=-999., float etacut = -999.);
  std::vector<snu::KMuon> GetMuons(BaseSelection::ID muid, float ptcut=-999., float etacut = -999.);
  std::vector<snu::KElectron> GetElectrons( BaseSelection::ID elid , float ptcut=-999., float etacut = -999.);

  float GetKFactor();

  bool ISCF(snu::KElectron el);
  bool TruthMatched(std::vector<snu::KElectron> el, bool tightdxy, bool allowCF);
  bool TruthMatched(snu::KElectron el, bool keepcf);


  snu::KJet GetCorrectedJetCloseToLepton(snu::KElectron lep,snu::KJet jet, bool usem=true);
  snu::KJet GetCorrectedJetCloseToLepton(snu::KMuon lep, snu::KJet jet);
  snu::KTruth GetTruthMatchedParticle(snu::KElectron el);
  
  
  float GetVirtualMassConv(int cmindex,int nconvindx);
  float GetVirtualMass(bool inph=false);

  float MassDrop(snu::KElectron electron, std::vector<snu::KJet> jets, bool usecorrectedpt=true);
  float MassDrop(snu::KMuon muon, std::vector<snu::KJet> jets, bool usecorrectedpt=true);

  float GetJetsCloseToLeptonPt(snu::KElectron electron, std::vector<snu::KJet> jets, bool usecorrectedpt=true);
  float GetJetsCloseToLeptonPt(snu::KMuon muon, std::vector<snu::KJet> jets, bool usecorrectedpt=true);
  float GetPtRelLepTJet(snu::KElectron electron, std::vector<snu::KJet> jets, bool usecorrectedpt=true);
  float GetPtRelLepTJet(snu::KMuon muon, std::vector<snu::KJet> jets, bool usecorrectedpt=true);


  void SetupLuminosityMap(bool initialsetup, TString forceperiod="");
  Int_t GetMCPeriod();
  Int_t GetDataPeriod();  
  int GetPeriod();
  int GetMCPeriodRandom();

  void  SetupID();
  void  SetupDDBkg();

  bool IsDiEl();

  bool  Check(float val);

  std::vector<snu::KMuon> GetMuons(BaseSelection::ID muid,bool keepfakes, float ptcut=-999., float etacut = -999.);
  std::vector<snu::KElectron> GetElectrons(bool keepcf, bool keepfake, BaseSelection::ID elid , float ptcut=-999., float etacut = -999.);

  std::vector<snu::KJet>  GetJets(TString jetid,  float ptcut=-999., float etacut = -999.);
  std::vector<snu::KFatJet>  GetFatJets(TString jetid,  float ptcut=-999., float etacut = -999.);
  std::vector<snu::KMuon> GetMuons(TString muid, float ptcut=-999., float etacut = -999.);
  std::vector<snu::KMuon> GetMuons(TString muid, bool keepfakes, float ptcut=-999., float etacut = -999.);
  std::vector<snu::KElectron> GetElectrons(bool keepcf, bool keepfake,   TString elid, float ptcut=-999., float etacut = -999.);
  std::vector<snu::KElectron> GetElectrons( TString elid , float ptcut=-999., float etacut = -999.);

  bool Is2015Analysis();
  void SetupSelectionMuon(std::string path_sel);
  void SetupSelectionJet(std::string path_sel);
  void SetupSelectionFatJet(std::string path_sel);
  void SetupSelectionElectron(std::string path_sel);

  

  void MakeBTagEfficiencyPlots();
  void GetJetTaggerEfficiences(TString taggerWP, snu::KJet::Tagger tag,  snu::KJet::WORKING_POINT wp);

  vector<TString >  GetHNDiLepElTriggers();

  bool CheckEventComparison(TString user, TString label);
  void FillEventComparisonFile(TString label);
  void FillCutFlow(TString cut, float weight);
  bool TriggerMatch(TString trigname, vector<snu::KMuon> mu);

  bool EtaRegion(TString reg,  std::vector<snu::KMuon> muons);
  bool EtaRegion(TString reg,  std::vector<snu::KElectron> electrons);

  void FillHistPerLumi(TString histname, float value, float w, float xmin, float xmax,int nbins, int nlumibins);

  bool HasCloseBJet(snu::KElectron el, snu::KJet::Tagger tag=snu::KJet::CSVv2 , snu::KJet::WORKING_POINT wp= snu::KJet::Medium, int mcperiod=-1);
  int NBJet(std::vector<snu::KJet> jets,  snu::KJet::Tagger tag=snu::KJet::CSVv2, snu::KJet::WORKING_POINT wp = snu::KJet::Medium, int mcperiod=-1);

  bool IsBTagged(snu::KJet jet,  snu::KJet::Tagger tag, snu::KJet::WORKING_POINT wp, int mcperiod=-1);
  float BTagScaleFactor_1a(std::vector<snu::KJet> jetColl, snu::KJet::Tagger tag, snu::KJet::WORKING_POINT wp, int mcperiod=-1);

  int AssignnNumberOfTruth();
  bool IsSignal();

  void ClassInfo();

  int VersionStamp(TString cversion);


  bool IsCF(snu::KElectron el);

  float GetDiLepMass(std::vector<snu::KMuon> muons);
  float GetDiLepMass(std::vector<snu::KElectron> electrons);


  float  JetResCorr(snu::KJet jet, std::vector<snu::KGenJet> genjets);
  float SumPt( std::vector<snu::KJet> particles);
  float SumPt( std::vector<snu::KFatJet> particles);
  bool isPrompt(long pdgid);
  void TruthPrintOut();

  bool FailHNDataSetCheck();
  bool PassJets(std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets );
  bool PassID(snu::KElectron electron, TString id);
  bool PassID(snu::KMuon mu, TString muid);

  bool IsTight(snu::KMuon muon);
  std::vector<snu::KElectron> GetTruePrompt(vector<snu::KElectron> electrons,  bool keep_chargeflip, bool keepfake);
  std::vector<snu::KMuon> GetTruePrompt(vector<snu::KMuon> muons,   bool keepfake);

  bool  MCIsCF(snu::KElectron el);
  bool Zcandidate(vector<snu::KMuon> muons, float interval, bool require_os=true);
  bool Zcandidate(vector<snu::KElectron> electrons, float interval, bool require_os=true);
  bool SameCharge(std::vector<snu::KMuon> muons);
  bool SameCharge(std::vector<snu::KElectron> electrons, bool runcf=false);
  bool OppositeCharge(std::vector<snu::KElectron> electrons, bool runcf=false);

  bool OppositeCharge(std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons);  
  float CorrectedMETRochester(std::vector<snu::KMuon> muons , bool updatemet);
  float CorrectedMETElectron(std::vector<snu::KElectron> electrons,  int syst=0);
  float CorrectedMETMuon(std::vector<snu::KMuon> muons ,int syst=0);


  void CorrectMuonMomentum(vector<snu::KMuon>& k_muons);
  void SetCorrectedMomentum(vector<snu::KMuon>& k_muons);

  std::vector<TLorentzVector> MakeTLorentz( std::vector<snu::KElectron> el);
  std::vector<TLorentzVector> MakeTLorentz( std::vector<snu::KMuon> mu);
  std::vector<TLorentzVector> MakeTLorentz( std::vector<snu::KJet> jet);
  std::vector<TLorentzVector> MakeTLorentz( std::vector<snu::KFatJet> jet);



  double MuonDYMassCorrection(std::vector<snu::KMuon> mu, double w);

  // enum for plotting functions/classes
  enum histtype {muhist, elhist, jethist, sighist_e,sighist_ee,sighist_eee,sighist_eeee, sighist_m,sighist_mm,sighist_mmm,sighist_mmmm, sighist_em, trilephist, hnpairmm, hntrilephist};
  
  
  //
  // Useful message function 
  //
  void Message(TString message, LQMsgType type=INFO);


  //
  //  Specify which triggers will be avaiable in KTrigger
  //
  void AddTriggerToList(TString triggername);
  
  /// Pileup Reweighting class
  static const Bool_t MC_pu = true;

  bool k_debugmode;
  MCDataCorrections* mcdata_correction;
  DataDrivenBackgrounds* m_datadriven_bkg;

  //// Event base pointer. Used to get all objects for analysis
  EventBase* eventbase;
  
  UInt_t numberVertices;
  Bool_t *goodVerticiesB;

  bool k_usegenjet;
  bool k_usetruth;
  bool k_usephotons;
  bool k_usefatjet;

  TDirectory *Dir;
  map<TString, TH1*> maphist;
  map<TString, TH2*> maphist2D;
  map<TString, TH3*> maphist3D;
  map<TString, TNtupleD*> mapntp;
  map<int,int> compmap;

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

  /// Event weights
  Double_t MCweight, weight;


  bool reset_lumi_mask;
  bool changed_target_lumi;

  // used to get trigger prescale
  
  bool k_onlyfromtaus;
  bool  k_reset_period;
  int a_mcperiod;
  bool IDSetup;
  bool setupDDBkg;

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
  map<TString, HNTriLeptonPlots*> mapCLhistHNTriLep;
  
  float WeightByTrigger(TString triggername, float tlumi);
  float WeightByTrigger(vector<TString> triggername, float tlumi);  

  //
  // Function that closes rootfile
  //
  void CloseFiles();
  

  //
  // Make Histograms and fill maphist
  //
  void MakeHistograms();
  void MakeHistograms(TString hname, int nbins, float xmin, float xmax, TString label="");
  void MakeHistograms(TString hname, int nbins, float xbins[], TString label="");
  void MakeHistograms2D(TString hname, int nbinsx, float xbins[], int nbinsy, float ybins[], TString label="");
  void MakeHistograms2D(TString hname, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax, TString label="");
  void MakeHistograms3D(TString hname, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax, int nbinsz, float zmin, float zmax,TString label="");
    //
    // Makes temporary dir
    //
  TDirectory* GetTemporaryDirectory(void) const;                                                                                                                                 
  //
  // Checks if a file exists
  //
  
  //// Plotting 
  TH1* GetHist(TString hname);
  TH2* GetHist2D(TString hname);
  TH3* GetHist3D(TString hname);

  /// Changed  Default json file

  /// Fills hist in maphist
  void FillHist(TString histname, float value, float w , TString label="");
  void FillHist(TString histname, float value, float w , float xmin, float xmax, int nbins=0 , TString label="");
  void FillHist(TString histname, float value, float w , float xmin[], int nbins=0 , TString label="");
  void FillHist(TString histname, float value1, float value2, float w , float x[], int nbinsx, float y[], int nbinsy , TString label="");
  void FillHist(TString histname, float value1,  float value2, float w , float xmin, float xmax, int nbinsx,  float ymin, float ymax, int nbinsy , TString label="");
  void FillHist(TString histname, float value1,  float value2, float value3, float w , float xmin, float xmax, int nbinsx,  float ymin, float ymax, int nbinsy , float zmin, float zmax, int nbinz,TString label="");

  /// Fills clever hists
  void FillCLHist(histtype type, TString hist, snu::KEvent ev,vector<snu::KMuon> muons, vector<snu::KElectron> electrons, vector<snu::KJet> jets,double weight);
  void FillCLHist(histtype type, TString hist, snu::KEvent ev,vector<snu::KMuon> muons, vector<snu::KElectron> electrons, vector<snu::KJet> jets, vector<snu::KFatJet> fjets,double weight);
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
  bool  PassTrigger(TString list);
  bool PassTriggerOR(vector<TString> list);
  bool  PassTrigger(std::vector<std::pair<TString,TString> > list);

  void ListTriggersAvailable();
  bool PassMETFilter();

  std::map<TString,BTagSFUtil*> MapBTagSF;

  //==== (Trilepton) HeavyN stuffs

  HNGenMatching *m_HNgenmatch;

  void PutNuPz(TLorentzVector *nu, double Pz);
  void PutNuPz(snu::KParticle *nu, double Pz);
  double solveqdeq(double W_mass, TLorentzVector l1l2l3, double MET, double METphi, TString pm);
  int find_mlmet_closest_to_W(snu::KParticle  lep[], snu::KParticle  MET, int n_lep=3);
  double MT(TLorentzVector a, TLorentzVector b);
  bool GenMatching(snu::KParticle a, snu::KParticle b, double maxDeltaR, double maxPtDiff);
  std::vector<snu::KMuon> GetHNTriMuonsByLooseRelIso(double LooseRelIsoMax, bool keepfake);
  void PrintTruth();
  std::vector<snu::KMuon> sort_muons_ptorder(std::vector<snu::KMuon> muons);

  //private methods
  static Bool_t Compare_Jet_Pt(const snu::KJet& jet0, const snu::KJet& jet1);
  static Bool_t Compare_Jet_CSV(const snu::KJet& jet0, const snu::KJet& jet1);

  Double_t Distance(const snu::KTruth& truth, const snu::KJet& jet);

  Bool_t Search_Truth_Value(std::vector<snu::KTruth>& gen_truth_coll, snu::KTruth gen_quark[4], snu::KTruth& gen_neutrino, snu::KTruth& gen_lepton);
  enum TRUTH_TABLE {TOP, A_TOP, BOTTOM, D_0, D_0_0, D_0_1, A_BOTTOM, D_1, D_1_0, D_1_1} Truth_Table;

  Bool_t Truth_Jet_Match(const snu::KTruth gen_quark[], const vector<snu::KJet>& jet_vector, Int_t permutation_truth[4]);
};
#endif
