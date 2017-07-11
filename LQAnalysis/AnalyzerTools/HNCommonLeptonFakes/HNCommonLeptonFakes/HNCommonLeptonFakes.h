#ifndef HNCommonLeptonFakes_H
#define HNCommonLeptonFakes_H

// STL include(s):                                                                                                                                                                   
#include <vector>
#include <map>
#include <string>

// Forward declaration(s):                                                                                                                                                           

#include "TString.h"
#include "TString.h"
#include "TH1.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TLorentzVector.h"

class HNCommonLeptonFakes {

 public:
  /// Destructor
  ~HNCommonLeptonFakes();

  /// constructor
  HNCommonLeptonFakes(std::string path = "../share/", bool usegev=false);

//Declaration of the function

  /// Function to return event weight with error

  /** get temporary directory                                                                                                                                                        
   *                                                                                                                                                                                
   */
  TDirectory*
    getTemporaryDirectory(void) const;



  float get_dilepton_mm_eventweight(TString fakerates, bool geterr, std::vector<TLorentzVector> muons, bool isT1, bool isT2, TString ID);
  float get_dilepton_mm_eventweight(TString fakerates, bool geterr, std::vector<TLorentzVector> muons, bool isT1, bool isT2, TString ID1,TString ID2, bool c1, bool c2, TString method, float ptcorr1, float ptcorr2,bool checkcj);
							 

  float get_dilepton_em_eventweight(bool geterr,std::vector<TLorentzVector> muons, std::vector<TLorentzVector> electron, bool ismu1tight, bool isel1tight);

  float get_dilepton_ee_eventweight(bool geterr, std::vector<TLorentzVector> electrons, bool isel1tight, bool isel2tight, TString cut);
  float get_dilepton_ee_eventweight(bool geterr, std::vector<TLorentzVector> electrons, bool isel1tight, bool isel2tight, std::vector<TString> cut);
  float get_dilepton_ee_eventweight(bool geterr, std::vector<TLorentzVector> electrons, bool isel1tight, bool isel2tight, TString cut1, TString cut2);
  

  float getFakeRate_muon(int sys, float pt,  float eta, TString ID);
  float getEfficiency_muon(int sys,float pt, float eta);

  float getFakeRate_electron(int sys, float pt,  float eta);
  float getEfficiency_electron(int sys,float pt, float eta);
  float getPromptRate_electron(int sys,float pt, float eta, TString cut);

  float getFakeRate_mc_muon(TString tag, float pt, float eta);

  float getFakeRate_electron(int sys, float pt,  float ht, TString cut);
  float getFakeRate_electronEta(int sys, float pt,  float eta, TString cut);
  float getEfficiency_electron(int sys,float pt, float eta, TString cut);
 


  void APPLY_scalefactor(bool applysf) ;


  /// Function to calculate the event weight from the muon efficiencies
  float CalculateDiLepMMWeight(float r1, float f1, float r2, float f2, bool t1, bool t2 );
  float CalculateDiLepMMWeight(float r1, float f1, float r2, float f2, bool t1, bool t2 , int eventtype);
  
  //Incase we want single muon event weight
  double lepton_weight(bool loose, double r, double f);
  

  void CheckFile(TFile* file);
  /// Function to reset all variable in the class
  void ResetValues();
  void NullTotals();

  /// Function called to close root file used to read histograms
  void closefile();

  /// Masp for storing histograms in the rootfile
  std::map<TString,TH1*> _EfficiencyMap;
  std::map<TString,TH2F*> _2DEfficiencyMap;
  std::map<TString,TH2D*> _2DEfficiencyMap_Double;
  std::map<TString,TEfficiency*> _TEEfficiencyMap;
  std::map<TString,TH1*> _EfficiencyMapHST;
  std::map<TString,TH2F*> _2DEfficiencyMapHST;


  /// Function to fill maps , null variable and set up conditions for running algorithm
  void InitialiseFake();

  /// Function to exit code if tryng to access a  non existing histogram
  void NoHist(TString hist, TString function = "");

  /// Initialises histograms filled by FillHists
  void SetHists();
  void SetHists(TString hist);
 
  /// Fills histgrams with TLorentz vectors passed to tool and event weights calculated
  void prepFillHists(float w, float mu1pt, float mu2pt, bool tight1, bool tight2);
  void FillHists(TString cut);

  // Function to apply a systematic shift in fakes rates
  void ApplySF(int type);

  // Function to return histograms filled in tool
  TH1* GetHist(bool weighted, TString cut);
 
  /// Adds event by event info to countes
  void AddToTotals(float w,std::pair<float,float> err, bool mu1tight, bool mu2tight);

  //=============
  //==== Trilep
  //=============

  //==== dXYSigMin and LooseRelIsoMax to string
  TString DoubleToTString(double this_dXYSig, double this_RelIso);
  void SetTrilepWP(double this_dXYSig, double this_RelIso);
  void SetUseQCDFake(bool useit);
  void SetDataPeriod(TString period);
  void SetNJet(int nj);
  void SetNBJet(int nbj);
  //==== get PR/FR
  float getTrilepFakeRate_muon(bool geterr, float pt,  float eta);
  float getTrilepPromptRate_muon(bool geterr, float pt, float eta);
  //==== get weight
  float get_dilepton_mm_eventweight(bool geterr, std::vector<TLorentzVector> muons, bool isT1, bool isT2);
  float get_trilepton_mmm_eventweight(bool geterr, std::vector<TLorentzVector> muons, bool isT1, bool isT2, bool isT3);
  float get_eventweight(bool geterr, std::vector<TLorentzVector> muons, TString muid, std::vector<TLorentzVector> electrons, std::vector<TString> elid, std::vector<bool> isT, int HalfSampleErrorDir=0);

  //==== Large dXYSig working poins
  std::vector<double> GetdXYMins();
  std::vector<double> GetRelIsoMaxs();
  //==== After runing get_eventweight, we have # of Loose but not Tight
  int GetNLooseNotTight();

  void SetUsePtCone(bool b);

 private:
  /// vector for storing FakeCR strings
  std::vector<TString> regions;
  std::vector<TString> ch_regions;
  
  /// path to root file
  std::string path_to_root_files;  
  
  /// class string for fake CR used as the central value
  TString defCR;

  /// bool to set output messsage level
  bool m_debug;

  /// bool to correctly set pt/ht to GeV needed to acces correct bin in root file
  bool useGeV;
  
  int fsys;
  int rsys;

  bool apply_scalefactor;

  bool dosystematic;

  /// Event counter
  int nevent;

  /// Map for calculating systematics
  std::map<TString,float> SystMap;
  
  ///Map for plotting cuts
  std::map<TString, std::pair<TH1*,TH1*> > CutHists;

  /// variables used to fill hists
  float _percentageFsys;
  float _percentageRsys;

  //==== Trilep fake
  double Current_dXYSig, Current_RelIso;
  TString DataPeriod;
  bool UseQCDFake;
  int n_jet, n_bjet;
  std::vector<double> dXYMins, RelIsoMaxs;
  int n_Loose_not_Tight;
  bool UsePtCone;

};
#endif
