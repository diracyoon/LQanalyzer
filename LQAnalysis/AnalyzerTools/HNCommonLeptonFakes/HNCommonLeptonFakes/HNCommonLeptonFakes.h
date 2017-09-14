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
  HNCommonLeptonFakes(TString path, std::map<TString, std::pair<std::pair<TString,TString> ,std::pair<float,TString> > > fake_hist_config);
  HNCommonLeptonFakes();

//Declaration of the function

  /// Function to return event weight with error

  /** get temporary directory                                                                                                                                                        
   *                                                                                                                                                                                
   */
  TDirectory*
    getTemporaryDirectory(void) const;


  void IsData(bool isdata);
  void SetupFake(TString path, std::map<TString, std::pair<std::pair<TString,TString>  ,std::pair<float,TString> > >fake_hist_config);
  void SetupFake();
  void SetMinPtMuon(float minpt);
  void SetMinPtElectron(float minpt);
  bool IsInitialised();

  float get_dilepton_mm_eventweight(bool geterr, std::vector<TLorentzVector> muons, bool isT1, bool isT2, TString key, float value1, float value2, int vbkg, TString prkey="");
  float get_dilepton_em_eventweight(bool geterr,std::vector<TLorentzVector> muons, std::vector<TLorentzVector> electron, bool isT1, bool isT2, TString key1,TString key2, float value1, float value2, int vbkg);
  float get_dilepton_ee_eventweight(bool geterr, std::vector<TLorentzVector> electrons, bool isT1, bool isT2, TString key, float value1, float value2, int vbkg, TString prkey="");

  float getFakeRate(float pt,  float eta, TString tag, float ptmin, bool geterr=false);
  float getPromptRate(float pt, float eta,TString tag);
  
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
  std::map<TString,TH2D*> _2DEfficiencyMap_Double;
  std::map<TString,TEfficiency*> _TEEfficiencyMap;

  std::map<TString,float> _2DPtMaxMap_Double;


  /// Function to fill maps , null variable and set up conditions for running algorithm
  void InitialiseFake(TString path, std::map<TString, std::pair<std::pair<TString,TString> ,std::pair<float,TString> > > fake_hist_config);

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


 private:

  /// path to root file
  std::string path_to_root_files;  
  
  /// bool to set output messsage level
  bool m_debug;


  bool _isdata;
  bool _isinitialised;


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

  TString DataPeriod;


  float _minpt_muon;
  float _minpt_electron;

};
#endif
