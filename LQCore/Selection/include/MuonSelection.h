#ifndef MuonSelection_h
#define MuonSelection_h

#include "BaseSelection.h"
#include "KMuon.h"

class MuonSelection : public BaseSelection {
  Double_t ECalDeposit_max, HCalDeposit_max, ECalDeposit_min, HCalDeposit_min;

  static const Double_t AreaTrackerMu[5];
  static const Double_t AreaEcalMu[5];
  static const Double_t AreaHcalMu[5];

 public:

  Bool_t MuonID;

  MuonSelection(LQEvent ev);
  ~MuonSelection();
  
  MuonSelection& operator= (const MuonSelection& obj);
  MuonSelection(const MuonSelection& bs);

  void SelectMuons(std::vector<snu::KMuon>& leptonColl, ID muid,   float ptcut=-999., float etacut=-999.);
  void SelectMuons(std::vector<snu::KMuon>& leptonColl, TString muid, float ptcut=-999., float etacut=-999.);
  void SelectMuons(std::vector<snu::KMuon>& leptonColl, TString muid,vector<pair<TString, TString> > vids, vector<pair<TString, float> > vidf, float ptcut=-999., float etacut=-999.);
  void SelectMuons(std::vector<snu::KMuon>& leptonColl, ID muid,vector<pair<TString, TString> > vids, vector<pair<TString, float> > vidf, float ptcut=-999., float etacut=-999.);

  bool MuonPass(snu::KMuon muon, TString muid,  float ptcut=-999., float etacut=-999.);
  bool MuonPass(snu::KMuon muon, TString muid, vector<pair<TString, TString> > vids, vector<pair<TString, float> > vidf,  float ptcut=-999., float etacut=-999.);

  float IsoCutValue(snu::KMuon muon,TString muid);
  bool PassUserID(TString id, snu::KMuon mu);
  bool PassUserID(TString id, snu::KMuon mu,vector<pair<TString, TString> > vids, vector<pair<TString, float> > vidf);
  bool PassID(TString id, snu::KMuon mu, bool checkdxy=true, bool checkdz=true, bool checkchi2=true, bool m_debug = false);

  
  //// General Selection
  void Selection(std::vector<snu::KMuon>& leptonColl, bool applyrochester=true, TString Option="");

  /// selection for treemaker
  void BasicSelection(std::vector<snu::KMuon>& leptonColl, bool m_debug = false);

  void SkimSelection( std::vector<snu::KMuon>& leptonColl, bool m_debug = false); 

  //// HN analysis selection
  bool HNIsTight(snu::KMuon muon,  bool m_debug = false);
  bool POGID(snu::KMuon muon, TString ID);

  bool HNTightMuonSelection(snu::KMuon mu);
  bool HNLooseMuonSelection(snu::KMuon mu);
  bool HNVetoMuonSelection(snu::KMuon mu);

  //// Top analysis selection
  bool TopTightMuonSelection(snu::KMuon mu);
  bool TopLooseMuonSelection(snu::KMuon mu);
  bool TopVetoMuonSelection(snu::KMuon mu);


  
  //// General 
  void SetDeposits(Double_t ECalDeposit, Double_t HCalDeposit);
  void SetDeposits(Double_t ECalDeposit1 , Double_t HCalDeposit1, Double_t ECalDeposit2 , Double_t HCalDeposit2);
  
};

#endif
