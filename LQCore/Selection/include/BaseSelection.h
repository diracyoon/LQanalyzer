#ifndef BaseSelection_h
#define BaseSelection_h

#include <iostream>
using namespace std;

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "AnalysisBase.h"
#include "LQEvent.h"

class BaseSelection {

 public:


  enum ID {
    ELECTRON_POG_VETO              = 0,
    ELECTRON_POG_LOOSE             = 1,
    ELECTRON_POG_MEDIUM            = 2,
    ELECTRON_POG_TIGHT             = 3,
    ELECTRON_POG_MVATrig           = 4,
    ELECTRON_POG_MVANonTrig        = 5,
    ELECTRON_POG_MVA_WP80          = 6,
    ELECTRON_POG_MVA_WP90          = 7,
    ELECTRON_ECAL_FIDUCIAL         = 8,
    ELECTRON_HN_VETO               = 9,
    ELECTRON_HN_TIGHT              = 10,
    ELECTRON_HN_FAKELOOSE          = 11,
    ELECTRON_HN_FAKELOOSE_NOD0     = 12,
    ELECTRON_HN_MVA_LOOSE          = 13,
    ELECTRON_HN_MVA_TIGHT          = 14,
    ELECTRON_TOP_VETO              = 15,
    ELECTRON_TOP_LOOSE             = 16,
    ELECTRON_TOP_TIGHT             = 17,
    ELECTRON_PTETA                 = 18,
    ELECTRON_NOCUT                 = 19,
    MUON_POG_LOOSE                 = 20,
    MUON_POG_MEDIUM                = 21,
    MUON_POG_TIGHT                 = 22,
    MUON_HN_VETO                   = 23,
    MUON_HN_FAKELOOSE              = 24,
    MUON_HN_TIGHT                  = 25,
    MUON_FAKELOOSE                 = 26,
    MUON_TOP_VETO                  = 27,
    MUON_TOP_LOOSE                 = 28,
    MUON_TOP_TIGHT                 = 29,
    MUON_PTETA                     = 30,
    MUON_NOCUT                     = 31,
    PFJET_LOOSE                    = 32,
    PFJET_MEDIUM                   = 33,
    PFJET_TIGHT                    = 34,
    JET_HN                         = 35,
    JET_HN_TChannel                = 36,
    JET_NOLEPTONVETO               = 37,
    JET_LOOSE                      = 38,
    JET_TIGHT                      = 39,
    PHOTON_POG_LOOSE               = 40, 
    PHOTON_POG_MEDIUM              = 41, 
    PHOTON_POG_TIGHT               = 42,
  };

  Int_t ifid;
  
  //// Selection cuts
  Bool_t etaPt,RelIsod0Chi2,DepositVeto,individual,RelIsod0;//fiducial;//muonid,pTcut,isIso;

  //// variables to make cuts
  Double_t LeptonRelIso,dxy,dz,D0,D0Error,D0Significance,Vxy, Vz;
  Double_t LeptonchiNdof;
  Int_t numlep;
  Bool_t jetIsOK;

  /// variables to set cuts
  Double_t pt_cut_min, pt_cut_max, eta_cut_min, eta_cut, 
           jpt_cut_min, jpt_cut_max, jeta_cut_min, jeta_cut,
           relIso_cut, relIsoMIN_cut, chiNdof_cut, chiNdofMIN_cut, 
           dxy_cut, dxyMIN_cut, dz_cut,
           dxySig_min, dxySig_max;
  Double_t relIsoBarrel_max, relIsoBarrel_min, dxyBarrel_max, dzBarrel_max,
           relIsoEndcap_max, relIsoEndcap_min, dxyEndcap_max, dzEndcap_max;

  TString  RelIsoType, HLTSafeLevel, PUJetIDWP;

  Int_t casediscriminator,simpleselection;

  /// LQEvent class object. Stores all vectors
  LQEvent k_lqevent;


  TString AccessStringMap(TString label, TString id);
  float AccessFloatMap(TString label, TString id);
  bool CheckCutFloat(TString label, TString id);
  bool CheckCutString(TString label, TString id);
  TString GetCutString(TString label, TString id);
  TString GetString(ID id);

  vector<pair<TString,TString> > GetStringList(TString id);
  vector<pair<TString,float> > GetFloatList(TString id);

  map<TString,vector<pair<TString,TString> > > k_stringmap;
  map<TString,vector<pair<TString,float> > > k_floatmap;

  void SetIDSMap(std::map<TString, vector<pair<TString,TString> > > smap);
  void SetIDFMap(std::map<TString, vector<pair<TString,float> > > smap);

  bool Check(float val);
  void PrintSIDMap();
  void PrintFIDMap();


  BaseSelection();
  BaseSelection& operator= (const BaseSelection& obj);
  BaseSelection(const BaseSelection& bs);

  ~BaseSelection();

  /// bools to tell selector to apply cuts
  Bool_t apply_ptcut,apply_etacut, apply_jptcut,apply_jetacut, apply_relisocut, apply_chi2cut, apply_dxycut, apply_dzcut, apply_general, apply_deposit;
  Bool_t apply_ID, apply_convcut, apply_chargeconst, apply_HLTSafeCut, applypileuptool;
  Bool_t apply_dxysigmin, apply_dxysigmax;
  Bool_t apply_BESepCut, apply_BETrRegIncl;

  
  ID k_id;
  
  /// SetCut functions
  void reset();
  void SetPt(Double_t minPt, Double_t maxPt);
  void SetPt(Double_t minPt);
  void SetJetPt(Double_t minPt);
  void SetEta(Double_t Eta);
  void SetUseJetPileUp(bool use);
  void SetPileUpJetID(bool use, TString pujetidwp="Loose");
  void SetJetEta(Double_t Eta);
  void SetEta(Double_t minEta, Double_t Eta);
  void SetRelIso(Double_t RelIso);
  void SetRelIso(Double_t RelIsoMIN, Double_t RelIso);
  void SetChiNdof(Double_t ChiNdof);
  void SetChiNdof(Double_t ChiNdofMIN, Double_t ChiNdof);
  void SetBSdxy(Double_t dxy);
  void SetBSdxy(Double_t dxyMIN, Double_t dxy);
  void SetBSdz(Double_t dz);
  void SetID(ID  id);
  void SetHLTSafeCut(TString level);
  void SetCheckCharge(bool check);
  void SetApplyConvVeto(bool apply);

  void SetBETrRegIncl(bool include=false);
  void SetRelIsoBEMax(Double_t BarrelIso, Double_t EndcapIso);
  void SetRelIsoBEMin(Double_t BarrelIso, Double_t EndcapIso);
  void SetdxyBEMax(Double_t Barreldxy, Double_t Endcapdxy);
  void SetdzBEMax(Double_t Barreldz, Double_t Endcapdz);
  void SetRelIsoType(TString type);
  void SetdxySigMin(Double_t dxySigMin);
  void SetdxySigMax(Double_t dxySigMax);

};

#endif
