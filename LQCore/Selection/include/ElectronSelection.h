#ifndef ElectronSelection_h
#define ElectronSelection_h

#include "KElectron.h"
#include "LQEvent.h"
#include "BaseSelection.h" 

class ElectronSelection : public BaseSelection {

 public:
  
  Bool_t ElectronID;


  ElectronSelection(LQEvent ev);
  ~ElectronSelection();

  ElectronSelection& operator= (const ElectronSelection& obj);
  ElectronSelection(const ElectronSelection& bs);
  
  void Selection(std::vector<snu::KElectron>& leptonColl, bool m_debug= false);		 
  void BasicSelection(std::vector<snu::KElectron>& leptonColl, bool m_debug= false);
  void SkimSelection(std::vector<snu::KElectron>& leptonColl, bool m_debug= false);

  float IsoCutValue(snu::KElectron electron,TString elid);
  // POG 
  void PogID(std::vector<snu::KElectron>& leptonColl, ID elid);


  void SelectElectrons(std::vector<snu::KElectron>& leptonColl, ID id,  TString elid, bool check_cc,bool check_cv,double  isomax_b,double isomax03_e,double dxymax_b,double dxymax_e,double dzmax_b,double dzmax_e, double dxysigmax,double dxysigmin,double ptcut=-999., double etacut=-999.);
  void SelectElectrons(std::vector<snu::KElectron>& leptonColl, TString elid, TString el_id,bool check_cc,bool check_cv,double  isomax_b,double isomax03_e,double dxymax_b,double dxymax_e,double dzmax_b,double dzmax_e, double dxysigmax,double dxysigmin,double ptcut=-999., double etacut=-999.);
  void SelectElectrons(std::vector<snu::KElectron>& leptonColl,TString elid, vector<pair<TString, TString> > vids, vector<pair<TString, float> > vidf, double ptcut=-999., double etacut=-999.);
  
  bool ElectronPass(snu::KElectron el, TString elid, TString el_id,bool check_cc,bool check_cv,double  isomax_b,double isomax03_e,double dxymax_b,double dxymax_e,double dzmax_b,double dzmax_e, double dxysigmax,double dxysigmin,double ptcut=-999., double etacut=-999.);
  bool ElectronPass(snu::KElectron el, TString elid, double ptcut=-999., double etacut=-999.);
  bool ElectronPass(snu::KElectron el, TString elid,vector<pair<TString, TString> > vids, vector<pair<TString, float> > vidf, double ptcut=-999., double etacut=-999.);
  

  //// IDS
  bool PassUserID(TString id,snu::KElectron el);
  bool PassUserID(TString id,snu::KElectron el, vector<pair<TString, TString> > vids, vector<pair<TString, float> > vidf);
  
  bool PassUserID(snu::KElectron el,TString id, TString el_id, bool check_cc,bool check_cv,double  isomax_b,double isomax03_e,double dxymax_b,double dxymax_e,double dzmax_b, double dzmax_e, double dxysigmax,double dxysigmin);

  bool PassID(snu::KElectron el, ID id);
    
};
#endif
