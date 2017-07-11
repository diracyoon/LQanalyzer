#include "ElectronSelection.h"


using namespace snu;

ElectronSelection::ElectronSelection(LQEvent ev) :  BaseSelection() {
  k_lqevent = ev;
  ElectronID = ELECTRON_POG_TIGHT;

};

ElectronSelection::~ElectronSelection() {};


void ElectronSelection::BasicSelection(std::vector<KElectron>& leptonColl , bool m_debug) {
  
  /// For filling SKTrees
  std::vector<KElectron> allelectrons = k_lqevent.GetElectrons();

  for (std::vector<KElectron>::iterator el = allelectrons.begin(); el!=allelectrons.end(); el++){
    
    if ( fabs(el->SCEta()) < eta_cut && el->Pt() >= pt_cut_min ){
      leptonColl.push_back(*el);
    }
    else if(m_debug) {
      if( fabs(el->SCEta()) >= eta_cut )  cout <<"BasicSelection::Fail Eta Cut" <<endl;
      if( el->Pt() < pt_cut_min )  cout <<"BasicSelection::Fail Pt Cut" <<endl;
    }
  }
  
}

void ElectronSelection::SkimSelection(std::vector<KElectron>& leptonColl, bool m_debug) {
  
  std::vector<KElectron> allelectrons = k_lqevent.GetElectrons();
  
  for (std::vector<KElectron>::iterator el = allelectrons.begin(); el!=allelectrons.end(); el++){


    if ( fabs(el->SCEta()) < eta_cut && el->Pt() >= pt_cut_min){
      leptonColl.push_back(*el);
    }
    else if(m_debug) {
      if( fabs(el->SCEta()) >= eta_cut )  cout <<"SKimSelection::Fail Eta Cut" <<endl;
      if( el->Pt() < pt_cut_min )  cout <<"SkimSelection::Fail Pt Cut" <<endl;
    }
  }
 
}



bool ElectronSelection::ElectronPass(snu::KElectron el, TString elid, double ptcut, double etacut){

  //// This is not used now by GetElectron function in AnalyzerCore. Only used by DataDriven class which has no access to selection map

  bool pass_selection = true;

  if (ptcut == -999.) ptcut = AccessFloatMap("ptmin",elid);
  if (etacut == -999.) etacut = AccessFloatMap("|etamax|",elid);
  if ( fabs(el.SCEta())>1.4442 && fabs(el.SCEta())<1.566 ) return false;;

  ElectronID = PassUserID(elid, el);
  if(!ElectronID)  pass_selection = false;
  
  if(!(fabs(el.SCEta()) < etacut)) pass_selection = false;
  if((el.Pt() < ptcut))   pass_selection = false;

  return pass_selection;
}


bool ElectronSelection::ElectronPass(snu::KElectron el, TString elid, TString el_id, bool check_cc,bool check_cv,double  isomax_b,double isomax03_e,double dxymax_b,double dxymax_e,double dzmax_b,double dzmax_e, double dxysigmax,double dxysigmin,double ptcut, double etacut){

  //// This is not used by any code but can be used by user IF they prefer
  /// This is Ferdinando style

  bool pass_selection = true;
  if ( fabs(el.SCEta())>1.4442 && fabs(el.SCEta())<1.566 ) return false;;

  ElectronID = PassUserID(el ,elid, el_id, check_cc,check_cv, isomax_b,isomax03_e,dxymax_b,dxymax_e,dzmax_b,dzmax_e, dxysigmax,dxysigmin);
  if(!ElectronID)  pass_selection = false;

  if(!(fabs(el.SCEta()) < etacut)) pass_selection = false;
  if((el.Pt() < ptcut))   pass_selection = false;

  return pass_selection;
}


bool ElectronSelection::ElectronPass(snu::KElectron el, TString elid, vector<pair<TString, TString> > vids, vector<pair<TString, float> > vidf, double ptcut, double etacut){

  //// This is new function that inputs vector instead of searching in map
  ///  This saves time vs previous functions
  int icut(0);
  /// in case user has overwridden pt/eta cuts
  if (ptcut == -999. || etacut == -999.){
    for(unsigned int iv=0; iv < vidf.size(); iv++){
      if(!Check(vidf[iv].second)) continue;
      if (vidf[iv].first =="ptmin") { icut++; if(ptcut == -999.)ptcut=vidf[iv].second;}
      if (vidf[iv].first =="|etamax|") {icut++;  if (etacut == -999.)etacut=vidf[iv].second;}
      if(icut ==2) break;
    }
  }

  bool pass_selection = true;
  if ( fabs(el.SCEta())>1.4442 && fabs(el.SCEta())<1.566 ) return false;;

  ElectronID = PassUserID(elid, el ,vids, vidf);
  if(!ElectronID)  pass_selection = false;

  if(!(fabs(el.SCEta()) < etacut)) pass_selection = false;

  if((el.Pt() < ptcut))   pass_selection = false;
  
  return pass_selection;
}




void ElectronSelection::SelectElectrons(std::vector<KElectron>& leptonColl, ID elid,TString el_id,bool check_cc,bool check_cv,double  isomax_b,double isomax03_e,double dxymax_b,double dxymax_e,double dzmax_b,double dzmax_e, double dxysigmax,double dxysigmin,double ptcut, double etacut){
  return SelectElectrons(leptonColl,GetString(elid),  el_id,check_cc,check_cv, isomax_b,isomax03_e,dxymax_b,dxymax_e,dzmax_b,dzmax_e, dxysigmax,dxysigmin, ptcut,etacut);
}
  
void ElectronSelection::SelectElectrons(std::vector<KElectron>& leptonColl, TString elid, TString el_id,bool check_cc,bool check_cv,double  isomax_b,double isomax03_e,double dxymax_b,double dxymax_e,double dzmax_b,double dzmax_e, double dxysigmax,double dxysigmin,double ptcut, double etacut){
  std::vector<KElectron> allelectrons = k_lqevent.GetElectrons();

  for (std::vector<KElectron>::iterator el = allelectrons.begin(); el!=allelectrons.end(); el++){
    
    if(ElectronPass(*el, elid, el_id,check_cc,check_cv, isomax_b,isomax03_e,dxymax_b,dxymax_e,dzmax_b,dzmax_e, dxysigmax,dxysigmin, ptcut,etacut) ) leptonColl.push_back(*el);


  }// end of el loop

  return;
}


void ElectronSelection::SelectElectrons(std::vector<KElectron>& leptonColl, TString elid,vector<pair<TString, TString> > vids, vector<pair<TString, float> > vidf, double ptcut, double etacut){
  
  std::vector<KElectron> allelectrons = k_lqevent.GetElectrons();
  
  for (std::vector<KElectron>::iterator el = allelectrons.begin(); el!=allelectrons.end(); el++){
    
    if(ElectronPass(*el, elid, vids,vidf, ptcut, etacut)) leptonColl.push_back(*el);

  }// end of el loop                                                                                                                                                                                                                                     

  return;
}




void ElectronSelection::Selection(std::vector<KElectron>& leptonColl , bool m_debug) {

  std::vector<KElectron> allelectrons = k_lqevent.GetElectrons();

  for (std::vector<KElectron>::iterator el = allelectrons.begin(); el!=allelectrons.end(); el++){

    bool pass_selection = true;

    // ID cut
    if(apply_ID){
      if     (GetString(k_id).Contains("POG"))             ElectronID = PassID(*el, k_id);
      else if(GetString(k_id).Contains("ELECTRON_HN_MVA")) ElectronID = PassID(*el, k_id);
      else ElectronID = PassUserID( *el,GetString(k_id) , GetString(k_id), apply_chargeconst, apply_convcut,relIsoBarrel_max,relIsoEndcap_max,dxyBarrel_max,dxyEndcap_max,dzBarrel_max,dzEndcap_max, 999.,999.);

      if(!ElectronID) {
        pass_selection = false;
        if(m_debug)cout << "Selection: Fail ID Cut" << endl;
      }
    }

    // HLT Emulation Cuts 
    if(apply_HLTSafeCut){
      if     ( HLTSafeLevel.Contains("CaloIdL_TrackIdL_IsoVL") && (!el->IsTrigMVAValid()) ) pass_selection=false;
      else if( HLTSafeLevel.Contains("WPLoose")                && (!el->PassHLTID())      ) pass_selection=false;
    }

    // SetRelIso. Default: PFRelIso03
    float reliso(0.);
    if     (apply_relisocut && RelIsoType.Contains("Default"))    reliso=el->PFRelIso(0.3); 
    else if(apply_relisocut && RelIsoType.Contains("PFRelIso04")) reliso=el->PFRelIso(0.4); 
    else                                                          reliso=el->PFRelIso(0.3); 

    if(apply_ptcut && !(el->Pt() >= pt_cut_min && el->Pt() < pt_cut_max)){
      pass_selection = false;
      if(m_debug)cout << "Selection: Fail Pt Cut" << endl;
    }
    if(apply_etacut && !(fabs(el->SCEta()) < eta_cut)) {
      pass_selection = false;
      if(m_debug)cout << "Selection: Fail Eta Cut" << endl;
    }

    //Whether to include EE-EBtransition region(Default: not)
    if(!apply_BETrRegIncl){
      if(fabs(el->SCEta())>1.4442 && fabs(el->SCEta())<1.566) pass_selection = false;
    }
    // Check charge consistancy between different detectors
    if(apply_chargeconst && !el->GsfCtfScPixChargeConsistency()){
      pass_selection = false;
      if(m_debug)cout << "Selection: Fail charge Cut" << endl;
    }
    // extra cut to reduce conversions
    // https://twiki.cern.ch/twiki/bin/view/CMS/ConversionTools
    if(apply_convcut && (!el->PassesConvVeto()) ) {
      pass_selection = false;
      if(m_debug)cout << "Selection: Fail Conversion Cut" << endl;
    }

    //d0 Significance
    if(apply_dxysigmin && !(fabs(el->dxySig()) >= dxySig_min)) pass_selection = false;
    if(apply_dxysigmax && !(fabs(el->dxySig()) <  dxySig_max)) pass_selection = false;


    if(apply_BESepCut){
      if(fabs(el->SCEta())<1.479){
        if(apply_relisocut && !(reliso < relIsoBarrel_max && reliso >= relIsoBarrel_min)) pass_selection=false;
        if(apply_dxycut    && !(fabs(el->dxy()) < dxyBarrel_max)) pass_selection=false;
        if(apply_dzcut     && !(fabs(el->dz()) < dzBarrel_max))   pass_selection=false;
      }
      else{
        if(apply_relisocut && !(reliso < relIsoEndcap_max && reliso >= relIsoEndcap_min)) pass_selection=false;
        if(apply_dxycut    && !(fabs(el->dxy()) < dxyEndcap_max)) pass_selection=false;
        if(apply_dzcut     && !(fabs(el->dz()) < dzEndcap_max))   pass_selection=false;
      }
    }
    else{
      if(apply_relisocut   && !(reliso < relIso_cut && reliso >= relIsoMIN_cut)){
        pass_selection = false;
        if(m_debug)cout << "Selection: Fail Isolation Cut" << endl;
      }
      if(apply_dzcut       && !(fabs(el->dz())<  dz_cut )) {
        pass_selection = false;
        if(m_debug)cout << "Selection: Fail dZ Cut" << endl;
      }
      if(apply_dxycut      && !( fabs(el->dxy())< dxy_cut )) {
        pass_selection = false;
        if(m_debug)cout << "Selection: Fail dxy Cut" << endl;
      }
    }


    if(pass_selection){
      leptonColl.push_back(*el);
    }

  }// end of el loop

  BaseSelection::reset();

  return;
}




bool ElectronSelection::PassUserID(TString id, snu::KElectron el, vector<pair<TString, TString> > vids , vector<pair<TString, float> > vidf ) {
  
  /// New function in 806 to check if electron passes ID. 
  /// now input vector of pairs

  int snuid = el.SNUID();
  bool pass_veto_noiso = false;
  bool pass_loose_noiso = false;
  bool pass_medium_noiso = false;
  bool pass_tight_noiso = false;
  if(snuid >= 1000){
    pass_tight_noiso = true;
    snuid = snuid - 1000;
  }
  if(snuid >= 100){
    pass_medium_noiso= true;
    snuid = snuid-100;
  }
  if(snuid >= 10){
    pass_loose_noiso=true;
    snuid = snuid-10;
  }
  if(snuid >= 1){
    pass_veto_noiso=true;
  }
  bool debug=false;
  
  bool pass_trigger_emulation=true;
  if(el.Pt() < 15.){
    /// special case that mva id not used as 15 GeV cut
    /// use single electron HLT cuts which are tighter than mva for 10-15 GeV bin
    if(!el.PassHLTID()) pass_trigger_emulation=false;
  }
  else {
    /// safe for DoubleEG triggers with CaloIdL_TrackIdL_IsoVL
    //https://twiki.cern.ch/twiki/bin/view/CMS/ChangesEGMHLTAlgo2014
    if(!el.IsTrigMVAValid()) pass_trigger_emulation=false;
  }

  for(unsigned int idel =0; idel < vids.size(); idel++){
    /// all cuts are not "false"
    if(vids[idel].second == "false") continue;

    if(vids[idel].first == "IsTight(POG)")  { 
      if(!pass_tight_noiso) {if(debug){ cout << "Fail tight" << endl;} return false;}
    }
    if(vids[idel].first == "IsMedium(POG)"){
      if(!pass_medium_noiso)  {if(debug){ cout << "Fail medium" << endl;} return false;}
    }
    if(vids[idel].first == "IsLoose(POG)"){
      if(!pass_loose_noiso) {if(debug){ cout << "Fail loose" << endl;} return false;}
    }
    if(vids[idel].first == "IsVeto(POG)"){
      if(!pass_veto_noiso)   {if(debug){ cout << "Fail veto" << endl;} return false;}
    }
    if(vids[idel].first == "GsfCtfScPix") {
      if( !el.GsfCtfScPixChargeConsistency()) {if(debug){ cout << "Fail charge" << endl;}  return false;}
    }
    if(vids[idel].first == "convveto")  {
      if(!el.PassesConvVeto()) {if(debug){ cout << "Fail convveto" << endl;}  return false;}
    }
    if(vids[idel].first == "IsTight(MVA)"){
      if(!pass_trigger_emulation)  {if(debug){ cout << "Fail MVA tight" << endl;} return false;}
      if(!el.PassTrigMVATight()){if(debug){ cout << "Fail MVA tight" << endl;} return false;}
      
    }
    if(vids[idel].first == "IsMedium(MVA)"){
      if(!pass_trigger_emulation)  {if(debug){ cout << "Fail MVA medium" << endl;} return false;}
      if(!el.PassTrigMVAMedium()) {if(debug){ cout << "Fail MVA medium" << endl;} return false;}
    }
    if(vids[idel].first == "IsHNLoose(MVA)"){
      if(!pass_trigger_emulation)  {if(debug){ cout << "Fail MVA medium" << endl;} return false;}
      if(!el.PassTrigMVAHNLoose()){if(debug){ cout << "Fail MVA tight" << endl;} return false;}
    }
    if(vids[idel].first == "IsHNTight(MVA)"){
      if(!pass_trigger_emulation)  {if(debug){ cout << "Fail MVA medium" << endl;} return false;}
      if(!el.PassTrigMVAHNTightv4()){if(debug){ cout << "Fail MVA tight" << endl;} return false;}
    }
    if(vids[idel].first == "IsGENTTight"){
      if(!pass_trigger_emulation)  {if(debug){ cout << "Fail MVA medium" << endl;} return false;}
      if(!el.PassTrigMVAGENTTight()){if(debug){ cout << "Fail MVA tight" << endl;} return false;}
    }
    
    if(vids[idel].first == "IsZZ(MVA)"){
      //if(!el.IsTrigMVAValid())  {if(debug){ cout << "Fail MVA medium" << endl;} return false;}
      if(!el.PassMVAZZ()) {if(debug){ cout << "Fail MVA medium" << endl;} return false;}
    }
  }

  

  LeptonRelIso = el.PFRelIso(0.3);

  bool checkUseMiniIso=false;
  if(id.Contains("miniiso")) checkUseMiniIso=true;
  if(checkUseMiniIso){
    LeptonRelIso = el.PFRelMiniIso();
  }

  for(unsigned int idel =0; idel < vidf.size(); idel++){
    if(!Check(vidf[idel].second)) continue;
    if(vidf[idel].first == "isomax03_b") {
      if(fabs(el.SCEta())<1.479 ){
	if((LeptonRelIso > vidf[idel].second))  {if(debug){ cout << "Fail iso: " << LeptonRelIso << " " << vidf[idel].second << endl;} return false;}
      }
    }
    if(vidf[idel].first == "isomax03_e") {
      if(fabs(el.SCEta())>=1.479 ){
	if((LeptonRelIso > vidf[idel].second))  {if(debug){ cout << "Fail iso" << LeptonRelIso << " " << vidf[idel].second << endl;}  return false;}
      }
    }
    if(vidf[idel].first == "|dxymax_b|") {
      if(fabs(el.SCEta())<1.479 ){
	if(fabs(el.dxy()) > vidf[idel].second) {if(debug){ cout << "Faildxy " << endl;} return false;}	
      }
    }
    if(vidf[idel].first == "|dxymax_e|") {
      if(fabs(el.SCEta())>=1.479 ){
	if(fabs(el.dxy()) > vidf[idel].second) {if(debug){ cout << "Fail dxy" << endl;} return false;}
      }
    }
    if(vidf[idel].first == "|dzmax_b|") {
      if(fabs(el.SCEta())<1.479 ){
	if(fabs(el.dz()) > vidf[idel].second) {if(debug){ cout << "Fail dz" << endl;} return false;}
      }
    }
    if(vidf[idel].first == "|dzmax_e|") {
      if(fabs(el.SCEta())>=1.479 ){
	if(fabs(el.dz()) > vidf[idel].second) {if(debug){ cout << "Fail dz" << endl;} return false;}
      }
    }
    if(vidf[idel].first == "|dxysigmin|") {
      if(fabs(el.dxySig()) < vidf[idel].second) { if(debug){ cout << "Fail dsximin"  << endl;} return false;}
    }
    if(vidf[idel].first == "|dxysigmax|") {
      if(fabs(el.dxySig()) >   vidf[idel].second) { if(debug){ cout << "Fail dsigmax"  << endl;} return false;}
    }
  }
  return true;
}


bool ElectronSelection::PassUserID(TString id, snu::KElectron el){


  float isomax_b = AccessFloatMap("isomax03_b",id);
  float isomax_e = AccessFloatMap("isomax03_e",id);
  float dxymax_b = AccessFloatMap("|dxymax_b|",id);
  float dxymax_e = AccessFloatMap("|dxymax_e|",id);
  float dzmax_b = AccessFloatMap("|dzmax_b|",id);
  float dzmax_e = AccessFloatMap("|dzmax_e|",id);

  float dxysigmax = AccessFloatMap("|dxysigmax|",id);
  float dxysigmin = AccessFloatMap("|dxysigmin|",id);

  bool checkisloose= (CheckCutString("IsLoose(POG)",id));
  bool checkisveto = (CheckCutString("IsVeto(POG)",id));
  bool checkismedium = (CheckCutString("IsMedium(POG)",id));
  bool checkistight  = (CheckCutString("IsTight(POG)",id));
  bool checkisMVAHNLoose = (CheckCutString("IsHNLoose(MVA)",id));
  bool checkisMVAHNTight = (CheckCutString("IsHNTight(MVA)",id));
  bool checkisMVAGENTTight = (CheckCutString("IsGENTTight",id));
  bool checkisMVAtight  = (CheckCutString("IsTight(MVA)",id));


  bool checkchargeconsy = (CheckCutString("GsfCtfScPix",id));
  bool convveto = (CheckCutString("convveto",id));
  bool checkdxysigmin  = CheckCutFloat("|dxysigmin|",id);
  bool checkdxysigmax  = CheckCutFloat("|dxysigmax|",id);
  
  LeptonRelIso = el.PFRelIso(0.3);
  bool checkUseMiniIso=false;
  if(id.Contains("miniiso")) checkUseMiniIso=true;

  if(checkUseMiniIso){
    LeptonRelIso = el.PFRelMiniIso();
  }
  bool pass_selection=true;

  int snuid = el.SNUID();
  bool pass_veto_noiso = false;
  bool pass_loose_noiso = false;
  bool pass_medium_noiso = false;
  bool pass_tight_noiso = false;
  if(snuid >= 1000){
    pass_tight_noiso = true;
    snuid = snuid - 1000;
  }
  if(snuid >= 100){
    pass_medium_noiso= true;
    snuid = snuid-100;
  }
  if(snuid >= 10){
    pass_loose_noiso=true;
    snuid = snuid-10;
  }    
  if(snuid >= 1){
    pass_veto_noiso=true;
  }


  bool debug=false;

  bool pass_trigger_emulation=true;
  if(el.Pt() < 15.){
    /// special case that mva id not used as 15 GeV cut                                                                                                         
    /// use single electron HLT cuts which are tighter than mva for 10-15 GeV bin                                                                               
    if(!el.PassHLTID()) pass_trigger_emulation=false;
  }
  else {
    /// safe for DoubleEG triggers with CaloIdL_TrackIdL_IsoVL                                                                                                  
    //https://twiki.cern.ch/twiki/bin/view/CMS/ChangesEGMHLTAlgo2014                                                                                            
    if(!el.IsTrigMVAValid()) pass_trigger_emulation=false;
  }


  if(checkisveto && !pass_veto_noiso)  {pass_selection = false;if(debug){ cout << "Failveto " << endl;}}
  if(checkisloose && !pass_loose_noiso)  {pass_selection = false;if(debug){ cout << "Failloose " << endl;}}
  if(checkismedium && !pass_medium_noiso)  {pass_selection = false;if(debug){ cout << "Fail medium" << endl;}}
  if(checkistight && !pass_tight_noiso)  {pass_selection = false;if(debug){ cout << "Fail tight" << endl;}}

  if(checkisMVAtight && !pass_trigger_emulation)  {pass_selection = false;if(debug){ cout << "Fail MVA tight" << endl;}}
  if(checkisMVAtight && !el.PassTrigMVATight()){pass_selection = false;if(debug){ cout << "Fail MVA tight" << endl;}}
  if(checkisMVAHNLoose && !pass_trigger_emulation)  {pass_selection = false;if(debug){ cout << "Fail MVA medium" << endl;}}
  if(checkisMVAHNLoose && !el.PassTrigMVAHNLoose()){pass_selection = false;if(debug){ cout << "Fail MVA medium" << endl;}}
  if(checkisMVAHNTight && !pass_trigger_emulation)  {pass_selection = false;if(debug){ cout << "Fail MVA medium" << endl;}}
  if(checkisMVAHNTight && !el.PassTrigMVAHNTightv4()){pass_selection = false;if(debug){ cout << "Fail MVA medium" << endl;}}
  if(checkisMVAGENTTight && !pass_trigger_emulation)  {pass_selection = false;if(debug){ cout << "Fail MVA medium" << endl;}}
  if(checkisMVAGENTTight && !el.PassTrigMVAGENTTight()){pass_selection = false;if(debug){ cout << "Fail MVA medium" << endl;}}
  


  if(convveto&& (!el.PassesConvVeto()) ){pass_selection = false;if(debug){ cout << "Fail convveto" << endl;}}
  if(checkchargeconsy &&  !el.GsfCtfScPixChargeConsistency()) {pass_selection = false;if(debug){ cout << "Fail charge" << endl;}}

  if(checkdxysigmin &&(fabs(el.dxySig()) < dxysigmin)) { pass_selection = false;if(debug){ cout << "Fail dsximin"  << endl;}}
  if(checkdxysigmax &&(fabs(el.dxySig()) > dxysigmax)) { pass_selection = false;if(debug){ cout << "Fail dsigmax"  << endl;}}

  if(fabs(el.SCEta())<1.479 ){  

    if((LeptonRelIso > isomax_b))  {pass_selection = false;if(debug){ cout << "Fail iso: " << LeptonRelIso << " " << isomax_b << endl;}}
    if(fabs(el.dxy()) > dxymax_b) {pass_selection = false;if(debug){ cout << "Faildxy " << endl;}}
    if(fabs(el.dz()) > dzmax_b) {pass_selection = false;if(debug){ cout << "Fail dz" << endl;}}
  }
  else{
    if((LeptonRelIso > isomax_e))  {pass_selection = false;if(debug){ cout << "Fail iso" << LeptonRelIso << " " << isomax_b << endl;}}
    if(fabs(el.dxy()) > dxymax_e) {pass_selection = false;if(debug){ cout << "Fail dxy" << endl;}}
    if(fabs(el.dz()) > dzmax_e) {pass_selection = false;if(debug){ cout << "Fail dz" << endl;}}
  }
  



  vector<pair<TString, TString> > vids =GetStringList(id);
  vector<pair<TString, float> > vidf = GetFloatList(id);
  bool check2 = PassUserID(id, el, vids, vidf);
  
  if(check2 != pass_selection) cout << "EL  IDs are not the same" << endl;

  return pass_selection;
}

bool ElectronSelection::PassUserID(snu::KElectron el, TString id, TString el_id, bool check_cc,bool check_cv,double  isomax_b,double isomax03_e,double dxymax_b,double dxymax_e,double dzmax_b, double dzmax_e, double dxysigmax,double dxysigmin){

  int snuid = el.SNUID();
  bool pass_veto_noiso = false;
  bool pass_loose_noiso = false;
  bool pass_medium_noiso = false;
  bool pass_tight_noiso = false;
  if(snuid >= 1000){
    pass_tight_noiso = true;
    snuid = snuid - 1000;
  }
  if(snuid >= 100){
    pass_medium_noiso= true;
    snuid = snuid-100;
  }
  if(snuid >= 10){
    pass_loose_noiso=true;
    snuid = snuid-10;
  }
  if(snuid >= 1){
    pass_veto_noiso=true;
  }


  bool debug=false;
  


  bool pass_trigger_emulation=true;
  if(el.Pt() < 15.){
    /// special case that mva id not used as 15 GeV cut                                                                                                         
    /// use single electron HLT cuts which are tighter than mva for 10-15 GeV bin                                                                               
    if(!el.PassHLTID()) pass_trigger_emulation=false;
  }
  else {
    /// safe for DoubleEG triggers with CaloIdL_TrackIdL_IsoVL                                                                                                  
    //https://twiki.cern.ch/twiki/bin/view/CMS/ChangesEGMHLTAlgo2014                                                                                            
    if(!el.IsTrigMVAValid()) pass_trigger_emulation=false;
  }




  if((el_id ==	"IsTight(POG)") &&  !pass_tight_noiso)  {if(debug){ cout << "Fail tight" << endl;}return false;}
  if((el_id ==	"IsMedium(POG)") && !pass_medium_noiso)  {if(debug){ cout << "Fail medium" << endl;}return false;}
  if((el_id ==	"IsLoose(POG)") && !pass_loose_noiso)  {if(debug){ cout << "Failloose " << endl;}return false;}
  if((el_id ==	"IsVeto(POG)") && !pass_veto_noiso)  {if(debug){ cout << "Failveto " << endl;}return false;}
  if((el_id ==	"IsMedium(MVA)")){
    if(!pass_trigger_emulation)  {if(debug){cout << "Fail MVA medium" << endl;} return false;}
    if(!el.PassTrigMVAMedium()) {if(debug){cout << "Fail MVA medium" << endl;} return false;}
  }
  if(el_id == "IsTight(MVA)") {
    if(!pass_trigger_emulation)  {if(debug){cout << "Fail MVA tight" << endl;} return false;}
    if(!el.PassTrigMVATight()) {if(debug){cout << "Fail MVA tight" << endl;} return false;}
  }
  if(el_id == "IsHNTight(MVA)") {
    if(!pass_trigger_emulation)  {if(debug){cout << "Fail MVA tight" << endl;} return false;}
    if(!el.PassTrigMVAHNTightv4()) {if(debug){cout << "Fail MVA tight" << endl;} return false;}
  }
  if(el_id == "IsHNLoose(MVA)") {
    if(!pass_trigger_emulation)  {if(debug){cout << "Fail MVA tight" << endl;} return false;}
    if(!el.PassTrigMVAHNLoose()) {if(debug){cout << "Fail MVA tight" << endl;} return false;}
  }
  if(el_id == "IsGENTTight") {
    if(!pass_trigger_emulation)  {if(debug){cout << "Fail MVA tight" << endl;} return false;}
    if(!el.PassTrigMVAGENTTight()) {if(debug){cout << "Fail MVA tight" << endl;} return false;}
  }

  
  if( check_cc &&  !el.GsfCtfScPixChargeConsistency()) { if(debug){ cout << "Fail charge" << endl;} return false;}
  if( check_cv &&  (!el.PassesConvVeto()) ){ if(debug){ cout << "Fail convveto" << endl;}return false;}

  
  LeptonRelIso = el.PFRelIso(0.3);
  
  if(id.Contains("miniaod")){
    LeptonRelIso = el.PFRelMiniIso();
  }
  
  
  if(fabs(el.SCEta())<1.479 ){
    if(Check(isomax_b)){
      if((LeptonRelIso > isomax_b))  {if(debug){ cout << "Fail iso: " << LeptonRelIso << " " << isomax_b << endl;} return false;}
    }
    if(Check(dxymax_b)){
      if(fabs(el.dxy()) > dxymax_b) {if(debug){ cout << "Faildxy " << endl;}return false;}
    }
    if(Check(dzmax_b)){
      if(fabs(el.dz()) > dzmax_b) {if(debug){ cout << "Fail dz" << endl;}return false;}    
    }
  }
  else{
    if(Check(isomax03_e)){
      if((LeptonRelIso > isomax03_e))  {if(debug){ cout << "Fail iso" << LeptonRelIso << " " << isomax03_e << endl;} return false;}
    }
    if(Check(dxymax_e)){
      if(fabs(el.dxy()) > dxymax_e) {if(debug){ cout << "Fail dxy" << endl;}return false;}
    }
    if(Check(dzmax_e)){
      if(fabs(el.dz()) > dzmax_e) {if(debug){ cout << "Fail dz" << endl;}return false;}
    }
  }
  if(Check(dxysigmax)){
    if((fabs(el.dxySig()) > dxysigmax)) { if(debug){ cout << "Fail dsigmax"  << endl;}return false;}     
  }
  if(Check(dxysigmin)){
    if((fabs(el.dxySig()) < dxysigmin)) { if(debug){ cout << "Fail dsximin"  << endl;}return false;}
  }

  return true;
}


ElectronSelection& ElectronSelection::operator= (const ElectronSelection& ms) {
  if(this != &ms){    
    BaseSelection::operator = (ms); 
    k_lqevent = ms.k_lqevent;  
    ElectronID = ms.ElectronID;

  }
  return *this;
};

ElectronSelection::ElectronSelection(const ElectronSelection& ms): 
  BaseSelection(ms){
  k_lqevent = ms.k_lqevent;
  ElectronID = ms.ElectronID;

};

bool ElectronSelection::PassID(snu::KElectron el, ID id){

  bool pass_selection=true;  bool debug=false;

  //SNUID : CutbasedPOGID without Isolation requirements

  if     ( id == ELECTRON_POG_VETO   && el.SNUID()<1    ){pass_selection = false; if(debug){ cout << "Failveto "  << endl;}}
  else if( id == ELECTRON_POG_LOOSE  && el.SNUID()<10   ){pass_selection = false; if(debug){ cout << "Failloose " << endl;}}
  else if( id == ELECTRON_POG_MEDIUM && el.SNUID()<100  ){pass_selection = false; if(debug){ cout << "Fail medium"<< endl;}}
  else if( id == ELECTRON_POG_TIGHT  && el.SNUID()<1000 ){pass_selection = false; if(debug){ cout << "Fail tight" << endl;}}
  else if( id == ELECTRON_POG_MVA_WP90 && (!el.PassNotrigMVAMedium()) ){ pass_selection = false; }
  else if( id == ELECTRON_POG_MVA_WP80 && (!el.PassNotrigMVATight())  ){ pass_selection = false; }
  else if( id == ELECTRON_HN_MVA_LOOSE && (!el.PassTrigMVAHNLoose())  ){ pass_selection = false; }
  else if( id == ELECTRON_HN_MVA_TIGHT && (!el.PassTrigMVAHNTightv4())  ){ pass_selection = false; }

  return pass_selection;

}
