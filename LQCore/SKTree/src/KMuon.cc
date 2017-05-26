#include "KMuon.h"

using namespace snu;

ClassImp(KMuon)

/**
 *Default constructor.
 */
KMuon::KMuon() :
  KParticle(),
  k_dz(-999.),
  k_dxy(-999.),
  k_dxy_sig(-999.),
  k_globmuon_chi2(-999.),
  k_muonVtx(-999.),
  k_muonVty(-999.),
  k_muonVtz(-999.),
  
  k_muon_valid_hits(-999), 
  k_muon_valid_pixhits(-999),
  k_muon_valid_stations(-999),
  k_muon_layer_with_meas(-999),
  k_muon_ispf(0),
  k_muon_isglobal(0),
  k_muon_istracker(0),

  muon_pt_up(-999.), 
  muon_pt_down(-999.),
  muon_maod_pt(-999.),
  k_muon_reliso03(-999.),
  k_muon_reliso04(-999.),
  k_muon_relminiiso(-999.),
  k_muon_maod_reliso03(-999.),
  k_muon_maod_reliso04(-999.),

  k_roch_pt(-999.),
  k_roch_eta(-999.),
  k_roch_phi(-999.),
  k_roch_m(-999.),
  k_roch_e(-999.),

  k_isloose(0),
  k_istight(0),
  k_matched(0),
  k_is_cf(0),
  k_is_conv(0),
  k_is_fromtau(0),
  k_ismedium(0),
  k_issoft(0),
  

  k_mother_pdgid(-1),
  k_mc_pdgid(-1),
  k_mother_index(-1),
  k_mc_index(-1),
  
  k_trig_match(""),
  k_corrected_rc(false),
  k_mctype(-999),
  k_isprompt(false)

{
  //Reset();
}

/**
 * Copy constructor.
 */
KMuon::KMuon(const KMuon& muon) :
  KParticle(muon),
  k_dz(muon.k_dz),
  k_dxy(muon.k_dxy),
  k_dxy_sig(muon.k_dxy_sig),
  k_globmuon_chi2(muon.k_globmuon_chi2),
  k_muonVtx(muon.k_muonVtx),
  k_muonVty(muon.k_muonVty),
  k_muonVtz(muon.k_muonVtz),

  k_muon_valid_hits(muon.k_muon_valid_hits),
  k_muon_valid_pixhits(muon.k_muon_valid_pixhits),
  k_muon_valid_stations(muon.k_muon_valid_stations),
  k_muon_layer_with_meas(muon.k_muon_layer_with_meas),

  k_muon_ispf(muon.k_muon_ispf),
  k_muon_isglobal(muon.k_muon_isglobal),
  k_muon_istracker(muon.k_muon_istracker),

  muon_pt_up(muon.muon_pt_up),
  muon_pt_down(muon.muon_pt_down),
  muon_maod_pt(muon.muon_maod_pt),
  k_muon_reliso03(muon.k_muon_reliso03),
  k_muon_reliso04(muon.k_muon_reliso04),
  k_muon_relminiiso(muon.k_muon_relminiiso),
  k_muon_maod_reliso03(muon.k_muon_maod_reliso03),
  k_muon_maod_reliso04(muon.k_muon_maod_reliso04),


  k_roch_pt(muon.k_roch_pt),
  k_roch_eta(muon.k_roch_eta),
  k_roch_phi(muon.k_roch_phi),
  k_roch_m(muon.k_roch_m),
  k_roch_e(muon.k_roch_e),

  k_isloose(muon.k_isloose),
  k_istight(muon.k_istight), 
  k_matched(muon.k_matched), 
  k_is_cf(muon.k_is_cf),
  k_is_conv(muon.k_is_conv),
  k_is_fromtau(muon.k_is_fromtau),
  k_ismedium(muon.k_ismedium),
  k_issoft(muon.k_issoft),

  k_mother_pdgid(muon.k_mother_pdgid),
  k_mc_pdgid(muon.k_mc_pdgid),
  k_mother_index(muon.k_mother_index),
  k_mc_index(muon.k_mc_index),

  k_trig_match(muon.k_trig_match),
  k_corrected_rc(muon.k_corrected_rc),
  k_mctype(muon.k_mctype),
  k_isprompt(muon.k_isprompt)
{
}
  


KMuon::~KMuon()
{
}

void KMuon::Reset()
{
  KParticle::Reset();
  k_dz=-999.;
  k_dxy=-999.;
  k_dxy_sig=-999.;
  k_globmuon_chi2=-999;
  k_muonVtx=-999.;
  k_muonVty=-999.;
  k_muonVtz=-999.;
  
  k_muon_valid_hits=-999;
  k_muon_valid_pixhits=-999;
  k_muon_valid_stations=-999;
  k_muon_layer_with_meas=-999;

  
  k_muon_ispf=false;
  k_muon_isglobal=false;
  k_muon_istracker=false;
  
  muon_pt_up=0.;
  muon_pt_down=0.;
  muon_maod_pt=0.;
  k_muon_reliso03=-999.;
  k_muon_reliso04=-999.;
  k_muon_relminiiso=-999.;
  k_muon_maod_reliso03=-999.;
  k_muon_maod_reliso04=-999.;
  k_roch_pt=-999.;
  k_roch_eta=-999.;
  k_roch_phi=-999.;
  k_roch_m=-999.;
  k_roch_e=-999.;



  k_isloose=false;
  k_istight=false;
  k_matched=false;
  k_is_cf=false;
  k_is_conv=false;
  k_is_fromtau=false;
  k_ismedium=false;
  k_issoft=false;

  k_mother_pdgid=-1;
  k_mc_pdgid=1;
  k_mother_index=-1;
  k_mc_index=-1;

  k_trig_match= "";
  k_corrected_rc=false;
  k_mctype=-999;
  k_isprompt=false;
}



KMuon& KMuon::operator= (const KMuon& p)
{
    if (this != &p) {
        KParticle::operator=(p);
	k_muon_ispf = p.IsPF();
	k_muon_isglobal = p.IsGlobal();
	k_muon_istracker = p.IsTracker();
	k_dz=p.dZ();
	k_dxy=p.dXY();
	k_dxy_sig=p.dXYSig();
	k_globmuon_chi2=p.GlobalChi2();
	k_muon_valid_hits=p.validHits();
	k_muon_valid_pixhits=p.validPixHits();
	k_muon_valid_stations=p.validStations();
	k_muon_layer_with_meas=p.ActiveLayer();
	//i_muonVtx = p.MuonVertexIndex();       
	muon_pt_up = p.PtShiftedUp();
	muon_pt_down = p.PtShiftedDown();
	k_isloose = p.IsLoose();
	k_istight = p.IsTight();
	k_issoft = p.IsSoft();
        k_ismedium = p.IsMedium();
	k_matched = p.MCMatched();
	k_is_cf=p.MCIsCF();
	k_is_conv=p.MCIsFromConversion();

	k_is_fromtau=p.MCFromTau();
	k_mother_pdgid=p.MotherPdgId();
	k_mc_pdgid=p.MCMatchedPdgId();
	k_mother_index=p.MotherTruthIndex();
	k_mc_index=p.MCTruthIndex();

	muon_maod_pt = p.MiniAODPt();

	k_muon_reliso03 = p.RelIso03();
	k_muon_reliso04 = p.RelIso04();
	k_muon_relminiiso = p.RelMiniIso();
	
	k_muon_maod_reliso03 = p.RelMiniAODIso03();
        k_muon_maod_reliso04 = p.RelMiniAODIso04();
	
	k_roch_pt=p.RochPt();
	k_roch_eta=p.RochEta();
	k_roch_phi=p.RochPhi();
	k_roch_m=p.RochM();
	k_roch_e=p.RochE();

	k_muonVtx=p.muonVtx();
	k_muonVty=p.muonVty();
	k_muonVtz=p.muonVtz();
	k_trig_match= p.TrigMatch();
	k_corrected_rc=p.IsRochesterCorrected();
	k_mctype = p.GetType();
	k_isprompt=p.IsPromptFlag();
    }
    
    return *this;
}

//// SET CLASS VARIBALES



std::string KMuon::Type() const
{
    return "KMuon";
}

void KMuon::SetIsPromptFlag(bool pflag){

  k_isprompt=pflag;
}

bool KMuon::TriggerMatched(TString path){
  TString trig = k_trig_match;
  if(trig.Contains(path)) return true;
  return false;
}

void KMuon::SetType(int type){
  k_mctype= type;
}
void KMuon::SetIsFromTau(bool istau){
  k_is_fromtau=istau;
}

void KMuon::SetIsRochesterCorrected(bool rccorr){
  k_corrected_rc=rccorr;
}

void KMuon::SetIsChargeFlip(bool iscf){
  k_is_cf=iscf;
}

void KMuon::SetIsPhotonConversion(bool isconv){
  k_is_conv=isconv;
}


void KMuon::SetRochM(double m){
  
  k_roch_m=m;
}

void KMuon::SetRochPt(double pt){

  k_roch_pt=pt;
}

void KMuon::SetRochEta(double eta){

  k_roch_eta=eta;
}

void KMuon::SetRochPhi(double phi){

  k_roch_phi=phi;
}

void KMuon::SetRochE(double e){

  k_roch_e=e;
}


void KMuon::SetMotherPdgId(int type){
  k_mother_pdgid=type;
}

void KMuon::SetMCMatchedPdgId(int type){
  k_mc_pdgid=type;
}

void KMuon::SetMotherTruthIndex(int mindex){
  k_mother_index=mindex;
}

void KMuon::SetMCTruthIndex(int tindex){
  k_mc_index=tindex;
}


//void KMuon::SetMuonVtxIndex(int ivertex){
//  i_muonVtx = ivertex;
//}


void KMuon::Setdz(double dz){

  k_dz = dz;
}

void KMuon::Setdxy(double dxy){

  k_dxy = dxy;
}

void KMuon::Setdxy_sig(double dxysig){

  k_dxy_sig = dxysig;
}


void KMuon::SetTrackVx(double vtx){
  
  k_muonVtx = vtx;
}

void KMuon::SetTrackVy(double vty){

  k_muonVty = vty;
}

void KMuon::SetTrackVz(double vtz){

  k_muonVtz = vtz;
}

void KMuon::SetGlobalchi2(double glob_chi2){

  k_globmuon_chi2 = glob_chi2;
}

void KMuon::SetValidHits(int validhits){

  k_muon_valid_hits = validhits;
}

void KMuon::SetPixelValidHits(int valid_pix_hits){

  k_muon_valid_pixhits = valid_pix_hits;
}

void KMuon::SetValidStations(int validstations){

  k_muon_valid_stations = validstations;
}

void KMuon::SetLayersWithMeasurement(int layer_with_meas){

  k_muon_layer_with_meas = layer_with_meas;
}



void KMuon::SetISPF(bool ispf){

  k_muon_ispf = ispf;
}

void KMuon::SetIsGlobal(bool isglobal){

  k_muon_isglobal =isglobal;
}

void KMuon::SetIsTracker(bool istracker){

  k_muon_istracker =istracker;
}


void KMuon::SetShiftedEUp(double pt_up){
  muon_pt_up= pt_up;
}

void KMuon::SetShiftedEDown(double pt_down){
  muon_pt_down= pt_down;
}


void KMuon::SetIsLoose(bool isloose){
  k_isloose= isloose;
}
void KMuon::SetIsTight(bool istight){
  k_istight= istight;
}

void KMuon::SetIsMedium(bool ismedium){
  k_ismedium= ismedium;
}
void KMuon::SetIsSoft(bool issoft){
  k_issoft= issoft;
}


void KMuon::SetMCMatched(bool ismatched){
  k_matched= ismatched;
}


void KMuon::SetRelIso(double cone, double iso){
  if(cone == 0.3) k_muon_reliso03= iso;
  else if(cone ==  0.4) k_muon_reliso04= iso;

}

void KMuon::SetRelMiniIso( double iso){
  k_muon_relminiiso = iso;
}


void KMuon::SetMiniAODRelIso (double cone, double iso){
  if(cone == 0.3) k_muon_maod_reliso03= iso;
  else if(cone ==  0.4) k_muon_maod_reliso04= iso;

}
void KMuon::SetMiniAODPt(double maodpt){
  muon_maod_pt=maodpt;
}



void KMuon::SetTrigMatch(TString match){
  k_trig_match= match;
}
