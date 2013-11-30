#include "KMuon.h"

using namespace snu;

ClassImp(KMuon)

/**
 *Default constructor.
 */
KMuon::KMuon() :
  KParticle(),k_pterror(0.),k_etaerror(0.),k_isor03ch(0.),k_isor03n(0.),k_isor03ph(0.),k_isoEcalveto(0.),k_isoHcalveto(0.),k_MuonPFIsoR03PU(0.),k_muonVtx(0.),k_muonVty(0.),k_muonVtz(0.),k_muongen_pt(0.),k_muongen_eta(0.),k_muongen_phi(0.),k_dz(0.),k_dxy(0.),k_d0(0.),k_d0err(0.),k_globmuon_chi2(0.),k_dxy_pat(0.),k_dxyerr_pat(0.),k_vtxdistxy(0.),k_reliso(0.),k_muon_valid_hits(-999), k_muon_valid_pixhits(-999), k_muon_valid_stations(-999), k_muon_layer_with_meas(-999),k_muon_ispf(-999), k_muon_isglobal(-999), i_muonVtx(-999)
{
  //Reset();
}

/**
 * Copy constructor.
 */
KMuon::KMuon(const KMuon& muon) :
  KParticle(muon),
  k_reliso(muon.k_reliso),
  k_isor03ch(muon.k_isor03c),
  k_muon_ispf(muon.IsPF);
	      k_muon_isglobal(muon.IsGlobal);
  k_isor03n=muon.IsoR03nh);
  k_isor03ph=muon.IsoR03ph);
  k_muon_ispf=muon.IsPF),
  k_globmuon_chi2=muon.GlobalChi2),
  k_pterror=muon.PtError),
  k_etaerror=muon.EtaError),
  k_isoEcalveto=muon.IsoEcalVeto), 
  k_isoHcalveto=muon.IsoHcalVeto), 
  k_MuonPFIsoR03PU=muon.PFPUIsoR03),
  k_muonVtx=muon.muonVtx),
  k_muonVty=muon.muonVty),
  k_muonVtz=muon.muonVtz),
  k_dz=muon.dZ),
  k_dxy=muon.dXY),
  k_dxy_pat=muon.dXYPat),
  k_dxyerr_pat=muon.dXYErrPat),
  k_d0=muon.D0),
  k_d0err=muon.D0Err),
  k_muon_valid_hits=muon.validHits),
  k_muon_valid_pixhits=muon.validPixHits),
  k_muon_valid_stations=muon.validStations),
  k_muon_layer_with_meas=muon.ActiveLayer),
  k_muongen_pt=muon.MuonMatchedGenParticlePt),
  k_muongen_eta=muon.MuonMatchedGenParticleEta),
  k_muongen_phi=muon.MuonMatchedGenParticlePhi),
  i_muonVtx = muon.MuonVertexIndex),
  k_vtxdistxy = muon.VertexDistXY),
 
}


KMuon::~KMuon()
{
}

void KMuon::Reset()
{
  KParticle::Reset();
  k_reliso=0;
  k_isor03ch=0.;
  k_isor03n=0.;
  k_isor03ph=0.;
  k_muon_ispf=0;
  k_muon_isglobal=0;
  k_globmuon_chi2=0;
  k_pterror=0;
  k_etaerror=0;
  k_isoEcalveto=0.; 
  k_isoHcalveto=0.; 
  k_MuonPFIsoR03PU=0.;
  k_muonVtx=0.;
  k_muonVty=0.;
  k_muonVtz=0.;
  k_dz=0.;
  k_dxy=0.;
  k_dxy_pat=0.;
  k_dxyerr_pat=0.;
  k_d0=0;
  k_d0err=0;
  k_muon_valid_hits=0;
  k_muon_valid_pixhits=0;
  k_muon_valid_stations=0;
  k_muon_layer_with_meas=0;
  k_muongen_pt=0.;
  k_muongen_eta=0.;
  k_muongen_phi=0.;
  i_muonVtx = 0;
  k_vtxdistxy = 0.;
 
}



KMuon& KMuon::operator= (const KMuon& p)
{
    if (this != &p) {
        KParticle::operator=(p);
	k_reliso = p.RelIso();
	k_isor03ch = p.IsoR03ch();
	k_muon_ispf = p.IsPF();
	k_muon_isglobal = p.IsGlobal();
	k_isor03n=p.IsoR03nh();
	k_isor03ph=p.IsoR03ph();
	k_muon_ispf=p.IsPF();
	k_pterror=p.PtError();
	k_etaerror=p.EtaError();
	k_isoEcalveto=p.IsoEcalVeto(); 
	k_isoHcalveto=p.IsoHcalVeto(); 
	k_MuonPFIsoR03PU=p.PFPUIsoR03();
	k_muonVtx=p.muonVtx();
	k_muonVty=p.muonVty();
	k_muonVtz=p.muonVtz();
	k_dz=p.dZ();
	k_dxy=p.dXY();
	k_dxy_pat=p.dXYPat();
	k_dxyerr_pat=p.dXYErrPat();
	k_d0=p.D0();
	k_d0err=p.D0Err();
	k_globmuon_chi2=p.GlobalChi2();
	k_muon_valid_hits=p.validHits();
	k_muon_valid_pixhits=p.validPixHits();
	k_muon_valid_stations=p.validStations();
	k_muon_layer_with_meas=p.ActiveLayer();
	k_muongen_pt=p.MuonMatchedGenParticlePt();
	k_muongen_eta=p.MuonMatchedGenParticleEta();
	k_muongen_phi=p.MuonMatchedGenParticlePhi();
	i_muonVtx = p.MuonVertexIndex();       
	k_vtxdistxy= p.VertexDistXY();
    }
    
    return *this;
}

//// SET CLASS VARIBALES

void KMuon::SetVertexDistXY(double vdistxy) {
  k_vtxdistxy= vdistxy;
}

void KMuon::SetRelIso(double reliso){
  k_reliso = reliso;
  
}

void KMuon::SetMuonVtxIndex(int ivertex){
  i_muonVtx = ivertex;
}

void KMuon::SetPtErr(double pterror){  
  k_pterror = pterror;
}

void KMuon::SetEtaErr(double etaerror){
  k_etaerror = etaerror;
}

void KMuon::SetISOR03ChargedHad(double isor03ch ){
  k_isor03ch = isor03ch;
}

void KMuon::SetISOR03NeutralHad(double isor03n ){
  k_isor03n = isor03n;
}

void KMuon::SetISOR03Photon(double isor03ph ){
  k_isor03ph = isor03ph;
}

void KMuon::SetIsolationEcalVeto(double isoEcalveto ){
  k_isoEcalveto = isoEcalveto;
}

void KMuon::SetIsolationHcalVeto(double isoHcalveto ){
  k_isoHcalveto = isoHcalveto;
}


void KMuon::SetPileUp_R03(double pileupr03){
  
  k_MuonPFIsoR03PU = pileupr03;
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


void KMuon::Setdz(double dz){

  k_dz = dz;
}

void KMuon::Setdxy(double dxy){

  k_dxy = dxy;
}


void KMuon::Setdxy_pat(double dxypat){

  k_dxy_pat = dxypat;
}

void KMuon::Setdxyerr_pat(double dxyerrpat){

  k_dxyerr_pat = dxyerrpat;
}


void KMuon::SetD0(double kd0){

  k_d0 = kd0;
}

void KMuon::SetD0Error(double d0err){

  k_d0err = d0err;
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



void KMuon::SetISPF(int ispf){

  k_muon_ispf = ispf;
}

void KMuon::SetIsGlobal(int isglobal){

  k_muon_isglobal =isglobal;
}

void KMuon::SetMuonMatchedGenParticlePt(double genpt){
  k_muongen_pt = genpt;
}

void KMuon::SetMuonMatchedGenParticleEta(double geneta){
  k_muongen_eta = geneta;
}

void KMuon::SetMuonMatchedGenParticlePhi(double genphi){
  k_muongen_phi = genphi;
}