#include "KElectron.h"

#include <iomanip>
#include <iostream>

using namespace std;
using namespace snu;

ClassImp(KElectron)

/**
 *Default constructor.
 */
KElectron::KElectron() :
KParticle()
{
  k_sceta=-999;
  k_dxy=-999;
  k_dxy_sig=-999;
  k_dz=-999;
  k_gsf_ctscpix_charge=false;
  k_hasmatchconvphot=false; 
  k_pf_chargedhad_iso03=-999;
  k_pf_photon_iso03=-999;
  k_pf_neutral_iso03=-999;
  k_pf_chargedhad_iso04=-999;
  k_pf_photon_iso04=-999;
  k_pf_neutral_iso04=-999;
  k_rel_iso03=-999;
  k_rel_iso04=-999;
  k_rel_miniiso=-999;
  k_abs_iso03=-999;
  k_abs_iso04=-999;
  k_pt_shifted_up=-999;
  k_pt_shifted_down=-999;
  pass_veto=false;
  pass_loose=false;
  pass_medium=false;
  pass_tight=false;
  pass_hltid=false;
  pass_heep=false;
  pass_trigmva_medium=false;
  pass_trigmva_tight=false;
  pass_notrigmva_medium=false;
  pass_notrigmva_tight=false;
  pass_notrigmva_zz=false;
  k_mc_matched=false;
  k_is_cf=false;
  k_is_conv=false;
  k_is_fromtau=false;
  k_mother_pdgid=-1;
  k_mc_pdgid=-1;
  k_mother_index=-1;
  k_mc_index=-1;
  k_isPF=false;
  k_istrigmvavalid=false;
  snu_id = -999;
  k_trkvx=-999;
  k_trkvy=-999;
  k_trkvz=-999;
  k_gsf_ctscpix_charge=false;
  k_trig_match="";
  k_eltype=-999;
  k_mva=-999;
  k_zzmva=-999;
  k_missing_hits=-999;
  k_smearfactor=-999.;
  k_in_conv=false;
  k_isprompt=false;
  
}

/**
 * Copy constructor.
 */
KElectron::KElectron(const KElectron& el) :
    KParticle(el)
{
  
  k_sceta= el.SCEta();
  k_dxy= el.dxy();
  k_dxy_sig= el.dxySig();
  k_dz= el.dz();
  k_gsf_ctscpix_charge= el.GsfCtfScPixChargeConsistency();
  k_hasmatchconvphot= el.PassesConvVeto();
  k_pf_chargedhad_iso03= el.PFChargedHadronIso(0.3);
  k_pf_photon_iso03 = el.PFPhotonIso(0.3);
  k_pf_neutral_iso03= el.PFNeutralHadronIso(0.3);
  k_pf_chargedhad_iso04= el.PFChargedHadronIso(0.4);
  k_pf_photon_iso04= el.PFPhotonIso(0.4);
  k_pf_neutral_iso04= el.PFNeutralHadronIso(0.4);
  k_rel_iso03=el.PFRelIso(0.3);
  k_rel_iso04=el.PFRelIso(0.4);
  k_rel_miniiso=el.PFRelMiniIso();
  k_abs_iso03=el.PFAbsIso(0.3);
  k_abs_iso04=el.PFAbsIso(0.4);
  k_pt_shifted_up=el.PtShiftedUp();
  k_pt_shifted_down=el.PtShiftedDown();
  pass_veto=el.PassVeto();
  pass_loose=el.PassLoose();
  pass_medium=el.PassMedium();
  pass_tight=el.PassTight();
  pass_hltid=el.PassHLTID();
  pass_heep=el.PassHEEP();
  pass_trigmva_medium=el.PassTrigMVAMedium();
  pass_trigmva_tight=el.PassTrigMVATight(); 
  pass_notrigmva_medium=el.PassNotrigMVAMedium();
  pass_notrigmva_tight=el.PassNotrigMVATight();
  pass_notrigmva_zz=el.PassMVAZZ();
  k_mc_matched=el.MCMatched();
  k_is_cf=el.MCIsCF();
  k_is_conv=el.MCIsFromConversion();
  k_is_fromtau=el.MCFromTau();
  k_mother_pdgid=el.MotherPdgId();
  k_mc_pdgid=el.MCMatchedPdgId();
  k_mother_index=el.MotherTruthIndex();
  k_mc_index=el.MCTruthIndex();
  k_isPF=el.IsPF();
  k_istrigmvavalid=el.IsTrigMVAValid();
  snu_id = el.SNUID();
  k_trkvx= el.TrkVx();
  k_trkvy= el.TrkVy();
  k_trkvz= el.TrkVz();
  k_gsf_ctscpix_charge= el.GsfCtfScPixChargeConsistency();
  k_trig_match= el.TrigMatch();
  k_eltype=el.GetType();
  k_mva=el.MVA();
  k_zzmva=el.ZZMVA();
  k_missing_hits=el.MissingHits();
  k_smearfactor=el.SmearFactor();
  k_in_conv=el.MCIsExternalConversion();
  k_isprompt=el.IsPromptFlag();

}



KElectron::~KElectron()
{
}

void KElectron::Reset()
{
  KParticle::Reset();
  k_sceta=-999;
  k_dxy=-999;
  k_dxy_sig=-999;
  k_dz=0;
  k_gsf_ctscpix_charge=false;
  k_hasmatchconvphot=false;
  k_pf_chargedhad_iso03=-999;
  k_pf_photon_iso03=-999;
  k_pf_neutral_iso03=-999;
  k_pf_chargedhad_iso04=-999;
  k_pf_photon_iso04=-999;
  k_pf_neutral_iso04=-999;
  k_rel_iso03=-999;
  k_rel_iso04=-999;
  k_pt_shifted_up=-999;
  k_pt_shifted_down=-999;
  pass_veto=false;
  pass_loose=false;
  pass_medium=false;
  pass_tight=false;
  pass_hltid=false;
  pass_heep=false;
  pass_trigmva_medium=false;
  pass_trigmva_tight=false;
  pass_notrigmva_medium=false;
  pass_notrigmva_tight=false;
  pass_notrigmva_zz=false;
  k_mc_matched=false;
  k_is_cf=false;
  k_is_conv=false;
  k_is_fromtau=false;
  k_mother_pdgid=-1;
  k_mc_pdgid=1;
  k_mother_index=-1;
  k_mc_index=-1;
  k_isPF=false;
  k_istrigmvavalid=false;
  snu_id = -999;
  k_trkvx=-999;
  k_trkvy=-999;
  k_trkvz=-999;
  k_abs_iso03=-999;
  k_abs_iso04=-999;
  k_gsf_ctscpix_charge=false;
  k_trig_match="";
  k_eltype=-999;
  k_mva=-999;
  k_zzmva=-999;
  k_missing_hits=-999;
  k_smearfactor=-999.;
  k_in_conv=false;
  k_isprompt=false;
  
}



KElectron& KElectron::operator= (const KElectron& p)
{
  if (this != &p) {
    KParticle::operator=(p);
    k_sceta= p.SCEta();
    k_dxy= p.dxy();
    k_dxy_sig= p.dxySig();
    k_dz= p.dz();
    k_gsf_ctscpix_charge= p.GsfCtfScPixChargeConsistency();
    k_hasmatchconvphot= p.PassesConvVeto();
    k_pf_chargedhad_iso03= p.PFChargedHadronIso(0.3);
    k_pf_photon_iso03 = p.PFPhotonIso(0.3);
    k_pf_neutral_iso03= p.PFNeutralHadronIso(0.3);
    k_pf_chargedhad_iso04= p.PFChargedHadronIso(0.4);
    k_pf_photon_iso04= p.PFPhotonIso(0.4);
    k_pf_neutral_iso04= p.PFNeutralHadronIso(0.4);
    k_rel_iso03=p.PFRelIso(0.3);
    k_rel_iso04=p.PFRelIso(0.4);
    k_rel_miniiso=p.PFRelMiniIso();
    k_pt_shifted_up=p.PtShiftedUp();
    k_pt_shifted_down=p.PtShiftedDown();
    pass_veto=p.PassVeto();
    pass_loose=p.PassLoose();
    pass_medium=p.PassMedium();
    pass_tight=p.PassTight();
    pass_hltid=p.PassHLTID();
    pass_heep=p.PassHEEP();
    pass_trigmva_medium=p.PassTrigMVAMedium();
    pass_trigmva_tight=p.PassTrigMVATight();
    pass_notrigmva_medium=p.PassNotrigMVAMedium();
    pass_notrigmva_tight=p.PassNotrigMVATight();
    pass_notrigmva_zz=p.PassMVAZZ();
    k_mc_matched=p.MCMatched();
    k_is_cf =p.MCIsCF();
    k_is_conv =p.MCIsFromConversion();
    k_is_fromtau=p.MCFromTau();
    k_mother_pdgid=p.MotherPdgId();
    k_mc_pdgid=p.MCMatchedPdgId();
    k_mother_index=p.MotherTruthIndex();
    k_mc_index=p.MCTruthIndex();
    k_isPF=p.IsPF();
    k_istrigmvavalid=p.IsTrigMVAValid();
    snu_id = p.SNUID();
    k_trkvx= p.TrkVx();
    k_trkvy= p.TrkVy();
    k_trkvz= p.TrkVz();
    k_abs_iso03=p.PFAbsIso(0.3);
    k_abs_iso04=p.PFAbsIso(0.4);
    k_trig_match= p.TrigMatch();
    k_eltype=p.GetType();
    k_mva=p.MVA();
    k_zzmva=p.ZZMVA();
    k_missing_hits=p.MissingHits();
    k_smearfactor=p.SmearFactor();
    k_in_conv=p.MCIsExternalConversion();
    k_isprompt=p.IsPromptFlag();
  }
  
  return *this;
}


void KElectron::SetIsMCExternalConversion(bool isconv){
  k_in_conv=isconv;
}


void KElectron::SetIsPromptFlag(bool ispromptflag){
  k_isprompt=ispromptflag;
}

void KElectron::SetSmearFactor(double smeare){
  
  k_smearfactor=smeare;

}

void KElectron::SetMissingHits(int mhits){
  k_missing_hits=mhits;
}

void KElectron::SetMVA(double mva){
  k_mva=mva;
}

void KElectron::SetZZMVA(double zzmva){
  k_zzmva=zzmva;
}


void  KElectron::SetType(int type){
  k_eltype= type;
}

bool KElectron::TriggerMatched(TString path){
  TString trig = k_trig_match;
  if(trig.Contains(path)) return true;
  return false;
}


Bool_t KElectron::PassTrigMVAHNTight() const{
  if(this->Pt() < 50.) {
    if(k_mva > 0.95) return true;
  }
  else {
    if((fabs(this->SCEta()) < 0.8) && k_mva > 0.76) return true;
    if((fabs(this->SCEta())  > 0.8) &&(fabs(this->SCEta())  < 1.479)  && k_mva > 0.72) return true;
    if((fabs(this->SCEta())  < 2.5) &&(fabs(this->SCEta())  > 1.479) && k_mva > 0.70) return true;
  }
  return false;

}


Bool_t KElectron::PassTrigMVAGENTTight() const{
  if(this->Pt() < 15.) {
    if((fabs(this->SCEta()) < 0.8) && k_mva > 0.77) return true;
    if((fabs(this->SCEta())  > 0.8) &&(fabs(this->SCEta())  < 1.479)  && k_mva > 0.56) return true;
    if((fabs(this->SCEta())  < 2.5) &&(fabs(this->SCEta())  > 1.479) && k_mva > 0.48) return true;

  }
  else {
    if((fabs(this->SCEta()) < 0.8) && k_mva > 0.52) return true;
    if((fabs(this->SCEta())  > 0.8) &&(fabs(this->SCEta())  < 1.479)  && k_mva > 0.11) return true;
    if((fabs(this->SCEta())  < 2.5) &&(fabs(this->SCEta())  > 1.479) && k_mva > -0.01) return true;
  }
  return false;

}

Bool_t KElectron::PassTrigMVAHNLoose() const{
  if((fabs(this->SCEta()) < 0.8) && k_mva > -0.02) return true;
  if((fabs(this->SCEta())  > 0.8) &&(fabs(this->SCEta())  < 1.479)  && k_mva > -0.52) return true;
  if((fabs(this->SCEta())  < 2.5) &&(fabs(this->SCEta())  > 1.479) && k_mva > -0.52) return true;
  
  return false;
}



float KElectron::ScaleFactor(const std::string& name, int sign) const {
  
  float fsign =1.;
  if(sign == 0) fsign =0.;
  if(sign == -1) fsign =-1.;


  //need to fill errors. 
  float eta = fabs(this->Eta());
  if (name == "mvaEleID-Spring15-25ns-Trig-V1-wp90"){
    return 1.;
  }
  /// ID SF from  https://arun-afs.web.cern.ch/arun-afs/Fits_Data_ID_76X_AbsEta/CutBasedID_TightWP_76X_18Feb.txt
  else if (name == "cutBasedElectronID-Spring15-25ns-V1-standalone-loose"){
    if (this->Pt()>10. && this->Pt() <= 20.){
      if      ( eta <=  0.800)    return (0.661/0.655)  + fsign*(0.02);
      else if ( eta <=  1.444)    return (0.684/0.635)  + fsign*(0.06);
      else if ( eta <=  1.566)    return (0.480/0.448)  + fsign*(0.06);
      else if ( eta <=  2.000)    return (0.612/0.596)  + fsign*(0.06);
      else if ( eta <=  2.500)    return (0.667/0.642)  + fsign*(0.06);
      else return 1.;
    }
    else  if (this->Pt()>20. && this->Pt() <= 30.){
      if      ( eta <=  0.800)    return (0.779/0.787)  + fsign*(0.02);
      else if ( eta <=  1.444)    return (0.762/0.769)  + fsign*(0.06);
      else if ( eta <=  1.566)    return (0.502/0.514)  + fsign*(0.06);
      else if ( eta <=  2.000)    return (0.731/0.747)  + fsign*(0.06);
      else if ( eta <=  2.500)    return (0.764/0.760)  + fsign*(0.06);

    }   
    else  if (this->Pt()>30. && this->Pt() <= 40.){
      if      ( eta <=  0.800)    return (0.868/0.877)  + fsign*(0.02);
      else if ( eta <=  1.444)    return (0.862/0.870)  + fsign*(0.06);
      else if ( eta <=  1.566)    return (0.668/0.679)  + fsign*(0.06);
      else if ( eta <=  2.000)    return (0.836/0.845)  + fsign*(0.06);
      else if ( eta <=  2.500)    return (0.837/0.834)  + fsign*(0.06);
      else return 1.;
    }   
    else  if (this->Pt()>40. && this->Pt() <= 50.){
      if      ( eta <=  0.800)    return (0.913/0.922)  + fsign*(0.02);
      else if ( eta <=  1.444)    return (0.912/0.921)  + fsign*(0.06);
      else if ( eta <=  1.566)    return (0.794/0.809)  + fsign*(0.06);
      else if ( eta <=  2.000)    return (0.894/0.889)  + fsign*(0.06);
      else if ( eta <=  2.500)    return (0.877/0.865)  + fsign*(0.06);
      else return 1.;
    }   
    else  if (this->Pt()>50. && this->Pt() <= 200.){
      if      ( eta <=  0.800)    return (0.927/0.937)  + fsign*(0.02);
      else if ( eta <=  1.444)    return (0.925/0.935)  + fsign*(0.06);
      else if ( eta <=  1.566)    return (0.806/0.819)  + fsign*(0.06);
      else if ( eta <=  2.000)    return (0.914/0.901)  + fsign*(0.06);
      else if ( eta <=  2.500)    return (0.890/0.876)  + fsign*(0.06);
      else return 1.;
    }   
  }
  else if (name == "cutBasedElectronID-Spring15-25ns-V1-standalone-medium"){
    if (this->Pt()>10. && this->Pt() <= 20.){
      if      ( eta <=  0.800)    return (0.555/0.551)  + fsign*(0.02);
      else if ( eta <=  1.444)    return (0.556/0.510)  + fsign*(0.06);
      else if ( eta <=  1.566)    return (0.352/0.324)  + fsign*(0.06);
      else if ( eta <=  2.000)    return (0.443/0.450)  + fsign*(0.06);
      else if ( eta <=  2.500)    return (0.524/0.506)  + fsign*(0.06);
      else return 1.;
    }
    else  if (this->Pt()>20. && this->Pt() <= 30.){
      if      ( eta <=  0.800)    return (0.674/0.694)  + fsign*(0.02);
      else if ( eta <=  1.444)    return (0.650/0.661)  + fsign*(0.06);
      else if ( eta <=  1.566)    return (0.391/0.406)  + fsign*(0.06);
      else if ( eta <=  2.000)    return (0.593/0.633)  + fsign*(0.06);
      else if ( eta <=  2.500)    return (0.655/0.664)  + fsign*(0.06);

    }
    else  if (this->Pt()>30. && this->Pt() <= 40.){
      if      ( eta <=  0.800)    return (0.783/0.795)  + fsign*(0.02);
      else if ( eta <=  1.444)    return (0.770/0.780)  + fsign*(0.06);
      else if ( eta <=  1.566)    return (0.541/0.570)  + fsign*(0.06);
      else if ( eta <=  2.000)    return (0.743/0.762)  + fsign*(0.06);
      else if ( eta <=  2.500)    return (0.736/0.764)  + fsign*(0.06);
      else return 1.;
    }
    else  if (this->Pt()>40. && this->Pt() <= 50.){
      if      ( eta <=  0.800)    return (0.839/0.851)  + fsign*(0.02);
      else if ( eta <=  1.444)    return (0.830/0.841)  + fsign*(0.06);
      else if ( eta <=  1.566)    return (0.694/0.707)  + fsign*(0.06);
      else if ( eta <=  2.000)    return (0.828/0.834)  + fsign*(0.06);
      else if ( eta <=  2.500)    return (0.823/0.818)  + fsign*(0.06);
      else return 1.;
    }
    else  if (this->Pt()>50. && this->Pt() <= 200.){
      if      ( eta <=  0.800)    return (0.867/0.877)  + fsign*(0.02);
      else if ( eta <=  1.444)    return (0.855/0.867)  + fsign*(0.06);
      else if ( eta <=  1.566)    return (0.727/0.729)  + fsign*(0.06);
      else if ( eta <=  2.000)    return (0.869/0.864)  + fsign*(0.06);
      else if ( eta <=  2.500)    return (0.851/0.843)  + fsign*(0.06);
      else return 1.;
    }
  }
  else if (name == "cutBasedElectronID-Spring15-25ns-V1-standalone-tight"){
    if (this->Pt()>10. && this->Pt() <= 20.){
      if      ( eta <=  0.800)    return (0.421/0.419)  + fsign*(0.02);
      else if ( eta <=  1.444)    return (0.412/0.383)  + fsign*(0.06);
      else if ( eta <=  1.566)    return (0.279/0.244)  + fsign*(0.06);
      else if ( eta <=  2.000)    return (0.347/0.344)  + fsign*(0.06);
      else if ( eta <=  2.500)    return (0.401/0.379)  + fsign*(0.06);
      else return 1.;
    }
    else  if (this->Pt()>20. && this->Pt() <= 30.){
      if      ( eta <=  0.800)    return (0.547/0.553)  + fsign*(0.02);
      else if ( eta <=  1.444)    return (0.515/0.521)  + fsign*(0.06);
      else if ( eta <=  1.566)    return (0.321/0.322)  + fsign*(0.06);
      else if ( eta <=  2.000)    return (0.487/0.512)  + fsign*(0.06);
      else if ( eta <=  2.500)    return (0.538/0.538)  + fsign*(0.06);

    }
    else  if (this->Pt()>30. && this->Pt() <= 40.){
      if      ( eta <=  0.800)    return (0.655/0.667)  + fsign*(0.02);
      else if ( eta <=  1.444)    return (0.643/0.649)  + fsign*(0.06);
      else if ( eta <=  1.566)    return (0.461/0.474)  + fsign*(0.06);
      else if ( eta <=  2.000)    return (0.625/0.644)  + fsign*(0.06);
      else if ( eta <=  2.500)    return (0.645/0.647)  + fsign*(0.06);
      else return 1.;
    }
    else  if (this->Pt()>40. && this->Pt() <= 50.){
      if      ( eta <=  0.800)    return (0.731/0.745)  + fsign*(0.02);
      else if ( eta <=  1.444)    return (0.724/0.731)  + fsign*(0.06);
      else if ( eta <=  1.566)    return (0.604/0.608)  + fsign*(0.06);
      else if ( eta <=  2.000)    return (0.718/0.723)  + fsign*(0.06);
      else if ( eta <=  2.500)    return (0.719/0.712)  + fsign*(0.06);
      else return 1.;
    }
    else  if (this->Pt()>50. && this->Pt() <= 200.){
      if      ( eta <=  0.800)    return (0.782/0.795)  + fsign*(0.02);
      else if ( eta <=  1.444)    return (0.769/0.784)  + fsign*(0.06);
      else if ( eta <=  1.566)    return (0.649/0.645)  + fsign*(0.06);
      else if ( eta <=  2.000)    return (0.772/0.768)  + fsign*(0.06);
      else if ( eta <=  2.500)    return (0.771/0.753)  + fsign*(0.06);
      else return 1.;
    }
  }


  else {cout << "The ID " << name << " has no SFs" << endl; }
  
  return 1.;
}



///// SETTING CLASS VARIABLES

void KElectron::SetIsFromTau(bool istau){
  k_is_fromtau=istau;
}

void KElectron::SetMotherPdgId(int type){
  k_mother_pdgid=type;
}

void KElectron::SetMCMatchedPdgId(int type){
  k_mc_pdgid=type;
}

void KElectron::SetMotherTruthIndex(int mindex){
  k_mother_index=mindex;
}

void KElectron::SetMCTruthIndex(int tindex){
  k_mc_index=tindex;
}



void KElectron::SetIsChargeFlip(bool iscf){
  k_is_cf=iscf;
}


void KElectron::SetIsPhotonConversion(bool isconv){
  k_is_conv=isconv;
}

void KElectron::SetSCEta(Double_t sceta){
  k_sceta = sceta;
}


void KElectron::Setdz(double d_z){ 
  k_dz = d_z;
}

void KElectron::Setdxy(double d_xy){ 
  k_dxy = d_xy;
}
void KElectron::Setdxy_sig(double d_xysig){
  k_dxy_sig = d_xysig;
}

 
void KElectron::SetPFChargedHadronIso(Double_t cone, Double_t pf_ch){
   if(cone == 0.3)  k_pf_chargedhad_iso03= pf_ch;
   else  if(cone == 0.4)  k_pf_chargedhad_iso04= pf_ch;
}

void KElectron::SetPFPhotonIso(Double_t cone,Double_t pf_ph){
  if(cone == 0.3)   k_pf_photon_iso03 = pf_ph;
  else  if(cone == 0.4) k_pf_photon_iso04 = pf_ph; 
}

void KElectron::SetPFNeutralHadronIso(Double_t cone,Double_t pf_ne){
  if(cone == 0.3)  k_pf_neutral_iso03 = pf_ne;
  else  if(cone == 0.4)k_pf_neutral_iso04 = pf_ne;

}


void KElectron::SetPFRelIso(Double_t cone,Double_t reliso){
  if(cone == 0.3)   k_rel_iso03=reliso;
  else  if(cone == 0.4) k_rel_iso04=reliso;
}

void KElectron::SetPFRelMiniIso(Double_t reliso){
  k_rel_miniiso=reliso;
}

void KElectron::SetPFAbsIso(Double_t cone,Double_t absiso){
  if(cone == 0.3) k_abs_iso03=absiso;
  else  if(cone == 0.4) k_abs_iso04=absiso;
}




void KElectron::SetGsfCtfScPixCharge(bool gsfctfscpix_ch){
  k_gsf_ctscpix_charge = gsfctfscpix_ch;
}


void KElectron::SetHasMatchedConvPhot(Bool_t hasmatchConvPhot){  
  k_hasmatchconvphot = hasmatchConvPhot;
}

void KElectron::SetShiftedEUp(Double_t Eup){
  k_pt_shifted_up= Eup;
}

void KElectron::SetShiftedEDown(Double_t Edown){
  k_pt_shifted_down= Edown;
}

void KElectron::SetSNUID(Int_t snuid){
  snu_id=snuid;
}

void KElectron::SetPassVeto(Bool_t veto){
  pass_veto= veto;
}

void KElectron::SetPassLoose(Bool_t loose){
  pass_loose= loose;
}

void KElectron::SetPassMedium(Bool_t medium){
  pass_medium= medium;
}

void KElectron::SetPassTight(Bool_t tight){
  pass_tight= tight;
}
void KElectron::SetPassHLT(Bool_t hlt){
  pass_hltid= hlt;
}



void KElectron::SetPassHEEP(Bool_t pass){
  pass_heep= pass;
}

void KElectron::SetPassMVATrigMedium(Bool_t pass){
  pass_trigmva_medium= pass;
}

void KElectron::SetPassMVATrigTight(Bool_t pass){
  pass_trigmva_tight= pass;
}

void KElectron::SetPassMVANoTrigMedium(Bool_t pass){
  pass_notrigmva_medium= pass;
}

void KElectron::SetPassMVANoTrigTight(Bool_t pass){
  pass_notrigmva_tight= pass;
}

void KElectron::SetPassMVAZZ(Bool_t pass){
  pass_notrigmva_zz= pass;
}


void KElectron::SetIsPF(Bool_t isPF){
  k_isPF=isPF;
}

void KElectron::SetIsTrigMVAValid(Bool_t b){
  k_istrigmvavalid= b;
}

void KElectron::SetIsMCMatched(Bool_t matched){
  k_mc_matched = matched;
}

void KElectron::SetTrkVx(Double_t trkvx){  
  k_trkvx = trkvx;
}
void KElectron::SetTrkVy(Double_t trkvy){  
  k_trkvy = trkvy;
}
void KElectron::SetTrkVz(Double_t trkvz){  
  k_trkvz = trkvz;
}

void KElectron::SetTrigMatch(TString match){
  k_trig_match = match;
}
