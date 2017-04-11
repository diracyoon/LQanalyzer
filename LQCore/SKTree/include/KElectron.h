#ifndef _SKTree_KElectron_H__
#define _SKTree_KElectron_H__

/// Local includes
#include "KParticle.h"

#include <iosfwd>
#include <string>
#include "TLorentzVector.h"

namespace snu {
  
  class KElectron : public KParticle {
  public:


    enum ElectronType{PROMPT=0,
		      FAKE=1,
		      PHOTONFAKE=2,
		      CONV_CF=3,
		      CONV_NONECF=4,
		      CF=5,
    };
    enum ElectronMotherType{none=0,
			    Z=1,
			    W=2,
			    ZorW=3,
			    pion=4,
    };

    KElectron();
    
    ///Copy constructor
    KElectron(const KElectron& el);
    
    ///Destructor    
    virtual ~KElectron() ;

    KElectron& operator= (const KElectron& obj);
    
    
    float ScaleFactor(const std::string& name, int sign) const ;


    // set kinematic variables
    void SetSCEta(Double_t sceta);
    
    //##### NOTE charge/pt/eta/phi use tlv class
    
    /// Trigger matching
    
    
    //// set   vertex variables
    void SetType(int eltype);
    void Setdxy(Double_t d_xy);
    void Setdxy_sig(Double_t d_xysig);
    void Setdz(Double_t d_z);
    
    void SetSNUID(int id);
    void SetPassVeto(Bool_t pass);
    void SetPassLoose(Bool_t pass);
    void SetPassMedium(Bool_t pass);
    void SetPassTight(Bool_t pass);
    void SetPassHEEP(Bool_t pass);

    

    void SetPassMVATrigMedium(Bool_t pass);
    void SetPassMVATrigTight(Bool_t pass);
    void SetPassMVANoTrigMedium(Bool_t pass);
    void SetPassMVANoTrigTight(Bool_t pass);
    
    void SetIsPF(Bool_t ispf);
    void SetIsChargeFlip(Bool_t iscf);
    void SetIsPhotonConversion(Bool_t isconv);
    void SetIsFromTau(Bool_t istau);
    void SetIsMCMatched(Bool_t ismatch);
    
    void SetMCMatchedPdgId(Int_t pg);
    void SetMotherPdgId(Int_t pg);
    void SetMotherTruthIndex(Int_t mindex);
    void SetMCTruthIndex(Int_t t_index);
    /// set ISO variables
    void SetPFChargedHadronIso(Double_t cone,Double_t pf_ch);
    void SetPFPhotonIso(Double_t cone,Double_t pf_ph);
    void SetPFNeutralHadronIso(Double_t cone,Double_t pf_ne);

    
    void SetPFRelIso(Double_t cone, Double_t pf_rel);
    void SetPFRelMiniIso( Double_t pf_rel);
    void SetPFAbsIso(Double_t cone, Double_t pf_abs);

    
    // set charge variables
    
    void SetGsfCtfScPixCharge(bool gsfctfscpix_ch);
    
    /// set conversion variables
    void SetHasMatchedConvPhot(Bool_t hasmatchConvPhot);

    
    void SetShiftedEUp(Double_t Eup);
    void SetShiftedEDown(Double_t Edown);

    void SetTrkVx(Double_t trkvx);
    void SetTrkVy(Double_t trkvy);
    void SetTrkVz(Double_t trkvz);

    void SetTrigMatch(TString match);
    void SetIsTrigMVAValid(bool b);

    bool TriggerMatched(TString path);

    ///// Functions to call class variables
    
    inline Int_t GetType()  const {
      if(k_eltype < 32) return k_eltype;
      else {
	if(fabs(k_mc_pdgid) == 211) return 32;
	if(fabs(k_mc_pdgid) == 310) return 33;
	if(fabs(k_mc_pdgid) == 431) return 34;
	if(fabs(k_mc_pdgid) == 13) return 35;
	return 36;
      }
    }

    inline Bool_t  IsEBFiducial() {return bool (fabs(SCEta()) < 1.442);}
    inline Bool_t  IsEEFiducial() {return bool (fabs(SCEta()) > 1.560 && fabs(SCEta()) < 2.50);}
      
    /// // Kinematic variables
    inline Double_t  SCEta() const {return k_sceta;}
    
    
    inline Double_t PtShiftedUp() const{ return k_pt_shifted_up;}
    inline Double_t PtShiftedDown() const{ return k_pt_shifted_down;}

    /// Trigger matching
    
    // ID variables
    inline Bool_t PassVeto() const{return pass_veto;}
    inline Bool_t PassLoose() const{return pass_loose;}
    inline Bool_t PassMedium() const{return pass_medium;}
    inline Bool_t PassTight() const{return pass_tight;}

    // HEEP ID 6.0
    inline Bool_t PassHEEP() const{return pass_heep;}

    // MVA ID
    inline Bool_t PassTrigMVAMedium() const{return pass_trigmva_medium;}
    inline Bool_t PassTrigMVATight() const{return pass_trigmva_tight;}
    inline Bool_t PassNotrigMVAMedium() const{return pass_notrigmva_medium;}
    inline Bool_t PassNotrigMVATight() const{return pass_notrigmva_tight;}
    
    
    inline Bool_t MCMatched() const{return k_mc_matched;}
    inline Bool_t IsPF() const{return k_isPF;}
    inline Bool_t MCIsPrompt() const{return k_mc_matched;}
    inline Bool_t MCIsCF() const{return k_is_cf;}
    inline Bool_t MCIsFromConversion() const{return k_is_conv;}
    inline Bool_t MCFromTau() const{return k_is_fromtau;}
    inline Int_t MCMatchedPdgId() const{return k_mc_pdgid;}
    inline Int_t MotherPdgId() const{return k_mother_pdgid;}
    inline Int_t MotherTruthIndex() const{return k_mother_index;}
    inline Int_t MCTruthIndex() const{return k_mc_index;}
    
    inline Int_t SNUID() const{return snu_id;}

    inline KElectron::ElectronType GetParticleType() const{ 
      if(k_is_conv&&k_is_cf) return KElectron::CONV_CF;
      if(k_is_conv&&!k_is_cf)   return KElectron::CONV_NONECF; 
      if(k_is_cf)  return KElectron::CF;
      if(k_mc_matched) return KElectron::PROMPT;
      if(k_mc_pdgid==22) return KElectron::PHOTONFAKE;
      return KElectron::FAKE;

    }

    inline KElectron::ElectronMotherType GetMotherType() const{
      if(k_mother_pdgid == 23) return KElectron::Z;
      if(fabs(k_mother_pdgid) == 24) return KElectron::W;
      if(k_mother_pdgid == -99999 ) return KElectron::ZorW;
      return  KElectron::pion;

    }
    // charge variables
    
    inline Bool_t GsfCtfScPixChargeConsistency()  const {return k_gsf_ctscpix_charge;}
    
    // Conversion variables
    inline Bool_t PassesConvVeto() const {return k_hasmatchconvphot;}
    

    inline Bool_t IsTrigMVAValid() const{return k_istrigmvavalid;}
    
    // Isolation Variables
    inline Double_t PFChargedHadronIso(double cone) const {
      if(cone == 0.3)   return k_pf_chargedhad_iso03;
      else  if(cone == 0.4)  return k_pf_chargedhad_iso04;
      else return -999.;
    }
    
    inline Double_t PFPhotonIso(double cone) const {
      if(cone == 0.3)   return k_pf_photon_iso03;
      else  if(cone == 0.4)  return k_pf_photon_iso04;
      else return -999.;
    }
    
    inline Double_t PFNeutralHadronIso(double cone) const {
      if(cone == 0.3)   return k_pf_neutral_iso03;
      else  if(cone == 0.4)  return k_pf_neutral_iso04;
      else return -999.;
    }

		
    inline Double_t PFRelIso(double cone) const {
      if(cone == 0.3)   return k_rel_iso03;
      else  if(cone == 0.4)   return k_rel_iso04;
      else return -999.;
    }

    inline Double_t PFRelMiniIso() const { return k_rel_miniiso; }


    inline Double_t PFAbsIso(double cone) const {
      if(cone == 0.3)   return k_abs_iso03;
      else  if(cone == 0.4)   return k_abs_iso04;
      else return -999.;
    }

    /// VtxDist with vertex chosen to be primary   
    inline Double_t  dxy() const {return  k_dxy;}
    inline Double_t  dxySig() const {return  k_dxy_sig;}
    inline Double_t  dz() const {return  k_dz;}
    

    inline Double_t  TrkVx() const {return  k_trkvx;}
    inline Double_t  TrkVy() const {return  k_trkvy;}
    inline Double_t  TrkVz() const {return  k_trkvz;}

    inline TString TrigMatch() const{return k_trig_match;}

  protected:
    /// Reset function.                                                                  
    virtual void Reset();    
    
  private:
    /// decalre private functions

    Double_t k_pf_chargedhad_iso03, k_pf_photon_iso03, k_pf_neutral_iso03, k_pf_chargedhad_iso04, k_pf_photon_iso04, k_pf_neutral_iso04, k_rel_iso03, k_rel_iso04, k_rel_miniiso;
    Double_t k_abs_iso03, k_abs_iso04;
    Double_t k_dxy,k_dxy_sig, k_dz,k_trkvx,  k_trkvy,  k_trkvz;
    Double_t k_sceta;
    
    Bool_t k_gsf_ctscpix_charge,pass_tight, pass_veto, pass_medium, pass_loose, k_mc_matched,  k_is_cf,k_is_conv, k_is_fromtau,k_isPF,k_hasmatchconvphot, pass_heep, pass_trigmva_medium, pass_trigmva_tight, pass_notrigmva_medium, pass_notrigmva_tight, k_istrigmvavalid ;
    
    Double_t k_pt_shifted_up, k_pt_shifted_down;
    Int_t snu_id,k_mother_pdgid, k_mc_pdgid,k_mother_index, k_mc_index;
    TString k_trig_match;
    Int_t k_eltype;

    ClassDef(KElectron,25);
  }; 
  
}//namespace snu

#endif
