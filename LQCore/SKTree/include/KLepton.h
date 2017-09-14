#ifndef KLepton_h
#define KLepton_h

/// c++ std libraries
#include <iostream>
#include <vector>
#include <map>
using namespace std;

// ROOT include(s): 

/// SNU classes
#include "KParticle.h"
#include "KMuon.h"
#include "KElectron.h"

class KLepton: public snu::KParticle {

public:

  enum Flavour{NOTSET=0, MUON=1, ELECTRON=2};

  KLepton();
  KLepton(const KLepton& lep);
  KLepton(const snu::KMuon& muon);
  KLepton(const snu::KElectron& electron);
  ~KLepton();

  inline Flavour LeptonFlavour() const {return k_flavour;}
  inline Double_t dZ() const {return k_dz;}
  inline Double_t dXY() const {return k_dxy;}
  inline Double_t dXYSig() const {return k_dxy_sig;}
  inline Double_t RelIso() const {return k_reliso;}
  inline Int_t GetType() const {return k_leptype;}
  inline bool MCIsCF() const {return k_mciscf;}
  inline const snu::KMuon* GetMuonPtr() const {return MuonPtr;}
  inline const snu::KElectron* GetElectronPtr() const {return ElectronPtr;}
    
  //==== FIXME
  //template <typename T>
  //  T GetLeptonPtr();

  KLepton& operator= (const KLepton& obj);
  KLepton& operator= (const snu::KMuon& obj);
  KLepton& operator= (const snu::KElectron& obj);

private:

  Flavour k_flavour;
  Double_t k_dz, k_dxy, k_dxy_sig, k_reliso;
  Int_t k_leptype;
  bool k_mciscf;
  const snu::KMuon* MuonPtr;
  const snu::KElectron* ElectronPtr;

};

#endif
