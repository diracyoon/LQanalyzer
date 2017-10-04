#ifndef _SKTree_KEvent_H__
#define _SKTree_KEvent_H__

/// Local includes
#include "KParticle.h"

#include <iosfwd>
#include <string>
#include <vector>
#include "TObject.h"

namespace snu {
  
  class KEvent : public TObject{
  public:


    enum json{missing=0,
	      silver=1,
	      gold=2};

    enum met_syst{None=0,
		  MuonEn=1,
		  ElectronEn=2,
		  JetEn=3,
		  JetRes=4,
		  Unclustered=5};

    enum met_type{pfmet=0,
		  puppi=1,
		  nohf=2,
		  pfmva=3};
    
    enum syst_dir{central=0,
		  down=1,
		  up=2};

    KEvent();
    
    ///Copy constructor
    KEvent(const KEvent& el);
    
    ///Destructor    
    ~KEvent() ;

    KEvent& operator= (const KEvent& obj);
    
    /// MET variabkes
    void SetMET(met_type type, double met, double metphi, double sumet); ///defailt MET for analysis
    // MET systematics
    void SetPFMETShift(syst_dir dir, met_syst type, double val);
    void SetPFSumETShift(syst_dir dir, met_syst type, double val);

    void SetTopPtReweight(double tw);
    void SetPFMETx(double metx);
    void SetPFMETy(double mety);

    void SetPFMETType1x(double mett1x);
    void SetPFMETType1y(double mett1y);
    void SetPFMETType1SumEt(double mett1sumet);
    
    void SetPFMETType1Unsmearedx(double mett1);
    void SetPFMETType1Unsmearedy(double mett1);
    void SetPFMETType1xyUnsmearedx(double mett1);
    void SetPFMETType1xyUnsmearedy(double mett1);

    void SetPropagatedRochesterToMET(bool setpr);
    /// PDF
    void SetPDFWeights(std::vector<float> pdfw);
    void SetScaleWeights(std::vector<float> pdfw);
    
    /// Process ID
    void SetLumiSection(int ls);
    /// Vertex
    void SetNVertices(int nvert);
    void SetNGoodVertices(int nvert);
    void SetIsGoodEvent(int isgood);
    void SetVertexInfo(double vX , double vY, double vZ, double ndof);
    
    /// Event
    void SetWeight(double weight);
    void SetGenId(int  id1, int  id2);
    void SetLHEWeight(double  lhe);
    void SetGenX(double  x1, double  x2);
    void SetGenQ(double  Q);

    void SetEventNumber(int ev);
    void SetRunNumber(int rn);
    void SetIsData(bool isdata);
    
    
    /// MET filters
    void SetPassCSCHaloFilterTight(bool pass);
    void SetPassEcalDeadCellTriggerPrimitiveFilter(bool pass);
    void SetPassHBHENoiseFilter(bool pass);
    void SetPassHBHENoiseIsoFilter(bool pass);
    void SetPassBadEESupercrystalFilter(bool pass);
    void SetPassTightHalo2016Filter(bool pass); 
    void SetPassBadChargedCandidateFilter(bool pass);
    void SetPassBadPFMuonFilter(bool pass);
      
    
    /// PileUp reweighting (only in MC)
    void SetPileUpInteractionsTrue(double npu);
    void SetPUWeight(syst_dir sys, double puweight);
    void SetAltPUWeight(syst_dir sys, double puweight);
    void SetPeriodPileupWeight(double puweightB, double puweightC, double puweightD,double puweightE,double puweightF,double puweightG,double puweightH);
    void SetLumiMask(json type, int mask);


    void SetCatVersion(std::string cat);
    
    void SetRho(double rho);

    /// New for CAT v7-4-6 (silver/gold json files)
    Bool_t LumiMask();
    Double_t PileUpWeight(syst_dir dir=central);
    Double_t AltPileUpWeight(syst_dir dir=central);


    Double_t PFMETShifted (met_syst type,syst_dir dir) const;
    Double_t PFSumETShifted(met_syst type,syst_dir dir) const;
    Double_t MET(bool uncorrected=false) const;
    Double_t METPhi( bool uncorrected=false) const;
    Double_t SumET() const;
    
    inline Bool_t PropagatedRochesterToMET() const {return prop_metrc;}
    
    inline Double_t Rho() const {return k_rho;}
    // To make backward compatible
    inline Double_t PFMET(bool uncorrected=false) const {return MET(uncorrected);}
    inline Double_t PFMETPhi(bool uncorrected=false) const {return METPhi(uncorrected);}

    ///k_PF_MET_JetRes_up is same as MET in miniaod
    inline Double_t PFMETUnSmeared() const {return k_PF_MET_JetRes_up;}
    inline Double_t PFMETx() const {return k_PF_METx;}
    inline Double_t PFMETy() const {return k_PF_METy;}

    inline Double_t TopReweight() const{return sqrt(k_top_reweight);}

    inline Double_t PFMETType1() const {return sqrt(pow(k_PF_MET_Type1x,2.) + pow(k_PF_MET_Type1y,2.));}
    inline Double_t PFMETType1x() const {return k_PF_MET_Type1x;}
    inline Double_t PFMETType1y() const {return k_PF_MET_Type1y;}
    inline Double_t PFMETType1Phi() const {return TMath::ATan2(k_PF_MET_Type1y,k_PF_MET_Type1x);}
    inline Double_t PFMETType1SumEt() const {return k_PF_MET_Type1sumet;}

    inline Double_t PFMETUnSmearedType1Phi() const {return TMath::ATan2(k_PF_MET_UnSmearedType1y,k_PF_MET_UnSmearedType1x); }
    inline Double_t PFMETUnSmearedType1Pt() const {return sqrt(pow(k_PF_MET_UnSmearedType1x,2.) + pow(k_PF_MET_UnSmearedType1y,2.));}
    inline Double_t PFMETUnSmearedType1Px() const {return k_PF_MET_UnSmearedType1x;}
    inline Double_t PFMETUnSmearedType1Py() const {return k_PF_MET_UnSmearedType1y;}
    inline Double_t PFMETUnSmearedType1xyPhi() const {return TMath::ATan2(k_PF_MET_UnSmearedType1XYy,k_PF_MET_UnSmearedType1XYx);}
    inline Double_t PFMETUnSmearedType1xyPt() const {return  sqrt(pow(k_PF_MET_UnSmearedType1XYx,2.) + pow(k_PF_MET_UnSmearedType1XYy,2.));}
    inline Double_t PFMETUnSmearedType1xyPx() const {return k_PF_MET_UnSmearedType1XYx;}
    inline Double_t PFMETUnSmearedType1xyPy() const {return k_PF_MET_UnSmearedType1XYy;}


    //// Functions to call in analysis code/selection code

    
    inline std::string CatVersion() const {return k_catversion;}

    
    /// Lumi
    inline Int_t LumiSection()  const {return k_lumisec;}


    
    /// Event specific
    inline Double_t MCWeight() const{return k_mcweight;}
    inline Int_t IsData() const {return k_isData;}
    inline Int_t EventNumber() const {return k_EventNumber;}
    inline Int_t RunNumber() const {return k_RunNumber;}
    
    // Vertex
    inline Int_t nVertices()  const {return k_nvertices;}
    inline Int_t nGoodVertices()  const {return k_ngoodvertices;}

    inline Int_t LumiMaskSilver() const {return k_lumi_mask_silver;}
    inline Int_t LumiMaskGold() const {return k_lumi_mask_gold;}


    inline Bool_t HasGoodPrimaryVertex() const {return k_isgoodevent;}
    inline Double_t VertexX() const {return k_vertexX;}
    inline Double_t VertexY() const {return k_vertexY;}
    inline Double_t VertexZ() const {return k_vertexZ;}
    inline Double_t VertexNDOF() const {return k_vertexNDOF;}

    
    
    /// MET filters
    inline Bool_t   PassCSCHaloFilterTight() const {return   k_passCSCHaloFilterTight;}
    inline Bool_t   PassEcalDeadCellTriggerPrimitiveFilter() const {return   k_passEcalDeadCellTriggerPrimitiveFilter;}
    inline Bool_t   PassHBHENoiseFilter() const {return  k_passHBHENoiseFilter ;}
    inline Bool_t   PassHBHENoiseIsoFilter() const {return  k_passHBHENoiseIsoFilter ;}
    inline Bool_t   PassBadEESupercrystalFilter() const {return  k_passBadEESupercrystalFilter ;}
    inline Bool_t   PassTightHalo2016Filter() const {return  k_passTightHalo2016Filter; }
    inline Bool_t   PassBadChargedCandidateFilter() const {return k_passBadChargedCandFilter;}
    inline Bool_t   PassBadPFMuonFilter() const {return k_passBadpfMuonFilter;}

    //Pileup reweighting
    inline Double_t PileUpInteractionsTrue() const{ return k_PileUpInteractionsTrue;}

    inline Double_t PeriodPileUpWeight(int period) const{
      if (period == 1) return k_pu_gold_weightB;
      if (period == 2) return k_pu_gold_weightC;
      if (period == 3) return k_pu_gold_weightD;
      if (period == 4) return k_pu_gold_weightE;
      if (period == 5) return k_pu_gold_weightF;
      if (period == 6) return k_pu_gold_weightG;
      if (period == 7) return k_pu_gold_weightH;
      return -999.;
    }


    inline Double_t PileUpWeight_Gold(syst_dir sys) const{
      if(sys == central)return k_pu_gold_weight;
      else if(sys == up)return k_pu_gold_p_weight;
      else if(sys == down)return k_pu_gold_m_weight;
      else return -999.;
    }
    
    inline Double_t AltPileUpWeight_Gold(syst_dir sys) const{
      if(sys == central)return k_pu_gold_xs71000_weight;
      else if(sys == up)return k_pu_gold_xs71000_p_weight;
      else if(sys == down)return k_pu_gold_xs71000_m_weight;
      else return -999.;
    }

    inline Double_t LHEWeight() const{return k_lheweight;}
    inline Int_t Id1() const {return k_pdf_id1;}
    inline Int_t Id2() const {return k_pdf_id2;}
    inline Double_t Q() const {return k_pdf_q;}
    inline Double_t x1() const {return k_pdf_x1;}
    inline Double_t x2() const {return k_pdf_x2;}
    
    inline std::vector<Float_t> PdfWeights() const {return k_pdf_weights;}
    inline std::vector<Float_t> ScaleWeights() const {return k_scale_weights;}
    
    virtual void Reset();    
  protected:
    /// Reset function.                                                                  
    
  private:
    /// decalre private functions
    
    Int_t    k_EventNumber, k_RunNumber,k_nvertices,  k_lumisec, k_ngoodvertices,k_pdf_id1, k_pdf_id2, k_lumi_mask_silver, k_lumi_mask_gold;
    std::vector<Float_t> k_pdf_weights, k_scale_weights;
    Double_t k_vertexX,k_vertexY,k_vertexZ, k_vertexNDOF,  k_mcweight, k_lheweight, k_pdf_q, k_pdf_x1, k_pdf_x2;

    Double_t k_PF_MET, k_PF_METphi, k_PF_SumET ;
    Double_t k_PF_METx,k_PF_METy;

    Double_t k_PF_MET_MuonEn_up,k_PF_MET_MuonEn_down; 
    Double_t k_PF_MET_ElectronEn_up,k_PF_MET_ElectronEn_down; 
    Double_t k_PF_MET_JetEn_up,k_PF_MET_JetEn_down, k_PF_SumET_JetEn_up, k_PF_SumET_JetEn_down; 
    Double_t k_PF_MET_JetRes_up,k_PF_MET_JetRes_down, k_PF_SumET_JetRes_up, k_PF_SumET_JetRes_down; 
    Double_t k_PF_MET_unclusteredEn_up,k_PF_MET_unclusteredEn_down, k_PF_SumET_unclusteredEn_up, k_PF_SumET_unclusteredEn_down; 

    Double_t k_PF_MET_Type1, k_PF_MET_Type1x, k_PF_MET_Type1y,  k_PF_MET_Type1sumet;
    Double_t k_PF_MET_UnSmearedType1x, k_PF_MET_UnSmearedType1y, k_PF_MET_UnSmearedType1XYx, k_PF_MET_UnSmearedType1XYy;


    Bool_t k_isData, k_isgoodevent;

    Bool_t k_passBadEESupercrystalFilter,k_passCSCHaloFilterTight,k_passEcalDeadCellTriggerPrimitiveFilter,  k_passHBHENoiseFilter, k_passHBHENoiseIsoFilter, k_passTightHalo2016Filter, k_passBadChargedCandFilter, k_passBadpfMuonFilter;


    Double_t  k_PileUpInteractionsTrue, k_pu_gold_weight, k_pu_gold_p_weight, k_pu_gold_m_weight, k_pu_gold_xs71000_weight, k_pu_gold_xs71000_p_weight, k_pu_gold_xs71000_m_weight;
    Double_t k_pu_gold_weightB, k_pu_gold_weightC, k_pu_gold_weightD,k_pu_gold_weightE,k_pu_gold_weightF,k_pu_gold_weightG,k_pu_gold_weightH;
    Bool_t prop_metrc;

    std::string k_catversion;
    Double_t k_rho;
    Double_t k_top_reweight;

    ClassDef(KEvent,33);
  }; 
  
}//namespace snu

#endif
