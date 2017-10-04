#include "KEvent.h"

#include <iomanip>
#include <iostream>

using namespace std;
using namespace snu;

ClassImp(KEvent)

/**
 *Default constructor.
 */
KEvent::KEvent() :
  TObject(),
  k_EventNumber(0),
  k_RunNumber(0),
  k_nvertices(-999),
  k_lumisec(-999),
  k_ngoodvertices(-999),
  k_pdf_id1(-999),
  k_pdf_id2(-999),
  k_pdf_weights(),
  k_scale_weights(),
  k_vertexX(-999.),
  k_vertexY(-999.),
  k_vertexZ(-999.), 
  k_vertexNDOF(-999.), 
  k_mcweight(-999.),
  k_lheweight(-999.),
  k_pdf_q(-999.),
  k_pdf_x1(-999.),
  k_pdf_x2(-999.),

  k_PF_MET(-999.), 
  k_PF_METx(-999.), 
  k_PF_METy(-999.), 
  k_PF_METphi(-999.),
  k_PF_SumET(-999.), 
  k_PF_MET_MuonEn_up(-999.),
  k_PF_MET_MuonEn_down(-999.),
  k_PF_MET_ElectronEn_up(-999.),
  k_PF_MET_ElectronEn_down(-999.),
  k_PF_MET_JetEn_up(-999.),
  k_PF_MET_JetEn_down(-999.),
  k_PF_SumET_JetEn_up(-999.),
  k_PF_SumET_JetEn_down(-999.),
  k_PF_MET_JetRes_up(-999.),
  k_PF_MET_JetRes_down(-999.),
  k_PF_SumET_JetRes_up(-999.),
  k_PF_SumET_JetRes_down(-999.),
  k_PF_MET_unclusteredEn_up(-999.),
  k_PF_MET_unclusteredEn_down(-999.),
  k_PF_SumET_unclusteredEn_up(-999.),
  k_PF_SumET_unclusteredEn_down(-999.),

  k_PF_MET_Type1x(-999.), 
  k_PF_MET_Type1y(-999.), 
  k_PF_MET_Type1sumet(-999.),
  k_PF_MET_UnSmearedType1x(-999.),
  k_PF_MET_UnSmearedType1y(-999.),
  k_PF_MET_UnSmearedType1XYx(-999.), 
  k_PF_MET_UnSmearedType1XYy(-999.),


  k_isData(false), 
  k_isgoodevent(false),

  k_passBadEESupercrystalFilter(false),
  k_passCSCHaloFilterTight(false),
  k_passEcalDeadCellTriggerPrimitiveFilter(false),
  k_passHBHENoiseFilter(false),
  k_passHBHENoiseIsoFilter(false),
  k_passTightHalo2016Filter(false),
  k_passBadChargedCandFilter(false),
  k_passBadpfMuonFilter(false),
  k_lumi_mask_silver(-999),
  k_lumi_mask_gold(-999),
  k_PileUpInteractionsTrue(-999.),
  k_pu_gold_weight(-999.),
  k_pu_gold_p_weight(-999.),
  k_pu_gold_m_weight(-999.),
  k_pu_gold_xs71000_weight(-999.),
  k_pu_gold_xs71000_p_weight(-999.),
  k_pu_gold_xs71000_m_weight(-999.),

  k_pu_gold_weightB(-999.),
  k_pu_gold_weightC(-999.),
  k_pu_gold_weightD(-999.),
  k_pu_gold_weightE(-999.),
  k_pu_gold_weightF(-999.),
  k_pu_gold_weightG(-999.),
  k_pu_gold_weightH(-999.),
  k_catversion(""),
  k_rho(-999.),
  prop_metrc(false),
  k_top_reweight(1.)
{

}

/**
 * Copy constructor.
 */
KEvent::KEvent(const KEvent& ev) :
  TObject(),
  k_EventNumber(ev.k_EventNumber),
  k_RunNumber(ev.k_RunNumber),
  k_nvertices(ev.k_nvertices),
  k_lumisec(ev.k_lumisec),
  k_ngoodvertices(ev.k_ngoodvertices),
  k_pdf_id1(ev.k_pdf_id1),
  k_pdf_id2(ev.k_pdf_id2),
  k_pdf_weights(ev.k_pdf_weights),
  k_scale_weights(ev.k_scale_weights),

  k_vertexX(ev.k_vertexX),
  k_vertexY(ev.k_vertexY),
  k_vertexZ(ev.k_vertexZ),
  k_vertexNDOF(ev.k_vertexNDOF),
  k_mcweight(ev.k_mcweight),
  k_lheweight(ev.k_lheweight),
  k_pdf_q(ev.k_pdf_q),
  k_pdf_x1(ev.k_pdf_x1),
  k_pdf_x2(ev.k_pdf_x2),

  k_PF_MET(ev.k_PF_MET),
  k_PF_METx(ev.k_PF_METx),
  k_PF_METy(ev.k_PF_METy),
  k_PF_METphi(ev.k_PF_METphi),
  k_PF_SumET(ev.k_PF_SumET),
  k_PF_MET_MuonEn_up(ev.k_PF_MET_MuonEn_up),
  k_PF_MET_MuonEn_down(ev.k_PF_MET_MuonEn_down),
  k_PF_MET_ElectronEn_up(ev.k_PF_MET_ElectronEn_up),
  k_PF_MET_ElectronEn_down(ev.k_PF_MET_ElectronEn_down),
  k_PF_MET_JetEn_up(ev.k_PF_MET_JetEn_up),
  k_PF_MET_JetEn_down(ev.k_PF_MET_JetEn_down),
  k_PF_SumET_JetEn_up(ev.k_PF_SumET_JetEn_up),
  k_PF_SumET_JetEn_down(ev.k_PF_SumET_JetEn_down),
  k_PF_MET_JetRes_up(ev.k_PF_MET_JetRes_up),
  k_PF_MET_JetRes_down(ev.k_PF_MET_JetRes_down),
  k_PF_SumET_JetRes_up(ev.k_PF_SumET_JetRes_up),
  k_PF_SumET_JetRes_down(ev.k_PF_SumET_JetRes_down),
  k_PF_MET_unclusteredEn_up(ev.k_PF_MET_unclusteredEn_up),
  k_PF_MET_unclusteredEn_down(ev.k_PF_MET_unclusteredEn_down),
  k_PF_SumET_unclusteredEn_up(ev.k_PF_SumET_unclusteredEn_up),
  k_PF_SumET_unclusteredEn_down(ev.k_PF_SumET_unclusteredEn_down),

  k_PF_MET_Type1x(ev.k_PF_MET_Type1x),
  k_PF_MET_Type1y(ev.k_PF_MET_Type1y),
  k_PF_MET_Type1sumet(ev.k_PF_MET_Type1sumet),
  k_PF_MET_UnSmearedType1x(ev.k_PF_MET_UnSmearedType1x),
  k_PF_MET_UnSmearedType1y(ev.k_PF_MET_UnSmearedType1y),
  k_PF_MET_UnSmearedType1XYx(ev.k_PF_MET_UnSmearedType1XYx),
  k_PF_MET_UnSmearedType1XYy(ev.k_PF_MET_UnSmearedType1XYy),


  k_isData(ev.k_isData),
  k_isgoodevent(ev.k_isgoodevent),
  k_passBadEESupercrystalFilter(ev.k_passBadEESupercrystalFilter),
  k_passCSCHaloFilterTight(ev.k_passCSCHaloFilterTight),
  k_passEcalDeadCellTriggerPrimitiveFilter(ev.k_passEcalDeadCellTriggerPrimitiveFilter),
  k_passHBHENoiseFilter(ev.k_passHBHENoiseFilter),
  k_passHBHENoiseIsoFilter(ev.k_passHBHENoiseIsoFilter),
  k_passTightHalo2016Filter(ev.k_passTightHalo2016Filter),
  k_passBadChargedCandFilter(ev.k_passBadChargedCandFilter),
  k_passBadpfMuonFilter(ev.k_passBadpfMuonFilter),

  k_lumi_mask_silver(ev.k_lumi_mask_silver),
  k_lumi_mask_gold(ev.k_lumi_mask_gold),
  k_PileUpInteractionsTrue(ev.k_PileUpInteractionsTrue),

  k_pu_gold_weight(ev.k_pu_gold_weight),
  k_pu_gold_p_weight(ev.k_pu_gold_p_weight),
  k_pu_gold_m_weight(ev.k_pu_gold_m_weight),
  k_pu_gold_xs71000_weight(ev.k_pu_gold_xs71000_weight),
  k_pu_gold_xs71000_p_weight(ev.k_pu_gold_xs71000_p_weight),
  k_pu_gold_xs71000_m_weight(ev.k_pu_gold_xs71000_m_weight),

  k_pu_gold_weightB(ev.k_pu_gold_weightB),
  k_pu_gold_weightC(ev.k_pu_gold_weightC),
  k_pu_gold_weightD(ev.k_pu_gold_weightD),
  k_pu_gold_weightE(ev.k_pu_gold_weightE),
  k_pu_gold_weightF(ev.k_pu_gold_weightF),
  k_pu_gold_weightG(ev.k_pu_gold_weightG),
  k_pu_gold_weightH(ev.k_pu_gold_weightH),

  k_catversion(ev.k_catversion),
  prop_metrc(ev.prop_metrc)

{
}


KEvent::~KEvent()
{
}

void KEvent::Reset()
{
  
  k_EventNumber= -999;
  k_RunNumber= -999;
  k_nvertices= -999;
  k_lumisec=-999; 
  k_ngoodvertices= -999;
  k_pdf_id1=-999;
  k_pdf_id2=-999;
  k_vertexX= -999.;
  k_vertexY= -999.;
  k_vertexZ= -999.;
  k_vertexNDOF= -999.;
  k_mcweight= -999.;
  k_lheweight= -999.;
  k_pdf_q=-999.;
  k_pdf_x1= -999.;
  k_pdf_x2=-999.;
  k_pdf_weights.clear();
  
  k_scale_weights.clear();
  k_PF_MET= -999.;
  k_PF_METx= -999.;
  k_PF_METy= -999.;
  k_PF_SumET= -999.;
  k_PF_METphi= -999.;
  k_PF_MET_MuonEn_up = -999.;
  k_PF_MET_MuonEn_down = -999.;
  k_PF_MET_ElectronEn_up = -999.;
  k_PF_MET_ElectronEn_down = -999.;
  k_PF_MET_JetEn_up = -999.;
  k_PF_MET_JetEn_down = -999.;
  k_PF_SumET_JetEn_up = -999.;
  k_PF_SumET_JetEn_down = -999.;
  k_PF_MET_JetRes_up = -999.;
  k_PF_MET_JetRes_down = -999.;
  k_PF_SumET_JetRes_up = -999.;
  k_PF_SumET_JetRes_down = -999.;
  k_PF_MET_unclusteredEn_up = -999.;
  k_PF_MET_unclusteredEn_down = -999.;
  k_PF_SumET_unclusteredEn_up = -999.;
  k_PF_SumET_unclusteredEn_down = -999.;
  
  k_PF_MET_Type1x = -999.;
  k_PF_MET_Type1y = -999.;
  k_PF_MET_Type1sumet = -999.;
  k_PF_MET_UnSmearedType1x = -999.;
  k_PF_MET_UnSmearedType1y = -999.;
  k_PF_MET_UnSmearedType1XYx = -999.;
  k_PF_MET_UnSmearedType1XYy = -999.;



  k_isData= false;
  k_isgoodevent = false;
  k_passBadEESupercrystalFilter= false;
  k_passCSCHaloFilterTight= false;
  k_passEcalDeadCellTriggerPrimitiveFilter= false;
  k_passHBHENoiseFilter= false;
  k_passHBHENoiseIsoFilter= false;
  k_passTightHalo2016Filter= false;
  k_passBadChargedCandFilter=false;
  k_passBadpfMuonFilter=false;
  
  k_lumi_mask_silver=-999;
  k_lumi_mask_gold=-999;

  k_PileUpInteractionsTrue = -999.;
  k_pu_gold_weight = -999.;
  k_pu_gold_p_weight=-999.;
  k_pu_gold_m_weight = -999.;
  k_pu_gold_xs71000_weight = -999.;
  k_pu_gold_xs71000_p_weight=-999.;
  k_pu_gold_xs71000_m_weight = -999.;

  k_pu_gold_weightB = -999;
  k_pu_gold_weightC = -999;
  k_pu_gold_weightD = -999;
  k_pu_gold_weightE = -999;
  k_pu_gold_weightF = -999;
  k_pu_gold_weightG = -999;
  k_pu_gold_weightH = -999;

  k_catversion="";
  k_rho= -999.;

  prop_metrc=false;
  k_top_reweight=1.;
}



KEvent& KEvent::operator= (const KEvent& p)
{
    if (this != &p) {

      k_EventNumber=p.EventNumber();
      k_RunNumber=p.RunNumber();
      k_nvertices= p.nVertices();
      k_lumisec= p.LumiSection();
      k_ngoodvertices = p.nGoodVertices();
      k_pdf_id1 = p.Id1();
      k_pdf_id2 = p.Id2();
      k_pdf_weights = p.PdfWeights();
      k_scale_weights = p.ScaleWeights();
      
      k_vertexX=p.VertexX();
      k_vertexY=p.VertexY();
      k_vertexZ=p.VertexZ();
      k_vertexNDOF=p.VertexNDOF();
      k_mcweight = p.MCWeight();
      k_lheweight = p.LHEWeight();
      k_pdf_q = p.Q();
      k_pdf_x1 = p.x1();
      k_pdf_x2 = p.x2();
           

      
      k_PF_MET= p.MET();
      k_PF_METx= p.PFMETx();
      k_PF_METy= p.PFMETy();
      k_PF_METphi= p.METPhi();
      k_PF_SumET = p.SumET();
      k_PF_MET_MuonEn_up = p.PFMETShifted(MuonEn,up);
      k_PF_MET_MuonEn_down = p.PFMETShifted(MuonEn,down);
      k_PF_MET_ElectronEn_up = p.PFMETShifted(ElectronEn,up);
      k_PF_MET_ElectronEn_down = p.PFMETShifted(ElectronEn,down);
      k_PF_MET_JetEn_up = p.PFMETShifted(JetEn,up);
      k_PF_MET_JetEn_down = p.PFMETShifted(JetEn,down);
      k_PF_SumET_JetEn_up = p.PFSumETShifted(JetEn,up);
      k_PF_SumET_JetEn_down = p.PFSumETShifted(JetEn,down);
      k_PF_MET_JetRes_up = p.PFMETShifted(JetRes,up);
      k_PF_MET_JetRes_down = p.PFMETShifted(JetRes,down);
      k_PF_SumET_JetRes_up = p.PFSumETShifted(JetRes,up);
      k_PF_SumET_JetRes_down = p.PFSumETShifted(JetRes,down);
      k_PF_MET_unclusteredEn_up = p.PFMETShifted(Unclustered,up);
      k_PF_MET_unclusteredEn_down = p.PFMETShifted(Unclustered,down);
      k_PF_SumET_unclusteredEn_up = p.PFSumETShifted(Unclustered,up);
      k_PF_SumET_unclusteredEn_down = p.PFSumETShifted(Unclustered,down);
      
      k_PF_MET_Type1x = p.PFMETType1x();
      k_PF_MET_Type1y = p.PFMETType1y();
      k_PF_MET_Type1sumet = p.PFMETType1SumEt();
      k_PF_MET_UnSmearedType1x = p.PFMETUnSmearedType1Px();
      k_PF_MET_UnSmearedType1y = p.PFMETUnSmearedType1Py();
      k_PF_MET_UnSmearedType1XYx = p.PFMETUnSmearedType1xyPx();
      k_PF_MET_UnSmearedType1XYy = p.PFMETUnSmearedType1xyPy();

      k_isData=p.IsData();
      k_isgoodevent= p.HasGoodPrimaryVertex();

      k_passBadEESupercrystalFilter = p.PassBadEESupercrystalFilter();
      k_passCSCHaloFilterTight = p.PassCSCHaloFilterTight();
      k_passEcalDeadCellTriggerPrimitiveFilter = p.PassEcalDeadCellTriggerPrimitiveFilter();
      k_passHBHENoiseFilter = p.PassHBHENoiseFilter();
      k_passHBHENoiseIsoFilter= p.PassHBHENoiseIsoFilter();
      k_passTightHalo2016Filter= p.PassTightHalo2016Filter();
      k_passBadChargedCandFilter= p.PassBadChargedCandidateFilter();
      k_passBadpfMuonFilter=p.PassBadPFMuonFilter();


      k_PileUpInteractionsTrue = p.PileUpInteractionsTrue();
      k_lumi_mask_silver = p.LumiMaskSilver();
      k_lumi_mask_gold = p.LumiMaskGold();
      k_pu_gold_xs71000_weight = p.AltPileUpWeight_Gold(central);
      k_pu_gold_xs71000_p_weight = p.AltPileUpWeight_Gold(up);
      k_pu_gold_xs71000_m_weight = p.AltPileUpWeight_Gold(down);


      k_pu_gold_weight = p.PileUpWeight_Gold(central);
      k_pu_gold_p_weight= p.PileUpWeight_Gold(up);
      k_pu_gold_m_weight= p.PileUpWeight_Gold(down);

      k_pu_gold_weightB = p.PeriodPileUpWeight(1);
      k_pu_gold_weightC = p.PeriodPileUpWeight(2);
      k_pu_gold_weightD = p.PeriodPileUpWeight(3);
      k_pu_gold_weightE = p.PeriodPileUpWeight(4);
      k_pu_gold_weightF = p.PeriodPileUpWeight(5);
      k_pu_gold_weightG = p.PeriodPileUpWeight(6);
      k_pu_gold_weightH = p.PeriodPileUpWeight(7);

      k_catversion = p.CatVersion();
      k_rho= p.Rho();
	
      prop_metrc = p.PropagatedRochesterToMET();
      k_top_reweight=p.TopReweight();
    }
    
    return *this;
}

//// SET CLASS VARIBALES
void KEvent::SetCatVersion(std::string cat){
  k_catversion = cat;
}

void KEvent::SetRho(double rho){
  k_rho=rho;
}
void KEvent::SetTopPtReweight(double topw){
  k_top_reweight = topw;
}

void KEvent::SetPropagatedRochesterToMET(bool hasprop){
  prop_metrc = hasprop;
}

void KEvent::SetPassCSCHaloFilterTight(bool pass){
  k_passCSCHaloFilterTight = pass;
}

void KEvent::SetPassBadEESupercrystalFilter(bool pass){
  k_passBadEESupercrystalFilter = pass;
}

void KEvent::SetPassEcalDeadCellTriggerPrimitiveFilter(bool pass){
  k_passEcalDeadCellTriggerPrimitiveFilter = pass;
}

void KEvent::SetPassHBHENoiseFilter(bool pass){
  k_passHBHENoiseFilter= pass;
}
void KEvent::SetPassHBHENoiseIsoFilter(bool pass){
  k_passHBHENoiseIsoFilter= pass;
}

void KEvent::SetPassTightHalo2016Filter(bool pass){
  k_passTightHalo2016Filter= pass;
}

void KEvent::SetPassBadChargedCandidateFilter(bool pass){
  
  k_passBadChargedCandFilter= pass;
}

void KEvent::SetPassBadPFMuonFilter(bool pass){
  k_passBadpfMuonFilter=pass;
}

void KEvent::SetPileUpInteractionsTrue(double npu){
  k_PileUpInteractionsTrue = npu;
}

void KEvent::SetPUWeight( syst_dir sys, double puw){
  if(sys==central)  k_pu_gold_weight = puw;
  if(sys==up)  k_pu_gold_p_weight = puw;
  if(sys==down)  k_pu_gold_m_weight = puw;
  
}
void KEvent::SetPeriodPileupWeight( double puwB,double puwC,double puwD,double puwE,double puwF,double puwG,double puwH){

  k_pu_gold_weightB = puwB;
  k_pu_gold_weightC = puwC;
  k_pu_gold_weightD = puwD;
  k_pu_gold_weightE = puwE;
  k_pu_gold_weightF = puwF;
  k_pu_gold_weightG = puwG;
  k_pu_gold_weightH = puwH;
}

void KEvent::SetAltPUWeight(syst_dir sys, double puw){
  if(sys==central)  k_pu_gold_xs71000_weight = puw;
  if(sys==up)  k_pu_gold_xs71000_p_weight = puw;
  if(sys==down)  k_pu_gold_xs71000_m_weight = puw;
  
}

void KEvent::SetLumiMask(json type, int mask){
  if(type==silver) k_lumi_mask_silver=mask;
  else if(type==gold) k_lumi_mask_gold=mask;
  else {std::cout<< "LumiMask not set correctly" << std::endl; exit(1);}
}


Bool_t KEvent::LumiMask(){
  
  if(k_lumi_mask_gold==1) return true;
  else return false;

}




void KEvent::SetWeight(double mcweight){
  k_mcweight = mcweight;
}




void KEvent::SetEventNumber(int ev){
  k_EventNumber = ev;
}

void KEvent::SetRunNumber(int rn){
  k_RunNumber = rn;
}

void KEvent::SetIsData(bool isdata){
  k_isData = isdata;
}



//// Setting Missing ET variables

void KEvent::SetMET(met_type type, double met, double metphi, double sumet){
  if(type == pfmet){
    k_PF_MET = met;
    k_PF_METphi = metphi;
    k_PF_SumET =sumet;
  }


  else {cout << "Problem setting MET" << endl; exit(0) ;}
}

void KEvent::SetPFMETx( double metx){
  k_PF_METx= metx;
}

void KEvent::SetPFMETy(double mety){
  k_PF_METy= mety;
}




void KEvent::SetPFMETShift(syst_dir dir, met_syst type, double val){
  if(dir == up){
    if(type==MuonEn)  k_PF_MET_MuonEn_up = val;
    if(type==ElectronEn)  k_PF_MET_ElectronEn_up= val;
    if(type==JetEn)  k_PF_MET_JetEn_up= val;
    if(type==JetRes)  k_PF_MET_JetRes_up= val;
    if(type==Unclustered)  k_PF_MET_unclusteredEn_up= val;
  }
  else if (dir == down){
    if(type==MuonEn)  k_PF_MET_MuonEn_down= val;
    if(type==ElectronEn)  k_PF_MET_ElectronEn_down= val;
    if(type==JetEn)  k_PF_MET_JetEn_down= val;
    if(type==JetRes)  k_PF_MET_JetRes_down= val;
    if(type==Unclustered)  k_PF_MET_unclusteredEn_down= val;
  }
  else {cout << "Problem setting MET shifted values" << endl ; exit(0);}
}

void KEvent::SetPFSumETShift(syst_dir dir, met_syst type, double val){
  if(dir == up){
    if(type==JetEn)  k_PF_SumET_JetEn_up= val;
    if(type==JetRes)  k_PF_SumET_JetRes_up= val;
    if(type==Unclustered)  k_PF_SumET_unclusteredEn_up= val;
  }
  else if (dir == down){
    if(type==JetEn)  k_PF_SumET_JetEn_down= val;
    if(type==JetRes)  k_PF_SumET_JetRes_down= val;
    if(type==Unclustered)  k_PF_SumET_unclusteredEn_down= val;
  }
 
  else {cout << "Problem setting SumET shifted values" <<endl ; exit(0);}

}


void KEvent::SetPFMETType1x(double mett1x){
  k_PF_MET_Type1x = mett1x;
}

void KEvent::SetPFMETType1y(double mett1y){
  k_PF_MET_Type1y = mett1y;
}


void KEvent::SetPFMETType1SumEt(double mett1sumet){
  k_PF_MET_Type1sumet = mett1sumet;
}

void KEvent::SetPFMETType1Unsmearedx(double mett1){
  k_PF_MET_UnSmearedType1x = mett1;
}

void KEvent::SetPFMETType1Unsmearedy(double mett1){
  k_PF_MET_UnSmearedType1y = mett1;
}

void KEvent::SetPFMETType1xyUnsmearedx(double mett1){
  k_PF_MET_UnSmearedType1XYx = mett1;
}
void KEvent::SetPFMETType1xyUnsmearedy(double mett1){
  k_PF_MET_UnSmearedType1XYy= mett1;
}


/// Settizang Vertex Variables
void KEvent::SetNVertices(int nvert){
  k_nvertices = nvert;
}


void KEvent::SetLumiSection(int ls){
  k_lumisec = ls;
}



void KEvent::SetIsGoodEvent(int nvert){
  if(nvert > 0) k_isgoodevent = true;
  else k_isgoodevent = false;
}


Double_t KEvent::PileUpWeight(syst_dir sys){
  if(TString(k_catversion).Contains("v7-4-")) return -999;
  return PileUpWeight_Gold(sys);

}

Double_t KEvent::AltPileUpWeight(syst_dir sys){
  return AltPileUpWeight_Gold(sys);
}

///New forCAT v7-4-5 (MET systematics in one function)
Double_t KEvent::PFMETShifted ( met_syst type, syst_dir dir) const{
  
  if(TString(k_catversion).Contains("v7-4-4")) return k_PF_MET;
  if(dir == up){
    if(type==None) return k_PF_MET;
    else if(type==MuonEn) return k_PF_MET_MuonEn_up;
    else if(type==ElectronEn) return k_PF_MET_ElectronEn_up;
    else if(type==JetEn) return k_PF_MET_JetEn_up;
    else if(type==JetRes) return k_PF_MET_JetRes_up;
    else if(type==Unclustered) return k_PF_MET_unclusteredEn_up;
    else return -999.;
  }
  else  if(dir == down){
    if(type==None) return k_PF_MET;
    else if(type==MuonEn) return k_PF_MET_MuonEn_down;
    else if(type==ElectronEn) return k_PF_MET_ElectronEn_down;
    else if(type==JetEn) return k_PF_MET_JetEn_down;
    else if(type==JetRes) return k_PF_MET_JetRes_down;
    else if(type==Unclustered) return k_PF_MET_unclusteredEn_down;
    else return-999.;

  }
  return k_PF_MET;
}

Double_t KEvent::PFSumETShifted(met_syst type,syst_dir dir) const{
  if(TString(k_catversion).Contains("v7-4-4")) return k_PF_SumET;
 
  if(dir == up){
    if(type==None) return k_PF_SumET;
    if(type==JetEn) return k_PF_SumET_JetEn_up;
    if(type==JetRes) return k_PF_SumET_JetRes_up;
    if(type==Unclustered) return k_PF_SumET_unclusteredEn_up;
  }
  else if(dir == down){
    if(type==None) return k_PF_SumET;
    if(type==JetEn) return k_PF_SumET_JetEn_down;
    if(type==JetRes) return k_PF_SumET_JetRes_down;
    if(type==Unclustered) return k_PF_SumET_unclusteredEn_down;
  }
  return k_PF_SumET;
}


Double_t KEvent::MET(bool uncorrected) const{
  
  if(uncorrected) return {sqrt(pow(k_PF_METx,2.)+ pow(k_PF_METy,2.))};
  else  return k_PF_MET;
  
  return -999.;
}
 
Double_t KEvent::METPhi(bool uncorrected) const{
  
  if(uncorrected) return {TMath::ATan2(k_PF_METy,k_PF_METx)};
  else  return k_PF_METphi;

  return -999.;
}

Double_t KEvent::SumET() const{
  
 
  return k_PF_SumET;
  return -999.;
}



/// Settizang Vertex Variables
void KEvent::SetNGoodVertices(int nvert){
  k_ngoodvertices = nvert;
}


void KEvent::SetVertexInfo(double vX, double vY, double vZ, double NDOF){
  k_vertexX = vX;
  k_vertexY = vY;
  k_vertexZ = vZ;
  k_vertexNDOF = NDOF;
}

void KEvent::SetGenId(int id1, int id2){
  k_pdf_id1 = id1;
  k_pdf_id2 = id2;
}

void KEvent::SetPDFWeights(std::vector<float> pdfweights){
  k_pdf_weights = pdfweights;
}

void KEvent::SetScaleWeights(std::vector<float> scaleweights){
  k_scale_weights= scaleweights;
}



void KEvent::SetGenQ(double q ){
  k_pdf_q = q;
}


void KEvent::SetGenX(double x1, double x2){
  k_pdf_x1 = x1;
  k_pdf_x2 = x2;
}

void KEvent::SetLHEWeight(double lw){
  k_lheweight = lw;
}





