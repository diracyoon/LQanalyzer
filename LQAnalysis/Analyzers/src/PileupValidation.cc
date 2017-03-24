// $Id: ExampleAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQPileupValidation Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *d
 ***************************************************************************/

/// Local includes
#include "PileupValidation.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (PileupValidation);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
PileupValidation::PileupValidation() :  AnalyzerCore() {

  mapcounter.clear();

  functionality = PileupValidation::VALIDATION;


  // To have the correct name in the log:                                                                                                                            
  SetLogName("PileupValidation");

  Message("In PileupValidation constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void PileupValidation::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //

   Message("Making clever hists for Z ->ll test code", INFO);

   //// Initialise Plotting class functions
   /// MakeCleverHistograms ( type, "label")  type can be muhist/elhist/jethist/sighist


   return;
}


void PileupValidation::ExecuteEvents()throw( LQError ){
  
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
  
  if(!isData) weight*= MCweight;
  
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return;
  
  FillHist("MET/PFMet_uncleaned" , eventbase->GetEvent().PFMET(), MCweight,  0. , 300., 300);
  if(!PassMETFilter()) return;     /// Initial event cuts :                                                                                                                                                                                                                 

  FillHist("MET/PFMet_t1xy" , eventbase->GetEvent().PFMET(), weight,  0. , 300., 300);
  FillHist("MET/PFMetPhi_t1xy" , eventbase->GetEvent().METPhi(snu::KEvent::pfmet), weight,  -3. , 3., 100);
  FillHist("MET/PFMet_t1" , eventbase->GetEvent().PFMETType1(), weight,  0. , 300., 300);
  FillHist("MET/PFMetPhi_t1" , eventbase->GetEvent().PFMETType1Phi(), weight,  -3. , 3., 100);
  FillHist("MET/PFMet_unsmeared" , eventbase->GetEvent().PFMETUnSmearedType1Pt(), weight,  0. , 300., 300);
  FillHist("MET/PFMetPhi_unsmeared" , eventbase->GetEvent().PFMETUnSmearedType1Phi(), weight,  -3. , 3., 100);
  std::vector<snu::KJet> jets =  GetJets("JET_HN", 30., 2.5);
  std::vector<snu::KJet> jets_pu =  GetJets("JET_HN_PU", 30., 2.5);
  if(isData){
    // met
    if(GetDataPeriod() == 1)   FillHist("METPerPeriod/PFMet_t1xy_B" , eventbase->GetEvent().PFMET(), weight,  0. , 300., 300);
    if(GetDataPeriod() == 2)   FillHist("METPerPeriod/PFMet_t1xy_C" , eventbase->GetEvent().PFMET(), weight,  0. , 300., 300);
    if(GetDataPeriod() == 3)   FillHist("METPerPeriod/PFMet_t1xy_D" , eventbase->GetEvent().PFMET(), weight,  0. , 300., 300);
    if(GetDataPeriod() == 4)   FillHist("METPerPeriod/PFMet_t1xy_E" , eventbase->GetEvent().PFMET(), weight,  0. , 300., 300);
    if(GetDataPeriod() == 5)   FillHist("METPerPeriod/PFMet_t1xy_F" , eventbase->GetEvent().PFMET(), weight,  0. , 300., 300);
    if(GetDataPeriod() == 6)   FillHist("METPerPeriod/PFMet_t1xy_G" , eventbase->GetEvent().PFMET(), weight,  0. , 300., 300);
    if(GetDataPeriod() == 7)   FillHist("METPerPeriod/PFMet_t1xy_H" , eventbase->GetEvent().PFMET(), weight,  0. , 300., 300);
    /// nvtx
    if(GetDataPeriod() == 1)   FillHist("NvtxPerPeriod/PFMet_t1xy_B" , eventbase->GetEvent().nVertices() ,weight, 0. , 75., 75);
    if(GetDataPeriod() == 2)   FillHist("NvtxPerPeriod/PFMet_t1xy_C" , eventbase->GetEvent().nVertices() ,weight, 0. , 75., 75);
    if(GetDataPeriod() == 3)   FillHist("NvtxPerPeriod/PFMet_t1xy_D" , eventbase->GetEvent().nVertices() ,weight, 0. , 75., 75);
    if(GetDataPeriod() == 4)   FillHist("NvtxPerPeriod/PFMet_t1xy_E" , eventbase->GetEvent().nVertices() ,weight, 0. , 75., 75);
    if(GetDataPeriod() == 5)   FillHist("NvtxPerPeriod/PFMet_t1xy_F" , eventbase->GetEvent().nVertices() ,weight, 0. , 75., 75);
    if(GetDataPeriod() == 6)   FillHist("NvtxPerPeriod/PFMet_t1xy_G" , eventbase->GetEvent().nVertices() ,weight, 0. , 75., 75);
    if(GetDataPeriod() == 7)   FillHist("NvtxPerPeriod/PFMet_t1xy_H" , eventbase->GetEvent().nVertices() ,weight, 0. , 75., 75);
    // njets
    if(GetDataPeriod() == 1)   FillHist("NJetPerPeriod/PFMet_t1xy_B" , jets.size(), weight, 0., 10., 10);
    if(GetDataPeriod() == 2)   FillHist("NJetPerPeriod/PFMet_t1xy_C" , jets.size(), weight, 0., 10., 10);
    if(GetDataPeriod() == 3)   FillHist("NJetPerPeriod/PFMet_t1xy_D" , jets.size(), weight, 0., 10., 10);
    if(GetDataPeriod() == 4)   FillHist("NJetPerPeriod/PFMet_t1xy_E" , jets.size(), weight, 0., 10., 10);
    if(GetDataPeriod() == 5)   FillHist("NJetPerPeriod/PFMet_t1xy_F" , jets.size(), weight, 0., 10., 10);
    if(GetDataPeriod() == 6)   FillHist("NJetPerPeriod/PFMet_t1xy_G" , jets.size(), weight, 0., 10., 10);
    if(GetDataPeriod() == 7)   FillHist("NJetPerPeriod/PFMet_t1xy_H" , jets.size(), weight, 0., 10., 10);
    // pu cut jets
    if(GetDataPeriod() == 1)   FillHist("NJetPUPerPeriod/PFMet_t1xy_B" , jets_pu.size(), weight, 0., 10., 10);
    if(GetDataPeriod() == 2)   FillHist("NJetPUPerPeriod/PFMet_t1xy_C" , jets_pu.size(), weight, 0., 10., 10);
    if(GetDataPeriod() == 3)   FillHist("NJetPUPerPeriod/PFMet_t1xy_D" , jets_pu.size(), weight, 0., 10., 10);
    if(GetDataPeriod() == 4)   FillHist("NJetPUPerPeriod/PFMet_t1xy_E" , jets_pu.size(), weight, 0., 10., 10);
    if(GetDataPeriod() == 5)   FillHist("NJetPUPerPeriod/PFMet_t1xy_F" , jets_pu.size(), weight, 0., 10., 10);
    if(GetDataPeriod() == 6)   FillHist("NJetPUPerPeriod/PFMet_t1xy_G" , jets_pu.size(), weight, 0., 10., 10);
    if(GetDataPeriod() == 7)   FillHist("NJetPUPerPeriod/PFMet_t1xy_H" , jets_pu.size(), weight, 0., 10., 10);


  }


  if(false){
    TString analysis_trigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
    std::vector<TString> triggerslist;
    triggerslist.push_back(analysis_trigger);


    if(!isData ||  (isData&&PassTrigger(analysis_trigger))  ){

      std::vector<snu::KElectron> electronTightColl=GetElectrons(TString("ELECTRON_POG_TIGHT"),15., 2.5);
      
      if(!isData) weight*=  mcdata_correction->GetDoubleEGTriggerEff(electronTightColl);
      if(OppositeCharge(electronTightColl)){
	if(electronTightColl.at(0).Pt() > 26. && electronTightColl.at(1).Pt() > 22.){
	  if(GetDiLepMass(electronTightColl) < 120. && GetDiLepMass(electronTightColl)  > 60. ){
	    if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 75., 75);
	    else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 75., 75);
	    if(isData) FillHist("Nint_nocut_data", eventbase->GetEvent().PileUpInteractionsTrue(),weight, 0. , 75., 75);
	    if(isData) FillHist("Nint_nocut_mc", eventbase->GetEvent().PileUpInteractionsTrue(),weight, 0. , 75., 75);
	  }
	}
      }
    }
  }
  else{
    TString analysis_trigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
    std::vector<TString> triggerslist;
    triggerslist.push_back(analysis_trigger);

    if(!isData ||  (isData&&PassTrigger(analysis_trigger))  ){
      std::vector<snu::KElectron> electronTightColl=GetElectrons(TString("ELECTRON_POG_TIGHT"),15., 2.5);
      
      if(!isData) weight*=  mcdata_correction->GetDoubleEGTriggerEff(electronTightColl);
      
      if(OppositeCharge(electronTightColl)){
        if(electronTightColl.at(0).Pt() > 26. && electronTightColl.at(1).Pt() > 22.){
	  
	  float mc_corr(1.);
	  float weight2(1.);
	  float weight2d(1.);
	  float weight3(1.);
	  float weight4(1.);
	  float weight5(1.);
	  
	  float id_iso_sf(1.);
	  float trigger_ps(1.);
	  float reco_weight=1.;
	  float puweight(1.);

	  
	  if(!isData){
	    id_iso_sf=   mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronTightColl,0);
	    trigger_ps= WeightByTrigger(analysis_trigger, TargetLumi)  ;
	    weight*= trigger_ps;
	    reco_weight = mcdata_correction->ElectronRecoScaleFactor(electronTightColl);
	    /// standard CMS pu weight
	    weight2 = weight * eventbase->GetEvent().PileUpWeight();
	    weight2d = weight * eventbase->GetEvent().PileUpWeight(snu::KEvent::down);
	    /// reweight nvtx to DY
	    weight3 = weight * mcdata_correction->UserPileupWeight(eventbase->GetEvent());
	    /// standard CMS pu weight but different if -p X is set in submittion                                                                                                                                                                                                                                       
	    weight4 = weight * mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());
	    mc_corr =  id_iso_sf * reco_weight * puweight;

	    /// need to ad dZ on POG https://indico.cern.ch/event/604947/contributions/2467335/attachments/1409577/2158348/EleIDSFs_IPCuts_Update.pdf
	  }
	  
	  weight5 = weight3 * mc_corr;
	  weight2 = weight2 *  mc_corr;
	  weight2d = weight2d *  mc_corr;
	  
	  /// Z peak
	  if(GetDiLepMass(electronTightColl) < 120. && GetDiLepMass(electronTightColl)  > 60. ){
	    /// out of box MET
	    FillHist("ZMET_PFMet_t1xy" , eventbase->GetEvent().PFMET(), weight,  0. , 200., 200);
	    FillHist("ZMET_PFMet_t1" , eventbase->GetEvent().PFMETType1(), weight,  0. , 200., 200);
	    FillHist("ZMET_PFMet_t1_unsmeared" , eventbase->GetEvent().PFMETUnSmearedType1Pt(), weight,  0. , 200., 200);
	    FillHist("ZMET_NJet", jets.size(), weight, 0., 10., 10);
	    FillHist("ZMET_Nvtx",  eventbase->GetEvent().nVertices(), weight, 0., 75., 75);
	    FillHist("ZMET_Nint",  eventbase->GetEvent().PileUpInteractionsTrue(), weight, 0., 75., 75);

	    //_ weighted for 69.2 mb
	    FillHist("ZMETNintPU69p2_PFMet_t1xy" , eventbase->GetEvent().PFMET(), weight2,  0. , 200., 200);
            FillHist("ZMETNintPU69p2_PFMet_t1" , eventbase->GetEvent().PFMETType1(), weight2,  0. , 200., 200);
            FillHist("ZMETNintPU69p2_PFMet_t1_unsmeared" , eventbase->GetEvent().PFMETUnSmearedType1Pt(), weight2,  0. , 200., 200);
	    FillHist("ZMETNintPU69p2_NJet", jets.size(), weight2, 0., 10., 10);
	    FillHist("ZMETNintPU69p2_Nvtx",  eventbase->GetEvent().nVertices(), weight2, 0., 75., 75);
            FillHist("ZMETNintPU69p2_Nint",  eventbase->GetEvent().PileUpInteractionsTrue(), weight2, 0., 75., 75);

	    FillHist("ZMETNintPU69p2d_PFMet_t1xy" , eventbase->GetEvent().PFMET(), weight2d,  0. , 200., 200);
            FillHist("ZMETNintPU69p2d_PFMet_t1" , eventbase->GetEvent().PFMETType1(), weight2d,  0. , 200., 200);
            FillHist("ZMETNintPU69p2d_PFMet_t1_unsmeared" , eventbase->GetEvent().PFMETUnSmearedType1Pt(), weight2d,  0. , 200., 200);
            FillHist("ZMETNintPU69p2d_NJet", jets.size(), weight2d, 0., 10., 10);
            FillHist("ZMETNintPU69p2d_Nvtx",  eventbase->GetEvent().nVertices(), weight2d, 0., 75., 75);
            FillHist("ZMETNintPU69p2d_Nint",  eventbase->GetEvent().PileUpInteractionsTrue(), weight2d, 0., 75., 75);


	    /// reweight nvtx to DY
	    FillHist("ZMETNvtx_PFMet_t1xy" , eventbase->GetEvent().PFMET(), weight3,  0. , 200., 200);
            FillHist("ZMETNvtx_PFMet_t1" , eventbase->GetEvent().PFMETType1(), weight3,  0. , 200., 200);
            FillHist("ZMETNvtx_PFMet_t1_unsmeared" , eventbase->GetEvent().PFMETUnSmearedType1Pt(), weight3,  0. , 200., 200);
	    FillHist("ZMETNvtx_NJet", jets.size(), weight3, 0., 10., 10);
	    FillHist("ZMETNvtx_Nvtx",  eventbase->GetEvent().nVertices(), weight3, 0., 75., 75);
            FillHist("ZMETNvtx_Nint",  eventbase->GetEvent().PileUpInteractionsTrue(), weight3, 0., 75., 75);

	    /// per period
	    FillHist("ZMETNintPU69p2PP_PFMet_t1xy" , eventbase->GetEvent().PFMET(), weight4,  0. , 200., 200);
            FillHist("ZMETNintPU69p2PP_PFMet_t1" , eventbase->GetEvent().PFMETType1(), weight4,  0. , 200., 200);
            FillHist("ZMETNintPU69p2PP_PFMet_t1_unsmeared" , eventbase->GetEvent().PFMETUnSmearedType1Pt(), weight4,  0. , 200., 200);
	    FillHist("ZMETNintPU69p2PP_NJet", jets.size(), weight4, 0., 10., 10);
	    FillHist("ZMETNintPU69p2PP_Nvtx",  eventbase->GetEvent().nVertices(), weight4, 0., 75., 75);
            FillHist("ZMETNintPU69p2PP_Nint",  eventbase->GetEvent().PileUpInteractionsTrue(), weight4, 0., 75., 75);
	    
	    /// reweight nvtx to DY                                                                                                                                                                                                                                            
            FillHist("ZMETNvtxCorr_PFMet_t1xy" , eventbase->GetEvent().PFMET(), weight5,  0. , 200., 200);
            FillHist("ZMETNvtxCorr_PFMet_t1" , eventbase->GetEvent().PFMETType1(), weight5,  0. , 200., 200);
            FillHist("ZMETNvtxCorr_PFMet_t1_unsmeared" , eventbase->GetEvent().PFMETUnSmearedType1Pt(), weight5,  0. , 200., 200);
            FillHist("ZMETNvtxCorr_NJet", jets.size(), weight5, 0., 10., 10);
	    FillHist("ZMETNvtxCorr_Nvtx",  eventbase->GetEvent().nVertices(), weight5, 0., 75., 75);
            FillHist("ZMETNvtxCorr_Nint",  eventbase->GetEvent().PileUpInteractionsTrue(), weight5, 0., 75., 75);


	    
	    if(jets.size() == 0) {
	      FillHist("ZMETNvtxCorr_nj0_PFMet_t1xy" , eventbase->GetEvent().PFMET(), weight5,  0. , 200., 200);
	      FillHist("ZMETNvtxCorr_nj0_PFMet_t1" , eventbase->GetEvent().PFMETType1(), weight5,  0. , 200., 200);
	      FillHist("ZMETNvtxCorr_nj0_PFMet_t1_unsmeared" , eventbase->GetEvent().PFMETUnSmearedType1Pt(), weight5,  0. , 200., 200);
	      FillHist("ZMETNvtxCorr_nj0_NJet", jets.size(), weight5, 0., 10., 10);
	      FillHist("ZMETNvtxCorr_nj0_Nvtx",  eventbase->GetEvent().nVertices(), weight5, 0., 75., 75);
	      FillHist("ZMETNvtxCorr_nj0_Nint",  eventbase->GetEvent().PileUpInteractionsTrue(), weight5, 0., 75., 75);

	    }
	    else{
	      FillHist("ZMETNvtxCorr_nj1_PFMet_t1xy" , eventbase->GetEvent().PFMET(), weight5,  0. , 200., 200);
              FillHist("ZMETNvtxCorr_nj1_PFMet_t1" , eventbase->GetEvent().PFMETType1(), weight5,  0. , 200., 200);
              FillHist("ZMETNvtxCorr_nj1_PFMet_t1_unsmeared" , eventbase->GetEvent().PFMETUnSmearedType1Pt(), weight5,  0. , 200., 200);
              FillHist("ZMETNvtxCorr_nj1_NJet", jets.size(), weight5, 0., 10., 10);
	      FillHist("ZMETNvtxCorr_nj1_Nvtx",  eventbase->GetEvent().nVertices(), weight5, 0., 75., 75);
	      FillHist("ZMETNvtxCorr_nj1_Nint",  eventbase->GetEvent().PileUpInteractionsTrue(), weight5, 0., 75., 75);

	    }
	  }// Zpeak
	  if(NBJet(jets) > 1 && jets.size() > 2){

	    FillHist("TopMETNintPU69p2_PFMet_t1xy" , eventbase->GetEvent().PFMET(), weight2,  0. , 200., 200);
            FillHist("TopMETNintPU69p2_PFMet_t1" , eventbase->GetEvent().PFMETType1(), weight2,  0. , 200., 200);
            FillHist("TopMETNintPU69p2_PFMet_t1_unsmeared" , eventbase->GetEvent().PFMETUnSmearedType1Pt(), weight2,  0. , 200., 200);
            FillHist("TopMETNintPU69p2_NJet", jets.size(), weight2, 0., 10., 10);
            FillHist("TopMETNintPU69p2_Nvtx",  eventbase->GetEvent().nVertices(), weight2, 0., 75., 75);
            FillHist("TopMETNintPU69p2_Nint",  eventbase->GetEvent().PileUpInteractionsTrue(), weight2, 0., 75., 75);


	    FillHist("TopMETNvtxCorr_PFMet_t1xy" , eventbase->GetEvent().PFMET(), weight5,  0. , 200., 200);
            FillHist("TopMETNvtxCorr_PFMet_t1" , eventbase->GetEvent().PFMETType1(), weight5,  0. , 200., 200);
            FillHist("TopMETNvtxCorr_PFMet_t1_unsmeared" , eventbase->GetEvent().PFMETUnSmearedType1Pt(), weight5,  0. , 200., 200);
            FillHist("TopMETNvtxCorr_NJet", jets.size(), weight5, 0., 10., 10);
            FillHist("TopMETNvtxCorr_Nvtx",  eventbase->GetEvent().nVertices(), weight5, 0., 75., 75);
            FillHist("TopMETNvtxCorr_Nint",  eventbase->GetEvent().PileUpInteractionsTrue(), weight5, 0., 75., 75);

	  }
	  
	}// pt cuts
      } // l+l-
    }// trigger
  }/// not pileup weights mode
  
  return;
}



void PileupValidation::counter(TString cut, float w){
  w=1.;
  map<TString,float>::iterator itmapcounter = mapcounter.find(cut) ;
  if (itmapcounter == mapcounter.end()){
    mapcounter[cut] = w;
  }
  else{
    float sum = itmapcounter->second;
    mapcounter[cut] = sum+w;
  }

}


void PileupValidation::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger << DEBUG << "END OF CYCLE: isdata=" << isData <<  LQLogger::endmsg;

  for ( map<TString,float>::iterator itmapcounter = mapcounter.begin(); itmapcounter != mapcounter.end(); itmapcounter++){
    cout << itmapcounter->first << " has count "<< itmapcounter->second << endl;
  }

}


void PileupValidation::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  

  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //  DeclareVariable(out_muons, "Signal_Muons");

  
  return;
  
}

PileupValidation::~PileupValidation() {
  
  Message("In PileupValidation Destructor" , INFO);
  
}



void PileupValidation::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void PileupValidation::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this PileupValidationCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void PileupValidation::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //

}



