// $Id: ExampleAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQGetFakeSF Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *d
 ***************************************************************************/

/// Local includes
#include "GetFakeSF.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (GetFakeSF);


/**
 *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
 *
 */
GetFakeSF::GetFakeSF() :  AnalyzerCore() {


  // To have the correct name in the log:                                                                                                                            
  SetLogName("GetFakeSF");

  Message("In GetFakeSF constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void GetFakeSF::InitialiseAnalysis() throw( LQError ) {
  
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


void GetFakeSF::ExecuteEvents()throw( LQError ){
  


  if(!isData)weight*= MCweight;

  FillHist("NoCut" , 1., MCweight,  0. , 2., 2);
  FillHist("NoCut_w" , 1., weight,  0. , 2., 2);

  std::vector<snu::KElectron> electronLooose_highdxy = GetElectrons(false,true, "ELECTRON_HN_HIGHDXY_FAKELOOSE");
  std::vector<snu::KElectron> electronLooose_lowdxy = GetElectrons(false,true,"ELECTRON_HN_LOWDXY_FAKELOOSE");
  std::vector<snu::KElectron> electronTight_highdxy = GetElectrons(false,true, "ELECTRON_HN_HIGHDXY_TIGHT");
  std::vector<snu::KElectron> electronTight_lowdxy = GetElectrons(false,true, "ELECTRON_HN_LOWDXY_TIGHT");

  

  std::vector<snu::KJet> jetColl           = GetJets("JET_NOLEPTONVETO");
  std::vector<snu::KMuon> muonColl = GetMuons("MUON_HN_LOOSE");  // loose selection                                                                                                                                                                                  

  
  FillCLHist(sighist_ee, "ELECTRON_HN_HIGHDXY_FAKELOOSE", eventbase->GetEvent(), muonColl,electronLooose_highdxy,jetColl, weight);
  FillCLHist(sighist_ee, "ELECTRON_HN_LOWDXY_FAKELOOSE", eventbase->GetEvent(), muonColl,electronLooose_lowdxy,jetColl, weight);
  FillCLHist(sighist_ee, "ELECTRON_HN_HIGHDXY_TIGHT", eventbase->GetEvent(), muonColl,electronTight_highdxy,jetColl, weight);
  FillCLHist(sighist_ee, "ELECTRON_HN_LOWDXY_TIGHT",  eventbase->GetEvent(), muonColl,electronTight_lowdxy,jetColl, weight);

  Float_t ptbins[10] = { 10., 15.,20.,25.,30.,35.,45.,60.,100., 200.};
  Float_t etabins[5] = { 0.,0.8,  1.479, 2.,  2.5};

  for(unsigned int iel =0; iel < electronLooose_lowdxy.size(); iel++){
    FillHist(("LooseEl_lowdxy_pt_eta"), electronLooose_lowdxy[iel].Pt(), fabs(electronLooose_lowdxy[iel].Eta()),  weight, ptbins, 9 , etabins, 4);
    FillHist(("LooseEl_lowdxy_pt"), electronLooose_lowdxy[iel].Pt(),  weight, ptbins, 9 );
    FillHist(("LooseEl_lowdxy_eta"), fabs(electronLooose_lowdxy[iel].Eta()),  weight, etabins, 4);

    if(PassID(electronLooose_lowdxy[iel], "ELECTRON_HN_LOWDXY_TIGHT")){
      FillHist(("TightEl_lowdxy_pt_eta"), electronLooose_lowdxy[iel].Pt(), fabs(electronLooose_lowdxy[iel].Eta()),  weight, ptbins, 9 , etabins, 4);
      FillHist(("TightEl_lowdxy_pt"), electronLooose_lowdxy[iel].Pt(),  weight, ptbins, 9 );
      FillHist(("TightEl_lowdxy_eta"), fabs(electronLooose_lowdxy[iel].Eta()),  weight, etabins, 4);
    }
  }
  
  for(unsigned int iel =0; iel < electronLooose_highdxy.size(); iel++){
    FillHist(("LooseEl_highdxy_pt_eta"), electronLooose_highdxy[iel].Pt(), fabs(electronLooose_highdxy[iel].Eta()),  weight, ptbins, 9 , etabins, 4);
    FillHist(("LooseEl_highdxy_pt"), electronLooose_highdxy[iel].Pt(),  weight, ptbins, 9 );
    FillHist(("LooseEl_highdxy_eta"), fabs(electronLooose_highdxy[iel].Eta()),  weight, etabins, 4);

    if(PassID(electronLooose_highdxy[iel], "ELECTRON_HN_HIGHDXY_TIGHT")){
      FillHist(("TightEl_highdxy_pt_eta"), electronLooose_highdxy[iel].Pt(), fabs(electronLooose_highdxy[iel].Eta()),  weight, ptbins, 9 , etabins, 4);
      FillHist(("TightEl_highdxy_pt"), electronLooose_highdxy[iel].Pt(),  weight, ptbins, 9 );
      FillHist(("TightEl_highdxy_eta"), fabs(electronLooose_highdxy[iel].Eta()),  weight, etabins, 4);
    }
  }


  if(1){
    for(unsigned int iel =0; iel < electronLooose_lowdxy.size(); iel++){
      FillHist(("LooseEl_lowdxy_pt_eta_noweight"), electronLooose_lowdxy[iel].Pt(), fabs(electronLooose_lowdxy[iel].Eta()),  1., ptbins, 9 , etabins, 4);
      FillHist(("LooseEl_lowdxy_pt_noweight"), electronLooose_lowdxy[iel].Pt(),  1., ptbins, 9 );
      FillHist(("LooseEl_lowdxy_eta_noweight"), fabs(electronLooose_lowdxy[iel].Eta()),  1., etabins, 4);

      if(PassID(electronLooose_lowdxy[iel], "ELECTRON_HN_LOWDXY_TIGHT")){
	FillHist(("TightEl_lowdxy_pt_eta_noweight"), electronLooose_lowdxy[iel].Pt(), fabs(electronLooose_lowdxy[iel].Eta()),  1., ptbins, 9 , etabins, 4);
	FillHist(("TightEl_lowdxy_pt_noweight"), electronLooose_lowdxy[iel].Pt(),  1., ptbins, 9 );
	FillHist(("TightEl_lowdxy_eta_noweight"), fabs(electronLooose_lowdxy[iel].Eta()),  1., etabins, 4);
      }
    }

    for(unsigned int iel =0; iel < electronLooose_highdxy.size(); iel++){
      FillHist(("LooseEl_highdxy_pt_eta_noweight"), electronLooose_highdxy[iel].Pt(), fabs(electronLooose_highdxy[iel].Eta()),  1., ptbins, 9 , etabins, 4);
      FillHist(("LooseEl_highdxy_pt_noweight"), electronLooose_highdxy[iel].Pt(),  1., ptbins, 9 );
      FillHist(("LooseEl_highdxy_eta_noweight"), fabs(electronLooose_highdxy[iel].Eta()),  1., etabins, 4);

      if(PassID(electronLooose_highdxy[iel], "ELECTRON_HN_HIGHDXY_TIGHT")){
	FillHist(("TightEl_highdxy_pt_eta_noweight"), electronLooose_highdxy[iel].Pt(), fabs(electronLooose_highdxy[iel].Eta()),  1., ptbins, 9 , etabins, 4);
	FillHist(("TightEl_highdxy_pt_noweight"), electronLooose_highdxy[iel].Pt(),  1., ptbins, 9 );
	FillHist(("TightEl_highdxy_eta_noweight"), fabs(electronLooose_highdxy[iel].Eta()),  1., etabins, 4);
      }
    }

    
  }
  

  return;
}// End of execute event loop
  



void GetFakeSF::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);
  m_logger << DEBUG << "END OF CYCLE: isdata=" << isData <<  LQLogger::endmsg;

}


void GetFakeSF::BeginCycle() throw( LQError ){
  
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

GetFakeSF::~GetFakeSF() {
  
  Message("In GetFakeSF Destructor" , INFO);
  
}



void GetFakeSF::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


///############### THESE ARE FUNCTIONS SPECIFIC TO THIS CYCLE

void GetFakeSF::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this GetFakeSFCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void GetFakeSF::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
}



