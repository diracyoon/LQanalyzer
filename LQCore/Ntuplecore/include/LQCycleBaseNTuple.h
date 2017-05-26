#ifndef LQCycleBaseNTuple_H
#define LQCycleBaseNTuple_H

class TTree;
class TFile;
class TString;
class TDirectory;
class TBranch;

#include <list>
#include "LQError.h"
#include "LQCycleBaseBase.h"



// forward declaration
class TH1F;


class LQCycleBaseNTuple :  public virtual LQCycleBaseBase{

public:
  LQCycleBaseNTuple();
  virtual ~LQCycleBaseNTuple();

  TFile*                m_outputFile;
  std::vector< TTree* > m_outputTrees;
  std::list< TObject* >   m_inputVarPointers;
  std::list< void* > m_outputVarPointers;  
  

  virtual void CreateOutputTrees(TFile* outputFile, TString name, TString title);
  virtual void FillOutTree();
  virtual void SaveOutputTrees( TDirectory* /*output*/ );
  virtual TTree* GetOutputTree( const char* treeName ) const;
  virtual void MakeOutPutFile(TString path, TString name= "NOTREE");

  
  TFile* GetOutputFile();
  void DeleteInputVariables();
  void CloseFiles();
  void SetDataType(bool type);
  void SetMCPeriod(int period);
  void SetNPStatus( bool type);
  void SetCFStatus( bool type);
  void SetTauStatus( bool type);
  void SetSampleName( TString sample);
  void SetSkimName( TString sample);
  void SetAnalyzerClassName( TString name); 
  void SetDataChannel(TString channel);

  void SetNSampleEvents(double nev);
  void SetOutPutStep(int step);
  void SetNEventsToProcess(int nentries);
  void WriteCycleHists(TH1F* htime, TH1F* mem1, TH1F* mem2);
  
  /// Bool to tell if sample is data or MC in analysis code
  bool k_isdata;

  bool k_running_nonprompt;
  bool k_running_chargeflip;
  bool k_running_taudecays;
  TString k_sample_name;
  TString k_channel;
  TString k_skim;
  TString k_classname;
  TString muscale_fitParametersFile;
  
  // How many entries in the fullsample
  double sample_entries;
  int output_interval;
  int events_to_process;
  int k_mcperiod;
  /// Declare an output variable
  template< class T >
    TBranch* DeclareVariable( T& obj, const char* name,
			       const char* treeName = 0 ) throw( LQError );

 private:
  static const char* RootType( const char* typeid_type ) ;

  ClassDef( LQCycleBaseNTuple, 1 );
  
};

// Don't include the templated function(s) when we're generating
// a dictionary:
#ifndef __CINT__
#include "LQCycleBaseNTuple.icc"
#endif

#endif
