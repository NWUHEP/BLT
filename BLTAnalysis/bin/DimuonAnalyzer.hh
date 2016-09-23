// =============================================================================
// A simple analysis on Bacon ntuples
//
// Input arguments:
//   argv[1] => input bacon file name
//   argv[2] => number of entries
//   argv[3] => ...
//
// Users should inherit from BLTSelector and implement the three functions:
//   Begin()
//   Process()
//   Terminate()
// =============================================================================


#ifndef DEMOANALYZER_HH
#define DEMOANALYZER_HH

// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"
#include "BLT/BLTAnalysis/interface/Cuts.hh"
#include "BLT/BLTAnalysis/interface/TriggerSelector.hh"
#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"

#include "BLT/BLTAnalysis/interface/muresolution.h"
#include "BLT/BLTAnalysis/interface/rochcor2012wasym.h"

// ROOT headers
#include <TLorentzVector.h>
#include <TVector3.h>

// C++ headers
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <regex>


class DimuonAnalyzer: public BLTSelector {
public:
    DimuonAnalyzer();
    ~DimuonAnalyzer();

    void    Begin(TTree *tree);
    Bool_t  Process(Long64_t entry);
    void    Terminate();

    void    ReportPostBegin();
    void    ReportPostTerminate();

    TFile       *outFile;
    TTree       *outTree;

    // Lumi mask
    RunLumiRangeMap lumiMask;

    // rochester muon corrections
    rochcor2012 *muonCorr;

    // Params and cuts
    std::unique_ptr<Parameters>         params;
    std::unique_ptr<Cuts>               cuts;
    std::unique_ptr<ParticleSelector>   particleSelector;
    std::unique_ptr<baconhep::TTrigger> trigger;

    // Branches in the output file
    TLorentzVector muonOneP4, muonTwoP4, jetP4, bjetP4;
    Float_t muonOneIso, muonTwoIso;
    Int_t muonOneQ, muonTwoQ;
    Float_t jetD0, bjetD0;
    Float_t bjetTag;
    UInt_t nJets, nFwdJets, nBJets;
    Float_t met, metPhi;
    UInt_t runNumber, lumiSection;
    ULong64_t evtNumber;
    Bool_t triggerStatus;
    Float_t eventWeight;

    //ClassDef(DimuonAnalyzer,0);
};


#endif  // DEMOANALYZER_HH
