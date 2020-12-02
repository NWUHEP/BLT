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


#ifndef HLTPRINTOUTANALYZER_HH
#define HLTPRINTOUTANALYZER_HH

// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"
#include "BLT/BLTAnalysis/interface/TriggerSelector.hh"

// BaconAna class definitions (might need to add more)
#include "BaconAna/Utils/interface/TTrigger.hh"

// ROOT headers
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRandom3.h>

// C++ headers
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <iostream>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <regex>
#define MZ 91.1876
#define WZ 2.4952

class hltPrintoutAnalyzer: public BLTSelector {
public:
    hltPrintoutAnalyzer();
    ~hltPrintoutAnalyzer();

    void   Begin(TTree *tree);
    Bool_t Process(Long64_t entry);
    void   Terminate();
    void   ReportPostBegin();
    void   ReportPostTerminate();

    TFile *outFile;
    //std::map<string, TTree*> outTrees;
    //std::map<string, TH1D*> eventCounts;
    TTree *outTree;

    // Lumi mask
    RunLumiRangeMap lumiMask;

    // Params 
    std::unique_ptr<Parameters>         params;

    // Utilities
    std::unique_ptr<baconhep::TTrigger> trigger;
    std::vector<std::string> triggerNames;

    // Branches in the output file
    
    // event data
    UInt_t runNumber, lumiSection;
    ULong64_t evtNumber;

    //ClassDef(hltPrintoutAnalyzer,0);
};

#endif  // HLTPRINTOUTANALYZER_HH
