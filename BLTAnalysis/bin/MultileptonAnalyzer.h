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


#ifndef MULTILEPTONANALYZER_HH
#define MULTILEPTONANALYZER_HH

// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"
#include "BLT/BLTAnalysis/interface/Cuts.hh"
#include "BLT/BLTAnalysis/interface/TriggerSelector.hh"
#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"
#include "BLT/BLTAnalysis/interface/WeightUtils.h"

#include "BLT/BLTAnalysis/interface/muresolution.h"
#include "BLT/BLTAnalysis/interface/rochcor2012wasym.h"

// BaconAna class definitions (might need to add more)
#include "BaconAna/Utils/interface/TTrigger.hh"

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


class MultileptonAnalyzer: public BLTSelector {
public:
    MultileptonAnalyzer();
    ~MultileptonAnalyzer();

    void   Begin(TTree *tree);
    Bool_t Process(Long64_t entry);
    void   Terminate();
    void   ReportPostBegin();
    void   ReportPostTerminate();

    TFile *outFile;
    TTree *outTree;

    // Lumi mask
    RunLumiRangeMap lumiMask;

    // rochester muon corrections
    rochcor2012 *muonCorr;

    // Params and cuts
    std::unique_ptr<Parameters>         params;
    std::unique_ptr<Cuts>               cuts;

    // Utilities
    std::unique_ptr<ParticleSelector>   particleSelector;
    std::unique_ptr<baconhep::TTrigger> trigger;
    std::unique_ptr<WeightUtils>        weights;

    std::vector<string> triggerNames;

    // Branches in the output file
    UInt_t runNumber, lumiSection, nPU, nPartons;
    ULong64_t evtNumber;
    Bool_t triggerStatus;
    Float_t eventWeight;
    UInt_t nJets, nFwdJets, nBJets, nMuons, nElectrons;

    TLorentzVector leptonOneP4, leptonTwoP4, jetP4, bjetP4, genBJetP4;
    Float_t leptonOneIso, leptonTwoIso;
    Int_t leptonOneQ, leptonTwoQ;
    Int_t leptonOneFlavor, leptonTwoFlavor;
    Bool_t leptonOneTrigger, leptonTwoTrigger;

    Float_t jetD0, bjetD0;
    Float_t bjetTag, jetTag, genBJetTag;
    Float_t bjetPUID, jetPUID;
    Float_t bjetFlavor, jetFlavor;
    Float_t met, metPhi;

    //ClassDef(MultileptonAnalyzer,0);
};


#endif  // MULTILEPTONANALYZER_HH
