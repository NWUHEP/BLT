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


#ifndef DIMUONANALYZER_HH
#define DIMUONANALYZER_HH

// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"
#include "BLT/BLTAnalysis/interface/Cuts.hh"
#include "BLT/BLTAnalysis/interface/TriggerSelector.hh"
#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"

#include "BLT/BLTAnalysis/interface/RoccoR.h"
#include "BLT/BLTAnalysis/interface/rochcor2016.h"

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
        std::string  outFileName;
        std::string  outTreeName;

        // Params and cuts
        std::unique_ptr<Parameters>         params;
        std::unique_ptr<Cuts>               cuts;
        std::unique_ptr<ParticleSelector>   particleSelector;
        std::unique_ptr<baconhep::TTrigger> trigger;

        // Lumi mask
        RunLumiRangeMap lumiMask;

        // rochester muon corrections
        rochcor2016 *muonCorr;

        // Branches in the output file
        TLorentzVector muonOneP4, muonTwoP4, jetP4, bjetP4;
        Float_t muonOneIso, muonTwoIso;
        Float_t met, met_phi;
        Float_t jetD0, bjetD0;
        UInt_t runNumber, lumiSection;
        ULong64_t evtNumber;
        Bool_t triggerStatus;
        UInt_t nJets, nFwdJets, nBJets;

        //ClassDef(DimuonAnalyzer,0);
};


#endif  // DIMUONANALYZER_HH
