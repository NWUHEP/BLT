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
#include "TH1.h"


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

    // 

    // Params and cuts
    std::unique_ptr<Parameters>         params;
    std::unique_ptr<Cuts>               cuts;
    std::unique_ptr<ParticleSelector>   particleSelector;
    std::unique_ptr<baconhep::TTrigger> trigger;

    // Branches in the output file
    TLorentzVector muonOneP4, muonTwoP4, jetP4, bjetP4, dimuonP4;
    TLorentzVector GenZd, GenMuonLead, GenMuonTrail, GenBQuark, GenQPrime;
    Float_t muonOneIso, muonTwoIso;
    Float_t met, met_phi;
    Float_t jetD0, bjetD0;
    Float_t muonDeltaPhi, muonDeltaEta, muonDeltaR;
    Double_t costheta;
    UInt_t runNumber, evtNumber, lumiSection;
    UInt_t nJets, nBJets;

    //ClassDef(DimuonAnalyzer,0);
};


#endif  // DEMOANALYZER_HH
