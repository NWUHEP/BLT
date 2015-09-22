// =============================================================================
// A simple analysis on Bacon ntuples
//
// Input arguments:
//   argv[1] => input bacon file name
//   argv[2] => output bacon bits file name
//
// =============================================================================


#ifndef DEMOANALYZER_HH
#define DEMOANALYZER_HH


// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/AnalysisParameters.hh"
#include "BLT/BLTAnalysis/interface/AnalysisCuts.hh"
#include "BLT/BLTAnalysis/interface/TriggerSelector.hh"
#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"

// ROOT headers
#include <TLorentzVector.h>
#include <TMath.h>

// C++ headers
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <memory>


class DemoAnalyzer: public BLTSelector {
public:
    DemoAnalyzer();
    ~DemoAnalyzer();

    void    Begin(TTree *tree);
    Bool_t  Process(Long64_t entry);
    void    Terminate();

    TFile       *outFile;
    TTree       *outTree;
    std::string  outFileName;
    std::string  outTreeName;

    // Params and cuts
    std::unique_ptr<AnalysisParameters> params;
    std::unique_ptr<AnalysisCuts>       cuts;
    std::unique_ptr<TriggerSelector>    triggerSelector;
    std::unique_ptr<ParticleSelector>   particleSelector;

    // Branches in the output file
    TLorentzVector muonOne;
    TLorentzVector muonTwo;
    TLorentzVector dimuon;

    TLorentzVector genMuonOne;
    TLorentzVector genMuonTwo;
    TLorentzVector genZ;

private:
};


#endif  // DEMOANALYZER_HH
