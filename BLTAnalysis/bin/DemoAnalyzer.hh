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
#include <memory>


class DemoAnalyzer: public BLTSelector {
public:
    DemoAnalyzer();
    ~DemoAnalyzer();

    void    Begin(TTree *tree);
    Bool_t  Process(Long64_t entry);
    void    Terminate();

private:
};


#endif  // DEMOANALYZER_HH
