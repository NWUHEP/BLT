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




#ifndef TauSelector_HH
#define TauSelector_HH

// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"
#include "BLT/BLTAnalysis/interface/Cuts.hh"
#include "BLT/BLTAnalysis/interface/TriggerSelector.hh"
#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"
#include "BLT/BLTAnalysis/interface/WeightUtils.h"

#include "BLT/BLTAnalysis/interface/RoccoR.h"

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


class TauSelector: public BLTSelector {
public:
    TauSelector();
    ~TauSelector();

    void   Begin(TTree *tree);
    Bool_t Process(Long64_t entry);
    void   Terminate();
    void   ReportPostBegin();
    void   ReportPostTerminate();

    TFile *outFile;
    std::map<string, TTree*> outTrees;
    std::map<string, TH1D*> eventCounts;

    // Lumi mask
    RunLumiRangeMap lumiMask;
    
    // rochester muon corrections
    RoccoR *muonCorr;
    TRandom3 *rng;


    // Params and cuts
    std::unique_ptr<Parameters>         params;
    std::unique_ptr<Cuts>               cuts;

    // Utilities
    std::unique_ptr<ParticleSelector>   particleSelector;
    std::unique_ptr<baconhep::TTrigger> trigger;
    std::unique_ptr<WeightUtils>        weights;

    std::vector<string> triggerNames;

    // Branches in the output file
    // event data
    UInt_t runNumber, lumiSection;
    ULong64_t evtNumber;
    UInt_t nBJets,nJets,nMuons,nElectrons,nTaus,nMistaus;

    // weights and uncertainties
    Float_t eventWeight;

    // physics object Lorentz vectors
    // TLorentzVector tauOneP4, tauTwoP4;


    Float_t tauOnePt, tauTwoPt, tauOneEta, tauTwoEta;

    // Additional lepton data
    Int_t tauOneFlavor, tauTwoFlavor, tauOneIso, tauTwoIso, tauOneDecayMode, tauTwoDecayMode;
    Float_t tauOneIsoMVA,tauOnePuppiChHadIso,tauOnePuppiGammaIso,tauOnePuppiNeuHadIso;
    Float_t tauTwoIsoMVA,tauTwoPuppiChHadIso,tauTwoPuppiGammaIso,tauTwoPuppiNeuHadIso;


    float GetMuonIsolation(const baconhep::TMuon*);
    float GetElectronIsolation(const baconhep::TElectron*, float);
};


#endif  // TauSelector_HH
