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


#ifndef ZEEAnalyzer_HH
#define ZEEAnalyzer_HH


// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"
#include "BLT/BLTAnalysis/interface/Cuts.hh"
#include "BLT/BLTAnalysis/interface/TriggerSelector.hh"
#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"
#include "BLT/BLTAnalysis/interface/WeightUtils.h"
#include "BLT/BLTAnalysis/interface/ElectronCorrector.h"

#include "BLT/BLTAnalysis/interface/RoccoR.h"

// BaconAna class definitions (might need to add more)
#include "BaconAna/Utils/interface/TTrigger.hh"
#include "BaconAna/DataFormats/interface/TLHEWeight.hh"

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


class ZEEAnalyzer: public BLTSelector {
public:
    ZEEAnalyzer();
    ~ZEEAnalyzer();

    void   Begin(TTree *tree);
    Bool_t Process(Long64_t entry);
    void   Terminate();
    void   ReportPostBegin();
    void   ReportPostTerminate();

    TFile *outFile;
    std::map<string, TTree*> outTrees;
    std::map<string, TH1D*> eventCounts;
    TH1D *hGenCat;

    // Lumi mask
    RunLumiRangeMap lumiMask;

    // rochester muon corrections
    RoccoR *muonCorr;
    TRandom3 *rng;

    // electron scale and smear corrector (trash)
    EnergyScaleCorrection *electronScaler; 

    // Params and cuts
    std::unique_ptr<Parameters>         params;
    std::unique_ptr<Cuts>               cuts;

    // Utilities
    std::unique_ptr<ParticleSelector>   particleSelector;
    std::unique_ptr<baconhep::TTrigger> trigger;
    std::unique_ptr<WeightUtils>        weights;

    std::vector<string> triggerNames;


    ////////////////////////////////
    // Branches in the output file
    // event data
    UInt_t runNumber, lumiSection, nPV;
    ULong64_t evtNumber;
    Float_t nPU;

    // generator level data
    UInt_t nPartons;
    Float_t genWeight;

    // weights and uncertainties
    Float_t eventWeight;

    // lepton variable

    Float_t leptonOnePt,  leptonTwoPt;
    Float_t leptonOneEta, leptonTwoEta;
    Float_t leptonOnePhi, leptonTwoPhi;
    Float_t leptonOneIso, leptonTwoIso;
    Int_t leptonOneFlavor, leptonTwoFlavor;
    Int_t leptonOnePassTriggerTest, leptonTwoPassTriggerTest;
    Int_t leptonOneNotInGap, leptonTwoNotInGap;

    Float_t dileptonMass, dileptonPt, dileptonEta, dileptonPhi;

    

    // met
    Float_t met, metPhi;

    // object counters
    UInt_t nElectrons;

    ////////////////////////////////


private:


    // weights and uncertainties
    Float_t puWeight, puVar;
    Float_t leptonOneIDWeight, leptonOneIDVar;
    Float_t leptonTwoIDWeight, leptonTwoIDVar;
    Float_t leptonOneRecoWeight, leptonOneRecoVar;
    Float_t leptonTwoRecoWeight, leptonTwoRecoVar;



    // helper functions
    float GetMuonIsolation(const baconhep::TMuon*);
    float GetElectronIsolation(const baconhep::TElectron*, float);

};

#endif  // ZEEAnalyzer_HH
