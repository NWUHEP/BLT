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


#ifndef MultilepAnalyzer_HH
#define MultilepAnalyzer_HH


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


class MultilepAnalyzer: public BLTSelector {
public:
    MultilepAnalyzer();
    ~MultilepAnalyzer();

    void   Begin(TTree *tree);
    Bool_t Process(Long64_t entry);
    void   Terminate();
    void   ReportPostBegin();
    void   ReportPostTerminate();

    TFile *outFile;
    std::map<string, TTree*> outTrees;
    std::map<string, TH1D*> eventCounts;
    TH1D *hGenCat;

    std::string singleElectronTriggerName = "";

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
    UInt_t runNumber, lumiSection, nPV, triggerLeptonStatus, mcEra;
    ULong64_t evtNumber;
    Float_t nPU;

    // generator level data
    UInt_t nPartons;
    Float_t genWeight;
    UInt_t genCategory;
    UInt_t genTauOneDaughters, genTauTwoDaughters;

    // weights and uncertainties
    Float_t eventWeight;

    // lepton variable
    TLorentzVector leptonOneP4, leptonTwoP4;
    Float_t leptonOneIso, leptonTwoIso;
    Int_t leptonOneFlavor, leptonTwoFlavor;

    // tau variable
    Int_t tauDecayMode, tauGenFlavor, tauGenFlavorHad;
    Float_t tauMVA, tauVetoedJetPt, tauVetoedJetPtUnc;

    // ht
    Float_t htSum, ht, htPhi;

    // met
    Float_t met, metPhi;

    // object counters
    UInt_t nMuons, nElectrons, nTaus, nJets, nBJets;
    UInt_t nFailMuons, nFailElectrons;

    ////////////////////////////////


private:


    // weights and uncertainties
    Float_t triggerWeight, puWeight, topPtWeight, zPtWeight,leptonOneIDWeight, leptonTwoIDWeight, leptonOneRecoWeight, leptonTwoRecoWeight;
    Float_t triggerVar, puVar, topPtVar, zPtVar, leptonOneIDVar, leptonTwoIDVar, leptonOneRecoVar, leptonTwoRecoVar;
    Float_t eleTriggerVarTagSyst, eleTriggerVarProbeSyst;


    // modified multiplicities for jet related uncertainties
    unsigned nJetsCut, nBJetsCut;
    unsigned nJetsJESUp,  nJetsJESDown,  nJetsJERUp,  nJetsJERDown;
    unsigned nBJetsJESUp, nBJetsJESDown, nBJetsJERUp, nBJetsJERDown;
    unsigned nBJetsBTagUp, nBJetsBTagDown, nBJetsMistagUp, nBJetsMistagDown;

    // helper functions
    float GetMuonIsolation(const baconhep::TMuon*);
    float GetElectronIsolation(const baconhep::TElectron*, float);
    int GetTauGenFlavor(TLorentzVector, vector<baconhep::TGenParticle*>,  vector<TGenParticle*>, vector<TGenParticle*>, vector<baconhep::TJet*>, bool );
    pair<float, float> GetTauVetoedJetPt(TLorentzVector, vector<TJet*> );
    float GetTriggerSF(EfficiencyContainer, EfficiencyContainer);
    float GetTriggerSFError(EfficiencyContainer, EfficiencyContainer);
    void ResetJetCounters();
    void JetCounting(TJet* jet, float , float);

};

#endif  // MultilepAnalyzer_HH
