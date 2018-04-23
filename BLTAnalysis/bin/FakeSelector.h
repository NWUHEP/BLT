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


#ifndef FAKESELECTOR_HH
#define FAKESELECTOR_HH

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


class FakeSelector: public BLTSelector {
public:
    FakeSelector();
    ~FakeSelector();

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
    UInt_t runNumber, lumiSection, nPV, nPartons;
    ULong64_t evtNumber;
    Bool_t triggerStatus;
    Float_t nPU;
    TVector3 rPV;
    UInt_t nJets, nFwdJets, nBJets, nMuons, nElectrons, nPhotons;

    // weights and uncertainties
    Float_t eventWeight;

    // physics object Lorentz vectors
    TLorentzVector leptonOneP4, leptonTwoP4, leptonThreeP4;

    // Additional lepton data
    Float_t leptonOneIso, leptonTwoIso, leptonThreeIso;
    Int_t leptonOneMother, leptonTwoMother, leptonThreeMother; 
    Int_t leptonOneFlavor, leptonTwoFlavor, leptonThreeFlavor;
    Float_t leptonOneD0, leptonTwoD0, leptonThreeD0;
    Float_t leptonOneDZ, leptonTwoDZ, leptonThreeDZ;


    // jet data
    TLorentzVector jetOneP4, jetTwoP4, jetThreeP4, jetFourP4;
    Float_t jetOneTag, jetTwoTag, jetThreeTag, jetFourTag;
    Int_t jetOneFlavor, jetTwoFlavor, jetThreeFlavor, jetFourFlavor;
    Float_t met, metPhi, ht, htPhi, htSum;

    // generator level data
    Int_t genOneId, genTwoId, genOneMother, genTwoMother, genCategory;
    TLorentzVector genOneP4, genTwoP4;

    float MetKluge(float);
    float GetMuonIsolation(const baconhep::TMuon*);
    float GetElectronIsolation(const baconhep::TElectron*, float);
    float GetTriggerSF(EfficiencyContainer, EfficiencyContainer);
    float GetTriggerSFError(EfficiencyContainer, EfficiencyContainer);
    vector<unsigned> PairDileptonToZ(vector<TLorentzVector>);
    int GetGenMotherId(vector<baconhep::TGenParticle*>, TLorentzVector);
    vector<baconhep::TJet*> KinematicTopTag(vector<baconhep::TJet*>, TVector2, TLorentzVector);

    //ClassDef(FakeSelector,0);
};


#endif  // FAKESELECTOR_HH
