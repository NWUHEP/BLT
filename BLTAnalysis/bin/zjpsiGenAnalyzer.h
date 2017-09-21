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


#ifndef ZJPSIGENANALYZER_HH
#define ZJPSIGENANALYZER_HH

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
#include <map>
#include <memory>
#include <iostream>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <regex>


class zjpsiGenAnalyzer: public BLTSelector {
public:
    zjpsiGenAnalyzer();
    ~zjpsiGenAnalyzer();

    void   Begin(TTree *tree);
    Bool_t Process(Long64_t entry);
    void   Terminate();
    void   ReportPostBegin();
    void   ReportPostTerminate();

    TFile *outFile;
    TTree *outTree;

    // Stuff for binned accepances and efficiencies
    std::map<std::pair<int, int>, std::vector<int>> cutMap;
    FILE *csvOut;

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
    UInt_t runNumber, lumiSection, nPV, nPartons;
    ULong64_t evtNumber;
    Bool_t triggerStatus;
    Float_t eventWeight, nPU;

    TVector3 rPV;
    TVector3 rDimuon;
    TVector3 rZCand;
    TVector3 rJpsiCand;

    TVector3 rPVErr;
    TVector3 rDimuonErr;
    TVector3 rZCandErr;
    TVector3 rJpsiCandErr;
    

    Float_t rPVChi2, rPVNdof;
    Float_t rDimuonChi2, rDimuonNdof;
    Float_t rZCandChi2, rZCandNdof;
    Float_t rZCandProb, rZCandRxy, rZCandRxyErr;
    Float_t rJpsiCandChi2, rJpsiCandNdof;
    Float_t rJpsiCandProb, rJpsiCandRxy, rJpsiCandRxyErr;

    TLorentzVector leptonOneP4, leptonTwoP4, leptonThreeP4, leptonFourP4, jetP4, bjetP4, genJetP4, genBJetP4, photonOneP4;
    Float_t leptonOneIso, leptonTwoIso, leptonThreeIso, leptonFourIso;
    Float_t leptonOneTkIsoNR, leptonTwoTkIsoNR, leptonThreeTkIsoNR, leptonFourTkIsoNR;
    Float_t leptonOneTkIso, leptonTwoTkIso, leptonThreeTkIso, leptonFourTkIso;
    Float_t leptonOnePFIso, leptonTwoPFIso, leptonThreePFIso, leptonFourPFIso;
    Int_t leptonOneQ, leptonTwoQ, leptonThreeQ, leptonFourQ;
    Int_t leptonOneFlavor, leptonTwoFlavor, leptonThreeFlavor, leptonFourFlavor, jetFlavor, bjetFlavor;
    Bool_t leptonOneTrigger, leptonTwoTrigger, leptonThreeTrigger, leptonFourTrigger;
    Bool_t leptonOneIsTight, leptonTwoIsTight, leptonThreeIsTight, leptonFourIsTight;

    Float_t leptonOneD0, leptonTwoD0, leptonThreeD0, leptonFourD0;
    Float_t leptonOneDZ, leptonTwoDZ, leptonThreeDZ, leptonFourDZ;

    Float_t jetD0, bjetD0;
    Float_t bjetPUID, bjetTag, jetPUID, jetTag, genJetTag, genBJetTag;
    Float_t met, metPhi;

    UInt_t nJets, nFwdJets, nBJets, nMuons, nElectrons, nPhotons;

    // MET kluge 
    float MetKluge(float);

    //ClassDef(zjpsiGenAnalyzer,0);
};


#endif  // ZJPSIGENANALYZER_HH
