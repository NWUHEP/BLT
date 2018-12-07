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

    // Branches in the output file
    // event data
    UInt_t runNumber, lumiSection, nPV, nPartons,triggerLeptonStatus;
    ULong64_t evtNumber;
    Bool_t triggerStatus;
    Float_t nPU;
    TVector3 rPV;
    UInt_t nJets, nFwdJets, nBJets, nMuons, nElectrons, nTaus, nPhotons;

    // modified multiplicities for jet related uncertainties
    unsigned nJetsCut,nJetsCut2, nJetsJESUp, nJetsJESDown, nJetsJERUp, nJetsJERDown;
    unsigned nBJetsCut, nBJetsJESUp, nBJetsJESDown, nBJetsJERUp, nBJetsJERDown;
    unsigned nBJetsBTagUp, nBJetsBTagDown, nBJetsMistagUp, nBJetsMistagDown;

    // weights and uncertainties
    Float_t eventWeight, triggerWeight, puWeight, topPtWeight, genWeight;
    Float_t leptonOneRecoWeight, leptonTwoRecoWeight, leptonThreeRecoWeight, leptonFourRecoWeight;
    Float_t triggerVar, puVar, topPtVar, leptonOneRecoVar, leptonTwoRecoVar;

    vector<float> qcdWeights;
    float pdfWeight, alphaS;

    // physics object Lorentz vectors
    TLorentzVector leptonOneP4, leptonTwoP4, leptonThreeP4, leptonFourP4, tauP4;

    // Additional lepton data
    Float_t leptonOneIso, leptonTwoIso, leptonThreeIso, leptonFourIso;
    Int_t leptonOneMother, leptonTwoMother, leptonThreeMother, leptonFourMother; 
    Int_t leptonOneFlavor, leptonTwoFlavor, leptonThreeFlavor, leptonFourFlavor;
    Float_t leptonOneD0, leptonTwoD0, leptonThreeD0, leptonFourD0;
    Float_t leptonOneDZ, leptonTwoDZ, leptonThreeDZ, leptonFourDZ;

    // tau variables
    Int_t tauDecayMode, tauGenFlavor, tauGenFlavorHad;
    Float_t tauMVA, tauPuppiChHadIso, tauPuppiGammaIso, tauPuppiNeuHadIso, tauVetoedJetPT, tauVetoedJetPTUnc;

    // dimuon vertex data
    TVector3 dileptonVertexOne, dileptonVertexTwo, dileptonVertexErrOne, dileptonVertexErrTwo;
    Float_t dileptonVertexChi2One, dileptonVertexDOFOne;
    Float_t dileptonVertexChi2Two, dileptonVertexDOFTwo;

    // jet data
    TLorentzVector jetOneP4, jetTwoP4, jetThreeP4, jetFourP4;
    Float_t jetOneTag, jetTwoTag, jetThreeTag, jetFourTag;
    Int_t jetOneFlavor, jetTwoFlavor, jetThreeFlavor, jetFourFlavor;
    Float_t met, metPhi, ht, htPhi, htSum;

    // generator level data
    Int_t genOneId, genTwoId, genOneMother, genTwoMother, genCategory;
    TLorentzVector genOneP4, genTwoP4;

    // helper functions
    float MetKluge(float);
    float GetMuonIsolation(const baconhep::TMuon*);
    float GetElectronIsolation(const baconhep::TElectron*, float);
    float GetTriggerSF(EfficiencyContainer, EfficiencyContainer);
    float GetTriggerSFError(EfficiencyContainer, EfficiencyContainer);
    vector<unsigned> PairDileptonToZ(vector<TLorentzVector>);
    int GetGenMotherId(vector<baconhep::TGenParticle*>, TLorentzVector);
    int GetTauGenFlavor(TLorentzVector, vector<baconhep::TGenParticle*>, vector<baconhep::TJet*>, bool );
    pair<float, float> GetTauVetoedJetPT(TLorentzVector, vector<TJet*> );
    vector<baconhep::TJet*> KinematicTopTag(vector<baconhep::TJet*>, TVector2, TLorentzVector);
    void ResetJetCounters();
    void JetCounting(TJet* jet, float , float);

    //ClassDef(MultileptonAnalyzer,0);
};

#endif  // MULTILEPTONANALYZER_HH
