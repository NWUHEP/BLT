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


class MuonTagAndProbe: public BLTSelector {
public:
    MuonTagAndProbe();
    ~MuonTagAndProbe();

    void   Begin(TTree *tree);
    Bool_t Process(Long64_t entry);
    void   Terminate();
    void   ReportPostBegin();
    void   ReportPostTerminate();

    TFile *outFile;
    TTree *evtTree, *muTree;

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

    Float_t  mu_pt, mu_eta, mu_phi, mu_ptErr;                   // kinematics
    Float_t  mu_genPt, mu_genEta, mu_genPhi;                   // kinematics
    Float_t  mu_staPt, mu_staEta, mu_staPhi;                 // STA track kinematics
    Float_t  mu_pfPt, mu_pfEta, mu_pfPhi;                    // matched PFCandidate
    Float_t  mu_trkIso, mu_ecalIso, mu_hcalIso;              // detector isolation (R=0.3)
    Float_t  mu_chHadIso, mu_gammaIso, mu_neuHadIso, mu_puIso;  // PF isolation variables (R=0.4)
    Float_t  mu_puppiChHadIso, mu_puppiGammaIso, mu_puppiNeuHadIso;  // Puppi Isolation R=0.4
    Float_t  mu_puppiChHadIsoNoLep, mu_puppiGammaIsoNoLep, mu_puppiNeuHadIsoNoLep; // Puppi Isolation R=0.4 no lep
    Float_t  mu_d0, mu_dz, mu_sip3d;                         // impact parameter
    Float_t  mu_tkNchi2, mu_muNchi2;                      // track fit normalized chi-square
    Float_t  mu_trkKink, mu_glbKink;                      // track kink
    Float_t  mu_trkHitFrac;                            // fraction of valid tracker hits
    Float_t  mu_chi2LocPos;                            // TRK-STA position match
    Float_t  mu_segComp;                               // compatibility of tracker track with muon segment
    Float_t  mu_caloComp;                              // muon hypothesis compatibility with calo energy
    Int_t    mu_q;                                     // charge
    Int_t    mu_nValidHits;                            // number of valid muon hits in global fit
    UInt_t   mu_nTkHits, mu_nPixHits;                     // number of hits in tracker
    UInt_t   mu_nTkLayers, mu_nPixLayers;                 // number of hit layers in tracker
    UInt_t   mu_nMatchStn;                             // number of stations with muon segments
    UInt_t   mu_typeBits;                              // muon type bits
    UInt_t   mu_selectorBits;                          // MuonSelector bits
    UInt_t   mu_pogIDBits;                             // POG muon IDs from CMSSW
    Bool_t   mu_genMatched;

    Float_t met, metPhi;
    UInt_t nJets, nFwdJets, nBJets, nMuons, nElectrons;

    // MET kluge 
    float MetKluge(float);
    void  FillMuonTree(baconhep::TMuon&, TLorentzVector&, bool);

    //ClassDef(MuonTagAndProbe,0);
};


#endif  // MULTILEPTONANALYZER_HH
