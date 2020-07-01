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


#ifndef ZTauTauAnalyzer_HH
#define ZTauTauAnalyzer_HH


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

// SVFit headers
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

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


class ZTauTauAnalyzer: public BLTSelector {
public:
    ZTauTauAnalyzer();
    ~ZTauTauAnalyzer();

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
    UInt_t runNumber, lumiSection, nPV, triggerLeptonStatus, mcEra;
    ULong64_t evtNumber;
    Float_t nPU;

    // generator level data
    UInt_t nPartons;
    Float_t genWeight;
    UInt_t genCategory;
    UInt_t nGenElectrons, nGenMuons, nGenTausHad, nGenTausLep;
    UInt_t nGenPromptElectrons, nGenPromptMuons, nGenPromptTaus;

    // weights and uncertainties
    Float_t eventWeight;
    Float_t genTauFlavorWeight;

    //for jets study channel
    const static UInt_t kMaxJets = 20; //maximum of 20 jets saved ever
    Float_t jetsPt[kMaxJets], jetsEta[kMaxJets], jetsPhi[kMaxJets], jetsM[kMaxJets], jetsbMVA[kMaxJets], jetsGenFlavor[kMaxJets];
    
    // lepton variable
    TLorentzVector leptonOneP4, leptonTwoP4, photonOneP4, jetOneP4, jetTwoP4, tauOneP4;    
    TLorentzVector genLeptonOneP4, genLeptonTwoP4;
    Float_t leptonOneIso, leptonTwoIso, photonMVA;
    Float_t leptonOneD0, leptonTwoD0;
    Int_t leptonOneFlavor, leptonTwoFlavor;

    // tau variable
    Int_t tauDecayMode, tauGenFlavor, tauGenFlavorHad, tauFlavor;
    Float_t tauGenPt, tauGenEta;
    Float_t tauPt, tauEta, tauMVA, tauVetoedJetPt, tauVetoedJetPtUnc;
    unsigned long tauHPSDisc; //tau ids
    
    // ht
    Float_t htSum, ht, htPhi;
    Float_t htFwdSum, htFwd, htFwdPhi;

    // met
    Float_t met, metPhi;
    Float_t covMet00, covMet01, covMet11;
    // Particle flow met
    Float_t pfMET, pfMETphi, pfMETCov00, pfMETCov01, pfMETCov11;
    Float_t pfMETC, pfMETCphi, pfMETCCov00, pfMETCCov01, pfMETCCov11;
    //  PUPPI met
    Float_t puppMET, puppMETphi, puppMETCov00, puppMETCov01, puppMETCov11;
    Float_t puppMETC, puppMETCphi, puppMETCCov00, puppMETCCov01, puppMETCCov11;
    // Alapaca met
    Float_t alpacaMET, alpacaMETphi;
    // Alapaca + PUPPI met
    Float_t pcpMET, pcpMETphi;
    // Trk met
    Float_t trkMET, trkMETphi;

    // Correction to MET due to Pt adjustments of reconstructed objects
    Float_t metCorr, metCorrPhi;
    
    // object counters
    UInt_t nMuons, nElectrons, nTaus,nPhotons, nJets, nFwdJets, nBJets, nBJetsM, nBJetsL;
    UInt_t nJets25, nFwdJets25, nBJets25, nBJets25M, nBJets25L;
    UInt_t nJets20, nFwdJets20, nBJets20, nBJets20M, nBJets20L;

    UInt_t nFailMuons, nFailElectrons;
    UInt_t nLowPtElectrons;
    
    //SVFit variables
    Float_t massSVFit, massErrSVFit;
    Int_t svFitStatus;
    TLorentzVector leptonOneSVP4, leptonTwoSVP4; //corrected lorentz vectors    
    ////////////////////////////////

    bool saveExtraInfo; //for debugging purposes
    bool dollgStudy; //for saving info for llg study
  
private:

    //deleting svfit histograms
    Int_t trkhistos = 0;
    // weights and uncertainties
    Float_t triggerWeight, puWeight, topPtWeight, leptonOneIDWeight, leptonTwoIDWeight, leptonOneRecoWeight, leptonTwoRecoWeight;
    Float_t triggerVar, puVar, topPtVar, zPtVar, leptonOneIDVar, leptonTwoIDVar, leptonOneRecoVar, leptonTwoRecoVar;

    // for Drell-Yan pT
    Float_t zPt, zDiLeptonPt, zPtWeight;

    // modified multiplicities for jet related uncertainties
    unsigned nJetsCut, nBJetsCut;
    unsigned nJetsJESUp,  nJetsJESDown,  nJetsJERUp,  nJetsJERDown;
    unsigned nBJetsJESUp, nBJetsJESDown, nBJetsJERUp, nBJetsJERDown;
    unsigned nBJetsBTagUp, nBJetsBTagDown, nBJetsMistagUp, nBJetsMistagDown;

    // helper functions
    float GetMuonIsolation(const baconhep::TMuon*);
    float GetElectronIsolation(const baconhep::TElectron*, float);
    pair<int,TLorentzVector> GetTauGenFlavor(TLorentzVector, vector<baconhep::TGenParticle*>,  vector<TGenParticle*>, vector<TGenParticle*>, vector<baconhep::TJet*>, bool );
    double FakeTauSF(int pdgid, double pt, double eta, bool isVeryTight);
    pair<float, float> GetTauVetoedJetPt(TLorentzVector, vector<TJet*> );
    float GetZPtWeight(float zpt);
    float GetTriggerSF(EfficiencyContainer, EfficiencyContainer);
    float GetTriggerSFError(EfficiencyContainer, EfficiencyContainer);
    void ResetJetCounters();
    bool JetCounting(TJet* jet, float , float);

};

#endif  // ZTauTauAnalyzer_HH
