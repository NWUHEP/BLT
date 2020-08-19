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


#ifndef HZGANALYZER_HH
#define HZGANALYZER_HH

// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"
#include "BLT/BLTAnalysis/interface/TriggerSelector.hh"
#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"
#include "BLT/BLTAnalysis/interface/WeightUtils.h"

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
#define MZ 91.1876
#define WZ 2.4952

class hzgAnalyzer: public BLTSelector {
public:
    hzgAnalyzer();
    ~hzgAnalyzer();

    void   Begin(TTree *tree);
    Bool_t Process(Long64_t entry);
    void   Terminate();
    void   ReportPostBegin();
    void   ReportPostTerminate();

    TFile *outFile;
    std::map<string, TTree*> outTrees;
    std::map<string, TH1D*> eventCounts;
    //TTree *outTree;

    // Lumi mask
    RunLumiRangeMap lumiMask;

    TRandom3 *rng;

    // Params 
    std::unique_ptr<Parameters>         params;

    // Utilities
    std::unique_ptr<ParticleSelector>   particleSelector;
    std::unique_ptr<baconhep::TTrigger> trigger;
    std::unique_ptr<WeightUtils>        weights;

    std::vector<string> triggerNames;

    // Branches in the output file
    
    // event data
    UInt_t runNumber, lumiSection, nPV, nPartons;
    ULong64_t evtNumber;
    Float_t nPU;
    Float_t xPV, yPV, zPV;
    UInt_t nJets, nCentralJets, nFwdJets, nBJets, nMuons, nElectrons, nTaus, nPhotons;
    Float_t rho;
   
    // weights
    Int_t genWeight;
    Float_t eventWeight, triggerWeight;
    Float_t puWeight, puWeightUp, puWeightDown;
    Float_t prefWeight, prefWeightUp, prefWeightDown;
    Float_t leptonOneRecoWeight, leptonTwoRecoWeight;
    Float_t leptonOneRecoErr, leptonTwoRecoErr;
    Float_t leptonOneIDWeight, leptonTwoIDWeight;
    Float_t leptonOneIDErr, leptonTwoIDErr;
    Float_t leptonOneRecoIDWeight, leptonTwoRecoIDWeight;
    Float_t trigOneWeight, trigTwoWeight;
    Float_t trigOneErr, trigTwoErr;
    Float_t photonCSEVWeight, photonIDWeight;
    Float_t photonCSEVErr, photonIDErr;
    Float_t photonR9Weight;
    Float_t addLeptonRecoWeight, addLeptonIDWeight, addLeptonRecoIDWeight;
    Int_t addLeptonFlavor;

    // physics object Lorentz vectors
    Float_t leptonOnePt, leptonOneEta, leptonOnePhi;
    Float_t leptonTwoPt, leptonTwoEta, leptonTwoPhi;
    Float_t leptonOnePtKin, leptonTwoPtKin;
    Float_t leptonOnePtKinErr, leptonTwoPtKinErr;
    Float_t leptonOnePtKinRes, leptonTwoPtKinRes;

    // Additional lepton data
    Float_t leptonOneIso, leptonTwoIso;
    Int_t leptonOneFlavor, leptonTwoFlavor;
    Float_t leptonOneD0, leptonTwoD0;
    Float_t leptonOneDZ, leptonTwoDZ;
    
    // tau data
    Int_t tauDecayMode;
    Float_t tauMVA; 
    //UInt_t tauPhotonMult, tauChHadMult;

    // photon data
    Float_t photonPt, photonEta, photonPhi;
    Float_t photonRawPt;
    Float_t photonSCPt, photonSCEta, photonSCPhi;
    Float_t photonR9, photonRawR9;
    Float_t photonMVA;
    Float_t photonERes;
    Float_t photonE, photonErrE;
    //Float_t photonWorstChIso; // on hold
    Bool_t passElectronVeto;

    Bool_t isLeptonTag;
    Bool_t isDijetTag;
    Bool_t isTightDijetTag;

    // jet data
    Float_t jetOnePt, jetOneEta, jetOnePhi, jetOneM;
    Float_t jetOneL1Corr, jetOneL2Corr, jetOneL3Corr, jetOneL4Corr;
    Float_t jetOneArea, jetOneRawPt;
    Float_t jetTwoPt, jetTwoEta, jetTwoPhi, jetTwoM;
    Float_t jetTwoL1Corr, jetTwoL2Corr, jetTwoL3Corr, jetTwoL4Corr;
    Float_t jetTwoArea, jetTwoRawPt;
    Float_t jetOneTag, jetTwoTag, jetThreeTag, jetFourTag;
    Float_t met, metPhi, metNC, metPhiNC, ht, htPhi, htSum;

    // generator level data

    UInt_t nGenPhotons, nGenMuons, nGenElectrons, nGenTaus;
    Int_t genLeptonOneId, genLeptonTwoId;
    Float_t genLeptonOnePt, genLeptonOneEta, genLeptonOnePhi;
    Float_t genLeptonTwoPt, genLeptonTwoEta, genLeptonTwoPhi;
    Float_t genPhotonPt, genPhotonEta, genPhotonPhi;
    Bool_t vetoDY;

    // dilepton data
    Float_t dileptonPt, dileptonEta, dileptonPhi, dileptonM;
    Float_t dileptonDEta, dileptonDPhi, dileptonDR;
    Float_t dileptonMKin;
     
    // dijet data
    Float_t dijetPt, dijetEta, dijetPhi, dijetM;
    Float_t dijetDEta, dijetDPhi, dijetDR;

    // jet, lepton data
    Float_t l1j1DEta, l1j1DPhi, l1j1DR;
    Float_t l1j2DEta, l1j2DPhi, l1j2DR;
    Float_t l2j1DEta, l2j1DPhi, l2j1DR;
    Float_t l2j2DEta, l2j2DPhi, l2j2DR;

    // jet, photon data
    Float_t j1PhotonDEta, j1PhotonDPhi, j1PhotonDR;
    Float_t j2PhotonDEta, j2PhotonDPhi, j2PhotonDR;
    Float_t jPhotonDRMax, jPhotonDRMin;

    // three body
    Float_t llgPt, llgEta, llgPhi, llgM, llgPtOverM;
    Float_t llgMKin;
    Float_t l1PhotonDEta, l1PhotonDPhi, l1PhotonDR;
    Float_t l2PhotonDEta, l2PhotonDPhi, l2PhotonDR;
    Float_t lPhotonDRMax, lPhotonDRMin;
    Float_t dileptonPhotonDEta, dileptonPhotonDPhi, dileptonPhotonDR;
    Float_t ptt;

    // angles
    Float_t zgBigTheta, zgLittleTheta, zgPhi;
    Float_t zgBigThetaMY, zgLittleThetaMY, zgPhiMY;
    Float_t zgBigThetaJames, zgLittleThetaJames, zgPhiJames;
    Float_t genBigTheta, genLittleTheta, genPhi;

    // other
    Float_t llgJJDEta, llgJJDPhi, llgJJDR;
    Float_t zepp, photonZepp;
    Float_t vbfPtBalance;
    Bool_t jetOneMatched, jetTwoMatched;
    Bool_t leptonOneMatched, leptonTwoMatched;
    Bool_t photonMatched;

    //ClassDef(hzgAnalyzer,0);
};

#endif  // HZGANALYZER_HH
