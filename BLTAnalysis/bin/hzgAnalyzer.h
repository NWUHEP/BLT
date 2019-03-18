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
    Float_t xPV, yPV, zPV;
    UInt_t nJets, nCentralJets, nFwdJets, nBJets, nMuons, nElectrons, nTaus, nPhotons;
   
    // weights
    Int_t genWeight;
    Float_t eventWeight, triggerWeight, puWeight;
    Float_t elIDWeightOne, elIDWeightTwo;
    Float_t elTrigWeightOne, elTrigWeightTwo;
    Float_t muonIDWeightOne, muonIDWeightTwo;
    Float_t muonTrigWeightOne, muonTrigWeightTwo;
    Float_t photonIDWeight;

    // physics object Lorentz vectors
    Float_t leptonOnePt, leptonOneEta, leptonOnePhi;
    Float_t leptonTwoPt, leptonTwoEta, leptonTwoPhi;

    Float_t leptonOnePtKin, leptonTwoPtKin;
    Float_t leptonOnePtKinJames, leptonTwoPtKinJames;

    // Additional lepton data
    Float_t leptonOneIso, leptonTwoIso;
    Int_t leptonOneMother, leptonTwoMother;
    Int_t leptonOneFlavor, leptonTwoFlavor;
    Float_t leptonOneD0, leptonTwoD0;
    Float_t leptonOneDZ, leptonTwoDZ;
    Float_t leptonOneRecoWeight, leptonTwoRecoWeight;

    Bool_t leptonOneECALDriven, leptonTwoECALDriven;
    
    // tau data
    Int_t tauDecayMode;
    Float_t tauMVA; 
    //UInt_t tauPhotonMult, tauChHadMult;

    // photon data
    Float_t photonOnePt, photonOneEta, photonOnePhi;
    Float_t photonOneR9;
    Float_t photonOneMVA;
    Float_t photonOneERes;
    Bool_t passElectronVeto;

    Bool_t isLeptonTag;
    Bool_t isDijetTag;
    Bool_t isTightDijetTag;

    // jet data
    Float_t jetOnePt, jetOneEta, jetOnePhi, jetOneM;
    Float_t jetTwoPt, jetTwoEta, jetTwoPhi, jetTwoM;
    Float_t jetOneTag, jetTwoTag, jetThreeTag, jetFourTag;
    Float_t met, metPhi, metNC, metPhiNC, ht, htPhi, htSum;

    // generator level data
    Float_t genLeptonOnePt, genLeptonOneEta, genLeptonOnePhi;
    Float_t genLeptonTwoPt, genLeptonTwoEta, genLeptonTwoPhi;
    Float_t genPhotonPt, genPhotonEta, genPhotonPhi;
    Int_t genLeptonOneId, genLeptonTwoId;
    Bool_t genPhotonFHPFS, genPhotonIPFS;
    
    //Int_t genOneId, genTwoId, genOneMother, genTwoMother, genCategory;
    //TLorentzVector genOneP4, genTwoP4;
    //Bool_t fromHardProcessFinalState, isPromptFinalState, hasPhotonMatch;
    Bool_t vetoDY;

    // dilepton data
    Float_t dileptonPt, dileptonEta, dileptonPhi, dileptonM;
    Float_t dileptonDEta, dileptonDPhi, dileptonDR;
    Float_t dileptonMKin;
    Float_t dileptonMKinJames;
    
    // dilepton vertex data
    //Float_t dileptonVertexOneX, dileptonVertexOneY, dileptonVertexOneZ;
    //Float_t dileptonVertexTwoX, dileptonVertexTwoY, dileptonVertexTwoZ;
    //Float_t dileptonVertexOneXErr, dileptonVertexOneYErr, dileptonVertexOneZErr;
    //Float_t dileptonVertexTwoXErr, dileptonVertexTwoYErr, dileptonVertexTwoZErr;
    //Float_t dileptonVertexChi2One, dileptonVertexDOFOne;
    //Float_t dileptonVertexChi2Two, dileptonVertexDOFTwo;
    
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
    Float_t llgMKinJames;
    Float_t l1PhotonDEta, l1PhotonDPhi, l1PhotonDR;
    Float_t l2PhotonDEta, l2PhotonDPhi, l2PhotonDR;
    Float_t lPhotonDRMax, lPhotonDRMin;
    Float_t dileptonPhotonDEta, dileptonPhotonDPhi, dileptonPhotonDR;
    Float_t ptt;
    Float_t zgBigTheta, zgLittleTheta, zgPhi;
    Float_t zgLittleThetaMY;
    Float_t zgBigThetaJames, zgLittleThetaJames, zgPhiJames;
    Float_t genBigTheta, genLittleTheta, genPhi;

    // other
    Float_t llgJJDEta, llgJJDPhi, llgJJDR;
    Float_t zepp;

    float GetMuonIsolation(const baconhep::TMuon*);
    float GetElectronIsolation(const baconhep::TElectron*, float);
    float GetPhotonIsolation(const baconhep::TPhoton*, float);

    void EvalMuonEnergyResolution(std::map<string, float>, std::map<string, int>, float&, float&, float&, float&, float&, float&);
    void EvalElectronEnergyResolution(std::map<string, float>, float&, float&, float&, float&, float&, float&);
    void find_optimized(double*, double&, double&);

    //ClassDef(hzgAnalyzer,0);
};

double GetDoubleSidedCB(double x, double mean, double sigma, double alphaL,
                                          double powerL, double alphaR, double powerR)
{
   // Returns value of double-sided Crystal Ball function for given parameters.
   //
   // Negative powerL and/or powerR means an infinite value.
    double a = (x - mean) / sigma;
    // left power-law or exponential tail
   if (a < -alphaL) {
      if (powerL < 0) // infinite powerL
         return exp(0.5 * alphaL*alphaL + alphaL*a);
       double b = powerL/alphaL;
      return exp(-0.5 * alphaL*alphaL) * pow(b/(b - alphaL - a), powerL);
   }
    // Gaussian core
   if (a <= alphaR)
      return exp(-0.5*a*a);
    // right exponential tail
   if (powerR < 0) // infinite powerR
      return exp(0.5 * alphaR*alphaR - alphaR*a);
    // right power-law tail
   double b = powerR/alphaR;
   return exp(-0.5 * alphaR*alphaR) * pow(b/(b - alphaR + a), powerR);
}
 Double_t NegativeProbability(Double_t* x, Double_t* p)
{
   // Negative probability density function to minimize.
    // lepton energies to optimize
   double e1 = x[0];
   double e2 = x[1];
    double m02 = p[0];       // lepton mass squared
   double cosAlpha = p[1];  // cosine of the angle between lepton directions
    // two-lepton invariant mass squared
   double mz2 = 2 * (e1*e2 + m02
                     - sqrt(e1*e1 - m02) * sqrt(e2*e2 - m02) * cosAlpha);
    // relativistic Breit-Wigner
   double a = mz2 - MZ*MZ;
   double probZ = 1/(a*a + mz2*mz2 * WZ*WZ/(MZ*MZ));
    // p[2] = first reconstructed energy
   // p[3] = second reconstructed energy
   // p[4], ... = parameters of energy resolution functions
    double prob1 = GetDoubleSidedCB(e1/p[2], p[4], p[5], p[6], p[7], p[8], p[9]);
   double prob2 = GetDoubleSidedCB(e2/p[3], p[10], p[11], p[12], p[13], p[14], p[15]);
    return -probZ * prob1 * prob2;
}

#endif  // HZGANALYZER_HH
