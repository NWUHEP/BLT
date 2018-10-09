#ifndef PARTICLESELECTOR_HH
#define PARTICLESELECTOR_HH


// Bacon header files
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"
#include "BLT/BLTAnalysis/interface/Cuts.hh"

#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TRandom3.h"

#include <string>
#include <vector>
#include <memory>
#include <cassert>

//CMSSW libraries
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

using namespace std;

class ParticleSelector {
public:
    ParticleSelector(const Parameters& parameters, const Cuts& cuts);
    ~ParticleSelector() {}

    // Setters
    void SetRealData(bool isRealData)       { _isRealData = isRealData; }
    void SetPV(const TVector3& pv)          { _pv = pv; }
    void SetNPV(int npv)                    { _npv = npv; }
    void SetRho(float rhoFactor)            { _rhoFactor = rhoFactor; }

    // Muons
    bool PassMuonID(const baconhep::TMuon* mu, const Cuts::muIDCuts& cutLevel) const;
    bool PassMuonIso(const baconhep::TMuon* mu, const Cuts::muIsoCuts& cutLevel) const;
    bool PassMuonIso(const baconhep::TMuon* mu, const Cuts::muDetIsoCuts& cutLevel) const;

    // Electrons
    bool PassElectronID(const baconhep::TElectron* el, const Cuts::elIDCuts& cutLevel) const;
    bool PassElectronMVA(const baconhep::TElectron* el, const Cuts::elMVACuts& cutLevel) const;
    bool PassElectronIso(const baconhep::TElectron* el, const Cuts::elIsoCuts& cutLevel, float EAEl[7]) const;

    // Photons
    bool PassPhotonID(const baconhep::TPhoton* ph, const Cuts::phIDCuts& cutLevel) const;
    bool PassPhotonMVA(const baconhep::TPhoton* ph, const Cuts::phMVACuts& cutLevel) const;
    bool PassPhotonIso(const baconhep::TPhoton* ph, const Cuts::phIsoCuts& cutLevel, float EAPho[7][3]) const;

    // Jets
    bool PassJetID(const baconhep::TJet* jet, const Cuts::jetIDCuts& cutLevel) const;
    bool PassJetPUID(const baconhep::TJet* jet) const;
    bool BTagModifier(const baconhep::TJet* jet, string, int, int, float) const;
    double JetCorrector(const baconhep::TJet* jet, string) const;
    double JetUncertainty(const baconhep::TJet* jet, string) const;
    pair<float, float> JetResolutionAndSF(const baconhep::TJet* jet, int) const;
    vector<string> GetJECSourceNames() { return this->_jecNames;}

private:
    Parameters _parameters;
    Cuts       _cuts;
    bool       _isRealData;
    TVector3   _pv;
    int        _npv;
    float      _rhoFactor;
    TRandom3*  _rng;

    // For offline jet corrections
    FactorizedJetCorrector* _jetCorrector;
    map<string, JetCorrectionUncertainty*> _jecUncertaintyMap;
    JME::JetResolution jetResolution;
    JME::JetResolutionScaleFactor jetResolutionSF;

    // jec uncertainty names
    vector<string> _jecNames {
                              "AbsoluteScale", "AbsoluteStat",
                              "Fragmentation", "AbsoluteMPFBias", //"CorrelationGroupIntercalibration", 
                              "PileUpDataMC", "PileUpPtBB", "PileUpPtEC1", "PileUpPtRef",
                              "RelativeBal", "RelativeJEREC1",
                              "RelativePtBB","RelativePtEC1",
                              "RelativeStatFSR", "RelativeStatEC",
                              "SinglePionECAL", "SinglePionHCAL",
                              "TimePtEta", "FlavorQCD",
                              "Total"
                             };
};

#endif  // PARTICLESELECTOR_HH
