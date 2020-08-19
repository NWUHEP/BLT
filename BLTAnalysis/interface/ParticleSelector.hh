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
#include "BLT/BLTAnalysis/interface/RoccoR.h"

#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TRandom3.h"

#include <string>
#include <vector>
#include <memory>
#include <cassert>
#include <map>

//CMSSW libraries
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

using namespace std;

class ParticleSelector {
public:
    ParticleSelector(const Parameters& parameters);
    ~ParticleSelector() {}

    // Setters
    void SetRho(float rhoFactor)            { _rhoFactor = rhoFactor; }

    // Muons
    bool PassMuonID(const baconhep::TMuon* mu, string) const;
    float GetMuonIsolation(const baconhep::TMuon* mu) const;
    void ApplyMuonMomentumCorrection(baconhep::TMuon* mu, bool isData) const;

    // Electrons
    bool PassElectronMVA(const baconhep::TElectron* el, string) const;
    float GetElectronIsolation(const baconhep::TElectron* el, const float rho) const;
    pair<float, float> GetElectronScaleErr(const baconhep::TElectron* el) const;

    // Photons
    bool PassPhotonMVA(const baconhep::TPhoton* ph, string) const;
    float GetPhotonIsolation(const baconhep::TPhoton* pho, const float rho) const;
    //float GetPhotonWorstChargedIsolation(baconhep::TPhoton* pho, TClonesArray* pfCands, TClonesArray* vertices, baconhep::TVertex* pv) const; // on hold

    // Taus
    void ApplyTauEnergyScaleCorrection(baconhep::TTau* tau) const;



    // Jets
    bool PassJetID(const baconhep::TJet* jet, string) const;
    bool PassJetPUID(const baconhep::TJet* jet) const;
    bool BTagModifier(const baconhep::TJet* jet, string, int, int) const;
    double JetCorrector(const baconhep::TJet* jet, string) const;
    std::vector<float> GetJetSubcorrections(const baconhep::TJet* jet, string) const;
    double JetUncertainty(const baconhep::TJet* jet) const;
    //pair<float, float> JetResolutionAndSF(const baconhep::TJet* jet, int) const;
    pair<float, float> JetResolutionAndSF(const baconhep::TJet* jet, int);

    vector<float> GetCorrectedPhotonMVA(const baconhep::TPhoton* ph, bool isData) const;

private:
    Parameters _parameters;
    float      _rhoFactor;
    TRandom3*  _rng;

    // For offline jet corrections
    FactorizedJetCorrector* _jetCorrector;
    JetCorrectionUncertainty* _jecUncertainty;
    //JME::JetResolution jetResolution;
    //JME::JetResolutionScaleFactor jetResolutionSF;
    std::map<std::string, JME::JetResolution> jetResolution;
    std::map<std::string, JME::JetResolutionScaleFactor> jetResolutionSF;

    // For offline Rochester muon corrections
    RoccoR *muonCorr;
};

#endif  // PARTICLESELECTOR_HH
