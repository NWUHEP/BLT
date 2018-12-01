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

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

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
    bool BTagModifier(const baconhep::TJet* jet, string, string, float) const;
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

    // For getting b tag scale factors and their uncertainties
    BTagCalibration *btagCalibrator, *mistagCalibrator;
    BTagCalibrationReader *btagReader, *mistagReader;
    vector<string> _btagUncSources {
                                    "up", "down",
                                    "up_bfragmentation", "up_btempcorr", "up_cb",
                                    "up_cfragmentation", "up_dmux", "up_gluonsplitting",
                                    "up_jes", "up_jetaway", "up_ksl", "up_l2c", 
                                    "up_ltothers", "up_mudr",
                                    "up_mupt", "up_sampledependence", "up_pileup",
                                    "up_ptrel", "up_statistic",
                                    "down_bfragmentation", "down_btempcorr", "down_cb",
                                    "down_cfragmentation", "down_dmux", "down_gluonsplitting",
                                    "down_jes", "down_jetaway", "down_ksl", "down_l2c", 
                                    "down_ltothers", "down_mudr", 
                                    "down_mupt", "down_sampledependence", "down_pileup", 
                                    "down_ptrel", "down_statistic"
                                   };
};

#endif  // PARTICLESELECTOR_HH
