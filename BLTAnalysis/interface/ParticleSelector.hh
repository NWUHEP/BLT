#ifndef PARTICLESELECTOR_HH
#define PARTICLESELECTOR_HH

#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"

#include "BLT/BLTAnalysis/interface/Parameters.hh"
#include "BLT/BLTAnalysis/interface/Cuts.hh"

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMath.h>


class ParticleSelector {
public:
    ParticleSelector(const Parameters& parameters, const Cuts& cuts);
    ~ParticleSelector() {}

    // Setters
    void SetRealData(bool isRealData)       { _isRealData = isRealData; }
    void SetPV(const TVector3& pv)          { _pv = pv; }
    void SetRho(float rho)                  { _rho = rho; }
    void SetRunNumber(unsigned int run)     { _runNumber = run; }
    void SetLumiSection(unsigned int lumi)  { _lumiSection = lumi; }
    void SetEventNumber(unsigned int event) { _evtNumber = event; }

    // Identifiers
    bool PassMuonID(const baconhep::TMuon* mu, const Cuts::muIDCuts& cutLevel) const;
    bool PassMuonIso(const baconhep::TMuon* mu, const Cuts::muIsoCuts& cutLevel) const;
    bool PassElectronID(const baconhep::TElectron* el, const Cuts::elIDCuts& cutLevel) const;
    bool PassElectronMVA(const baconhep::TElectron* el, const Cuts::elMVACuts& cutLevel) const;
    bool PassElectronIso(const baconhep::TElectron* el, const Cuts::elIsoCuts& cutLevel) const;
    bool PassPhotonID(const baconhep::TPhoton* ph, const Cuts::phIDCuts& cutLevel) const;
    bool PassPhotonMVA(const baconhep::TPhoton* ph, const Cuts::phMVACuts& cutLevel) const;
    bool PassPhotonIso(const baconhep::TPhoton* ph, const Cuts::phIsoCuts& cutLevel) const;
    bool PassJetID(const baconhep::TJet* jet, const Cuts::jetIDCuts& cutLevel) const;

    // Finders
    bool FindGoodDiMuons(const std::vector<baconhep::TMuon*>& muons, int& index1, int& index2) const;
    bool FindGoodDiElectrons(const std::vector<baconhep::TElectron*>& electrons, int& index1, int& index2) const;


private:
    Parameters      _parameters;
    Cuts            _cuts;
    bool            _isRealData;
    unsigned int    _runNumber;
    unsigned int    _lumiSection;
    unsigned int    _evtNumber;
    TVector3        _pv;
    float           _rho;
};

#endif  // PARTICLESELECTOR_HH
