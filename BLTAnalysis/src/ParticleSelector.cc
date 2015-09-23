#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"

//FIXME: implement these

ParticleSelector::ParticleSelector(const Parameters& parameters, const Cuts& cuts) {

}

bool ParticleSelector::PassMuonID(const baconhep::TMuon* mu, const Cuts::muIDCuts& cutLevel) const {
    return 1;
}

bool ParticleSelector::PassMuonIso(const baconhep::TMuon* mu, const Cuts::muIsoCuts& cutLevel) const {
    return 1;
}

bool ParticleSelector::PassElectronID(const baconhep::TElectron* el, const Cuts::elIDCuts& cutLevel) const {
    return 1;
}

bool ParticleSelector::PassElectronMVA(const baconhep::TElectron* el, const Cuts::elMVACuts& cutLevel) const {
    return 1;
}

bool ParticleSelector::PassElectronIso(const baconhep::TElectron* el, const Cuts::elIsoCuts& cutLevel) const {
    return 1;
}

bool ParticleSelector::PassPhotonID(const baconhep::TPhoton* ph, const Cuts::phIDCuts& cutLevel) const {
    return 1;
}

bool ParticleSelector::PassPhotonMVA(const baconhep::TPhoton* ph, const Cuts::phMVACuts& cutLevel) const {
    return 1;
}

bool ParticleSelector::PassPhotonIso(const baconhep::TPhoton* ph, const Cuts::phIsoCuts& cutLevel) const {
    return 1;
}

bool ParticleSelector::PassJetID(const baconhep::TJet* jet, const Cuts::jetIDCuts& cutLevel) const {
    return 1;
}

bool ParticleSelector::FindGoodDiMuons(const std::vector<baconhep::TMuon*>& muons, int& index1, int& index2) const {
    return 1;
}

bool ParticleSelector::FindGoodDiElectrons(const std::vector<baconhep::TElectron*>& electrons, int& index1, int& index2) const {
    return 1;
}
