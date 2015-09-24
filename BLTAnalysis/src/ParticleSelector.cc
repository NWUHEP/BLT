#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"

#include <iostream>

//FIXME: implement these

using namespace baconhep;
using namespace std;

bool test_bits(unsigned int bits, unsigned int test) {
    return !!(bits & test);
}

ParticleSelector::ParticleSelector(const Parameters& parameters, const Cuts& cuts) {

}

bool ParticleSelector::PassMuonID(const baconhep::TMuon* mu, const Cuts::muIDCuts& cutLevel) const {
    bool muPass = false;

    if (cutLevel.cutName == "tightMuID") {
        if (
            //fabs(mu->eta) < 2.4  // uncomment to apply eta requirement
            fabs(mu->eta) < 99.
            && test_bits(mu->typeBits, baconhep::kPFMuon) == cutLevel.IsPF
            && test_bits(mu->typeBits, baconhep::kGlobal) == cutLevel.IsGLB
            && mu->muNchi2    < cutLevel.NormalizedChi2
            && mu->nValidHits > cutLevel.NumberOfValidMuonHits
            && mu->nMatchStn  > cutLevel.NumberOfMatchedStations
            && mu->nPixHits   > cutLevel.NumberOfValidPixelHits
            && mu->nTkLayers  > cutLevel.TrackLayersWithMeasurement
            && fabs(mu->d0)   < cutLevel.dxy
            && fabs(mu->dz)   < cutLevel.dz
            ) muPass = true;

        if (test_bits(mu->pogIDBits, baconhep::kPOGTightMuon) != muPass)
            std::cout << warning() << "BLT tightMuID does not agree with the CMSSW version" << std::endl;
    }
    return muPass;
}

bool ParticleSelector::PassMuonIso(const baconhep::TMuon* mu, const Cuts::muIsoCuts& cutLevel) const {
    bool isoPass = false;
    if (cutLevel.cutName == "tightMuIso" || cutLevel.cutName == "looseMuIso") {
        float combIso = (mu->chHadIso + max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
        if (combIso/mu->pt < cutLevel.relCombIso04) isoPass = true;
    }
    return isoPass;
}

bool ParticleSelector::PassElectronID(const baconhep::TElectron* el, const Cuts::elIDCuts& cutLevel) const {
    bool elPass = false;
    //if (fabs(el->scEta) > 2.5) return elPass;  // uncomment to apply eta requirement
    if (cutLevel.cutName == "mvaPreElID"){
        if (fabs(el->scEta) < 1.479) {
            if (
                el->sieie                            < cutLevel.sigmaIetaIeta[0]
                //&& el->hadronicOverEm                < cutLevel.HadOverEm[0]  //FIXME
                && el->trkIso/el->pt                 < cutLevel.dr03TkSumPt[0]
                && el->ecalIso/el->pt                < cutLevel.dr03EcalRecHitSumEt[0]
                && el->hcalIso/el->pt                < cutLevel.dr03HcalTowerSumEt[0]
                //&& el->gsfNumberOfLostHits          == cutLevel.numberOfLostHits[0]  //FIXME
            ) elPass = true;
        } else {
            if (
                el->sieie                            < cutLevel.sigmaIetaIeta[1]
                //&& el->hadronicOverEm                < cutLevel.HadOverEm[1]  //FIXME
                && el->trkIso/el->pt                 < cutLevel.dr03TkSumPt[1]
                && el->ecalIso/el->pt                < cutLevel.dr03EcalRecHitSumEt[1]
                && el->hcalIso/el->pt                < cutLevel.dr03HcalTowerSumEt[1]
                //&& el->gsfNumberOfLostHits          == cutLevel.numberOfLostHits[1]  //FIXME
            ) elPass = true;
        }

    } else {
        if (fabs(el->scEta) > 1.4446 && fabs(el->scEta) < 1.566)  return elPass;
        if (
            (fabs(el->scEta) < 1.479
            //&& fabs(el->deltaEtaSeedClusterTrackAtCalo) < cutLevel.dEtaIn[0]  //FIXME
            //&& fabs(el->deltaPhiSeedClusterTrackAtCalo) < cutLevel.dPhiIn[0]  //FIXME
            && el->sieie                                < cutLevel.sigmaIetaIeta[0]
            //&& el->hadronicOverEm                       < cutLevel.HadOverEm[0]  //FIXME
            && fabs(el->d0)                             < cutLevel.dxy[0]
            && fabs(el->dz)                             < cutLevel.dz[0]
            //&& el->inverseEnergyMomentumDiff            < cutLevel.fabsEPDiff[0]  //FIXME
            //&& el->conversionMissHits                  <= cutLevel.ConversionMissHits[0]  //FIXME
            && (!el->isConv)                           == cutLevel.PassedConversionProb[0]
            ) ||
            (fabs(el->scEta) > 1.479
            //&& fabs(el->deltaEtaSeedClusterTrackAtCalo) < cutLevel.dEtaIn[1]  //FIXME
            //&& fabs(el->deltaPhiSeedClusterTrackAtCalo) < cutLevel.dPhiIn[1]  //FIXME
            && el->sieie                                < cutLevel.sigmaIetaIeta[1]
            //&& el->hadronicOverEm                       < cutLevel.HadOverEm[1]  //FIXME
            && fabs(el->d0)                             < cutLevel.dxy[1]
            && fabs(el->dz)                             < cutLevel.dz[1]
            //&& el->inverseEnergyMomentumDiff            < cutLevel.fabsEPDiff[1]  //FIXME
            //&& el->ConversionMissHits                  <= cutLevel.ConversionMissHits[1]  //FIXME
            && (!el->isConv)                           == cutLevel.PassedConversionProb[1]
            )
        ) elPass = true;
    }
    return elPass;
}

bool ParticleSelector::PassElectronMVA(const baconhep::TElectron* el, const Cuts::elMVACuts& cutLevel) const {
    bool elPass = false;
    //FIXME
    return elPass;
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

void ParticleSelector::FindGoodDiMuons(const std::vector<baconhep::TMuon*>& muons,
                         TLorentzVector& mu1, TLorentzVector& mu2, int& index1, int& index2) const {

}

void ParticleSelector::FindGoodDiElectrons(const std::vector<baconhep::TElectron*>& electrons,
                             TLorentzVector& el1, TLorentzVector& el2, int& index1, int& index2) const {

}

void ParticleSelector::FindGenDYToLL(const std::vector<baconhep::TGenParticle*>& genParticles,
                       TLorentzVector& genZ, TLorentzVector& genLep1, TLorentzVector& genLep2,
                       int& indexZ, int& indexLep1, int& indexLep2) const {

}
