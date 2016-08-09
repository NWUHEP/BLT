#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"

#include <iostream>

//FIXME: implement these

using namespace baconhep;
using namespace std;

bool test_bits(unsigned int bits, unsigned int test) {
    return (bits & test) == test;
}

ParticleSelector::ParticleSelector(const Parameters& parameters, const Cuts& cuts) {
    this->_parameters = parameters;
    this->_cuts = cuts;
}

bool ParticleSelector::PassMuonID(const baconhep::TMuon* mu, const Cuts::muIDCuts& cutLevel) const {
    bool muPass = false;
    if (cutLevel.cutName == "tightMuID") {
        if (this->_parameters.period == "2012") {
            if (
                    mu->muNchi2    < cutLevel.NormalizedChi2
                    && mu->nValidHits > cutLevel.NumberOfValidMuonHits
                    && mu->nMatchStn  > cutLevel.NumberOfMatchedStations
                    && mu->nPixHits   > cutLevel.NumberOfValidPixelHits
                    && mu->nTkLayers  > cutLevel.TrackLayersWithMeasurement
                    && fabs(mu->d0)   < cutLevel.dxy
                    && fabs(mu->dz)   < cutLevel.dz
                    && test_bits(mu->typeBits, baconhep::kPFMuon) == cutLevel.IsPF
                    && test_bits(mu->typeBits, baconhep::kGlobal) == cutLevel.IsGLB
               ) muPass = true;

        } else {
            muPass = test_bits(mu->pogIDBits, baconhep::kPOGTightMuon);
        }
    }
    return muPass;
}

bool ParticleSelector::PassMuonIso(const baconhep::TMuon* mu, const Cuts::muIsoCuts& cutLevel) const {
    bool isoPass = false;
    if (cutLevel.cutName == "tightMuIso" || cutLevel.cutName == "looseMuIso") {
        float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
        if (combIso/mu->pt < cutLevel.relCombIso04) 
            isoPass = true;

    } else if (cutLevel.cutName == "amumuDetIso") {
        if (mu->trkIso/mu->pt > 0.1) 
            isoPass = true;
    }
    return isoPass;
}

bool ParticleSelector::PassElectronID(const baconhep::TElectron* el, const Cuts::elIDCuts& cutLevel) const {
    bool elPass = false;
    std::cout << warning() << "BLT electron ID is missing certain variables." << std::endl;

    //if (fabs(el->scEta) > 2.5) return elPass;  // uncomment to apply eta requirement
    if (cutLevel.cutName == "mvaPreElID"){
        if (fabs(el->scEta) < 1.479) {
            if (
                el->sieie                            < cutLevel.sigmaIetaIeta[0]
                && el->hovere                        < cutLevel.HadOverEm[0]
                && el->trkIso/el->pt                 < cutLevel.dr03TkSumPt[0]
                && el->ecalIso/el->pt                < cutLevel.dr03EcalRecHitSumEt[0]
                && el->hcalIso/el->pt                < cutLevel.dr03HcalTowerSumEt[0]
                //&& el->nMissingLostHits             == cutLevel.numberOfLostHits[0]  //FIXME
            ) elPass = true;
        } else {
            if (
                el->sieie                            < cutLevel.sigmaIetaIeta[1]
                && el->hovere                        < cutLevel.HadOverEm[1]
                && el->trkIso/el->pt                 < cutLevel.dr03TkSumPt[1]
                && el->ecalIso/el->pt                < cutLevel.dr03EcalRecHitSumEt[1]
                && el->hcalIso/el->pt                < cutLevel.dr03HcalTowerSumEt[1]
                //&& el->nMissingLostHits             == cutLevel.numberOfLostHits[1]  //FIXME
            ) elPass = true;
        }

    } else {
        if (fabs(el->scEta) > 1.4446 && fabs(el->scEta) < 1.566)  return elPass;
        if (
            (
                fabs(el->scEta) < 1.479
                && fabs(el->dEtaIn)                         < cutLevel.dEtaIn[0]
                && fabs(el->dPhiIn)                         < cutLevel.dPhiIn[0]
                && el->sieie                                < cutLevel.sigmaIetaIeta[0]
                && el->hovere                               < cutLevel.HadOverEm[0]
                && fabs(el->d0)                             < cutLevel.dxy[0]
                && fabs(el->dz)                             < cutLevel.dz[0]
                //&& el->inverseEnergyMomentumDiff            < cutLevel.fabsEPDiff[0]  //FIXME
                && el->nMissingHits                        <= cutLevel.ConversionMissHits[0]
                && (!el->isConv)                           == cutLevel.PassedConversionProb[0]
            ) || (
                fabs(el->scEta) > 1.479
                && fabs(el->dEtaIn)                         < cutLevel.dEtaIn[1]
                && fabs(el->dPhiIn)                         < cutLevel.dPhiIn[1]
                && el->sieie                                < cutLevel.sigmaIetaIeta[1]
                && el->hovere                               < cutLevel.HadOverEm[1]
                && fabs(el->d0)                             < cutLevel.dxy[1]
                && fabs(el->dz)                             < cutLevel.dz[1]
                //&& el->inverseEnergyMomentumDiff            < cutLevel.fabsEPDiff[1]  //FIXME
                && el->nMissingHits                        <= cutLevel.ConversionMissHits[1]
                && (!el->isConv)                           == cutLevel.PassedConversionProb[1]
            )
        ) elPass = true;
    }
    return elPass;
}

bool ParticleSelector::PassElectronMVA(const baconhep::TElectron* el, const Cuts::elMVACuts& cutLevel) const {
    bool elPass = false;
    std::cout << warning() << "BLT electron MVA ID is missing certain variables." << std::endl;

    if (cutLevel.cutName == "hzgMVAID") {
        if (el->pt > cutLevel.pt[0] && el->pt < cutLevel.pt[1]) {
            //if (el->mvaOld > cutLevel.mvaVal[0]) elPass = true;  //FIXME
        } else if (el->pt > cutLevel.pt[1]) {
            //if (el->mvaOld > cutLevel.mvaVal[1]) elPass = true;  //FIXME
        }

    } else if (cutLevel.cutName == "hzzMVAID") {
        if (el->pt > cutLevel.pt[0] && el->pt < cutLevel.pt[1]) {
            if (fabs(el->eta) < cutLevel.eta[0]) {
                if (
                    //(el->mvaHZZ                 > cutLevel.mvaVal[0])  //FIXME
                    //&& (el->nMissingHits        < cutLevel.missHits[0])  //FIXME
                    1
                    && (el->sip3d               < cutLevel.sip[0])
                ) elPass = true;
            } else if (fabs(el->eta) < cutLevel.eta[1]) {
                if (
                    //(el->mvaHZZ                 > cutLevel.mvaVal[1])  //FIXME
                    //&& (el->nMissingHits        < cutLevel.missHits[1])  //FIXME
                    1
                    && (el->sip3d               < cutLevel.sip[1])
                ) elPass = true;
            } else if (fabs(el->eta) < cutLevel.eta[2]) {
                if (
                    //(el->mvaHZZ                 > cutLevel.mvaVal[2])  //FIXME
                    //&& (el->nMissingHits        < cutLevel.missHits[2])  //FIXME
                    1
                    && (el->sip3d               < cutLevel.sip[2])
                ) elPass = true;
            }
        } else if (el->pt > cutLevel.pt[1]) {
            if (fabs(el->eta) < cutLevel.eta[0]) {
                if (
                    //(el->mvaHZZ                 > cutLevel.mvaVal[3])  //FIXME
                    //&& (el->nMissingHits        < cutLevel.missHits[3])  //FIXME
                    1
                    && (el->sip3d               < cutLevel.sip[3])
                ) elPass = true;
            } else if (fabs(el->eta) < cutLevel.eta[1]) {
                if (
                    //(el->mvaHZZ                 > cutLevel.mvaVal[4])  //FIXME
                    //&& (el->nMissingHits        < cutLevel.missHits[4])  //FIXME
                    1
                    && (el->sip3d               < cutLevel.sip[4])
                ) elPass = true;
            } else if (fabs(el->eta) < cutLevel.eta[2]) {
                if (
                    //(el->mvaHZZ                 > cutLevel.mvaVal[5])  //FIXME
                    //&& (el->nMissingHits        < cutLevel.missHits[5])  //FIXME
                    1
                    && (el->sip3d               < cutLevel.sip[5])
                ) elPass = true;
            }
        }

    } else if (cutLevel.cutName == "hwwMVAID") {
        if (el->pt > cutLevel.pt[0] && el->pt < cutLevel.pt[1]) {
            if (fabs(el->eta) < cutLevel.eta[0]) {
                if (
                    //(el->mvaHZZ                 > cutLevel.mvaVal[0])  //FIXME
                    //&& (el->nMissingHits       == cutLevel.missHits[0])  //FIXME
                    1
                    && ((!el->isConv)          == cutLevel.conversionVeto[0])
                ) elPass = true;
            } else if (fabs(el->eta) < cutLevel.eta[1]) {
                if (
                    //(el->mvaHZZ                 > cutLevel.mvaVal[1])  //FIXME
                    //&& (el->nMissingHits       == cutLevel.missHits[1])  //FIXME
                    1
                    && ((!el->isConv)          == cutLevel.conversionVeto[1])
                ) elPass = true;
            } else if (fabs(el->eta) < cutLevel.eta[2]) {
                if (
                    //(el->mvaHZZ                 > cutLevel.mvaVal[2])  //FIXME
                    //&& (el->nMissingHits       == cutLevel.missHits[2])  //FIXME
                    1
                    && ((!el->isConv)          == cutLevel.conversionVeto[2])
                ) elPass = true;
            }
        } else if (el->pt > cutLevel.pt[1]) {
            if (fabs(el->eta) < cutLevel.eta[0]) {
                if (
                    //(el->mvaHZZ                  > cutLevel.mvaVal[3])  //FIXME
                    //&& (el->nMissingHits        == cutLevel.missHits[3])  //FIXME
                    1
                    && ((!el->isConv)           == cutLevel.conversionVeto[3])
                ) elPass = true;
            } else if (fabs(el->eta) < cutLevel.eta[1]) {
                if (
                    //(el->mvaHZZ                  > cutLevel.mvaVal[4])  //FIXME
                    //&& (el->nMissingHits        == cutLevel.missHits[4])  //FIXME
                    1
                    && ((!el->isConv)           == cutLevel.conversionVeto[4])
                ) elPass = true;
            } else if (fabs(el->eta) < cutLevel.eta[2]) {
                if (
                    //(el->mvaHZZ                  > cutLevel.mvaVal[5])  //FIXME
                    //&& (el->nMissingHits        == cutLevel.missHits[5])  //FIXME
                    1
                    && ((!el->isConv)           == cutLevel.conversionVeto[5])
                ) elPass = true;
            }
        }
    }
    return elPass;
}

bool ParticleSelector::PassElectronIso(const baconhep::TElectron* el, const Cuts::elIsoCuts& cutLevel) const {
    bool isoPass = false;
    std::cout << warning() << "BLT electron iso is missing certain variables." << std::endl;

    float effArea = 0;
    //float effArea = el->effArea;  //FIXME

    float combIso = (el->chHadIso + std::max(0.,(double)el->neuHadIso + el->gammaIso - _rhoFactor*effArea));

    if (cutLevel.cutName == "mediumElIso") {
        if (el->pt < 20) {
            if (combIso/el->pt < 0.10) isoPass = true;
        } else {
            if (combIso/el->pt < 0.15) isoPass = true;
        }

    } else {
        if (combIso/el->pt < cutLevel.relCombIso04) isoPass = true;
    }
    return isoPass;
}

bool ParticleSelector::PassPhotonID(const baconhep::TPhoton* ph, const Cuts::phIDCuts& cutLevel) const {
    bool phoPass = false;
    bool phoPass1 = false;
    bool phoPass2 = false;
    std::cout << warning() << "BLT photon ID is missing certain variables." << std::endl;

    //if (fabs(ph->scEta) > 2.5) return phoPass;  // uncomment to apply eta requirement
    if (cutLevel.cutName == "preSelPhID") {
        if (fabs(ph->scEta) > 1.4442 && fabs(ph->scEta) < 1.566) return phoPass;
        if (fabs(ph->scEta) < 1.479) {
            if (
                (!ph->isConv)                   == cutLevel.PassedEleSafeVeto[0]
                && ph->sieie                     < cutLevel.sigmaIetaIeta[0]
            ) {
                if (ph->r9 > 0.9) {
                    if (ph->hovere               < cutLevel.HadOverEm[0]) phoPass1 = true;
                } else {
                    if (ph->hovere               < cutLevel.HadOverEm[1]) phoPass1 = true;
                }
            }
        } else {
            if (
                (!ph->isConv)                   == cutLevel.PassedEleSafeVeto[1]
                && ph->sieie                     < cutLevel.sigmaIetaIeta[1]
                && ph->hovere                    < cutLevel.HadOverEm[1]
            ) phoPass1 = true;
        }

        if (ph->r9 > 0.9) {
            if (
                //ph->hcalIso03 - 0.005*ph->pt     < cutLevel.HcalIso[0]  //FIXME
                //&& ph->trkIso03 - 0.002*ph->pt   < cutLevel.TrkIso[0]   //FIXME
                //&& ph->cicPF4chgpfIso02          < cutLevel.ChPfIso[0]  //FIXME
                1
            ) phoPass2 = true;
        } else {
            if (
                //ph->hcalIso03 - 0.005*ph->pt     < cutLevel.HcalIso[1]  //FIXME
                //&& ph->trkIso03 - 0.002*ph->pt   < cutLevel.TrkIso[1]   //FIXME
                //&& ph->cicPF4chgpfIso02          < cutLevel.ChPfIso[1]  //FIXME
                1
            ) phoPass2 = true;
        }
        if (phoPass1 && phoPass2) phoPass = true;

    } else if (cutLevel.cutName == "medPhID"){
        if (fabs(ph->scEta) > 1.4442 && fabs(ph->scEta) < 1.566) return phoPass;
        if (
            (
                fabs(ph->scEta)  < 1.4442
                && (!ph->isConv)                == cutLevel.PassedEleSafeVeto[0]
                && ph->hovere                    < cutLevel.HadOverEm[0]
                && ph->sieie                     < cutLevel.sigmaIetaIeta[0]
            ) || (
                fabs(ph->scEta)  > 1.566
                && (!ph->isConv)                == cutLevel.PassedEleSafeVeto[1]
                && ph->hovere                    < cutLevel.HadOverEm[1]
                && ph->sieie                     < cutLevel.sigmaIetaIeta[1]
            )
        ) phoPass = true;
    }
    return phoPass;
}

bool ParticleSelector::PassPhotonMVA(const baconhep::TPhoton* ph, const Cuts::phMVACuts& cutLevel) const {
    bool phoPass = false;
    std::cout << warning() << "BLT photon MVA ID is missing certain variables." << std::endl;

    //FIXME

    return phoPass;
}

bool ParticleSelector::PassPhotonIso(const baconhep::TPhoton* ph, const Cuts::phIsoCuts& cutLevel, float EAPho[7][3]) const {
    bool isoPass = false;
    std::cout << warning() << "BLT photon iso is missing certain variables." << std::endl;

    //if (fabs(ph->scEta) > 2.5) return isoPass;  // uncomment to apply eta requirement

    float chEA, nhEA, phEA, chIsoCor, nhIsoCor, phIsoCor;
    if (fabs(ph->scEta) < 1.0) {
        chEA = EAPho[0][0];
        nhEA = EAPho[0][1];
        phEA = EAPho[0][2];
    } else if (fabs(ph->scEta) < 1.479) {
        chEA = EAPho[1][0];
        nhEA = EAPho[1][1];
        phEA = EAPho[1][2];
    } else if (fabs(ph->scEta) < 2.0) {
        chEA = EAPho[2][0];
        nhEA = EAPho[2][1];
        phEA = EAPho[2][2];
    } else if (fabs(ph->scEta) < 2.2) {
        chEA = EAPho[3][0];
        nhEA = EAPho[3][1];
        phEA = EAPho[3][2];
    } else if (fabs(ph->scEta) < 2.3) {
        chEA = EAPho[4][0];
        nhEA = EAPho[4][1];
        phEA = EAPho[4][2];
    } else if (fabs(ph->scEta) < 2.4) {
        chEA = EAPho[5][0];
        nhEA = EAPho[5][1];
        phEA = EAPho[5][2];
    } else {
        chEA = EAPho[6][0];
        nhEA = EAPho[6][1];
        phEA = EAPho[6][2];
    }

    chIsoCor = ph->chHadIso - _rhoFactor*chEA;
    nhIsoCor = ph->neuHadIso - _rhoFactor*nhEA;
    phIsoCor = ph->gammaIso -_rhoFactor*phEA;

    if (cutLevel.cutName == "loosePhIso"){
        if (
            (
                fabs(ph->scEta) < 1.4442
                && max((double)chIsoCor,0.)      < cutLevel.chIso03[0]
                && max((double)nhIsoCor,0.)      < cutLevel.nhIso03[0] + 0.04*ph->pt
                && max((double)phIsoCor,0.)      < cutLevel.phIso03[0] + 0.005*ph->pt
            ) || (
                fabs(ph->scEta) > 1.566
                && max((double)chIsoCor,0.)      < cutLevel.chIso03[1]
                && max((double)nhIsoCor,0.)      < cutLevel.nhIso03[1] + 0.04*ph->pt
            )
        ) isoPass = true;
    } else {
        if (
            (
                fabs(ph->scEta) < 1.4442
                && max((double)chIsoCor,0.)      < cutLevel.chIso03[0]
                && max((double)nhIsoCor,0.)      < cutLevel.nhIso03[0] + 0.04*ph->pt
                && max((double)phIsoCor,0.)      < cutLevel.phIso03[0] + 0.005*ph->pt
            ) || (
                fabs(ph->scEta) > 1.566
                && max((double)chIsoCor,0.)      < cutLevel.chIso03[1]
                && max((double)nhIsoCor,0.)      < cutLevel.nhIso03[1] + 0.04*ph->pt
                && max((double)phIsoCor,0.)      < cutLevel.phIso03[1] + 0.005*ph->pt
            )
        ) isoPass = true;
    }
    return isoPass;
}

bool ParticleSelector::PassVBFJetID(const baconhep::TJet* jet, const Cuts::vbfJetIDCuts& cutLevel) const {
    bool jetPass = false;
    std::cout << warning() << "BLT VBF jet ID is missing certain variables." << std::endl;

    //if (fabs(jet->eta) > 4.7) return isoPass;  // uncomment to apply eta requirement
    if (cutLevel.cutName == "vbfJetID"){
        if (fabs(jet->eta) < 2.5) {
            if (
                jet->betaStar/log(_npv-0.64)     < cutLevel.betaStarC[0]
                && jet->dR2Mean                  < cutLevel.dR2Mean[0]
            ) jetPass = true;
        } else if (fabs(jet->eta) < 2.75) {
            if(
                jet->betaStar/log(_npv-0.64)     < cutLevel.betaStarC[1]
                && jet->dR2Mean                  < cutLevel.dR2Mean[1]
            ) jetPass = true;
        } else if (fabs(jet->eta) < 3.0) {
            if(
                jet->dR2Mean                     < cutLevel.dR2Mean[2]
            ) jetPass = true;
        } else {
            if(
                jet->dR2Mean                     < cutLevel.dR2Mean[3]
            ) jetPass = true;
        }
    }
    return jetPass;
}

bool ParticleSelector::PassJetID(const baconhep::TJet* jet, const Cuts::jetIDCuts& cutLevel) const {
    bool jetPass = false;
    //std::cout << warning() << "BLT jet ID is not doing cross cleaning." << std::endl;

    //if (fabs(jet->eta) > 4.7) return isoPass;  // uncomment to apply eta requirement
    if (cutLevel.cutName == "looseJetID") {
        if (fabs(jet->eta) <= 2.4) {
            if (
                jet->neuHadFrac       < cutLevel.NHF
                && jet->neuEmFrac     < cutLevel.NEMF
                && jet->nParticles    > cutLevel.NumConst
                && jet->chHadFrac     > cutLevel.CHF
                && jet->nCharged      > cutLevel.CHM
                && jet->chEmFrac      < cutLevel.CEMF
                && jet->csv           > cutLevel.CSV
            ) jetPass = true;
        } else {
            if (
                jet->neuHadFrac       < cutLevel.NHF
                && jet->neuEmFrac     < cutLevel.NEMF
                && jet->nParticles    > cutLevel.NumConst
            ) jetPass = true;
        }

    } else {
        if (fabs(jet->eta) <= 2.4) {
            if (
                jet->neuHadFrac       < cutLevel.NHF
                && jet->neuEmFrac     < cutLevel.NEMF
                && jet->nParticles    > cutLevel.NumConst
                && jet->chHadFrac     > cutLevel.CHF
                && jet->nCharged      > cutLevel.CHM
                && jet->chEmFrac      < cutLevel.CEMF
                && jet->csv           > cutLevel.CSV
            ) jetPass = true;
        }
    }
    return jetPass;
}

bool ParticleSelector::FindGoodDiMuons(const std::vector<baconhep::TMuon*>& muons,
                         TLorentzVector& mu1, TLorentzVector& mu2, int& index1, int& index2) const {
    bool goodZ = false;
    TLorentzVector tmpZ;
    index1 = -1;
    index2 = -1;

    for (unsigned i=0; i<muons.size(); i++) {
        const TMuon* imuon = muons.at(i);
        if (imuon->pt > _cuts.leadMuPt) {
            for (unsigned j=i+1; j<muons.size(); j++) {
                const TMuon* jmuon = muons.at(j);
                if (jmuon->pt > _cuts.trailMuPt && imuon->q != jmuon->q) {
                    copy_p4(imuon, MUON_MASS, mu1);
                    copy_p4(jmuon, MUON_MASS, mu2);
                    index1 = i;
                    index2 = j;
                    tmpZ = mu1 + mu2;

                    if (0<tmpZ.M() && tmpZ.M()<999999)  // Not yet imposed Z mass cut
                        goodZ = true;
                }
                if (goodZ)  break;
            }
        }
        if (goodZ) break;
    }

    if (!goodZ) {
        mu1 = TLorentzVector();
        mu2 = TLorentzVector();
        index1 = -1;
        index2 = -1;
    }
    return goodZ;
}

bool ParticleSelector::FindGoodDiElectrons(const std::vector<baconhep::TElectron*>& electrons,
                             TLorentzVector& el1, TLorentzVector& el2, int& index1, int& index2) const {
    bool goodZ = false;
    TLorentzVector tmpZ;
    index1 = -1;
    index2 = -1;

    for (unsigned i=0; i<electrons.size(); i++) {
        const TElectron* ielectron = electrons.at(i);
        if (ielectron->pt > _cuts.leadElPt) {
            for (unsigned j=i+1; j<electrons.size(); j++) {
                const TElectron* jelectron = electrons.at(j);
                if (jelectron->pt > _cuts.trailElPt && ielectron->q != jelectron->q) {
                    copy_p4(ielectron, ELE_MASS, el1);
                    copy_p4(jelectron, ELE_MASS, el2);
                    index1 = i;
                    index2 = j;
                    tmpZ = el1 + el2;

                    if (0<tmpZ.M() && tmpZ.M()<999999)  // Not yet imposed Z mass cut
                        goodZ = true;
                }
                if (goodZ)  break;
            }
        }
        if (goodZ) break;
    }

    if (!goodZ) {
        el1 = TLorentzVector();
        el2 = TLorentzVector();
        index1 = -1;
        index2 = -1;
    }
    return goodZ;
}

bool ParticleSelector::FindGenZToLL(const TClonesArray* genParticles,
                      TLorentzVector& genZ, TLorentzVector& genLep1, TLorentzVector& genLep2,
                      int& indexZ, int& indexLep1, int& indexLep2) const {
    bool goodZ = false;
    indexZ = -1;
    indexLep1 = -1;
    indexLep2 = -1;

    for (int i=0; i<genParticles->GetEntries(); i++) {
        const TGenParticle* particle = (TGenParticle*) genParticles->At(i);
        assert(particle);

        if (std::abs(particle->pdgId) == Z_PDGID && particle->status == 22) {
            goodZ = true;
            indexZ = i;
            copy_p4(particle, particle->mass, genZ);
        }
    }
    if (!goodZ) {
        indexZ = -1;
        genZ = TLorentzVector();
    }

    bool goodZToLL = false;
    for (int i=0; i<genParticles->GetEntries(); i++) {
        const TGenParticle* particle = (TGenParticle*) genParticles->At(i);
        assert(particle);

        if (
            (std::abs(particle->pdgId) == ELE_PDGID || std::abs(particle->pdgId) == MUON_PDGID)  // not including tau
            && particle->status == 1
        ) {
            // Check parent
            const TGenParticle* particleParent = 0;
            if (particle->parent >= 0) {
                particleParent = (TGenParticle*) genParticles->At(particle->parent);
                while (particleParent->pdgId == particle->pdgId && particleParent->parent >= 0) {
                    particleParent = (TGenParticle*) genParticles->At(particleParent->parent);
                }
            }

            if (!particleParent || particleParent->pdgId != Z_PDGID)
                continue;

            if (indexLep1 == -1) {
                indexLep1 = i;
                copy_p4(particle, particle->mass, genLep1);
            } else if (indexLep2 == -1) {
                indexLep2 = i;
                copy_p4(particle, particle->mass, genLep2);
                goodZToLL = true;
            } else {
                // More than 2 leptons
                goodZToLL = false;
            }
        }
    }

    if (!goodZToLL) {
        indexLep1 = -1;
        indexLep2 = -1;
        genLep1 = TLorentzVector();
        genLep2 = TLorentzVector();
    }

    return goodZToLL;
}
