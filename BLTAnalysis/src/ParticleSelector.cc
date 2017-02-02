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

    _rng = new TRandom3(1337);
}

bool ParticleSelector::PassMuonID(const baconhep::TMuon* mu, const Cuts::muIDCuts& cutLevel) const {
    bool muPass = false;
    if (cutLevel.cutName == "tightMuID") {
        if (this->_parameters.period == "2012") {
            if (
                    mu->muNchi2       < cutLevel.NormalizedChi2
                    && mu->nValidHits > cutLevel.NumberOfValidMuonHits
                    && mu->nMatchStn  > cutLevel.NumberOfMatchedStations
                    && mu->nPixHits   > cutLevel.NumberOfValidPixelHits
                    && mu->nTkLayers  > cutLevel.TrackLayersWithMeasurement
                    && fabs(mu->d0)   < cutLevel.dxy
                    && fabs(mu->dz)   < cutLevel.dz
                    && test_bits(mu->typeBits, baconhep::kPFMuon) == cutLevel.IsPF
                    && test_bits(mu->typeBits, baconhep::kGlobal) == cutLevel.IsGLB
               ) muPass = true;

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

    }     return isoPass;
}

bool ParticleSelector::PassMuonIso(const baconhep::TMuon* mu, const Cuts::muDetIsoCuts& cutLevel) const {
    bool isoPass = false;
    if (
            mu->trkIso/mu->pt < cutLevel.trkIso
            && mu->hcalIso/mu->pt < cutLevel.hcalIso
            && mu->ecalIso/mu->pt < cutLevel.ecalIso
       ) 
        isoPass = true;
    return isoPass;
}

bool ParticleSelector::PassElectronID(const baconhep::TElectron* el, const Cuts::elIDCuts& cutLevel) const 
{
    bool elPass = false;
    float energyInverse = fabs(1. - el->eoverp)/el->ecalEnergy;

    if (fabs(el->scEta) < 1.479) 
    { // barrel
        if (
                fabs(el->dEtaIn)    < cutLevel.dEtaIn[0]
                && fabs(el->dPhiIn) < cutLevel.dPhiIn[0]
                && el->sieie        < cutLevel.sigmaIetaIeta[0]
                && el->hovere       < cutLevel.HadOverEm[0]
                && fabs(el->d0)     < cutLevel.dxy[0]
                && fabs(el->dz)     < cutLevel.dz[0]
                && energyInverse   < cutLevel.fabsEPDiff[0]
                && el->nMissingHits <= cutLevel.ConversionMissHits[0]
                && !el->isConv
           ) elPass = true;
    } 
    else if (fabs(el->scEta) > 1.4446 && fabs(el->scEta) < 1.566) 
    { // transition
        return elPass;
    } 
    else if (fabs(el->scEta) > 1.566) 
    { // endcap
        if (
                fabs(el->dEtaIn)    < cutLevel.dEtaIn[1]
                && fabs(el->dPhiIn) < cutLevel.dPhiIn[1]
                && el->sieie        < cutLevel.sigmaIetaIeta[1]
                && el->hovere       < cutLevel.HadOverEm[1]
                && fabs(el->d0)     < cutLevel.dxy[1]
                && fabs(el->dz)     < cutLevel.dz[1]
                && energyInverse   < cutLevel.fabsEPDiff[1]
                && el->nMissingHits <= cutLevel.ConversionMissHits[1]
                && !el->isConv    
           ) elPass = true;
    }    
    return elPass;
}

bool ParticleSelector::PassElectronMVA(const baconhep::TElectron* el, const Cuts::elMVACuts& cutLevel) const {
    bool elPass = false;

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

bool ParticleSelector::PassElectronIso(const baconhep::TElectron* el, const Cuts::elIsoCuts& cutLevel, float EAEl[7]) const 
{
    int iEta = 0;
    float etaBins[8] = {0., 1., 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};
    for (unsigned i = 0; i < 8; ++i) {
        if (fabs(el->scEta) > etaBins[i] && fabs(el->scEta) < etaBins[i+1]) {
            iEta = i;
            break;
        }
    }

    float combIso = el->chHadIso
                    + std::max(0.,(double)el->neuHadIso 
                    + el->gammaIso 
                    - _rhoFactor*EAEl[iEta]);

    bool isoPass = false;
    if (cutLevel.cutName == "mediumElIso") {
        if (el->pt < 20) {
            if (combIso/el->pt < 0.10) isoPass = true;
        } else {
            if (combIso/el->pt < 0.15) isoPass = true;
        }

    } else {
        if (combIso/el->pt < cutLevel.relCombIso) isoPass = true;
    }
    return isoPass;
}

bool ParticleSelector::PassPhotonID(const baconhep::TPhoton* ph, const Cuts::phIDCuts& cutLevel) const {
    bool phoPass = false;
    bool phoPass1 = false;
    bool phoPass2 = false;

    //if (fabs(ph->scEta) > 2.5) return phoPass;  // uncomment to apply eta requirement
    if (cutLevel.cutName == "preSelPhID") {
        if (fabs(ph->scEta) > 1.4442 && fabs(ph->scEta) < 1.566) return phoPass;
        if (fabs(ph->scEta) < 1.479) {
            if (
                    (!ph->isConv) == cutLevel.PassedEleSafeVeto[0] 
                    && ph->sieie < cutLevel.sigmaIetaIeta[0]
               ) {
                if (ph->r9 > 0.9) {
                    if (ph->hovere < cutLevel.HadOverEm[0]) phoPass1 = true;
                } else {
                    if (ph->hovere < cutLevel.HadOverEm[1]) phoPass1 = true;
                }
            }
        } else {
            if (
                    (!ph->isConv) == cutLevel.PassedEleSafeVeto[1]
                    && ph->sieie  < cutLevel.sigmaIetaIeta[1]
                    && ph->hovere < cutLevel.HadOverEm[1]
               ) phoPass1 = true;
        }

        if (ph->r9 > 0.9) {
            if (
                    //ph->hcalIso - 0.005*ph->pt     < cutLevel.HcalIso[0]  //FIXME
                    //&& ph->trkIso - 0.002*ph->pt   < cutLevel.TrkIso[0]   //FIXME
                    //&& ph->cicPF4chgpfIso02          < cutLevel.ChPfIso[0]  //FIXME
                    1
               ) phoPass2 = true;
        } else {
            if (
                    //ph->hcalIso - 0.005*ph->pt     < cutLevel.HcalIso[1]  //FIXME
                    //&& ph->trkIso - 0.002*ph->pt   < cutLevel.TrkIso[1]   //FIXME
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
    return phoPass;
}

bool ParticleSelector::PassPhotonIso(const baconhep::TPhoton* ph, const Cuts::phIsoCuts& cutLevel, float EAPho[7][3]) const {
    bool isoPass = false;

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
                 && max((double)chIsoCor,0.)      < cutLevel.chIso[0]
                 && max((double)nhIsoCor,0.)      < cutLevel.nhIso[0] + 0.04*ph->pt
                 && max((double)phIsoCor,0.)      < cutLevel.phIso[0] + 0.005*ph->pt
                ) || (
                    fabs(ph->scEta) > 1.566
                    && max((double)chIsoCor,0.)      < cutLevel.chIso[1]
                    && max((double)nhIsoCor,0.)      < cutLevel.nhIso[1] + 0.04*ph->pt
                    )
           ) isoPass = true;
    } else {
        if (
                (
                 fabs(ph->scEta) < 1.4442
                 && max((double)chIsoCor,0.)      < cutLevel.chIso[0]
                 && max((double)nhIsoCor,0.)      < cutLevel.nhIso[0] + 0.04*ph->pt
                 && max((double)phIsoCor,0.)      < cutLevel.phIso[0] + 0.005*ph->pt
                ) || (
                    fabs(ph->scEta) > 1.566
                    && max((double)chIsoCor,0.)      < cutLevel.chIso[1]
                    && max((double)nhIsoCor,0.)      < cutLevel.nhIso[1] + 0.04*ph->pt
                    && max((double)phIsoCor,0.)      < cutLevel.phIso[1] + 0.005*ph->pt
                    )
           ) isoPass = true;
    }
    return isoPass;
}


bool ParticleSelector::PassJetID(const baconhep::TJet* jet, const Cuts::jetIDCuts& cutLevel) const {
    bool jetPass = false;
    if (fabs(jet->eta) <= 2.7) {
        if (
                jet->neuHadFrac       < 0.99
                && jet->neuEmFrac     < 0.99
                && jet->nParticles    > 1
           ) {
            if (fabs(jet->eta) <= 2.4) {
                if (jet->chHadFrac > 0 && jet->nCharged > 0 && jet->chEmFrac < 0.99) 
                    jetPass = true;
            } else {
                jetPass = true;
            }
        }
    } else if (fabs(jet->eta) <= 3.0) { 
        if (jet->neuEmFrac > 0.01 && jet->neuHadFrac < 0.98 && jet->nNeutrals > 2) 
            jetPass = true;
    } else {
        if (jet->neuEmFrac < 0.9 && jet->nNeutrals > 10) 
            jetPass = true;
    }

    return jetPass;
}

bool ParticleSelector::PassJetPUID(const baconhep::TJet* jet) const {
    bool pass = false;
    if (0 <= fabs(jet->eta) && fabs(jet->eta) < 2.5) {
        if(jet->mva >= -0.89) 
            pass = true;
    } else {
        pass = true;
    }
    
    /*else if(2.5 <= fabs(jet->eta) && fabs(jet->eta) < 2.75) {
        if (jet->mva >= -0.52) 
            pass = true;
    } else if(2.75 <= fabs(jet->eta) && fabs(jet->eta) < 3) {
        if (jet->mva >= -0.38) 
            pass = true;
    } else if(3 <= fabs(jet->eta) && fabs(jet->eta) < 4.7) {
        if (jet->mva >= -0.30) 
            pass = true;
    }*/

    return pass;
}
bool ParticleSelector::BTagModifier(const baconhep::TJet* jet, string tagName) const
{
    bool  isBTagged = false;
    float jetPt     = jet->pt;
    //float jetEta    = jet->eta;
    float bTag      = jet->csv;
    int   jetFlavor = jet->hadronFlavor;

    // Get b tag efficiency and mistag scale factor
    float btagSF = 1.;
    float mistagSF = 1.;
    if (tagName == "CSVT") {
        // These SF are provided by the b tag POG 
        if (bTag > 0.935) isBTagged = true;
        btagSF   = 0.857294 + 3.75846e-05*jetPt; 
        mistagSF = 0.688619 + 260.84/(jetPt*jetPt); 
    } else if (tagName == "MVAT") {
        // These SF are provided by the b tag POG 
        if (bTag > 0.9432) isBTagged = true;
        btagSF   = 0.52032*((1.+(0.370141*jetPt))/(1.+(0.196566*jetPt))); 
        mistagSF = 0.985864+122.8/(jetPt*jetPt)+0.000416939*jetPt;
    } else if (tagName == "MVAM") {
        // These SF are provided by the b tag POG 
        if (bTag > 0.4432) isBTagged = true;
        btagSF   = 0.885562 - 5.67668e-05*jetPt; 
        mistagSF = 1.07108+-0.0015866*jetPt+3.06322e-06*jetPt*jetPt+-1.7438e-09*jetPt*jetPt*jetPt;
    }


    // Upgrade or downgrade jet
    float rNumber = _rng->Uniform(1.);
    if (abs(jetFlavor) == 5 || abs(jetFlavor) == 4) {
        float mcEff = 1.;
        if (abs(jetFlavor) == 4) 
            mcEff = 0.3;//_ctagEff->Eval(jetPt);
        else if (abs(jetFlavor) == 5) 
            mcEff = 0.6;//_btagEff->Eval(jetPt);

        if(btagSF >= 1){  // use this if SF>1
            if (!isBTagged) { //upgrade to tagged
                float mistagRate = (1. - btagSF) / (1. - 1./mcEff);
                if (rNumber < mistagRate) isBTagged = true;
            }
        } else if (btagSF < 1) { //downgrade tagged to untagged
            if(isBTagged && rNumber > btagSF) isBTagged = false;
        }

    } else {
        float mcEff = 0.003; //_misTagEff->Eval(jetPt);
        if(mistagSF > 1){  // use this if SF>1
            if (!isBTagged) { //upgrade to tagged
                float mistagPercent = (1. - mistagSF) / (1. - (1./mcEff));
                if (rNumber < mistagPercent) isBTagged = true;
            }
        } else if (mistagSF < 1) { //downgrade tagged to untagged
            if (isBTagged && rNumber > mistagSF) isBTagged = false;
        }
    }

    return isBTagged;
}
