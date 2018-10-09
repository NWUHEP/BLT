#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"

#include <iostream>

using namespace baconhep;
using namespace std;

// Why the fuck is this not part of the standard library!?
vector<string> split(const string &s, char delim) {
    stringstream ss(s);
    string item;
    vector<string> tokens;
    while (getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}

bool test_bits(unsigned int bits, unsigned int test) {
    return (bits & test) == test;
}

ParticleSelector::ParticleSelector(const Parameters& parameters, const Cuts& cuts) {
    this->_parameters = parameters;
    this->_cuts = cuts;

    this->_rng = new TRandom3(1337);

    // offline jet corrections on-the-fly
    vector<string> ds = split(parameters.datasetgroup, '_');
    string datasetname = ds[0];
    const std::string cmssw_base = getenv("CMSSW_BASE");
    if (datasetname == "muon" || datasetname == "electron") { // data 
        // jet energy corrections
        string runPeriod   = ds[1];
        cout << datasetname << " " << runPeriod << endl;
        if (runPeriod == "2016B" || runPeriod == "2016C" || runPeriod == "2016D") {
            runPeriod = "BCD";
        } else if (runPeriod == "2016E" || runPeriod == "2016F") {
            runPeriod = "EF";
        } else if (runPeriod == "2016G") {
            runPeriod = "G";
        } else if (runPeriod == "2016H") {
            runPeriod = "H";
        }

        std::string jecPath = cmssw_base + "/src/BLT/BLTAnalysis/data/Summer16_23Sep2016" + runPeriod + "V4_DATA/Summer16_23Sep2016" + runPeriod + "V4_DATA";
        std::cout << jecPath << std::endl;
        JetCorrectorParameters *ResJetPar = new JetCorrectorParameters(jecPath + "_L2L3Residual_AK4PFchs.txt"); 
        JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(jecPath + "_L3Absolute_AK4PFchs.txt");
        JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(jecPath + "_L2Relative_AK4PFchs.txt");
        JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(jecPath + "_L1FastJet_AK4PFchs.txt");

        vector<JetCorrectorParameters> vPar;
        vPar.push_back(*L1JetPar);
        vPar.push_back(*L2JetPar);
        vPar.push_back(*L3JetPar);
        vPar.push_back(*ResJetPar);

        _jetCorrector = new FactorizedJetCorrector(vPar);

        // jec uncertainties
        //JetCorrectorParameters *jecUnc  = new JetCorrectorParameters(jecPath + "_UncertaintySources_AK4PFchs.txt", "Total");
        //_jecUncertainty = new JetCorrectionUncertainty(*jecUnc);

    } else { // MC 
        // jet energy corrections
        std::string jecPath = cmssw_base + "/src/BLT/BLTAnalysis/data/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC";
        std::cout << jecPath << std::endl;
        JetCorrectorParameters *ResJetPar = new JetCorrectorParameters(jecPath + "_L2L3Residual_AK4PFchs.txt"); 
        JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(jecPath + "_L3Absolute_AK4PFchs.txt");
        JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(jecPath + "_L2Relative_AK4PFchs.txt");
        JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(jecPath + "_L1FastJet_AK4PFchs.txt");

        vector<JetCorrectorParameters> vPar;
        vPar.push_back(*L1JetPar);
        vPar.push_back(*L2JetPar);
        vPar.push_back(*L3JetPar);
        vPar.push_back(*ResJetPar);

        _jetCorrector = new FactorizedJetCorrector(vPar);

        for (const auto& name: _jecNames) {
           JetCorrectorParameters *jecUnc  = new JetCorrectorParameters(jecPath + "_UncertaintySources_AK4PFchs.txt", name);
           _jecUncertaintyMap[name] = new JetCorrectionUncertainty(*jecUnc);
        } 


        // jet energy resolution
        jetResolution   = JME::JetResolution(cmssw_base + "/src/BLT/BLTAnalysis/data/jet_pt_resolution.dat");
        jetResolutionSF = JME::JetResolutionScaleFactor(cmssw_base + "/src/BLT/BLTAnalysis/data/jet_resolution_scale_factors.dat");
    }
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

    if (fabs(el->scEta) < 1.479) { // barrel
        if (
                el->sieie           < 0.00998
                && fabs(el->dEtaIn) < 0.00308
                && fabs(el->dPhiIn) < 0.0816
                && el->hovere       < 0.0414
                && energyInverse    < 0.0129
                && el->nMissingHits <= 1
                && fabs(el->d0)     < 0.05
                && fabs(el->dz)     < 0.1
                && !el->isConv
           ) elPass = true;
    } else {
        if (
                el->sieie           < 0.0292
                && fabs(el->dEtaIn) < 0.00605
                && fabs(el->dPhiIn) < 0.0394
                && el->hovere       < 0.0641
                && energyInverse    < 0.0129
                && el->nMissingHits <= 1
                && fabs(el->d0)     < 0.1
                && fabs(el->dz)     < 0.2
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
    float effArea[8] = {0.1703, 0.1715, 0.1213, 0.1230, 0.1635, 0.1937, 0.2393};
    for (unsigned i = 0; i < 8; ++i) {
        if (fabs(el->scEta) > etaBins[i] && fabs(el->scEta) < etaBins[i+1]) {
            iEta = i;
            break;
        }
    }

    float combIso = el->chHadIso 
                    + std::max(0., (double)el->neuHadIso + el->gammaIso - _rhoFactor*effArea[iEta]);

    bool isoPass = false;
    if (fabs(el->scEta) <= 1.479) {
        if (combIso/el->pt < 0.0588) isoPass = true;
    } else {
        if (combIso/el->pt < 0.0571) isoPass = true;
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

bool ParticleSelector::BTagModifier(const baconhep::TJet* jet, string tagName, int btagSyst, int mistagSyst, float rNumber) const
{
    bool  isBTagged = false;
    float jetPt     = jet->pt;
    int   jetFlavor = jet->hadronFlavor;
    float bTag      = -1;

    float binningPt[] = {30, 50, 70, 100, 140, 200, 300, 600};
    int ptBin = 0;
    for (int i = 0; i < 7; ++i) {
        if (jetPt > binningPt[i] && jetPt <= binningPt[i+1]) {
            ptBin = i;
            break;
        }
    }

    // Get b tag efficiency and mistag scale factor
    float btagSF   = 1.;
    float mistagSF = 1.;
    float mcEff  = 1.;
    if (tagName == "CSVT") { // not up-to-date
        bTag = jet->csv;
        // These SF are provided by the b tag POG 
        if (bTag > 0.935) 
            isBTagged = true;
        btagSF   = 0.857294 + 3.75846e-05*jetPt; 
        mistagSF = 0.688619 + 260.84/(jetPt*jetPt); 

        if (abs(jetFlavor) == 5) {
            float bEff[] = {0.41637905, 0.45007627, 0.47419147, 0.48388148, 0.4745329, 0.45031636, 0.40974969};
            mcEff = bEff[ptBin];
        } else if (abs(jetFlavor) == 4) {
            mcEff = 0.03;
        } else {
            mcEff = 0.002;
        }
    } else if (tagName == "MVAT") {
        // These SF are provided by the b tag POG 
        bTag = jet->bmva;
        if (bTag > 0.9432) 
            isBTagged = true;

        btagSF   = 0.517971*(1.+0.332528*jetPt) / (1. + 0.174914*jetPt); 
        mistagSF = 0.985864 + 122.8/(jetPt*jetPt) + 0.000416939*jetPt;

        if (abs(jetFlavor) == 5) {
            float bEff[] = {0.41637905, 0.45007627, 0.47419147, 0.48388148, 0.4745329, 0.45031636, 0.40974969};
            mcEff = bEff[ptBin];

            float scale[] = {0.0661, 0.0513, 0.0477, 0.0453, 0.0575, 0.0802, 0.3285}; 
            if (btagSyst == 1) {
                btagSF += scale[ptBin];
            } else if (btagSyst == -1) {
                btagSF -= scale[ptBin];
            }
        } else if (abs(jetFlavor) == 4) {
            mcEff = 0.03;
            float scale[] = {0.01889, 0.01466, 0.01362, 0.0129, 0.0164, 0.0229, 0.0939}; 
            if (btagSyst == 1) {
                btagSF += scale[ptBin];
            } else if (btagSyst == -1) {
                btagSF -= scale[ptBin];
            }
        } else {
            mcEff = 0.002;
            if (mistagSyst == 1) {
                mistagSF *= (1 + (0.253674 - 0.000127486*jetPt + 8.91567e-08*jetPt*jetPt));
            } else if (mistagSyst == -1) {
                mistagSF *= (1 - (0.253674 - 0.000127486*jetPt + 8.91567e-08*jetPt*jetPt));
            }
        }
    } else if (tagName == "MVAM") { // not up-to-date
        // These SF are provided by the b tag POG 
        bTag = jet->bmva;
        if (bTag > 0.4432) 
            isBTagged = true;
        btagSF   = 0.600657*((1.+(0.753343*jetPt))/(1.+(0.472587*jetPt))); 
        mistagSF = 1.11046 - 0.00042021*jetPt + 1.48012e-06*jetPt*jetPt - 8.44735e-10*jetPt*jetPt*jetPt;

        if (abs(jetFlavor) == 5) {
            float bEff[] = {0.41637905, 0.45007627, 0.47419147, 0.48388148, 0.4745329, 0.45031636, 0.40974969};
            mcEff = bEff[ptBin];
        } else if (abs(jetFlavor) == 4) {
            mcEff = 0.03;
        } else {
            mcEff = 0.002;
        }
    }

    // Upgrade or downgrade jet
    if (abs(jetFlavor) == 5 || abs(jetFlavor) == 4) {
        if (btagSF > 1) {  // use this if SF>1
            if (!isBTagged) { // upgrade to b tagged
                float tagRate = (1. - btagSF) / (1. - 1./mcEff);
                if (rNumber < tagRate) 
                    isBTagged = true;
            }
        } else if (btagSF < 1) { // downgrade b tagged to untagged
            if (isBTagged) {
                float tagRate = 1. - btagSF; 
                if (rNumber < tagRate) 
                    isBTagged = false;
            }
        }
    } else {
        //cout << mistagSF << " " << mistagSyst << " " << isBTagged << " ";
        if (mistagSF > 1) {  // use this if SF>1
            if (!isBTagged) { //upgrade to tagged
                float tagRate = (1. - mistagSF) / (1. - 1./mcEff);
                if (rNumber < tagRate) 
                    isBTagged = true;
            }
        } else if (mistagSF < 1) { //downgrade tagged to untagged
            if (isBTagged) {
                float tagRate = 1. - btagSF; 
                if (rNumber < tagRate) 
                    isBTagged = false;
            }
        }
        //cout << isBTagged << endl;
    }

    return isBTagged;
}

double ParticleSelector::JetCorrector(const baconhep::TJet* jet, string tagName) const
{
    _jetCorrector->setJetEta(jet->eta);
    _jetCorrector->setJetPt(jet->ptRaw);
    _jetCorrector->setJetA(jet->area);
    _jetCorrector->setRho(_rhoFactor); 

    double correction = _jetCorrector->getCorrection();

    return correction;
}

double ParticleSelector::JetUncertainty(const baconhep::TJet* jet, string srcName) const
{
    JetCorrectionUncertainty* _jecUncertainty = _jecUncertaintyMap.at(srcName);
    _jecUncertainty->setJetEta(jet->eta);
    _jecUncertainty->setJetPt(jet->ptRaw);

    double uncertainty = _jecUncertainty->getUncertainty(true);

    return uncertainty;
}

pair<float, float> ParticleSelector::JetResolutionAndSF(const baconhep::TJet* jet, int syst) const
{
    JME::JetParameters jetParams;
    jetParams.setJetPt(jet->pt);
    jetParams.setJetEta(jet->eta);
    jetParams.setRho(_rhoFactor);

    float jres = jetResolution.getResolution(jetParams);
    float jsf = 1.;
    if (syst == 1) {
        jsf  = jetResolutionSF.getScaleFactor(jetParams, Variation::UP);
    } else if (syst == -1) {
        jsf  = jetResolutionSF.getScaleFactor(jetParams, Variation::DOWN);
    } else {
        jsf  = jetResolutionSF.getScaleFactor(jetParams);
    }

    return make_pair(jres, jsf);
}
