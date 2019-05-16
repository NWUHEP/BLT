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

    const std::string _cmssw_base = getenv("CMSSW_BASE");
    _rng = new TRandom3();

    std::string rcPath = _cmssw_base + "/src/BLT/BLTAnalysis/data/roccor/roccor_" + parameters.period + ".txt";
    muonCorr = new RoccoR(rcPath);

    // offline jet corrections on-the-fly
    vector<string> ds = split(parameters.datasetgroup, '_');
    string datasetname = ds[0];
    const std::string cmssw_base = getenv("CMSSW_BASE");
    if (datasetname == "muon" || datasetname == "electron") { // data 
        // jet energy corrections
        string runPeriod   = ds[1];
        cout << datasetname << " " << runPeriod << endl;
       
        std::string jecPath;

        if (parameters.period == "2016") {
            if (runPeriod == "2016B" || runPeriod == "2016C" || runPeriod == "2016D") {
                runPeriod = "BCD";
            } else if (runPeriod == "2016E" || runPeriod == "2016F") {
                runPeriod = "EF";
            } else if (runPeriod == "2016G") {
                runPeriod = "G";
            } else if (runPeriod == "2016H") {
                runPeriod = "H";
            } 
            jecPath = cmssw_base + "/src/BLT/BLTAnalysis/data/jec/Summer16_23Sep2016" + runPeriod + "V4_DATA/Summer16_23Sep2016" + runPeriod + "V4_DATA";
        }

        else if (parameters.period == "2017") {
            if (runPeriod == "2017B") {
                runPeriod = "B";
            } else if (runPeriod == "2017C") {
                runPeriod = "C";
            } else if (runPeriod == "2017D" || runPeriod == "2017E") {
                runPeriod = "DE";
            } else if (runPeriod == "2017F") {
                runPeriod = "F";
            } 
            jecPath = cmssw_base + "/src/BLT/BLTAnalysis/data/jec/Fall17_17Nov2017" + runPeriod + "_V32_DATA/Fall17_17Nov2017" + runPeriod + "_V32_DATA";
        }
        
        else if (parameters.period == "2018") {
            if (runPeriod == "2018A") {
                runPeriod = "A";
            } else if (runPeriod == "2018B") {
                runPeriod = "B";
            } else if (runPeriod == "2018C") {
                runPeriod = "C";
            } else if (runPeriod == "2018D") {
                runPeriod = "D";
            } 
            jecPath = cmssw_base + "/src/BLT/BLTAnalysis/data/jec/Autumn18_Run" + runPeriod + "_V8_DATA/Autumn18_Run" + runPeriod + "_V8_DATA";
        }

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
        JetCorrectorParameters *jecUnc  = new JetCorrectorParameters(jecPath + "_UncertaintySources_AK4PFchs.txt", "Total");
        _jecUncertainty = new JetCorrectionUncertainty(*jecUnc);

    } else { // MC 
        // jet energy corrections
        std::string jecPath;
        if (parameters.period == "2016") {
            jecPath = cmssw_base + "/src/BLT/BLTAnalysis/data/jec/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC";
        } 
        else if (parameters.period == "2017") {
            jecPath = cmssw_base + "/src/BLT/BLTAnalysis/data/jec/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC";
        } 
        else if (parameters.period == "2018") {
            jecPath = cmssw_base + "/src/BLT/BLTAnalysis/data/jec/Autumn18_V8_MC/Autumn18_V8_MC";
        } 

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
        JetCorrectorParameters *jecUnc  = new JetCorrectorParameters(jecPath + "_UncertaintySources_AK4PFchs.txt", "Total");
        _jecUncertainty = new JetCorrectionUncertainty(*jecUnc);

        // jet energy resolution
        jetResolution   = JME::JetResolution(cmssw_base + "/src/BLT/BLTAnalysis/data/jer/jet_pt_resolution.dat");
        jetResolutionSF = JME::JetResolutionScaleFactor(cmssw_base + "/src/BLT/BLTAnalysis/data/jer/jet_resolution_scale_factors.dat");
    }
}

bool ParticleSelector::PassMuonID(const baconhep::TMuon* mu, string idName) const
{
    bool muPass = false;

    if (idName == "tight") {
        if (
               (mu->typeBits & baconhep::kPFMuon) 
               && (mu->typeBits & baconhep::kGlobal) 
               && mu->muNchi2    < 10.
               && mu->nMatchStn  > 1
               && mu->nPixHits   > 0
               && fabs(mu->d0)   < 0.2
               && fabs(mu->dz)   < 0.5
               && mu->nTkLayers  > 5 
               && mu->nValidHits > 0
           ) {
            muPass = true;
        }
    }

    else if (idName == "HZZ") {
        if (
            fabs(mu->d0) < 0.5
            && fabs(mu->dz) < 1.0
            && (((mu->typeBits & baconhep::kGlobal) || 
               ((mu->typeBits & baconhep::kTracker) && mu->nMatchStn > 0)) &&
               (mu->btt != 2)) // Global muon or (arbitrated) tracker muon
            && fabs(mu->sip3d) < 4.0                 
           ) {
                // We now have h->ZZ->4l "loose" muons
                if (mu->pt < 200.0) {
                    if (mu->pogIDBits & baconhep::kPOGLooseMuon)
                        muPass = true;
                }
                else {
                    // need to pass the tracker high-pt ID
                    if (
                            (mu->pogIDBits & baconhep::kPOGLooseMuon) ||
                            (mu->nMatchStn > 1
                            && (mu->ptErr/mu->pt) < 0.3 
                            && fabs(mu->d0) < 0.2
                            && fabs(mu->dz) < 0.5
                            && mu->nPixHits > 0
                            && mu->nTkLayers > 5)
                       )
                        muPass = true;
                }
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
                && fabs(el->d0)     < 1.
                && fabs(el->dz)     < 1.
                && !el->isConv
           ) elPass = true;
    } else if (fabs(el->scEta) > 1.4446 && fabs(el->scEta) < 1.566) { // transition
        elPass = false;
    } else if (fabs(el->scEta) > 1.566) { // endcap
        if (
                el->sieie           < 0.0292
                && fabs(el->dEtaIn) < 0.00605
                && fabs(el->dPhiIn) < 0.0394
                && el->hovere       < 0.0641
                && energyInverse    < 0.0129
                && el->nMissingHits <= 1
                && fabs(el->d0)     < 1.
                && fabs(el->dz)     < 1.
                && !el->isConv    
           ) elPass = true;
    }    
    return elPass;
}

bool ParticleSelector::PassElectronMVA(const baconhep::TElectron* el, string idName) const {
    bool elPass = false;

    if (idName == "loose") 
    {
        if (fabs(el->eta) < 0.8) {
            if (el->mvaSpring16 > 0.837) 
                elPass = true;
        }
        else if (fabs(el->eta) >= 0.8 && fabs(el->eta) < 1.479) {
            if (el->mvaSpring16 > 0.715)
                elPass = true;
        }
        else {
            if (el->mvaSpring16 > 0.357)
                elPass = true;
        }
    } 
    
    else if (idName == "tight") 
    {
        if (fabs(el->eta) < 0.8) {
            if (el->mvaSpring16 > 0.941) 
                elPass = true;
        }
        else if (fabs(el->eta) >= 0.8 && fabs(el->eta) < 1.479) {
            if (el->mvaSpring16 > 0.899)
                elPass = true;
        }
        else {
            if (el->mvaSpring16 > 0.758)
                elPass = true;
        }
    } 
    
    else if (idName == "HZZ") 
    {
        if (_parameters.period == "2016") {
            if (el->pt > 5. && el->pt < 10.) {
                if (fabs(el->scEta) < 0.8) {
                    if (el->mvaSpring16HZZ > -0.211)
                        elPass = true;
                } else if (fabs(el->scEta) < 1.479) {
                    if (el->mvaSpring16HZZ > -0.396)  
                        elPass = true;
                } else if (fabs(el->scEta) < 2.5) {
                    if (el->mvaSpring16HZZ > -0.215)  
                       elPass = true;
                }
            } else if (el->pt > 10.) {
                if (fabs(el->scEta) < 0.8) {
                    if (el->mvaSpring16HZZ > -0.870)
                        elPass = true;
                } else if (fabs(el->scEta) < 1.479) {
                    if (el->mvaSpring16HZZ > -0.838)  
                        elPass = true;
                } else if (fabs(el->scEta) < 2.5) {
                    if (el->mvaSpring16HZZ > -0.763)
                        elPass = true;
                }
            }
        }
        else if (_parameters.period == "2017") {
            if (el->pt > 5. && el->pt < 10.) {
                if (fabs(el->scEta) < 0.8) {
                    if (el->mvaFall17V1Iso > -0.1)
                        elPass = true;
                } else if (fabs(el->scEta) < 1.479) {
                    if (el->mvaFall17V1Iso > -0.28)  
                        elPass = true;
                } else if (fabs(el->scEta) < 2.5) {
                    if (el->mvaFall17V1Iso > -0.05)  
                       elPass = true;
                }
            } else if (el->pt > 10.) {
                if (fabs(el->scEta) < 0.8) {
                    if (el->mvaFall17V1Iso > -0.83)
                        elPass = true;
                } else if (fabs(el->scEta) < 1.479) {
                    if (el->mvaFall17V1Iso > -0.77)  
                        elPass = true;
                } else if (fabs(el->scEta) < 2.5) {
                    if (el->mvaFall17V1Iso > -0.69)
                        elPass = true;
                }
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
                    + std::max(0.,(double)el->neuHadIso 
                            + el->gammaIso 
                            - _rhoFactor*effArea[iEta]);

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

bool ParticleSelector::PassPhotonMVA(const baconhep::TPhoton* ph, string idName) const 
{
    bool phoPass = false;
    if (_parameters.period == "2016") {
        if (idName == "loose") {
            if (ph->mvaSpring16 > 0.2) {
                phoPass = true;
            }
        }
        else if (idName == "tight") {
            if (fabs(ph->scEta) <= 1.479) { // barrel
                if (ph->mvaSpring16 > 0.68) {
                    phoPass = true;
                }
            }
            else { //endcap
                if (ph->mvaSpring16 > 0.60) {
                    phoPass = true;
                }
            }
        }
    }
    else { // 2017 or 2018 
        if (idName == "loose") {
            if (fabs(ph->scEta) <= 1.479) { // barrel
                if (ph->mvaFall17V2 > -0.02) {
                    phoPass = true;
                }
            }
            else { //endcap
                if (ph->mvaFall17V2 > -0.26) {
                    phoPass = true;
                }
            }
        }
        else if (idName == "tight") {
            if (fabs(ph->scEta) <= 1.479) { // barrel
                if (ph->mvaFall17V2 > 0.42) {
                    phoPass = true;
                }
            }
            else { //endcap
                if (ph->mvaFall17V2 > 0.14) {
                    phoPass = true;
                }
            }
        }
    }
    
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

bool ParticleSelector::BTagModifier(const baconhep::TJet* jet, string tagName, int btagSyst, int mistagSyst) const
{
    bool  isBTagged = false;
    float jetPt     = jet->pt;
    int   jetFlavor = jet->hadronFlavor;
    float bTag      = -1;

    float rNumber = _rng->Uniform(1.);

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
    if (tagName == "CSVT") {
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
    } else if (tagName == "MVAM") {
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
                float mistagRate = (1. - btagSF) / (1. - 1./mcEff);
                if (rNumber < mistagRate) 
                    isBTagged = true;
            }
        } else if (btagSF < 1) { // downgrade b tagged to untagged
            if (isBTagged && rNumber > btagSF) 
                isBTagged = false;
        }
    } else {
        //cout << mistagSF << " " << mistagSyst << " " << isBTagged << " ";
        if (mistagSF > 1) {  // use this if SF>1
            if (!isBTagged) { //upgrade to tagged
                float mistagRate = (1. - mistagSF) / (1. - 1./mcEff);
                if (rNumber < mistagRate) 
                    isBTagged = true;
            }
        } else if (mistagSF < 1) { //downgrade tagged to untagged
            if (isBTagged && rNumber > mistagSF) 
                isBTagged = false;
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

double ParticleSelector::JetUncertainty(const baconhep::TJet* jet) const
{
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

float ParticleSelector::GetMuonIsolation(const baconhep::TMuon* mu) const
{
    //float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    float combIso = (mu->chHadIso03 + std::max(0.,(double)mu->neuHadIso03 + mu->gammaIso03 - 0.5*mu->puIso03));
    return combIso;
}

float ParticleSelector::GetElectronIsolation(const baconhep::TElectron* el, const float rho) const
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

    float combIso = el->chHadIso + std::max(0., (double)el->neuHadIso + el->gammaIso - rho*effArea[iEta]);

    return combIso;
}

float ParticleSelector::GetPhotonIsolation(const baconhep::TPhoton* pho, const float rho) const
{
    int iEta = 0;
    float etaBins[8] = {0., 1., 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};
    float effArea[8] = {0.1703, 0.1715, 0.1213, 0.1230, 0.1635, 0.1937, 0.2393};
    for (unsigned i = 0; i < 8; ++i) {
        if (fabs(pho->scEta) > etaBins[i] && fabs(pho->scEta) < etaBins[i+1]) {
            iEta = i;
            break;
        }
    }

    float combIso = pho->chHadIso + std::max(0., (double)pho->neuHadIso + pho->gammaIso - rho*effArea[iEta]);

    return combIso;
}

void ParticleSelector::ApplyMuonMomentumCorrection(baconhep::TMuon* mu, bool isData) const
{
    double muonSF = 1.;
    if (isData)
        muonSF = muonCorr->kScaleDT(mu->q, mu->pt, mu->eta, mu->phi, 0, 0);
    else
        muonSF = muonCorr->kSmearMC(mu->q, mu->pt, mu->eta, mu->phi, mu->nTkLayers, _rng->Rndm(), 0, 0);
    mu->pt = muonSF*mu->pt;
}

void ParticleSelector::ApplyTauEnergyScaleCorrection(baconhep::TTau* tau) const
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Tau_energy_scale 
{
    if (tau->decaymode == 0) {
        tau->pt *= 0.995;
    } else if (tau->decaymode == 1) {
        tau->pt *= 1.01;
    } else if (tau->decaymode == 10) {
        tau->pt *= 1.006;
    }
}

