#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"

#include <iostream>
#include <math.h>
#include "TMVA/Reader.h"
#include "TFormula.h"

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

ParticleSelector::ParticleSelector(const Parameters& parameters) {
    this->_parameters = parameters;

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
            } else if (runPeriod == "2016G" || runPeriod == "2016H") {
                runPeriod = "GH";
            } 
            //jecPath = cmssw_base + "/src/BLT/BLTAnalysis/data/jec/Summer16_23Sep2016" + runPeriod + "V4_DATA/Summer16_23Sep2016" + runPeriod + "V4_DATA";
            jecPath = cmssw_base + "/src/BLT/BLTAnalysis/data/jec/Summer16_07Aug2017" + runPeriod + "_V11_DATA/Summer16_07Aug2017" + runPeriod + "_V11_DATA";
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
            //jecPath = cmssw_base + "/src/BLT/BLTAnalysis/data/jec/Summer16_23Sep2016V4_MC/Summer16_23Sep2016V4_MC";
            jecPath = cmssw_base + "/src/BLT/BLTAnalysis/data/jec/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC";
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
        //jetResolution   = JME::JetResolution(cmssw_base + "/src/BLT/BLTAnalysis/data/jer/jet_pt_resolution.dat");
        //jetResolutionSF = JME::JetResolutionScaleFactor(cmssw_base + "/src/BLT/BLTAnalysis/data/jer/jet_resolution_scale_factors.dat");
        
        std::string jerPath;
        if (parameters.period == "2016") {
            jerPath = cmssw_base + "/src/BLT/BLTAnalysis/data/jer/Summer16_25nsV1";
            jetResolution["2016"] = JME::JetResolution(jerPath + "_MC_PtResolution_AK4PFchs.txt");
            jetResolutionSF["2016"] = JME::JetResolutionScaleFactor(jerPath + "_MC_SF_AK4PFchs.txt");
        }
        else if (parameters.period == "2017") {
            jerPath = cmssw_base + "/src/BLT/BLTAnalysis/data/jer/Fall17_V3";
            jetResolution["2017"] = JME::JetResolution(jerPath + "_MC_PtResolution_AK4PFchs.txt");
            jetResolutionSF["2017"] = JME::JetResolutionScaleFactor(jerPath + "_MC_SF_AK4PFchs.txt");
        }
        else if (parameters.period == "2018") {
            jerPath = cmssw_base + "/src/BLT/BLTAnalysis/data/jer/Autumn18_RunABC_V7b";
            jetResolution["2018ABC"] = JME::JetResolution(jerPath + "_MC_PtResolution_AK4PFchs.txt");
            jetResolutionSF["2018ABC"] = JME::JetResolutionScaleFactor(jerPath + "_MC_SF_AK4PFchs.txt");
            jerPath = cmssw_base + "/src/BLT/BLTAnalysis/data/jer/Autumn18_RunD_V7b";
            jetResolution["2018D"] = JME::JetResolution(jerPath + "_MC_PtResolution_AK4PFchs.txt");
            jetResolutionSF["2018D"] = JME::JetResolutionScaleFactor(jerPath + "_MC_SF_AK4PFchs.txt");
            //jerPath = cmssw_base + "/src/BLT/BLTAnalysis/data/jer/Autumn18_V7b";
            //float rNumber = _rng->Uniform(1.);
            //if (rNumber > 0.543367) {
                //jerPath = cmssw_base + "/src/BLT/BLTAnalysis/data/jer/Autumn18_RunABC_V7b";
            //}
            //else {
                //jerPath = cmssw_base + "/src/BLT/BLTAnalysis/data/jer/Autumn18_RunD_V7b";
            //}

        }

        //std::cout << jerPath << std::endl;
        
        //jetResolution   = JME::JetResolution(jerPath + "_MC_PtResolution_AK4PFchs.txt");
        //jetResolutionSF = JME::JetResolutionScaleFactor(jerPath + "_MC_SF_AK4PFchs.txt");



        //jetResolution   = JME::JetResolution(cmssw_base + "/src/BLT/BLTAnalysis/data/jer/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt");
        //jetResolutionSF = JME::JetResolutionScaleFactor(cmssw_base + "/src/BLT/BLTAnalysis/data/jer/Summer16_25nsV1_MC_SF_AK4PFchs.txt");
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

    else if (idName == "looseFall17V2") {
        if (el->pt > 5. && el->pt < 10.) {
            if (fabs(el->scEta) < 0.8) {
                //if (el->mvaFall17V2Iso > 0.604775)
                if (el->mvaFall17V2Iso > 0.700642584415)
                    elPass = true;
            } else if (fabs(el->scEta) < 1.479) {
                //if (el->mvaFall17V2Iso > 0.628743)  
                if (el->mvaFall17V2Iso > 0.739335420875)  
                    elPass = true;
            } else if (fabs(el->scEta) < 2.5) {
                //if (el->mvaFall17V2Iso > 0.896462)  
                if (el->mvaFall17V2Iso > 1.45390456109)  
                   elPass = true;
            }
        } else if (el->pt > 10.) {
            if (fabs(el->scEta) < 0.8) {
                //if (el->mvaFall17V2Iso > -0.145237)
                if (el->mvaFall17V2Iso > -0.146270871164)
                    elPass = true;
            } else if (fabs(el->scEta) < 1.479) {
                //if (el->mvaFall17V2Iso > -0.0315746)  
                if (el->mvaFall17V2Iso > -0.0315850882679)  
                    elPass = true;
            } else if (fabs(el->scEta) < 2.5) {
                //if (el->mvaFall17V2Iso > -0.032173)
                if (el->mvaFall17V2Iso > -0.0321841194737)
                    elPass = true;
            }
        }
    }

    else if (idName == "wp90Fall17V2") {
        double rawMVAValue = 0.5*log((1+el->mvaFall17V2Iso)/(1-el->mvaFall17V2Iso));
        if (el->calibPt > 5. && el->calibPt < 10.) {
            if (fabs(el->scEta) < 0.8) {
                if (rawMVAValue > (2.84704783417 - exp(-1*el->calibPt / 3.32529515837) * 9.38050947827))
                    elPass = true;
            } else if (fabs(el->scEta) < 1.479) {
                if (rawMVAValue > (2.03833922005 - exp(-1*el->calibPt / 1.93288758682) *  15.364588247))
                    elPass = true;
            } else if (fabs(el->scEta) < 2.5) {
                if (rawMVAValue > (1.82704158461 - exp(-1*el->calibPt / 1.89796754399) *  19.1236071158))
                   elPass = true;
            }
        } else if (el->calibPt > 10.) {
            if (fabs(el->scEta) < 0.8) {
                if (rawMVAValue > (6.12931925263 - exp(-1*el->calibPt / 13.281753835) *  8.71138432196))
                    elPass = true;
            } else if (fabs(el->scEta) < 1.479) {
                if (rawMVAValue > (5.26289004857 - exp(-1*el->calibPt / 13.2154971491) *  8.0997882835))
                    elPass = true;
            } else if (fabs(el->scEta) < 2.5) {
                if (rawMVAValue > (4.37338792902 - exp(-1*el->calibPt / 14.0776094696) *  8.48513324496))
                    elPass = true;
            }
        }
    }
    
    else if (idName == "wp80Fall17V2") {
        double rawMVAValue = 0.5*log((1+el->mvaFall17V2Iso)/(1-el->mvaFall17V2Iso));
        if (el->calibPt > 5. && el->calibPt < 10.) {
            if (fabs(el->scEta) < 0.8) {
                if (rawMVAValue > (3.53495358797 - exp(-1*el->calibPt / 3.07272325141) * 9.94262764352))
                    elPass = true;
            } else if (fabs(el->scEta) < 1.479) {
                if (rawMVAValue > (3.06015605623 - exp(-1*el->calibPt / 1.95572234114) *  14.3091184421))
                    elPass = true;
            } else if (fabs(el->scEta) < 2.5) {
                if (rawMVAValue > (3.02052519639 - exp(-1*el->calibPt / 1.59784164742) *  28.719380105))
                   elPass = true;
            }
        } else if (el->calibPt > 10.) {
            if (fabs(el->scEta) < 0.8) {
                if (rawMVAValue > (7.35752275071 - exp(-1*el->calibPt / 15.87907864) *  7.61288809226))
                    elPass = true;
            } else if (fabs(el->scEta) < 1.479) {
                if (rawMVAValue > (6.41811074032 - exp(-1*el->calibPt / 14.730562874) *  6.96387331587))
                    elPass = true;
            } else if (fabs(el->scEta) < 2.5) {
                if (rawMVAValue > (5.64936312428 - exp(-1*el->calibPt / 16.3664949747) *  7.19607610311))
                    elPass = true;
            }
        }
    }

    else if (idName == "tightFall17V1NoIso") {
        double rawMVAValue = 0.5*log((1+el->mvaFall17V1NoIso)/(1-el->mvaFall17V1NoIso));
        if (el->calibPt > 5. && el->calibPt < 10.) {
            if (fabs(el->scEta) < 0.8) {
                if (rawMVAValue > (0.9530240956555949 - exp(-1*el->calibPt / 2.7591425841003647) *  0.4669644718545271))
                    elPass = true;
            } else if (fabs(el->scEta) < 1.479) {
                if (rawMVAValue > (0.9336564763961019 - exp(-1*el->calibPt / 2.709276284272272) *  0.33512286599215946))
                    elPass = true;
            } else if (fabs(el->scEta) < 2.5) {
                if (rawMVAValue > (0.9313133688365339 - exp(-1*el->calibPt / 1.5821934800715558) *  3.8889462619659265))
                   elPass = true;
            }
        } else if (el->calibPt > 10.) {
            if (fabs(el->scEta) < 0.8) {
                if (rawMVAValue > (0.9825268564943458 - exp(-1*el->calibPt / 8.702601455860762) *  1.1974861596609097))
                    elPass = true;
            } else if (fabs(el->scEta) < 1.479) {
                if (rawMVAValue > (0.9727509457929913 - exp(-1*el->calibPt / 8.179525631018565) *  1.7111755094657688))
                    elPass = true;
            } else if (fabs(el->scEta) < 2.5) {
                if (rawMVAValue > (0.9562619539540145 - exp(-1*el->calibPt / 8.109845366281608) *  3.013927699126942))
                    elPass = true;
            }
        }
    }

    return elPass;
}

bool ParticleSelector::PassPhotonMVA(const baconhep::TPhoton* ph, string idName) const 
{
    bool phoPass = false;

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
    else if (idName == "looseMingyan") {
        if (fabs(ph->scEta) <= 1.479) { // barrel
            if (ph->mvaFall17V2 > -0.4) {
                phoPass = true;
            }
        }
        else { //endcap
            if (ph->mvaFall17V2 > -0.59) {
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
    
    return phoPass;
}

bool ParticleSelector::PassJetID(const baconhep::TJet* jet, string idName) const {
    bool jetPass = false;
    if (idName == "loose") {
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
    }
    else if (idName == "tight") {
        if (_parameters.period == "2017") {
            if (fabs(jet->eta) <= 2.7) {
                if (
                        jet->neuHadFrac       < 0.90
                        && jet->neuEmFrac     < 0.90
                        && jet->nParticles    > 1
                   ) {
                    if (fabs(jet->eta) <= 2.4) {
                        if (jet->chHadFrac > 0 && jet->nCharged > 0) 
                            jetPass = true;
                    } else {
                        jetPass = true;
                    }
                }
            } else if (fabs(jet->eta) <= 3.0) { 
                if (jet->neuEmFrac > 0.02 && jet->neuEmFrac < 0.99 && jet->nNeutrals > 2) 
                    jetPass = true;
            } else {
                if (jet->neuEmFrac < 0.9 && jet->neuHadFrac > 0.02 && jet->nNeutrals > 10) 
                    jetPass = true;
            }
        }
        else if (_parameters.period == "2018") {
            if (fabs(jet->eta) <= 2.6) {
                if (
                        jet->neuHadFrac       < 0.90
                        && jet->neuEmFrac     < 0.90
                        && jet->nParticles    > 1
                        && jet->chHadFrac     > 0
                        && jet->nCharged      > 0
                   ) {
                    jetPass = true;
                }
            }
            else if (fabs(jet->eta) <= 2.7) {
                if (
                        jet->neuHadFrac       < 0.90
                        && jet->neuEmFrac     < 0.99
                        && jet->nCharged      > 0
                   ) {
                    jetPass = true;
                }
            }
            else if (fabs(jet->eta) <= 3.0) { 
                if (jet->neuEmFrac > 0.02 && jet->neuEmFrac < 0.99 && jet->nNeutrals > 2) 
                    jetPass = true;
            }
            else {
                if (jet->neuEmFrac < 0.9 && jet->neuHadFrac > 0.2 && jet->nNeutrals > 10) 
                    jetPass = true;
            }
        }
    }
    /*std::cout << "jet eta, neuHadFrac, neuEmFrac, nParticles, chHadFrac, nCharged, chEmFrac, nNeutrals, passID: " << 
                jet->eta << ", " << jet->neuHadFrac << ", " << jet->neuEmFrac << ", " << jet->nParticles << ", " << 
                jet->chHadFrac << ", " << jet->nCharged << ", " << jet->chEmFrac << ", " << jet->nNeutrals << ", " << jetPass << std::endl;*/

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

    /*std::vector<float> vv = _jetCorrector->getSubCorrections();
    for (unsigned int i = 0; i < vv.size(); ++i) {
        std::cout << "correction " << i << " = " << vv.at(i) << std::endl;
    }*/

    double correction = _jetCorrector->getCorrection();

    return correction;
}

std::vector<float> ParticleSelector::GetJetSubcorrections(const baconhep::TJet* jet, string tagName) const 
{
    _jetCorrector->setJetEta(jet->eta);
    _jetCorrector->setJetPt(jet->ptRaw);
    _jetCorrector->setJetA(jet->area);
    _jetCorrector->setRho(_rhoFactor); 

    return _jetCorrector->getSubCorrections();
}

    
    


double ParticleSelector::JetUncertainty(const baconhep::TJet* jet) const
{
    _jecUncertainty->setJetEta(jet->eta);
    _jecUncertainty->setJetPt(jet->ptRaw);

    double uncertainty = _jecUncertainty->getUncertainty(true);

    return uncertainty;
}

//pair<float, float> ParticleSelector::JetResolutionAndSF(const baconhep::TJet* jet, int syst) const
pair<float, float> ParticleSelector::JetResolutionAndSF(const baconhep::TJet* jet, int syst)
{
    JME::JetParameters jetParams;
    jetParams.setJetPt(jet->pt);
    jetParams.setJetEta(jet->eta);
    jetParams.setRho(_rhoFactor);

    JME::JetResolution thisJetResolution;
    JME::JetResolutionScaleFactor thisJetResolutionSF;

    if (_parameters.period == "2016") {
        thisJetResolution = jetResolution["2016"];
        thisJetResolutionSF = jetResolutionSF["2016"];
    }
    else if (_parameters.period == "2017") {
        thisJetResolution = jetResolution["2017"];
        thisJetResolutionSF = jetResolutionSF["2017"];
    }
    else if (_parameters.period == "2018") {
        float rNumber = _rng->Uniform(1.);
        if (rNumber > 0.543367) {
            thisJetResolution = jetResolution["2018ABC"];
            thisJetResolutionSF = jetResolutionSF["2018ABC"];
        }
        else {
            thisJetResolution = jetResolution["2018D"];
            thisJetResolutionSF = jetResolutionSF["2018D"];
        }
    }

    //float jres = jetResolution.getResolution(jetParams);
    //float jsf = 1.;
    //if (syst == 1) {
    //    jsf  = jetResolutionSF.getScaleFactor(jetParams, Variation::UP);
    //} else if (syst == -1) {
    //    jsf  = jetResolutionSF.getScaleFactor(jetParams, Variation::DOWN);
    //} else {
    //    jsf  = jetResolutionSF.getScaleFactor(jetParams);
    //}
    
    float jres = thisJetResolution.getResolution(jetParams);
    float jsf = 1.;
    if (syst == 1) {
        jsf  = thisJetResolutionSF.getScaleFactor(jetParams, Variation::UP);
    } else if (syst == -1) {
        jsf  = thisJetResolutionSF.getScaleFactor(jetParams, Variation::DOWN);
    } else {
        jsf  = thisJetResolutionSF.getScaleFactor(jetParams);
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

float ParticleSelector::GetPFPhotonIsolation(const baconhep::TPFPart* pf_pho, const float rho) const
{
    /*int iEta = 0;
    float etaBins[8] = {0., 1., 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};
    float effArea[8] = {0.1703, 0.1715, 0.1213, 0.1230, 0.1635, 0.1937, 0.2393};
    for (unsigned i = 0; i < 8; ++i) {
        if (fabs(pf_pho->eta) > etaBins[i] && fabs(pf_pho->eta) < etaBins[i+1]) {
            iEta = i;
            break;
        }
    }*/

    //float combIso = pf_pho->chHadIso03 + std::max(0., (double)pf_pho->neuHadIso03 + pf_pho->gammaIso03 - rho*effArea[iEta]);
    float combIso = pf_pho->chHadIso03 + pf_pho->neuHadIso03 + pf_pho->gammaIso03;

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

pair<float, float> ParticleSelector::GetElectronScaleErr(const baconhep::TElectron* el) const
{
    float energyNominal = el->calibE;

    float statUp = fabs(el->energyScaleStatUp - energyNominal)/energyNominal;
    float statDown = fabs(el->energyScaleStatDown - energyNominal)/energyNominal;
    float systUp = fabs(el->energyScaleSystUp - energyNominal)/energyNominal;
    float systDown = fabs(el->energyScaleSystDown - energyNominal)/energyNominal;
    float gainUp = fabs(el->energyScaleGainUp - energyNominal)/energyNominal;
    float gainDown = fabs(el->energyScaleGainDown - energyNominal)/energyNominal;

    float scaleUp = sqrt(pow(statUp, 2) + pow(systUp, 2) + pow(gainUp, 2));
    float scaleDown = sqrt(pow(statDown, 2) + pow(systDown, 2) + pow(gainDown, 2));

    float energyUp = energyNominal + scaleUp*energyNominal;
    float energyDown = energyNominal - scaleDown*energyNominal;

    return make_pair(energyUp, energyDown);
}

float get_w(float p_data, float p_mc) {
  return 1-p_data/p_mc;
}

float get_z(float p_mc_1, float p_data_1, float p_mc_2, float p_data_2) {
  return (p_data_1-p_mc_1)/(p_mc_2-p_data_2);  
}

vector<float> ParticleSelector::GetCorrectedPhotonMVA(const baconhep::TPhoton* ph, bool isData) const
{
  const std::string _cmssw_base = getenv("CMSSW_BASE");
  // this function would return the corrected variables and also the MVA scores
  // they are ordered as:
  // SigmaIEtaIEta
  // SigmaIEtaIPhi
  // SCEtaWidth
  // SCPhiWidth
  // R9 full5x5
  // S4
  // Photon Iso
  // Charged Hardon Iso
  // Worst Charged Hadron Iso
  // MVA score
  // corrected MVA score

  vector<float> corrValues;
  corrValues.clear();

  vector<float> rnd; 
  rnd.clear();
  for (int i=0; i<6; ++ i) rnd.push_back(_rng->Rndm());
  
  static TMVA::Reader* tmvaReaderEGMPhoID[2]        = {NULL, NULL};
  static TMVA::Reader* tmvaReaderSieie[2]           = {NULL, NULL};
  static TMVA::Reader* tmvaReaderSieip[2]           = {NULL, NULL};
  static TMVA::Reader* tmvaReaderEtaWidth[2]        = {NULL, NULL};
  static TMVA::Reader* tmvaReaderPhiWidth[2]        = {NULL, NULL};
  static TMVA::Reader* tmvaReaderR9[2]              = {NULL, NULL};
  static TMVA::Reader* tmvaReaderS4[2]              = {NULL, NULL};
  static TMVA::Reader* tmvaReaderPhoIsoTail[2]      = {NULL, NULL};
  static TMVA::Reader* tmvaReaderPhoIsoData[2]      = {NULL, NULL};
  static TMVA::Reader* tmvaReaderPhoIsoMC[2]        = {NULL, NULL};
  static TMVA::Reader* tmvaReaderPhoIsoMorph[2]     = {NULL, NULL};
  static TMVA::Reader* tmvaReaderChIsoTail[2]       = {NULL, NULL};
  static TMVA::Reader* tmvaReaderChIsoData[2]       = {NULL, NULL};
  static TMVA::Reader* tmvaReaderChIsoMC[2]         = {NULL, NULL};
  static TMVA::Reader* tmvaReaderChIsoMorph[2]      = {NULL, NULL};
  static TMVA::Reader* tmvaReaderWorstChIsoTail[2]  = {NULL, NULL};
  static TMVA::Reader* tmvaReaderWorstChIsoMorph[2] = {NULL, NULL};

  static TFormula *fSieie[2]      = {NULL, NULL};
  static TFormula *fSieip[2]      = {NULL, NULL};
  static TFormula *fEtaWidth[2]   = {NULL, NULL};
  static TFormula *fPhiWidth[2]   = {NULL, NULL};
  static TFormula *fR9[2]         = {NULL, NULL};
  static TFormula *fS4[2]         = {NULL, NULL};
  static TFormula *fPhoIso[2]     = {NULL, NULL};
  static TFormula *fChIso[2]      = {NULL, NULL};
  static TFormula *fWorstChIso[2] = {NULL, NULL};
  
  //Float_t* phoEt                   = ph->pt;
  //Float_t* phoEta                  = ph->eta;
  //Float_t* phoPhi                  = ph->phi;
  //Float_t* phoR9                   = ph->r9;
  //Float_t* phoSCEta                = ph->scEta;
  //Float_t* phoSCRawE               = ph->scRawE;
  //Float_t* phoSCEtaWidth           = ph->scEtaWidth;
  //Float_t* phoSCPhiWidth           = ph->scPhiWidth;
  //Float_t* phoSigmaIEtaIEtaFull5x5 = ph->sieie;
  //Float_t* phoSigmaIEtaIPhiFull5x5 = ph->sieip;
  //Float_t* phoR9Full5x5            = ph->r9_full5x5;
  //Float_t* phoE2x2Full5x5          = ph->e2x2;
  //Float_t* phoE5x5Full5x5          = ph->e5x5;
  //Float_t* phoPFPhoIso             = ph->phoPhIso; 
  //Float_t* phoPFChIso              = ph->phoChIso; 
  //Float_t* phoPFChWorstIso         = ph->phoWorstChIso; 
  ////Float_t* phoESEnP1               = data.GetPtrFloat("phoESEnP1");
  ////Float_t* phoESEnP2               = data.GetPtrFloat("phoESEnP2");
  //Float_t* phoESEffSigmaRR         = ph->srr;
  //Float_t  rho                     = _rhoFactor;
  
  // MVA variables
  //static float phoEt_, phoEta_, phoPhi_, phoR9_, phoR9Full5x5_;
  static float phoEt_, phoPhi_, phoR9Full5x5_;
  static float phoSCEtaWidth_, phoSCPhiWidth_, rho_;
  static float phoSCEta_, phoSCRawE_;
  static float phoPFPhoIso_, phoPFChIso_, phoPFChIsoWorst_;
  static float phoESEnToRawE_, phoESEffSigmaRR_;
  static float sieieFull5x5, sieipFull5x5, s4Full5x5;
  static float rndV1, rndV2, rndV3;
  
  // 0=ECAL barrel or 1=ECAL endcaps
  int iBE = (fabs(ph->scEta) < 1.479) ? 0 : 1;
  
  // one-time MVA initialization
  if (!tmvaReaderEGMPhoID[iBE]) {

    tmvaReaderEGMPhoID[iBE] = new TMVA::Reader("!Color:Silent");

    tmvaReaderEGMPhoID[iBE]->AddVariable("SCRawE", &phoSCRawE_);
    tmvaReaderEGMPhoID[iBE]->AddVariable("r9", &phoR9Full5x5_);
    tmvaReaderEGMPhoID[iBE]->AddVariable("sigmaIetaIeta", &sieieFull5x5);
    tmvaReaderEGMPhoID[iBE]->AddVariable("etaWidth", &phoSCEtaWidth_);
    tmvaReaderEGMPhoID[iBE]->AddVariable("phiWidth", &phoSCPhiWidth_);
    tmvaReaderEGMPhoID[iBE]->AddVariable("covIEtaIPhi", &sieipFull5x5);
    tmvaReaderEGMPhoID[iBE]->AddVariable("s4", &s4Full5x5); 
    tmvaReaderEGMPhoID[iBE]->AddVariable("phoIso03", &phoPFPhoIso_);
    tmvaReaderEGMPhoID[iBE]->AddVariable("chgIsoWrtChosenVtx", &phoPFChIso_);
    tmvaReaderEGMPhoID[iBE]->AddVariable("chgIsoWrtWorstVtx", &phoPFChIsoWorst_);
    tmvaReaderEGMPhoID[iBE]->AddVariable("scEta", &phoSCEta_);
    tmvaReaderEGMPhoID[iBE]->AddVariable("rho", &rho_);
    if (iBE == 1) {
      tmvaReaderEGMPhoID[iBE]->AddVariable("esEffSigmaRR", &phoESEffSigmaRR_);
      tmvaReaderEGMPhoID[iBE]->AddVariable("esEnergyOverRawE", &phoESEnToRawE_);
    }

    if (iBE == 0) tmvaReaderEGMPhoID[0]->BookMVA("BDT", _cmssw_base + "/src/BLT/BLTAnalysis/data/tmvaWeightFiles/EgammaPhoId_94X_barrel_BDT.weights.xml"); // FIX ME
    else tmvaReaderEGMPhoID[1]->BookMVA("BDT", _cmssw_base + "/src/BLT/BLTAnalysis/data/tmvaWeightFiles/EgammaPhoId_94X_endcap_BDT.weights.xml"); // FIX ME
  }
  
  // set MVA variables
  //phoPhi_          = phoPhi;
  //phoR9_           = phoR9;
  //phoSCEta_        = phoSCEta;
  //phoSCRawE_       = phoSCRawE;
  //phoSCEtaWidth_   = phoSCEtaWidth;
  //phoSCPhiWidth_   = phoSCPhiWidth;
  //rho_             = _rhoFactor;
  ////phoESEnToRawE_   = (phoESEnP1[i]+phoESEnP2[i])/phoSCRawE[i];
  //phoESEnToRawE_   = ph->scESEn/phoSCRawE;
  //phoESEffSigmaRR_ = phoESEffSigmaRR;
  //phoEt_           = phoEt;
  //phoEta_          = phoEta;
  //phoR9Full5x5_    = phoR9Full5x5;
  //sieieFull5x5     = phoSigmaIEtaIEtaFull5x5;
  //sieipFull5x5     = phoSigmaIEtaIPhiFull5x5;
  //s4Full5x5        = phoE2x2Full5x5/phoE5x5Full5x5;
  //phoPFPhoIso_     = phoPFPhoIso;
  //phoPFChIso_      = phoPFChIso; 
  //phoPFChIsoWorst_ = phoPFChWorstIso; 
  
  phoPhi_          = ph->phi;
  //phoR9_           = ph->r9;
  phoSCEta_        = ph->scEta;
  phoSCRawE_       = ph->scRawE;
  phoSCEtaWidth_   = ph->scEtaWidth;
  phoSCPhiWidth_   = ph->scPhiWidth;
  rho_             = _rhoFactor;
  //phoESEnToRawE_   = (phoESEnP1[i]+phoESEnP2[i])/phoSCRawE[i];
  phoESEnToRawE_   = ph->scESEn/ph->scRawE;
  phoESEffSigmaRR_ = ph->srr;
  phoEt_           = ph->pt;
  //phoEta_          = ph->eta;
  phoR9Full5x5_    = ph->r9_full5x5;
  sieieFull5x5     = ph->sieie;
  sieipFull5x5     = ph->sieip;
  s4Full5x5        = ph->e2x2/ph->e5x5;
  //phoPFPhoIso_     = ph->phoPhIso;
  phoPFPhoIso_     = ph->phoNeuHadIso; // THERE WAS A BUG AT NTUPLE LEVEL SWAPPING PH and NEU ISO 
  phoPFChIso_      = ph->phoChIso; 
  phoPFChIsoWorst_ = ph->phoWorstChIso; 

  //std::cout << "photon SCEta, Pt, sieie, sieip, scEtaWidth, scPhiWidth, R9, S4, phoIso, chIso, worstChIso, MVA" << std::endl;
  //std::cout << phoSCEta_ << ", " << phoEt_ << ", " << sieieFull5x5 << ", " << sieipFull5x5 << ", " << phoSCEtaWidth_ << ", " << 
  //             phoSCPhiWidth_ << ", " << phoR9Full5x5_ << ", " << s4Full5x5 << ", " << phoPFPhoIso_ << ", " << phoPFChIso_ << ", " << 
  //             phoPFChIsoWorst_ << ", " << ph->mvaFall17V2 << std::endl;

  //std::cout << "WHAT IS LISTED AS phoNeuHadIso = " << ph->phoNeuHadIso << std::endl;
  
  float mva = tmvaReaderEGMPhoID[iBE]->EvaluateMVA("BDT");
  
  if (!isData) { 
    
    if (!tmvaReaderSieie[iBE]) {
      
      tmvaReaderSieie[iBE] = new TMVA::Reader("!Color:Silent"); 

      if (_parameters.period == "2016") {	  
	tmvaReaderSieie[iBE]->AddVariable("f0", &phoEt_);
	tmvaReaderSieie[iBE]->AddVariable("f1", &phoSCEta_);
	tmvaReaderSieie[iBE]->AddVariable("f2", &phoPhi_);
	tmvaReaderSieie[iBE]->AddVariable("f3", &rho_);      
	tmvaReaderSieie[iBE]->AddVariable("f4", &phoSCPhiWidth_);
	tmvaReaderSieie[iBE]->AddVariable("f5", &sieipFull5x5);
	tmvaReaderSieie[iBE]->AddVariable("f6", &s4Full5x5);
	tmvaReaderSieie[iBE]->AddVariable("f7", &phoR9Full5x5_);
	tmvaReaderSieie[iBE]->AddVariable("f8", &sieieFull5x5);
	tmvaReaderSieie[iBE]->AddVariable("f9", &phoSCEtaWidth_);
      } else {
	tmvaReaderSieie[iBE]->AddVariable("f0", &phoEt_);
	tmvaReaderSieie[iBE]->AddVariable("f1", &phoSCEta_);
	tmvaReaderSieie[iBE]->AddVariable("f2", &phoPhi_);
	tmvaReaderSieie[iBE]->AddVariable("f3", &rho_);
	tmvaReaderSieie[iBE]->AddVariable("f4", &sieipFull5x5);
	tmvaReaderSieie[iBE]->AddVariable("f5", &s4Full5x5);
	tmvaReaderSieie[iBE]->AddVariable("f6", &phoR9Full5x5_);
	tmvaReaderSieie[iBE]->AddVariable("f7", &phoSCPhiWidth_);
	tmvaReaderSieie[iBE]->AddVariable("f8", &sieieFull5x5);
	tmvaReaderSieie[iBE]->AddVariable("f9", &phoSCEtaWidth_);
      }
       
      if (_parameters.period == "2016") {
	if (iBE == 0) {
	  tmvaReaderSieie[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EB_sieie.xml"); // FIX ME
	  fSieie[0] = new TFormula("", "x[0]*0.0001104662098613408-3.3204630331623054e-05");
	} else {
	  tmvaReaderSieie[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EE_sieie.xml"); // FIX ME
	  fSieie[1] = new TFormula("", "x[0]*0.00047618253683114065+0.0003225000376058697");
	}
      } else if (_parameters.period == "2017") {
	if (iBE == 0) {
	  tmvaReaderSieie[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EB_sieie.xml"); // FIX ME
	  fSieie[0] = new TFormula("", "x[0]*9.584385303930019e-05-5.1717290076120845e-05");
	} else {
	  tmvaReaderSieie[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EE_sieie.xml"); // FIX ME
	  fSieie[1] = new TFormula("", "x[0]*0.0004952022425841057+0.000232510955048542");
	}
      } else if (_parameters.period == "2018") {
	if (iBE == 0) {
	  tmvaReaderSieie[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EB_sieie.xml"); // FIX ME
	  fSieie[0] = new TFormula("", "x[0]*9.447610685463367e-05-1.338798598885145e-05");
	} else {
	  tmvaReaderSieie[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EE_sieie.xml"); // FIX ME
	  fSieie[1] = new TFormula("", "x[0]*0.00047765361707808197+0.00033791164231124736");
	}
      }    
    }
    
    if (!tmvaReaderSieip[iBE]) {
      
      tmvaReaderSieip[iBE]    = new TMVA::Reader("!Color:Silent"); 

      if (_parameters.period == "2016") {
	tmvaReaderSieip[iBE]->AddVariable("f0", &phoEt_);
	tmvaReaderSieip[iBE]->AddVariable("f1", &phoSCEta_);
	tmvaReaderSieip[iBE]->AddVariable("f2", &phoPhi_);
	tmvaReaderSieip[iBE]->AddVariable("f3", &rho_);
	tmvaReaderSieip[iBE]->AddVariable("f4", &phoSCPhiWidth_);
	tmvaReaderSieip[iBE]->AddVariable("f5", &sieipFull5x5);
	tmvaReaderSieip[iBE]->AddVariable("f6", &s4Full5x5);
	tmvaReaderSieip[iBE]->AddVariable("f7", &phoR9Full5x5_);
	tmvaReaderSieip[iBE]->AddVariable("f8", &sieieFull5x5);
	tmvaReaderSieip[iBE]->AddVariable("f9", &phoSCEtaWidth_);
      } else {
	tmvaReaderSieip[iBE]->AddVariable("f0", &phoEt_);
	tmvaReaderSieip[iBE]->AddVariable("f1", &phoSCEta_);
	tmvaReaderSieip[iBE]->AddVariable("f2", &phoPhi_);
	tmvaReaderSieip[iBE]->AddVariable("f3", &rho_);
	tmvaReaderSieip[iBE]->AddVariable("f4", &sieipFull5x5);
	tmvaReaderSieip[iBE]->AddVariable("f5", &s4Full5x5);
	tmvaReaderSieip[iBE]->AddVariable("f6", &phoR9Full5x5_);
	tmvaReaderSieip[iBE]->AddVariable("f7", &phoSCPhiWidth_);
	tmvaReaderSieip[iBE]->AddVariable("f8", &sieieFull5x5);
	tmvaReaderSieip[iBE]->AddVariable("f9", &phoSCEtaWidth_);
      }
      
      if (_parameters.period == "2016") {
	if (iBE == 0) {
	  tmvaReaderSieip[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EB_sieip.xml"); // FIX ME
	  fSieip[0] = new TFormula("", "x[0]*9.256491692147542e-07-3.2280634187410214e-08");
	} else {
	  tmvaReaderSieip[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EE_sieip.xml"); // FIX ME
	  fSieip[1] = new TFormula("", "x[0]*2.2162202063046954e-05-9.087623072804566e-06");
	}
      } else if (_parameters.period == "2017") {
	if (iBE == 0) {
	  tmvaReaderSieip[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EB_sieip.xml"); // FIX ME
	  fSieip[0] = new TFormula("", "x[0]*1.0669183816085082e-06-3.2259631995695154e-08");
	} else {
	  tmvaReaderSieip[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EE_sieip.xml"); // FIX ME
	  fSieip[1] = new TFormula("", "x[0]*1.9968203087899473e-05-1.1005423958127062e-05");
	}
      } else if (_parameters.period == "2018") {
	if (iBE == 0) {
	  tmvaReaderSieip[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EB_sieip.xml"); // FIX ME
	  fSieip[0] = new TFormula("", "x[0]*1.025769961188504e-06-2.1602032511843264e-08");
	} else {
	  tmvaReaderSieip[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EE_sieip.xml"); // FIX ME
	  fSieip[1] = new TFormula("", "x[0]*1.976322925152308e-05-1.1049352291879044e-05");
	}
      }    
    }
    
    if (!tmvaReaderEtaWidth[iBE]) {
      
      tmvaReaderEtaWidth[iBE] = new TMVA::Reader("!Color:Silent"); 

      if (_parameters.period == "2016") {
	tmvaReaderEtaWidth[iBE]->AddVariable("f0", &phoEt_);
	tmvaReaderEtaWidth[iBE]->AddVariable("f1", &phoSCEta_);
	tmvaReaderEtaWidth[iBE]->AddVariable("f2", &phoPhi_);
	tmvaReaderEtaWidth[iBE]->AddVariable("f3", &rho_);
	tmvaReaderEtaWidth[iBE]->AddVariable("f4", &phoSCPhiWidth_);
	tmvaReaderEtaWidth[iBE]->AddVariable("f5", &sieipFull5x5);
	tmvaReaderEtaWidth[iBE]->AddVariable("f6", &s4Full5x5);
	tmvaReaderEtaWidth[iBE]->AddVariable("f7", &phoR9Full5x5_);
	tmvaReaderEtaWidth[iBE]->AddVariable("f8", &sieieFull5x5);
	tmvaReaderEtaWidth[iBE]->AddVariable("f9", &phoSCEtaWidth_);
      } else {
	tmvaReaderEtaWidth[iBE]->AddVariable("f0", &phoEt_);
	tmvaReaderEtaWidth[iBE]->AddVariable("f1", &phoSCEta_);
	tmvaReaderEtaWidth[iBE]->AddVariable("f2", &phoPhi_);
	tmvaReaderEtaWidth[iBE]->AddVariable("f3", &rho_);
	tmvaReaderEtaWidth[iBE]->AddVariable("f4", &sieipFull5x5);
	tmvaReaderEtaWidth[iBE]->AddVariable("f5", &s4Full5x5);
	tmvaReaderEtaWidth[iBE]->AddVariable("f6", &phoR9Full5x5_);
	tmvaReaderEtaWidth[iBE]->AddVariable("f7", &phoSCPhiWidth_);
	tmvaReaderEtaWidth[iBE]->AddVariable("f8", &sieieFull5x5);
	tmvaReaderEtaWidth[iBE]->AddVariable("f9", &phoSCEtaWidth_);
      }
      
      if (_parameters.period == "2016") {
	if (iBE == 0) {
	  tmvaReaderEtaWidth[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EB_etaWidth.xml"); // FIX ME
	  fEtaWidth[0] = new TFormula("", "x[0]*0.0003420906834020758-1.2947742904550509e-05");
	} else {
	  tmvaReaderEtaWidth[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EE_etaWidth.xml"); // FIX ME
	  fEtaWidth[1] = new TFormula("", "x[0]*0.000527115361936753+0.00032721686342924994");
	}
      } else if (_parameters.period == "2017") {
	if (iBE == 0) {
	  tmvaReaderEtaWidth[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EB_etaWidth.xml"); // FIX ME
	  fEtaWidth[0] = new TFormula("", "x[0]*0.0003493050480966255-0.00016565428928989791");
	} else {
	  tmvaReaderEtaWidth[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EE_etaWidth.xml"); // FIX ME
	  fEtaWidth[1] = new TFormula("", "x[0]*0.0005480486104594688+0.0002232480231738592");
	}
      } else if (_parameters.period == "2018") {
	if (iBE == 0) {
	  tmvaReaderEtaWidth[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EB_etaWidth.xml"); // FIX ME
	  fEtaWidth[0] = new TFormula("", "x[0]*0.0003381714524934305-3.983385867640729e-06");
	} else {
	  tmvaReaderEtaWidth[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EE_etaWidth.xml"); // FIX ME
	  fEtaWidth[1] = new TFormula("", "x[0]*0.000547223078659973+0.0004153295725724329");
	}
      }    
    }
    
    if (!tmvaReaderPhiWidth[iBE]) {
      
      tmvaReaderPhiWidth[iBE] = new TMVA::Reader("!Color:Silent"); 

      if (_parameters.period == "2016") {
	tmvaReaderPhiWidth[iBE]->AddVariable("f0", &phoEt_);
	tmvaReaderPhiWidth[iBE]->AddVariable("f1", &phoSCEta_);
	tmvaReaderPhiWidth[iBE]->AddVariable("f2", &phoPhi_);
	tmvaReaderPhiWidth[iBE]->AddVariable("f3", &rho_);
	tmvaReaderPhiWidth[iBE]->AddVariable("f4", &phoSCPhiWidth_);
	tmvaReaderPhiWidth[iBE]->AddVariable("f5", &sieipFull5x5);
	tmvaReaderPhiWidth[iBE]->AddVariable("f6", &s4Full5x5);
	tmvaReaderPhiWidth[iBE]->AddVariable("f7", &phoR9Full5x5_);
	tmvaReaderPhiWidth[iBE]->AddVariable("f8", &sieieFull5x5);
	tmvaReaderPhiWidth[iBE]->AddVariable("f9", &phoSCEtaWidth_);
      } else {
	tmvaReaderPhiWidth[iBE]->AddVariable("f0", &phoEt_);
	tmvaReaderPhiWidth[iBE]->AddVariable("f1", &phoSCEta_);
	tmvaReaderPhiWidth[iBE]->AddVariable("f2", &phoPhi_);
	tmvaReaderPhiWidth[iBE]->AddVariable("f3", &rho_);
	tmvaReaderPhiWidth[iBE]->AddVariable("f4", &sieipFull5x5);
	tmvaReaderPhiWidth[iBE]->AddVariable("f5", &s4Full5x5);
	tmvaReaderPhiWidth[iBE]->AddVariable("f6", &phoR9Full5x5_);
	tmvaReaderPhiWidth[iBE]->AddVariable("f7", &phoSCPhiWidth_);
	tmvaReaderPhiWidth[iBE]->AddVariable("f8", &sieieFull5x5);
	tmvaReaderPhiWidth[iBE]->AddVariable("f9", &phoSCEtaWidth_);
      }
      
      if (_parameters.period == "2016") {
	if (iBE == 0) {
	  tmvaReaderPhiWidth[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EB_phiWidth.xml"); // FIX ME
	  fPhiWidth[0] = new TFormula("", "x[0]*0.0025889840810865723+0.0011585117927925426");
	} else {
	  tmvaReaderPhiWidth[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EE_phiWidth.xml"); // FIX ME
	  fPhiWidth[1] = new TFormula("", "x[0]*0.00262391873289608+0.0009113112561151201");
	}
      } else if (_parameters.period == "2017") {
	if (iBE == 0) {
	  tmvaReaderPhiWidth[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EB_phiWidth.xml"); // FIX ME
	  fPhiWidth[0] = new TFormula("", "x[0]*0.002104991113194962+0.0007377486601366596");
	} else {
	  tmvaReaderPhiWidth[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EE_phiWidth.xml"); // FIX ME
	  fPhiWidth[1] = new TFormula("", "x[0]*0.0019464228762895372+0.0009612037536670028");
	}
      } else if (_parameters.period == "2018") {
	if (iBE == 0) {
	  tmvaReaderPhiWidth[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EB_phiWidth.xml"); // FIX ME
	  fPhiWidth[0] = new TFormula("", "x[0]*0.0019421618661353145+0.0007529814748217659");
	} else {
	  tmvaReaderPhiWidth[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EE_phiWidth.xml"); // FIX ME
	  fPhiWidth[1] = new TFormula("", "x[0]*0.002128917659751403+0.0012794900435001803");
	}
      }    
    }
    
    if (!tmvaReaderR9[iBE]) {
      
      tmvaReaderR9[iBE] = new TMVA::Reader("!Color:Silent"); 

      if (_parameters.period == "2016") {
	tmvaReaderR9[iBE]->AddVariable("f0", &phoEt_);
	tmvaReaderR9[iBE]->AddVariable("f1", &phoSCEta_);
	tmvaReaderR9[iBE]->AddVariable("f2", &phoPhi_);
	tmvaReaderR9[iBE]->AddVariable("f3", &rho_);
	tmvaReaderR9[iBE]->AddVariable("f4", &phoSCPhiWidth_);
	tmvaReaderR9[iBE]->AddVariable("f5", &sieipFull5x5);
	tmvaReaderR9[iBE]->AddVariable("f6", &s4Full5x5);
	tmvaReaderR9[iBE]->AddVariable("f7", &phoR9Full5x5_);
	tmvaReaderR9[iBE]->AddVariable("f8", &sieieFull5x5);
	tmvaReaderR9[iBE]->AddVariable("f9", &phoSCEtaWidth_);
      } else {
	tmvaReaderR9[iBE]->AddVariable("f0", &phoEt_);
	tmvaReaderR9[iBE]->AddVariable("f1", &phoSCEta_);
	tmvaReaderR9[iBE]->AddVariable("f2", &phoPhi_);
	tmvaReaderR9[iBE]->AddVariable("f3", &rho_);
	tmvaReaderR9[iBE]->AddVariable("f4", &sieipFull5x5);
	tmvaReaderR9[iBE]->AddVariable("f5", &s4Full5x5);
	tmvaReaderR9[iBE]->AddVariable("f6", &phoR9Full5x5_);
	tmvaReaderR9[iBE]->AddVariable("f7", &phoSCPhiWidth_);
	tmvaReaderR9[iBE]->AddVariable("f8", &sieieFull5x5);
	tmvaReaderR9[iBE]->AddVariable("f9", &phoSCEtaWidth_);
      }
	
      if (_parameters.period == "2016") {
	if (iBE == 0) {
	  tmvaReaderR9[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EB_r9.xml"); // FIX ME
	  fR9[0] = new TFormula("", "x[0]*0.00915591682273384+0.00042335154161760036");
	} else {
	  tmvaReaderR9[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EE_r9.xml"); // FIX ME
	  fR9[1] = new TFormula("", "x[0]*0.010253359666186235-0.0031073467744174854");
	}
      } else if (_parameters.period == "2017") {
	if (iBE == 0) {
	  tmvaReaderR9[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EB_r9.xml"); // FIX ME
	  fR9[0] = new TFormula("", "x[0]*0.009489495430770628-0.002028252255487417");
	} else {
	  tmvaReaderR9[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EE_r9.xml"); // FIX ME
	  fR9[1] = new TFormula("", "x[0]*0.014832657700482338-0.007372621685483638");
	}
      } else if (_parameters.period == "2018") {
	if (iBE == 0) {
	  tmvaReaderR9[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EB_r9.xml"); // FIX ME 
	  fR9[0] = new TFormula("", "x[0]*0.010041591228904162-0.0034957641959179053");
	} else {
	  tmvaReaderR9[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EE_r9.xml"); // FIX ME
	  fR9[1] = new TFormula("", "x[0]*0.014926955853444557-0.00869450644995956");
	}
      }    
    }
    
    if (!tmvaReaderS4[iBE]) {
      
      tmvaReaderS4[iBE] = new TMVA::Reader("!Color:Silent"); 

      if (_parameters.period == "2016") {
	tmvaReaderS4[iBE]->AddVariable("f0", &phoEt_);
	tmvaReaderS4[iBE]->AddVariable("f1", &phoSCEta_);
	tmvaReaderS4[iBE]->AddVariable("f2", &phoPhi_);
	tmvaReaderS4[iBE]->AddVariable("f3", &rho_);
	tmvaReaderS4[iBE]->AddVariable("f4", &phoSCPhiWidth_);
	tmvaReaderS4[iBE]->AddVariable("f5", &sieipFull5x5);
	tmvaReaderS4[iBE]->AddVariable("f6", &s4Full5x5);
	tmvaReaderS4[iBE]->AddVariable("f7", &phoR9Full5x5_);
	tmvaReaderS4[iBE]->AddVariable("f8", &sieieFull5x5);
	tmvaReaderS4[iBE]->AddVariable("f9", &phoSCEtaWidth_);
      } else {
	tmvaReaderS4[iBE]->AddVariable("f0", &phoEt_);
	tmvaReaderS4[iBE]->AddVariable("f1", &phoSCEta_);
	tmvaReaderS4[iBE]->AddVariable("f2", &phoPhi_);
	tmvaReaderS4[iBE]->AddVariable("f3", &rho_);
	tmvaReaderS4[iBE]->AddVariable("f4", &sieipFull5x5);
	tmvaReaderS4[iBE]->AddVariable("f5", &s4Full5x5);
	tmvaReaderS4[iBE]->AddVariable("f6", &phoR9Full5x5_);
	tmvaReaderS4[iBE]->AddVariable("f7", &phoSCPhiWidth_);
	tmvaReaderS4[iBE]->AddVariable("f8", &sieieFull5x5);
	tmvaReaderS4[iBE]->AddVariable("f9", &phoSCEtaWidth_);
      }
      
      if (_parameters.period == "2016") {
	if (iBE == 0) {
	  tmvaReaderS4[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EB_s4.xml"); // FIX ME
	  fS4[0] = new TFormula("", "x[0]*0.007630080763117303+0.00038370660197839523");
	} else {
	  tmvaReaderS4[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EE_s4.xml"); // FIX ME
	  fS4[1] = new TFormula("", "x[0]*0.011291341541304956-0.0049337720663549245");
	}
      } else if (_parameters.period == "2017") {
	if (iBE == 0) {
	  tmvaReaderS4[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EB_s4.xml"); // FIX ME
	  fS4[0] = new TFormula("", "x[0]*0.0065310642236871275-0.0014767145572411322");
	} else {
	  tmvaReaderS4[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EE_s4.xml"); // FIX ME
	  fS4[1] = new TFormula("", "x[0]*0.014044106856233973-0.008718535579634146");
	}
      } else if (_parameters.period == "2018") {
	if (iBE == 0) {
	  tmvaReaderS4[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EB_s4.xml"); // FIX ME
	  fS4[0] = new TFormula("", "x[0]*0.00728814761466437-0.002791432255344395");
	} else {
	  tmvaReaderS4[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EE_s4.xml"); // FIX ME
	  fS4[1] = new TFormula("", "x[0]*0.014542940630291212-0.011083602033276074");
	}
      }    
    }
    
    if (!tmvaReaderPhoIsoTail[iBE]) {
      
      tmvaReaderPhoIsoTail[iBE] = new TMVA::Reader("!Color:Silent"); 
      
      tmvaReaderPhoIsoTail[iBE]->AddVariable("f0", &phoEt_);
      tmvaReaderPhoIsoTail[iBE]->AddVariable("f1", &phoSCEta_);
      tmvaReaderPhoIsoTail[iBE]->AddVariable("f2", &phoPhi_);
      tmvaReaderPhoIsoTail[iBE]->AddVariable("f3", &rho_);
      tmvaReaderPhoIsoTail[iBE]->AddVariable("f4", &rndV1);
      
      tmvaReaderPhoIsoData[iBE] = new TMVA::Reader("!Color:Silent"); 
      
      tmvaReaderPhoIsoData[iBE]->AddVariable("f0", &phoEt_);
      tmvaReaderPhoIsoData[iBE]->AddVariable("f1", &phoSCEta_);
      tmvaReaderPhoIsoData[iBE]->AddVariable("f2", &phoPhi_);
      tmvaReaderPhoIsoData[iBE]->AddVariable("f3", &rho_);
      
      tmvaReaderPhoIsoMC[iBE] = new TMVA::Reader("!Color:Silent"); 
      
      tmvaReaderPhoIsoMC[iBE]->AddVariable("f0", &phoEt_);
      tmvaReaderPhoIsoMC[iBE]->AddVariable("f1", &phoSCEta_);
      tmvaReaderPhoIsoMC[iBE]->AddVariable("f2", &phoPhi_);
      tmvaReaderPhoIsoMC[iBE]->AddVariable("f3", &rho_);
      
      tmvaReaderPhoIsoMorph[iBE] = new TMVA::Reader("!Color:Silent"); 
      
      tmvaReaderPhoIsoMorph[iBE]->AddVariable("f0", &phoEt_);
      tmvaReaderPhoIsoMorph[iBE]->AddVariable("f1", &phoSCEta_);
      tmvaReaderPhoIsoMorph[iBE]->AddVariable("f2", &phoPhi_);
      tmvaReaderPhoIsoMorph[iBE]->AddVariable("f3", &rho_);
      tmvaReaderPhoIsoMorph[iBE]->AddVariable("f4", &phoPFPhoIso_);
      
      if (_parameters.period == "2016") {
	if (iBE == 0) {
	  tmvaReaderPhoIsoTail[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalTailRegressor_EB_phoIso.xml"); // FIX ME
	  tmvaReaderPhoIsoData[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/data_clf_p2t_EB_phoIso.xml"); // FIX ME 
	  tmvaReaderPhoIsoMC[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/mc_clf_p2t_EB_phoIso.xml"); // FIX ME
	  tmvaReaderPhoIsoMorph[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EB_phoIso.xml"); // FIX ME
	  fPhoIso[0] = new TFormula("", "x[0]*0.0935077084495785+0.039768874381154895");
	} else {
	  tmvaReaderPhoIsoTail[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalTailRegressor_EE_phoIso.xml"); // FIX ME
	  tmvaReaderPhoIsoData[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/data_clf_p2t_EE_phoIso.xml"); // FIX ME
	  tmvaReaderPhoIsoMC[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/mc_clf_p2t_EE_phoIso.xml"); // FIX ME
	  tmvaReaderPhoIsoMorph[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EE_phoIso.xml"); // FIX ME
	  fPhoIso[1] = new TFormula("", "x[0]*0.08596048057765085-0.012660368261709964");
	}
      } else if (_parameters.period == "2017") {
	if (iBE == 0) { 
	  tmvaReaderPhoIsoTail[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalTailRegressor_EB_phoIso.xml"); // FIX ME
	  tmvaReaderPhoIsoData[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/data_clf_p2t_EB_phoIso.xml"); // FIX ME
	  tmvaReaderPhoIsoMC[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/mc_clf_p2t_EB_phoIso.xml"); // FIX ME
	  tmvaReaderPhoIsoMorph[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EB_phoIso.xml"); // FIX ME
	  fPhoIso[0] = new TFormula("", "x[0]*0.1256188820912958+0.06895645737583889");
	} else {
	  tmvaReaderPhoIsoTail[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalTailRegressor_EE_phoIso.xml"); // FIX ME
	  tmvaReaderPhoIsoData[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/data_clf_p2t_EE_phoIso.xml"); // FIX ME 
	  tmvaReaderPhoIsoMC[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/mc_clf_p2t_EE_phoIso.xml"); // FIX ME 
	  tmvaReaderPhoIsoMorph[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EE_phoIso.xml"); // FIX ME
	  fPhoIso[1] = new TFormula("", "x[0]*0.0886983137006268+0.025296074198263616");
	}
      } else if (_parameters.period == "2018") {
	if (iBE == 0) {
	  tmvaReaderPhoIsoTail[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalTailRegressor_EB_phoIso.xml"); // FIX ME
	  tmvaReaderPhoIsoData[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/data_clf_p2t_EB_phoIso.xml"); // FIX ME 
	  tmvaReaderPhoIsoMC[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/mc_clf_p2t_EB_phoIso.xml"); // FIX ME 
	  tmvaReaderPhoIsoMorph[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EB_phoIso.xml"); // FIX ME
	  fPhoIso[0] = new TFormula("", "x[0]*0.12040089375899635+0.07988597084627963");
	} else {
	  tmvaReaderPhoIsoTail[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalTailRegressor_EE_phoIso.xml"); // FIX ME
	  tmvaReaderPhoIsoData[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/data_clf_p2t_EE_phoIso.xml"); // FIX ME
	  tmvaReaderPhoIsoMC[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/mc_clf_p2t_EE_phoIso.xml"); // FIX ME
	  tmvaReaderPhoIsoMorph[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EE_phoIso.xml"); // FIX ME
	  fPhoIso[1] = new TFormula("", "x[0]*0.06640579530911996+0.005250037356742843");
	}
      }    
    }
    
    // charged hadron isolation and worst charged hadron isolation
    if (!tmvaReaderChIsoTail[iBE]) {
      
      tmvaReaderChIsoTail[iBE] = new TMVA::Reader("!Color:Silent"); 
      
      tmvaReaderChIsoTail[iBE]->AddVariable("f0", &phoEt_);
      tmvaReaderChIsoTail[iBE]->AddVariable("f1", &phoSCEta_);
      tmvaReaderChIsoTail[iBE]->AddVariable("f2", &phoPhi_);
      tmvaReaderChIsoTail[iBE]->AddVariable("f3", &rho_);
      tmvaReaderChIsoTail[iBE]->AddVariable("f4", &phoPFChIsoWorst_);
      tmvaReaderChIsoTail[iBE]->AddVariable("f5", &rndV2);
      
      tmvaReaderChIsoData[iBE] = new TMVA::Reader("!Color:Silent"); 
      
      tmvaReaderChIsoData[iBE]->AddVariable("f0", &phoEt_);
      tmvaReaderChIsoData[iBE]->AddVariable("f1", &phoSCEta_);
      tmvaReaderChIsoData[iBE]->AddVariable("f2", &phoPhi_);
      tmvaReaderChIsoData[iBE]->AddVariable("f3", &rho_);
      
      tmvaReaderChIsoMC[iBE] = new TMVA::Reader("!Color:Silent"); 
      
      tmvaReaderChIsoMC[iBE]->AddVariable("f0", &phoEt_);
      tmvaReaderChIsoMC[iBE]->AddVariable("f1", &phoSCEta_);
      tmvaReaderChIsoMC[iBE]->AddVariable("f2", &phoPhi_);
      tmvaReaderChIsoMC[iBE]->AddVariable("f3", &rho_);
      
      tmvaReaderChIsoMorph[iBE] = new TMVA::Reader("!Color:Silent"); 
      
      tmvaReaderChIsoMorph[iBE]->AddVariable("f0", &phoEt_);
      tmvaReaderChIsoMorph[iBE]->AddVariable("f1", &phoSCEta_);
      tmvaReaderChIsoMorph[iBE]->AddVariable("f2", &phoPhi_);
      tmvaReaderChIsoMorph[iBE]->AddVariable("f3", &rho_);
      tmvaReaderChIsoMorph[iBE]->AddVariable("f4", &phoPFChIso_);
      tmvaReaderChIsoMorph[iBE]->AddVariable("f5", &phoPFChIsoWorst_);
      
      tmvaReaderWorstChIsoTail[iBE] = new TMVA::Reader("!Color:Silent");
      
      tmvaReaderWorstChIsoTail[iBE]->AddVariable("f0", &phoEt_);
      tmvaReaderWorstChIsoTail[iBE]->AddVariable("f1", &phoSCEta_);
      tmvaReaderWorstChIsoTail[iBE]->AddVariable("f2", &phoPhi_);
      tmvaReaderWorstChIsoTail[iBE]->AddVariable("f3", &rho_);
      tmvaReaderWorstChIsoTail[iBE]->AddVariable("f4", &phoPFChIso_);
      tmvaReaderWorstChIsoTail[iBE]->AddVariable("f5", &rndV3);
      
      tmvaReaderWorstChIsoMorph[iBE] = new TMVA::Reader("!Color:Silent"); 
      
      tmvaReaderWorstChIsoMorph[iBE]->AddVariable("f0", &phoEt_);
      tmvaReaderWorstChIsoMorph[iBE]->AddVariable("f1", &phoSCEta_);
      tmvaReaderWorstChIsoMorph[iBE]->AddVariable("f2", &phoPhi_);
      tmvaReaderWorstChIsoMorph[iBE]->AddVariable("f3", &rho_);
      tmvaReaderWorstChIsoMorph[iBE]->AddVariable("f4", &phoPFChIso_);
      tmvaReaderWorstChIsoMorph[iBE]->AddVariable("f5", &phoPFChIsoWorst_);
      
      if (_parameters.period == "2016") {
	if (iBE == 0) {
	  tmvaReaderChIsoTail[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalTailRegressor_EB_chIso.xml"); // FIX ME
	  tmvaReaderChIsoData[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/data_clf_3Cat_EB_chIso_chIsoWorst.xml");  // FIX ME
	  tmvaReaderChIsoMC[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/mc_clf_3Cat_EB_chIso_chIsoWorst.xml"); // FIX ME
	  tmvaReaderChIsoMorph[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EB_chIso.xml"); // FIX ME 
	  tmvaReaderWorstChIsoTail[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalTailRegressor_EB_chIsoWorst.xml"); // FIX ME 
	  tmvaReaderWorstChIsoMorph[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EB_chIsoWorst.xml"); // FIX ME
	  fChIso[0] = new TFormula("", "x[0]*0.06188532145612949+0.02015489488479011");
	  fWorstChIso[0] = new TFormula("", "x[0]*0.05524731145218942+0.003043560613187335");
	} else { 
	  tmvaReaderChIsoTail[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalTailRegressor_EE_chIso.xml"); // FIX ME
	  tmvaReaderChIsoData[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/data_clf_3Cat_EE_chIso_chIsoWorst.xml"); // FIX ME
	  tmvaReaderChIsoMC[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/mc_clf_3Cat_EE_chIso_chIsoWorst.xml"); // FIX ME
	  tmvaReaderChIsoMorph[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EE_chIso.xml"); // FIX ME 
	  tmvaReaderWorstChIsoTail[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalTailRegressor_EE_chIsoWorst.xml"); // FIX ME
	  tmvaReaderWorstChIsoMorph[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2016/weights_finalRegressor_EE_chIsoWorst.xml"); // FIX ME
	  fChIso[1] = new TFormula("", "x[0]*0.09449158290398274+0.005233862600715372");
	  fWorstChIso[1] = new TFormula("", "x[0]*0.07908877024409315-0.016709346023939642");
	}
      } else if (_parameters.period == "2017") {
	if (iBE == 0) {
	  tmvaReaderChIsoTail[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalTailRegressor_EB_chIso.xml"); // FIX ME
	  tmvaReaderChIsoData[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/data_clf_3Cat_EB_chIso_chIsoWorst.xml"); // FIX ME 
	  tmvaReaderChIsoMC[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/mc_clf_3Cat_EB_chIso_chIsoWorst.xml"); // FIX ME 
	  tmvaReaderChIsoMorph[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EB_chIso.xml"); // FIX ME 
	  tmvaReaderWorstChIsoTail[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalTailRegressor_EB_chIsoWorst.xml"); // FIX ME 
	  tmvaReaderWorstChIsoMorph[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EB_chIsoWorst.xml"); // FIX ME
	  fChIso[0] = new TFormula("", "x[0]*0.07795545910678636+0.04056159366893253");
	  fWorstChIso[0] = new TFormula("", "x[0]*0.09341274835827507+0.09945467977615818");
	} else {
	  tmvaReaderChIsoTail[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalTailRegressor_EE_chIso.xml"); // FIX ME
	  tmvaReaderChIsoData[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/data_clf_3Cat_EE_chIso_chIsoWorst.xml"); // FIX ME
	  tmvaReaderChIsoMC[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/mc_clf_3Cat_EE_chIso_chIsoWorst.xml"); // FIX ME
	  tmvaReaderChIsoMorph[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EE_chIso.xml"); // FIX ME
	  tmvaReaderWorstChIsoTail[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalTailRegressor_EE_chIsoWorst.xml"); // FIX ME
	  tmvaReaderWorstChIsoMorph[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2017/weights_finalRegressor_EE_chIsoWorst.xml"); // FIX ME
	  fChIso[1] = new TFormula("", "x[0]*0.08121688032488544+0.012777520521935815");
	  fWorstChIso[1] = new TFormula("", "x[0]*0.08261618370104969-0.0002665169782493232");
	}
      } else if (_parameters.period == "2018") {
	if (iBE == 0) {
	  tmvaReaderChIsoTail[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalTailRegressor_EB_chIso.xml"); // FIX ME
	  tmvaReaderChIsoData[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/data_clf_3Cat_EB_chIso_chIsoWorst.xml"); // FIX ME
	  tmvaReaderChIsoMC[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/mc_clf_3Cat_EB_chIso_chIsoWorst.xml"); // FIX ME 
	  tmvaReaderChIsoMorph[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EB_chIso.xml"); // FIX ME
	  tmvaReaderWorstChIsoTail[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalTailRegressor_EB_chIsoWorst.xml"); // FIX ME
	  tmvaReaderWorstChIsoMorph[0]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EB_chIsoWorst.xml"); // FIX ME
	  fChIso[0] = new TFormula("", "x[0]*0.09683323426563749+0.05272420382216636");
	  fWorstChIso[0] = new TFormula("", "x[0]*0.09610418976923851+0.12885914702326673");
	} else {
	  tmvaReaderChIsoTail[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalTailRegressor_EE_chIso.xml"); // FIX ME
	  tmvaReaderChIsoData[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/data_clf_3Cat_EE_chIso_chIsoWorst.xml"); // FIX ME
	  tmvaReaderChIsoMC[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/mc_clf_3Cat_EE_chIso_chIsoWorst.xml"); // FIX ME
	  tmvaReaderChIsoMorph[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EE_chIso.xml"); // FIX ME
	  tmvaReaderWorstChIsoTail[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalTailRegressor_EE_chIsoWorst.xml"); // FIX ME
	  tmvaReaderWorstChIsoMorph[1]->BookMVA("BDT", _cmssw_base +"/src/BLT/BLTAnalysis/data/tmvaWeightFiles/2018/weights_finalRegressor_EE_chIsoWorst.xml"); // FIX ME
	  fChIso[1] = new TFormula("", "x[0]*0.07311271907865564+0.015334704146744094");
	  fWorstChIso[1] = new TFormula("", "x[0]*0.07549508562988608+0.02261090428387391");
	}
      }    
    }
        
    float corrsieieFull5x5 = sieieFull5x5+fSieie[iBE]->Eval(tmvaReaderSieie[iBE]->EvaluateRegression(0, "BDT"));
    float corrsieipFull5x5 = sieipFull5x5+fSieip[iBE]->Eval(tmvaReaderSieip[iBE]->EvaluateRegression(0, "BDT"));
    float corrEtaWidth     = phoSCEtaWidth_+fEtaWidth[iBE]->Eval(tmvaReaderEtaWidth[iBE]->EvaluateRegression(0, "BDT"));
    float corrPhiWidth     = phoSCPhiWidth_+fPhiWidth[iBE]->Eval(tmvaReaderPhiWidth[iBE]->EvaluateRegression(0, "BDT"));
    float corrR9           = phoR9Full5x5_+fR9[iBE]->Eval(tmvaReaderR9[iBE]->EvaluateRegression(0, "BDT"));
    float corrS4           = s4Full5x5+fS4[iBE]->Eval(tmvaReaderS4[iBE]->EvaluateRegression(0, "BDT"));
    
    // corrected photon isolation
    auto p_tail_data = 1./(1.+sqrt(2./(1.+tmvaReaderPhoIsoData[iBE]->EvaluateMVA("BDT"))-1.));
    auto p_tail_mc   = 1./(1.+sqrt(2./(1.+tmvaReaderPhoIsoMC[iBE]->EvaluateMVA("BDT"))-1.));
    auto p_peak_data = 1 - p_tail_data;
    auto p_peak_mc   = 1 - p_tail_mc;
    auto migration_rnd_value = rnd[0];
    
    double p_move_to_tail = (p_tail_data-p_tail_mc)/p_peak_mc;
    double p_move_to_peak = (p_peak_data-p_peak_mc)/p_tail_mc;
    
    float corrPhoIso = phoPFPhoIso_;
    if (phoPFPhoIso_ == 0 && p_tail_data > p_tail_mc && migration_rnd_value < p_move_to_tail) {
      rndV1 = rnd[1]*(0.99-0.01)+0.01;
      corrPhoIso = tmvaReaderPhoIsoTail[iBE]->EvaluateRegression(0, "BDT");
    } else if (phoPFPhoIso_ > 0 && p_peak_data > p_peak_mc && migration_rnd_value <= p_move_to_peak)
      corrPhoIso = 0.;
    
    if (corrPhoIso > 0.) corrPhoIso += fPhoIso[iBE]->Eval(tmvaReaderPhoIsoMorph[iBE]->EvaluateRegression(0, "BDT"));
    
    // charged hadron isolation and worst charged hadron isolation
    auto p_00_data = tmvaReaderChIsoData[iBE]->EvaluateMulticlass(0, "BDT"); 
    auto p_01_data = tmvaReaderChIsoData[iBE]->EvaluateMulticlass(1, "BDT"); 
    auto p_11_data = tmvaReaderChIsoData[iBE]->EvaluateMulticlass(2, "BDT"); 
    auto p_00_mc   = tmvaReaderChIsoMC[iBE]->EvaluateMulticlass(0, "BDT"); 
    auto p_01_mc   = tmvaReaderChIsoMC[iBE]->EvaluateMulticlass(1, "BDT"); 
    auto p_11_mc   = tmvaReaderChIsoMC[iBE]->EvaluateMulticlass(2, "BDT");
    migration_rnd_value = rnd[2];
    rndV2 = rnd[3]*(0.99-0.01)+0.01;
    rndV3 = rnd[4]*(0.99-0.01)+0.01;
    
    float corrChIso      = phoPFChIso_;
    float corrWorstChIso = phoPFChIsoWorst_;
    
    // 00
    if (corrChIso == 0 && corrWorstChIso == 0 && p_00_mc > p_00_data && migration_rnd_value <= get_w(p_00_data, p_00_mc)) {
      
      // 00->01
      if (p_01_mc < p_01_data && p_11_mc > p_11_data) {
	corrWorstChIso = tmvaReaderWorstChIsoTail[iBE]->EvaluateRegression(0, "BDT");
      }
      // 00->11
      else if (p_01_mc > p_01_data && p_11_mc < p_11_data) {
	corrChIso      = tmvaReaderChIsoTail[iBE]->EvaluateRegression(0, "BDT");
	corrWorstChIso = tmvaReaderWorstChIsoTail[iBE]->EvaluateRegression(0, "BDT");
      }
      // 00->either 
      else if (p_01_mc < p_01_data && p_11_mc < p_11_data) {
	migration_rnd_value = rnd[5];
	if (migration_rnd_value <= get_z(p_01_mc, p_01_data, p_00_mc, p_00_data))
	  corrWorstChIso = tmvaReaderWorstChIsoTail[iBE]->EvaluateRegression(0, "BDT");
	else {
	  corrChIso      = tmvaReaderChIsoTail[iBE]->EvaluateRegression(0, "BDT");
	  corrWorstChIso = tmvaReaderWorstChIsoTail[iBE]->EvaluateRegression(0, "BDT");
	}                   
      }
    }
    // 01
    else if (corrChIso == 0. && corrWorstChIso > 0. && p_01_mc > p_01_data && migration_rnd_value <= get_w(p_01_data, p_01_mc)) {
      // 01->00
      if (p_00_mc < p_00_data && p_11_mc > p_11_data) {
	corrChIso      = 0.;
	corrWorstChIso = 0.;
      }
      // 01->11
      else if (p_00_mc>p_00_data && p_11_mc<p_11_data)
	corrChIso = tmvaReaderChIsoTail[iBE]->EvaluateRegression(0, "BDT");
      // 01->either
      else if (p_00_mc < p_00_data && p_11_mc < p_11_data) {
	migration_rnd_value = rnd[5];
	if (migration_rnd_value <= get_z(p_00_mc, p_00_data, p_01_mc, p_01_data))
	  corrWorstChIso = 0.;
	else
	  corrChIso = tmvaReaderChIsoTail[iBE]->EvaluateRegression(0, "BDT");
      }
    }
    // 11
    else if (corrChIso > 0. && corrWorstChIso > 0. && p_11_mc > p_11_data && migration_rnd_value <= get_w(p_11_data, p_11_mc)) {                
      // 11->00
      if (p_00_mc < p_00_data && p_01_mc > p_01_data) {
	corrChIso      = 0.;
	corrWorstChIso = 0.;
      }
      // 11->01
      else if (p_00_mc > p_00_data && p_01_mc < p_01_data)
	corrChIso = 0.;
      // 11->either
      else if (p_00_mc < p_00_data && p_01_mc < p_01_data) {
	migration_rnd_value = rnd[5];
	if (migration_rnd_value <= get_z(p_00_mc, p_00_data, p_11_mc, p_11_data)) {
	  corrChIso      = 0.;
	  corrWorstChIso = 0.;
	}
	else
	  corrChIso = 0.;
      }
    }
    
    if (corrChIso > 0.) corrChIso += fChIso[iBE]->Eval(tmvaReaderChIsoMorph[iBE]->EvaluateRegression(0, "BDT"));
    if (corrWorstChIso > 0.) corrWorstChIso += fWorstChIso[iBE]->Eval(tmvaReaderWorstChIsoMorph[iBE]->EvaluateRegression(0, "BDT"));
    
    if (corrChIso > corrWorstChIso) swap(corrChIso, corrWorstChIso);
    /*
    cout<<"Sieie      : "<<sieieFull5x5<<" "<<corrsieieFull5x5<<endl;
    cout<<"Sieip      : "<<sieipFull5x5<<" "<<corrsieipFull5x5<<endl;
    cout<<"EtaWidth   : "<<phoSCEtaWidth_<<" "<<corrEtaWidth<<endl;
    cout<<"PhiWidth   : "<<phoSCPhiWidth_<<" "<<corrPhiWidth<<endl;
    cout<<"R9         : "<<phoR9Full5x5_<<" "<<corrR9<<endl;
    cout<<"S4         : "<<s4Full5x5<<" "<<corrS4<<endl;
    cout<<"PhoIso     : "<<phoPFPhoIso_<<" "<<corrPhoIso<<endl;
    cout<<"ChIso      : "<<phoPFChIso_<<" "<<corrChIso<<endl;
    cout<<"WorstChIso : "<<phoPFChIsoWorst_<<" "<<corrWorstChIso<<endl;
    */
    corrValues.push_back(corrsieieFull5x5);
    corrValues.push_back(corrsieipFull5x5);
    corrValues.push_back(corrEtaWidth);
    corrValues.push_back(corrPhiWidth);
    corrValues.push_back(corrR9);
    corrValues.push_back(corrS4);
    corrValues.push_back(corrPhoIso); 
    corrValues.push_back(corrChIso);
    corrValues.push_back(corrWorstChIso);
    
    sieieFull5x5     = corrsieieFull5x5;
    sieipFull5x5     = corrsieipFull5x5;
    phoSCEtaWidth_   = corrEtaWidth;
    phoSCPhiWidth_   = corrPhiWidth;
    phoR9Full5x5_    = corrR9;
    //  s4Full5x5        = corrS4;
    phoPFPhoIso_     = corrPhoIso;
    phoPFChIso_      = corrChIso;
    phoPFChIsoWorst_ = corrWorstChIso;

    
    float corrmva = tmvaReaderEGMPhoID[iBE]->EvaluateMVA("BDT");   
    
    //std::cout << "photon corr_sieie, corr_sieip, corr_scEtaWidth, corr_scPhiWidth, corr_R9, corr_S4, corr_phoIso, corr_ChIso, corr_worstChIso, mva, corr_MVA" << std::endl;
    //std::cout << sieieFull5x5 << ", " << sieipFull5x5 << ", " << phoSCEtaWidth_ << ", " << phoSCPhiWidth_ << ", " << phoR9Full5x5_ << ", " << s4Full5x5 << ", " << 
    //             phoPFPhoIso_ << ", " << phoPFChIso_ << ", " << phoPFChIsoWorst_ << ", " << mva << ", " << corrmva << std::endl;

    corrValues.push_back(mva);
    corrValues.push_back(corrmva);
    //cout<<"mva        : "<<mva<<" "<<corrmva<<endl;
    //
    //std::cout << "mva POG, preCorr, postCorrNoS4:" << std::endl;
    //std::cout << ph->mvaFall17V2 << ", " << mva << ", " << corrmva << std::endl;
    //s4Full5x5        = corrS4;
    //corrmva = tmvaReaderEGMPhoID[iBE]->EvaluateMVA("BDT");   
    //std::cout << corrmva << std::endl;

    
  } else {
    corrValues.push_back(0);
    corrValues.push_back(0);
    corrValues.push_back(0);
    corrValues.push_back(0);
    corrValues.push_back(0);
    corrValues.push_back(0);
    corrValues.push_back(0);
    corrValues.push_back(0);
    corrValues.push_back(0);
    corrValues.push_back(mva);
    corrValues.push_back(0);
  }
  
  return corrValues;
}
