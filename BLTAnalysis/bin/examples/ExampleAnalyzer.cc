#include "ExampleAnalyzer.h"
#include <map>


//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;


// bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 


ExampleAnalyzer::ExampleAnalyzer() : BLTSelector()
{

}

ExampleAnalyzer::~ExampleAnalyzer()
{

}

void ExampleAnalyzer::Begin(TTree *tree)

{   cout<< "begin"<<endl;
    rng = new TRandom3();


    // Parse command line option
    cout<< "Parse command line option"<<endl;
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1), std::sregex_token_iterator(), std::back_inserter(options));

    // Set the parameters
    cout<< "Set the parameters"<<endl;
    params.reset(new Parameters());
    params->setup(options);

    // Set the cuts
    cout<< "Set the cuts"<<endl;
    cuts.reset(new Cuts());
    particleSelector.reset(new ParticleSelector(*params, *cuts));

    // Trigger bits mapping file
    cout<< "Trigger bits mapping file"<<endl;
    const std::string cmssw_base = getenv("CMSSW_BASE");
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns";
    trigger.reset(new baconhep::TTrigger(trigfilename));

    if (params->selection == "single_lepton") {
        triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");

    } else if (params->selection == "single_muon") {
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");

    } else if (params->selection == "single_electron") {
        triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");
    }

    // Weight utility class
    cout<< "Weight utility class"<<endl;
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
    //string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt";
    lumiMask.AddJSONFile(jsonFileName);

    // muon momentum corrections
    cout<< "muon momentum corrections"<<endl;
    muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/rcdata.2016.v3");

    // electron scale corrections
    cout<< "electron scale corrections"<<endl;
    electronScaler = new EnergyScaleCorrection(cmssw_base + "/src/BLT/BLTAnalysis/data");

    // Prepare the output tree
    cout<< "Prepare the output tree"<<endl;
    string outFileName = params->get_output_filename("output");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();

    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);

    vector<std::string> channelNames = {"mumu", "ee", "emu", 
                                        "etau", "mutau",
                                        };

    for (unsigned i = 0; i < channelNames.size(); ++i) {
        string channel = channelNames[i];
        outFile->mkdir(channel.c_str());
        outFile->cd(channel.c_str());
        string treeName = "bltTree_" + params->datasetgroup;
        tree = new TTree(treeName.c_str(), treeName.c_str());

        // event data
        tree->Branch("runNumber", &runNumber);
        tree->Branch("evtNumber", &evtNumber, "eventNumber/l");
        tree->Branch("lumiSection", &lumiSection);
        tree->Branch("nPV", &nPV);
        tree->Branch("nPU", &nPU);
        tree->Branch("triggerLeptonStatus",&triggerLeptonStatus);

        // gen level objects
        tree->Branch("nPartons", &nPartons);
        tree->Branch("genWeight", &genWeight);
        
        // weights and their uncertainties
        tree->Branch("eventWeight", &eventWeight);

        // leptons
        tree->Branch("leptonOneP4" , &leptonOneP4);
        tree->Branch("leptonOneIso", &leptonOneIso);
        tree->Branch("leptonOneFlavor", &leptonOneFlavor);
        tree->Branch("leptonTwoP4" , &leptonTwoP4);
        tree->Branch("leptonTwoIso", &leptonTwoIso);
        tree->Branch("leptonTwoFlavor", &leptonTwoFlavor);

        // taus
        if (channel == "mutau" || channel == "etau" ) {
            tree->Branch("tauDecayMode",     &tauDecayMode);
            tree->Branch("tauMVA",           &tauMVA); 
            tree->Branch("tauGenFlavor",     &tauGenFlavor);
            tree->Branch("tauGenFlavorHad",  &tauGenFlavorHad);
            tree->Branch("tauVetoedJetPt",   &tauVetoedJetPt);
            tree->Branch("tauVetoedJetPtUnc",&tauVetoedJetPtUnc);
        }

        // met and ht
        tree->Branch("met", &met);
        tree->Branch("metPhi", &metPhi);

        // object counters
        tree->Branch("nMuons", &nMuons);
        tree->Branch("nElectrons", &nElectrons);
        tree->Branch("nTaus", &nTaus);
        tree->Branch("nJets", &nJets);
        tree->Branch("nBJets", &nBJets);

        outTrees[channel] = tree;

        // event counter
        string outHistName = params->get_output_treename("TotalEvents_" + channel);
        eventCounts[channel] = new TH1D(outHistName.c_str(),"ChannelCounts",10,0.5,10.5);
    }
    cout<< "Done with initialization"<<endl;
    ReportPostBegin();

}

Bool_t ExampleAnalyzer::Process(Long64_t entry)
{
    /* beginning of event loop */
    GetEntry(entry, 1);  // load all branches
    outFile->cd();
    eventWeight = 1.;
    this->totalEvents++;
    hTotalEvents->Fill(1);
    
    if (entry%10000==0) {
        std::cout << "... Processing event " << entry 
            << " Run: " << fInfo->runNum 
            << " Lumi: " << fInfo->lumiSec 
            << " Event: " << fInfo->evtNum 
            << std::endl;
    }


    const bool isData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isData);

    /* Apply lumi mask */
    if (isData) {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(2);

    /* Trigger selection */
    bool passTrigger = false;
    vector<string> passTriggerNames;
    for (unsigned i = 0; i < triggerNames.size(); ++i) {
        bool triggered = false;
        triggered = trigger->pass(triggerNames[i], fInfo->triggerBits);
        passTrigger |= triggered;

        if (triggered) {
            passTriggerNames.push_back(triggerNames[i]);
        }
    }

    bool muonTriggered = find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_IsoMu24_v*") != passTriggerNames.end();
    muonTriggered = muonTriggered || find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_IsoTkMu24_v*") != passTriggerNames.end();
    bool electronTriggered = find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_Ele27_WPTight_Gsf_v*") != passTriggerNames.end();


    if (!passTrigger)
        return kTRUE;
    hTotalEvents->Fill(3);




    /////////////////////
    // Fill event info //
    /////////////////////

    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
    nPV           = fPVArr->GetEntries();
    triggerStatus = passTrigger;


    ///////////////////////
    // Generator objects //
    ///////////////////////

    vector<TGenParticle*> genTaus, genTausLep, genTausHad;

    if (!isData) {

        // Set data period for 2016 MC scale factors
        if (rng->Rndm() < 0.468) {
            weights->SetDataPeriod("2016BtoF");    
        } else {
            weights->SetDataPeriod("2016GH");
        }

        // save gen weight for amc@nlo Drell-Yan sample
        genWeight = fGenEvtInfo->weight > 0 ? 1 : -1;
        if (genWeight < 0) {
            hTotalEvents->Fill(10);
        }


        // pileup reweighting
        nPU          = fInfo->nPUmean;
        puWeight     = weights->GetPUWeight(fInfo->nPUmean); 
        eventWeight *= puWeight;


        // loop over gen particle collection:
        //   * parton counting for combining inclusive and jet-binned Drell-Yan samples
        //   * save gen level tau and leptonic taus

        unsigned count = 0;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            // parton counting for jet-binned Drell-Yan samples
            if (
                    particle->status == 23 
                    && (abs(particle->pdgId) < 6 || particle->pdgId == 21) 
                    && particle->parent != -2
               ) {
                ++count;
            }

            // This saves gen level tau and leptonic taus
            if ( particle->parent != -2 ) {
                if ( abs(particle->pdgId) == 16) {
                    TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                    if (abs(mother->pdgId) == 15) { 
                        genTaus.push_back(mother); 
                    }
                }

                if ( abs(particle->pdgId) == 11 || abs(particle->pdgId) == 13 ) {
                    TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                    if (abs(mother->pdgId) == 15) {
                        genTausLep.push_back(mother);
                    }
                }
            }
        }
        // finish loop over gen particles. now conclude gen level infomation
        //  *  get gen level nPartons
        //  *  get hadronic taus 

        // get nPartons
        nPartons = count; 
        
        // genTausHad = genTaus - genTausLep
        for (unsigned i = 0; i < genTaus.size(); ++i) {
            
            bool isHad = true;

            TLorentzVector genTauP4;
            genTauP4.SetPtEtaPhiM(genTaus[i]->pt, genTaus[i]->eta, genTaus[i]->phi, genTaus[i]->mass); 

            for (unsigned j = 0; j < genTausLep.size(); ++j) {
                TLorentzVector genTauLepP4;
                genTauLepP4.SetPtEtaPhiM(genTausLep[j]->pt, genTausLep[j]->eta, genTausLep[j]->phi, genTausLep[j]->mass); 

                if (genTauP4.DeltaR(genTauLepP4) < 0.1) {
                    isHad = false;
                }
            }
            if (isHad) genTausHad.push_back(genTaus[i]);
        }

    } else {
        nPU = 0;
        nPartons = 0;
    }



    ///////////////////
    // Select objects//
    ///////////////////


    /* -------- Vertices ---------*/

    //cout<<fInfo->hasGoodPV<<endl;

    if (fInfo->hasGoodPV) {
        assert(fPVArr->GetEntries() != 0);
        TVector3 pv;
        copy_xyz((TVertex*) fPVArr->At(0), pv);
        particleSelector->SetPV(pv);
    } else {
        return kTRUE;
    }
    hTotalEvents->Fill(4);
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);


    
    /* -------- MUONS ---------*/
    vector<TMuon*> muons, fail_muons;
    vector<TLorentzVector> veto_muons;
    for (int i=0; i < fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);

        // Apply rochester muon momentum corrections
        double muonSF = 1.;
        if (isData) {
            muonSF = muonCorr->kScaleDT(muon->q, muon->pt, muon->eta, muon->phi, 0, 0);
        } else {
            muonSF = muonCorr->kScaleAndSmearMC(muon->q, muon->pt, muon->eta, muon->phi,
                    muon->nTkLayers, rng->Rndm(), rng->Rndm(), 
                    0, 0);
        }
        muon->pt = muonSF*muon->pt; 
        
        TLorentzVector muonP4;
        muonP4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, 0.1052);


        if (
                muonP4.Pt() > 10.
                && fabs(muonP4.Eta()) < 2.4

                // tight muon ID and ISO
                && (muon->typeBits & baconhep::kPFMuon) 
                && (muon->typeBits & baconhep::kGlobal) 
                && muon->muNchi2    < 10.
                && muon->nMatchStn  > 1
                && muon->nPixHits   > 0
                && muon->nTkLayers  > 5 
                && muon->nValidHits > 0
                && fabs(muon->d0)   < 0.2
                && fabs(muon->dz)   < 0.5
                && GetMuonIsolation(muon)/muonP4.Pt() <  0.15
           ) {
                muons.push_back(muon);
                veto_muons.push_back(muonP4);
            } 
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);
    nMuons = muons.size();



    
    /* -------- ELECTRONS ---------*/
    vector<TElectron*> electrons;
    vector<TLorentzVector> veto_electrons;

    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);
        
        // apply electron scale and smear
        if (isData) {
            scaleData sdata = electronScaler->GetScaleData(electron, runNumber);
            electron->pt *= sdata.scale;
        } else {
            float sFactor = electronScaler->GetSmearingFactor(electron, 0, 0);
            electron->pt *= rng->Gaus(1, sFactor);
        }

        TLorentzVector electronP4;
        electronP4.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, 511e-6);

        if (
                electron->pt > 10
                && fabs(electron->scEta) < 2.5
                && particleSelector->PassElectronID(electron, cuts->tightElID)
                // && (GetElectronIsolation(electron, fInfo->rhoJet)/electron->pt) < 0.1
                && particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
           ) {
                electrons.push_back(electron);
                veto_electrons.push_back(electronP4);
            }
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);
    nElectrons = electrons.size();

    
    /* -------- TAUS ---------*/
    vector<TTau*> taus;
    vector<TLorentzVector> veto_taus;
    for (int i=0; i < fTauArr->GetEntries(); i++) {
        TTau *tau = (TTau*) fTauArr->At(i);
        assert(tau);

        TLorentzVector tauP4; 
        tauP4.SetPtEtaPhiM(tau->pt, tau->eta, tau->phi, tau->m);

        // Prevent overlap of muons and jets
        bool muOverlap = false;
        for (const auto& mu: veto_muons) {
            if (tauP4.DeltaR(mu) < 0.3) {
                muOverlap = true;
                break;
            }
        }
        bool elOverlap = false;
        for (const auto& el: veto_electrons) {
            if (tauP4.DeltaR(el) < 0.3) {
                elOverlap = true;
                break;
            }
        }


        // apply tau energy scale correction (https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Tau_energy_scale)
        if (!isData) {

            // check if can be tagged as hadronic tau
            bool genTauHadOverlap = false;
            for (unsigned i = 0; i < genTausHad.size(); ++i) {
                TLorentzVector genP4;
                genP4.SetPtEtaPhiM(genTausHad[i]->pt, genTausHad[i]->eta, genTausHad[i]->phi, genTausHad[i]->mass); 
                if (tauP4.DeltaR(genP4) < 0.3) {
                    genTauHadOverlap = true;
                    break;
                }
            }

            // apply the tau energy scale correction only for true tau ID
            if (genTauHadOverlap) {
                if (tau->decaymode == 0) {
                    tau->pt *= 0.995;
                } else if (tau->decaymode == 1) {
                    tau->pt *= 1.01;
                } else if (tau->decaymode == 10) {
                    tau->pt *= 1.006;
                }
            }

        }

        if( 
                tau->pt > 20  
                && abs(tau->eta) < 2.3 
                && !muOverlap
                && !elOverlap
                && (tau->hpsDisc & baconhep::kByDecayModeFinding)
                && (tau->hpsDisc & baconhep::kByTightIsolationMVA3newDMwLT)
                && (tau->hpsDisc & baconhep::kByMVA6VTightElectronRejection)
                && (tau->hpsDisc & baconhep::kByTightMuonRejection3)
            ) {
            taus.push_back(tau);
            veto_taus.push_back(tauP4);
            }
    }
    sort(taus.begin(), taus.end(), sort_by_higher_pt<TTau>);
    nTaus = taus.size();
    

    
    /* -------- JETS ---------*/
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets,bjets,vetoedJets;

    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

        // apply JEC offline and get scale uncertainties
        double jec = particleSelector->JetCorrector(jet, "NONE");
        jet->pt = jet->ptRaw*jec;


        float gRand = 1.;
        float jerc = 1.;
        if (!isData) { // apply jet energy resolution corrections to simulation
            pair<float, float> resPair = particleSelector->JetResolutionAndSF(jet, 0);
            gRand = rng->Gaus(0, resPair.first);
            jerc = 1 + gRand*sqrt(std::max((double)resPair.second*resPair.second - 1, 0.));
            jet->pt = jet->pt*jerc;
        }

        // Prevent overlap of muons and jets
        TLorentzVector jetP4; 
        jetP4.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
        bool muOverlap = false;
        for (const auto& mu: veto_muons) {
            if (jetP4.DeltaR(mu) < 0.4) {
                muOverlap = true;
                break;
            }
        }

        bool elOverlap = false;
        for (const auto& el: veto_electrons) {
            if (jetP4.DeltaR(el) < 0.4) {
                elOverlap = true;
                break;
            }
        }

        bool tauOverlap = false;
        for (const auto& tau: veto_taus) {
            if (jetP4.DeltaR(tau) < 0.4) {
                tauOverlap = true;
                vetoedJets.push_back(jet);
                break;
            }
        }

        if (
                fabs(jet->eta) <= 2.4
                && jet->pt > 30
                && particleSelector->PassJetID(jet, cuts->looseJetID)
                && !muOverlap 
                && !elOverlap
                && !tauOverlap
            ) {

            jets.push_back(jet);
            if (jet->bmva > 0.9432) {
                bjets.push_back(jet);
            }
            
        }
    }
    std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);
    nJets  = jets.size();
    nBJets = bjets.size();


    /* -------- MET ---------*/
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;


    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////


    string channel = "";
    if        (electrons.size() == 0 && muons.size() == 2 && taus.size() == 0 ){ // mu+mu selection
        

        channel = "mumu";
        eventCounts[channel]->Fill(1);

        if (!muonTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        if (muons[0]->pt < 25. || muons[1]->pt < 10.)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        if (nJets < 2 || nBJets < 1)
            return kTRUE;
        eventCounts[channel]->Fill(4);


        // convert to TLorentzVectors
        TLorentzVector muonOneP4, muonTwoP4, dimuonP4;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, 0.1052);
        dimuonP4 = muonOneP4 + muonTwoP4;

        float muonOneIso = GetMuonIsolation(muons[0]);
        float muonTwoIso = GetMuonIsolation(muons[1]);


        leptonOneP4     = muonOneP4;
        leptonOneIso    = muonOneIso;
        leptonOneFlavor = muons[0]->q*13;
        leptonTwoP4     = muonTwoP4;
        leptonTwoIso    = muonTwoIso;
        leptonTwoFlavor = muons[1]->q*13;

        // test if selected lepton fire the trigger.
        bitset<2> triggered;
        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, muons[0]->hltMatchBits) && muonOneP4.Pt() > 25)
                triggered.set(0);
            if (trigger->passObj(name, 1, muons[1]->hltMatchBits) && muonTwoP4.Pt() > 25)
                triggered.set(1);
        }
        if (!triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 0;
        if ( triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 1;
        if (!triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 2;
        if ( triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 3;

        // correct for MC, including reconstruction and trigger
        if (!isData) {

            // correct for reconstruction weights
            EfficiencyContainer effCont1, effCont2;
            effCont1 = weights->GetMuonRecoEff(muonOneP4);
            effCont2 = weights->GetMuonRecoEff(muonTwoP4);


            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));


            effs = effCont2.GetEff();
            errs = effCont2.GetErr();
            leptonTwoRecoWeight = effs.first/effs.second;
            leptonTwoRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            
            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;

            // correct for trigger.

            // check if lepton could pass the trigger threshold and is matched
            // to a trigger object.  When both muons pass the trigger, use the
            // efficiency for detecting either
            

            if (triggered.all()) {
                effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
                effCont2      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
                triggerWeight = GetTriggerSF(effCont1, effCont2);
                triggerVar    = GetTriggerSFError(effCont1, effCont2);
            } else if (triggered.test(0)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else if (triggered.test(1)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }
            eventWeight *= triggerWeight;
        }
    } else if (electrons.size() == 2 && muons.size() == 0 && taus.size() == 0) { // e+e selection


        channel = "ee";
        eventCounts[channel]->Fill(1);


        if (electrons[0]->pt < 30 || electrons[1]->pt < 10 )
            return kTRUE;
        eventCounts[channel]->Fill(2);

        if (!electronTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        if (nJets < 2 || nBJets < 1)
            return kTRUE;
        eventCounts[channel]->Fill(4);


        TLorentzVector electronOneP4, electronTwoP4, dielectronP4;
        electronOneP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        electronTwoP4.SetPtEtaPhiM(electrons[1]->pt, electrons[1]->eta, electrons[1]->phi, 511e-6);
        dielectronP4 = electronOneP4 + electronTwoP4;


        leptonOneP4     = electronOneP4;
        leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonOneFlavor = 11*electrons[0]->q;
        leptonTwoP4     = electronTwoP4;
        leptonTwoIso    = GetElectronIsolation(electrons[1], fInfo->rhoJet);
        leptonTwoFlavor = 11*electrons[1]->q;


        // test if selected lepton fire the trigger.
        bitset<2> triggered;
        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, electrons[0]->hltMatchBits) && electronOneP4.Pt() > 30)
                triggered.set(0);
            if (trigger->passObj(name, 1, electrons[1]->hltMatchBits) && electronTwoP4.Pt() > 30)
                triggered.set(1);
        }
        if (!triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 0;
        if ( triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 1;
        if (!triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 2;
        if ( triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 3;


        // correct for MC, including reconstruction and trigger
        if (!isData) {

            // correct for reconstruction weights

            EfficiencyContainer effCont1, effCont2;
            effCont1 = weights->GetElectronRecoEff(electronOneP4);
            effCont2 = weights->GetElectronRecoEff(electronTwoP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            effs = effCont2.GetEff();
            errs = effCont2.GetErr();
            leptonTwoRecoWeight = effs.first/effs.second;
            leptonTwoRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;

            // correct for trigger.

            // check if lepton could pass the trigger threshold and is matched
            // to a trigger object.  When both muons pass the trigger, use the
            // efficiency for detecting either

            if (triggered.all()) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronOneP4);
                effCont2      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronTwoP4);
                triggerWeight = GetTriggerSF(effCont1, effCont2);
                triggerVar    = GetTriggerSFError(effCont1, effCont2);
            } else if (triggered.test(0)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronOneP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else if (triggered.test(1)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronTwoP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }
            eventWeight *= triggerWeight;

        }
    } else if (electrons.size() == 1 && muons.size() == 1 && taus.size() == 0) { // e+mu selection
        channel = "emu";
        eventCounts[channel]->Fill(1);

        if (muons[0]->pt < 10. || electrons[0]->pt < 20.)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        if (!(muonTriggered || electronTriggered))
            return kTRUE;
        eventCounts[channel]->Fill(3);

        if (nJets < 2 || nBJets < 1)
            return kTRUE;
        eventCounts[channel]->Fill(4);


        TLorentzVector muonP4, electronP4, dileptonP4;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        electronP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        dileptonP4 = muonP4 + electronP4;


        leptonOneP4     = muonP4;
        leptonOneIso    = GetMuonIsolation(muons[0]);
        leptonOneFlavor = 13*muons[0]->q;
        leptonTwoP4     = electronP4;
        leptonTwoIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonTwoFlavor = 11*electrons[0]->q;


        // test if selected lepton fire the trigger.
        bitset<2> triggered;
        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, muons[0]->hltMatchBits))
                triggered.set(0);
            if (trigger->passObj(name, 1, electrons[0]->hltMatchBits))
                triggered.set(1);
        }
        if (!triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 0;
        if ( triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 1;
        if (!triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 2;
        if ( triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 3;


        // correct for MC, including reconstruction and trigger
        if (!isData) {

            // reconstruction weights
            EfficiencyContainer effCont1, effCont2;
            effCont1 = weights->GetMuonRecoEff(muonP4);
            effCont2 = weights->GetElectronRecoEff(electronP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            effs = effCont2.GetEff();
            errs = effCont2.GetErr();
            leptonTwoRecoWeight = effs.first/effs.second;
            leptonTwoRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;


            // correct for trigger.

            // check if lepton could pass the trigger threshold and is matched
            // to a trigger object.  When both muons pass the trigger, use the
            // efficiency for detecting either

            if (triggered.all()) {
                effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
                effCont2      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronP4);
                triggerWeight = GetTriggerSF(effCont1, effCont2);
                triggerVar    = GetTriggerSFError(effCont1, effCont2);
            } else if (triggered.test(0)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else if (triggered.test(1)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }
            eventWeight *= triggerWeight;
        }
    } else if (electrons.size() == 1 && muons.size() == 0 && taus.size() == 1) { // e+tau selection

        channel = "etau";
        eventCounts[channel]->Fill(1);

        if (electrons[0]->pt < 30 )
            return kTRUE;
        eventCounts[channel]->Fill(2);

        if (!electronTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        if (nJets < 2)
            return kTRUE;
        eventCounts[channel]->Fill(4);

        TLorentzVector electronP4, tauP4;
        electronP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        tauP4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, taus[0]->m);


        leptonOneP4     = electronP4;
        leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonOneFlavor = 11*electrons[0]->q;
        leptonTwoP4     = tauP4;
        leptonTwoIso    = 0.;
        leptonTwoFlavor = 15*taus[0]->q;

        ///////////tau info///////////////////
        tauDecayMode      = taus[0]->decaymode;
        tauMVA            = taus[0]->rawIsoMVA3newDMwLT;

        pair <float, float> tauVetoedJetPtPair;
        tauVetoedJetPtPair= GetTauVetoedJetPt(tauP4, vetoedJets);
        tauVetoedJetPt    = tauVetoedJetPtPair.first;
        tauVetoedJetPtUnc = tauVetoedJetPtPair.second;
        //////////////////////////////////////


        // test if selected lepton fire the trigger.
        bool triggered = false;
        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, electrons[0]->hltMatchBits)) {
                triggered = true;
                break;
            }
        }
        if (!triggered) triggerLeptonStatus = 0;
        if ( triggered) triggerLeptonStatus = 1;

        // correct for MC, including reconstruction and trigger
        if (!isData) {

            tauGenFlavor    = GetTauGenFlavor(tauP4,genTausHad,vetoedJets, false); // useHadronFlavor = false
            tauGenFlavorHad = GetTauGenFlavor(tauP4,genTausHad,vetoedJets, true);  // useHadronFlavor = true

            // reconstruction weights
            EfficiencyContainer effCont1;
            effCont1 = weights->GetElectronRecoEff(electronP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            leptonTwoRecoWeight = 0.95;
            leptonTwoRecoVar    = 0.05;
            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;

            // correct for trigger.
            if (triggered) {
                effCont1 = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronP4);
                effs = effCont1.GetEff();
                errs = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(triggerWeight, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }
            eventWeight *= triggerWeight;
            
            

        }
    } else if (electrons.size() == 0 && muons.size() == 1 && taus.size() == 1) { // mu+tau selection

        channel = "mutau";
        eventCounts[channel]->Fill(1);

        if (muons[0]->pt < 25)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        if (!muonTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        if (nJets < 2)
            return kTRUE;
        eventCounts[channel]->Fill(4);

        TLorentzVector muonP4, tauP4;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        tauP4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, taus[0]->m);

        leptonOneP4     = muonP4;
        leptonOneIso    = GetMuonIsolation(muons[0]);
        leptonOneFlavor = 13*muons[0]->q;
        leptonTwoP4     = tauP4;
        leptonTwoIso    = 0.;
        leptonTwoFlavor = 15*taus[0]->q;


        ///////////tau info///////////////////
        tauDecayMode      = taus[0]->decaymode;
        tauMVA            = taus[0]->rawIsoMVA3newDMwLT;

        pair <float, float> tauVetoedJetPtPair;
        tauVetoedJetPtPair= GetTauVetoedJetPt(tauP4, vetoedJets);
        tauVetoedJetPt    = tauVetoedJetPtPair.first;
        tauVetoedJetPtUnc = tauVetoedJetPtPair.second;
        //////////////////////////////////////


        // test if selected lepton fire the trigger.
        bool triggered = false;
        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, muons[0]->hltMatchBits)) {
                triggered = true;
                break;
            }
        }
        if (!triggered) triggerLeptonStatus = 0;
        if ( triggered) triggerLeptonStatus = 1;


        // correct for MC, including reconstruction and trigger
        if (!isData) {
            tauGenFlavor    = GetTauGenFlavor(tauP4,genTausHad,vetoedJets, false); // useHadronFlavor = false
            tauGenFlavorHad = GetTauGenFlavor(tauP4,genTausHad,vetoedJets, true);  // useHadronFlavor = true

            // reconstruction weights
            EfficiencyContainer effCont1, effCont2;
            effCont1 = weights->GetMuonRecoEff(muonP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            leptonTwoRecoWeight = 0.95;
            leptonTwoRecoVar    = 0.05;
            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;

            // correct for trigger.
            if (triggered) {
                effCont1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
                effs = effCont1.GetEff();
                errs = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(triggerWeight, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }
            eventWeight *= triggerWeight;
            

        }
    } else {
        return kTRUE;
    }


    ///////////////////
    // Fill jet info //
    ///////////////////

    outFile->cd(channel.c_str());
    outTrees[channel]->Fill();
    this->passedEvents++;
    return kTRUE;
}

void ExampleAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void ExampleAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void ExampleAnalyzer::ReportPostTerminate()
{
    std::cout << "  ==== Terminate Job =========================================" << std::endl;
    std::cout << "  output   : " << params->get_output_filename("demoFile") << std::endl;
    std::cout << "           : Processed " << this->fileCount << " files with " << this->unskimmedEventCount << " unskimmed events." << std::endl;
    std::cout << "           : Selected " << this->passedEvents << " / " << this->totalEvents << " events." << std::endl;
    std::cout << "  ============================================================" << std::endl;
}


// _____________________________________________________________________________
// Main function

int main(int argc, char **argv)
{
    std::unique_ptr<ExampleAnalyzer> selector(new ExampleAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


float ExampleAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    return combIso;
}

float ExampleAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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




float ExampleAnalyzer::GetTriggerSF(EfficiencyContainer eff1, EfficiencyContainer eff2)
{
    pair<double, double> trigEff1, trigEff2;
    trigEff1 = eff1.GetEff();
    trigEff2 = eff2.GetEff();

    float sf = 1.;
    if (trigEff1.second > 0 && trigEff2.second <= 0){
        sf *= trigEff1.first/trigEff1.second;
    }
    
    if (trigEff1.second <= 0 && trigEff2.second > 0){
        sf *= trigEff2.first/trigEff2.second;
    }
    if (trigEff1.second > 0 && trigEff2.second > 0){
        sf *= trigEff1.first/trigEff1.second;
        sf *= trigEff2.first/trigEff2.second;
    }

    //if (trigEff1.second > 0 || trigEff2.second > 0) {
    //  sf = (1 - (1 - trigEff1.first)*(1 - trigEff2.first))/(1 - (1 - trigEff1.second)*(1 - trigEff2.second));
    //}


    return sf;
}

float ExampleAnalyzer::GetTriggerSFError(EfficiencyContainer eff1, EfficiencyContainer eff2)
{
    pair<double, double> trigEff1, trigEff2, trigErr1, trigErr2;
    trigEff1 = eff1.GetEff();
    trigEff2 = eff2.GetEff();
    trigErr1 = eff1.GetErr();
    trigErr2 = eff2.GetErr();

    //cout << trigEff1.first << " :: " << trigEff1.second << " :: "
    //     << trigEff2.first << " :: " << trigEff2.second << endl;

    //cout << trigErr1.first << " :: " << trigErr1.second << " :: "
    //     << trigErr2.first << " :: " << trigErr2.second << endl;

    float denom = 1 - (1 - trigEff1.second)*(1 - trigEff2.second);
    float dEffData1 = (1 - trigEff2.first)/denom;
    float dEffData2 = (1 - trigEff1.first)/denom;
    float dEffMC1   = -(1 - trigEff2.second)/pow(denom, 2);
    float dEffMC2   = -(1 - trigEff1.second)/pow(denom, 2);

    float sfVar = pow(trigErr1.first*dEffData1, 2)
        + pow(trigErr2.first*dEffData2, 2)
        + pow(trigErr1.second*dEffMC1, 2)
        + pow(trigErr2.second*dEffMC2, 2);
    return sfVar;
}


int ExampleAnalyzer::GetTauGenFlavor(TLorentzVector p4, vector<TGenParticle*> genTausHad, vector<TJet*> vetoedJets, bool useHadronFlavor )
{
    int flavor = 26;
    
    
    // check if can be tagged as hadronic tau
    for (unsigned i = 0; i < genTausHad.size(); ++i) {
        TLorentzVector genP4;
        genP4.SetPtEtaPhiM(genTausHad[i]->pt, genTausHad[i]->eta, genTausHad[i]->phi, genTausHad[i]->mass); 
        if (genP4.DeltaR(p4) < 0.3) {
            flavor = 15;
        }
    }

    // check if can be tagged by jet flavor
    if (flavor!= 15){
        float jetPtMax = - 1.0;
        for (unsigned i = 0; i < vetoedJets.size(); ++i) {


            TLorentzVector jetP4; 
            jetP4.SetPtEtaPhiM(vetoedJets[i]->pt, vetoedJets[i]->eta, vetoedJets[i]->phi, vetoedJets[i]->mass);
            if (jetP4.DeltaR(p4) < 0.4 && jetP4.Pt()>jetPtMax) {
                // cout<<"flavor Jet: "<< jetPtMax<<","<< jetP4.Pt()<< endl;
                if(useHadronFlavor) {
                    flavor = vetoedJets[i]->hadronFlavor;
                } else {
                    flavor = abs(vetoedJets[i]->partonFlavor);
                }
                jetPtMax = jetP4.Pt();
            }
        }
    }
    return flavor;
}



pair<float, float> ExampleAnalyzer::GetTauVetoedJetPt(TLorentzVector p4, vector<TJet*> vetoedJets)
{

    float jetPt = -1;
    float jetPtUnc = -1;
    

    // check if can be tagged by jet flavor

    float jetPtMax = - 1.0;
    for (unsigned i = 0; i < vetoedJets.size(); ++i) {

        TLorentzVector jetP4; 
        jetP4.SetPtEtaPhiM(vetoedJets[i]->pt, vetoedJets[i]->eta, vetoedJets[i]->phi, vetoedJets[i]->mass);

        if (jetP4.DeltaR(p4) < 0.4 && jetP4.Pt()>jetPtMax) {
            // cout<<"Energy Jet: "<<jetPtMax<<","<< jetP4.Pt() << endl;

            // save 
            jetPt       = jetP4.Pt();
            jetPtUnc    = jetPt * float(particleSelector->JetUncertainty(vetoedJets[i]));
            jetPtUnc    = abs(jetPtUnc);

            jetPtMax = jetP4.Pt();
        }    
    }
    
    return make_pair(jetPt,jetPtUnc);
}
