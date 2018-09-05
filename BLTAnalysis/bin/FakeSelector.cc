#include "FakeSelector.h"
#include <map>

//
// See header file for class documentation
//


using namespace std;
using namespace baconhep;

const float ZMASS = 91.19;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 
bool sort_by_btag(const baconhep::TJet* lhs, const baconhep::TJet* rhs) 
{
    return lhs->bmva > rhs->bmva;
}


FakeSelector::FakeSelector() : BLTSelector()
{

}

FakeSelector::~FakeSelector()
{

}

void FakeSelector::Begin(TTree *tree)
{
    rng = new TRandom3();

    // Parse command line option
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1), std::sregex_token_iterator(), std::back_inserter(options));

    // Set the parameters
    params.reset(new Parameters());
    params->setup(options);

    // Set the cuts
    cuts.reset(new Cuts());
    particleSelector.reset(new ParticleSelector(*params, *cuts));

    // Trigger bits mapping file
    const std::string cmssw_base = getenv("CMSSW_BASE");
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns";
    trigger.reset(new baconhep::TTrigger(trigfilename));
    if (params->selection == "single_lepton") {
        triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");
    } else {
        triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");
    } 

    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask

    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
    lumiMask.AddJSONFile(jsonFileName);

    // muon momentum corrections
    muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/rcdata.2016.v3");

    // Prepare the output tree
    string outFileName = params->get_output_filename("output");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();

    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);

    vector<std::string> channelNames = {"ee_mu","mumu_mu","ee_e","mumu_e","emu_tau"};

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
        tree->Branch("triggerStatus", &triggerStatus);

        tree->Branch("nPV", &nPV);
        tree->Branch("nPU", &nPU);
        tree->Branch("nPartons", &nPartons);
        tree->Branch("rPV", &rPV);
        tree->Branch("eventWeight", &eventWeight);

        // met and ht
        tree->Branch("met", &met);
        tree->Branch("metPhi", &metPhi);
        tree->Branch("ht", &ht);
        tree->Branch("htPhi", &htPhi);
        tree->Branch("htSum", &htSum);

        // leptons
        tree->Branch("leptonOneIso", &leptonOneIso);
        tree->Branch("leptonOneP4", &leptonOneP4);
        tree->Branch("leptonOneFlavor", &leptonOneFlavor);
        tree->Branch("leptonOneMother", &leptonOneMother);
        tree->Branch("leptonOneD0", &leptonOneD0);
        tree->Branch("leptonOneDZ", &leptonOneDZ);

        tree->Branch("leptonTwoP4", &leptonTwoP4);
        tree->Branch("leptonTwoIso", &leptonTwoIso);
        tree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
        tree->Branch("leptonTwoMother", &leptonTwoMother);
        tree->Branch("leptonTwoD0", &leptonTwoD0);
        tree->Branch("leptonTwoDZ", &leptonTwoDZ);

        tree->Branch("leptonThreeP4", &leptonThreeP4);
        tree->Branch("leptonThreeIso", &leptonThreeIso);
        tree->Branch("leptonThreeIsoPass", &leptonThreeIsoPass);
        tree->Branch("leptonThreeFlavor", &leptonThreeFlavor);
        tree->Branch("leptonThreeMother", &leptonThreeMother);
        tree->Branch("leptonThreeD0", &leptonThreeD0);
        tree->Branch("leptonThreeDZ", &leptonThreeDZ);


        if (channel == "emu_tau") {
            //tree->Branch("tauChHadMult",  &tauChHadMult);
            //tree->Branch("tauPhotonMult", &tauPhotonMult);
            tree->Branch("tauDecayMode",  &tauDecayMode);
            tree->Branch("tauMVA",        &tauMVA);
            tree->Branch("tauMVAOld",     &tauMVAOld);
            
            tree->Branch("taupuppiChHadIso",     &taupuppiChHadIso);
            tree->Branch("taupuppiGammaIso",     &taupuppiGammaIso);
            tree->Branch("taupuppiNeuHadIso",    &taupuppiNeuHadIso);
            tree->Branch("taupuppiChHadIsoNoLep",     &taupuppiChHadIsoNoLep);
            tree->Branch("taupuppiGammaIsoNoLep",     &taupuppiGammaIsoNoLep);
            tree->Branch("taupuppiNeuHadIsoNoLep",    &taupuppiNeuHadIsoNoLep);

        }

        // gen level objects
        tree->Branch("genCategory", &genCategory);


        // object counters
        tree->Branch("nMuons", &nMuons);
        tree->Branch("nElectrons", &nElectrons);
        tree->Branch("nPhotons", &nPhotons);
        tree->Branch("nJets", &nJets);
        tree->Branch("nFwdJets", &nFwdJets);
        tree->Branch("nBJets", &nBJets);
        
        outTrees[channel] = tree;

        // event counter
        string outHistName = params->get_output_treename("TotalEvents_" + channel);
        eventCounts[channel] = new TH1D(outHistName.c_str(),"ChannelCounts",10,0.5,10.5);
    }

    ReportPostBegin();
}

Bool_t FakeSelector::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
    outFile->cd();
    eventWeight = 1.;
    this->totalEvents++;
    hTotalEvents->Fill(1);

    if (entry%10000==0)  
        std::cout << "... Processing event " << entry 
            << " Run: " << fInfo->runNum 
            << " Lumi: " << fInfo->lumiSec 
            << " Event: " << fInfo->evtNum 
            << std::endl;

    const bool isData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isData);

    ///////////////////////
    // Generator objects //
    ///////////////////////

    vector<TGenParticle*> genParticles;
    vector<int> genMotherId;

    if (!isData) {

        // Set data period for 2016 MC scale factors
        if (rng->Rndm() < 0.468) {
            weights->SetDataPeriod("2016BtoF");    
        } else {
            weights->SetDataPeriod("2016GH");
        }

        float topSF = 1.;
        unsigned count = 0;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            if (
                    particle->status == 23 
                    && (abs(particle->pdgId) < 6 || particle->pdgId == 21) 
                    && particle->parent != -2
               ) {
                ++count;
            }

            // top pt reweighting: get the scale factor based on the top quark pt
            if (abs(particle->pdgId) == 6 && particle->status == 62) {
                topSF *= exp(0.0615 - 0.0005*particle->pt);
            }
        }
        nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY

        // Account for the top pt weights

        if (params->datasetgroup.substr(0, 5) == "ttbar") {
            float topPtWeight = sqrt(topSF);
            eventWeight *= topPtWeight;
        }
    } else {
        nPartons = 0;
    }

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

    if (!passTrigger)
        return kTRUE;
    hTotalEvents->Fill(3);


    /////////////////////
    // Fill event info //
    /////////////////////

    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
    triggerStatus = passTrigger;
    nPV           = fPVArr->GetEntries();
    if (!isData) {
        eventWeight *= weights->GetPUWeight(fInfo->nPUmean); // pileup reweighting
    } else {
        nPU = 0;

    }

    ///////////////////
    // Select objects//
    ///////////////////

    /* Vertices */
    if (fInfo->hasGoodPV) {
        assert(fPVArr->GetEntries() != 0);
        TVector3 pv;
        copy_xyz((TVertex*) fPVArr->At(0), pv);
        rPV = pv; 
        particleSelector->SetPV(pv);
    } else {
        return kTRUE;
    }
    hTotalEvents->Fill(4);
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);

    /* MUONS */
    vector<TMuon*> muons;
    vector<TLorentzVector> veto_muons;
    for (int i=0; i < fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);

        // Apply rochester muon momentum corrections
        TLorentzVector muonP4;
        copy_p4(muon, MUON_MASS, muonP4);
        double muonSF = 1.;
        if (isData) {
            muonSF = muonCorr->kScaleDT(muon->q, muon->pt, muon->eta, muon->phi, 0, 0);
        } else {
            muonSF = muonCorr->kScaleAndSmearMC(muon->q, muon->pt, muon->eta, muon->phi,
                    muon->nTkLayers, rng->Rndm(), rng->Rndm(), 
                    0, 0);
        }
        muon->pt = muonSF*muon->pt; 
        muonP4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);

        // Remove muons with very small deltaR
        float minDeltaR = 1e6;
        for (unsigned j=0; j < muons.size(); ++j) {
            TLorentzVector tmpMuonP4;
            tmpMuonP4.SetPtEtaPhiM(muons[j]->pt, muons[j]->eta, muons[j]->phi, 0.1051);
            float dr = muonP4.DeltaR(tmpMuonP4);
            if (dr < minDeltaR) {
                minDeltaR = dr;
            }
        }


        if (
                muonP4.Pt() > 10.
                && fabs(muonP4.Eta()) < 2.4

                // tight muon ID and ISO
                && (muon->typeBits & baconhep::kPFMuon) 
                && (muon->typeBits & baconhep::kGlobal) 
                && muon->muNchi2    < 10.
                && muon->nMatchStn  > 1
                && muon->nPixHits   > 0
                && fabs(muon->d0)   < 0.2
                && fabs(muon->dz)   < 0.5
                && muon->nTkLayers  > 5 
                && muon->nValidHits > 0
           ) {
            muons.push_back(muon);

            if (GetMuonIsolation(muon)/muonP4.Pt() < 0.15) {
                veto_muons.push_back(muonP4);
            } 
        }
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);

    /* ELECTRONS */
    vector<TElectron*> electrons;
    vector<TLorentzVector> veto_electrons;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        TLorentzVector electronP4;
        electronP4.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, 511e-6);
        //cout << electron->regscale << ", " << electron->regsmear << endl;

        if (
                electron->pt > 10
                && fabs(electron->eta) < 2.5
                && particleSelector->PassElectronID(electron, cuts->tightElID)
                
           ) {
            electrons.push_back(electron);

            if (particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)){
                veto_electrons.push_back(electronP4);
            } 
        }

    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);


    /* TAUS */
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

        if( 
                tau->pt > 20
                && abs(tau->eta) < 2.3 
                && !muOverlap
                && !elOverlap
                && (tau->hpsDisc & baconhep::kByDecayModeFinding)
                // && (tau->hpsDisc & baconhep::kByVTightIsolationMVA3newDMwLT)
                && (tau->hpsDisc & baconhep::kByMVA6VTightElectronRejection)
                && (tau->hpsDisc & baconhep::kByTightMuonRejection3)
          ) {
            taus.push_back(tau);
            if (tau->hpsDisc & baconhep::kByVTightIsolationMVA3newDMwLT){
                veto_taus.push_back(tauP4);
            }
            
        }
    }
    sort(taus.begin(), taus.end(), sort_by_higher_pt<TTau>);
    


    /* JETS */
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets;
    TLorentzVector hadronicP4;
    float sumJetPt = 0;

    nJets    = 0;
    nFwdJets = 0;
    nBJets   = 0;
    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

        if (isData) { // fix for broken bacon JEC
            double jec = particleSelector->JetCorrector(jet, "NONE");
            jet->pt = jet->ptRaw*jec;
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

        if (
                jet->pt > 30 
                && fabs(jet->eta) < 4.7
                && particleSelector->PassJetID(jet, cuts->looseJetID)
                && !muOverlap 
                && !elOverlap
           ) {

            if (fabs(jet->eta) <= 2.4) { 
                hadronicP4 += jetP4;
                sumJetPt += jetP4.Pt();
                jets.push_back(jet);

                if (isData) {
                    if (jet->bmva > 0.9432) { 
                        ++nBJets;
                    } else {
                        ++nJets;
                    }
                } else {
                    if (particleSelector->BTagModifier(jet, "MVAT", 0, 0, rng->Uniform(1.))) { 
                        ++nBJets;
                    } else {
                        ++nJets;
                    }
                }
            } else {
                if (fabs(jet->eta) > 2.5) {
                    hadronicP4 += jetP4;
                    sumJetPt += jetP4.Pt();
                    ++nFwdJets;
                }
            }
        }
    }
    std::sort(jets.begin(), jets.end(), sort_by_btag);

    /* MET */
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;

    TVector2 metP2;
    metP2.SetMagPhi(met, metPhi);

    /* HT */
    htSum = sumJetPt;
    ht    = hadronicP4.Pt();
    htPhi = hadronicP4.Phi();

    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    nMuons     = muons.size();
    nElectrons = electrons.size();
    nTaus      = taus.size();

    // trigger selections
    bool muonTriggered =    find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_IsoMu24_v*") != passTriggerNames.end() 
                        ||  find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_IsoTkMu24_v*") != passTriggerNames.end();
    bool electronTriggered = find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_Ele27_WPTight_Gsf_v*") != passTriggerNames.end();

    string channel = "";




    if (muons.size() >= 3 && electrons.size() == 0) { // mumu+mu selection
        channel = "mumu_mu";
        eventCounts[channel]->Fill(1);


        ////////////////////////////////////////////////////////////////////////////////////////
        // https://github.com/NWUHEP/BLT/blob/topic_wbranch/BLTAnalysis/bin/FakeSelector.cc#L452
        // pick muon pair that has a mass closest to the Z pole
        ////////////////////////////////////////////////////////////////////////////////////////
        
        float ix[] = {0, 1, 2};
        float massDiff = 1e9;
        vector<TMuon*> tmp_muons = muons;

        for (unsigned i = 0; i < muons.size() - 1; ++i) {
            for (unsigned j = i+1; j < muons.size(); ++j) {

                TLorentzVector tmpOne, tmpTwo, tmpDimuon;
                tmpOne.SetPtEtaPhiM(muons[i]->pt, muons[i]->eta, muons[i]->phi, 0.1052);
                tmpTwo.SetPtEtaPhiM(muons[j]->pt, muons[j]->eta, muons[j]->phi, 0.1052);
                tmpDimuon = tmpOne + tmpTwo;
                float dm = fabs(tmpDimuon.M() - 91);

                if (
                    muons[i]->pt > 25. 
                    && muons[j]->pt > 10.
                    && GetMuonIsolation(muons[i])/tmpOne.Pt() < 0.15 
                    && GetMuonIsolation(muons[j])/tmpTwo.Pt() < 0.15
                    && bool(muons[i]->q != muons[j]->q)
                    && dm < massDiff
                ) {

                    massDiff = dm;

                    ix[0] = i;
                    ix[1] = j;
                    tmp_muons[0] = muons[i];
                    tmp_muons[1] = muons[j];
                    
                }
            }
        }
        // find the probe muon (highest pt muon not in the Z pair)
        for (unsigned i = 0; i < muons.size(); ++i) {
            if (i != ix[0] && i != ix[1]) {
                ix[2] = i;
                tmp_muons[2] = muons[i];
                break;
            }
        }
        muons = tmp_muons;
        //////////////////////////////////////////////////////////////////////////
        // finish
        //////////////////////////////////////////////////////////////////////////




    
        if (   !muonTriggered
            || muons[0]->pt <25 || GetMuonIsolation(muons[0])/muons[0]->pt > 0.15 
            || muons[1]->pt <20 || GetMuonIsolation(muons[1])/muons[1]->pt > 0.15 
            ) {
            return kTRUE;
        }
        eventCounts[channel]->Fill(2);


        // apply preselection and set tree vars
        TLorentzVector muonOneP4, muonTwoP4, dimuonP4;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, 0.1052);

        dimuonP4 = muonOneP4 + muonTwoP4;
        float muonOneIso = GetMuonIsolation(muons[0]);
        float muonTwoIso = GetMuonIsolation(muons[1]);

        if ( fabs(dimuonP4.M() - 91) > 15 || muons[0]->q == muons[1]->q)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        leptonOneP4     = muonOneP4;
        leptonOneIso    = muonOneIso;
        leptonOneFlavor = muons[0]->q*13;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;

        leptonTwoP4     = muonTwoP4;
        leptonTwoIso    = muonTwoIso;
        leptonTwoFlavor = muons[1]->q*13;
        leptonTwoD0     = muons[1]->d0;
        leptonTwoDZ     = muons[1]->dz;

        TLorentzVector muonThreeP4;
        muonThreeP4.SetPtEtaPhiM(muons[2]->pt, muons[2]->eta, muons[2]->phi, 0.1052);
        leptonThreeP4     = muonThreeP4;
        leptonThreeIso    = GetMuonIsolation(muons[2]);
        leptonThreeFlavor = muons[2]->q*13;
        leptonThreeDZ     = muons[2]->dz;
        leptonThreeD0     = muons[2]->d0;

        if (!isData) {

            // reconstruction weights
            EfficiencyContainer effCont1, effCont2;
            effCont1 = weights->GetMuonRecoEff(muonOneP4);
            effCont2 = weights->GetMuonRecoEff(muonTwoP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            float leptonOneRecoWeight = effs.first/effs.second;

            effs = effCont2.GetEff();
            errs = effCont2.GetErr();
            float leptonTwoRecoWeight = effs.first/effs.second;

            // trigger weights with trigger matching
            bitset<2> triggered;
            for (const auto& name: passTriggerNames) {
                if (trigger->passObj(name, 1, muons[0]->hltMatchBits) && muonOneP4.Pt() > 25)
                    triggered.set(0);
                if (trigger->passObj(name, 1, muons[1]->hltMatchBits) && muonTwoP4.Pt() > 25)
                    triggered.set(1);
            }

            float triggerWeight = 1.;
            if (triggered.any()) {
                if (triggered.test(0)) {
                    effCont1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
                }
                if (triggered.test(1)) {
                    effCont2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
                }
                triggerWeight = GetTriggerSF(effCont1, effCont2);

            } else {
                return kTRUE;
            }

            // update the event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;
        }
    } else if (electrons.size() == 2 && muons.size() == 1) { // ee+mu selection
        channel = "ee_mu";
        eventCounts[channel]->Fill(1);

        // pick muon pair that has a mass closest to the Z pole
        if(    !electronTriggered
            || electrons[0]->pt < 30. || !particleSelector->PassElectronIso(electrons[0], cuts->tightElIso, cuts->EAEl)
            || electrons[1]->pt < 20. || !particleSelector->PassElectronIso(electrons[1], cuts->tightElIso, cuts->EAEl)
        
        )
            return kTRUE;
        eventCounts[channel]->Fill(2);

        // apply preselection and set tree vars
        TLorentzVector electronOneP4, electronTwoP4, dielectronP4;
        electronOneP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 5.11e-6);
        electronTwoP4.SetPtEtaPhiM(electrons[1]->pt, electrons[1]->eta, electrons[1]->phi, 5.11e-6);
        dielectronP4 = electronOneP4 + electronTwoP4;
        if ( fabs(dielectronP4.M() - 91) > 15 || electrons[0]->q == electrons[1]->q)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        leptonOneP4     = electronOneP4;
        leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonOneFlavor = electrons[0]->q*11;
        leptonOneDZ     = electrons[0]->dz;
        leptonOneD0     = electrons[0]->d0;

        leptonTwoP4     = electronTwoP4;
        leptonTwoIso    = GetElectronIsolation(electrons[1], fInfo->rhoJet);
        leptonTwoFlavor = electrons[1]->q*11;
        leptonTwoDZ     = electrons[1]->dz;
        leptonTwoD0     = electrons[1]->d0;

        TLorentzVector muonP4;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        leptonThreeP4     = muonP4;
        leptonThreeIso    = GetMuonIsolation(muons[0]);
        leptonThreeFlavor = muons[0]->q*13;
        leptonThreeDZ     = muons[0]->dz;
        leptonThreeD0     = muons[0]->d0;

        if (!isData) {

            EfficiencyContainer effCont1, effCont2;
            effCont1 = weights->GetElectronRecoEff(leptonOneP4);
            effCont2 = weights->GetElectronRecoEff(leptonTwoP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            float leptonOneRecoWeight = effs.first/effs.second;

            effs = effCont2.GetEff();
            errs = effCont2.GetErr();
            float leptonTwoRecoWeight = effs.first/effs.second;

            // trigger weights with trigger matching
            bitset<2> triggered;
            for (const auto& name: passTriggerNames) {
                if (trigger->passObj(name, 1, electrons[0]->hltMatchBits) && leptonOneP4.Pt() > 30)
                    triggered.set(0);
                if (trigger->passObj(name, 1, electrons[1]->hltMatchBits) && leptonTwoP4.Pt() > 30)
                    triggered.set(1);
            }

            float triggerWeight = 1.;
            if (triggered.all()) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", leptonOneP4);
                effCont2      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", leptonTwoP4);
                triggerWeight = GetTriggerSF(effCont1, effCont2);
            } else if (triggered.test(0)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", leptonOneP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
            } else if (triggered.test(1)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", leptonTwoP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
            } else {
                return kTRUE;
            }

            // update event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;
        }
    } else if (muons.size() == 2 && electrons.size() == 1) { // mumu+e selection
        channel = "mumu_e";
        eventCounts[channel]->Fill(1);

        if (   !muonTriggered
            || muons[0]->pt <25 || GetMuonIsolation(muons[0])/muons[0]->pt > 0.15 
            || muons[1]->pt <20 || GetMuonIsolation(muons[1])/muons[1]->pt > 0.15 )
        {
            return kTRUE;
        }
        eventCounts[channel]->Fill(2);


        // apply preselection and set tree vars
        TLorentzVector muonOneP4, muonTwoP4, dimuonP4;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, 0.1052);

        dimuonP4 = muonOneP4 + muonTwoP4;
        float muonOneIso = GetMuonIsolation(muons[0]);
        float muonTwoIso = GetMuonIsolation(muons[1]);

        if ( fabs(dimuonP4.M() - 91) > 15 || muons[0]->q == muons[1]->q)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        leptonOneP4     = muonOneP4;
        leptonOneIso    = muonOneIso;
        leptonOneFlavor = muons[0]->q*13;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;

        leptonTwoP4     = muonTwoP4;
        leptonTwoIso    = muonTwoIso;
        leptonTwoFlavor = muons[1]->q*13;
        leptonTwoD0     = muons[1]->d0;
        leptonTwoDZ     = muons[1]->dz;

        TLorentzVector electronOneP4;
        electronOneP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 5.11e-6);
        leptonThreeP4     = electronOneP4;
        leptonThreeIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonThreeFlavor = electrons[0]->q*11;
        leptonThreeDZ     = electrons[0]->dz;
        leptonThreeD0     = electrons[0]->d0;
        
        if (particleSelector->PassElectronIso(electrons[0], cuts->tightElIso, cuts->EAEl) )
            leptonThreeIsoPass = 1;
        else
            leptonThreeIsoPass = 0;

        if (!isData) {

            // reconstruction weights
            EfficiencyContainer effCont1, effCont2;
            effCont1 = weights->GetMuonRecoEff(muonOneP4);
            effCont2 = weights->GetMuonRecoEff(muonTwoP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            float leptonOneRecoWeight = effs.first/effs.second;

            effs = effCont2.GetEff();
            errs = effCont2.GetErr();
            float leptonTwoRecoWeight = effs.first/effs.second;

            // trigger weights with trigger matching
            bitset<2> triggered;
            for (const auto& name: passTriggerNames) {
                if (trigger->passObj(name, 1, muons[0]->hltMatchBits) && muonOneP4.Pt() > 25)
                    triggered.set(0);
                if (trigger->passObj(name, 1, muons[1]->hltMatchBits) && muonTwoP4.Pt() > 25)
                    triggered.set(1);
            }

            float triggerWeight = 1.;
            if (triggered.any()) {
                if (triggered.test(0)) {
                    effCont1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
                }
                if (triggered.test(1)) {
                    effCont2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
                }
                triggerWeight = GetTriggerSF(effCont1, effCont2);

            } else {
                return kTRUE;
            }

            // update the event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;
        }
    } else if (electrons.size() >= 3 && muons.size() == 0) { // ee+e selection
        channel = "ee_e";
        eventCounts[channel]->Fill(1);

        ////////////////////////////////////////////////////////////////////////////////////////
        // https://github.com/NWUHEP/BLT/blob/topic_wbranch/BLTAnalysis/bin/FakeSelector.cc#L452
        // pick muon pair that has a mass closest to the Z pole
        ////////////////////////////////////////////////////////////////////////////////////////
        
        float ix[] = {0, 1, 2};
        float massDiff = 1e9;
        vector<TElectron*> tmp_electrons = electrons;

        for (unsigned i = 0; i < electrons.size() - 1; ++i) {
            for (unsigned j = i+1; j < electrons.size(); ++j) {

                TLorentzVector tmpOne, tmpTwo, tmpDielectron;
                tmpOne.SetPtEtaPhiM(electrons[i]->pt, electrons[i]->eta, electrons[i]->phi, 5.11e-6);
                tmpTwo.SetPtEtaPhiM(electrons[j]->pt, electrons[j]->eta, electrons[j]->phi, 5.11e-6);
                tmpDielectron = tmpOne + tmpTwo;
                float dm = fabs(tmpDielectron.M() - 91);

                if (
                    electrons[i]->pt > 30. 
                    && electrons[j]->pt > 15.
                    && particleSelector->PassElectronIso(electrons[i], cuts->tightElIso, cuts->EAEl)
                    && particleSelector->PassElectronIso(electrons[j], cuts->tightElIso, cuts->EAEl)
                    && bool(electrons[i]->q != electrons[j]->q)
                    && dm < massDiff
                ) {

                    massDiff = dm;

                    ix[0] = i;
                    ix[1] = j;
                    tmp_electrons[0] = electrons[i];
                    tmp_electrons[1] = electrons[j];
                    
                }
            }
        }
        // find the probe electron (highest pt electron not in the Z pair)
        for (unsigned i = 0; i < electrons.size(); ++i) {
            if (i != ix[0] && i != ix[1]) {
                ix[2] = i;
                tmp_electrons[2] = electrons[i];
                break;
            }
        }
        electrons = tmp_electrons;

        //////////////////////////////////////////////////////////////////////////
        // finish
        //////////////////////////////////////////////////////////////////////////



        // pick muon pair that has a mass closest to the Z pole
        if(   !electronTriggered
            || electrons[0]->pt < 30. || !particleSelector->PassElectronIso(electrons[0], cuts->tightElIso, cuts->EAEl)
            || electrons[1]->pt < 20. || !particleSelector->PassElectronIso(electrons[1], cuts->tightElIso, cuts->EAEl)
        )
            return kTRUE;
        eventCounts[channel]->Fill(2);

        // apply preselection and set tree vars
        TLorentzVector electronOneP4, electronTwoP4, dielectronP4;
        electronOneP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 5.11e-6);
        electronTwoP4.SetPtEtaPhiM(electrons[1]->pt, electrons[1]->eta, electrons[1]->phi, 5.11e-6);
        dielectronP4 = electronOneP4 + electronTwoP4;
        if ( fabs(dielectronP4.M() - 91) > 15 || electrons[0]->q == electrons[1]->q)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        leptonOneP4     = electronOneP4;
        leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonOneFlavor = electrons[0]->q*11;
        leptonOneDZ     = electrons[0]->dz;
        leptonOneD0     = electrons[0]->d0;

        leptonTwoP4     = electronTwoP4;
        leptonTwoIso    = GetElectronIsolation(electrons[1], fInfo->rhoJet);
        leptonTwoFlavor = electrons[1]->q*11;
        leptonTwoDZ     = electrons[1]->dz;
        leptonTwoD0     = electrons[1]->d0;

        TLorentzVector electronThreeP4;
        electronThreeP4.SetPtEtaPhiM(electrons[2]->pt, electrons[2]->eta, electrons[2]->phi, 5.11e-6);
        leptonThreeP4     = electronThreeP4;
        leptonThreeIso    = GetElectronIsolation(electrons[2], fInfo->rhoJet);
        leptonThreeFlavor = electrons[2]->q*11;
        leptonThreeDZ     = electrons[2]->dz;
        leptonThreeD0     = electrons[2]->d0;

        if (particleSelector->PassElectronIso(electrons[2], cuts->tightElIso, cuts->EAEl) )
            leptonThreeIsoPass = 1;
        else
            leptonThreeIsoPass = 0;

        if (!isData) {

            EfficiencyContainer effCont1, effCont2;
            effCont1 = weights->GetElectronRecoEff(leptonOneP4);
            effCont2 = weights->GetElectronRecoEff(leptonTwoP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            float leptonOneRecoWeight = effs.first/effs.second;

            effs = effCont2.GetEff();
            errs = effCont2.GetErr();
            float leptonTwoRecoWeight = effs.first/effs.second;

            // trigger weights with trigger matching
            bitset<2> triggered;
            for (const auto& name: passTriggerNames) {
                if (trigger->passObj(name, 1, electrons[0]->hltMatchBits) && leptonOneP4.Pt() > 30)
                    triggered.set(0);
                if (trigger->passObj(name, 1, electrons[1]->hltMatchBits) && leptonTwoP4.Pt() > 30)
                    triggered.set(1);
            }

            float triggerWeight = 1.;
            if (triggered.all()) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", leptonOneP4);
                effCont2      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", leptonTwoP4);
                triggerWeight = GetTriggerSF(effCont1, effCont2);
            } else if (triggered.test(0)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", leptonOneP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
            } else if (triggered.test(1)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", leptonTwoP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
            } else {
                return kTRUE;
                
            }

            // update event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;
        }
    }  
    
    else if (electrons.size() == 1 && muons.size() == 1 && taus.size() == 1 ) { // emu+tau selection
        channel = "emu_tau";
        eventCounts[channel]->Fill(1);

        // pick muon pair that has a mass closest to the Z pole

        if (   !muonTriggered
            || muons[0]->pt <25 || GetMuonIsolation(muons[0])/muons[0]->pt > 0.15 
            || electrons[0]->pt < 15. || ! particleSelector->PassElectronIso(electrons[0], cuts->tightElIso, cuts->EAEl)
            )
        {
            return kTRUE;
        }
        eventCounts[channel]->Fill(2);

        // apply preselection and set tree vars
        TLorentzVector muonOneP4, electronOneP4;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        electronOneP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 5.11e-6);

        if (  muons[0]->q == electrons[0]->q)
            return kTRUE;
        eventCounts[channel]->Fill(3);


        leptonOneP4     = muonOneP4;
        leptonOneIso    = GetMuonIsolation(muons[0]);
        leptonOneFlavor = muons[0]->q*13;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;

        leptonTwoP4     = electronOneP4;
        leptonTwoIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonTwoFlavor = electrons[0]->q*11;
        leptonTwoDZ     = electrons[0]->dz;
        leptonTwoD0     = electrons[0]->d0;

        TLorentzVector tauOneP4;
        tauOneP4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, 1.77682);
        int tauISO = 0;

        if (taus[0]->hpsDisc & baconhep::kByLooseIsolationMVA3newDMwLT)  tauISO+=1;
        if (taus[0]->hpsDisc & baconhep::kByMediumIsolationMVA3newDMwLT) tauISO+=1;
        if (taus[0]->hpsDisc & baconhep::kByTightIsolationMVA3newDMwLT)  tauISO+=1;
        if (taus[0]->hpsDisc & baconhep::kByVTightIsolationMVA3newDMwLT) tauISO+=1;

        tauDecayMode  = taus[0]->decaymode;

        tauMVA        = taus[0]->rawIsoMVA3newDMwLT;
        tauMVAOld     = taus[0]->rawIsoMVA3oldDMwLT;
        taupuppiChHadIso = taus[0]->puppiChHadIso;
        taupuppiGammaIso = taus[0]->puppiGammaIso;
        taupuppiNeuHadIso = taus[0]->puppiNeuHadIso;
        taupuppiChHadIsoNoLep = taus[0]->puppiChHadIsoNoLep;
        taupuppiGammaIsoNoLep = taus[0]->puppiGammaIsoNoLep;
        taupuppiNeuHadIsoNoLep = taus[0]->puppiNeuHadIsoNoLep;


        leptonThreeP4     = tauOneP4;
        leptonThreeIso    = tauISO;
        leptonThreeFlavor = taus[0]->q*15;

        if (!isData) {
            // reconstruction weights
            EfficiencyContainer effCont1, effCont2;
            effCont1 = weights->GetMuonRecoEff(muonOneP4);
            effCont2 = weights->GetElectronRecoEff(electronOneP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            float leptonOneRecoWeight = effs.first/effs.second;

            effs = effCont2.GetEff();
            errs = effCont2.GetErr();
            float leptonTwoRecoWeight = effs.first/effs.second;



            float triggerWeight = 1.;
            effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
            effs          = effCont1.GetEff();
            errs          = effCont1.GetErr();
            triggerWeight = effs.first/effs.second;
            // triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));


            // update the event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;
        }


    } 

    
    
    else {
        return kTRUE;
    }

    outFile->cd(channel.c_str());
    outTrees[channel]->Fill();
    this->passedEvents++;
    return kTRUE;

}

void FakeSelector::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void FakeSelector::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void FakeSelector::ReportPostTerminate()
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
    std::unique_ptr<FakeSelector> selector(new FakeSelector());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

float FakeSelector::GetMuonIsolation(const baconhep::TMuon* mu)
{
    float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    return combIso;
}

float FakeSelector::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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


float FakeSelector::GetTriggerSF(EfficiencyContainer eff1, EfficiencyContainer eff2)
{
    pair<double, double> trigEff1, trigEff2;
    trigEff1 = eff1.GetEff();
    trigEff2 = eff2.GetEff();

    float sf = 1.;
    if (trigEff1.second > 0 || trigEff2.second > 0) {
        sf = (1 - (1 - trigEff1.first)*(1 - trigEff2.first))/(1 - (1 - trigEff1.second)*(1 - trigEff2.second));
    }

    return sf;
}

float FakeSelector::GetTriggerSFError(EfficiencyContainer eff1, EfficiencyContainer eff2)
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
