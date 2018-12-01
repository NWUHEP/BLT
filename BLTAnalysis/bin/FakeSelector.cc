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
    triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");
    triggerNames.push_back("HLT_IsoMu24_v*");
    triggerNames.push_back("HLT_IsoTkMu24_v*");

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

    string channel = "muonFakes";
    outFile->mkdir(channel.c_str());
    outFile->cd(channel.c_str());
    string treeName = "bltTree_" + params->datasetgroup;
    outTree = new TTree(treeName.c_str(), treeName.c_str());

    // event data
    outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber, "eventNumber/l");
    outTree->Branch("lumiSection", &lumiSection);
    outTree->Branch("triggerStatus", &triggerStatus);
    outTree->Branch("nPV", &nPV);
    outTree->Branch("nPU", &nPU);
    outTree->Branch("nPartons", &nPartons);
    outTree->Branch("rPV", &rPV);
    outTree->Branch("eventWeight", &eventWeight);

    // met and ht
    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);
    outTree->Branch("ht", &ht);
    outTree->Branch("htPhi", &htPhi);
    outTree->Branch("htSum", &htSum);

    // leptons
    outTree->Branch("leptonOneIso", &leptonOneIso);
    outTree->Branch("leptonOneP4", &leptonOneP4);
    outTree->Branch("leptonOneFlavor", &leptonOneFlavor);
    outTree->Branch("leptonOneMother", &leptonOneMother);
    outTree->Branch("leptonOneD0", &leptonOneD0);
    outTree->Branch("leptonOneDZ", &leptonOneDZ);

    outTree->Branch("leptonTwoP4", &leptonTwoP4);
    outTree->Branch("leptonTwoIso", &leptonTwoIso);
    outTree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
    outTree->Branch("leptonTwoMother", &leptonTwoMother);
    outTree->Branch("leptonTwoD0", &leptonTwoD0);
    outTree->Branch("leptonTwoDZ", &leptonTwoDZ);

    outTree->Branch("leptonThreeP4", &leptonThreeP4);
    outTree->Branch("leptonThreeIso", &leptonThreeIso);
    outTree->Branch("leptonThreeFlavor", &leptonThreeFlavor);
    outTree->Branch("leptonThreeMother", &leptonThreeMother);
    outTree->Branch("leptonThreeD0", &leptonThreeD0);
    outTree->Branch("leptonThreeDZ", &leptonThreeDZ);

    // gen level objects
    outTree->Branch("genCategory", &genCategory);
    outTree->Branch("genOneP4", &genOneP4);
    outTree->Branch("genOneId", &genOneId);
    //outTree->Branch("genOneMother", &genOneMother);
    outTree->Branch("genTwoP4", &genTwoP4);
    outTree->Branch("genTwoId", &genTwoId);
    //outTree->Branch("genTwoMother", &genTwoMother);

    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nPhotons", &nPhotons);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nFwdJets", &nFwdJets);
    outTree->Branch("nBJets", &nBJets);

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


            //cout << i  << ", " << particle->parent << ", " << particle->pdgId << ", " << particle->status;
            //cout << "\t" << particle->pt << ", " << particle->eta;
            //cout << endl;

            // parton counting for jet-binned Drell-Yan samples
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
                && particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
           ) {
            electrons.push_back(electron);
            veto_electrons.push_back(electronP4);
        }
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);

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
                    if (particleSelector->BTagModifier(jet, "MVAT", "central", rng->Uniform(1.))) { 
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

    if (muons.size() >= 3 && electrons.size() == 0) { // mumu+mu selection
        hTotalEvents->Fill(5);

        // pick muon pair that has a mass closest to the Z pole
        float ix[] = {0, 1, 2};
        float massDiff = 1e9;
        vector<TMuon*> tmp_muons = muons;
        for (unsigned i = 0; i < muons.size() - 1; ++i) {
            for (unsigned j = i+1; j < muons.size(); ++j) {

                TLorentzVector tmpOne, tmpTwo, tmpDimuon;
                tmpOne.SetPtEtaPhiM(muons[i]->pt, muons[i]->eta, muons[i]->phi, 0.1052);
                tmpTwo.SetPtEtaPhiM(muons[j]->pt, muons[j]->eta, muons[j]->phi, 0.1052);
                float muonOneIso = GetMuonIsolation(muons[i]);
                float muonTwoIso = GetMuonIsolation(muons[j]);
                tmpDimuon = tmpOne + tmpTwo;

                if (
                        muons[i]->pt > 25. 
                        && muons[j]->pt > 10.
                        && muonOneIso/tmpOne.Pt() < 0.15 
                        && muonTwoIso/tmpTwo.Pt() < 0.15
                        && bool(muons[i]->q != muons[j]->q)
                   ) {

                    float dm = fabs(tmpDimuon.M() - 91);
                    if (dm < massDiff) {
                        massDiff = dm;

                        ix[0] = i;
                        ix[1] = j;
                        tmp_muons[0] = muons[i];
                        tmp_muons[1] = muons[j];
                    }
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

        // apply preselection and set tree vars
        TLorentzVector muonOneP4, muonTwoP4, dimuonP4;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, 0.1052);
        dimuonP4 = muonOneP4 + muonTwoP4;
        float muonOneIso = GetMuonIsolation(muons[0]);
        float muonTwoIso = GetMuonIsolation(muons[1]);
        if (
                muons[0]->pt < 25. 
                || muons[1]->pt < 10.
                || muonOneIso/muonOneP4.Pt() > 0.15 
                || muonTwoIso/muonTwoP4.Pt() > 0.15 
                || fabs(dimuonP4.M() - 91) > 15
                || muons[0]->q == muons[1]->q
           )
            return kTRUE;
        hTotalEvents->Fill(6);

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
            //leptonOneMother = GetGenMotherId(genParticles, muonOneP4);
            //leptonTwoMother = GetGenMotherId(genParticles, muonTwoP4);

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
        hTotalEvents->Fill(5);

        // pick muon pair that has a mass closest to the Z pole


        // apply preselection and set tree vars
        TLorentzVector electronOneP4, electronTwoP4, dielectronP4;
        electronOneP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 5.11e-6);
        electronTwoP4.SetPtEtaPhiM(electrons[1]->pt, electrons[1]->eta, electrons[1]->phi, 5.11e-6);
        dielectronP4 = electronOneP4 + electronTwoP4;
        if (
                electrons[0]->pt < 30. 
                || electrons[1]->pt < 10.
                || fabs(dielectronP4.M() - 91) > 15
                || electrons[0]->q == electrons[1]->q
           )
            return kTRUE;
        hTotalEvents->Fill(6);

        leptonOneP4     = electronOneP4;
        leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonOneFlavor = 11*electrons[0]->q;
        leptonOneDZ     = electrons[0]->dz;
        leptonOneD0     = electrons[0]->d0;

        leptonTwoP4     = electronTwoP4;
        leptonTwoIso    = GetElectronIsolation(electrons[1], fInfo->rhoJet);
        leptonTwoFlavor = 11*electrons[1]->q;
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
    } else {
        return kTRUE;
    }


    outTree->Fill();
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

vector<unsigned> FakeSelector::PairDileptonToZ(vector<TLorentzVector> leptons)
{
    unsigned pairings[3][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}};
    float minMassDiff = 1e6;
    unsigned bestPairing[4] = {0, 1, 2, 3};
    for (unsigned i = 0; i < 3; ++i) {
        unsigned *p = pairings[i];
        TLorentzVector dilepton1 = leptons[p[0]] + leptons[p[1]];
        TLorentzVector dilepton2 = leptons[p[2]] + leptons[p[3]];

        if (fabs(dilepton1.M() - 91.2) < minMassDiff) {
            minMassDiff = fabs(dilepton1.M() - 91.2);

            // This makes me loathe c++
            bestPairing[0] = p[0];
            bestPairing[1] = p[1];
            bestPairing[2] = p[2];
            bestPairing[3] = p[3];
        }
        if (fabs(dilepton2.M() - 91.2) < minMassDiff) {
            minMassDiff = fabs(dilepton2.M() - 91.2);
            bestPairing[0] = p[2];
            bestPairing[1] = p[3];
            bestPairing[2] = p[0];
            bestPairing[3] = p[1];
        }
    }

    // this too
    vector<unsigned> outPairing(bestPairing, bestPairing + sizeof(bestPairing)/sizeof(unsigned));
    return outPairing;
}

int FakeSelector::GetGenMotherId(vector<TGenParticle*> particles, TLorentzVector p4)
{
    int motherId = 0;
    for (unsigned i = 0; i < particles.size(); ++i) {
        TLorentzVector genP4;
        genP4.SetPtEtaPhiM(particles[i]->pt, particles[i]->eta, particles[i]->phi, particles[i]->mass); 
        if (genP4.DeltaR(p4) < 0.3) {
            TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particles[i]->parent);
            motherId = mother->pdgId;
        }
    }
    return motherId;
}

vector<TJet*> FakeSelector::KinematicTopTag(vector<TJet*> jets, TVector2 metP2, TLorentzVector leptonP4)
{
    float topMass = 172.4;
    float wMass   = 80.385;
    float minChi2 = 1e9;

    vector<int> ix(jets.size());
    std::iota(std::begin(ix), std::end(ix), 0);
    vector<int> optimalOrder = ix;
    do {
        if (ix[0] > ix[1]) continue;

        vector<TLorentzVector> jetP4;
        for (unsigned i = 0; i < 4; ++i) {
            TLorentzVector p;
            p.SetPtEtaPhiM(jets[ix[i]]->pt, jets[ix[i]]->eta, jets[ix[i]]->phi, jets[ix[i]]->mass);
            jetP4.push_back(p);
            //cout << ix[i];
        }

        TLorentzVector wCandidate    = jetP4[0] + jetP4[1];
        TLorentzVector hTopCandidate = wCandidate + jetP4[2];
        TLorentzVector lTopCandidate = leptonP4 + jetP4[3];

        // hadronic top mass tag
        float chi2 = pow(wCandidate.M() - wMass, 2) + pow(hTopCandidate.M() - topMass, 2);

        // leptonic top mass tag

        // b jet assignment

        // angular correlations

        if (chi2 < minChi2) {
            minChi2 = chi2;
            optimalOrder = ix;
        }
        //cout << " " << chi2 << endl;
        std::reverse(ix.begin() + 4, ix.end());
    } while (next_permutation(ix.begin(), ix.end()));
    //cout << minChi2 << "\n\n";

    vector<TJet*> newJets;
    for (unsigned i = 0; i < 4; ++i) {
        newJets.push_back(jets[optimalOrder[i]]);
    }
    return newJets;
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
