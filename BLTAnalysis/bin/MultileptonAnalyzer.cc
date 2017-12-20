#include "MultileptonAnalyzer.h"
#include <map>

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

const float ZMASS = 91.19;
const float PSIMASS = 3.096;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 
bool sort_by_btag(const baconhep::TJet* lhs, const baconhep::TJet* rhs) 
{
    return lhs->bmva > rhs->bmva;
}


MultileptonAnalyzer::MultileptonAnalyzer() : BLTSelector()
{

}

MultileptonAnalyzer::~MultileptonAnalyzer()
{

}

void MultileptonAnalyzer::Begin(TTree *tree)
{

    // Parse command line option
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1),
            std::sregex_token_iterator(), std::back_inserter(options));

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

    vector<std::string> channelNames = {"mumu", "emu", "mutau", "mu4j"};
    if (params->selection == "single_muon") {
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");

    } else if (params->selection == "single electron") {
        triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");

    } else if (params->selection == "double muon") {
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
        //triggerNames.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
        //triggerNames.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
        //triggerNames.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*");
        //triggerNames.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*");
        //
    } else if (params->selection == "double electron") {
        triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");

    } else if (params->selection == "mueg") {
        triggerNames.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
        triggerNames.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
        triggerNames.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*");
        triggerNames.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*");
    }

    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
    lumiMask.AddJSONFile(jsonFileName);

    // muon momentum corrections
    muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/rcdata.2016.v3");
    rng = new TRandom3();

    // Prepare the output tree
    string outFileName = params->get_output_filename("output");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();

    for (unsigned i = 0; i < channelNames.size(); ++i) {
        string channel = channelNames[i];
        string treeName = params->get_output_treename("bltTree_" + channel);
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

        // weights and their uncertainties
        tree->Branch("eventWeight", &eventWeight);
        tree->Branch("leptonOneRecoWeight", &leptonOneRecoWeight);
        tree->Branch("leptonTwoRecoWeight", &leptonTwoRecoWeight);
        tree->Branch("topPtWeight", &topPtWeight);
        tree->Branch("puWeight", &puWeight);
        tree->Branch("triggerWeight", &triggerWeight);

        tree->Branch("met", &met);
        tree->Branch("metPhi", &metPhi);
        tree->Branch("ht", &ht);
        tree->Branch("htPhi", &htPhi);
        tree->Branch("htSum", &htSum);

        // leptons
        tree->Branch("leptonOneP4", &leptonOneP4);
        tree->Branch("leptonOneIso", &leptonOneIso);
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


        if (channel == "mutau") {
            //tree->Branch("tauChHadMult",  &tauChHadMult);
            //tree->Branch("tauPhotonMult", &tauPhotonMult);
            tree->Branch("tauDecayMode",  &tauDecayMode);
            tree->Branch("tauMVA",        &tauMVA);
        }


        // jets
        tree->Branch("jetOneP4", &jetOneP4);
        tree->Branch("jetOneTag", &jetOneTag);
        tree->Branch("jetTwoP4", &jetTwoP4);
        tree->Branch("jetTwoTag", &jetTwoTag);

        if (channel == "mu4j") {
            tree->Branch("jetThreeP4", &jetThreeP4);
            tree->Branch("jetThreeTag", &jetThreeTag);
            tree->Branch("jetFourP4", &jetFourP4);
            tree->Branch("jetFourTag", &jetFourTag);
        }

        // gen level objects
        tree->Branch("genCategory", &genCategory);
        tree->Branch("genOneP4", &genOneP4);
        tree->Branch("genOneId", &genOneId);
        //tree->Branch("genOneMother", &genOneMother);
        tree->Branch("genTwoP4", &genTwoP4);
        tree->Branch("genTwoId", &genTwoId);
        //tree->Branch("genTwoMother", &genTwoMother);

        // object counters
        tree->Branch("nMuons", &nMuons);
        tree->Branch("nElectrons", &nElectrons);
        tree->Branch("nTaus", &nTaus);
        tree->Branch("nPhotons", &nPhotons);
        tree->Branch("nJets", &nJets);
        tree->Branch("nFwdJets", &nFwdJets);
        tree->Branch("nBJets", &nBJets);

        outTrees[channel] = tree;

        // event counter
        string outHistName = params->get_output_treename("TotalEvents_" + channel);
        eventCounts[channel] = new TH1D(outHistName.c_str(),"ChannelCounts",10,0.5,10.5);
    }

    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);

    outHistName = params->get_output_treename("GenCategory");
    hGenCat = new TH1D(outHistName.c_str(), "WW decay modes",30,0.5,30.5);

    ReportPostBegin();
}

Bool_t MultileptonAnalyzer::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
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

    unsigned wCount = 0;
    unsigned tauCount = 0;
    bitset<6> wDecay;
    bitset<4> tauDecay;
    if (!isData) {
        unsigned count = 0;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            if (abs(particle->pdgId) == 6 || particle->parent == -2) continue;

            //cout << i  << ", " << particle->parent << ", " << particle->pdgId << ", " << particle->status;
            //cout << "\t" << particle->pt << ", " << particle->eta;
            //cout << endl;

            if (
                    particle->status == 23 
                    && (abs(particle->pdgId) < 6 || particle->pdgId == 21) 
                    && particle->parent != -2
               ) {
                ++count;
            }

            // Tag events based on neutrino flavors
            unsigned flavor = abs(particle->pdgId);
            if (
                    (flavor == 12 || flavor == 14 || flavor == 16)
                    && particle->parent != -2
               ) {
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                //cout << i << ", " << particle->pdgId << ", " << mother->pdgId;
                if (abs(mother->pdgId) == 24) {
                    wDecay.set(3*wCount + (flavor - 12)/2);
                    ++wCount;
                } 

                if (abs(mother->pdgId) == 15 && flavor != 16 && mother->parent != -2) {
                    TGenParticle* gmother = (TGenParticle*) fGenParticleArr->At(mother->parent);
                    if (gmother->parent == -2 || tauCount == 2) continue;
                    tauDecay.set(2*tauCount + (flavor - 12)/2);
                    ++tauCount;
                } 
            }

            // This will find final state leptons (fix this to discriminate between decays)
            if (
                    (abs(particle->pdgId) == 11 || abs(particle->pdgId) == 13) 
                    && particle->status == 1 
                    && particle->parent != -2
               ) {

                // Find if the lepton comes from a top quark
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                //cout << mother->pdgId << ", ";
                int origin = abs(mother->pdgId);
                int intermediary = origin;
                while (origin != 6 && mother->parent != -2) {
                    mother = (TGenParticle*) fGenParticleArr->At(mother->parent);
                    origin = abs(mother->pdgId);
                    if (origin <= 5) { // remove leptons that have been radiated from light quarks
                        intermediary = -1; 
                        break;
                    } else if (origin == 15) {
                        intermediary = 15;
                    } else if (intermediary != 15 && origin == 24) {
                        intermediary = 24;
                    }
                }

                if (origin == 6 && (intermediary == 24 || intermediary == 15)) {
                    genParticles.push_back(particle);
                    genMotherId.push_back(intermediary);
                }
            }
        }
        nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY
        // categorize events

        genCategory = 0;
        if (wCount == 2) { // fully leptonic
            if (wDecay.test(0) && wDecay.test(3)) {// W->e, W->e
                hGenCat->Fill(1); 
                genCategory = 1;
            } else if (wDecay.test(1) && wDecay.test(4)) {// W->mu, W->mu
                hGenCat->Fill(2); 
                genCategory = 2;
            } else if ((wDecay.test(0) && wDecay.test(4)) || (wDecay.test(1) && wDecay.test(3))) {// W->e, W->mu
                hGenCat->Fill(3); 
                genCategory = 3;
            } else if (wDecay.test(2) && wDecay.test(5)) {// W->tau, W->tau
                if (tauDecay.test(0) && tauDecay.test(2)) {// tau->e, tau->e
                    hGenCat->Fill(4); 
                    genCategory = 4;
                } else if (tauDecay.test(1) && tauDecay.test(3)) {// tau->mu, tau->mu
                    hGenCat->Fill(5); 
                    genCategory = 5;
                } else if ((tauDecay.test(1) && tauDecay.test(2)) || (tauDecay.test(0) && tauDecay.test(3))) {// tau->e, tau->mu
                    hGenCat->Fill(6); 
                    genCategory = 6;
                } else if (tauDecay.test(0) && tauCount == 1) {// tau->e, tau->h
                    hGenCat->Fill(7); 
                    genCategory = 7;
                } else if (tauDecay.test(1) && tauCount == 1) {// tau->mu, tau->h
                    hGenCat->Fill(8); 
                    genCategory = 8;
                } else if (tauCount == 0) {// tau->h, tau->h
                    hGenCat->Fill(9); 
                    genCategory = 9;
                }
            } else if ((wDecay.test(0) && wDecay.test(5)) || (wDecay.test(2) && wDecay.test(3))) {// W->e, W->tau
                if (tauDecay.test(0)) {// tau->e
                    hGenCat->Fill(10);
                    genCategory = 10;
                } else if (tauDecay.test(1)) {// tau->mu
                    hGenCat->Fill(11);
                    genCategory = 11;
                } else if (tauDecay.none()) {// tau->h
                    hGenCat->Fill(12);
                    genCategory = 12;
                }
            } else if ((wDecay.test(1) && wDecay.test(5)) || (wDecay.test(2) && wDecay.test(4))) {// W->mu, W->tau
                if (tauDecay.test(0)) {// tau->e
                    hGenCat->Fill(13);
                    genCategory = 13;
                } else if (tauDecay.test(1)) {// tau->mu
                    hGenCat->Fill(14);
                    genCategory = 14;
                } else if (tauDecay.none()) {// tau->h
                    hGenCat->Fill(15);
                    genCategory = 15;
                }
            } 
        } else if (wCount == 1) { // semileptonic
            if (wDecay.test(0)) {// W->e, W->h
                hGenCat->Fill(16); 
                genCategory = 16;
            } else if (wDecay.test(1)) {// W->mu, W->h
                hGenCat->Fill(17); 
                genCategory = 17;
            } else if (wDecay.test(2)) {// W->tau, W->h
                if (tauDecay.test(0)) {// tau->e
                    hGenCat->Fill(18);
                    genCategory = 18;
                } else if (tauDecay.test(1)) {// tau->mu
                    hGenCat->Fill(19);
                    genCategory = 19;
                } else if (tauDecay.none()) {// tau->h
                    hGenCat->Fill(20);
                    genCategory = 20;
                }
            }
        } else { // hadronic
            hGenCat->Fill(21); 
            genCategory = 21;
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

    if (!passTrigger)// && isData)
        return kTRUE;
    hTotalEvents->Fill(3);


    /////////////////////
    // Fill event info //
    /////////////////////

    eventWeight   = 1;
    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
    triggerStatus = passTrigger;
    nPV           = fPVArr->GetEntries();
    if (!isData) {
        nPU = fInfo->nPUmean;
        puWeight = weights->GetPUWeight(fInfo->nPUmean); // pileup reweighting
        eventWeight *= puWeight;
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
                muonP4.Pt() > 3.
                && fabs(muonP4.Eta()) < 2.4

                // Loose muon ID
                //&& (muon->typeBits & baconhep::kPFMuon)
                //&& ((muon->typeBits & baconhep::kGlobal) || (muon->selectorBits & baconhep::kTrackerMuonArbitrated))
                //&& fabs(muon->d0)  < 0.5
                //&& fabs(muon->dz)  < 20.0
                //&& minDeltaR > 0.02

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
                && GetMuonIsolation(muon)/muonP4.Pt() < 0.15
                ) {
                    muons.push_back(muon);
                }

        // muons for jet veto
        if (
                muonP4.Pt() > 10
                && fabs(muonP4.Eta()) < 2.4
                // tight muon ID and ISO
                && (muon->typeBits & baconhep::kPOGTightMuon)
                && GetMuonIsolation(muon)/muonP4.Pt() < 0.15
           ) {
            veto_muons.push_back(muonP4);
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

        if( // Selection adopted from Stany's tau counting selection for ttDM
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

        bool tauOverlap = false;
        for (const auto& tau: veto_taus) {
            if (jetP4.DeltaR(tau) < 0.4) {
                tauOverlap = true;
                break;
            }
        }

        if (
                jet->pt > 30 
                && fabs(jet->eta) < 4.7
                && particleSelector->PassJetID(jet, cuts->looseJetID)
                && !muOverlap 
                && !elOverlap
                && !tauOverlap
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
                    if (particleSelector->BTagModifier(jet, "MVAT")) { 
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

    string channel = "";
    if (muons.size() >= 2 && electrons.size() == 0 && taus.size() == 0) { // mu+mu selection
        channel = "mumu";

        // convert to TLorentzVectors
        TLorentzVector muonOneP4, muonTwoP4, dimuonP4;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, 0.1052);
        dimuonP4 = muonOneP4 + muonTwoP4;

        float muonOneIso = GetMuonIsolation(muons[0]);
        float muonTwoIso = GetMuonIsolation(muons[1]);
        if (
                !(muons[0]->typeBits & baconhep::kPOGTightMuon)
                || !(muons[1]->typeBits & baconhep::kPOGTightMuon)
                || muonOneIso/muonOneP4.Pt() > 0.15 
                || muonTwoIso/muonTwoP4.Pt() > 0.15 
           )
            return kTRUE;
        eventCounts[channel]->Fill(1);

        if (muons[0]->pt < 25. || muons[1]->pt < 3.)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        if (nJets + nBJets < 2 || nBJets < 1)
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
        leptonTwoDZ     = muons[1]->dz;
        leptonTwoD0     = muons[1]->d0;

        if (!isData) {
            leptonOneMother = GetGenMotherId(genParticles, muonOneP4);
            leptonTwoMother = GetGenMotherId(genParticles, muonTwoP4);

            leptonOneRecoWeight = weights->GetMuonIDEff(muonOneP4);
            leptonOneRecoWeight *= weights->GetMuonISOEff(muonOneP4);
            leptonTwoRecoWeight = weights->GetMuonIDEff(muonTwoP4);
            leptonTwoRecoWeight *= weights->GetMuonISOEff(muonTwoP4);

            // trigger weight
            pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
            pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
            triggerWeight = (1 - (1 - trigEff1.first)*(1 - trigEff2.first))/(1 - (1 - trigEff1.second)*(1 - trigEff2.second));

            eventWeight *= triggerWeight*leptonOneRecoWeight*leptonTwoRecoWeight;
        }
    }  else if (electrons.size() >= 2 && muons.size() == 0 && taus.size() == 0) { // e+e selection
        channel = "ee";

        return kTRUE;

        if (electrons[0]->pt < 30 || electrons[1]->pt < 10)
            return kTRUE;
        eventCounts[channel]->Fill(1);

        TLorentzVector electronOneP4, electronTwoP4, dielectronP4;
        electronOneP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        electronTwoP4.SetPtEtaPhiM(electrons[1]->pt, electrons[1]->eta, electrons[1]->phi, 511e-6);
        dielectronP4 = electronOneP4 + electronTwoP4;
        if (dielectronP4.M() < 12.)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        if (nJets + nBJets < 2)
            return kTRUE;
        eventCounts[channel]->Fill(3);

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

        if (!isData) {
            leptonOneMother = GetGenMotherId(genParticles, electronOneP4);
            leptonTwoMother = GetGenMotherId(genParticles, electronTwoP4);

            leptonOneRecoWeight = weights->GetElectronRecoIdEff(electronOneP4);
            leptonTwoRecoWeight = weights->GetElectronRecoIdEff(electronTwoP4);

            // trigger weight
            //pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
            //pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
            //eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);

            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;
        }
    } else if (muons.size() == 1 && electrons.size() == 1 && taus.size() == 0) { // e+mu selection
        channel = "emu";

        // trigger matching for thresholds
        float muPtThreshold = 25.;
        //if (trigger->passObj("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", 1, muons[0]->hltMatchBits)
        //        || trigger->passObj("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*", 1, muons[0]->hltMatchBits)
        //   ) {
        //    muPtThreshold = 10;
        //} else if (trigger->passObj("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", 1, muons[0]->hltMatchBits)
        //        || trigger->passObj("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", 1, muons[0]->hltMatchBits)
        //        ) {
        //    muPtThreshold = 25;
        //}

        float elPtThreshold = 10.;
        //if (trigger->passObj("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", 2, electrons[0]->hltMatchBits)
        //        || trigger->passObj("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", 2, electrons[0]->hltMatchBits)
        //   ) {
        //    elPtThreshold = 15;
        //} else if (trigger->passObj("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", 2, electrons[0]->hltMatchBits)
        //        || trigger->passObj("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*", 2, electrons[0]->hltMatchBits)
        //        ) {
        //    elPtThreshold = 25;
        //}

        if (muons[0]->pt < muPtThreshold || electrons[0]->pt < elPtThreshold)
            return kTRUE;
        eventCounts[channel]->Fill(1);

        TLorentzVector muonP4, electronP4, dilepton;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        electronP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        dilepton = muonP4 + electronP4;
        if (dilepton.M() < 12)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        if (nJets + nBJets < 2)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        leptonOneP4     = muonP4;
        leptonOneIso    = GetMuonIsolation(muons[0]);
        leptonOneFlavor = 13*muons[0]->q;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;

        leptonTwoP4     = electronP4;
        leptonTwoIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonTwoFlavor = 11*electrons[0]->q;
        leptonTwoDZ     = electrons[0]->dz;
        leptonTwoD0     = electrons[0]->d0;

        if (!isData) {
            leptonOneMother = GetGenMotherId(genParticles, muonP4);
            leptonTwoMother = GetGenMotherId(genParticles, electronP4);

            leptonOneRecoWeight = weights->GetMuonIDEff(muonP4);
            leptonOneRecoWeight *= weights->GetMuonISOEff(muonP4);
            leptonTwoRecoWeight = weights->GetElectronRecoIdEff(electronP4);
            //cout << leptonTwoRecoWeight << endl;

            // trigger weight
            pair<float, float> trigEff = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
            triggerWeight = trigEff.first/trigEff.second;

            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight*triggerWeight;

        }
    } else if (muons.size() == 1 && taus.size() == 1 && electrons.size() == 0) { // mu+tau selection
        channel = "mutau";

        if (muons[0]->pt < 25)
            return kTRUE;
        eventCounts[channel]->Fill(1);

        TLorentzVector muonP4, tauP4, dilepton;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        tauP4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, 1.776);
        dilepton = muonP4 + tauP4;
        if (dilepton.M() < 12)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        if (nJets + nBJets < 2 || nBJets < 1)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        leptonOneP4     = muonP4;
        leptonOneIso    = GetMuonIsolation(muons[0]);
        leptonOneFlavor = 13*muons[0]->q;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;

        leptonTwoP4     = tauP4;
        leptonTwoIso    = 0.;
        leptonTwoFlavor = 15*taus[0]->q;
        leptonTwoDZ     = taus[0]->dzLeadChHad;
        leptonTwoD0     = taus[0]->d0LeadChHad;

        //tauChHadMult  = taus[0]->nSignalChHad;
        //tauPhotonMult = taus[0]->nSignalGamma;
        tauDecayMode  = taus[0]->decaymode;
        tauMVA        = taus[0]->rawIsoMVA3newDMwLT;

        if (!isData) {
            leptonOneMother = GetGenMotherId(genParticles, muonP4);

            leptonOneRecoWeight = weights->GetMuonIDEff(muonP4);
            leptonOneRecoWeight *= weights->GetMuonISOEff(muonP4);
            //leptonTwoRecoWeight = weights->GetElectronRecoIdEff(electronP4);

            // trigger weight
            pair<float, float> trigEff = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
            triggerWeight = trigEff.first/trigEff.second;

            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight*triggerWeight;

        }
    } else if (muons.size() == 1 && electrons.size() == 0 && taus.size() == 0) { // mu+h selection
        channel = "mu4j";

        // convert to TLorentzVectors
        TLorentzVector muonP4;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        float muonIso = GetMuonIsolation(muons[0]);
        if (muonIso/muonP4.Pt() > 0.15 )
            return kTRUE;
        eventCounts[channel]->Fill(1);

        if (muons[0]->pt < 25.)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        if (nJets + nBJets < 4 || nBJets < 1)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        leptonOneP4     = muonP4;
        leptonOneIso    = muonIso;
        leptonOneFlavor = muons[0]->q*13;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;


        // Collect the highest pt jets in the event
        std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);
        //jets = KinematicTopTag(jets, metP2, muonP4);
        jetOneP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
        jetOneTag = jets[0]->csv;
        jetTwoP4.SetPtEtaPhiM(jets[1]->pt, jets[1]->eta, jets[1]->phi, jets[1]->mass);
        jetTwoTag = jets[1]->csv;
        jetThreeP4.SetPtEtaPhiM(jets[2]->pt, jets[2]->eta, jets[2]->phi, jets[2]->mass);
        jetThreeTag = jets[2]->csv;
        jetFourP4.SetPtEtaPhiM(jets[3]->pt, jets[3]->eta, jets[3]->phi, jets[3]->mass);
        jetFourTag = jets[3]->csv;

        if (!isData) {
            leptonOneMother = GetGenMotherId(genParticles, muonP4);

            leptonOneRecoWeight = weights->GetMuonIDEff(muonP4);
            leptonOneRecoWeight *= weights->GetMuonISOEff(muonP4);

            // trigger weight
            pair<float, float> trigEff = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
            triggerWeight = trigEff.first/trigEff.second;

            eventWeight *= triggerWeight*leptonOneRecoWeight;
        }
    } else {
        return kTRUE;
    }


    ///////////////////
    // Fill jet info //
    ///////////////////

    if (params->selection != "mu4j") { // jets are handled differently for the lepton + jet selection
        if (jets.size() > 0) {
            jetOneP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
            jetOneTag    = jets[0]->csv;
        } else {
            jetOneP4.SetPtEtaPhiM(0., 0., 0., 0.);
            jetOneTag    = 0.;
        }

        if (jets.size() > 1) {
            jetTwoP4.SetPtEtaPhiM(jets[1]->pt, jets[1]->eta, jets[1]->phi, jets[1]->mass);
            jetTwoTag    = jets[1]->csv;
        } else {
            jetTwoP4.SetPtEtaPhiM(0., 0., 0., 0.);
            jetTwoTag    = 0.;
        } 
    }

    if (!isData && genParticles.size() == 2) {
        genOneId = genParticles[0]->pdgId;
        genOneP4.SetPtEtaPhiM(genParticles[0]->pt, genParticles[0]->eta, genParticles[0]->phi, genParticles[0]->mass); 
        genTwoId = genParticles[1]->pdgId;
        genTwoP4.SetPtEtaPhiM(genParticles[1]->pt, genParticles[1]->eta, genParticles[1]->phi, genParticles[1]->mass); 
    } else {
        genOneP4.SetPtEtaPhiM(0., 0., 0., 0.); 
        genTwoP4.SetPtEtaPhiM(0., 0., 0., 0.); 
    }

    outTrees[channel]->Fill();
    this->passedEvents++;
    return kTRUE;
}

void MultileptonAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void MultileptonAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void MultileptonAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<MultileptonAnalyzer> selector(new MultileptonAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

float MultileptonAnalyzer::MetKluge(float met)
{
    if (met > 500) {
        return 1.;
    }

    float bins[]        = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 250., 500.};
    float corrections[] = {1., 0.99, 0.96, 0.95, 0.948, 0.947, 0.95, 0.957, 0.966, 0.97, 0.971, 0.98};

    int bin = 0;
    for (int i = 0; i < 12; ++i) {
        if (met > bins[i] && met <= bins[i+1]) {
            bin = i;
            break;
        }
    }

    return corrections[bin];
}

float MultileptonAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    return combIso;
}

float MultileptonAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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

vector<unsigned> MultileptonAnalyzer::PairDileptonToZ(vector<TLorentzVector> leptons)
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

int MultileptonAnalyzer::GetGenMotherId(vector<TGenParticle*> particles, TLorentzVector p4)
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

vector<TJet*> MultileptonAnalyzer::KinematicTopTag(vector<TJet*> jets, TVector2 metP2, TLorentzVector leptonP4)
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

