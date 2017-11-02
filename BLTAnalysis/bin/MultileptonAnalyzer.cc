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

    if (params->selection == "4l") {
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
        triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
        //triggerNames.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
        //triggerNames.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
        //triggerNames.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*");
        //triggerNames.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*");
    } else if (params->selection == "mumu" || params->selection == "mu4j") {
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");
    } else if (params->selection == "ee") {
        triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");
    } else if (params->selection == "emu") {
        //triggerNames.push_back("HLT_IsoMu24_v*");
        //triggerNames.push_back("HLT_IsoTkMu24_v*");
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
    string outTreeName = params->get_output_treename("tree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "bltTree");

    // event data
    outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber, "eventNumber/l");
    outTree->Branch("lumiSection", &lumiSection);
    outTree->Branch("triggerStatus", &triggerStatus);
    outTree->Branch("eventWeight", &eventWeight);
    outTree->Branch("puWeight", &puWeight);
    outTree->Branch("triggerWeight", &triggerWeight);
    outTree->Branch("nPV", &nPV);
    outTree->Branch("nPU", &nPU);
    outTree->Branch("nPartons", &nPartons);
    outTree->Branch("rPV", &rPV);

    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);
    outTree->Branch("ht", &ht);
    outTree->Branch("htPhi", &htPhi);
    outTree->Branch("htSum", &htSum);

    // leptons
    outTree->Branch("leptonOneP4", &leptonOneP4);
    outTree->Branch("leptonOneIso", &leptonOneIso);
    outTree->Branch("leptonOneFlavor", &leptonOneFlavor);
    outTree->Branch("leptonOneMother", &leptonOneMother);
    outTree->Branch("leptonOneD0", &leptonOneD0);
    outTree->Branch("leptonOneDZ", &leptonOneDZ);
    outTree->Branch("leptonOneRecoWeight", &leptonOneRecoWeight);

    outTree->Branch("leptonTwoP4", &leptonTwoP4);
    outTree->Branch("leptonTwoIso", &leptonTwoIso);
    outTree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
    outTree->Branch("leptonTwoMother", &leptonTwoMother);
    outTree->Branch("leptonTwoD0", &leptonTwoD0);
    outTree->Branch("leptonTwoDZ", &leptonTwoDZ);
    outTree->Branch("leptonTwoRecoWeight", &leptonTwoRecoWeight);

    // lepton vertices
    outTree->Branch("dileptonVertexOne", &dileptonVertexOne);
    outTree->Branch("dileptonVertexErrOne", &dileptonVertexErrOne);
    outTree->Branch("dileptonVertexChi2One", &dileptonVertexChi2One);
    outTree->Branch("dileptonVertexDOFOne", &dileptonVertexDOFOne);

    if (params->selection == "4l") {
        outTree->Branch("leptonThreeP4", &leptonThreeP4);
        outTree->Branch("leptonThreeIso", &leptonThreeIso);
        outTree->Branch("leptonThreeFlavor", &leptonThreeFlavor);
        outTree->Branch("leptonThreeD0", &leptonThreeD0);
        outTree->Branch("leptonThreeDZ", &leptonThreeDZ);
        outTree->Branch("leptonThreeRecoWeight", &leptonThreeRecoWeight);

        outTree->Branch("leptonFourP4", &leptonFourP4);
        outTree->Branch("leptonFourIso", &leptonFourIso);
        outTree->Branch("leptonFourFlavor", &leptonFourFlavor);
        outTree->Branch("leptonFourD0", &leptonFourD0);
        outTree->Branch("leptonFourDZ", &leptonFourDZ);
        outTree->Branch("leptonFourRecoWeight", &leptonFourRecoWeight);

        outTree->Branch("dileptonVertexTwo", &dileptonVertexTwo);
        outTree->Branch("dileptonVertexErrTwo", &dileptonVertexErrTwo);
        outTree->Branch("dileptonVertexChi2Two", &dileptonVertexChi2Two);
        outTree->Branch("dileptonVertexDOFTwo", &dileptonVertexDOFTwo);
    }

    // jets
    outTree->Branch("jetOneP4", &jetOneP4);
    outTree->Branch("jetOneTag", &jetOneTag);
    outTree->Branch("jetTwoP4", &jetTwoP4);
    outTree->Branch("jetTwoTag", &jetTwoTag);

    if (params->selection == "mu4j") {
        outTree->Branch("jetThreeP4", &jetThreeP4);
        outTree->Branch("jetThreeTag", &jetThreeTag);
        outTree->Branch("jetFourP4", &jetFourP4);
        outTree->Branch("jetFourTag", &jetFourTag);
    }

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
    outTree->Branch("nTaus", &nTaus);
    outTree->Branch("nPhotons", &nPhotons);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nFwdJets", &nFwdJets);
    outTree->Branch("nBJets", &nBJets);

    // event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",30,0.5,30.5);

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
    if (!isData) {
        unsigned count = 0;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
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

            // This will save Z and J/psi in 4l events
            if (
                    params->selection == "4l"
                    && (abs(particle->pdgId) == 443 || abs(particle->pdgId) == 23) 
                    && particle->status == 62 
               ) {
                genParticles.push_back(particle);
            } 

            // This will save leptons (mostly for use with dilepton analyses)
            if ((abs(particle->pdgId) == 11 || abs(particle->pdgId) == 13) && particle->status == 1) {
                if (particle->parent != -2) {
                    TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                    //cout << i << ", " << particle->pdgId  << ", " << mother->pdgId << ", " << particle->status << endl;;
                    if (fabs(mother->pdgId) == 24 || abs(mother->pdgId) == 15) {
                        //cout << i << ", " << particle->pdgId  << ", " << mother->pdgId;
                        //cout << "\t" << particle->pt << ", " << particle->eta;
                        //cout << endl;

                        if (mother->parent != -2) {
                            genParticles.push_back(particle);
                            genMotherId.push_back(mother->pdgId);

                            //TGenParticle* gmother = (TGenParticle*) fGenParticleArr->At(mother->parent);
                            //cout << gmother->pdgId << endl;
                        }

                    }
                }
            }

            nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY
        }
        //cout << genParticles.size() << endl;

        // categorize events based on decay type
        if (genParticles.size() == 2) {
            //cout << abs(genParticles[0]->pdgId) << ", " << abs(genMotherId[0]) << ", " << abs(genParticles[1]->pdgId) << ", " << abs(genMotherId[1]) << endl; 
            if (abs(genParticles[0]->pdgId) == 11 && abs(genMotherId[0]) == 24 && abs(genParticles[1]->pdgId) == 11 && abs(genMotherId[1]) == 24) { 
                // W->e + W->e
                genCategory = 0;
                hTotalEvents->Fill(10);
            } else if (abs(genParticles[0]->pdgId) == 13 && abs(genMotherId[0]) == 24 && abs(genParticles[1]->pdgId) == 13 && abs(genMotherId[1]) == 24) { 
                // W->mu + W->mu
                genCategory = 1;
                hTotalEvents->Fill(11);
            } else if (
                    (abs(genParticles[0]->pdgId) == 11 && abs(genMotherId[0]) == 24 && abs(genParticles[1]->pdgId) == 13 && abs(genMotherId[1]) == 24)
                    || (abs(genParticles[1]->pdgId) == 11 && abs(genMotherId[1]) == 24 && abs(genParticles[0]->pdgId) == 13 && abs(genMotherId[0]) == 24)
                    ) { 
                // W->e + W->mu
                genCategory = 2;
                hTotalEvents->Fill(12);
            } else if (
                    (abs(genParticles[0]->pdgId) == 13 && abs(genMotherId[0]) == 15 && abs(genParticles[1]->pdgId) == 13 && abs(genMotherId[1]) == 24)
                    || (abs(genParticles[1]->pdgId) == 13 && abs(genMotherId[1]) == 15 && abs(genParticles[0]->pdgId) == 13 && abs(genMotherId[0]) == 24)
                    ) { 
                // tau->mu + W->mu
                genCategory = 3;
                hTotalEvents->Fill(13);
            } else if (
                    (abs(genParticles[0]->pdgId) == 11 && abs(genMotherId[0]) == 15 && abs(genParticles[1]->pdgId) == 11 && abs(genMotherId[1]) == 24)
                    || (abs(genParticles[1]->pdgId) == 11 && abs(genMotherId[1]) == 15 && abs(genParticles[0]->pdgId) == 11 && abs(genMotherId[0]) == 24)
                    ) { 
                // tau->e + W->e
                genCategory = 4;
                hTotalEvents->Fill(14);
            } else if (
                    (abs(genParticles[0]->pdgId) == 13 && abs(genMotherId[0]) == 15 && abs(genParticles[1]->pdgId) == 11 && abs(genMotherId[1]) == 24)
                    || (abs(genParticles[1]->pdgId) == 13 && abs(genMotherId[1]) == 15 && abs(genParticles[0]->pdgId) == 11 && abs(genMotherId[0]) == 24)
                    ) { 
                // tau->mu + W->e
                genCategory = 5;
                hTotalEvents->Fill(15);
            } else if (
                    (abs(genParticles[0]->pdgId) == 11 && abs(genMotherId[0]) == 15 && abs(genParticles[1]->pdgId) == 13 && abs(genMotherId[1]) == 24)
                    || (abs(genParticles[1]->pdgId) == 11 && abs(genMotherId[1]) == 15 && abs(genParticles[0]->pdgId) == 13 && abs(genMotherId[0]) == 24)
                    ) { 
                // tau->e + W->mu
                genCategory = 6;
                hTotalEvents->Fill(16);
            } else if (abs(genParticles[0]->pdgId) == 13 && abs(genMotherId[0]) == 15 && abs(genParticles[1]->pdgId) == 13 && abs(genMotherId[1]) == 15) { 
                // tau->mu + tau->mu
                genCategory = 7;
                hTotalEvents->Fill(17);
            } else if (abs(genParticles[0]->pdgId) == 11 && abs(genMotherId[0]) == 15 && abs(genParticles[1]->pdgId) == 11 && abs(genMotherId[1]) == 15) { 
                // tau->e + tau->e
                genCategory = 8;
                hTotalEvents->Fill(18);
            } else if (
                    (abs(genParticles[0]->pdgId) == 11 && abs(genMotherId[0]) == 15 && abs(genParticles[1]->pdgId) == 13 && abs(genMotherId[1]) == 15)
                    || (abs(genParticles[1]->pdgId) == 11 && abs(genMotherId[1]) == 15 && abs(genParticles[0]->pdgId) == 13 && abs(genMotherId[0]) == 15)
                    ) { // tau->e + tau->mu

                genCategory = 9;
                hTotalEvents->Fill(19);
            } else {
                genCategory = 10;
                hTotalEvents->Fill(20);
            }
        } else if (genParticles.size() == 1) {
            if (abs(genParticles[0]->pdgId) == 11 && abs(genMotherId[0]) == 24) { 
                // W->e + tau->h
                genCategory = 11;
                hTotalEvents->Fill(21);
            } else if (abs(genParticles[0]->pdgId) == 13 && abs(genMotherId[0]) == 24) { 
                // W->mu + tau->h
                genCategory = 12;
                hTotalEvents->Fill(22);
            } else if (abs(genParticles[0]->pdgId) == 11 && abs(genMotherId[0]) == 15) { 
                // tau->e + tau->h
                genCategory = 13;
                hTotalEvents->Fill(23);
            } else if (abs(genParticles[0]->pdgId) == 13 && abs(genMotherId[0]) == 15) { 
                // tau->mu + tau->h
                genCategory = 14;
                hTotalEvents->Fill(24);
            }
        } else {
            genCategory = 15;
            hTotalEvents->Fill(25);
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
        puWeight *= weights->GetPUWeight(fInfo->nPUmean); // pileup reweighting
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
                && GetMuonIsolation(muon)/muonP4.Pt() < 0.25
                ) {
                    muons.push_back(muon);
                }

        // muons for jet veto
        if (
                muonP4.Pt() > 10
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
    for (int i=0; i < fTauArr->GetEntries(); i++) {
         TTau *tau = (TTau*) fTauArr->At(i);
        assert(tau);

        TLorentzVector tauP4; 
        tauP4.SetPtEtaPhiM(tau->pt, tau->eta, tau->phi, tau->m);

        // Prevent overlap of muons and jets
        bool muOverlap = false;
        for (const auto& mu: veto_muons) {
            if (tauP4.DeltaR(mu) < 0.4) {
                muOverlap = true;
                break;
            }
        }
        bool elOverlap = false;
        for (const auto& el: veto_electrons) {
            if (tauP4.DeltaR(el) < 0.4) {
                elOverlap = true;
                break;
            }
        }

        if( // Selection adopted from Stany's tau counting selection for ttDM
                tau->pt > 18  
                && abs(tau->eta) < 2.3 
                && !muOverlap
                && !elOverlap
                && (tau->hpsDisc & baconhep::kByDecayModeFinding)
                && (tau->hpsDisc & baconhep::kByLooseCombinedIsolationDBSumPtCorr3Hits)
          ) {
            taus.push_back(tau);
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

    if (params->selection == "mu4j") {

        if (muons.size() != 1 || electrons.size() != 0) 
            return kTRUE;
        hTotalEvents->Fill(5);

        // convert to TLorentzVectors
        TLorentzVector muonP4;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        float muonIso = GetMuonIsolation(muons[0]);
        if ( !(muons[0]->typeBits & baconhep::kPOGTightMuon) || muonIso/muonP4.Pt() > 0.15 )
            return kTRUE;
        hTotalEvents->Fill(6);

        if (muons[0]->pt < 25.)
            return kTRUE;
        hTotalEvents->Fill(7);

        if (nJets + nBJets < 4 || nBJets < 1)
            return kTRUE;
        hTotalEvents->Fill(8);

        leptonOneP4     = muonP4;
        leptonOneIso    = muonIso;
        leptonOneFlavor = muons[0]->q*13;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;

        // Collect the highest pt jets in the event
        //jets = KinematicTopTag(muonP4, metP2, jets); // do this here probably
        std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);
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
            triggerWeight *= trigEff.first/trigEff.second;

            eventWeight = triggerWeight*leptonOneRecoWeight;
        }

    } if (params->selection == "mumu") {

        if (muons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

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
        hTotalEvents->Fill(6);

        if (muons[0]->pt < 25. || muons[1]->pt < 3.)
            return kTRUE;
        hTotalEvents->Fill(7);

        if (nJets + nBJets < 2)
            return kTRUE;
        hTotalEvents->Fill(8);

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
            triggerWeight *= (1 - (1 - trigEff1.first)*(1 - trigEff2.first))/(1 - (1 - trigEff1.second)*(1 - trigEff2.second));

            eventWeight = triggerWeight*leptonOneRecoWeight*leptonTwoRecoWeight;
        }
    } else if (params->selection == "ee") {

        if (electrons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        if (electrons[0]->pt < 30 || electrons[1]->pt < 10)
            return kTRUE;
        hTotalEvents->Fill(6);

        TLorentzVector electronOneP4, electronTwoP4, dielectronP4;
        electronOneP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        electronTwoP4.SetPtEtaPhiM(electrons[1]->pt, electrons[1]->eta, electrons[1]->phi, 511e-6);
        dielectronP4 = electronOneP4 + electronTwoP4;
        if (dielectronP4.M() < 12.)
            return kTRUE;
        hTotalEvents->Fill(7);

        if (nJets + nBJets < 2)
            return kTRUE;
        hTotalEvents->Fill(8);

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

            leptonOneRecoWeight *= weights->GetElectronRecoIdEff(electronOneP4);
            leptonTwoRecoWeight *= weights->GetElectronRecoIdEff(electronTwoP4);

            // trigger weight
            //pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
            //pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
            //eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);

            eventWeight = leptonOneRecoWeight*leptonTwoRecoWeight;
        }

    } else if (params->selection == "emu") {

        if (muons.size() != 1 || electrons.size() != 1)
            return kTRUE;
        hTotalEvents->Fill(5);

        // trigger matching for thresholds
        float muPtThreshold = 25.;
        if (trigger->passObj("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", 1, muons[0]->hltMatchBits)
                || trigger->passObj("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*", 1, muons[0]->hltMatchBits)
           ) {
            muPtThreshold = 10;
        } else if (trigger->passObj("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", 1, muons[0]->hltMatchBits)
                || trigger->passObj("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", 1, muons[0]->hltMatchBits)
                ) {
            muPtThreshold = 25;
        }

        float elPtThreshold = 10.;
        if (trigger->passObj("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", 2, electrons[0]->hltMatchBits)
                || trigger->passObj("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*", 2, electrons[0]->hltMatchBits)
           ) {
            elPtThreshold = 15;
        } else if (trigger->passObj("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*", 2, electrons[0]->hltMatchBits)
                || trigger->passObj("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*", 2, electrons[0]->hltMatchBits)
                ) {
            elPtThreshold = 25;
        }

        if (muons[0]->pt < muPtThreshold || electrons[0]->pt < elPtThreshold)
            return kTRUE;
        hTotalEvents->Fill(6);

        TLorentzVector muonP4, electronP4, dilepton;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        electronP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        dilepton = muonP4 + electronP4;
        if (dilepton.M() < 12)
            return kTRUE;
        hTotalEvents->Fill(7);

        if (nJets + nBJets < 2)
            return kTRUE;
        hTotalEvents->Fill(8);

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

            leptonOneRecoWeight *= weights->GetMuonIDEff(muonP4);
            leptonOneRecoWeight *= weights->GetMuonISOEff(muonP4);
            leptonTwoRecoWeight *= weights->GetElectronRecoIdEff(electronP4);

            // trigger weight
            pair<float, float> trigEff = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
            triggerWeight *= trigEff.first/trigEff.second;

            eventWeight = leptonOneRecoWeight*leptonTwoRecoWeight*triggerWeight;
        }
    } else if (params->selection == "4l") {

        if (muons.size() >= 4) {

            if (muons[0]->pt < 20. || muons[1]->pt < 10.)
                return kTRUE;
            hTotalEvents->Fill(6);

            vector<TLorentzVector> muonP4(4, TLorentzVector());
            muonP4[0].SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
            muonP4[1].SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, 0.1052);
            muonP4[2].SetPtEtaPhiM(muons[2]->pt, muons[2]->eta, muons[2]->phi, 0.1052);
            muonP4[3].SetPtEtaPhiM(muons[3]->pt, muons[3]->eta, muons[3]->phi, 0.1052);

            // This will reorder the lepton collection so that the lead two
            // muons have invariant mass closest to a Z
            vector<unsigned> bestPair = PairDileptonToZ(muonP4); 
            vector<TMuon*> sortedMuons;
            vector<TLorentzVector> sortedMuonP4;
            for (unsigned i = 0; i < bestPair.size(); ++i) {
                sortedMuons.push_back(muons[bestPair[i]]);
                sortedMuonP4.push_back(muonP4[bestPair[i]]);
            }
            muons  = sortedMuons;
            muonP4 = sortedMuonP4;
            //--------------------------------------------------//

            leptonOneP4     = muonP4[0];
            leptonOneIso    = GetMuonIsolation(muons[0]);
            leptonOneFlavor = muons[0]->q*13;
            leptonOneDZ     = muons[0]->dz;
            leptonOneD0     = muons[0]->d0;

            leptonTwoP4     = muonP4[1];
            leptonTwoIso    = GetMuonIsolation(muons[1]);
            leptonTwoFlavor = muons[1]->q*13;
            leptonTwoDZ     = muons[1]->dz;
            leptonTwoD0     = muons[1]->d0;

            leptonThreeP4     = muonP4[2];
            leptonThreeIso    = GetMuonIsolation(muons[2]);
            leptonThreeFlavor = muons[2]->q*13;
            leptonThreeDZ     = muons[2]->dz;
            leptonThreeD0     = muons[2]->d0;

            leptonFourP4     = muonP4[3];
            leptonFourIso    = GetMuonIsolation(muons[3]);
            leptonFourFlavor = muons[3]->q*13;
            leptonFourDZ     = muons[3]->dz;
            leptonFourD0     = muons[3]->d0;

            // Apply James' Z+J/psi selection
            if (
                    !(muons[0]->typeBits & baconhep::kPOGTightMuon)
                    || !(muons[1]->typeBits & baconhep::kPOGTightMuon)
                    || leptonOneIso/leptonOneP4.Pt() > 0.15 
                    || leptonTwoIso/leptonTwoP4.Pt() > 0.15 
                    || leptonThreeIso/leptonThreeP4.Pt() > 0.25
               )
                return kTRUE;
            hTotalEvents->Fill(6);

            if (!isData) {
                eventWeight *= weights->GetMuonIDEff(muonP4[0]);
                eventWeight *= weights->GetMuonISOEff(muonP4[0]);
                eventWeight *= weights->GetMuonIDEff(muonP4[1]);
                eventWeight *= weights->GetMuonISOEff(muonP4[1]);

                eventWeight *= weights->GetMuonIDEff(muonP4[2]);
                //eventWeight *= weights->GetMuonISOEff(muonP4[2]);
                eventWeight *= weights->GetMuonIDEff(muonP4[3]);
                //eventWeight *= weights->GetMuonISOEff(muonP4[3]);

                // trigger weight
                //pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
                //pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
                //eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);

            }

        } else if (muons.size() == 2 and electrons.size() == 2) {

            if (electrons[0]->pt < 25. || electrons[1]->pt < 15.)
                return kTRUE;
            hTotalEvents->Fill(6);

            TLorentzVector lepP4;
            lepP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 0.1052);
            leptonOneP4     = lepP4;
            leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
            leptonOneFlavor = electrons[0]->q*13;
            leptonOneDZ     = electrons[0]->dz;
            leptonOneD0     = electrons[0]->d0;

            lepP4.SetPtEtaPhiM(electrons[1]->pt, electrons[1]->eta, electrons[1]->phi, 0.1052);
            leptonTwoP4     = lepP4;
            leptonTwoIso    = GetElectronIsolation(electrons[1], fInfo->rhoJet);
            leptonTwoFlavor = electrons[1]->q*13;
            leptonTwoDZ     = electrons[1]->dz;
            leptonTwoD0     = electrons[1]->d0;

            lepP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
            leptonThreeP4     = lepP4;
            leptonThreeIso    = GetMuonIsolation(muons[0]);
            leptonThreeFlavor = muons[0]->q*13;
            leptonThreeDZ     = muons[0]->dz;
            leptonThreeD0     = muons[0]->d0;

            lepP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, 0.1052);
            leptonFourP4     = lepP4;
            leptonFourIso    = GetMuonIsolation(muons[1]);
            leptonFourFlavor = muons[1]->q*13;
            leptonFourDZ     = muons[1]->dz;
            leptonFourD0     = muons[1]->d0;

            if (!isData) {
                // figure this out
                //eventWeight *= weights->GetMuonIDEff(muonP4[2]);
                //eventWeight *= weights->GetMuonISOEff(muonP4[2]);
                //eventWeight *= weights->GetMuonIDEff(muonP4[3]);
                //eventWeight *= weights->GetMuonISOEff(muonP4[3]);

                // trigger weight
                //pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
                //pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
                //eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);

            }

        } else {
            hTotalEvents->Fill(5);
            return kTRUE;
        }

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

    outTree->Fill();
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

//vector<TJet*> MultileptonAnalyzer::KinematicTopTag(vector<TJet*> jets, TVector2 metP2, TLorentzVector leptonP4)
//{
//    vector<TJet*> newJets;
//    float topMass = 172.4;
//    float wMass   = 80.385;
//    do {
//        for (unsigned = i = 0; i < n; ++i) {
//            if (v[i]) {
//                cout << i+1 << " ";
//            }
//            cout << endl;
//        } 
//    } while (prev_permutation(jets.being(), 
//            
//
//    return newJets;
//}
//
//// taken from https://stackoverflow.com/questions/12991758/creating-all-possible-k-combinations-of-n-items-in-c
//void comb(int N, int K)
//{
//    std::string bitmask(K, 1); // K leading 1's
//    bitmask.resize(N, 0); // N-K trailing 0's
//
//    // print integers and permute bitmask
//    vector<vector<unsigned> > combinations;
//    do {
//        vector<unsigned> thisCombination;
//        for (int i = 0; i < N; ++i) // [0..N-1] integers
//        {
//            if (bitmask[i]) thisCombination.push_back(i)
//        }
//        combination.push_back(thisCombination);
//    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
//    return combinations;
//}

