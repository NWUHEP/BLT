#include "hzgAnalyzer.h"
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


hzgAnalyzer::hzgAnalyzer() : BLTSelector()
{

}

hzgAnalyzer::~hzgAnalyzer()
{

}

void hzgAnalyzer::Begin(TTree *tree)
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

    if (params->selection == "mumu" || params->selection == "mumug") {
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
    }
    else if (params->selection == "elelg") {
        triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
        triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
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
    outTree->Branch("genWeight", &genWeight);
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

    // photons
    outTree->Branch("photonOneP4", &photonOneP4);
    outTree->Branch("photonOneR9", &photonOneR9);

    // lepton vertices
    outTree->Branch("dileptonVertexOne", &dileptonVertexOne);
    outTree->Branch("dileptonVertexErrOne", &dileptonVertexErrOne);
    outTree->Branch("dileptonVertexChi2One", &dileptonVertexChi2One);
    outTree->Branch("dileptonVertexDOFOne", &dileptonVertexDOFOne);

    // jets
    outTree->Branch("jetOneP4", &jetOneP4);
    outTree->Branch("jetOneTag", &jetOneTag);
    outTree->Branch("jetTwoP4", &jetTwoP4);
    outTree->Branch("jetTwoTag", &jetTwoTag);

//    if (params->selection == "mu4j") {
//        outTree->Branch("jetThreeP4", &jetThreeP4);
//        outTree->Branch("jetThreeTag", &jetThreeTag);
//        outTree->Branch("jetFourP4", &jetFourP4);
//        outTree->Branch("jetFourTag", &jetFourTag);
//    }

    // gen level objects
    outTree->Branch("genCategory", &genCategory);
    outTree->Branch("genOneP4", &genOneP4);
    outTree->Branch("genOneId", &genOneId);
    //outTree->Branch("genOneMother", &genOneMother);
    outTree->Branch("genTwoP4", &genTwoP4);
    outTree->Branch("genTwoId", &genTwoId);
    //outTree->Branch("genTwoMother", &genTwoMother);
    outTree->Branch("fromHardProcessFinalState", &fromHardProcessFinalState);
    outTree->Branch("isPromptFinalState", &isPromptFinalState);
    outTree->Branch("hasPhotonMatch", &hasPhotonMatch);
    outTree->Branch("vetoDY", &vetoDY);

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

Bool_t hzgAnalyzer::Process(Long64_t entry)
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

    vetoDY = false; 

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

            // This will save leptons 
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

            // DY photon veto
            if (abs(particle->pdgId) == 22 && particle->parent > 0) {
                genParticles.push_back(particle);
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                if (fabs(mother->pdgId) < 22 && particle->pt > 5) {
                    TLorentzVector particleP4;
                    TLorentzVector motherP4;
                    particleP4.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, particle->mass); 
                    motherP4.SetPtEtaPhiM(mother->pt, mother->eta, mother->phi, mother->mass); 
                    if (particleP4.DeltaR(motherP4) > 0.3) 
                        vetoDY = true;
                }                    
            }

            nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY
        }
        //cout << genParticles.size() << endl;

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
        puWeight = weights->GetPUWeight(nPU); // pileup reweighting
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
               //// tight muon ID and ISO
               //muonP4.Pt() > 5.
               //&& fabs(muonP4.Eta()) < 2.4
               //&& (muon->typeBits & baconhep::kPFMuon) 
               //&& (muon->typeBits & baconhep::kGlobal) 
               //&& muon->muNchi2    < 10.
               //&& muon->nMatchStn  > 1
               //&& muon->nPixHits   > 0
               //&& fabs(muon->d0)   < 0.2
               //&& fabs(muon->dz)   < 0.5
               //&& muon->nTkLayers  > 5 
               //&& muon->nValidHits > 0
               //&& GetMuonIsolation(muon)/muonP4.Pt() < 0.25
           //){
           // muons.push_back(muon);
       // }

                // h->ZZ->4l "tight" ID and pf_isorel < 0.35
                muonP4.Pt() > 5.
                && fabs(muonP4.Eta()) < 2.4
                && fabs(muon->d0) < 0.5
                && fabs(muon->dz) < 1.0
                && ((muon->typeBits & baconhep::kGlobal) || 
                   ((muon->typeBits & baconhep::kTracker) && muon->nMatchStn > 0)) // Global muon or (arbitrated) tracker muon
                //&& !(muon->selectorBits & baconhep::kAllStandAloneMuons) // no standalone muons
                //&& !(muon->typeBits & baconhep::kStandalone) // no standalone muons
                && muon->sip3d < 4.0                 
                && GetMuonIsolation(muon)/muonP4.Pt() < 0.35

                ) {
                    // We now have h->ZZ->4l "loose" muons
                    if (muonP4.Pt() < 200.0) {
                        if (muon->pogIDBits & baconhep::kPOGLooseMuon)
                            muons.push_back(muon);
                    }
                    else {
                        // need to pass the tracker high-pt ID
                        if (
                                muon->nMatchStn > 2
                                && (muonP4.Pt()/muon->ptErr) < 0.3 
                                && fabs(muon->d0) < 0.2
                                && fabs(muon->dz) < 0.5
                                && muon->nPixHits > 0
                                && muon->nTkLayers > 5
                           )
                            muons.push_back(muon);
                    }
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
        electronP4.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, ELE_MASS);

        if (
                electron->pt > 10
                && fabs(electron->eta) < 2.5
                //&& particleSelector->PassElectronID(electron, cuts->tightElID)
                && particleSelector->PassElectronID(electron, cuts->tightMVAElID)
                //&& particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
                && GetElectronIsolation(electron, fInfo->rhoJet)/electronP4.Pt() < 0.35
                && fabs(electron->d0) < 0.5
                && fabs(electron->dz) < 1.0
                && electron->sip3d < 4.0 
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

    /* PHOTONS */
    vector <TPhoton*> photons;
    for (int i=0; i<fPhotonArr->GetEntries(); i++) {
        TPhoton* photon = (TPhoton*) fPhotonArr->At(i);
        assert(photon);

    if (
            // ID conditions
            photon->pt > 10
            && fabs(photon->eta) < 2.5 
            && particleSelector->PassPhotonMVA(photon, cuts->looseMVAPhID)
            && photon->passElectronVeto
        )
        photons.push_back(photon);
    } 
    sort(photons.begin(), photons.end(), sort_by_higher_pt<TPhoton>);

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
    nPhotons   = photons.size();

    if (params->selection == "mumu") {
        if (muons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);
        TLorentzVector muonOneP4, muonTwoP4;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, MUON_MASS);
        muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, MUON_MASS);
        if (muonOneP4.Pt() < 20.0) 
            return kTRUE;
        hTotalEvents->Fill(6);
        if (muonTwoP4.Pt() < 10.0)
            return kTRUE;
        hTotalEvents->Fill(7);
        TLorentzVector dimuon = muonOneP4 + muonTwoP4;
        if (muons[0]->q == muons[1]->q)
            return kTRUE;
        hTotalEvents->Fill(8);
        if (dimuon.M() < 80.0 || dimuon.M() > 100.0)
            return kTRUE;
        hTotalEvents->Fill(9);

        int isGlobalOne = (muons[0]->typeBits & baconhep::kGlobal) > 0;
        int isGlobalTwo = (muons[1]->typeBits & baconhep::kGlobal) > 0;
        int isTrackerOne = (muons[0]->typeBits & baconhep::kTracker) > 0;
        int isTrackerTwo = (muons[1]->typeBits & baconhep::kTracker) > 0;
     
        bool printForSync = true;
        
        if (printForSync)
            cout << runNumber << "," << lumiSection << "," << evtNumber << ","  << 
                    muonOneP4.Pt() << "," << muonOneP4.Eta() << "," << 
                    muonOneP4.Phi() << "," <<
                    isGlobalOne << "," << isTrackerOne << "," <<  
                    muons[0]->nMatchStn << "," << muons[0]->d0 << "," << 
                    muons[0]->dz << "," << muons[0]->sip3d << "," << 
                    muons[0]->nPixHits << "," << muons[0]->nTkLayers << "," << 
                    muons[0]->pt << "," << muons[0]->ptErr << "," << 
                    GetMuonIsolation(muons[0]) << "," <<
                    muonTwoP4.Pt() << "," << muonTwoP4.Eta() << "," << 
                    muonTwoP4.Phi() << "," << 
                    isGlobalTwo << "," << isTrackerTwo << "," <<
                    muons[1]->nMatchStn << "," << muons[1]->d0 << "," << 
                    muons[1]->dz << "," << muons[1]->sip3d << "," << 
                    muons[1]->nPixHits << "," << muons[1]->nTkLayers << "," << 
                    muons[1]->pt << "," << muons[1]->ptErr << "," << 
                    GetMuonIsolation(muons[1]) << endl;
        return kTRUE;
    }

    else if (params->selection == "mumug") {
        if (muons.size() < 2) 
            return kTRUE;
        hTotalEvents->Fill(5);
        if (photons.size() < 1)
            return kTRUE;
        hTotalEvents->Fill(6);
        
        unsigned int muonOneIndex = 0;
        unsigned int muonTwoIndex = 1;
        bool hasValidPair = false;
        float zMassDiff = 100.;
        for (unsigned int i = 0; i < muons.size(); ++i) {
            for (unsigned int j = i+1; j < muons.size(); ++j) {
                TLorentzVector tempMuonOne, tempMuonTwo;
                tempMuonOne.SetPtEtaPhiM(muons[i]->pt, muons[i]->eta, muons[i]->phi, MUON_MASS);
                tempMuonTwo.SetPtEtaPhiM(muons[j]->pt, muons[j]->eta, muons[j]->phi, MUON_MASS);
                float thisMass = (tempMuonOne + tempMuonTwo).M();
                if (
                        muons[i]->q != muons[j]->q
                        && muons[i]->pt > 20.0 
                        && muons[j]->pt > 10.0
                        && thisMass > 50.0
                   ) {
                    if (hasValidPair) {
                        if (fabs(thisMass - ZMASS) < zMassDiff) {
                            zMassDiff = fabs(thisMass - ZMASS);
                            leptonOneP4 = tempMuonOne;
                            leptonTwoP4 = tempMuonTwo;
                            muonOneIndex = i;
                            muonTwoIndex = j;
                        }
                    }
                    else {
                        zMassDiff = fabs(thisMass - ZMASS);
                        leptonOneP4 = tempMuonOne;
                        leptonTwoP4 = tempMuonTwo;
                        muonOneIndex = i;
                        muonTwoIndex = j;
                        hasValidPair = true;
                    }
                }
            }
        }

        if (!hasValidPair)
            return kTRUE;
        hTotalEvents->Fill(7);

        photonOneP4.SetPtEtaPhiM(photons[0]->pt, photons[0]->eta, photons[0]->phi, 0.);
        if (photonOneP4.Pt() < 15.0)
            return kTRUE;
        hTotalEvents->Fill(8);

        //float m_ll = (leptonOneP4 + leptonTwoP4).M();
        //float m_llg = (leptonOneP4 + leptonTwoP4 + photonOneP4).M();
        //if (photonOneP4.Et() < 0.14*m_llg)
        //    return kTRUE;
        //hTotalEvents->Fill(9);
        //if (leptonOneP4.DeltaR(photonOneP4) < 0.4)
        //    return kTRUE;
        //hTotalEvents->Fill(10);
        //if (leptonTwoP4.DeltaR(photonOneP4) < 0.4)
        //    return kTRUE;
        //hTotalEvents->Fill(11);
        //if (m_ll + m_llg < 185.0)
        //    return kTRUE;
        //hTotalEvents->Fill(12);
        //if (m_ll < 100.0 || m_ll > 180.0)
        //    return kTRUE;
        //hTotalEvents->Fill(13);
 
        leptonOneIso    = GetMuonIsolation(muons[muonOneIndex]);
        leptonOneFlavor = muons[muonOneIndex]->q*13;
        leptonOneDZ     = muons[muonOneIndex]->dz;
        leptonOneD0     = muons[muonOneIndex]->d0;
            
        leptonTwoIso    = GetMuonIsolation(muons[muonTwoIndex]);
        leptonTwoFlavor = muons[muonTwoIndex]->q*13;
        leptonTwoDZ     = muons[muonTwoIndex]->dz;
        leptonTwoD0     = muons[muonTwoIndex]->d0;

        photonOneR9 = photons[0]->r9;

        if (!isData) {

            float photon_gen_dr = 100.;
            TGenParticle *photon_match = new TGenParticle;
            for (unsigned int ip = 0; ip < genParticles.size(); ++ip) {
                TGenParticle *genPart = genParticles.at(ip);
                if (genPart->pdgId == 22) {
                    TLorentzVector genP4;
                    genP4.SetPtEtaPhiM(genPart->pt, genPart->eta, genPart->phi, genPart->mass); 
                    if (genP4.DeltaR(photonOneP4) < photon_gen_dr) {
                        *photon_match = *genPart;
                        photon_gen_dr = genP4.DeltaR(photonOneP4);
                    }
                }
            }
            
            if (photon_gen_dr < 0.3) {
                hasPhotonMatch            = true;
                fromHardProcessFinalState = photon_match->fromHardProcessFinalState;
                isPromptFinalState        = photon_match->isPromptFinalState;
            }
            else {
                hasPhotonMatch            = false;
                fromHardProcessFinalState = false;
                isPromptFinalState        = false;
            }

            delete photon_match;

            genWeight = fGenEvtInfo->weight;

            //eventWeight *= weights->GetMuonIDEff(leptonOneP4); // Fix for h->zz->4l id
            //eventWeight *= weights->GetMuonISOEff(leptonOneP4);
            //eventWeight *= weights->GetMuonIDEff(leptonTwoP4); // Fix for h->zz->4l id
            //eventWeight *= weights->GetMuonISOEff(leptonTwoP4);
        } 
    } // end mumug selection
    
    else if (params->selection == "elelg") {
        if (electrons.size() < 2) 
            return kTRUE;
        hTotalEvents->Fill(5);
        if (photons.size() < 1)
            return kTRUE;
        hTotalEvents->Fill(6);
        
        unsigned int electronOneIndex = 0;
        unsigned int electronTwoIndex = 1;
        bool hasValidPair = false;
        float zMassDiff = 100.;
        for (unsigned int i = 0; i < electrons.size(); ++i) {
            for (unsigned int j = i+1; j < electrons.size(); ++j) {
                TLorentzVector tempElectronOne, tempElectronTwo;
                tempElectronOne.SetPtEtaPhiM(electrons[i]->pt, electrons[i]->eta, electrons[i]->phi, ELE_MASS);
                tempElectronTwo.SetPtEtaPhiM(electrons[j]->pt, electrons[j]->eta, electrons[j]->phi, ELE_MASS);
                float thisMass = (tempElectronOne + tempElectronTwo).M();
                if (
                        electrons[i]->q != electrons[j]->q
                        && electrons[i]->pt > 25.0 
                        && electrons[j]->pt > 15.0
                        && thisMass > 50.0
                   ) {
                    if (hasValidPair) {
                        if (fabs(thisMass - ZMASS) < zMassDiff) {
                            zMassDiff = fabs(thisMass - ZMASS);
                            leptonOneP4 = tempElectronOne;
                            leptonTwoP4 = tempElectronTwo;
                            electronOneIndex = i;
                            electronTwoIndex = j;
                        }
                    }
                    else {
                        zMassDiff = fabs(thisMass - ZMASS);
                        leptonOneP4 = tempElectronOne;
                        leptonTwoP4 = tempElectronTwo;
                        electronOneIndex = i;
                        electronTwoIndex = j;
                        hasValidPair = true;
                    }
                }
            }
        }

        if (!hasValidPair)
            return kTRUE;
        hTotalEvents->Fill(7);

        photonOneP4.SetPtEtaPhiM(photons[0]->pt, photons[0]->eta, photons[0]->phi, 0.);
        if (photonOneP4.Pt() < 15.0)
            return kTRUE;
        hTotalEvents->Fill(8);

        //float m_ll = (leptonOneP4 + leptonTwoP4).M();
        //float m_llg = (leptonOneP4 + leptonTwoP4 + photonOneP4).M();
        //if (photonOneP4.Et() < 0.14*m_llg)
        //    return kTRUE;
        //hTotalEvents->Fill(9);
        //if (leptonOneP4.DeltaR(photonOneP4) < 0.4)
        //    return kTRUE;
        //hTotalEvents->Fill(10);
        //if (leptonTwoP4.DeltaR(photonOneP4) < 0.4)
        //    return kTRUE;
        //hTotalEvents->Fill(11);
        //if (m_ll + m_llg < 185.0)
        //    return kTRUE;
        //hTotalEvents->Fill(12);
        //if (m_ll < 100.0 || m_ll > 180.0)
        //    return kTRUE;
        //hTotalEvents->Fill(13);
 
        leptonOneIso    = GetElectronIsolation(electrons[electronOneIndex], fInfo->rhoJet);
        leptonOneFlavor = electrons[electronOneIndex]->q*13;
        leptonOneDZ     = electrons[electronOneIndex]->dz;
        leptonOneD0     = electrons[electronOneIndex]->d0;
            
        leptonTwoIso    = GetElectronIsolation(electrons[electronTwoIndex], fInfo->rhoJet);
        leptonTwoFlavor = electrons[electronTwoIndex]->q*13;
        leptonTwoDZ     = electrons[electronTwoIndex]->dz;
        leptonTwoD0     = electrons[electronTwoIndex]->d0;

        photonOneR9 = photons[0]->r9;

        if (!isData) {

            float photon_gen_dr = 100.;
            TGenParticle *photon_match = new TGenParticle;
            for (unsigned int ip = 0; ip < genParticles.size(); ++ip) {
                TGenParticle *genPart = genParticles.at(ip);
                if (genPart->pdgId == 22) {
                    TLorentzVector genP4;
                    genP4.SetPtEtaPhiM(genPart->pt, genPart->eta, genPart->phi, genPart->mass); 
                    if (genP4.DeltaR(photonOneP4) < photon_gen_dr) {
                        *photon_match = *genPart;
                        photon_gen_dr = genP4.DeltaR(photonOneP4);
                    }
                }
            }
            
            if (photon_gen_dr < 0.3) {
                hasPhotonMatch            = true;
                fromHardProcessFinalState = photon_match->fromHardProcessFinalState;
                isPromptFinalState        = photon_match->isPromptFinalState;
            }
            else {
                hasPhotonMatch            = false;
                fromHardProcessFinalState = false;
                isPromptFinalState        = false;
            }

            delete photon_match;

            genWeight = fGenEvtInfo->weight;

            //eventWeight *= weights->GetMuonIDEff(leptonOneP4); // Fix for h->zz->4l id
            //eventWeight *= weights->GetMuonISOEff(leptonOneP4);
            //eventWeight *= weights->GetMuonIDEff(leptonTwoP4); // Fix for h->zz->4l id
            //eventWeight *= weights->GetMuonISOEff(leptonTwoP4);
        } 
    } // end elelg selection

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

void hzgAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void hzgAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void hzgAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<hzgAnalyzer> selector(new hzgAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

float hzgAnalyzer::MetKluge(float met)
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

float hzgAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    return combIso;
}

float hzgAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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

vector<unsigned> hzgAnalyzer::PairDileptonToZ(vector<TLorentzVector> leptons)
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

int hzgAnalyzer::GetGenMotherId(vector<TGenParticle*> particles, TLorentzVector p4)
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

//vector<TJet*> hzgAnalyzer::KinematicTopTag(vector<TJet*> jets, TVector2 metP2, TLorentzVector leptonP4)
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

