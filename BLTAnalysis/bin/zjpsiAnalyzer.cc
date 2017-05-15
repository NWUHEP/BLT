#include "zjpsiAnalyzer.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 

zjpsiAnalyzer::zjpsiAnalyzer() : BLTSelector()
{

}

zjpsiAnalyzer::~zjpsiAnalyzer()
{

}

void zjpsiAnalyzer::Begin(TTree *tree)
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

    if (params->selection == "mumu" || params->selection == "emu" || params->selection == "4l" || params->selection == "4mu") {
        //triggerNames.push_back("HLT_IsoMu22_v*");
        //triggerNames.push_back("HLT_IsoTkMu22_v*");
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");

    } else if (params->selection == "ee") {
        triggerNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
    }

    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    if (true) { // this will need to be turned off for MC
        string jsonFileName = cmssw_base + 
            "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
        //"/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt";
        lumiMask.AddJSONFile(jsonFileName);
    }

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
    outTree->Branch("nPV", &nPV);
    outTree->Branch("nPU", &nPU);
    outTree->Branch("nPartons", &nPartons);
    
    outTree->Branch("rPV", &rPV);

    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);

    // leptons
    outTree->Branch("leptonOneP4", &leptonOneP4);
    outTree->Branch("leptonOneIso", &leptonOneIso);
    outTree->Branch("leptonOneQ", &leptonOneQ);
    outTree->Branch("leptonOneFlavor", &leptonOneFlavor);
    outTree->Branch("leptonOneTrigger", &leptonOneTrigger);

    outTree->Branch("leptonOneD0", &leptonOneD0);
    outTree->Branch("leptonOneDZ", &leptonOneDZ);

    outTree->Branch("leptonTwoP4", &leptonTwoP4);
    outTree->Branch("leptonTwoIso", &leptonTwoIso);
    outTree->Branch("leptonTwoQ", &leptonTwoQ);
    outTree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
    outTree->Branch("leptonTwoTrigger", &leptonTwoTrigger);
    
    outTree->Branch("leptonTwoD0", &leptonTwoD0);
    outTree->Branch("leptonTwoDZ", &leptonTwoDZ);
    
    outTree->Branch("leptonThreeP4", &leptonThreeP4);
    outTree->Branch("leptonThreeIso", &leptonThreeIso);
    outTree->Branch("leptonThreeQ", &leptonThreeQ);
    outTree->Branch("leptonThreeFlavor", &leptonThreeFlavor);
    outTree->Branch("leptonThreeTrigger", &leptonThreeTrigger);
    outTree->Branch("leptonFourP4", &leptonFourP4);
    outTree->Branch("leptonFourIso", &leptonFourIso);
    outTree->Branch("leptonFourQ", &leptonFourQ);
    outTree->Branch("leptonFourFlavor", &leptonFourFlavor);
    outTree->Branch("leptonFourTrigger", &leptonFourTrigger);

    // jets
    outTree->Branch("jetP4", &jetP4);
    outTree->Branch("jetD0", &jetD0);
    outTree->Branch("jetTag", &jetTag);
    outTree->Branch("jetPUID", &jetPUID);
    outTree->Branch("jetFlavor", &jetFlavor);

    outTree->Branch("bjetP4", &bjetP4);
    outTree->Branch("bjetD0", &bjetD0);
    outTree->Branch("bjetTag", &bjetTag);
    outTree->Branch("bjetPUID", &bjetPUID);
    outTree->Branch("bjetFlavor", &bjetFlavor);

    outTree->Branch("genBJetP4", &genBJetP4);
    outTree->Branch("genBJetTag", &genBJetTag);
    outTree->Branch("genJetP4", &genJetP4);
    outTree->Branch("genJetTag", &genJetTag);

    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nFwdJets", &nFwdJets);
    outTree->Branch("nBJets", &nBJets);

    // event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);

    ReportPostBegin();
}

Bool_t zjpsiAnalyzer::Process(Long64_t entry)
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

    // Apply lumi mask
    if (isData) {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(2);

    /* Trigger selection */
    bool passTrigger = false;
    for (unsigned i = 0; i < triggerNames.size(); ++i) {
        passTrigger |= trigger->pass(triggerNames[i], fInfo->triggerBits);
    }

    if (!passTrigger && isData)
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
        eventWeight *= weights->GetPUWeight(fInfo->nPUmean); // pileup reweighting
    } else {
        nPU = 0;
    }

    ///////////////////////
    // Generator objects //
    ///////////////////////

    if (!isData) {
        unsigned count = 0;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            //cout << particle->status << ", "
            //     << particle->pdgId  << ", "
            //     << particle->parent
            //     << endl;

            if (
                    particle->status == 23 
                    && (abs(particle->pdgId) < 6 || particle->pdgId == 21) 
                    && particle->parent != -2
               ) {
                ++count;

            }
        }
        nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY
        //cout << nPartons << "\n" << endl;
    } else {
        nPartons = 0;
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
    vector<TMuon*> muons_type1;
    vector<TMuon*> muons_type2;
    vector<TMuon*> muons_type3;

    /* Apply Rochester Muon Corrections */
    for (int i=0; i < fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);
        
        double muonSF = 1.;
        if (isData) {
            muonSF = muonCorr->kScaleDT(muon->q, muon->pt, muon->eta, muon->phi, 0, 0);
        } 
        else {
            muonSF = muonCorr->kScaleAndSmearMC(muon->q, muon->pt, muon->eta, muon->phi,
                                                muon->nTkLayers, rng->Rndm(), rng->Rndm(), 
                                                0, 0);
        }
        muon->pt = muonSF*muon->pt;

        /* pfiso variable */
        float pfiso = muon->chHadIso + max(0., muon->neuHadIso + muon->gammaIso -0.5*(muon->puIso));

        /* Type 1 Muons */
        if (
                muon->pt > 20
                && fabs(muon->eta) < 2.1
                
                // tight muon ID
                && (muon->typeBits & baconhep::kPFMuon) 
                && (muon->typeBits & baconhep::kGlobal) 
                && muon->muNchi2    < 10.
                && muon->nMatchStn  > 1
                && muon->nPixHits   > 0
                && fabs(muon->d0)   < 0.2
                && fabs(muon->dz)   < 0.5
                && muon->nTkLayers  > 5 
                && muon->nValidHits > 0

                // tight PF iso
                && pfiso/muon->pt < 0.15
           ) {
                muons_type1.push_back(muon);
        }

        /* Type 2 Muons */
        else if (
                    muon->pt > 4
                    && fabs(muon->eta) < 2.1 

                    // loose muon ID
                    && (muon->typeBits & baconhep::kPFMuon)
                    && ((muon->typeBits & baconhep::kGlobal) || (muon->typeBits & baconhep::kTracker))

                    // loose PF iso
                    && pfiso/muon->pt < 0.25
                ){
                    muons_type2.push_back(muon);
        }

        /* Type 3 Muons */
        else if (
                    muon->pt > 3
                    && fabs(muon->eta) < 2.1

                    // loose muon ID
                    && (muon->typeBits & baconhep::kPFMuon)
                    && ((muon->typeBits & baconhep::kGlobal) || (muon->typeBits & baconhep::kTracker))
                ){
                    muons_type3.push_back(muon);
        }
                    
                    
    }

    sort(muons_type1.begin(), muons_type1.end(), sort_by_higher_pt<TMuon>);
    sort(muons_type2.begin(), muons_type2.end(), sort_by_higher_pt<TMuon>);
    sort(muons_type3.begin(), muons_type3.end(), sort_by_higher_pt<TMuon>);

    // Second pass
    vector<TLorentzVector> muons_P4_type1;
    vector<TLorentzVector> muons_P4_type2;
    vector<TLorentzVector> muons_P4_type3;
    
    vector<float> muons_iso_type1;
    vector<float> muons_iso_type2;
    vector<float> muons_iso_type3;
    
    vector<float> muons_q_type1;
    vector<float> muons_q_type2;
    vector<float> muons_q_type3;
    
    vector<bool> muons_trigger_type1;
    vector<bool> muons_trigger_type2;
    vector<bool> muons_trigger_type3;
    
    vector<float> muons_d0_type1;
    vector<float> muons_d0_type2;
    vector<float> muons_d0_type3;
    
    vector<float> muons_dz_type1;
    vector<float> muons_dz_type2;
    vector<float> muons_dz_type3;

    for (unsigned i = 0; i < muons_type1.size(); i++) {
        TMuon* muon = muons_type1[i];
        TLorentzVector muonP4;
        copy_p4(muons_type1[i], MUON_MASS, muonP4);

        // Fill containers

        muons_P4_type1.push_back(muonP4);
        muons_iso_type1.push_back(muon->trkIso);
        muons_q_type1.push_back(muon->q);
        muons_d0_type1.push_back(muon->d0);
        muons_dz_type1.push_back(muon->dz);

        // trigger matching
        bool triggered = false;
        for (unsigned i = 0; i < triggerNames.size(); ++i) {
            triggered |= trigger->passObj(triggerNames[i], 1, muon->hltMatchBits);
        }
        muons_trigger_type1.push_back(triggered);
    }
    std::sort(muons_P4_type1.begin(), muons_P4_type1.end(), P4SortCondition);
    
    for (unsigned i = 0; i < muons_type2.size(); i++) {
        TMuon* muon = muons_type2[i];
        TLorentzVector muonP4;
        copy_p4(muons_type2[i], MUON_MASS, muonP4);

        // Fill containers

        muons_P4_type2.push_back(muonP4);
        muons_iso_type2.push_back(muon->trkIso);
        muons_q_type2.push_back(muon->q);
        muons_d0_type2.push_back(muon->d0);
        muons_dz_type2.push_back(muon->dz);

        // trigger matching
        bool triggered = false;
        for (unsigned i = 0; i < triggerNames.size(); ++i) {
            triggered |= trigger->passObj(triggerNames[i], 1, muon->hltMatchBits);
        }
        muons_trigger_type2.push_back(triggered);
    }
    std::sort(muons_P4_type2.begin(), muons_P4_type2.end(), P4SortCondition);
    
    for (unsigned i = 0; i < muons_type3.size(); i++) {
        TMuon* muon = muons_type3[i];
        TLorentzVector muonP4;
        copy_p4(muons_type3[i], MUON_MASS, muonP4);

        // Fill containers

        muons_P4_type3.push_back(muonP4);
        muons_iso_type3.push_back(muon->trkIso);
        muons_q_type3.push_back(muon->q);
        muons_d0_type3.push_back(muon->d0);
        muons_dz_type3.push_back(muon->dz);

        // trigger matching
        bool triggered = false;
        for (unsigned i = 0; i < triggerNames.size(); ++i) {
            triggered |= trigger->passObj(triggerNames[i], 1, muon->hltMatchBits);
        }
        muons_trigger_type3.push_back(triggered);
    }
    std::sort(muons_P4_type3.begin(), muons_P4_type3.end(), P4SortCondition);

    /* ELECTRONS */
    std::vector<TLorentzVector> electrons;
    vector<float> electrons_iso;
    vector<float> electrons_q;
    vector<bool> electrons_trigger;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        if (
                electron->pt > 20 
                && fabs(electron->eta) < 2.5
                && particleSelector->PassElectronID(electron, cuts->tightElID)
                && particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
           ) {
            TLorentzVector electronP4;
            copy_p4(electron, ELE_MASS, electronP4);
            electrons.push_back(electronP4);
            electrons_iso.push_back(0.);
            electrons_q.push_back(electron->q);

            // trigger matching
            bool triggered = false;
            for (unsigned i = 0; i < triggerNames.size(); ++i) {
                triggered |= trigger->passObj(triggerNames[i], 1, electron->hltMatchBits);
            }
            electrons_trigger.push_back(triggered);
        }
    }

    std::sort(electrons.begin(), electrons.end(), P4SortCondition);


    /* JETS */
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets;
    std::vector<TJet*> fwdjets;
    std::vector<TJet*> bjets;
    std::vector<TJet*> genbjets;
    std::vector<TJet*> genjets;
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
       /* TLorentzVector vJet; 
        vJet.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
        bool muOverlap = false;
        for (const auto& mu: veto_muons) {
            if (vJet.DeltaR(mu) < 0.5) {
                muOverlap = true;
                break;
            }
        }*/
        //bool elOverlap = false;
        //for (const auto& el: electrons) {
        //    if (vJet.DeltaR(el) < 0.5) {
        //        elOverlap = true;
        //        break;
        //    }
        //}

        //cout << jet->genpt << ", " << jet->genm << ", " 
        //     << jet->geneta << ", " << jet->genphi << ", " 
        //     << jet->partonFlavor << ", " << jet->hadronFlavor << ", " 
        //     << endl;


        if (!isData) {
            if (abs(jet->hadronFlavor) == 5) {
                genbjets.push_back(jet);
            } else {
                genjets.push_back(jet);
            }
        }

        if (
                jet->pt > 30 
                && fabs(jet->eta) < 4.7
                && particleSelector->PassJetID(jet, cuts->looseJetID)
           ) {

            if (fabs(jet->eta) <= 2.4) { 
                if (
                        jet->pt > 30 
                        && jet->mva > -0.89
                        //&& !muOverlap 
                        //&& !elOverlap
                   ) { 
                    if (isData) {
                        if (jet->bmva > 0.9432) { 
                            bjets.push_back(jet);
                            ++nBJets;
                        } else {
                            jets.push_back(jet);
                            ++nJets;
                        }
                    } else {
                        if (particleSelector->BTagModifier(jet, "MVAT")) { 
                            bjets.push_back(jet);
                            ++nBJets;
                        } else {
                            jets.push_back(jet);
                            ++nJets;
                        }
                    }
                }
            } else {
                if ((fabs(jet->eta) < 2.5 && jet->mva > -0.89) || fabs(jet->eta) > 2.5) {
                    fwdjets.push_back(jet);
                    ++nFwdJets;
                }
            }
        }
    }

    std::sort(fwdjets.begin(), fwdjets.end(), sort_by_higher_pt<TJet>);
    std::sort(bjets.begin(), bjets.end(), sort_by_higher_pt<TJet>);
    std::sort(genjets.begin(), genjets.end(), sort_by_higher_pt<TJet>);
    std::sort(genbjets.begin(), genbjets.end(), sort_by_higher_pt<TJet>);

    // Add additional b jets to the central jet collection
    if (bjets.size() > 1) {
        jets.insert(jets.end(), bjets.begin()+1, bjets.end());
    }
    std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);

    /* MET */
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;

    if (!isData) {
        //met = MetKluge(met)*met;
        met = 0.96*met;
    }


    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    //nMuons     = muons.size();
    nMuons     = muons_P4_type1.size() + muons_P4_type2.size() + muons_P4_type3.size();
    nElectrons = electrons.size();

    if (params->selection == "mumu") {
        if (muons_P4_type1.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        // Find leading positive and negatively charged muons and convert to TLorentzVectors
        TLorentzVector muonOneP4, muonTwoP4;
        unsigned muonTwoIndex = 1;
        muonOneP4 = muons_P4_type1[0];
        for (unsigned i = 1; i < muons_P4_type1.size(); ++i) {
            if (muons_q_type1[0] != muons_q_type1[i]) {
                muonTwoP4 = muons_P4_type1[i];
                muonTwoIndex = i;
                break;
            }
        }

        if (muonOneP4.Pt() < 25)
            return kTRUE;
        hTotalEvents->Fill(6);

        leptonOneP4      = muonOneP4;
        leptonOneIso     = muons_iso_type1[0];
        leptonOneQ       = muons_q_type1[0];
        leptonOneTrigger = muons_trigger_type1[0];
        leptonOneFlavor  = 1;

        leptonOneD0        = muons_d0_type1[0];
        leptonOneDZ        = muons_dz_type1[0];
        
        leptonTwoP4      = muonTwoP4;
        leptonTwoIso     = muons_iso_type1[muonTwoIndex];
        leptonTwoQ       = muons_q_type1[muonTwoIndex];
        leptonTwoTrigger = muons_trigger_type1[muonTwoIndex];
        leptonTwoFlavor  = 1;
            
        leptonTwoD0        = muons_d0_type1[muonTwoIndex];
        leptonTwoDZ        = muons_dz_type1[muonTwoIndex];

        if (!isData) {
            eventWeight *= weights->GetMuonIDEff(muonOneP4);
            eventWeight *= weights->GetMuonISOEff(muonOneP4);
            eventWeight *= weights->GetMuonIDEff(muonTwoP4);
            eventWeight *= weights->GetMuonISOEff(muonTwoP4);

            // trigger weight
            pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
            pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
            eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);
        }

    }

    else if (params->selection == "4mu") {
        if (nMuons < 4)
            return kTRUE;
        hTotalEvents->Fill(5);
        if (muons_P4_type1.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(6);

        // Identify the Z candidate muons
        TLorentzVector muonOneP4, muonTwoP4;
        unsigned muonTwoIndex = 1;
        muonOneP4 = muons_P4_type1[0];
        for (unsigned i = 1; i < muons_P4_type1.size(); ++i) {
            if (muons_q_type1[0] != muons_q_type1[i]) {
                muonTwoP4 = muons_P4_type1[i];
                muonTwoIndex = i;
                break;
            }
        }

        if (muonOneP4.Pt() < 25)
            return kTRUE;
        hTotalEvents->Fill(7);
        
        leptonOneP4      = muonOneP4;
        leptonOneIso     = muons_iso_type1[0];
        leptonOneQ       = muons_q_type1[0];
        leptonOneTrigger = muons_trigger_type1[0];
        leptonOneFlavor  = 1;

        leptonOneD0        = muons_d0_type1[0];
        leptonOneDZ        = muons_dz_type1[0];
        
        leptonTwoP4      = muonTwoP4;
        leptonTwoIso     = muons_iso_type1[muonTwoIndex];
        leptonTwoQ       = muons_q_type1[muonTwoIndex];
        leptonTwoTrigger = muons_trigger_type1[muonTwoIndex];
        leptonTwoFlavor  = 1;
            
        leptonTwoD0        = muons_d0_type1[muonTwoIndex];
        leptonTwoDZ        = muons_dz_type1[muonTwoIndex];

        // Identify the J/Psi candidate muons
        // Separate cases depending on the number of type 1, 2, and 3 muons
        if (muons_P4_type1.size() >= 4) {
            TLorentzVector muonThreeP4, muonFourP4;
            unsigned muonThreeIndex = 1;
            unsigned muonFourIndex = 1;
            for (unsigned i = 1; i < muons_P4_type1.size(); ++i) {
                if (i != muonTwoIndex) { 
                muonThreeP4 = muons_P4_type1[i];
                muonThreeIndex = i;
                break;
                }
            }
            for (unsigned i = 1; i < muons_P4_type1.size(); ++i) {
                if (i != muonTwoIndex && i != muonThreeIndex) { 
                muonFourP4 = muons_P4_type1[i];
                muonFourIndex = i;
                break;
                }
            }
            leptonThreeP4      = muonThreeP4;
            leptonThreeIso     = muons_iso_type1[muonThreeIndex];
            leptonThreeQ       = muons_q_type1[muonThreeIndex];
            leptonThreeTrigger = muons_trigger_type1[muonThreeIndex];
            leptonThreeFlavor  = 1;

            //leptonThreeD0        = muons_d0_type1[muonThreeIndex];
            //leptonThreeDZ        = muons_dz_type1[muonThreeIndex];
            
            leptonFourP4      = muonFourP4;
            leptonFourIso     = muons_iso_type1[muonFourIndex];
            leptonFourQ       = muons_q_type1[muonFourIndex];
            leptonFourTrigger = muons_trigger_type1[muonFourIndex];
            leptonFourFlavor  = 1;
                
            //leptonFourD0        = muons_d0_type1[muonFourIndex];
            //leptonFourDZ        = muons_dz_type1[muonFourIndex];
        }

        else if (muons_P4_type1.size() == 3) {
            TLorentzVector muonThreeP4;
            unsigned muonThreeIndex = 1;
            for (unsigned i = 1; i < muons_P4_type1.size(); ++i) {
                if (i != muonTwoIndex) { 
                muonThreeP4 = muons_P4_type1[i];
                muonThreeIndex = i;
                break;
                }
            }
            leptonThreeP4      = muonThreeP4;
            leptonThreeIso     = muons_iso_type1[muonThreeIndex];
            leptonThreeQ       = muons_q_type1[muonThreeIndex];
            leptonThreeTrigger = muons_trigger_type1[muonThreeIndex];
            leptonThreeFlavor  = 1;

            //leptonThreeD0        = muons_d0_type1[muonThreeIndex];
            //leptonThreeDZ        = muons_dz_type1[muonThreeIndex];

            if (muons_P4_type2.size() >= 1) {
                leptonFourP4      = muons_P4_type2[0];
                leptonFourIso     = muons_iso_type2[0];
                leptonFourQ       = muons_q_type2[0];
                leptonFourTrigger = muons_trigger_type2[0];
                leptonFourFlavor  = 1;
                    
                //leptonFourD0        = muons_d0_type2[0];
                //leptonFourDZ        = muons_dz_type2[0];
            }
            else {
                leptonFourP4      = muons_P4_type3[0];
                leptonFourIso     = muons_iso_type3[0];
                leptonFourQ       = muons_q_type3[0];
                leptonFourTrigger = muons_trigger_type3[0];
                leptonFourFlavor  = 1;
                    
                //leptonFourD0        = muons_d0_type3[0];
                //leptonFourDZ        = muons_dz_type3[0];
            }
        }

        else {
            if (muons_P4_type2.size() == 0)
                return kTRUE;
            else if (muons_P4_type2.size() == 1) {
                leptonThreeP4      = muons_P4_type2[0];
                leptonThreeIso     = muons_iso_type2[0];
                leptonThreeQ       = muons_q_type2[0];
                leptonThreeTrigger = muons_trigger_type2[0];
                leptonThreeFlavor  = 1;

                //leptonThreeD0        = muons_d0_type2[0];
                //leptonThreeDZ        = muons_dz_type2[0];
                
                leptonFourP4      = muons_P4_type3[0];
                leptonFourIso     = muons_iso_type3[0];
                leptonFourQ       = muons_q_type3[0];
                leptonFourTrigger = muons_trigger_type3[0];
                leptonFourFlavor  = 1;
                    
                //leptonFourD0        = muons_d0_type3[0];
                //leptonFourDZ        = muons_dz_type3[0];
            }
            else {
                leptonThreeP4      = muons_P4_type2[0];
                leptonThreeIso     = muons_iso_type2[0];
                leptonThreeQ       = muons_q_type2[0];
                leptonThreeTrigger = muons_trigger_type2[0];
                leptonThreeFlavor  = 1;

                //leptonThreeD0        = muons_d0_type2[0];
                //leptonThreeDZ        = muons_dz_type2[0];
                
                leptonFourP4      = muons_P4_type2[1];
                leptonFourIso     = muons_iso_type2[1];
                leptonFourQ       = muons_q_type2[1];
                leptonFourTrigger = muons_trigger_type2[1];
                leptonFourFlavor  = 1;
                    
                //leptonFourD0        = muons_d0_type2[1];
                //leptonFourDZ        = muons_dz_type2[1];
            }
        }
    hTotalEvents->Fill(8);
    }

    else if (params->selection == "ee") {

        if (electrons.size() != 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        TLorentzVector dielectron;
        dielectron = electrons[0] + electrons[1];
        if (dielectron.M() < 12. || dielectron.M() > 70.)
            return kTRUE;
        hTotalEvents->Fill(6);

        leptonOneP4      = electrons[0];
        leptonOneIso     = electrons_iso[0];
        leptonOneQ       = electrons_q[0];
        leptonOneTrigger = electrons_trigger[0];
        leptonOneFlavor  = 0;

        leptonTwoP4      = electrons[1];
        leptonTwoIso     = electrons_iso[1];
        leptonTwoQ       = electrons_q[1];
        leptonTwoTrigger = electrons_trigger[1];
        leptonTwoFlavor  = 0;
    }

    else if (params->selection == "4l") {

        /*if (nMuons >= 4) { 

            if (muons_P4_type1[0].Pt() < 25)
                return kTRUE;
            hTotalEvents->Fill(6);

            leptonOneP4        = muons[0];
            leptonOneIso       = muons_iso[0];
            leptonOneQ         = muons_q[0];
            leptonOneTrigger   = muons_trigger[0];
            leptonOneFlavor    = 0;

            leptonOneD0        = muons_d0[0];
            leptonOneDZ        = muons_dz[0];

            leptonTwoP4        = muons[1];
            leptonTwoIso       = muons_iso[1];
            leptonTwoQ         = muons_q[1];
            leptonTwoTrigger   = muons_trigger[1];
            leptonTwoFlavor    = 0;
            
            leptonTwoD0        = muons_d0[1];
            leptonTwoDZ        = muons_dz[1];
            
            leptonThreeP4      = muons[2];
            leptonThreeIso     = muons_iso[2];
            leptonThreeQ       = muons_q[2];
            leptonThreeTrigger = muons_trigger[2];
            leptonThreeFlavor  = 0;
            leptonFourP4       = muons[3];
            leptonFourIso      = muons_iso[3];
            leptonFourQ        = muons_q[3];
            leptonFourTrigger  = muons_trigger[3];
            leptonFourFlavor   = 0;

        } else if (muons.size() >= 2 and electrons.size() >= 2) {

            if (muons[0].Pt() < 25)
                return kTRUE;
            hTotalEvents->Fill(6);
        
            leptonOneP4      = muons[0];
            leptonOneIso     = muons_iso[0];
            leptonOneQ       = muons_q[0];
            leptonOneTrigger = muons_trigger[0];
            leptonOneFlavor  = 0;
            
            leptonOneD0        = muons_d0[0];
            leptonOneDZ        = muons_dz[0];
           
            leptonTwoP4      = muons[1];
            leptonTwoIso     = muons_iso[1];
            leptonTwoQ       = muons_q[1];
            leptonTwoTrigger = muons_trigger[1];
            leptonTwoFlavor  = 0;
            
            leptonTwoD0        = muons_d0[1];
            leptonTwoDZ        = muons_dz[1];

            leptonThreeP4      = electrons[0];
            leptonThreeIso     = electrons_iso[0];
            leptonThreeQ       = electrons_q[0];
            leptonThreeTrigger = electrons_trigger[0];
            leptonThreeFlavor  = 1;
            leptonFourP4       = electrons[1];
            leptonFourIso      = electrons_iso[1];
            leptonFourQ        = electrons_q[1];
            leptonFourTrigger  = electrons_trigger[1];
            leptonFourFlavor   = 1;
        } else {
            hTotalEvents->Fill(5);
            return kTRUE;
        }*/
    }


    ///////////////////
    // Fill jet info //
    ///////////////////

    if (bjets.size() > 0) {
        bjetP4.SetPtEtaPhiM(bjets[0]->pt, bjets[0]->eta, bjets[0]->phi, bjets[0]->mass);
        bjetD0     = bjets[0]->d0;
        bjetTag    = bjets[0]->csv;
        bjetPUID   = bjets[0]->mva;
        bjetFlavor = bjets[0]->hadronFlavor;
    } else {
        bjetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        bjetD0     = 0.;
        bjetTag    = 0.;
        bjetPUID   = 0.;
        bjetFlavor = 0.;
    }

    if (fwdjets.size() > 0) {
        jetP4.SetPtEtaPhiM(fwdjets[0]->pt, fwdjets[0]->eta, fwdjets[0]->phi, fwdjets[0]->mass);
        jetD0     = fwdjets[0]->d0;
        jetTag    = 0.;
        jetPUID   = fwdjets[0]->mva;
        jetFlavor = fwdjets[0]->hadronFlavor;
    } else if (jets.size() > 0) {
        jetP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
        jetD0     = jets[0]->d0;
        jetTag    = jets[0]->csv;
        jetPUID   = jets[0]->mva;
        jetFlavor = jets[0]->hadronFlavor;
    } else {
        jetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        jetD0     = 0.;
        jetTag    = 0.;
        jetPUID   = 0.;
        jetFlavor = 0.;
    } 

    if (genbjets.size() > 0 && !isData) {
        genBJetP4.SetPtEtaPhiM(genbjets[0]->genpt, genbjets[0]->geneta, genbjets[0]->genphi, genbjets[0]->genm);
        genBJetTag = genbjets[0]->csv;
    } else {
        genBJetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        genBJetTag = 0;
    }

    if (genjets.size() > 0 && !isData) {
        genJetP4.SetPtEtaPhiM(genjets[0]->genpt, genjets[0]->geneta, genjets[0]->genphi, genjets[0]->genm);
        genJetTag = genjets[0]->csv;
    } else {
        genJetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        genJetTag = 0;
    }

    outTree->Fill();
    this->passedEvents++;
    return kTRUE;
}

void zjpsiAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void zjpsiAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void zjpsiAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<zjpsiAnalyzer> selector(new zjpsiAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

float zjpsiAnalyzer::MetKluge(float met)
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
