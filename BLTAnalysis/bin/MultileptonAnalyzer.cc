#include "MultileptonAnalyzer.h"
#include <map>

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool sync_print = false;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 

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

    if (params->selection == "mumu" || params->selection == "emu") {
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");

    } else if (params->selection == "ee") {
        triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
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

    // MET
    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);

    // leptons
    outTree->Branch("leptonOneP4", &leptonOneP4);
    outTree->Branch("leptonTwoP4", &leptonTwoP4);
    outTree->Branch("leptonOneQ", &leptonOneQ);
    outTree->Branch("leptonTwoQ", &leptonTwoQ);

    // jets
    outTree->Branch("bjetOneP4", &bjetOneP4);
    outTree->Branch("bjetTwoP4", &bjetTwoP4);

    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nBJets", &nBJets);
    outTree->Branch("nGenElectrons",&nGenElectrons);
    outTree->Branch("nGenMuons",&nGenMuons);
    outTree->Branch("nGenTaus",&nGenTaus);
    

    // event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);

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

    if (sync_print) {
        cout << "trigger status: " << passTrigger << "\n" << endl;
    }

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
        unsigned genelectronscount = 0;
        unsigned genmuoncount = 0;
        unsigned gentaucount = 0;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
            if(particle->parent <0) continue;
            TGenParticle* particleparent = (TGenParticle*) fGenParticleArr->At(particle->parent);
            //if(particleparent->parent <0) continue;
            //TGenParticle* particlegrandparent = (TGenParticle*) fGenParticleArr->At(particleparent->parent); && abs(particlegrandparent->pdgId) == 24

            if ( abs(particle->pdgId) == 11 && abs(particleparent->pdgId) == 24 ) ++genelectronscount;
            if ( abs(particle->pdgId) == 13 && abs(particleparent->pdgId) == 24 ) ++genmuoncount;
            if ( abs(particle->pdgId) == 15 && abs(particleparent->pdgId) == 24 ) ++gentaucount;
            
        }
        nGenElectrons = genelectronscount;
        nGenMuons = genmuoncount;
        nGenTaus = gentaucount;
        cout << "elec:" << nGenElectrons << ", muon:" << nGenMuons << ", tau:"<< nGenTaus << endl;
        
    } 
    


    ///////////////////
    // Select objects//
    ///////////////////

    /* 0. Vertices */
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

    /* 1. MUONS */
    vector<TMuon*> muons;
    vector<TLorentzVector> veto_muons; // for vetoing jets that overlap with muons (not saved!!!)
    for (int i=0; i < fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);

        TLorentzVector muonP4;
        copy_p4(muon, MUON_MASS, muonP4);

        // Apply rochester muon momentum corrections
        double muonSF = 1.;
        if (isData) {
            muonSF = muonCorr->kScaleDT(muon->q, muon->pt, muon->eta, muon->phi, 0, 0);
        } else {
            muonSF = muonCorr->kScaleAndSmearMC(muon->q, muon->pt, muon->eta, muon->phi,
                    muon->nTkLayers, rng->Rndm(), rng->Rndm(), 
                    0, 0);
        }
        muonP4.SetPtEtaPhiM(muonSF*muon->pt, muon->eta, muon->phi, MUON_MASS);

        // Fill containers
        if (
                muon->pt > 10
                && fabs(muon->eta) < 2.4
                && particleSelector->PassMuonID(muon, cuts->tightMuID)
                && muon->trkIso/muon->pt < 0.1
            ) {
            muons.push_back(muon);
            //muons for jet veto
            //if (muonP4.Pt() > 20) {
            veto_muons.push_back(muonP4);
            //}   
        }
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);

    /* 2. ELECTRONS */
    std::vector<TElectron*> electrons;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        if (
                electron->pt > 10 
                && fabs(electron->eta) < 2.5
                && particleSelector->PassElectronID(electron, cuts->tightElID)
                && particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
           ) {
            TLorentzVector electronP4;
            copy_p4(electron, ELE_MASS, electronP4);
            electrons.push_back(electron);
        }
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);

    /* 3.JETS */
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets;
    std::vector<TJet*> bjets;
    nJets    = 0;
    nBJets   = 0;
    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

        // Prevent overlap of muons and jets
        TLorentzVector vJet; 
        vJet.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
        bool muOverlap = false;
        for (const auto& mu: veto_muons) {
            if (vJet.DeltaR(mu) < 0.5) {
                muOverlap = true;
                break;
            }
        }

        if (
                jet->pt > 30 
                && fabs(jet->eta) < 2.4
                && particleSelector->PassJetID(jet, cuts->looseJetID)
                && !muOverlap 
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
                if (particleSelector->BTagModifier(jet, "MVAT")) { // accounts for b jet efficiency (don't worry about this for now)
                    bjets.push_back(jet);
                    ++nBJets;
                } else {
                    jets.push_back(jet);
                    ++nJets;
                }
            }
        }
    }
    std::sort(bjets.begin(), bjets.end(), sort_by_higher_pt<TJet>);
    std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);

    /* 4.MET */
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;

    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    nMuons     = muons.size();
    nElectrons = electrons.size();

    if (params->selection == "mumu") {
        if (muons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        if (muons[0]->pt < 25 || muons[1]->pt < 20)
            return kTRUE;
        hTotalEvents->Fill(6);

        //if (muons[0]->q != muons[1]->q) // remove opposite sign muons
        //if (muons[0]->q == muons[1]->q)  // remove same sign muons
        //    return kTRUE;
        hTotalEvents->Fill(7);

        if (bjets.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(8);

        TLorentzVector muonOneP4, muonTwoP4;
        copy_p4(muons[0], MUON_MASS, muonOneP4);
        copy_p4(muons[1], MUON_MASS, muonTwoP4);
        leptonOneP4 = muonOneP4;
        leptonTwoP4 = muonTwoP4;
        leptonOneQ = muons[0]->q;
        leptonTwoQ = muons[1]->q;

        // fill b jets
        bjetOneP4.SetPtEtaPhiM(bjets[0]->pt, bjets[0]->eta, bjets[0]->phi, bjets[0]->mass);
        bjetTwoP4.SetPtEtaPhiM(bjets[1]->pt, bjets[1]->eta, bjets[1]->phi, bjets[1]->mass);

        if (!isData) {
            eventWeight *= weights->GetMuonIDEff(muonOneP4);
            eventWeight *= weights->GetMuonISOEff(muonOneP4);
            eventWeight *= weights->GetMuonIDEff(muonTwoP4);
            eventWeight *= weights->GetMuonISOEff(muonTwoP4);

            // trigger weight
            //pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
            //pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
            //eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);
        }
    } 
    if (params->selection == "ee") {
        if (electrons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        if (electrons[0]->pt < 25 || electrons[1]->pt < 20)
            return kTRUE;
        hTotalEvents->Fill(6);

        //if (electrons[0]->q == electrons[1]->q)  // remove same sign muons
        //    return kTRUE;
        hTotalEvents->Fill(7);

        if (bjets.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(8);

        TLorentzVector electronOneP4, electronTwoP4;
        copy_p4(electrons[0], ELE_MASS, electronOneP4);
        copy_p4(electrons[1], ELE_MASS, electronTwoP4);
        leptonOneP4 = electronOneP4;
        leptonTwoP4 = electronTwoP4;
        leptonOneQ = electrons[0]->q;
        leptonTwoQ = electrons[1]->q;
        // fill b jets
        bjetOneP4.SetPtEtaPhiM(bjets[0]->pt, bjets[0]->eta, bjets[0]->phi, bjets[0]->mass);
        bjetTwoP4.SetPtEtaPhiM(bjets[1]->pt, bjets[1]->eta, bjets[1]->phi, bjets[1]->mass);

        if (!isData) {
            eventWeight *= weights->GetElectronRecoIdEff(electronOneP4);
            eventWeight *= weights->GetElectronRecoIdEff(electronTwoP4);

            // trigger weight
            //pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", electronOneP4);
            //pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", electronTwoP4);
            //eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);
        }
    }
    if (params->selection == "emu") {
        if (electrons.size() < 1 || muons.size()<1 )
            return kTRUE;
        hTotalEvents->Fill(5);

        if (muons[0]->pt < 25 || electrons[0]->pt < 20) 
            return kTRUE;
        hTotalEvents->Fill(6);

        //if (electrons[0]->q == electrons[1]->q)  // remove same sign muons
        //    return kTRUE;
        hTotalEvents->Fill(7);

        if (bjets.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(8);


        TLorentzVector muonP4, electronP4;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        electronP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        leptonOneP4 = muonP4;
        leptonTwoP4 = electronP4;
        leptonOneQ = muons[0]->q;
        leptonTwoQ = electrons[0]->q;

        // fill b jets
        bjetOneP4.SetPtEtaPhiM(bjets[0]->pt, bjets[0]->eta, bjets[0]->phi, bjets[0]->mass);
        bjetTwoP4.SetPtEtaPhiM(bjets[1]->pt, bjets[1]->eta, bjets[1]->phi, bjets[1]->mass);

        if (!isData) {
            eventWeight *= weights->GetMuonIDEff(muonP4);
            eventWeight *= weights->GetElectronRecoIdEff(electronP4);
        }
    }
    
    


    outTree->Fill();
    this->passedEvents++; // controller variables
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

