#include "MultileptonAnalyzer.h"
#include <map>

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool sync_print = false;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 
bool sort_by_btag(const baconhep::TJet* lhs, const baconhep::TJet* rhs) {return lhs->bmva > rhs->bmva;}

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

    if (params->selection == "mumu" || params->selection == "emu" || params->selection == "mutau" ) {
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");

    } else if (params->selection == "ee" || params->selection == "etau") {
        //triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
        triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");

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
    outTree->Branch("eventWeightPu", &eventWeightPu);
    outTree->Branch("eventWeightTrigger", &eventWeightTrigger);
    outTree->Branch("eventWeightId", &eventWeightId);
    outTree->Branch("eventWeightIso", &eventWeightIso);
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
    outTree->Branch("leptonOneMother", &leptonOneMother);
    outTree->Branch("leptonTwoMother", &leptonTwoMother);
    outTree->Branch("leptonOneD0", &leptonOneD0);
    outTree->Branch("leptonOneDZ", &leptonOneDZ);
    outTree->Branch("leptonTwoD0", &leptonTwoD0);
    outTree->Branch("leptonTwoDZ", &leptonTwoDZ);

    outTree->Branch("leptonOneSip3d", &leptonOneSip3d);
    outTree->Branch("leptonTwoSip3d", &leptonTwoSip3d);

    

    // jets
    outTree->Branch("bjetOneP4", &bjetOneP4);
    outTree->Branch("bjetTwoP4", &bjetTwoP4);
    outTree->Branch("bjetOneBmva", &bjetOneBmva);
    outTree->Branch("bjetTwoBmva", &bjetTwoBmva);

    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nBJets", &nBJets);
    outTree->Branch("nGenElectrons",&nGenElectrons);
    outTree->Branch("nGenMuons",&nGenMuons);
    outTree->Branch("nGenTaus",&nGenTaus);
    outTree->Branch("nGenTause",&nGenTause);
    outTree->Branch("nGenTausmu",&nGenTausmu);
    outTree->Branch("nGenTaush",&nGenTaush);
    

    // event counter
    string outHistName1 = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName1.c_str(),"TotalEvents",20,0.5,20.5);
    string outHistName2 = params->get_output_treename("GenEvents");
    hGenEvents = new TH1D(outHistName2.c_str(),"GenEvents",500,0.5,500.5);


    ReportPostBegin();
}

Bool_t MultileptonAnalyzer::Process(Long64_t entry)
{

    ////////////////////////////////////////
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

    eventWeightPu      = 1;
    eventWeightTrigger = 1;
    eventWeightIso     = 1;
    eventWeightId      = 1;
    eventWeight        = 1;
    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
    triggerStatus = passTrigger;
    nPV           = fPVArr->GetEntries();
    if (!isData) {
        nPU = fInfo->nPUmean;
        eventWeightPu *= weights->GetPUWeight(fInfo->nPUmean); // pileup reweighting
    } else {
        nPU = 0;
    }

    ///////////////////////
    // Generator objects //
    ///////////////////////
    vector<TGenParticle*> genParticles;

    if (!isData) {
        unsigned genelectronscount = 0;
        unsigned genmuoncount = 0;
        unsigned gentaucount = 0;
        unsigned gentauecount=0;
        unsigned gentaumucount=0;

        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
            if(particle->parent <0) continue;
            TGenParticle* particleparent = (TGenParticle*) fGenParticleArr->At(particle->parent);

            if ( abs(particle->pdgId) == 11 && abs(particleparent->pdgId) == 24 ) ++genelectronscount;
            if ( abs(particle->pdgId) == 13 && abs(particleparent->pdgId) == 24 ) ++genmuoncount;
            if ( abs(particle->pdgId) == 15 && abs(particleparent->pdgId) == 24 ) ++gentaucount;
            if ( abs(particle->pdgId) == 12 && abs(particleparent->pdgId) == 15 ) {
                TGenParticle* temptau = particleparent;
                if(temptau->parent <0) continue;
                TGenParticle* temptaumother = (TGenParticle*) fGenParticleArr->At(temptau->parent);
                while (abs( temptaumother->pdgId)==15) {
                    temptau = temptaumother;
                    if(temptau->parent >=0) temptaumother = (TGenParticle*) fGenParticleArr->At(temptau->parent);
                    else break;
                    }
                if (abs(temptaumother->pdgId) == 24) ++gentauecount;
                }
            if ( abs(particle->pdgId) == 14 && abs(particleparent->pdgId) == 15 ){
                TGenParticle* temptau = particleparent;
                if(temptau->parent <0) continue;
                TGenParticle* temptaumother = (TGenParticle*) fGenParticleArr->At(temptau->parent);
                while (abs( temptaumother->pdgId)==15) {
                    temptau = temptaumother;
                    if(temptau->parent >=0) temptaumother = (TGenParticle*) fGenParticleArr->At(temptau->parent);
                    else break;
                    }
                if (abs(temptaumother->pdgId) == 24) ++gentaumucount;
                }

            if ( abs(particle->pdgId) == 11 && ((abs(particleparent->pdgId) == 24)|| (abs(particleparent->pdgId) == 15)) ) genParticles.push_back(particle); 
            if ( abs(particle->pdgId) == 13 && ((abs(particleparent->pdgId) == 24)|| (abs(particleparent->pdgId) == 15)) ) genParticles.push_back(particle); 
        }
        nGenElectrons = genelectronscount;
        nGenMuons = genmuoncount;
        nGenTaus = gentaucount;
        nGenTause = gentauecount;
        nGenTausmu= gentaumucount;
        nGenTaush = gentaucount - gentauecount - gentaumucount;
    } 


    unsigned nGenState = 300;// others or data
    if (nGenTaush<0) {nGenState = 0;}
    else if (nGenElectrons==2 && nGenMuons==0 && nGenTaus==0) nGenState = 10; //ee
    else if (nGenElectrons==1 && nGenMuons==1 && nGenTaus==0) nGenState = 20; //emu
    else if (nGenElectrons==0 && nGenMuons==2 && nGenTaus==0) nGenState = 30; //mumu
    else if (nGenElectrons==1 && nGenMuons==0 && nGenTause ==1) nGenState = 40; //etau_e
    else if (nGenElectrons==1 && nGenMuons==0 && nGenTausmu==1) nGenState = 50; //etau_mu
    else if (nGenElectrons==1 && nGenMuons==0 && nGenTaush ==1) nGenState = 60; //etau_h
    else if (nGenElectrons==0 && nGenMuons==1 && nGenTause ==1) nGenState = 70; //mutau_e
    else if (nGenElectrons==0 && nGenMuons==1 && nGenTausmu==1) nGenState = 80; //mutau_mu
    else if (nGenElectrons==0 && nGenMuons==1 && nGenTaush ==1) nGenState = 90; //mutau_h

    else if (nGenElectrons==0 && nGenMuons==0 && nGenTause==2) nGenState = 100; // tau_e + tau_e
    else if (nGenElectrons==0 && nGenMuons==0 && nGenTause==1 && nGenTausmu==1) nGenState = 110; // tau_e + tau_mu
    else if (nGenElectrons==0 && nGenMuons==0 && nGenTausmu==2) nGenState = 120; // tau_mu + tau_mu
    else if (nGenElectrons==0 && nGenMuons==0 && nGenTause==1 && nGenTaush==1) nGenState = 130; // tau_e + tau_h
    else if (nGenElectrons==0 && nGenMuons==0 && nGenTausmu==1 && nGenTaush==1) nGenState = 140; // tau_mu + tau_h
    else if (nGenElectrons==0 && nGenMuons==0 && nGenTaush==2) nGenState = 150; // tau_h + tau_h


    else if (nGenElectrons + nGenMuons + nGenTaus ==1 && nGenElectrons==1) nGenState = 200; // eh
    else if (nGenElectrons + nGenMuons + nGenTaus ==1 && nGenMuons==1)     nGenState = 210; // muh
    else if (nGenElectrons + nGenMuons + nGenTaus ==1 && nGenTaus==1 && nGenTause ==1) nGenState = 220; // tau_e h
    else if (nGenElectrons + nGenMuons + nGenTaus ==1 && nGenTaus==1 && nGenTausmu==1) nGenState = 230; // tau_mu h
    else if (nGenElectrons + nGenMuons + nGenTaus ==1 && nGenTaus==1 && nGenTaush ==1) nGenState = 240; // tau_h h

    else if (nGenElectrons + nGenMuons + nGenTaus ==0) nGenState = 250; // hadronic
    hGenEvents->Fill(nGenState+1);
    hGenEvents->Fill(nGenState+2); // for MC, this does not matter
    hGenEvents->Fill(nGenState+3); // for MC, this does not matter

    //cout<<nGenTaus<<","<<nGenTause<<","<<nGenTausmu<<","<<nGenTaush<<","<<nGenState<<endl;
     


    


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
    hGenEvents->Fill(nGenState+4);

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
                muon->pt > 3
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
    vector<TElectron*> electrons;
    vector<TLorentzVector> veto_electrons; // for vetoing jets that overlap with electrons (not saved!!!)
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        if (
                electron->pt > 3
                && fabs(electron->eta) < 2.5
                && particleSelector->PassElectronID(electron, cuts->tightElID)
                && particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
           ) {
            TLorentzVector electronP4;
            copy_p4(electron, ELE_MASS, electronP4);
            electrons.push_back(electron);
            veto_electrons.push_back(electronP4);
        }
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);

    /* 2.5 TAUS */
    vector<TTau*> taus;
    vector<TLorentzVector> veto_taus;
    for (int i=0; i<fTauArr->GetEntries(); i++){
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
        if (elOverlap ) hTotalEvents->Fill(11);
        if (muOverlap ) hTotalEvents->Fill(12);
        if( // Selection adopted from Stany's tau counting selection for ttDM
                tau->pt > 18  
                && abs(tau->eta) < 2.3 
                && !muOverlap
                && !elOverlap
                && (tau->hpsDisc & baconhep::kByDecayModeFinding)
                && (tau->hpsDisc & baconhep::kByLooseCombinedIsolationDBSumPtCorr3Hits) //kByTightIsolationMVA3oldDMwLT)
                && (tau->hpsDisc & baconhep::kByMVA6VTightElectronRejection)
                && (tau->hpsDisc & baconhep::kByTightMuonRejection3)
          ) {
            taus.push_back(tau);
            veto_taus.push_back(tauP4);
        }
    }
    sort(taus.begin(), taus.end(), sort_by_higher_pt<TTau>);

    

    /* 3.JETS */
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets;
    std::vector<TJet*> bjets;
    nJets    = 0; //light jet
    nBJets   = 0;
    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

        // Prevent overlap of muons and jets
        TLorentzVector vJet; 
        vJet.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
        bool muOverlap = false;
        for (const auto& mu: veto_muons) {
            if (vJet.DeltaR(mu) < 0.4) {
                muOverlap = true;
                break;
            }
        }
        bool elOverlap = false;
        for (const auto& ele: veto_electrons) {
            if (vJet.DeltaR(ele) < 0.4) {
                elOverlap = true;
                break;
            }
        }

        bool tauOverlap = false;
        for (const auto& tau: veto_taus) {
            if (vJet.DeltaR(tau) < 0.4) {
                tauOverlap = true;
                break;
            }
        }
        if (elOverlap ) hTotalEvents->Fill(15);
        if (muOverlap ) hTotalEvents->Fill(16);
        if (tauOverlap) hTotalEvents->Fill(17);

        if (
                jet->pt > 30 
                && fabs(jet->eta) < 2.4
                && particleSelector->PassJetID(jet, cuts->looseJetID)
                && !muOverlap 
                && !elOverlap
                && !tauOverlap
           ) {
            jets.push_back(jet); 
            ++nJets;
            if (isData) {
                if (jet->bmva > 0.9432) { 
                    ++nBJets;
                } 
            } else {
                if (particleSelector->BTagModifier(jet, "MVAT")) { // accounts for b jet efficiency (don't worry about this for now)
                    ++nBJets;
                }
            }
        }
    }
    std::sort(jets.begin(), jets.end(), sort_by_btag);  


    /* 4.MET */
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;

    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    nMuons     = muons.size();
    nElectrons = electrons.size();

    if (params->selection == "mumu") {
        if ((muons.size()!=2) || (nElectrons!=0) )
            return kTRUE;
        hTotalEvents->Fill(5);
        hGenEvents->Fill(nGenState+5);

        if (muons[0]->pt < 25 || muons[1]->pt < 3)
            return kTRUE;
        hTotalEvents->Fill(6);
        hGenEvents->Fill(nGenState+6);

        if (muons[0]->q == muons[1]->q)  // remove same sign muons
            return kTRUE;
        hTotalEvents->Fill(7);
        hGenEvents->Fill(nGenState+7);

        if ( (nJets< 2) || (nBJets<1) ) 
            return kTRUE;
        hTotalEvents->Fill(8);
        hGenEvents->Fill(nGenState+8);

        TLorentzVector muonOneP4, muonTwoP4;
        copy_p4(muons[0], MUON_MASS, muonOneP4);
        copy_p4(muons[1], MUON_MASS, muonTwoP4);
        leptonOneP4 = muonOneP4;
        leptonTwoP4 = muonTwoP4;
        leptonOneQ = muons[0]->q;
        leptonTwoQ = muons[1]->q;
        leptonOneDZ  = muons[0]->dz;
        leptonOneD0  = muons[0]->d0;
        leptonTwoDZ  = muons[1]->dz;
        leptonTwoD0  = muons[1]->d0;
        leptonOneSip3d = muons[0]->sip3d;
        leptonTwoSip3d = muons[1]->sip3d;
    

        // fill b jets
        bjetOneP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);   
        bjetTwoP4.SetPtEtaPhiM(jets[1]->pt, jets[1]->eta, jets[1]->phi, jets[1]->mass);  
        bjetOneBmva = jets[0]->bmva;
        bjetTwoBmva = jets[1]->bmva; 

        if (!isData) {

            leptonOneMother = GetGenMotherId(genParticles, leptonOneP4);
            leptonTwoMother = GetGenMotherId(genParticles, leptonTwoP4);


            eventWeightId  *= weights->GetMuonIDEff(muonOneP4);
            eventWeightIso *= weights->GetMuonISOEff(muonOneP4);
            eventWeightId  *= weights->GetMuonIDEff(muonTwoP4);
            eventWeightIso *= weights->GetMuonISOEff(muonTwoP4);

            // trigger weight
            pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
            pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
            eventWeightTrigger *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);
        }
    } 

    if (params->selection == "mutau") {
        //cout<<taus.size() <<endl;
        if ( muons.size() !=1 || taus.size() !=1 || nElectrons!=0)
            return kTRUE;
        hTotalEvents->Fill(5);
        hGenEvents->Fill(nGenState+5);

        if (muons[0]->pt < 25 || taus[0]->pt < 18)
            return kTRUE;
        hTotalEvents->Fill(6);
        hGenEvents->Fill(nGenState+6);

        if (muons[0]->q == taus[0]->q)
            return kTRUE;
        hTotalEvents->Fill(7);
        hGenEvents->Fill(nGenState+7);

        if ( (nJets< 2) )// || (nBJets<0.1)) 
            return kTRUE;
        hTotalEvents->Fill(8);
        hGenEvents->Fill(nGenState+8);

        TLorentzVector muonOneP4, tauOneP4;
        copy_p4(muons[0], MUON_MASS, muonOneP4);
        copy_p4(taus[0], taus[0]->m, tauOneP4);
        leptonOneP4 = muonOneP4;
        leptonTwoP4 = tauOneP4;

        leptonOneQ = muons[0]->q;
        leptonOneDZ  = muons[0]->dz;
        leptonOneD0  = muons[0]->d0;
        leptonOneSip3d = muons[0]->sip3d;

        leptonTwoQ = taus[0]->q;
        //leptonTwoDZ  = taus[0]->dz;
        //leptonTwoD0  = taus[0]->d0;
        //leptonTwoSip3d = taus[0]->sip3d;


        // fill b jets
        bjetOneP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);   
        bjetTwoP4.SetPtEtaPhiM(jets[1]->pt, jets[1]->eta, jets[1]->phi, jets[1]->mass);  
        bjetOneBmva = jets[0]->bmva;
        bjetTwoBmva = jets[1]->bmva; 

        if (!isData) {

            leptonOneMother = GetGenMotherId(genParticles, leptonOneP4);
            leptonTwoMother = GetGenMotherId(genParticles, leptonTwoP4);

            eventWeightId  *= weights->GetMuonIDEff(muonOneP4);
            eventWeightIso *= weights->GetMuonISOEff(muonOneP4);

            // trigger weight
            pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
            eventWeightTrigger *= trigEff1.first;
        }
    } 

    if (params->selection == "ee") {

        if ( (nMuons!=0) || (nElectrons!=2) )
            return kTRUE;
        hTotalEvents->Fill(5);
        hGenEvents->Fill(nGenState+5);

        if (electrons[0]->pt < 28 || electrons[1]->pt < 3)
            return kTRUE;
        hTotalEvents->Fill(6);
        hGenEvents->Fill(nGenState+6);

        if (electrons[0]->q == electrons[1]->q)  // remove same sign muons
            return kTRUE;
        hTotalEvents->Fill(7);
        hGenEvents->Fill(nGenState+7);

        if ( (nJets< 2) || (nBJets<1) ) 
            return kTRUE;
        hTotalEvents->Fill(8);
        hGenEvents->Fill(nGenState+8);

        TLorentzVector electronOneP4, electronTwoP4;
        copy_p4(electrons[0], ELE_MASS, electronOneP4);
        copy_p4(electrons[1], ELE_MASS, electronTwoP4);
        leptonOneP4 = electronOneP4;
        leptonTwoP4 = electronTwoP4;
        leptonOneQ = electrons[0]->q;
        leptonTwoQ = electrons[1]->q;
        // fill b jets
        bjetOneP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);   
        bjetTwoP4.SetPtEtaPhiM(jets[1]->pt, jets[1]->eta, jets[1]->phi, jets[1]->mass);  
        bjetOneBmva = jets[0]->bmva;
        bjetTwoBmva = jets[1]->bmva; 
        

        if (!isData) {
            leptonOneMother = GetGenMotherId(genParticles, leptonOneP4);
            leptonTwoMother = GetGenMotherId(genParticles, leptonTwoP4);

            eventWeight *= weights->GetElectronRecoIdEff(electronOneP4);
            eventWeight *= weights->GetElectronRecoIdEff(electronTwoP4);

            // trigger weight
            //pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", electronOneP4);
            //pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", electronTwoP4);
            //eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);
        }
    }
    
    if (params->selection == "etau"){
    if ( muons.size() !=0 || taus.size() !=1 || nElectrons!=1 )
            return kTRUE;
    hTotalEvents->Fill(5);
    hGenEvents->Fill(nGenState+5);

    if (electrons[0]->pt < 28 || taus[0]->pt < 18)
        return kTRUE;
    hTotalEvents->Fill(6);
    hGenEvents->Fill(nGenState+6);
    cout<<"start charge"<<endl;
    if (electrons[0]->q == taus[0]->q)
        return kTRUE;
    hTotalEvents->Fill(7);
    hGenEvents->Fill(nGenState+7);

    if ( (nJets< 2)  || (nBJets<1)) 
        return kTRUE;
    hTotalEvents->Fill(8);
    hGenEvents->Fill(nGenState+8);

    TLorentzVector electronOneP4, tauOneP4;
    copy_p4(electrons[0], ELE_MASS, electronOneP4);
    copy_p4(taus[0], taus[0]->m, tauOneP4);
    leptonOneP4 = electronOneP4;
    leptonTwoP4 = tauOneP4;
    leptonOneQ = electrons[0]->q;
    leptonTwoQ = taus[0]->q;

    // fill b jets
    bjetOneP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);   
    bjetTwoP4.SetPtEtaPhiM(jets[1]->pt, jets[1]->eta, jets[1]->phi, jets[1]->mass);  
    bjetOneBmva = jets[0]->bmva;
    bjetTwoBmva = jets[1]->bmva; 


    if (!isData) {
        leptonOneMother = GetGenMotherId(genParticles, leptonOneP4);
        leptonTwoMother = GetGenMotherId(genParticles, leptonTwoP4);

        eventWeight *= weights->GetElectronRecoIdEff(electronOneP4);
        // trigger weight
        //pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", electronOneP4);
        //pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*", electronTwoP4);
        //eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);
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


    //bjets.push_back(jet); //************************
    //jets.push_back(jet);  //************************
    //bjets.push_back(jet); //************************
    //jets.push_back(jet);  //************************
    //std::sort(bjets.begin(), bjets.end(), sort_by_higher_pt<TJet>);  //************************
    //std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);  //************************
    //if (bjets.size() < 2) //************************
    //************************
    //bjetOneP4.SetPtEtaPhiM(bjets[0]->pt, bjets[0]->eta, bjets[0]->phi, bjets[0]->mass); 
    //bjetTwoP4.SetPtEtaPhiM(bjets[1]->pt, bjets[1]->eta, bjets[1]->phi, bjets[1]->mass); 
    //bjetOneBmva = bjets[0]->bmva;
    //bjetTwoBmva = bjets[1]->bmva;
    //************************