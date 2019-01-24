#include "TauSelector.h"
#include <map>

//
// See header file for class documentation
//


using namespace std;
using namespace baconhep;



bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 

bool sort_by_btag(const baconhep::TJet* lhs, const baconhep::TJet* rhs) 
{
    return lhs->bmva > rhs->bmva;
}


TauSelector::TauSelector() : BLTSelector()
{

}

TauSelector::~TauSelector()
{

}

void TauSelector::Begin(TTree *tree)
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

    vector<std::string> channelNames = {"tauID"};

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

        tree->Branch("eventWeight", &eventWeight);

        // taus
        // tree->Branch("tauOneP4",         &tauOneP4);
        tree->Branch("tauOnePt",         &tauOnePt);
        tree->Branch("tauOneEta",        &tauOneEta);
        tree->Branch("tauOneFlavor",     &tauOneFlavor);
        tree->Branch("tauOneIso",        &tauOneIso);
        tree->Branch("tauOneDecayMode",  &tauOneDecayMode);
        tree->Branch("tauOneIsoMVA",     &tauOneIsoMVA);
        tree->Branch("tauOnePuppiChHadIso",     &tauOnePuppiChHadIso);
        tree->Branch("tauOnePuppiGammaIso",     &tauOnePuppiGammaIso);
        tree->Branch("tauOnePuppiNeuHadIso",    &tauOnePuppiNeuHadIso);

        // mistaus
        // tree->Branch("tauTwoP4",         &tauTwoP4);
        tree->Branch("tauTwoPt",         &tauTwoPt);
        tree->Branch("tauTwoEta",        &tauTwoEta);
        tree->Branch("tauTwoFlavor",     &tauTwoFlavor);
        tree->Branch("tauTwoIso",        &tauTwoIso);
        tree->Branch("tauTwoDecayMode",  &tauTwoDecayMode);
        tree->Branch("tauTwoIsoMVA",     &tauTwoIsoMVA);
        tree->Branch("tauTwoPuppiChHadIso",     &tauTwoPuppiChHadIso);
        tree->Branch("tauTwoPuppiGammaIso",     &tauTwoPuppiGammaIso);
        tree->Branch("tauTwoPuppiNeuHadIso",    &tauTwoPuppiNeuHadIso);

        // object counters
        tree->Branch("nJets", &nJets);
        tree->Branch("nBJets", &nBJets);
        tree->Branch("nMuons", &nMuons);
        tree->Branch("nElectrons", &nElectrons);
        tree->Branch("nTaus", &nTaus);
        tree->Branch("nMistaus", &nMistaus);

        outTrees[channel] = tree;

        // event counter
        string outHistName = params->get_output_treename("TotalEvents_" + channel);
        eventCounts[channel] = new TH1D(outHistName.c_str(),"ChannelCounts",10,0.5,10.5);
    }

    ReportPostBegin();
}

Bool_t TauSelector::Process(Long64_t entry)
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

    if (!isData) {
        // tauh+tauh or tauh+(e,mu,h)
        // check no leptonic tau

        bool hasPromptTau = false;
        bool hasLeptonicTau = false;
        
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            if (particle->parent > 0) {
                TGenParticle* mother   = (TGenParticle*) fGenParticleArr->At(particle->parent);

                if (    abs(particle->pdgId) == 16
                        && abs(mother->pdgId) == 24
                ) {
                    hasPromptTau = true;
                }
                
                if (    (abs(particle->pdgId) == 12 || abs(particle->pdgId) == 14)
                        && abs(mother->pdgId) == 15
                ) {
                    hasLeptonicTau = true;
                }
            }
        }

        if (hasPromptTau && (!hasLeptonicTau)){

            for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
                TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

                if (particle->parent > 0) {
                    TGenParticle* mother   = (TGenParticle*) fGenParticleArr->At(particle->parent);

                    if (    abs(particle->pdgId) == 16 
                            && abs(mother->pdgId) == 15
                    ) {
                        genParticles.push_back(mother);
                    }
                }
            }

        } else {
            return kTRUE;
        }
    }
    hTotalEvents->Fill(2);

    /* Apply lumi mask */
    if (isData) {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(3);

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

    // if (!passTrigger)
    //     return kTRUE;
    hTotalEvents->Fill(4);


    /////////////////////
    // Fill event info //
    /////////////////////

    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
   

    ///////////////////
    // Select objects//
    ///////////////////

    /* Vertices */
    if (fInfo->hasGoodPV) {
        assert(fPVArr->GetEntries() != 0);
        TVector3 pv;
        copy_xyz((TVertex*) fPVArr->At(0), pv);
        particleSelector->SetPV(pv);
    } else {
        return kTRUE;
    }
    hTotalEvents->Fill(5);
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);

    /* MUONS */
    vector<TMuon*> muons;
    vector<TLorentzVector> veto_muons;
    nMuons = 0;
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
            muonP4.Pt() > 20.
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

            && GetMuonIsolation(muon)/muonP4.Pt() < 0.15

            ) {
                
            muons.push_back(muon);
            veto_muons.push_back(muonP4);
            nMuons ++;
        }
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);

    

    /* ELECTRONS */
    vector<TElectron*> electrons;
    vector<TLorentzVector> veto_electrons;
    nElectrons = 0;
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
            electronP4.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, 511e-6);

            electrons.push_back(electron);
            veto_electrons.push_back(electronP4);
            nElectrons++;
        }

    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);



    /* TAUS */
    vector<TTau*> taus, mistaus;
    nTaus = 0;
    nMistaus = 0;
    for (int i=0; i < fTauArr->GetEntries(); i++) {
        TTau *tau = (TTau*) fTauArr->At(i);
        assert(tau);


        TLorentzVector tauP4; 
        tauP4.SetPtEtaPhiM(tau->pt, tau->eta, tau->phi, 1.77682);

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

            // match tauP4 with gen level hadronic taus
            bool genOverlap = false;
            for (unsigned i = 0; i < genParticles.size(); ++i) {
                TLorentzVector genP4;
                genP4.SetPtEtaPhiM(genParticles[i]->pt, genParticles[i]->eta, genParticles[i]->phi, genParticles[i]->mass); 
                
                if (genP4.DeltaR(tauP4) < 0.1) {
                    genOverlap = true;

                }
            }

            if (genOverlap){
                taus.push_back(tau);
                nTaus ++;
            } else {
                mistaus.push_back(tau);
                nMistaus ++;
            }

        }
    }
    sort(taus.begin(), taus.end(), sort_by_higher_pt<TTau>);

    /* JETS */
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;
    nJets = 0;
    nBJets= 0;
    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

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
            && fabs(jet->eta) <= 2.4
            && particleSelector->PassJetID(jet, cuts->looseJetID)
            && !muOverlap 
            && !elOverlap
            
            ) {

            nJets ++;
            if (particleSelector->BTagModifier(jet, "MVAT", 0, 0, rng->Uniform(1.))) { 
                ++nBJets;
            }

        }
    }

    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    string channel = "";


    if (nMistaus>=1 && nTaus==1  ) { // tauh+jets selection && (nElectrons==1 || nMuons==1)
        channel = "tauID";
        eventCounts[channel]->Fill(1);


        TLorentzVector tau1P4;
        tau1P4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, 1.77682);

        int tau1ISO = 0;
        if (taus[0]->hpsDisc & baconhep::kByLooseIsolationMVA3newDMwLT)  tau1ISO+=1;
        if (taus[0]->hpsDisc & baconhep::kByMediumIsolationMVA3newDMwLT) tau1ISO+=1;
        if (taus[0]->hpsDisc & baconhep::kByTightIsolationMVA3newDMwLT)  tau1ISO+=1;
        if (taus[0]->hpsDisc & baconhep::kByVTightIsolationMVA3newDMwLT) tau1ISO+=1;

        tauOnePt     = taus[0]->pt;
        tauOneEta    = taus[0]->eta;
        tauOneIso    = tau1ISO;
        tauOneFlavor = taus[0]->q*15;
        tauOneDecayMode      = taus[0]->decaymode;
        tauOneIsoMVA         = taus[0]->rawIsoMVA3newDMwLT;
        tauOnePuppiChHadIso  = taus[0]->puppiChHadIso;
        tauOnePuppiGammaIso  = taus[0]->puppiGammaIso;
        tauOnePuppiNeuHadIso = taus[0]->puppiNeuHadIso;


        TLorentzVector tau2P4;
        tau2P4.SetPtEtaPhiM(mistaus[0]->pt, mistaus[0]->eta, mistaus[0]->phi, 1.77682);

        int tau2ISO = 0;
        if (mistaus[0]->hpsDisc & baconhep::kByLooseIsolationMVA3newDMwLT)  tau2ISO+=1;
        if (mistaus[0]->hpsDisc & baconhep::kByMediumIsolationMVA3newDMwLT) tau2ISO+=1;
        if (mistaus[0]->hpsDisc & baconhep::kByTightIsolationMVA3newDMwLT)  tau2ISO+=1;
        if (mistaus[0]->hpsDisc & baconhep::kByVTightIsolationMVA3newDMwLT) tau2ISO+=1;

        // tauTwoP4     = tau2P4;
        tauTwoPt     = mistaus[0]->pt;
        tauTwoEta    = mistaus[0]->eta;
        tauTwoIso    = tau2ISO;
        tauTwoFlavor = mistaus[0]->q*15;
        tauTwoDecayMode      = mistaus[0]->decaymode;
        tauTwoIsoMVA         = mistaus[0]->rawIsoMVA3newDMwLT;
        tauTwoPuppiChHadIso  = mistaus[0]->puppiChHadIso;
        tauTwoPuppiGammaIso  = mistaus[0]->puppiGammaIso;
        tauTwoPuppiNeuHadIso = mistaus[0]->puppiNeuHadIso;

    } else {
        return kTRUE;
    }

    outFile->cd(channel.c_str());
    outTrees[channel]->Fill();
    this->passedEvents++;
    return kTRUE;

}

void TauSelector::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void TauSelector::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void TauSelector::ReportPostTerminate()
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
    std::unique_ptr<TauSelector> selector(new TauSelector());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}




float TauSelector::GetMuonIsolation(const baconhep::TMuon* mu)
{
    float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    return combIso;
}

float TauSelector::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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
