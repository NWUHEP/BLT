#include "ThreeMuonAnalyzer.h"
#include <map>

//

// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool sync_print = false;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 
bool sort_by_btag(const baconhep::TJet* lhs, const baconhep::TJet* rhs) {return lhs->bmva > rhs->bmva;}

ThreeMuonAnalyzer::ThreeMuonAnalyzer() : BLTSelector()
{

}

ThreeMuonAnalyzer::~ThreeMuonAnalyzer()
{

}

void ThreeMuonAnalyzer::Begin(TTree *tree)
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


    if (params->selection == "3mu" ) {
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");
    } else if (params->selection == "ee" ) {
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
    outTree->Branch("leptonThreeP4", &leptonThreeP4);

    outTree->Branch("leptonThreeISO", &leptonThreeISO);



    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nBJets", &nBJets);

    

    // event counter
    string outHistName1 = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName1.c_str(),"TotalEvents",20,0.5,20.5);


    ReportPostBegin();
}

Bool_t ThreeMuonAnalyzer::Process(Long64_t entry)
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
    // nothing to do

     


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
        muon->pt = muonSF*muon->pt; 
        muonP4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);

        // Fill containers
        if (
                muon->pt > 10
                && fabs(muon->eta) < 2.4
                // tight muon ID
                && (muon->typeBits & baconhep::kPOGTightMuon)
                // no ISO required
                // && GetMuonIsolation(muon)/muonP4.Pt() < 1.0  
            ) {
            muons.push_back(muon);
            //muons for jet veto
            if (muonP4.Pt() > 20
                // tight ISO
                && GetMuonIsolation(muon)/muonP4.Pt() < 0.15
                ) veto_muons.push_back(muonP4);
        }
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);

    /* 2. ELECTRONS */
    vector<TElectron*> electrons;
    vector<TLorentzVector> veto_electrons; // for vetoing jets that overlap with electrons (not saved!!!)
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);
        //cout<<electron->mva <<endl;

        if (
                electron->pt > 10
                && fabs(electron->eta) < 2.5
                // tight ID and tight ISO
                && particleSelector->PassElectronID(electron, cuts->tightElID) // cut-based ID in barrel+endcup //&& particleSelector->PassElectronMVA(electron, cuts->tightMVAElID)
                && particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl) //PF ISO
            ) {
            electrons.push_back(electron);

            TLorentzVector electronP4;
            copy_p4(electron, ELE_MASS, electronP4);
            if (electronP4.Pt() > 20) veto_electrons.push_back(electronP4);
        }
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);

    

    /* 3.JETS */
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    nJets    = 0; 
    nBJets   = 0;

    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

        if (isData) { // fix for broken bacon JEC
            double jec = particleSelector->JetCorrector(jet, "NONE");
            jet->pt = jet->ptRaw*jec;
        }


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

        if (
                jet->pt > 30 
                && fabs(jet->eta) < 2.4
                && particleSelector->PassJetID(jet, cuts->looseJetID)
                && !muOverlap 
                && !elOverlap
            ) {
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


    /* 4.MET */
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;

    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    nMuons     = muons.size();
    nElectrons = electrons.size();


    if (params->selection == "3mu") {

        
        if (  muons.size() < 3 || muons[0]->pt <25 || muons[1]->pt <20) 
            return kTRUE;
        hTotalEvents->Fill(5);

        if ( GetMuonIsolation(muons[0])/muons[0]->pt > 0.15 || GetMuonIsolation(muons[1])/muons[1]->pt > 0.15  )
            return kTRUE;
        hTotalEvents->Fill(6);          


        TLorentzVector muonOneP4, muonTwoP4, muonThreeP4;

        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, 0.1052);
        muonThreeP4.SetPtEtaPhiM(muons[2]->pt, muons[2]->eta, muons[1]->phi, 0.1052);

        TLorentzVector muonOneTwoP4 = muonOneP4 + muonTwoP4;

        // remove leading and trailing come from Z
        if ( (muonOneTwoP4.M() > 105) || (muonOneTwoP4.M() < 75) || (muons[0]->q == muons[1]->q) )  
            return kTRUE;
        hTotalEvents->Fill(7);

        if (met >30)
             return kTRUE;
        hTotalEvents->Fill(8);

        if (muonOneTwoP4.DeltaPhi(muonThreeP4) < 0.5)
             return kTRUE;
        hTotalEvents->Fill(9);

        
        leptonOneP4 = muonOneP4;
        leptonTwoP4 = muonTwoP4;
        leptonThreeP4 = muonThreeP4;

        leptonThreeISO = GetMuonIsolation(muons[2])/muons[2]->pt ;

    

        if (!isData) {

            eventWeightId  *= weights->GetMuonIDEff(muonOneP4);
            eventWeightIso *= weights->GetMuonISOEff(muonOneP4);
            eventWeightId  *= weights->GetMuonIDEff(muonTwoP4);
            eventWeightIso *= weights->GetMuonISOEff(muonTwoP4);
            
            eventWeightId  *= weights->GetMuonIDEff(muonThreeP4);
            eventWeightIso *= weights->GetMuonISOEff(muonThreeP4);

            // trigger weight
            pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
            pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
            pair<float, float> trigEff3 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonThreeP4);
            eventWeightTrigger *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first)*(1 - trigEff3.first);
        }
    } 
    

    outTree->Fill();
    this->passedEvents++; // controller variables
    return kTRUE;
}

void ThreeMuonAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void ThreeMuonAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void ThreeMuonAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<ThreeMuonAnalyzer> selector(new ThreeMuonAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}



float ThreeMuonAnalyzer::GetMuonIsolation( baconhep::TMuon* mu)
{
    float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    return combIso;
}