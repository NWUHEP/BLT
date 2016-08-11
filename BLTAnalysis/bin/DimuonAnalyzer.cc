#include "DimuonAnalyzer.hh"

//
// See header file for class documentation
//

using namespace baconhep;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 

DimuonAnalyzer::DimuonAnalyzer() : BLTSelector()
{

}

DimuonAnalyzer::~DimuonAnalyzer()
{

}

void DimuonAnalyzer::Begin(TTree *tree)
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
    std::string trigfilename = cmssw_base + std::string("/src/BaconAna/DataFormats/data/HLTFile_v2");
    trigger.reset(new baconhep::TTrigger(trigfilename));

    // Prepare the output tree
    outFileName = params->get_output_filename("output");
    outTreeName = params->get_output_treename("data");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "demoTree");

    outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber);
    outTree->Branch("lumiSection", &lumiSection);
    outTree->Branch("muonOneP4", &muonOneP4);
    outTree->Branch("muonTwoP4", &muonTwoP4);
    outTree->Branch("jetP4", &jetP4);
    outTree->Branch("bjetP4", &bjetP4);
    outTree->Branch("met", &met);
    outTree->Branch("met_phi", &met_phi);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nBJets", &nBJets);

    ReportPostBegin();
}

Bool_t DimuonAnalyzer::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;

    //if (entry%1==0)  std::cout << "... Processing event: " << entry << "." << std::endl;
    if (entry%10000==0)  std::cout << "... Processing event: " << entry << " Run: " << fInfo->runNum << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum << "." << std::endl;

    const bool isRealData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isRealData);

    /* Trigger selection */
    bool passTrigger;
    passTrigger= trigger->pass("HLT_IsoMu24_eta2p1_v*", fInfo->triggerBits);

    if (!passTrigger)
        return kTRUE;

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
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);

    /* MUONS */
    std::vector<TLorentzVector> muons;
    std::vector<TLorentzVector> veto_muons;
    std::vector<unsigned> muon_q;
    for (int i=0; i<fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);

        if (
                muon->pt > 25 
                && std::abs(muon->eta) < 2.1
                && particleSelector->PassMuonID(muon, cuts->tightMuID)
                && particleSelector->PassMuonIso(muon, cuts->amumuMuDetIso)
           ) {
            TLorentzVector muonP4;
            copy_p4(muon, MUON_MASS, muonP4);
            muons.push_back(muonP4);
            muon_q.push_back(muon->q);
        } else if (
                muon->pt > 5 
                && std::abs(muon->eta) < 2.4
                && particleSelector->PassMuonID(muon, cuts->tightMuID)
                ) {
            TLorentzVector muonP4;
            copy_p4(muon, MUON_MASS, muonP4);
            veto_muons.push_back(muonP4);
        }

    }

    std::sort(muons.begin(), muons.end(), P4SortCondition);

    /* JETS */
    TClonesArray* jetCollection;
    if (params->period == "2012") 
        jetCollection = fAK5Arr;
    else 
        jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets;
    std::vector<TJet*> fwdjets;
    std::vector<TJet*> bjets;
    nJets  = 0;
    nBJets = 0;
    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

        if (
                jet->pt > 30 
                && particleSelector->PassJetID(jet, cuts->looseJetID)
                && particleSelector->PassJetPUID(jet, cuts->looseJetID)
           ) {

            // Prevent overlap of muons and jets
            TLorentzVector vJet; 
            vJet.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
            bool muOverlap = false;
            for (const auto& mu: muons) {
                if (vJet.DeltaR(mu) < 0.5) {
                    muOverlap = true;
                    break;
                }
            }
            if (muOverlap) continue;

            if (fabs(jet->eta) <= 2.4 && particleSelector->PassJetID(jet, cuts->bJetID)) {
                bjets.push_back(jet);
                ++nBJets;
            } else {
                if (fabs(jet->eta) <= 2.4) {
                    ++nJets;
                    jets.push_back(jet);
                } else {
                    fwdjets.push_back(jet);
                }
            }
        }
    }
    std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);
    std::sort(fwdjets.begin(), fwdjets.end(), sort_by_higher_pt<TJet>);
    std::sort(bjets.begin(), bjets.end(), sort_by_higher_pt<TJet>);

    /* MET */
    met = fInfo->pfMET;
    met_phi = fInfo->pfMETphi;


    ////////////////////////////
    /* Apply dimuon selection */
    ////////////////////////////

    if (muons.size() != 2) 
        return kTRUE;

    if (muon_q[0] == muon_q[1]) 
        return kTRUE;

    TLorentzVector dimuon;
    dimuon = muons[0] + muons[1];
    if (dimuon.M() < 12. || dimuon.M() > 70.)
        return kTRUE;

    //////////
    // Fill //
    //////////

    runNumber = fInfo->runNum;
    evtNumber = fInfo->evtNum;
    lumiSection = fInfo->lumiSec;

    muonOneP4 = muons[0];
    muonTwoP4 = muons[1];
    dimuonP4  = dimuon;

    if (fwdjets.size() > 0)
        jetP4.SetPtEtaPhiM(fwdjets[0]->pt, fwdjets[0]->eta, fwdjets[0]->phi, fwdjets[0]->mass);
    else
        jetP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);

    bjetP4.SetPtEtaPhiM(bjets[0]->pt, bjets[0]->eta, bjets[0]->phi, bjets[0]->mass);


    outTree->Fill();
    this->passedEvents++;

    return kTRUE;
}

void DimuonAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void DimuonAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void DimuonAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<DimuonAnalyzer> selector(new DimuonAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
