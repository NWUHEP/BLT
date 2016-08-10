#include "DimuonAnalyzer.hh"

//
// See header file for class documentation
//


using namespace baconhep;

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
    triggerSelector.reset(new TriggerSelector());

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
    triggerSelector->SetRealData(isRealData);

    /* Trigger selection */

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
    std::vector<TMuon*> muons;
    for (int i=0; i<fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);

        if (
            muon->pt > 20 
            && std::abs(muon->eta) < 2.4
            && particleSelector->PassMuonID(muon, cuts->tightMuID)
            && particleSelector->PassMuonIso(muon, cuts->tightMuIso)
           ) {
            muons.push_back(muon);
        }
    }

    std::sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);

    /* JETS */
    TClonesArray* jetCollection;
    if (params->period == "2012") 
        jetCollection = fAK5Arr;
    else 
        jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets;
    std::vector<TJet*> bjets;
    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);
        if (
                jet->pt > 25
                && particleSelector->PassJetID(jet, cuts->looseJetID)
           ) {

            if (particleSelector->PassJetID(jet, cuts->bJetID))
                bjets.push_back(jet);
            else 
                jets.push_back(jet);
        }
    }
    std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);
    std::sort(bjets.begin(), bjets.end(), sort_by_higher_pt<TJet>);

    /* MET */
    met = fInfo->pfMET;
    met_phi = fInfo->pfMETphi;

    ////////////////////////////
    /* Apply dimuon selection */
    ////////////////////////////
    
    if (muons.size() != 2) 
        return kTRUE;

    if (muons[0]->q == muons[1]->q)
        return kTRUE;

    if (bjets.size() < 1)
        return kTRUE;

    if (jets.size() < 1)
        return kTRUE;

    /* Prepare Lorentz vectors */
    TLorentzVector muon1, muon2, dimuon;
    copy_p4(muons[0], MUON_MASS, muon1);
    copy_p4(muons[1], MUON_MASS, muon2);
    dimuon = muon1 + muon2;

    if (dimuon.M() < 12.)
        return kTRUE;

    //////////
    // Fill //
    //////////
      
    runNumber = fInfo->runNum;
    evtNumber = fInfo->evtNum;
    lumiSection = fInfo->lumiSec;

    muonOneP4 = muon1;
    muonTwoP4 = muon2;
    dimuonP4  = dimuon;

    jetP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
    bjetP4.SetPtEtaPhiM(bjets[0]->pt, bjets[0]->eta, bjets[0]->phi, bjets[0]->mass);

    outTree->Fill();
    this->passedEvents++;

    delete pfMET;

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
