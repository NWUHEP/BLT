#include "DemoAnalyzer.hh"

// _____________________________________________________________________________
// See header file for class documentation

DemoAnalyzer::DemoAnalyzer() : BLTSelector()
{

}

DemoAnalyzer::~DemoAnalyzer()
{

}

void DemoAnalyzer::Begin(TTree *tree)
{
    params.reset(new AnalysisParameters());
    cuts.reset(new AnalysisCuts());
    triggerSelector.reset(new TriggerSelector());
    particleSelector.reset(new ParticleSelector());

    std::string outFileName = "demoFile_"+params->dataname+"_"+params->selection+"_"+params->jobCount+".root";
    std::string outTreeName = "demoTree_"+params->suffix;

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "demoTree");

    outTree->Branch("muonOne", &muonOne);
    outTree->Branch("muonTwo", &muonTwo);
    outTree->Branch("dimuon", &dimuon);
    outTree->Branch("genMuonOne", &genMuonOne);
    outTree->Branch("genMuonTwo", &genMuonTwo);
    outTree->Branch("genZ", &genZ);
}

Bool_t DemoAnalyzer::Process(Long64_t entry)
{
    GetEntry(entry, 1);

    return kTRUE;
}

void DemoAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();
}


// _____________________________________________________________________________
// Main function

int main(int argc, char **argv)
{
    using namespace std;

    string arguments = "";
    for (int i=1; i<argc; ++i)
        arguments = arguments + " " + argv[i];

    TChain* fChain = new TChain("Events");

    unique_ptr<DemoAnalyzer> selector(new DemoAnalyzer());
    assert(selector);

    selector->SetOption(arguments.c_str());
    selector->Begin(fChain);
    selector->Init(fChain);   // Init() follows Begin() in TSelector
    selector->Notify();       // Needed explicitly for the first tree; handled
                              // by TChain for the following trees

    Long64_t firstentry = 0, nentries = fChain->GetEntriesFast(),
             entryNumber = 0, localEntry = 0;
    Long64_t totalEvents = 0, passedEvents = 0;
    for (Long64_t entry=firstentry; entry<firstentry+nentries; entry++) {
        entryNumber = fChain->GetEntryNumber(entry);
        if (entryNumber < 0) break;
        localEntry = fChain->LoadTree(entryNumber);
        if (localEntry < 0) break;

        Bool_t pass = selector->Process(localEntry);
        ++totalEvents;
        if (pass)
            ++passedEvents;
    }

    selector->Terminate();

    return EXIT_SUCCESS;
}
