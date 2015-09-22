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
    TString option = GetOption();

    params.reset(new AnalysisParameters());
    cuts.reset(new AnalysisCuts());
    triggerSelector.reset(new TriggerSelector());
    particleSelector.reset(new ParticleSelector());

    outFileName = "demoFile_"+params->dataname+"_"+params->selection+"_"+params->jobCount+".root";
    outTreeName = "demoTree_"+params->suffix;

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

    outTree->Fill();
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
    std::unique_ptr<DemoAnalyzer> selector(new DemoAnalyzer());
    int ret = selector->MakeMeSandwich(argc, argv);
    return ret;
}
