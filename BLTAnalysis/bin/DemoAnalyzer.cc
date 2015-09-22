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
    // Parse command line option
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1),
              std::sregex_token_iterator(), std::back_inserter(options));
    assert (options.size() == 7);

    // Set the parameters
    params.reset(new AnalysisParameters());
    params->suffix    = options[0];
    params->abcd      = options[1];
    params->selection = options[2];
    params->period    = options[3];
    params->dataname  = options[4];
    params->jobcount  = options[5];
    params->pileup    = options[6];

    // Set the cuts
    cuts.reset(new AnalysisCuts());
    triggerSelector.reset(new TriggerSelector());
    particleSelector.reset(new ParticleSelector());

    // Prepare the output tree
    outFileName = params->get_output_filename("demoFile");
    outTreeName = params->get_output_treename("demoTree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "demoTree");

    outTree->Branch("muonOne", &muonOne);
    outTree->Branch("muonTwo", &muonTwo);
    outTree->Branch("dimuon", &dimuon);
    outTree->Branch("genMuonOne", &genMuonOne);
    outTree->Branch("genMuonTwo", &genMuonTwo);
    outTree->Branch("genZ", &genZ);

    ReportPostBegin();
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

    ReportPostTerminate();
}

void DemoAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job ====" << std::endl;
    std::cout << *params << std::endl;
}

void DemoAnalyzer::ReportPostTerminate()
{
    std::cout << "  ==== Terminate Job ====" << std::endl;
    std::cout << "  output   : " << params->get_output_filename("demoFile") << std::endl;
}


// _____________________________________________________________________________
// Main function

int main(int argc, char **argv)
{
    std::unique_ptr<DemoAnalyzer> selector(new DemoAnalyzer());
    int ret = selector->MakeMeSandwich(argc, argv);
    return ret;
}
