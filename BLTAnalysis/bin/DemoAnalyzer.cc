#include "DemoAnalyzer.hh"

//
// See header file for class documentation
//

// _____________________________________________________________________________
// Useful functions

namespace baconhep {
    class TMET : public TObject
    {
    public:
        TMET():  pt(0), phi(0) {}
        ~TMET() {}
        float pt, phi;
    };
}

std::ostream& operator<<(std::ostream& os, const baconhep::TMuon* p) {
    return os << "pt: " << p->pt << " eta: " << p->eta << " phi: " << p->phi << " q: " << p->q;
}

std::ostream& operator<<(std::ostream& os, const baconhep::TElectron* p) {
    return os << "pt: " << p->pt << " eta: " << p->eta << " phi: " << p->phi << " q: " << p->q;
}

std::ostream& operator<<(std::ostream& os, const baconhep::TPhoton* p) {
    return os << "pt: " << p->pt << " eta: " << p->eta << " phi: " << p->phi;
}

std::ostream& operator<<(std::ostream& os, const baconhep::TTau* p) {
    return os << "pt: " << p->pt << " eta: " << p->eta << " phi: " << p->phi << " q: " << p->q;
}

std::ostream& operator<<(std::ostream& os, const baconhep::TJet* p) {
    return os << "pt: " << p->pt << " eta: " << p->eta << " phi: " << p->phi << " mass: " << p->mass;
}

std::ostream& operator<<(std::ostream& os, const baconhep::TMET* p) {
    return os << "pt: " << p->pt << " phi: " << p->phi;
}


// _____________________________________________________________________________
// DemoAnalyzer implementations

using namespace baconhep;

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

    GetEntry(entry, 1);  // load all branches

    if (entry%1==0)  std::cout << "... Processing event: " << entry << "." << std::endl;

    ///////////
    // Muons //
    ///////////

    for(int i=0; i<fMuonArr->GetEntries(); i++) {
        const TMuon* thisMuon = (TMuon*) fMuonArr->At(i);
        assert(thisMuon);

        std::cout << "Muon " << i << ": " << thisMuon << std::endl;
    }

    ///////////////
    // Electrons //
    ///////////////

    for(int i=0; i<fElectronArr->GetEntries(); i++) {
        const TElectron* thisElectron = (TElectron*) fElectronArr->At(i);
        assert(thisElectron);

        std::cout << "Electron " << i << ": " << thisElectron << std::endl;
    }

    /////////////
    // Photons //
    /////////////

    for(int i=0; i<fPhotonArr->GetEntries(); i++) {
        const TPhoton* thisPhoton = (TPhoton*) fPhotonArr->At(i);
        assert(thisPhoton);

        std::cout << "Photon " << i << ": " << thisPhoton << std::endl;
    }

    //////////
    // Taus //
    //////////

    for (int i=0; i<fTauArr->GetEntries(); i++) {
        const TTau* thisTau = (TTau*) fTauArr->At(i);
        assert(thisTau);

        std::cout << "Tau " << i << ": " << thisTau << std::endl;
    }

    //////////
    // Jets //
    //////////

    for(int i=0; i<fAK4CHSArr->GetEntries(); i++) {
        const TJet* thisJet = (TJet*) fAK4CHSArr->At(i);

        std::cout << "Jet " << i << ": " << thisJet << std::endl;
    }

    /////////
    // MET //
    /////////

    TMET* pfMET = new TMET();
    pfMET->pt = fInfo->pfMET;
    pfMET->phi = fInfo->pfMETphi;

    TMET* caloMET = new TMET();
    caloMET->pt = fInfo->caloMET;
    caloMET->phi = fInfo->caloMETphi;

    std::cout << "MET " << "(PF) " << ": " << pfMET << std::endl;
    std::cout << "MET " << "(C)  " << ": " << caloMET << std::endl;


    //////////
    // Fill //
    //////////

    outTree->Fill();

    delete pfMET;
    delete caloMET;

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
    std::cout << "  ==== Begin Job ==========" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  =========================" << std::endl;
}

void DemoAnalyzer::ReportPostTerminate()
{
    std::cout << "  ==== Terminate Job ======" << std::endl;
    std::cout << "  output   : " << params->get_output_filename("demoFile") << std::endl;
    std::cout << "  =========================" << std::endl;
}


// _____________________________________________________________________________
// Main function

int main(int argc, char **argv)
{
    std::unique_ptr<DemoAnalyzer> selector(new DemoAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
