#include "BTagSelector.h"
#include <map>

//
// See header file for class documentation
//

using namespace std;
using namespace baconhep;

const float ZMASS = 91.19;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 
bool sort_by_btag(const baconhep::TJet* lhs, const baconhep::TJet* rhs) 
{
    return lhs->bmva > rhs->bmva;
}


BTagSelector::BTagSelector() : BLTSelector()
{

}

BTagSelector::~BTagSelector()
{

}

void BTagSelector::Begin(TTree *tree)
{

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
    //triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");
    triggerNames.push_back("HLT_IsoMu24_v*");
    triggerNames.push_back("HLT_IsoTkMu24_v*");

    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask

    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
    lumiMask.AddJSONFile(jsonFileName);

    // Prepare the output tree
    string outFileName = params->get_output_filename("output");
    outFile = new TFile(outFileName.c_str(),"RECREATE");

    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);

    outFile->cd();
    string treeName = "bltTree_" + params->datasetgroup;
    outTree = new TTree(treeName.c_str(), treeName.c_str());

    // event data
    outTree->Branch("jetP4", &jetP4);
    outTree->Branch("jetTag", &jetTag);
    outTree->Branch("jetFlavor", &jetFlavor);

    ReportPostBegin();
}

Bool_t BTagSelector::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
    outFile->cd();
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

    hTotalEvents->Fill(4);
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);

    /* JETS */
    for (int i=0; i < fAK4CHSArr->GetEntries(); i++) {
        TJet* jet = (TJet*) fAK4CHSArr->At(i);
        assert(jet);

        if (
                jet->pt > 30 
                && fabs(jet->eta) <= 2.4
                && particleSelector->PassJetID(jet, cuts->looseJetID)
           ) {

            TLorentzVector jetP4; 
            jetP4.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
            jetTag    = jet->csv;
            jetFlavor = jet->hadronFlavor;
            ++this->passedEvents;
            outTree->Fill();

        }
    }

    return kTRUE;
}

void BTagSelector::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void BTagSelector::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void BTagSelector::ReportPostTerminate()
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
    std::unique_ptr<BTagSelector> selector(new BTagSelector());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
