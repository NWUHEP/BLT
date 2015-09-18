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

}

Bool_t DemoAnalyzer::Process(Long64_t entry)
{
   return kTRUE;
}

void DemoAnalyzer::Terminate()
{

}


// _____________________________________________________________________________
// Main function

int main(int argc, char **argv)
{
    using namespace std;

    string option = "";
    for (int i=1; i<argc; ++i)
        option = option + " " + argv[i];

    TChain* fChain = new TChain("Events");

    unique_ptr<DemoAnalyzer> selector(new DemoAnalyzer());
    assert(selector);

    selector->SetOption(option.c_str());
    selector->Init(fChain);
    selector->Begin(fChain);
    selector->Notify();

    Long64_t firstentry = 0, nentries = fChain->GetEntriesFast(),
             entryNumber = 0, localEntry = 0;
    for (Long64_t entry=firstentry; entry<firstentry+nentries; entry++) {
        entryNumber = fChain->GetEntryNumber(entry);
        if (entryNumber < 0) break;
        localEntry = fChain->LoadTree(entryNumber);
        if (localEntry < 0) break;

        selector->Process(localEntry);
    }

    selector->Terminate();

    return EXIT_SUCCESS;
}
