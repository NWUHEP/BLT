#include "BLT/BLTAnalysis/interface/BLTSelector.hh"

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.

void BLTSelector::Begin(TTree *tree)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    //TString option = GetOption();

}

Bool_t BLTSelector::Process(Long64_t entry)
{
    // The Process() function is called for each entry in the tree (or possibly
    // keyed object in the case of PROOF) to be processed. The entry argument
    // specifies which entry in the currently loaded tree is to be processed.
    // It can be passed to either BLTSelector::GetEntry() or TBranch::GetEntry()
    // to read either all or the required parts of the data. When processing
    // keyed objects with PROOF, the object is already loaded and is available
    // via the fObject pointer.
    //
    // This function should contain the "body" of the analysis. It can contain
    // simple or elaborate selection criteria, run algorithms on the data
    // of the event and typically fill histograms.
    //
    // The processing can be stopped by calling Abort().
    //
    // Use fStatus to set the return value of TTree::Process().
    //
    // The return value is currently not used.


    return kTRUE;
}

void BLTSelector::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

}

// _____________________________________________________________________________
// MakeMeSandwich() acts like the "main" function

#include <TStopwatch.h>

#include <iostream>
#include <string>

Int_t BLTSelector::MakeMeSandwich(int argc, char **argv) {
    int verbose = 1;

    int showTiming = 1;
    TStopwatch timer;

    // _________________________________________________________________________
    // Get arguments

    std::string option = "";
    std::string inFileName = "";
    for (int i=1; i<argc; ++i) {
        std::string argi = argv[i];
        if (argi == "-v") {
            verbose = 1;
        } else if (argi == "-vv") {
            verbose = 2;
        } else if (argi == "-vvv") {
            verbose = 3;
        } else {
            if (inFileName == "") {
                inFileName = argi;
            } else if (option == "") {
                option = argi;
            } else {
                option = option + " " + argi;
            }
        }
    }

    // _________________________________________________________________________
    // Handle input files

    TChain* fChain = new TChain("Events");
    std::string inFileExtension = get_file_extension(inFileName);

    if (inFileExtension == "root") {
        if (fChain->Add(inFileName.c_str())) {
            if (verbose >= 1)  std::cout << Info() << "Successfully opened " << inFileName << "." << std::endl;
        } else {
            std::cout << Error() << "Failed to open" << inFileName << std::endl;
            return EXIT_FAILURE;
        }
    } else if (inFileExtension == "txt") {
        TFileCollection fc("fileinfolist", "", inFileName.c_str());
        if (fChain->AddFileInfoList((TCollection*) fc.GetList()) ) {
            if (verbose >= 1)  std::cout << Info() << "Successfully opened " << inFileName << "." << std::endl;
        } else {
            std::cout << Error() << "Failed to open" << inFileName << std::endl;
            return EXIT_FAILURE;
        }
    } else {
        std::cout << Error() << "Unrecognized file extension: " << inFileExtension << std::endl;
        return EXIT_FAILURE;
    }

    if (fChain->LoadTree(0) < 0) {  // needed because TChain =/= TTree
        std::cout << Error() << "Failed to load the entries in the file." << std::endl;
        return EXIT_FAILURE;
    }

    // _________________________________________________________________________
    // Set up the selector

    if (showTiming)
        timer.Start();

    this->SetOption(option.c_str());
    this->Begin(fChain);
    this->Init(fChain); // Init() follows Begin() in TSelector
    if (verbose >= 2)  std::cout << Debug() << "DemoAnalyzer is initialized." << std::endl;

    this->Notify();     // Notify() needs to be called explicitly for the first tree;
                        // handled by TChain for the following trees
    if (verbose >= 2)  std::cout << Debug() << "DemoAnalyzer has opened the first input file." << std::endl;

    // _________________________________________________________________________
    // Loop over events

    Long64_t firstentry = 0, nentries = fChain->GetEntriesFast(),
             entryNumber = 0, localEntry = 0;
    Long64_t totalEvents = 0, passedEvents = 0;

    for (Long64_t entry=firstentry; entry<firstentry+nentries; ++entry) {
        entryNumber = fChain->GetEntryNumber(entry);
        if (entryNumber < 0) break;
        localEntry = fChain->LoadTree(entryNumber);
        if (localEntry < 0) break;

        // The real work is being done here
        Bool_t pass = this->Process(localEntry);

        ++totalEvents;
        if (pass)
            ++passedEvents;
    }
    if (verbose >= 1)  std::cout << Info() << "Selected " << passedEvents << "/" << totalEvents << " events." << std::endl;

    // _________________________________________________________________________
    // Write the output

    this->Terminate();
    if (verbose >= 2)  std::cout << Debug() << "DemoAnalyzer has finished processing." << std::endl;

    // Done
    if (showTiming) {
        timer.Stop();

        std::cout << "  ==== Timing ====" << std::endl;
        std::cout << "  CPU  time: " << timer.CpuTime() << std::endl;
        std::cout << "  Real time: " << timer.RealTime() << std::endl;
    }

    return EXIT_SUCCESS;
}
