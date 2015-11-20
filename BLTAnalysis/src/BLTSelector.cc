#define BLTSelector_cxx
// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//

#include "BLT/BLTAnalysis/interface/BLTSelector.hh"

#include <TStopwatch.h>

#include <iostream>
#include <string>
#include <stdexcept>


void BLTSelector::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

}

void BLTSelector::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

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

void BLTSelector::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.

}

void BLTSelector::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

}

// _____________________________________________________________________________
// Notify() keeps track of the files that are processed

Bool_t BLTSelector::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    // Keep traack of number of files and number of events
    fileCount += 1;
    fCurrentFile = fChain->GetCurrentFile();
    if (!fCurrentFile) {
        throw std::runtime_error("Cannot get current file");
    }
    TH1* h1 = (TH1*) fCurrentFile->Get("TotalEvents");
    if (!h1) {
        throw std::runtime_error("Cannot get 'TotalEvents' histogram");
    }
    hTotalEvents->Add(h1);

    long int eventCount = h1->GetBinContent(1);
    unskimmedEventCount += eventCount;

    std::cout << info() << "Processing " << fCurrentFile->GetName() << " with " << eventCount << " unskimmed events." << std::endl;

    return kTRUE;
}

// _____________________________________________________________________________
// MakeMeSandwich() acts like the "main" function

Int_t BLTSelector::MakeMeSandwich(int argc, char **argv) {

    int verbose = 1;

    int showTiming = 1;
    TStopwatch timer;

    gROOT->SetBatch();


    // _________________________________________________________________________
    // Get arguments

    std::string option = "";
    std::string inputFiles = "";
    std::string sMaxEvents = "";
    for (int i=1; i<argc; ++i) {
        std::string argi = argv[i];
        if (argi == "-v") {
            verbose = 1;
        } else if (argi == "-vv") {
            verbose = 2;
        } else if (argi == "-vvv") {
            verbose = 3;
        } else {
            if (inputFiles == "") {
                inputFiles = argi;
            } else if (sMaxEvents == "") {
                sMaxEvents = argi;
            } else if (option == "") {
                option = argi;
            } else {
                option = option + " " + argi;
            }
        }
    }

    long int maxEvents = std::stol(sMaxEvents);

    // _________________________________________________________________________
    // Handle input files

    TChain* myChain = new TChain("Events");
    std::string inputFileExt = get_file_extension(inputFiles);

    if (inputFileExt == "root") {
        if (myChain->AddFile(inputFiles.c_str())) {
            if (verbose >= 1)  std::cout << info() << "Trying to open " << inputFiles << "." << std::endl;
        } else {
            std::cout << error() << "Failed to open " << inputFiles << std::endl;
            throw std::runtime_error("bad input");
        }
    } else if (inputFileExt == "txt") {
        TFileCollection fc;
        if (fc.AddFromFile(inputFiles.c_str()) && myChain->AddFileInfoList((TCollection*) fc.GetList()) ) {
            if (verbose >= 1)  std::cout << info() << "Trying to open " << inputFiles << "." << std::endl;
        } else {
            std::cout << error() << "Failed to open " << inputFiles << std::endl;
            throw std::runtime_error("bad input");
        }
    } else {
        std::cout << error() << "Unrecognized file extension: " << inputFileExt << std::endl;
        throw std::runtime_error("bad input");
    }

    if (myChain->LoadTree(0) < 0) {  // needed because TChain =/= TTree
        std::cout << error() << "Failed to load the entries in the file." << std::endl;
        throw std::runtime_error("bad input");
    }


    // _________________________________________________________________________
    // Run

    if (showTiming) {
        timer.Start();
    }

    if (maxEvents < 0)
        maxEvents = myChain->GetEntriesFast();
    myChain->Process(this, option.c_str(), maxEvents, 0);

    // Done
    if (showTiming) {
        timer.Stop();
        std::cout << info() << "CPU  time: " << timer.CpuTime() << std::endl;
        std::cout << info() << "Real time: " << timer.RealTime() << std::endl;
    }

    return 0;
}
