// =============================================================================
// This class is based on the automatically generated codes from calling
//   Events->MakeSelector("BLTSelector")
// on a Bacon ntuple by ROOT version 6.02/05. But this class is not
// inherited from TSelector.
// =============================================================================


#ifndef BLTSELECTOR_HH
#define BLTSELECTOR_HH

// ROOT headers
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TFileCollection.h>

// Header file for the classes stored in the TTree if any.
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TTau.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"

#include "BLTHelper.hh"


class BLTSelector
{
public:
    TTree              *fChain;   //!pointer to the analyzed TTree or TChain
    Long64_t            fStatus;
    TString             fOption;

    long int            fileCount;  // keep track of number of files opened
    long int            unskimmedEventCount;  // keep track of number of events processed before skimming
    long int            totalEvents;
    long int            passedEvents;

    BLTSelector() : fChain(0), fStatus(0), fOption(""),
                    fileCount(0), unskimmedEventCount(0), totalEvents(0), passedEvents(0) { }
    virtual             ~BLTSelector() { }
    virtual void        Init(TTree *tree);
    virtual Bool_t      Notify();
    virtual Long64_t    GetStatus() const { return fStatus; }
    virtual const char *GetOption() const { return fOption; }
    virtual Int_t       GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
    virtual void        SetOption(const char *option) { fOption = option; }
    virtual void        Begin(TTree *tree);
    virtual Bool_t      Process(Long64_t entry);
    virtual void        Terminate();

    Int_t               MakeMeSandwich(int argc, char **argv);  // DO NOT OVERRIDE!

    baconhep::TEventInfo    *fInfo;
    baconhep::TGenEventInfo *fGenEvtInfo;
    TClonesArray            *fGenParticleArr;
    TClonesArray            *fElectronArr;
    TClonesArray            *fMuonArr;
    TClonesArray            *fTauArr;
    TClonesArray            *fPhotonArr;
    TClonesArray            *fPVArr;
    TClonesArray            *fAK4CHSArr;
    TClonesArray            *fAK8CHSArr;
    TClonesArray            *fAddAK8CHSArr;
    TClonesArray            *fCA15CHSArr;
    TClonesArray            *fAddCA15CHSArr;
    TClonesArray            *fAK4PuppiArr;
    TClonesArray            *fCA8PuppiArr;
    TClonesArray            *fAddCA8PuppiArr;
    TClonesArray            *fCA15PuppiArr;
    TClonesArray            *fAddCA15PuppiArr;

    TBranch                 *b_Info;
    TBranch                 *b_GenEvtInfo;
    TBranch                 *b_GenParticleArr;
    TBranch                 *b_ElectronArr;
    TBranch                 *b_MuonArr;
    TBranch                 *b_TauArr;
    TBranch                 *b_PhotonArr;
    TBranch                 *b_PVArr;
    TBranch                 *b_AK4CHSArr;
    TBranch                 *b_AK8CHSArr;
    TBranch                 *b_AddAK8CHSArr;
    TBranch                 *b_CA15CHSArr;
    TBranch                 *b_AddCA15CHSArr;
    TBranch                 *b_AK4PuppiArr;
    TBranch                 *b_CA8PuppiArr;
    TBranch                 *b_AddCA8PuppiArr;
    TBranch                 *b_CA15PuppiArr;
    TBranch                 *b_AddCA15PuppiArr;

    TFile                   *fCurrentFile;
    TH1D                    *hTotalEvents;
};

void BLTSelector::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrentFile = tree->GetCurrentFile();
    //fChain->SetMakeClass(1);

    hTotalEvents = new TH1D("TotalEvents","TotalEvents",1,-10,10);

    fInfo                    = 0;
    fGenEvtInfo              = 0;
    fGenParticleArr          = 0;
    fElectronArr             = 0;
    fMuonArr                 = 0;
    fTauArr                  = 0;
    fPhotonArr               = 0;
    fPVArr                   = 0;
    fAK4CHSArr               = 0;
    fAK8CHSArr               = 0;
    fAddAK8CHSArr            = 0;
    fCA15CHSArr              = 0;
    fAddCA15CHSArr           = 0;
    fAK4PuppiArr             = 0;
    fCA8PuppiArr             = 0;
    fAddCA8PuppiArr          = 0;
    fCA15PuppiArr            = 0;
    fAddCA15PuppiArr         = 0;

    fChain->SetBranchAddress("Info", &fInfo, &b_Info);
    fChain->SetBranchAddress("GenEvtInfo", &fGenEvtInfo, &b_GenEvtInfo);
    fChain->SetBranchAddress("GenParticle", &fGenParticleArr, &b_GenParticleArr);
    fChain->SetBranchAddress("Electron", &fElectronArr, &b_ElectronArr);
    fChain->SetBranchAddress("Muon", &fMuonArr, &b_MuonArr);
    fChain->SetBranchAddress("Tau", &fTauArr, &b_TauArr);
    fChain->SetBranchAddress("Photon", &fPhotonArr, &b_PhotonArr);
    fChain->SetBranchAddress("PV", &fPVArr, &b_PVArr);
    fChain->SetBranchAddress("AK4CHS", &fAK4CHSArr, &b_AK4CHSArr);
    fChain->SetBranchAddress("AK8CHS", &fAK8CHSArr, &b_AK8CHSArr);
    fChain->SetBranchAddress("AddAK8CHS", &fAddAK8CHSArr, &b_AddAK8CHSArr);
    fChain->SetBranchAddress("CA15CHS", &fCA15CHSArr, &b_CA15CHSArr);
    fChain->SetBranchAddress("AddCA15CHS", &fAddCA15CHSArr, &b_AddCA15CHSArr);
    fChain->SetBranchAddress("AK4Puppi", &fAK4PuppiArr, &b_AK4PuppiArr);
    fChain->SetBranchAddress("CA8Puppi", &fCA8PuppiArr, &b_CA8PuppiArr);
    fChain->SetBranchAddress("AddCA8Puppi", &fAddCA8PuppiArr, &b_AddCA8PuppiArr);
    fChain->SetBranchAddress("CA15Puppi", &fCA15PuppiArr, &b_CA15PuppiArr);
    fChain->SetBranchAddress("AddCA15Puppi", &fAddCA15PuppiArr, &b_AddCA15PuppiArr);
}

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
    TH1* h1 = (TH1*) fCurrentFile->Get("TotalEvents");
    hTotalEvents->Add(h1);
    unskimmedEventCount += h1->GetBinContent(1);

    return kTRUE;
}

#endif  // BLTSELECTOR_HH
