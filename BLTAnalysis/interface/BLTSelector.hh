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

class BLTSelector
{
public:
    TTree              *fChain;   //!pointer to the analyzed TTree or TChain
    Long64_t            fStatus;
    TString             fOption;

    BLTSelector() : fChain(0), fStatus(0) { }
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

    baconhep::TEventInfo    *fInfo;
    baconhep::TGenEventInfo *fGenEvtInfo;
    TClonesArray            *fGenParticle;
    TClonesArray            *fElectron;
    TClonesArray            *fMuon;
    TClonesArray            *fTau;
    TClonesArray            *fPhoton;
    TClonesArray            *fPV;
    TClonesArray            *fAK4CHS;
    TClonesArray            *fAK8CHS;
    TClonesArray            *fAddAK8CHS;
    TClonesArray            *fCA15CHS;
    TClonesArray            *fAddCA15CHS;
    TClonesArray            *fAK4Puppi;
    TClonesArray            *fCA8Puppi;
    TClonesArray            *fAddCA8Puppi;
    TClonesArray            *fCA15Puppi;
    TClonesArray            *fAddCA15Puppi;

    TBranch                 *b_Info;
    TBranch                 *b_GenEvtInfo;
    TBranch                 *b_GenParticle;
    TBranch                 *b_Electron;
    TBranch                 *b_Muon;
    TBranch                 *b_Tau;
    TBranch                 *b_Photon;
    TBranch                 *b_PV;
    TBranch                 *b_AK4CHS;
    TBranch                 *b_AK8CHS;
    TBranch                 *b_AddAK8CHS;
    TBranch                 *b_CA15CHS;
    TBranch                 *b_AddCA15CHS;
    TBranch                 *b_AK4Puppi;
    TBranch                 *b_CA8Puppi;
    TBranch                 *b_AddCA8Puppi;
    TBranch                 *b_CA15Puppi;
    TBranch                 *b_AddCA15Puppi;
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
    //fChain->SetMakeClass(1);

    fInfo                    = 0;
    fGenEvtInfo              = 0;
    fGenParticle             = 0;
    fElectron                = 0;
    fMuon                    = 0;
    fTau                     = 0;
    fPhoton                  = 0;
    fPV                      = 0;
    fAK4CHS                  = 0;
    fAK8CHS                  = 0;
    fAddAK8CHS               = 0;
    fCA15CHS                 = 0;
    fAddCA15CHS              = 0;
    fAK4Puppi                = 0;
    fCA8Puppi                = 0;
    fAddCA8Puppi             = 0;
    fCA15Puppi               = 0;
    fAddCA15Puppi            = 0;

    fChain->SetBranchAddress("Info", &fInfo, &b_Info);
    fChain->SetBranchAddress("GenEvtInfo", &fGenEvtInfo, &b_GenEvtInfo);
    fChain->SetBranchAddress("GenParticle", &fGenParticle, &b_GenParticle);
    fChain->SetBranchAddress("Electron", &fElectron, &b_Electron);
    fChain->SetBranchAddress("Muon", &fMuon, &b_Muon);
    fChain->SetBranchAddress("Tau", &fTau, &b_Tau);
    fChain->SetBranchAddress("Photon", &fPhoton, &b_Photon);
    fChain->SetBranchAddress("PV", &fPV, &b_PV);
    fChain->SetBranchAddress("AK4CHS", &fAK4CHS, &b_AK4CHS);
    fChain->SetBranchAddress("AK8CHS", &fAK8CHS, &b_AK8CHS);
    fChain->SetBranchAddress("AddAK8CHS", &fAddAK8CHS, &b_AddAK8CHS);
    fChain->SetBranchAddress("CA15CHS", &fCA15CHS, &b_CA15CHS);
    fChain->SetBranchAddress("AddCA15CHS", &fAddCA15CHS, &b_AddCA15CHS);
    fChain->SetBranchAddress("AK4Puppi", &fAK4Puppi, &b_AK4Puppi);
    fChain->SetBranchAddress("CA8Puppi", &fCA8Puppi, &b_CA8Puppi);
    fChain->SetBranchAddress("AddCA8Puppi", &fAddCA8Puppi, &b_AddCA8Puppi);
    fChain->SetBranchAddress("CA15Puppi", &fCA15Puppi, &b_CA15Puppi);
    fChain->SetBranchAddress("AddCA15Puppi", &fAddCA15Puppi, &b_AddCA15Puppi);
}

Bool_t BLTSelector::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

#endif  // BLTSELECTOR_HH
