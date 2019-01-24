#include "BaconSkimmer.hh"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 

BaconSkimmer::BaconSkimmer() : BLTSelector()
{

}

BaconSkimmer::~BaconSkimmer()
{

}

void BaconSkimmer::Begin(TTree *tree)
{
    // Parse command line option
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1),
            std::sregex_token_iterator(), std::back_inserter(options));

    // Set the parameters
    params.reset(new Parameters());
    params->setup(options);

    // Set the cuts
    cuts.reset(new Cuts());
    particleSelector.reset(new ParticleSelector(*params, *cuts));

    // Trigger bits mapping file
    const std::string cmssw_base = getenv("CMSSW_BASE");
    std::string trigfilename = cmssw_base + std::string("/src/BaconAna/DataFormats/data/HLTFile_25ns");
    trigger.reset(new baconhep::TTrigger(trigfilename));

    // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary

    lumiMask =  RunLumiRangeMap();
    if (true) { // this will need to be turned off for MC
        std::string jsonFileName = cmssw_base + std::string("/src/BLT/BLTAnalyzer/test/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt");
        lumiMask.AddJSONFile(jsonFileName);
    }

    // Prepare the output tree
    outFileName = params->get_output_filename("output");
    outTreeName = params->get_output_treename("data");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "dataTree");

    outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber);
    outTree->Branch("lumiSection", &lumiSection);
    outTree->Branch("triggerStatus", &triggerStatus);

    outTree->Branch("leptons", &leptons);
    outTree->Branch("jets", &jets);
    outTree->Branch("bjets", &bjets);

    outTree->Branch("met", &met);
    outTree->Branch("met_phi", &met_phi);

    ReportPostBegin();
}

Bool_t BaconSkimmer::Process(Long64_t entry)
{
    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;

    if (entry%10000==0)  std::cout << "... Processing event: " << entry << " Run: " << fInfo->runNum << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum << "." << std::endl;

    const bool isRealData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isRealData);

    /* Lumi mask */
    //if (isRealData) {
    //    RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
    //    if(!lumiMask.HasRunLumi(rl)) 
    //        return kTRUE;
    //}
    hTotalEvents->Fill(2);

    /* Trigger selection */
    bool passTrigger;
    passTrigger = trigger->pass("HLT_IsoMu22_v*", fInfo->triggerBits);
    triggerStatus = passTrigger;
    if (!passTrigger)
        return kTRUE;

    hTotalEvents->Fill(3);

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

    /* MUONS */
    std::vector<TLorentzVector> muons;
    std::vector<TLorentzVector> veto_muons;
    std::vector<float> muons_iso;
    std::vector<float> muons_q;
    bool triggerMuon = false;
    for (int i=0; i<fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);

        if (
                muon->pt > 5 
                && fabs(muon->et
                // tight muon ID
                && (muon->typeBits & baconhep::kPFMuon) 
                && (muon->typeBits & baconhep::kGlobal) 
                && muon->muNchi2    < 10.
                && muon->nMatchStn  > 1
                && muon->nPixHits   > 0
                && fabs(muon->d0)   < 0.2
                && fabs(muon->dz)   < 0.5
                && muon->nTkLayers  > 5 
                && muon->nValidHits > 0
                // loose muon detector ISO
                && muon->trkIso/muon->pt < 0.1
           ) {
            TLorentzVector muonP4;
            copy_p4(muon, MUON_MASS, muonP4);
            muons.push_back(muonP4);
            muons_iso.push_back(muon->trkIso);
            muons_q.push_back(muon->q);

            if (muon->pt > 25 && fabs(muon->eta) < 2.1) {
                triggerMuon = true;
            }
        }
    }
    std::sort(muons.begin(), muons.end(), P4SortCondition);

    if (!triggerMuon)
        return kTRUE;

    /* ELECTRONS */
    std::vector<TLorentzVector> electrons;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        if (
                electron->pt > 20 
                && fabs(electron->eta) < 2.5
                && particleSelector->PassElectronID(electron, cuts->looseElID)
                && particleSelector->PassElectronIso(electron, cuts->looseElIso)
           ) {
            TLorentzVector electronP4;
            copy_p4(electron, ELE_MASS, electronP4);
            electrons.push_back(electronP4);
        }
    }

    std::sort(electrons.begin(), electrons.end(), P4SortCondition);

    /* PHOTONS */
    /* Don't need these just now
    std::vector<TLorentzVector> photons;
    for (int i=0; i<fPhotonArr->GetEntries(); i++) {
        TPhoton* photon = (TPhoton*) fPhotonArr->At(i);
        assert(photon);

        if (
                photon->pt > 10 
                && fabs(photon->eta) < 2.5
                && particleSelector->PassPhotonID(photon, cuts->loosePhID)
                && particleSelector->PassPhotonIso(photon, cuts->loosePhIso, cuts->EAPho)
           ) {
            TLorentzVector photonP4;
            copy_p4(photon, 0., photonP4);
            photons.push_back(photonP4);
        }
    }

    std::sort(photons.begin(), photons.end(), P4SortCondition);
    */

    /* JETS */
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets;
    std::vector<TJet*> fwdjets;
    std::vector<TJet*> bjets;
    nJets  = 0;
    nBJets = 0;
    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

        // Prevent overlap of muons and jets
        TLorentzVector vJet; 
        vJet.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
        bool muOverlap = false;
        for (const auto& mu: veto_muons) {
            if (vJet.DeltaR(mu) < 0.5) {
                muOverlap = true;
                break;
            }
        }
        bool elOverlap = false;
        for (const auto& el: electrons) {
            if (vJet.DeltaR(el) < 0.5) {
                elOverlap = true;
                break;
            }
        }
        if (muOverlap || elOverlap) continue;

        if (
                jet->pt > 30 
                && fabs(jet->eta < 4.7)
                && particleSelector->PassJetID(jet, cuts->looseJetID)
           ) {

            if (fabs(jet->eta) <= 2.4) { 
                if (
                        jet->pt > 30 
                        //&& particleSelector->PassJetPUID(jet, cuts->looseJetID)
                        && !muOverlap 
                        && !elOverlap
                   ) { 
                    if (jet->csv > 0.935) {
                        bjets.push_back(jet);
                        ++nBJets;
                    } else {
                        jets.push_back(jet);
                        ++nJets;
                    }
                }
            } else {
                if (jet->pt > 30) {
                    fwdjets.push_back(jet);
                }
            }
        }
    }

    std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);
    std::sort(fwdjets.begin(), fwdjets.end(), sort_by_higher_pt<TJet>);
    std::sort(bjets.begin(), bjets.end(), sort_by_higher_pt<TJet>);

    /* MET */
    met     = fInfo->pfMET;
    met_phi = fInfo->pfMETphi;

    ////////////////////////////
    /* Apply dimuon selection */
    ////////////////////////////

    if (muons.size() < 2)
        return kTRUE;
    hTotalEvents->Fill(5);

    if (muons_q[0] == muons_q[1])
        return kTRUE;

    TLorentzVector dimuon;
    dimuon = muons[0] + muons[1];
    if (dimuon.M() > 70.)
        return kTRUE;
    hTotalEvents->Fill(6);

    if (fwdjets.size() + jets.size() < 1)
        return kTRUE;
    hTotalEvents->Fill(7);

    if (bjets.size() < 1)
        return kTRUE;
    hTotalEvents->Fill(8);

    //////////
    // Fill //
    //////////

    runNumber = fInfo->runNum;
    evtNumber = fInfo->evtNum;
    lumiSection = fInfo->lumiSec;

    muonOneP4  = muons[0];
    muonOneIso = muons_iso[0];
    muonTwoP4  = muons[1];
    muonTwoIso = muons_iso[1];
    dimuonP4   = dimuon;

    bjetP4.SetPtEtaPhiM(bjets[0]->pt, bjets[0]->eta, bjets[0]->phi, bjets[0]->mass);
    bjetD0   = bjets[0]->d0;
    if (fwdjets.size() > 0) {
        jetP4.SetPtEtaPhiM(fwdjets[0]->pt, fwdjets[0]->eta, fwdjets[0]->phi, fwdjets[0]->mass);
        jetD0   = fwdjets[0]->d0;
    } else {
        jetP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
        jetD0   = jets[0]->d0;
    }

    outTree->Fill();
    this->passedEvents++;
    return kTRUE;
}

void BaconSkimmer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void BaconSkimmer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void BaconSkimmer::ReportPostTerminate()
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
    std::unique_ptr<BaconSkimmer> selector(new BaconSkimmer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
