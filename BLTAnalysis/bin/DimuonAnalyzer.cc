#include "DimuonAnalyzer.hh"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 

DimuonAnalyzer::DimuonAnalyzer() : BLTSelector()
{

}

DimuonAnalyzer::~DimuonAnalyzer()
{

}

void DimuonAnalyzer::Begin(TTree *tree)
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
    std::string trigfilename = cmssw_base + std::string("/src/BaconAna/DataFormats/data/HLTFile_v2");
    trigger.reset(new baconhep::TTrigger(trigfilename));

    // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask.reset(new RunLumiRangeMap());
    if (true) { // this will need to be turned off for MC
        string jsonFileName = "src/BLT/BLTAnalyzer/test/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt";
        lumiMask->AddJSONFile(jsonFileName);
    }

    // Prepare the output tree
    outFileName = params->get_output_filename("output");
    outTreeName = params->get_output_treename("data");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "demoTree");

    outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber);
    outTree->Branch("lumiSection", &lumiSection);

    outTree->Branch("muonOneP4", &muonOneP4);
    outTree->Branch("muonOneIso", &muonOneIso);
    outTree->Branch("muonTwoP4", &muonTwoP4);
    outTree->Branch("muonTwoIso", &muonTwoIso);

    outTree->Branch("jetP4", &jetP4);
    outTree->Branch("jet_d0", &jetD0);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("bjetP4", &bjetP4);
    outTree->Branch("bjet_d0", &bjetD0);
    outTree->Branch("nBJets", &nBJets);

    outTree->Branch("met", &met);
    outTree->Branch("met_phi", &met_phi);

    ReportPostBegin();
}

Bool_t DimuonAnalyzer::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;

    //if (entry%1==0)  std::cout << "... Processing event: " << entry << "." << std::endl;
    if (entry%10000==0)  std::cout << "... Processing event: " << entry << " Run: " << fInfo->runNum << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum << "." << std::endl;

    const bool isRealData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isRealData);

    // Apply lumi mask
    if (isRealData) {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask->HasRunLumi(rl)) 
            return kTRUE;
    }

    /* Trigger selection */
    bool passTrigger;
    passTrigger= trigger->pass("HLT_IsoMu24_eta2p1_v*", fInfo->triggerBits);

    if (!passTrigger)
        return kTRUE;

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
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);

    /* MUONS */
    std::vector<TLorentzVector> muons;
    std::vector<float> muons_iso;
    std::vector<TLorentzVector> veto_muons;
    int leadMuonQ = 0;
    for (int i=0; i<fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);

        if (
                muon->pt > 20 
                && fabs(muon->eta) < 2.4
                // tight muon ID 2012
                //&& (muon->typeBits & baconhep::kPFMuon) 
                && (muon->typeBits & baconhep::kGlobal) 
                && muon->muNchi2    < 10.
                && muon->nMatchStn  > 1
                && muon->nPixHits   > 0
                && fabs(muon->d0)   < 0.2
                && fabs(muon->dz)   < 0.5
                && muon->nTkLayers  > 5 
                && muon->nValidHits > 0
                // tight muon detector ISO 2012
                && muon->trkIso03/muon->pt < 0.1
           ) {
            TLorentzVector muonP4;
            copy_p4(muon, MUON_MASS, muonP4);
            veto_muons.push_back(muonP4);

            if (muon->pt > 25 && fabs(muon->eta) < 2.1) {
                if (muons.size() == 0) { // get the highest pt muon
                    muons.push_back(muonP4);
                    muons_iso.push_back(muon->trkIso03);
                    leadMuonQ = muon->q;
                } else if (muons.size() == 1) { // require second muon to be os
                    if (muon->q != leadMuonQ) {
                        muons.push_back(muonP4);
                        muons_iso.push_back(muon->trkIso03);
                    }
                } 
            }
        }
    }
    std::sort(muons.begin(), muons.end(), P4SortCondition);

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
    if (params->period == "2012") 
        jetCollection = fAK5Arr;
    else 
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
        //if (muOverlap || elOverlap) continue;

        if (
                jet->pt > 30 
                && fabs(jet->eta < 4.7)
                && particleSelector->PassJetID(jet, cuts->looseJetID)
           ) {

            if (fabs(jet->eta) <= 2.4) { 
                if (
                        jet->pt > 30 
                        && particleSelector->PassJetPUID(jet, cuts->looseJetID)
                        && !muOverlap 
                        && !elOverlap
                   ) { 
                    if (jet->csv > 0.898) {
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

    TLorentzVector dimuon;
    dimuon = muons[0] + muons[1];
    if (dimuon.M() < 12. || dimuon.M() > 70.)
        return kTRUE;

    if (fwdjets.size() + jets.size() < 1)
        return kTRUE;

    if (bjets.size() < 1)
        return kTRUE;

    //////////
    // Fill //
    //////////

    runNumber = fInfo->runNum;
    evtNumber = fInfo->evtNum;
    lumiSection = fInfo->lumiSec;

    muonOneP4 = muons[0];
    muonOneIso = muons_iso[0];
    muonTwoP4 = muons[1];
    muonTwoIso = muons_iso[1];
    dimuonP4  = dimuon;

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

void DimuonAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void DimuonAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void DimuonAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<DimuonAnalyzer> selector(new DimuonAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
