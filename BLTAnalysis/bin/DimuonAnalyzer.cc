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
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_v2";
    trigger.reset(new baconhep::TTrigger(trigfilename));

    // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    if (true) { // this will need to be turned off for MC
        string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/test/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt";
        lumiMask.AddJSONFile(jsonFileName);
    }

    // muon momentum corrections
    muonCorr = new rochcor2012();

    // Prepare the output tree
    string outFileName = params->get_output_filename("output");
    string outTreeName = params->get_output_treename("tree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "bltTree");

    outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber, "eventNumber/l");
    outTree->Branch("lumiSection", &lumiSection);
    outTree->Branch("triggerStatus", &triggerStatus);
    outTree->Branch("eventWeight", &eventWeight);

    outTree->Branch("muonOneP4", &muonOneP4);
    outTree->Branch("muonOneIso", &muonOneIso);
    outTree->Branch("muonOneQ", &muonOneQ);
    outTree->Branch("muonTwoP4", &muonTwoP4);
    outTree->Branch("muonTwoIso", &muonTwoIso);
    outTree->Branch("muonTwoQ", &muonTwoQ);

    outTree->Branch("jetP4", &jetP4);
    outTree->Branch("jetD0", &jetD0);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nFwdJets", &nFwdJets);

    outTree->Branch("bjetP4", &bjetP4);
    outTree->Branch("bjetTag", &bjetTag);
    outTree->Branch("bjetD0", &bjetD0);
    outTree->Branch("nBJets", &nBJets);

    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);

    // Event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);
    
    ReportPostBegin();
}

Bool_t DimuonAnalyzer::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;
    hTotalEvents->Fill(1);

    //if (entry%1==0)  std::cout << "... Processing event: " << entry << "." << std::endl;
    if (entry%10000==0)  std::cout << "... Processing event: " << entry << " Run: " << fInfo->runNum << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum << "." << std::endl;

    const bool isRealData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isRealData);

    // Apply lumi mask
    if (isRealData) {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(2);

    /* Trigger selection */
    bool passTrigger;
    passTrigger= trigger->pass("HLT_IsoMu24_eta2p1_v*", fInfo->triggerBits);
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
    /* Apply a preselection so we can make a collection of muons to clean against */
    vector<TMuon*> tmp_muons;
    for (int i=0; i<fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);

        if (
                muon->pt > 20 
                && fabs(muon->eta) < 2.4
                // tight muon ID
                //&& (muon->typeBits & baconhep::kPFMuon) 
                && (muon->typeBits & baconhep::kGlobal) 
                && muon->muNchi2    < 10.
                && muon->nMatchStn  > 1
                && muon->nPixHits   > 0
                && fabs(muon->d0)   < 0.2
                && fabs(muon->dz)   < 0.5
                && muon->nTkLayers  > 5 
                && muon->nValidHits > 0
           ) {
            tmp_muons.push_back(muon);
        }
    }
    sort(tmp_muons.begin(), tmp_muons.end(), sort_by_higher_pt<TMuon>);

    // Second pass
    vector<TLorentzVector> muons;
    vector<TLorentzVector> veto_muons;
    vector<float> muons_iso;
    vector<float> muons_q;
    for (unsigned i = 0; i < tmp_muons.size(); i++) {
        TMuon* muon = tmp_muons[i];

        TLorentzVector muonP4;
        copy_p4(tmp_muons[i], MUON_MASS, muonP4);

        // Remove muon track pt from muon track isolation variable
        for (unsigned j = i+1; j < tmp_muons.size(); j++) {
            TLorentzVector muon_j;
            copy_p4(tmp_muons[j], MUON_MASS, muon_j);

            if (muonP4.DeltaR(muon_j) < 0.3) {
                muon->trkIso03 = max(0., muon->trkIso03 - muon_j.Pt());
                tmp_muons[j]->trkIso03 = max(0., tmp_muons[j]->trkIso03 - muonP4.Pt());
            }
        }

        // Apply rochester muon momentum corrections
        float qter = 1.;
        if (isRealData) {
            muonCorr->momcor_data(muonP4, muon->q, 0, qter);
        } else {
            muonCorr->momcor_mc(muonP4, muon->q, 0, qter);
        }
        if (muon->trkIso03/muonP4.Pt() < 0.1) {
            veto_muons.push_back(muonP4);

            if (muonP4.Pt() > 25 && fabs(muonP4.Eta()) < 2.1) {
                if (muons.size() == 0) {
                    muons.push_back(muonP4);
                    muons_iso.push_back(muon->trkIso03);
                    muons_q.push_back(muon->q);
                } else if (muons.size() == 1) {
                    if (muon->q != muons_q[0]) {
                        muons.push_back(muonP4);
                        muons_iso.push_back(muon->trkIso03);
                        muons_q.push_back(muon->q);
                    }
                }
            }
        }
    }

    /* ELECTRONS */
    /*std::vector<TLorentzVector> electrons;
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
    */

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
    jetCollection = fAK5Arr;

    std::vector<TJet*> jets;
    std::vector<TJet*> fwdjets;
    std::vector<TJet*> bjets;
    nJets    = 0;
    nFwdJets = 0;
    nBJets   = 0;
    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

        // Prevent overlap of muons and jets
        TLorentzVector vJet; 
        vJet.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
        bool muOverlap = false;
        for (const auto& mu: muons) {
            if (vJet.DeltaR(mu) < 0.5) {
                muOverlap = true;
                break;
            }
        }
        bool elOverlap = false;
        //for (const auto& el: electrons) {
        //    if (vJet.DeltaR(el) < 0.5) {
        //        elOverlap = true;
        //        break;
        //    }
        //}

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
                    ++nFwdJets;
                }
            }
        }
    }

    std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);
    std::sort(fwdjets.begin(), fwdjets.end(), sort_by_higher_pt<TJet>);
    std::sort(bjets.begin(), bjets.end(), sort_by_higher_pt<TJet>);

    /* MET */
    met    = fInfo->pfMET;
    metPhi = fInfo->pfMETphi;

    ////////////////////////////
    /* Apply dimuon selection */
    ////////////////////////////

    if (muons.size() < 2)
        return kTRUE;
    hTotalEvents->Fill(5);

    TLorentzVector dimuon;
    dimuon = muons[0] + muons[1];
    if (dimuon.M() < 12. || dimuon.M() > 70.)
        return kTRUE;
    hTotalEvents->Fill(6);

    ////////////////////////////
    // Get event weights here //
    ////////////////////////////

    /*HERE!!!*/

    //////////
    // Fill //
    //////////

    runNumber   = fInfo->runNum;
    evtNumber   = fInfo->evtNum;
    lumiSection = fInfo->lumiSec;
    eventWeight = 1;

    muonOneP4  = muons[0];
    muonOneIso = muons_iso[0];
    muonOneQ   = muons_q[0];
    muonTwoP4  = muons[1];
    muonTwoIso = muons_iso[1];
    muonTwoQ   = muons_q[1];

    if (bjets.size() > 0) {
        bjetP4.SetPtEtaPhiM(bjets[0]->pt, bjets[0]->eta, bjets[0]->phi, bjets[0]->mass);
        bjetD0   = bjets[0]->d0;
        bjetTag  = bjets[0]->csv;
    } else {
        bjetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        bjetD0   = 0.;
        bjetTag  = 0.;
    }

    if (fwdjets.size() > 0) {
        jetP4.SetPtEtaPhiM(fwdjets[0]->pt, fwdjets[0]->eta, fwdjets[0]->phi, fwdjets[0]->mass);
        jetD0   = fwdjets[0]->d0;
    } else if (jets.size() > 0) {
        jetP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
        jetD0   = jets[0]->d0;
    } else {
        jetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        jetD0   = 0.;
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
