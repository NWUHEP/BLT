#include "MuonTagAndProbe.h"
#include <map>

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 

MuonTagAndProbe::MuonTagAndProbe() : BLTSelector()
{

}

MuonTagAndProbe::~MuonTagAndProbe()
{

}

void MuonTagAndProbe::Begin(TTree *tree)
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
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns";
    trigger.reset(new baconhep::TTrigger(trigfilename));

    if (params->selection == "mumu" || params->selection == "emu") {
        //triggerNames.push_back("HLT_IsoMu22_v*");
        //triggerNames.push_back("HLT_IsoTkMu22_v*");
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        //triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");

    } else if (params->selection == "ee") {
        triggerNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
    }

    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    if (true) { // this will need to be turned off for MC
        string jsonFileName = cmssw_base + 
            "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
        //"/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt";
        lumiMask.AddJSONFile(jsonFileName);
    }

    // muon momentum corrections
    muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/rcdata.2016.v3");
    rng = new TRandom3();

    // Prepare the output tree
    string outFileName = params->get_output_filename("output");
    string treeName    = params->get_output_treename("tree");
    string muTreeName  = "muons_" + treeName;
    string evtTreeName = "evt_" + treeName;

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    muTree  = new TTree(muTreeName.c_str(), "muTree");
    evtTree = new TTree(evtTreeName.c_str(), "evtTree");

    // event data
    evtTree->Branch("runNumber", &runNumber);
    evtTree->Branch("evtNumber", &evtNumber, "eventNumber/l");
    evtTree->Branch("lumiSection", &lumiSection);
    evtTree->Branch("triggerStatus", &triggerStatus);
    evtTree->Branch("eventWeight", &eventWeight);
    evtTree->Branch("nPV", &nPV);
    evtTree->Branch("nPU", &nPU);
    evtTree->Branch("nPartons", &nPartons);

    evtTree->Branch("met", &met);
    evtTree->Branch("metPhi", &metPhi);
    evtTree->Branch("nMuons", &nMuons);

    evtTree->Branch("nElectrons", &nElectrons);
    evtTree->Branch("nJets", &nJets);
    evtTree->Branch("nFwdJets", &nFwdJets);
    evtTree->Branch("nBJets", &nBJets);

    // probe muon characteristics
    muTree->Branch("pt",                  &mu_pt);
    muTree->Branch("eta",                 &mu_eta);
    muTree->Branch("phi",                 &mu_phi);
    muTree->Branch("ptErr",               &mu_ptErr);
    muTree->Branch("staPt",               &mu_staPt);
    muTree->Branch("staEta",              &mu_staEta);
    muTree->Branch("staPhi",              &mu_staPhi);
    muTree->Branch("pfPt",                &mu_pfPt);
    muTree->Branch("pfEta",               &mu_pfEta);
    muTree->Branch("pfPhi",               &mu_pfPhi);
    muTree->Branch("trkIso",              &mu_trkIso);
    muTree->Branch("ecalIso",             &mu_ecalIso);
    muTree->Branch("hcalIso",             &mu_hcalIso);
    muTree->Branch("chHadIso",            &mu_chHadIso);
    muTree->Branch("gammaIso",            &mu_gammaIso);
    muTree->Branch("neuHadIso",           &mu_neuHadIso);
    muTree->Branch("puIso",               &mu_puIso);
    muTree->Branch("puppiChHadIso",       &mu_puppiChHadIso);
    muTree->Branch("puppiGammaIso",       &mu_puppiGammaIso);
    muTree->Branch("puppiNeuHadIso",      &mu_puppiNeuHadIso);
    muTree->Branch("puppiChHadIsoNoLep",  &mu_puppiChHadIsoNoLep);
    muTree->Branch("puppiGammaIsoNoLep",  &mu_puppiGammaIsoNoLep);
    muTree->Branch("puppiNeuHadIsoNoLep", &mu_puppiNeuHadIsoNoLep);
    muTree->Branch("d0",                  &mu_d0);
    muTree->Branch("dz",                  &mu_dz);
    muTree->Branch("sip3d",               &mu_sip3d);
    muTree->Branch("tkNchi2",             &mu_tkNchi2);
    muTree->Branch("muNchi2",             &mu_muNchi2);
    muTree->Branch("trkKink",             &mu_trkKink);
    muTree->Branch("glbKink",             &mu_glbKink);
    muTree->Branch("trkHitFrac",          &mu_trkHitFrac);
    muTree->Branch("chi2LocPos",          &mu_chi2LocPos);
    muTree->Branch("segComp",             &mu_segComp);
    muTree->Branch("caloComp",            &mu_caloComp);
    muTree->Branch("q",                   &mu_q);
    muTree->Branch("nValidHits",          &mu_nValidHits);
    muTree->Branch("nTkHits",             &mu_nTkHits);
    muTree->Branch("nPixHits",            &mu_nPixHits);
    muTree->Branch("nTkLayers",           &mu_nTkLayers);
    muTree->Branch("nPixLayers",          &mu_nPixLayers);
    muTree->Branch("nMatchStn",           &mu_nMatchStn);
    muTree->Branch("typeBits",            &mu_typeBits);
    muTree->Branch("selectorBits",        &mu_selectorBits);
    muTree->Branch("pogIDBits",           &mu_pogIDBits);
    muTree->Branch("genPt",               &mu_genPt);
    muTree->Branch("genEta",              &mu_genEta);
    muTree->Branch("genPhi",              &mu_genPhi);
    muTree->Branch("genMatched",          &mu_genMatched);

    // event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);

    ReportPostBegin();
}

Bool_t MuonTagAndProbe::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;
    hTotalEvents->Fill(1);

    if (entry%10000==0) { 
        std::cout << "... Processing event " << entry 
            << " Run: " << fInfo->runNum 
            << " Lumi: " << fInfo->lumiSec 
            << " Event: " << fInfo->evtNum 
            << std::endl;
    }


    const bool isData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isData);

    // Apply lumi mask
    if (isData) {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(2);

    /* Trigger selection */
    bool passTrigger = false;
    for (unsigned i = 0; i < triggerNames.size(); ++i) {
        passTrigger |= trigger->pass(triggerNames[i], fInfo->triggerBits);
    }

    //if (!passTrigger && isData)
    //return kTRUE;

    hTotalEvents->Fill(3);

    /////////////////////
    // Fill event info //
    /////////////////////

    eventWeight   = 1;
    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
    triggerStatus = passTrigger;
    nPV           = fPVArr->GetEntries();
    if (!isData) {
        nPU = fInfo->nPUmean;
        eventWeight *= weights->GetPUWeight(fInfo->nPUmean); // pileup reweighting
    } else {
        nPU = 0;
    }

    ///////////////////////
    // Generator objects //
    ///////////////////////

    vector<TGenParticle> genMuons;
    if (!isData) {
        unsigned count = 0;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            //cout << i << ", "
            //    << particle->status << ", "
            //    << particle->pdgId  << ", "
            //    << particle->parent
            //    << endl;

            if (
                    particle->status == 23 
                    && (abs(particle->pdgId) < 6 || particle->pdgId == 21) 
                    && particle->parent != -1
               ) {
                ++count;
            }

            if (
                    (abs(particle->pdgId) == 13 || abs(particle->pdgId) == 5)
                    && particle->parent > 0
               ) {
                TGenParticle* parent = (TGenParticle*) fGenParticleArr->At(particle->parent);

                if (
                        true 
                        //(abs(parent->pdgId) == 23 || abs(parent->pdgId) == 24)
                        //&& abs(parent->pdgId) != 13 
                        //&& abs(parent->pdgId) != 15
                   ) {
                    genMuons.push_back(*particle);

                    cout << particle->pt << ", " 
                        << particle->eta << ", " 
                        << particle->phi << ", " 
                        << particle->mass << ", " 
                        << particle->pdgId << ", " 
                        << parent->pdgId << ", " 
                        << parent->status
                        << endl;
                }
            }
        }
        nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY
        //cout << nPartons << "\n" << endl;
    } else {
        nPartons = 0;
    }

    cout << "-------------------" << endl;

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
    hTotalEvents->Fill(4);
    }
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);

    /* MUONS */

    //int trigger_muon_index = -1;
    vector<TMuon> muons;
    for (int i=0; i < fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);

        // Apply rochester muon momentum corrections
        double muonSF = 1.;
        if (isData) {
            muonSF = muonCorr->kScaleDT(muon->q, muon->pt, muon->eta, muon->phi, 0, 0);
        } else {
            muonSF = muonCorr->kScaleAndSmearMC(muon->q, muon->pt, muon->eta, muon->phi,
                    muon->nTkLayers, rng->Rndm(), rng->Rndm(), 0, 0);
        }
        muon->pt = muonSF*muon->pt;

        if (
                muon->pt > 3 
                && fabs(muon->eta) < 2.4

                // Loose muon ID
                && (muon->typeBits & baconhep::kPFMuon)
                && ((muon->typeBits & baconhep::kGlobal) || (muon->typeBits & baconhep::kTracker))

                // tight muon ID
                //&& (muon->typeBits & baconhep::kPFMuon) 
                //&& (muon->typeBits & baconhep::kGlobal) 
                //&& muon->muNchi2    < 10.
                //&& muon->nMatchStn  > 1
                //&& muon->nPixHits   > 0
                //&& fabs(muon->d0)   < 0.2
                //&& fabs(muon->dz)   < 0.5
                //&& muon->nTkLayers  > 5 
                //&& muon->nValidHits > 0
           ) {
            muons.push_back(*muon);
        }

        // trigger matching
        //bool triggered = false;
        //for (unsigned i = 0; i < triggerNames.size(); ++i) {
        //    triggered |= trigger->passObj(triggerNames[i], 1, muon->hltMatchBits);
        //}
        //muons_trigger.push_back(triggered);
    }
    //sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);

    /* ELECTRONS */
    std::vector<TLorentzVector> electrons;
    vector<float> electrons_iso;
    vector<float> electrons_q;
    vector<bool> electrons_trigger;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        if (
                electron->pt > 20 
                && fabs(electron->eta) < 2.5
                && particleSelector->PassElectronID(electron, cuts->tightElID)
                && particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
           ) {
            TLorentzVector electronP4;
            copy_p4(electron, ELE_MASS, electronP4);
            electrons.push_back(electronP4);
            electrons_iso.push_back(0.);
            electrons_q.push_back(electron->q);

            // trigger matching
            bool triggered = false;
            for (unsigned i = 0; i < triggerNames.size(); ++i) {
                triggered |= trigger->passObj(triggerNames[i], 1, electron->hltMatchBits);
            }
            electrons_trigger.push_back(triggered);
        }
    }

    std::sort(electrons.begin(), electrons.end(), P4SortCondition);


    /* JETS */
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets;
    std::vector<TJet*> fwdjets;
    std::vector<TJet*> bjets;
    std::vector<TJet*> genbjets;
    std::vector<TJet*> genjets;
    nJets    = 0;
    nFwdJets = 0;
    nBJets   = 0;

    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

        //if (isData) { // fix for broken bacon JEC
        //    double jec = particleSelector->JetCorrector(jet, "NONE");
        //    jet->pt = jet->ptRaw*jec;
        //}

        // Prevent overlap of muons and jets
        TLorentzVector vJet; 
        vJet.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
        bool muOverlap = false;
        //for (const auto& mu: veto_muons) {
        //    if (vJet.DeltaR(mu) < 0.5) {
        //        muOverlap = true;
        //        break;
        //    }
        //}
        //bool elOverlap = false;
        //for (const auto& el: electrons) {
        //    if (vJet.DeltaR(el) < 0.5) {
        //        elOverlap = true;
        //        break;
        //    }
        //}

        //cout << jet->genpt << ", " << jet->genm << ", " 
        //     << jet->geneta << ", " << jet->genphi << ", " 
        //     << jet->partonFlavor << ", " << jet->hadronFlavor << ", " 
        //     << endl;

        if (!isData) {
            if (abs(jet->hadronFlavor) == 5) {
                genbjets.push_back(jet);
            } else {
                genjets.push_back(jet);
            }
        }

        if (
                jet->pt > 30 
                && fabs(jet->eta) < 4.7
                && particleSelector->PassJetID(jet, cuts->looseJetID)
           ) {

            if (fabs(jet->eta) <= 2.4) { 
                if (
                        jet->pt > 30 
                        && jet->mva > -0.89
                        && !muOverlap 
                        //&& !elOverlap
                   ) { 
                    if (isData) {
                        if (jet->bmva > 0.9432) { 
                            bjets.push_back(jet);
                            ++nBJets;
                        } else {
                            jets.push_back(jet);
                            ++nJets;
                        }
                    } else {
                        if (particleSelector->BTagModifier(jet, "MVAT")) { 
                            bjets.push_back(jet);
                            ++nBJets;
                        } else {
                            jets.push_back(jet);
                            ++nJets;
                        }
                    }
                }
            } else {
                if ((fabs(jet->eta) < 2.5 && jet->mva > -0.89) || fabs(jet->eta) > 2.5) {
                    fwdjets.push_back(jet);
                    ++nFwdJets;
                }
            }
        }
    }

    std::sort(fwdjets.begin(), fwdjets.end(), sort_by_higher_pt<TJet>);
    std::sort(bjets.begin(), bjets.end(), sort_by_higher_pt<TJet>);
    std::sort(genjets.begin(), genjets.end(), sort_by_higher_pt<TJet>);
    std::sort(genbjets.begin(), genbjets.end(), sort_by_higher_pt<TJet>);

    // Add additional b jets to the central jet collection
    if (bjets.size() > 1) {
        jets.insert(jets.end(), bjets.begin()+1, bjets.end());
    }
    std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);

    /* MET */
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;

    if (!isData) {
        //met = MetKluge(met)*met;
        met = 0.96*met;
    }

    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    nMuons     = muons.size();
    nElectrons = electrons.size();
    if (params->selection == "mc") {

        // Match selected muons to truth muons
        for (unsigned i = 0; i < muons.size(); ++i) {
            TLorentzVector muon(1., 1., 1., 1.); 
            muon.SetPtEtaPhiM(muons[i].pt, muons[i].eta, muons[i].phi, 0.1057); 

            bool matched     = false;
            bool antimatched = true;
            TLorentzVector genMuonP4(1., 1., 1., 1.);
            for (unsigned j = 0; j < genMuons.size(); ++j) {
                TLorentzVector gmuon(1., 1., 1., 1.); 
                gmuon.SetPtEtaPhiM(genMuons[j].pt, genMuons[j].eta, genMuons[j].phi, 0.1057); 

                if (muon.DeltaR(gmuon) < 0.1) {
                    genMuonP4   = gmuon;
                    matched     = true;
                    antimatched = false;
                    break;
                } else if (muon.DeltaR(gmuon) < 0.3) {
                    antimatched = false;
                } 
            }

            if (matched) {
                FillMuonTree(muons[i], genMuonP4, true);
            } else if (antimatched) {
                FillMuonTree(muons[i], genMuonP4, false);
            }
        }

    }

    ///////////////////
    // Fill jet info //
    ///////////////////


    evtTree->Fill();
    this->passedEvents++;
    return kTRUE;
}

void MuonTagAndProbe::FillMuonTree(TMuon& muon, TLorentzVector& genMuonP4, bool isReal)
{
    // reco level quantities
    mu_pt                  = muon.pt;
    mu_eta                 = muon.eta;
    mu_phi                 = muon.phi;
    mu_ptErr               = muon.ptErr;
    mu_staPt               = muon.staPt;
    mu_staEta              = muon.staEta;
    mu_staPhi              = muon.staPhi;
    mu_pfPt                = muon.pfPt;
    mu_pfEta               = muon.pfEta;
    mu_pfPhi               = muon.pfPhi;
    mu_trkIso              = muon.trkIso;
    mu_ecalIso             = muon.ecalIso;
    mu_hcalIso             = muon.hcalIso;
    mu_chHadIso            = muon.chHadIso;
    mu_gammaIso            = muon.gammaIso;
    mu_neuHadIso           = muon.neuHadIso;
    mu_puIso               = muon.puIso;
    mu_puppiChHadIso       = muon.puppiChHadIso;
    mu_puppiGammaIso       = muon.puppiGammaIso;
    mu_puppiNeuHadIso      = muon.puppiNeuHadIso;
    mu_puppiChHadIsoNoLep  = muon.puppiChHadIsoNoLep;
    mu_puppiGammaIsoNoLep  = muon.puppiGammaIsoNoLep;
    mu_puppiNeuHadIsoNoLep = muon.puppiNeuHadIsoNoLep;
    mu_d0                  = muon.d0;
    mu_dz                  = muon.dz;
    mu_sip3d               = muon.sip3d;
    mu_tkNchi2             = muon.tkNchi2;
    mu_muNchi2             = muon.muNchi2;
    mu_trkKink             = muon.trkKink;
    mu_glbKink             = muon.glbKink;
    mu_trkHitFrac          = muon.trkHitFrac;
    mu_chi2LocPos          = muon.chi2LocPos;
    mu_segComp             = muon.segComp;
    mu_caloComp            = muon.caloComp;
    mu_q                   = muon.q;
    mu_nValidHits          = muon.nValidHits;
    mu_nTkHits             = muon.nTkHits;
    mu_nPixHits            = muon.nPixHits;
    mu_nTkLayers           = muon.nTkLayers;
    mu_nPixLayers          = muon.nPixLayers;
    mu_nMatchStn           = muon.nMatchStn;
    mu_typeBits            = muon.typeBits;
    mu_selectorBits        = muon.selectorBits;
    mu_pogIDBits           = muon.pogIDBits;

    // gen level quantities
    mu_genMatched          = isReal;
    mu_pfPt                = genMuonP4.Pt();
    mu_pfEta               = genMuonP4.Eta();
    mu_pfPhi               = genMuonP4.Phi();

    muTree->Fill();
    return;
}

void MuonTagAndProbe::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void MuonTagAndProbe::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void MuonTagAndProbe::ReportPostTerminate()
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
    std::unique_ptr<MuonTagAndProbe> selector(new MuonTagAndProbe());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

