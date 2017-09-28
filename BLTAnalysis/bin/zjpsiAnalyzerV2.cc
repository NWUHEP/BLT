#include "zjpsiAnalyzerV2.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 

zjpsiAnalyzerV2::zjpsiAnalyzerV2() : BLTSelector()
{

}

zjpsiAnalyzerV2::~zjpsiAnalyzerV2()
{

}

void zjpsiAnalyzerV2::Begin(TTree *tree)
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

    if (params->selection == "mumu" || params->selection == "emu" || params->selection == "4l" || params->selection == "4mu" || params->selection == "4mu_otman") {
        //triggerNames.push_back("HLT_IsoMu22_v*");
        //triggerNames.push_back("HLT_IsoTkMu22_v*");
        //triggerNames.push_back("HLT_IsoMu24_v*");
        //triggerNames.push_back("HLT_IsoTkMu24_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");

    } else if (params->selection == "ee") {
        triggerNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
    
    } else if (params->selection == "jpsi_control") {
        triggerNames.push_back("HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing_v*");
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
    string outTreeName = params->get_output_treename("tree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "bltTree");

    // event data
    outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber, "eventNumber/l");
    outTree->Branch("lumiSection", &lumiSection);
    outTree->Branch("triggerStatus", &triggerStatus);
    outTree->Branch("eventWeight", &eventWeight);
    outTree->Branch("nPV", &nPV);
    outTree->Branch("nPU", &nPU);
    outTree->Branch("nPartons", &nPartons);
    
    outTree->Branch("rPV", &rPV);
    outTree->Branch("rPVErr", &rPVErr);
    outTree->Branch("rPVChi2", &rPVChi2);
    outTree->Branch("rPVNdof", &rPVNdof);
    
    outTree->Branch("rDimuon", &rDimuon);
    outTree->Branch("rDimuonErr", &rDimuonErr);
    outTree->Branch("rDimuonChi2", &rDimuonChi2);
    outTree->Branch("rDimuonNdof", &rDimuonNdof);
    
    outTree->Branch("rZCand", &rZCand);
    outTree->Branch("rZCandErr", &rZCandErr);
    outTree->Branch("rZCandChi2", &rZCandChi2);
    outTree->Branch("rZCandNdof", &rZCandNdof);
    outTree->Branch("rZCandProb", &rZCandProb);
    outTree->Branch("rZCandRxy", &rZCandRxy);
    outTree->Branch("rZCandRxyErr", &rZCandRxyErr);
    
    outTree->Branch("rJpsiCand", &rJpsiCand);
    outTree->Branch("rJpsiCandErr", &rJpsiCandErr);
    outTree->Branch("rJpsiCandChi2", &rJpsiCandChi2);
    outTree->Branch("rJpsiCandNdof", &rJpsiCandNdof);
    outTree->Branch("rJpsiCandProb", &rJpsiCandProb);
    outTree->Branch("rJpsiCandRxy", &rJpsiCandRxy);
    outTree->Branch("rJpsiCandRxyErr", &rJpsiCandRxyErr);

    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);

    // leptons
    outTree->Branch("leptonOneP4", &leptonOneP4);
    outTree->Branch("leptonOneIso", &leptonOneIso);
    outTree->Branch("leptonOneTkIsoNR", &leptonOneTkIsoNR);
    outTree->Branch("leptonOneTkIso", &leptonOneTkIso);
    outTree->Branch("leptonOnePFIso", &leptonOnePFIso);
    outTree->Branch("leptonOneQ", &leptonOneQ);
    outTree->Branch("leptonOneFlavor", &leptonOneFlavor);
    outTree->Branch("leptonOneTrigger", &leptonOneTrigger);
    outTree->Branch("leptonOneIsTight", &leptonOneIsTight);

    outTree->Branch("leptonOneD0", &leptonOneD0);
    outTree->Branch("leptonOneDZ", &leptonOneDZ);

    outTree->Branch("leptonTwoP4", &leptonTwoP4);
    outTree->Branch("leptonTwoIso", &leptonTwoIso);
    outTree->Branch("leptonTwoTkIsoNR", &leptonTwoTkIsoNR);
    outTree->Branch("leptonTwoTkIso", &leptonTwoTkIso);
    outTree->Branch("leptonTwoPFIso", &leptonTwoPFIso);
    outTree->Branch("leptonTwoQ", &leptonTwoQ);
    outTree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
    outTree->Branch("leptonTwoTrigger", &leptonTwoTrigger);
    outTree->Branch("leptonTwoIsTight", &leptonTwoIsTight);
    
    outTree->Branch("leptonTwoD0", &leptonTwoD0);
    outTree->Branch("leptonTwoDZ", &leptonTwoDZ);
    
    outTree->Branch("leptonThreeP4", &leptonThreeP4);
    outTree->Branch("leptonThreeIso", &leptonThreeIso);
    outTree->Branch("leptonThreeTkIsoNR", &leptonThreeTkIsoNR);
    outTree->Branch("leptonThreeTkIso", &leptonThreeTkIso);
    outTree->Branch("leptonThreePFIso", &leptonThreePFIso);
    outTree->Branch("leptonThreeQ", &leptonThreeQ);
    outTree->Branch("leptonThreeFlavor", &leptonThreeFlavor);
    outTree->Branch("leptonThreeTrigger", &leptonThreeTrigger);
    outTree->Branch("leptonThreeIsTight", &leptonThreeIsTight);

    outTree->Branch("leptonThreeD0", &leptonThreeD0);
    outTree->Branch("leptonThreeDZ", &leptonThreeDZ);

    outTree->Branch("leptonFourP4", &leptonFourP4);
    outTree->Branch("leptonFourIso", &leptonFourIso);
    outTree->Branch("leptonFourTkIsoNR", &leptonFourTkIsoNR);
    outTree->Branch("leptonFourTkIso", &leptonFourTkIso);
    outTree->Branch("leptonFourPFIso", &leptonFourPFIso);
    outTree->Branch("leptonFourQ", &leptonFourQ);
    outTree->Branch("leptonFourFlavor", &leptonFourFlavor);
    outTree->Branch("leptonFourTrigger", &leptonFourTrigger);
    outTree->Branch("leptonFourIsTight", &leptonFourIsTight);
    
    outTree->Branch("leptonFourD0", &leptonFourD0);
    outTree->Branch("leptonFourDZ", &leptonFourDZ);

    // jets
    outTree->Branch("jetP4", &jetP4);
    outTree->Branch("jetD0", &jetD0);
    outTree->Branch("jetTag", &jetTag);
    outTree->Branch("jetPUID", &jetPUID);
    outTree->Branch("jetFlavor", &jetFlavor);

    outTree->Branch("bjetP4", &bjetP4);
    outTree->Branch("bjetD0", &bjetD0);
    outTree->Branch("bjetTag", &bjetTag);
    outTree->Branch("bjetPUID", &bjetPUID);
    outTree->Branch("bjetFlavor", &bjetFlavor);

    outTree->Branch("genBJetP4", &genBJetP4);
    outTree->Branch("genBJetTag", &genBJetTag);
    outTree->Branch("genJetP4", &genJetP4);
    outTree->Branch("genJetTag", &genJetTag);

    outTree->Branch("photonOneP4", &photonOneP4);

    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nFwdJets", &nFwdJets);
    outTree->Branch("nBJets", &nBJets);
    outTree->Branch("nPhotons", &nPhotons);

    // event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",16,0.5,16.5);

    ReportPostBegin();
}

Bool_t zjpsiAnalyzerV2::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
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

    if (!passTrigger && isData)
        return kTRUE;

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

    if (!isData) {
        unsigned count = 0;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            //cout << particle->status << ", "
            //     << particle->pdgId  << ", "
            //     << particle->parent
            //     << endl;

            if (
                    particle->status == 23 
                    && (abs(particle->pdgId) < 6 || particle->pdgId == 21) 
                    && particle->parent != -2
               ) {
                ++count;

            }
        }
        nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY
        //cout << nPartons << "\n" << endl;
    } else {
        nPartons = 0;
    }

    ///////////////////
    // Select objects//
    ///////////////////

    /* Vertices */
    if (fInfo->hasGoodPV) {
        assert(fPVArr->GetEntries() != 0);
        TVertex* primary = (TVertex*) fPVArr->At(0);
        TVector3 pv;
        TVector3 pv_err;  
        pv.SetXYZ(primary->x, primary->y, primary->z);
        pv_err.SetXYZ(primary->xerr, primary->yerr, primary->zerr);
        rPV = pv;
        rPVErr = pv_err;
        rPVChi2 = primary->chi2;
        rPVNdof = primary->ndof;
        particleSelector->SetPV(pv);
    } else {
        return kTRUE;
    }
    hTotalEvents->Fill(4);
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);

    /* MUONS */
    vector<TMuon*> muons;

    /* Apply Rochester Muon Corrections */
    for (int i=0; i < fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);
        
        double muonSF = 1.;
        if (isData) {
            muonSF = muonCorr->kScaleDT(muon->q, muon->pt, muon->eta, muon->phi, 0, 0);
        } 
        else {
            muonSF = muonCorr->kScaleAndSmearMC(muon->q, muon->pt, muon->eta, muon->phi,
                                                muon->nTkLayers, rng->Rndm(), rng->Rndm(), 
                                                0, 0);
        }
        muon->pt = muonSF*muon->pt;

        if (
                    muon->pt > 3.0
                    && fabs(muon->eta) < 2.4

                    // loose muon ID
                    && (muon->typeBits & baconhep::kPFMuon)
                    && ((muon->typeBits & baconhep::kGlobal) || (muon->selectorBits & baconhep::kTrackerMuonArbitrated))
                    //&& fabs(muon->d0) < 0.5
                    //&& fabs(muon->dz) < 1.0
                    && fabs(muon->d0) < 0.5
                    && fabs(muon->dz) < 20.0
                    //&& muon->sip3d < 4.0
                ){
                    muons.push_back(muon);
        }
                                     
    }

    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);

    // Getting muon info
    vector<TLorentzVector> muons_P4;
    vector<float> muons_iso;
    vector<float> muons_TkIsoNR;
    vector<float> muons_TkIso;
    vector<float> muons_PFIso;
    vector<bool> muons_isTight;
    vector<float> muons_q;
    vector<bool> muons_trigger;
    vector<float> muons_d0;
    vector<float> muons_dz;
    vector<unsigned int> muons_index;
    //vector<map<int, TVector3>> muons_vertices;

    for (unsigned i = 0; i < muons.size(); i++) {
        TMuon* muon = muons[i];
        TLorentzVector muonP4;
        copy_p4(muons[i], MUON_MASS, muonP4);

        // Fill containers

        muons_P4.push_back(muonP4);
        muons_q.push_back(muon->q);
        muons_d0.push_back(muon->d0);
        muons_dz.push_back(muon->dz);
        muons_index.push_back(muon->muIndex);
     //   muons_vertices.push_back(muon->dimuonVertex);

        // ID and isolation
        if ( 
                // tight muon ID
                (muon->typeBits & baconhep::kPFMuon) 
                && (muon->typeBits & baconhep::kGlobal) 
                && muon->muNchi2    < 10.
                && muon->nMatchStn  > 1
                && muon->nPixHits   > 0
                && fabs(muon->d0)   < 0.2
                && fabs(muon->dz)   < 0.5
                && muon->nTkLayers  > 5 
                && muon->nValidHits > 0

           ) {
                muons_isTight.push_back(true);
        } 
        else 
            muons_isTight.push_back(false);

        float pfiso = muon->chHadIso + max(0., muon->neuHadIso + muon->gammaIso -0.5*(muon->puIso));
        muons_iso.push_back(pfiso); 
        muons_PFIso.push_back(pfiso);
        muons_TkIsoNR.push_back(muon->trkIso);

        // Tracker isolation Pt removal 
        // Remove muon track pt from muon track isolation variable
        //for (unsigned j = i+1; j < muons.size(); j++) {
        for (unsigned j = 0; j < muons.size(); j++) {
            if (i!=j) {
                TLorentzVector muon_j;
                copy_p4(muons[j], MUON_MASS, muon_j);

                if (muonP4.DeltaR(muon_j) < 0.3) {
                    muon->trkIso = max(0., muon->trkIso - muon_j.Pt());
                    //muons[j]->trkIso = max(0., muons[j]->trkIso - muonP4.Pt());
                }
            }
        }

        muons_TkIso.push_back(muon->trkIso);

        // trigger matching
        bool triggered = false;
        for (unsigned i = 0; i < triggerNames.size(); ++i) {
            triggered |= trigger->passObj(triggerNames[i], 1, muon->hltMatchBits);
        }
        muons_trigger.push_back(triggered);
    }
    std::sort(muons_P4.begin(), muons_P4.end(), P4SortCondition);
   

    /* ELECTRONS */
    std::vector<TLorentzVector> electrons;
    vector<float> electrons_iso;
    vector<float> electrons_q;
    vector<bool> electrons_trigger;
    vector<float> electrons_d0;
    vector<float> electrons_dz;
    vector<unsigned int> electrons_index;
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
            electrons_d0.push_back(electron->d0);
            electrons_dz.push_back(electron->dz);
            electrons_index.push_back(electron->eleIndex);

            // trigger matching
            bool triggered = false;
            for (unsigned i = 0; i < triggerNames.size(); ++i) {
                triggered |= trigger->passObj(triggerNames[i], 1, electron->hltMatchBits);
            }
            electrons_trigger.push_back(triggered);
        }
    }

    std::sort(electrons.begin(), electrons.end(), P4SortCondition);

    
    /* PHOTONS */
    std::vector<TLorentzVector> photons;
    for (int i=0; i<fPhotonArr->GetEntries(); i++) {
        TPhoton* photon = (TPhoton*) fPhotonArr->At(i);
        assert(photon);

        if (
                photon->pt > 10
                && fabs(photon->eta) < 2.5
                //&& particleSelector->PassPhotonID(photon, cuts->loosePhID)
                && particleSelector->PassPhotonID(photon, cuts->preSelPhID)
                && particleSelector->PassPhotonIso(photon, cuts->loosePhIso, cuts->EAPho)
           ) {
            TLorentzVector photonP4;
            copy_p4(photon, 0., photonP4);
            photons.push_back(photonP4);
        }
    }

    std::sort(photons.begin(), photons.end(), P4SortCondition);


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

        if (isData) { // fix for broken bacon JEC
            double jec = particleSelector->JetCorrector(jet, "NONE");
            jet->pt = jet->ptRaw*jec;
        }

        // Prevent overlap of muons and jets
       /* TLorentzVector vJet; 
        vJet.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
        bool muOverlap = false;
        for (const auto& mu: veto_muons) {
            if (vJet.DeltaR(mu) < 0.5) {
                muOverlap = true;
                break;
            }
        }*/
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
                        //&& !muOverlap 
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
    nElectrons  = electrons.size();
    nPhotons = photons.size();

    if (params->selection == "mumu") {
        if (muons_P4.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        // Select the leading muon
        TLorentzVector muonOneP4;
        unsigned muonOneIndex = 0;
        bool validMuon1 = false;
        for (unsigned i = 0; i < (muons_P4.size() - 1); ++i) {
            if (
                    //muons_isTight[i] 
                    //&& (muons_iso[i]/muons_P4[i].Pt()) < 0.15
                    /*&& */muons_P4[i].Pt() > 4.0
                )
            {
                muonOneP4 = muons_P4[i];
                muonOneIndex = i;
                validMuon1 = true;
                break;
            }
        }
        
        if (!validMuon1)
            return kTRUE;
        hTotalEvents->Fill(6);

        // Select the oppositely charged subleading muon
        TLorentzVector muonTwoP4;
        unsigned muonTwoIndex = muonOneIndex + 1;
        bool validMuon2 = false;
        for (unsigned i = muonOneIndex+1; i < muons_P4.size(); ++i) {
            if (
                    muons_q[muonOneIndex] != muons_q[i]
                    //&& muons_isTight[i]
                    //&& (muons_iso[i]/muons_P4[i].Pt()) < 0.15
                    && muons_P4[i].Pt() > 4.0
               )
            {
                muonTwoP4 = muons_P4[i];
                muonTwoIndex = i;
                validMuon2 = true;
                break;
            }
        }

        if (!validMuon2)
            return kTRUE;

        hTotalEvents->Fill(7);

        leptonOneP4      = muonOneP4;
        //leptonOneIso     = muons_iso[muonOneIndex];
        leptonOneQ       = muons_q[muonOneIndex];
        leptonOneTrigger = muons_trigger[muonOneIndex];
        leptonOneIsTight = muons_isTight[muonOneIndex];
        leptonOneFlavor  = 1;

        leptonOneD0        = muons_d0[muonOneIndex];
        leptonOneDZ      = muons_dz[muonOneIndex];
        leptonTwoP4      = muonTwoP4;
        //leptonTwoIso     = muons_iso[muonTwoIndex];
        leptonTwoQ       = muons_q[muonTwoIndex];
        leptonTwoTrigger = muons_trigger[muonTwoIndex];
        leptonTwoIsTight = muons_isTight[muonTwoIndex];
        leptonTwoFlavor  = 1;
            
        leptonTwoD0        = muons_d0[muonTwoIndex];
        leptonTwoDZ        = muons_dz[muonTwoIndex];

        // Check for dimuon vertex
        unsigned int leptonOneIndex = muons_index[muonOneIndex];
        unsigned int leptonTwoIndex = muons_index[muonTwoIndex];

        bool hasValidVertex = false;
         
        for (int i = 0; i < fDimuonVertexArr->GetEntries(); ++i) {
            TVertex* dimuonVert = (TVertex*) fDimuonVertexArr->At(i);
            if (
                    ((dimuonVert->index1 == leptonOneIndex && dimuonVert->index2 == leptonTwoIndex) ||
                    (dimuonVert->index1 == leptonTwoIndex && dimuonVert->index2 == leptonOneIndex)) &&
                    (dimuonVert->isValid)
               ) {
                    TVector3 vertPos;
                    TVector3 vertErr;  
                    vertPos.SetXYZ(dimuonVert->x, dimuonVert->y, dimuonVert->z);
                    vertErr.SetXYZ(dimuonVert->xerr, dimuonVert->yerr, dimuonVert->zerr);
                    rDimuon = vertPos; 
                    rDimuonErr = vertErr;
                    rDimuonChi2 = dimuonVert->chi2;
                    rDimuonNdof = dimuonVert->ndof;
                    hasValidVertex = true; 
                    break;
            }
        }

        if (!hasValidVertex)
            return kTRUE; 
 
        hTotalEvents->Fill(8);

        if (photons.size() > 0)
            photonOneP4 = photons[0];

        if (!isData) {
            eventWeight *= weights->GetMuonIDEff(muonOneP4);
            eventWeight *= weights->GetMuonISOEff(muonOneP4);
            eventWeight *= weights->GetMuonIDEff(muonTwoP4);
            eventWeight *= weights->GetMuonISOEff(muonTwoP4);

            // trigger weight
            pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
            pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
            eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);
        }

    }

    else if (params->selection == "4mu") {
        if (muons_P4.size() < 4)
            return kTRUE;
        hTotalEvents->Fill(5);
        
        TLorentzVector muonOneP4;
        unsigned muonOneIndex = 0;
        bool validMuon1 = false;
        TLorentzVector muonTwoP4;
        unsigned muonTwoIndex = muonOneIndex + 1;
        bool validMuon2 = false;

        unsigned int pairing_algo = 1;

        if (pairing_algo == 0) {
        
            // Select the leading muon
            for (unsigned i = 0; i < (muons_P4.size() - 1); ++i) {
                if (
                        muons_isTight[i] 
                        && (muons_iso[i]/muons_P4[i].Pt()) < 0.15
                       // && (muons_TkIso[i]/muons_P4[i].Pt()) < 0.05
                        && muons_P4[i].Pt() > 22.0
                    )
                {
                    muonOneP4 = muons_P4[i];
                    muonOneIndex = i;
                    validMuon1 = true;
                    break;
                }
            }

            // Select the oppositely charged second muon
            for (unsigned i = muonOneIndex+1; i < muons_P4.size(); ++i) {
                if (
                        muons_q[muonOneIndex] != muons_q[i]
                        && muons_isTight[i]
                        && (muons_iso[i]/muons_P4[i].Pt()) < 0.15
                        //&& (muons_TkIso[i]/muons_P4[i].Pt()) < 0.05
                        && muons_P4[i].Pt() > 10.0
                   )
                {
                    muonTwoP4 = muons_P4[i];
                    muonTwoIndex = i;
                    validMuon2 = true;
                    break;
                }
            }
        }

        else if (pairing_algo == 1) {
            float zMass_pdg = 91.188;
            float zMassDiff = 100.;
            for (unsigned i = 0; i < muons_P4.size(); ++i) {
                for (unsigned j = i + 1; j < muons_P4.size(); ++j) {
                    float thisMass = (muons_P4[i] + muons_P4[j]).M();
                    if (
                            fabs(thisMass - zMass_pdg) < zMassDiff
                            && muons_isTight[i]
                            && (muons_iso[i]/muons_P4[i].Pt()) < 0.15
                            //&& (muons_TkIso[i]/muons_P4[i].Pt()) < 0.05
                            && muons_P4[i].Pt() > 22.0
                            && muons_isTight[j]
                            && (muons_iso[j]/muons_P4[j].Pt()) < 0.15
                            //&& (muons_TkIso[i]/muons_P4[i].Pt()) < 0.05
                            && muons_P4[j].Pt() > 10.0
                            && muons_q[i] != muons_q[j]
                       ) 
                    {
                        muonOneIndex = i;
                        muonTwoIndex = j;
                        muonOneP4 = muons_P4[i];
                        muonTwoP4 = muons_P4[j];
                        validMuon1 = true;
                        validMuon2 = true;
                        zMassDiff = fabs(thisMass - zMass_pdg);
                    }
                }
            }
        }

        if (!validMuon1)
            return kTRUE;

        hTotalEvents->Fill(6);

        if (!validMuon2)
            return kTRUE;

        hTotalEvents->Fill(7);

        leptonOneP4      = muonOneP4;
        leptonOneIso     = muons_iso[muonOneIndex];
        leptonOneTkIsoNR     = muons_TkIsoNR[muonOneIndex];
        leptonOneTkIso     = muons_TkIso[muonOneIndex];
        leptonOnePFIso     = muons_PFIso[muonOneIndex];
        leptonOneQ       = muons_q[muonOneIndex];
        leptonOneTrigger = muons_trigger[muonOneIndex];
        leptonOneIsTight = muons_isTight[muonOneIndex];
        leptonOneFlavor  = 1;

        leptonOneD0        = muons_d0[muonOneIndex];
        leptonOneDZ      = muons_dz[muonOneIndex];
        leptonTwoP4      = muonTwoP4;
        leptonTwoIso     = muons_iso[muonTwoIndex];
        leptonTwoTkIsoNR     = muons_TkIsoNR[muonTwoIndex];
        leptonTwoTkIso     = muons_TkIso[muonTwoIndex];
        leptonTwoPFIso     = muons_PFIso[muonTwoIndex];
        leptonTwoQ       = muons_q[muonTwoIndex];
        leptonTwoTrigger = muons_trigger[muonTwoIndex];
        leptonTwoIsTight = muons_isTight[muonTwoIndex];
        leptonTwoFlavor  = 1;
            
        leptonTwoD0        = muons_d0[muonTwoIndex];
        leptonTwoDZ        = muons_dz[muonTwoIndex];

        float thisZMass = (muonOneP4 + muonTwoP4).M();
        if (
                thisZMass < 2.0 ||
                (thisZMass > 4.0 && thisZMass < 75.) ||
                thisZMass > 107.
           )
            return kTRUE;
        hTotalEvents->Fill(8);

        // Check for dimuon vertex 
        unsigned int leptonOneIndex = muons_index[muonOneIndex];
        unsigned int leptonTwoIndex = muons_index[muonTwoIndex];

        bool hasValidZVertex = false;
        //bool hasGoodZVertex = false;
         
        for (int i = 0; i < fDimuonVertexArr->GetEntries(); ++i) {
            TVertex* dimuonVert = (TVertex*) fDimuonVertexArr->At(i);
            if (
                    ((dimuonVert->index1 == leptonOneIndex && dimuonVert->index2 == leptonTwoIndex) ||
                    (dimuonVert->index1 == leptonTwoIndex && dimuonVert->index2 == leptonOneIndex)) &&
                    (dimuonVert->isValid)
               ) {
                    TVector3 vertPos;
                    TVector3 vertErr;  
                    vertPos.SetXYZ(dimuonVert->x, dimuonVert->y, dimuonVert->z);
                    vertErr.SetXYZ(dimuonVert->xerr, dimuonVert->yerr, dimuonVert->zerr);
                    rZCand = vertPos; 
                    rZCandErr = vertErr;
                    rZCandChi2 = dimuonVert->chi2;
                    rZCandNdof = dimuonVert->ndof;
                    rZCandProb = dimuonVert->prob;
                    rZCandRxy = dimuonVert->rxy;
                    rZCandRxyErr = dimuonVert->rxy_err;
                    hasValidZVertex = true; 
                    //if (rZCandProb > 0.005)
                    //    hasGoodZVertex = true;
                    break;
            }
        }
        
        if (!hasValidZVertex)
            return kTRUE; 
        
        hTotalEvents->Fill(9);

        //if (!hasGoodZVertex)
        //    return kTRUE;

        //hTotalEvents->Fill(10);

        TLorentzVector muonThreeP4;
        unsigned muonThreeIndex = 0;
        bool validMuon3 = false;
        TLorentzVector muonFourP4;
        unsigned muonFourIndex = muonThreeIndex + 1;
        bool validMuon4 = false;

        if (pairing_algo == 0) {

            // Select the J/Psi candidate leading muon
            for (unsigned i = 0; i < (muons_P4.size() - 1); ++i) {
                if (
                        i != muonOneIndex 
                        && i != muonTwoIndex
                        && (muons_iso[i]/muons_P4[i].Pt()) < 0.25
                        //&& (muons_TkIso[i]/muons_P4[i].Pt()) < 0.10
                        && muons_P4[i].Pt() > 4.0
                    )
                {
                    muonThreeP4 = muons_P4[i];
                    muonThreeIndex = i;
                    validMuon3 = true;
                    break;
                }
            }
            
            // Select the J/Psi candidate subleading muon
            for (unsigned i = muonThreeIndex+1; i < muons_P4.size(); ++i) {
                if (
                        i != muonOneIndex 
                        && i != muonTwoIndex
                        && muons_q[muonThreeIndex] != muons_q[i]
                   )
                {
                    muonFourP4 = muons_P4[i];
                    muonFourIndex = i;
                    validMuon4 = true;
                    break;
                }
            }
        }

        else if (pairing_algo == 1) {
            float jpsiMass_pdg = 3.097;
            float jpsiMassDiff = 100.;
            for (unsigned i = 0; i < muons_P4.size(); ++i) {
                for (unsigned j = i + 1; j < muons_P4.size(); ++j) {
                    float thisMass = (muons_P4[i] + muons_P4[j]).M();
                    if (
                            fabs(thisMass - jpsiMass_pdg) < jpsiMassDiff
                            && !(i == muonOneIndex || i == muonTwoIndex)
                            && !(j == muonOneIndex || j == muonTwoIndex)
                            && (muons_iso[i]/muons_P4[i].Pt()) < 0.25
                            //&& (muons_TkIso[i]/muons_P4[i].Pt()) < 0.10
                            && muons_P4[i].Pt() > 4.0
                            && muons_q[i] != muons_q[j]
                       ) 
                    {
                        muonThreeIndex = i;
                        muonFourIndex = j;
                        muonThreeP4 = muons_P4[i];
                        muonFourP4 = muons_P4[j];
                        validMuon3 = true;
                        validMuon4 = true;
                        jpsiMassDiff = fabs(thisMass - jpsiMass_pdg);
                    }
                }
            }
        }
        
        if (!validMuon3)
            return kTRUE;
       
        hTotalEvents->Fill(10);
        
        if (!validMuon4)
            return kTRUE;

        hTotalEvents->Fill(11);
        
        leptonThreeP4      = muonThreeP4;
        leptonThreeIso     = muons_iso[muonThreeIndex];
        leptonThreeTkIsoNR     = muons_TkIsoNR[muonThreeIndex];
        leptonThreeTkIso     = muons_TkIso[muonThreeIndex];
        leptonThreePFIso     = muons_PFIso[muonThreeIndex];
        leptonThreeQ       = muons_q[muonThreeIndex];
        leptonThreeTrigger = muons_trigger[muonThreeIndex];
        leptonThreeIsTight = muons_isTight[muonThreeIndex];
        leptonThreeFlavor  = 1;

        leptonThreeD0        = muons_d0[muonThreeIndex];
        leptonThreeDZ        = muons_dz[muonThreeIndex];

        leptonFourP4      = muonFourP4;
        leptonFourIso     = muons_iso[muonFourIndex];
        leptonFourTkIsoNR     = muons_TkIsoNR[muonFourIndex];
        leptonFourTkIso     = muons_TkIso[muonFourIndex];
        leptonFourPFIso     = muons_PFIso[muonFourIndex];
        leptonFourQ       = muons_q[muonFourIndex];
        leptonFourTrigger = muons_trigger[muonFourIndex];
        leptonFourIsTight = muons_isTight[muonFourIndex];
        leptonFourFlavor  = 1;
            
        leptonFourD0        = muons_d0[muonFourIndex];
        leptonFourDZ        = muons_dz[muonFourIndex];
        
        float thisJpsiMass = (muonThreeP4 + muonFourP4).M();
        if (
                thisJpsiMass < 2.0 ||
                (thisJpsiMass > 4.0 && thisJpsiMass < 75.) ||
                thisJpsiMass > 107.
           )
            return kTRUE;
        hTotalEvents->Fill(12);
       
        // Check for dimuon vertex 
        unsigned int leptonThreeIndex = muons_index[muonThreeIndex];
        unsigned int leptonFourIndex = muons_index[muonFourIndex];

        bool hasValidJpsiVertex = false;
        //bool hasGoodJpsiVertex = false;
         
        for (int i = 0; i < fDimuonVertexArr->GetEntries(); ++i) {
            TVertex* dimuonVert = (TVertex*) fDimuonVertexArr->At(i);
            if (
                    ((dimuonVert->index1 == leptonThreeIndex && dimuonVert->index2 == leptonFourIndex) ||
                    (dimuonVert->index1 == leptonFourIndex && dimuonVert->index2 == leptonThreeIndex)) &&
                    (dimuonVert->isValid)
               ) {
                    TVector3 vertPos;
                    TVector3 vertErr;  
                    vertPos.SetXYZ(dimuonVert->x, dimuonVert->y, dimuonVert->z);
                    vertErr.SetXYZ(dimuonVert->xerr, dimuonVert->yerr, dimuonVert->zerr);
                    rJpsiCand = vertPos; 
                    rJpsiCandErr = vertErr;
                    rJpsiCandChi2 = dimuonVert->chi2;
                    rJpsiCandNdof = dimuonVert->ndof;
                    rJpsiCandProb = dimuonVert->prob;
                    rJpsiCandRxy = dimuonVert->rxy;
                    rJpsiCandRxyErr = dimuonVert->rxy_err;
                    hasValidJpsiVertex = true; 
                    //if (rJpsiCandProb > 0.005)
                    //    hasGoodJpsiVertex = true;
                    break;
            }
        }
        
        if (!hasValidJpsiVertex)
            return kTRUE;
 
        hTotalEvents->Fill(13);

        //if (!hasGoodJpsiVertex)
        //    return kTRUE;

        hTotalEvents->Fill(14);

        //if (fabs(rJpsiCand.Z() - rZCand.Z()) > 1.)
        //    return kTRUE;

        //hTotalEvents->Fill(16);
        
        if (photons.size() > 0)
            photonOneP4 = photons[0];
    }

    else if (params->selection == "4mu_otman") {
        vector<pair<TLorentzVector, unsigned int>> tightMuons;
        vector<pair<TLorentzVector, unsigned int>> tightIsoMuons;
        vector<pair<TLorentzVector, unsigned int>> looseMuons;
        vector<pair<TLorentzVector, unsigned int>> looseIsoMuons;
        TLorentzVector muonOneP4, muonTwoP4, muonThreeP4, muonFourP4;
        unsigned int muonOneIndex = 0;
        unsigned int muonTwoIndex = 0;
        //unsigned int muonThreeIndex = 0;
        //unsigned int muonFourIndex = 0;
        for (unsigned i = 0; i < muons_P4.size(); ++i) {
            looseMuons.push_back(make_pair(muons_P4[i], i));
            if (muons_iso[i]/muons_P4[i].Pt() < 0.25) 
                looseIsoMuons.push_back(make_pair(muons_P4[i], i));
            if (muons_isTight[i]) {
                tightMuons.push_back(make_pair(muons_P4[i], i));
                if (muons_iso[i]/muons_P4[i].Pt() < 0.15)
                    tightIsoMuons.push_back(make_pair(muons_P4[i], i));
            }
        }
        if (tightMuons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);
        if (tightIsoMuons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(6);
        float zMass_pdg = 91.188;
        bool inZWindow = false;
        for (unsigned i = 0; i < tightIsoMuons.size(); ++i) {
            for (unsigned j = i + 1; j < tightIsoMuons.size(); ++j) {
                float thisDimuonMass = (tightIsoMuons[i].first + tightIsoMuons[j].first).M();
                if (fabs(thisDimuonMass - zMass_pdg) < 5.0) {
                    inZWindow = true;
                    muonOneP4 = tightIsoMuons[i].first;
                    muonOneIndex = tightIsoMuons[i].second;
                    muonTwoP4 = tightIsoMuons[j].first;
                    muonTwoIndex = tightIsoMuons[j].second;
                } 
            }
        }
        if (!inZWindow)
            return kTRUE;
        hTotalEvents->Fill(7);
        if (looseMuons.size() < 4)
            return kTRUE;
        hTotalEvents->Fill(8);
        if (looseIsoMuons.size() < 4)
            return kTRUE;
        hTotalEvents->Fill(9);
        float jpsiMass_pdg = 3.097;
        bool inJpsiWindow = false;
        for (unsigned i = 0; i < looseIsoMuons.size(); ++i) {
            for (unsigned j = i + 1; j < looseIsoMuons.size(); ++j) {
                float thisDimuonMass = (looseIsoMuons[i].first + looseIsoMuons[j].first).M();
                if (
                        !(looseIsoMuons[i].second == muonOneIndex || looseIsoMuons[i].second == muonTwoIndex) &&
                        !(looseIsoMuons[j].second == muonOneIndex || looseIsoMuons[j].second == muonTwoIndex) &&
                        fabs(thisDimuonMass - jpsiMass_pdg) < 0.5
                   ) {
                    inJpsiWindow = true;
                    muonThreeP4 = looseIsoMuons[i].first;
                    //muonThreeIndex = looseIsoMuons[i].second;
                    muonFourP4 = looseIsoMuons[j].first;
                    //muonFourIndex = looseIsoMuons[j].second;
                }
            }
        }
        if (!inJpsiWindow)
            return kTRUE;
        hTotalEvents->Fill(10);
        if (
                muonOneP4.Pt() < 22.0 ||
                muonTwoP4.Pt() < 10.0 ||
                muonThreeP4.Pt() < 4.0 ||
                muonFourP4.Pt() < 3.0
           )
            return kTRUE;
        hTotalEvents->Fill(11);
    }
    
    else if (params->selection == "jpsi_control") {
        if (muons_P4.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);
        
        TLorentzVector muonOneP4;
        unsigned muonOneIndex = 0;
        bool validMuon1 = false;
        TLorentzVector muonTwoP4;
        unsigned muonTwoIndex = muonOneIndex + 1;
        bool validMuon2 = false;

        unsigned int pairing_algo = 1;
        
        if (pairing_algo == 0) {

            // Select the J/Psi candidate leading muon
            for (unsigned i = 0; i < (muons_P4.size() - 1); ++i) {
                if (
                        (muons_iso[i]/muons_P4[i].Pt()) < 0.25
                        && muons_P4[i].Pt() > 4.0
                    )
                {
                    muonOneP4 = muons_P4[i];
                    muonOneIndex = i;
                    validMuon1 = true;
                    break;
                }
            }
            
            // Select the J/Psi candidate subleading muon
            for (unsigned i = muonOneIndex+1; i < muons_P4.size(); ++i) {
                if (
                        muons_q[muonOneIndex] != muons_q[i]
                   )
                {
                    muonTwoP4 = muons_P4[i];
                    muonTwoIndex = i;
                    validMuon2 = true;
                    break;
                }
            }
        }

        else if (pairing_algo == 1) {
            float jpsiMass_pdg = 3.097;
            float jpsiMassDiff = 100.;
            for (unsigned i = 0; i < muons_P4.size(); ++i) {
                for (unsigned j = i + 1; j < muons_P4.size(); ++j) {
                    float thisMass = (muons_P4[i] + muons_P4[j]).M();
                    if (
                            fabs(thisMass - jpsiMass_pdg) < jpsiMassDiff
                            && (muons_iso[i]/muons_P4[i].Pt()) < 0.25
                            && muons_P4[i].Pt() > 4.0
                            && muons_q[i] != muons_q[j]
                       ) 
                    {
                        muonOneIndex = i;
                        muonTwoIndex = j;
                        muonOneP4 = muons_P4[i];
                        muonTwoP4 = muons_P4[j];
                        validMuon1 = true;
                        validMuon2 = true;
                        jpsiMassDiff = fabs(thisMass - jpsiMass_pdg);
                    }
                }
            }
        }


        leptonOneP4      = muonOneP4;
        leptonOneIso     = muons_iso[muonOneIndex];
        leptonOneTkIsoNR     = muons_TkIsoNR[muonOneIndex];
        leptonOneTkIso     = muons_TkIso[muonOneIndex];
        leptonOnePFIso     = muons_PFIso[muonOneIndex];
        leptonOneQ       = muons_q[muonOneIndex];
        leptonOneTrigger = muons_trigger[muonOneIndex];
        leptonOneIsTight = muons_isTight[muonOneIndex];
        leptonOneFlavor  = 1;

        leptonOneD0        = muons_d0[muonOneIndex];
        leptonOneDZ      = muons_dz[muonOneIndex];
        leptonTwoP4      = muonTwoP4;
        leptonTwoIso     = muons_iso[muonTwoIndex];
        leptonTwoTkIsoNR     = muons_TkIsoNR[muonTwoIndex];
        leptonTwoTkIso     = muons_TkIso[muonTwoIndex];
        leptonTwoPFIso     = muons_PFIso[muonTwoIndex];
        leptonTwoQ       = muons_q[muonTwoIndex];
        leptonTwoTrigger = muons_trigger[muonTwoIndex];
        leptonTwoIsTight = muons_isTight[muonTwoIndex];
        leptonTwoFlavor  = 1;
            
        leptonTwoD0        = muons_d0[muonTwoIndex];
        leptonTwoDZ        = muons_dz[muonTwoIndex];

        
        if (!validMuon1)
            return kTRUE;
       
        hTotalEvents->Fill(6);
        
        if (!validMuon2)
            return kTRUE;

        hTotalEvents->Fill(7);
        
        leptonOneP4      = muonOneP4;
        leptonOneIso     = muons_iso[muonOneIndex];
        leptonOneTkIsoNR     = muons_TkIsoNR[muonOneIndex];
        leptonOneTkIso     = muons_TkIso[muonOneIndex];
        leptonOnePFIso     = muons_PFIso[muonOneIndex];
        leptonOneQ       = muons_q[muonOneIndex];
        leptonOneTrigger = muons_trigger[muonOneIndex];
        leptonOneIsTight = muons_isTight[muonOneIndex];
        leptonOneFlavor  = 1;

        leptonOneD0        = muons_d0[muonOneIndex];
        leptonOneDZ        = muons_dz[muonOneIndex];

        leptonTwoP4      = muonTwoP4;
        leptonTwoIso     = muons_iso[muonTwoIndex];
        leptonTwoTkIsoNR     = muons_TkIsoNR[muonTwoIndex];
        leptonTwoTkIso     = muons_TkIso[muonTwoIndex];
        leptonTwoPFIso     = muons_PFIso[muonTwoIndex];
        leptonTwoQ       = muons_q[muonTwoIndex];
        leptonTwoTrigger = muons_trigger[muonTwoIndex];
        leptonTwoIsTight = muons_isTight[muonTwoIndex];
        leptonTwoFlavor  = 1;
            
        leptonTwoD0        = muons_d0[muonTwoIndex];
        leptonTwoDZ        = muons_dz[muonTwoIndex];
        
        float thisJpsiMass = (muonOneP4 + muonTwoP4).M();
        if (
                thisJpsiMass < 2.0 ||
                (thisJpsiMass > 4.0 && thisJpsiMass < 75.) ||
                thisJpsiMass > 107.
           )
            return kTRUE;
        hTotalEvents->Fill(8);
       
        // Check for dimuon vertex 
        unsigned int leptonOneIndex = muons_index[muonOneIndex];
        unsigned int leptonTwoIndex = muons_index[muonTwoIndex];

        bool hasValidJpsiVertex = false;
         
        for (int i = 0; i < fDimuonVertexArr->GetEntries(); ++i) {
            TVertex* dimuonVert = (TVertex*) fDimuonVertexArr->At(i);
            if (
                    ((dimuonVert->index1 == leptonOneIndex && dimuonVert->index2 == leptonTwoIndex) ||
                    (dimuonVert->index1 == leptonTwoIndex && dimuonVert->index2 == leptonOneIndex)) &&
                    (dimuonVert->isValid)
               ) {
                    TVector3 vertPos;
                    TVector3 vertErr;  
                    vertPos.SetXYZ(dimuonVert->x, dimuonVert->y, dimuonVert->z);
                    vertErr.SetXYZ(dimuonVert->xerr, dimuonVert->yerr, dimuonVert->zerr);
                    rJpsiCand = vertPos; 
                    rJpsiCandErr = vertErr;
                    rJpsiCandChi2 = dimuonVert->chi2;
                    rJpsiCandNdof = dimuonVert->ndof;
                    hasValidJpsiVertex = true; 
                    break;
            }
        }
        
        if (!hasValidJpsiVertex)
            return kTRUE;
 
        hTotalEvents->Fill(9);
        
        if (photons.size() > 0)
            photonOneP4 = photons[0];
    }
    
    else if (params->selection == "2e2mu") {
        if ((muons_P4.size() < 2) || (electrons.size() < 2))
            return kTRUE;
        hTotalEvents->Fill(5);
        
        // Select the leading electron
        TLorentzVector electronOneP4;
        unsigned electronOneIndex = 0;
        bool validElectron1 = false;
        for (unsigned i = 0; i < (electrons.size() - 1); ++i) {
            if ( electrons[i].Pt() > 25.0 )
            {
                electronOneP4 = electrons[i];
                electronOneIndex = i;
                validElectron1 = true;
                break;
            }
        }

        if (!validElectron1)
            return kTRUE;

        hTotalEvents->Fill(6);

        // Select the oppositely charged second muon
        TLorentzVector electronTwoP4;
        unsigned electronTwoIndex = electronOneIndex + 1;
        bool validElectron2 = false;
        for (unsigned i = electronOneIndex+1; i < electrons.size(); ++i) {
            if (
                    electrons_q[electronOneIndex] != electrons_q[i]
                    && electrons[i].Pt() > 20.0
               )
            {
                electronTwoP4 = electrons[i];
                electronTwoIndex = i;
                validElectron2 = true;
                break;
            }
        }

        if (!validElectron2)
            return kTRUE;

        hTotalEvents->Fill(7);

        leptonOneP4      = electronOneP4;
        //leptonOneIso     = electrons_iso[electronOneIndex];
        leptonOneQ       = electrons_q[electronOneIndex];
        leptonOneTrigger = electrons_trigger[electronOneIndex];
        leptonOneFlavor  = 1;

        leptonOneD0        = electrons_d0[electronOneIndex];
        leptonOneDZ      = electrons_dz[electronOneIndex];
        leptonTwoP4      = electronTwoP4;
        //leptonTwoIso     = electrons_iso[electronTwoIndex];
        leptonTwoQ       = electrons_q[electronTwoIndex];
        leptonTwoTrigger = electrons_trigger[electronTwoIndex];
        leptonTwoFlavor  = 1;
            
        leptonTwoD0        = electrons_d0[electronTwoIndex];
        leptonTwoDZ        = electrons_dz[electronTwoIndex];
    
        // Check for dielectron vertex 
        unsigned int leptonOneIndex = electrons_index[electronOneIndex];
        unsigned int leptonTwoIndex = electrons_index[electronTwoIndex];

        bool hasValidZVertex = false;
         
        for (int i = 0; i < fDielectronVertexArr->GetEntries(); ++i) {
            TVertex* dielectronVert = (TVertex*) fDielectronVertexArr->At(i);
            if (
                    ((dielectronVert->index1 == leptonOneIndex && dielectronVert->index2 == leptonTwoIndex) ||
                    (dielectronVert->index1 == leptonTwoIndex && dielectronVert->index2 == leptonOneIndex)) &&
                    (dielectronVert->isValid)
               ) {
                    TVector3 vertPos;
                    TVector3 vertErr;  
                    vertPos.SetXYZ(dielectronVert->x, dielectronVert->y, dielectronVert->z);
                    vertErr.SetXYZ(dielectronVert->xerr, dielectronVert->yerr, dielectronVert->zerr);
                    rZCand = vertPos; 
                    rZCandErr = vertErr;
                    rZCandChi2 = dielectronVert->chi2;
                    rZCandNdof = dielectronVert->ndof;
                    hasValidZVertex = true; 
                    break;
            }
        }
        
        if (!hasValidZVertex)
            return kTRUE; 
        
        hTotalEvents->Fill(8);

        // Select the J/Psi candidate leading muon
        TLorentzVector muonThreeP4;
        unsigned muonThreeIndex = 0;
        bool validMuon3 = false;
        for (unsigned i = 0; i < (muons_P4.size() - 1); ++i) {
            if (
                    (muons_iso[i]/muons_P4[i].Pt()) < 0.25
                    && muons_P4[i].Pt() > 4.0
                )
            {
                muonThreeP4 = muons_P4[i];
                muonThreeIndex = i;
                validMuon3 = true;
                break;
            }
        }
        
        if (!validMuon3)
            return kTRUE;
       
        hTotalEvents->Fill(9);

        // Select the J/Psi candidate subleading muon
        TLorentzVector muonFourP4;
        unsigned muonFourIndex = muonThreeIndex + 1;
        bool validMuon4 = false;
        for (unsigned i = muonThreeIndex+1; i < muons_P4.size(); ++i) {
            if (
                    muons_q[muonThreeIndex] != muons_q[i]
               )
            {
                muonFourP4 = muons_P4[i];
                muonFourIndex = i;
                validMuon4 = true;
                break;
            }
        }
        
        if (!validMuon4)
            return kTRUE;

        hTotalEvents->Fill(10);
        
        leptonThreeP4      = muonThreeP4;
        leptonThreeIso     = muons_iso[muonThreeIndex];
        leptonThreeQ       = muons_q[muonThreeIndex];
        leptonThreeTrigger = muons_trigger[muonThreeIndex];
        leptonThreeIsTight = muons_isTight[muonThreeIndex];
        leptonThreeFlavor  = 1;

        leptonThreeD0        = muons_d0[muonThreeIndex];
        leptonThreeDZ        = muons_dz[muonThreeIndex];

        leptonFourP4      = muonFourP4;
        leptonFourIso     = muons_iso[muonFourIndex];
        leptonFourQ       = muons_q[muonFourIndex];
        leptonFourTrigger = muons_trigger[muonFourIndex];
        leptonFourIsTight = muons_isTight[muonFourIndex];
        leptonFourFlavor  = 1;
            
        leptonFourD0        = muons_d0[muonFourIndex];
        leptonFourDZ        = muons_dz[muonFourIndex];
       
        // Check for dimuon vertex 
        unsigned int leptonThreeIndex = muons_index[muonThreeIndex];
        unsigned int leptonFourIndex = muons_index[muonFourIndex];

        bool hasValidJpsiVertex = false;
         
        for (int i = 0; i < fDimuonVertexArr->GetEntries(); ++i) {
            TVertex* dimuonVert = (TVertex*) fDimuonVertexArr->At(i);
            if (
                    ((dimuonVert->index1 == leptonThreeIndex && dimuonVert->index2 == leptonFourIndex) ||
                    (dimuonVert->index1 == leptonFourIndex && dimuonVert->index2 == leptonThreeIndex)) &&
                    (dimuonVert->isValid)
               ) {
                    TVector3 vertPos;
                    TVector3 vertErr;  
                    vertPos.SetXYZ(dimuonVert->x, dimuonVert->y, dimuonVert->z);
                    vertErr.SetXYZ(dimuonVert->xerr, dimuonVert->yerr, dimuonVert->zerr);
                    rJpsiCand = vertPos; 
                    rJpsiCandErr = vertErr;
                    rJpsiCandChi2 = dimuonVert->chi2;
                    rJpsiCandNdof = dimuonVert->ndof;
                    hasValidJpsiVertex = true; 
                    break;
            }
        }
        
        if (!hasValidJpsiVertex)
            return kTRUE;
 
        hTotalEvents->Fill(11);
        
        if (photons.size() > 0)
            photonOneP4 = photons[0];
    }

    else if (params->selection == "ee") {

        if (electrons.size() != 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        TLorentzVector dielectron;
        dielectron = electrons[0] + electrons[1];
        if (dielectron.M() < 12. || dielectron.M() > 70.)
            return kTRUE;
        hTotalEvents->Fill(6);

        leptonOneP4      = electrons[0];
        //leptonOneIso     = electrons_iso[0];
        leptonOneQ       = electrons_q[0];
        leptonOneTrigger = electrons_trigger[0];
        leptonOneFlavor  = 0;

        leptonTwoP4      = electrons[1];
        //leptonTwoIso     = electrons_iso[1];
        leptonTwoQ       = electrons_q[1];
        leptonTwoTrigger = electrons_trigger[1];
        leptonTwoFlavor  = 0;
    }

    else if (params->selection == "4l") {

        /*if (nMuons >= 4) { 

            if (muons_P4_type1[0].Pt() < 25)
                return kTRUE;
            hTotalEvents->Fill(6);

            leptonOneP4        = muons[0];
            leptonOneIso       = muons_iso[0];
            leptonOneQ         = muons_q[0];
            leptonOneTrigger   = muons_trigger[0];
            leptonOneFlavor    = 0;

            leptonOneD0        = muons_d0[0];
            leptonOneDZ        = muons_dz[0];

            leptonTwoP4        = muons[1];
            leptonTwoIso       = muons_iso[1];
            leptonTwoQ         = muons_q[1];
            leptonTwoTrigger   = muons_trigger[1];
            leptonTwoFlavor    = 0;
            
            leptonTwoD0        = muons_d0[1];
            leptonTwoDZ        = muons_dz[1];
            
            leptonThreeP4      = muons[2];
            leptonThreeIso     = muons_iso[2];
            leptonThreeQ       = muons_q[2];
            leptonThreeTrigger = muons_trigger[2];
            leptonThreeFlavor  = 0;
            leptonFourP4       = muons[3];
            leptonFourIso      = muons_iso[3];
            leptonFourQ        = muons_q[3];
            leptonFourTrigger  = muons_trigger[3];
            leptonFourFlavor   = 0;

        } else if (muons.size() >= 2 and electrons.size() >= 2) {

            if (muons[0].Pt() < 25)
                return kTRUE;
            hTotalEvents->Fill(6);
        
            leptonOneP4      = muons[0];
            leptonOneIso     = muons_iso[0];
            leptonOneQ       = muons_q[0];
            leptonOneTrigger = muons_trigger[0];
            leptonOneFlavor  = 0;
            
            leptonOneD0        = muons_d0[0];
            leptonOneDZ        = muons_dz[0];
           
            leptonTwoP4      = muons[1];
            leptonTwoIso     = muons_iso[1];
            leptonTwoQ       = muons_q[1];
            leptonTwoTrigger = muons_trigger[1];
            leptonTwoFlavor  = 0;
            
            leptonTwoD0        = muons_d0[1];
            leptonTwoDZ        = muons_dz[1];

            leptonThreeP4      = electrons[0];
            leptonThreeIso     = electrons_iso[0];
            leptonThreeQ       = electrons_q[0];
            leptonThreeTrigger = electrons_trigger[0];
            leptonThreeFlavor  = 1;
            leptonFourP4       = electrons[1];
            leptonFourIso      = electrons_iso[1];
            leptonFourQ        = electrons_q[1];
            leptonFourTrigger  = electrons_trigger[1];
            leptonFourFlavor   = 1;
        } else {
            hTotalEvents->Fill(5);
            return kTRUE;
        }*/
    }


    ///////////////////
    // Fill jet info //
    ///////////////////

    if (bjets.size() > 0) {
        bjetP4.SetPtEtaPhiM(bjets[0]->pt, bjets[0]->eta, bjets[0]->phi, bjets[0]->mass);
        bjetD0     = bjets[0]->d0;
        bjetTag    = bjets[0]->csv;
        bjetPUID   = bjets[0]->mva;
        bjetFlavor = bjets[0]->hadronFlavor;
    } else {
        bjetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        bjetD0     = 0.;
        bjetTag    = 0.;
        bjetPUID   = 0.;
        bjetFlavor = 0.;
    }

    if (fwdjets.size() > 0) {
        jetP4.SetPtEtaPhiM(fwdjets[0]->pt, fwdjets[0]->eta, fwdjets[0]->phi, fwdjets[0]->mass);
        jetD0     = fwdjets[0]->d0;
        jetTag    = 0.;
        jetPUID   = fwdjets[0]->mva;
        jetFlavor = fwdjets[0]->hadronFlavor;
    } else if (jets.size() > 0) {
        jetP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
        jetD0     = jets[0]->d0;
        jetTag    = jets[0]->csv;
        jetPUID   = jets[0]->mva;
        jetFlavor = jets[0]->hadronFlavor;
    } else {
        jetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        jetD0     = 0.;
        jetTag    = 0.;
        jetPUID   = 0.;
        jetFlavor = 0.;
    } 

    if (genbjets.size() > 0 && !isData) {
        genBJetP4.SetPtEtaPhiM(genbjets[0]->genpt, genbjets[0]->geneta, genbjets[0]->genphi, genbjets[0]->genm);
        genBJetTag = genbjets[0]->csv;
    } else {
        genBJetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        genBJetTag = 0;
    }

    if (genjets.size() > 0 && !isData) {
        genJetP4.SetPtEtaPhiM(genjets[0]->genpt, genjets[0]->geneta, genjets[0]->genphi, genjets[0]->genm);
        genJetTag = genjets[0]->csv;
    } else {
        genJetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        genJetTag = 0;
    }

    outTree->Fill();
    this->passedEvents++;
    return kTRUE;
}

void zjpsiAnalyzerV2::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void zjpsiAnalyzerV2::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void zjpsiAnalyzerV2::ReportPostTerminate()
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
    std::unique_ptr<zjpsiAnalyzerV2> selector(new zjpsiAnalyzerV2());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

float zjpsiAnalyzerV2::MetKluge(float met)
{
    if (met > 500) {
        return 1.;
    }

    float bins[]        = {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 250., 500.};
    float corrections[] = {1., 0.99, 0.96, 0.95, 0.948, 0.947, 0.95, 0.957, 0.966, 0.97, 0.971, 0.98};

    int bin = 0;
    for (int i = 0; i < 12; ++i) {
        if (met > bins[i] && met <= bins[i+1]) {
            bin = i;
            break;
        }
    }
    
    return corrections[bin];
}
