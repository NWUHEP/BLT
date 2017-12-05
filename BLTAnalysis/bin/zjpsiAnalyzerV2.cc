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

    } else if (params->selection == "2e2mu") {
        triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
        triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*");

    } else if (params->selection == "ee") {
        triggerNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
    
    } else if (params->selection == "jpsi_control") {
        triggerNames.push_back("HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing_v*");
    }

    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
    lumiMask.AddJSONFile(jsonFileName);

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
    outTree->Branch("genWeight", &genWeight);
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
    outTree->Branch("rZValid", &rZValid);
    
    outTree->Branch("rJpsiCand", &rJpsiCand);
    outTree->Branch("rJpsiCandErr", &rJpsiCandErr);
    outTree->Branch("rJpsiCandChi2", &rJpsiCandChi2);
    outTree->Branch("rJpsiCandNdof", &rJpsiCandNdof);
    outTree->Branch("rJpsiCandProb", &rJpsiCandProb);
    outTree->Branch("rJpsiCandRxy", &rJpsiCandRxy);
    outTree->Branch("rJpsiCandRxyErr", &rJpsiCandRxyErr);
    outTree->Branch("rJpsiValid", &rJpsiValid);

    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);
    outTree->Branch("ht", &ht);
    outTree->Branch("htPhi", &htPhi);

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
    outTree->Branch("leptonOneMother", &leptonOneMother);

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
    outTree->Branch("leptonTwoMother", &leptonTwoMother);
    
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
    outTree->Branch("leptonThreeMother", &leptonThreeMother);

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
    outTree->Branch("leptonFourMother", &leptonFourMother);

    // jets
    //outTree->Branch("jetP4", &jetP4);
    //outTree->Branch("jetD0", &jetD0);
    //outTree->Branch("jetTag", &jetTag);
    //outTree->Branch("jetPUID", &jetPUID);
    //outTree->Branch("jetFlavor", &jetFlavor);

    //outTree->Branch("bjetP4", &bjetP4);
    //outTree->Branch("bjetD0", &bjetD0);
    //outTree->Branch("bjetTag", &bjetTag);
    //outTree->Branch("bjetPUID", &bjetPUID);
    //outTree->Branch("bjetFlavor", &bjetFlavor);

    //outTree->Branch("genBJetP4", &genBJetP4);
    //outTree->Branch("genBJetTag", &genBJetTag);
    //outTree->Branch("genJetP4", &genJetP4);
    //outTree->Branch("genJetTag", &genJetTag);
    
    outTree->Branch("jetOneP4", &jetOneP4);
    outTree->Branch("jetOneTag", &jetOneTag);

    outTree->Branch("jetTwoP4", &jetTwoP4);
    outTree->Branch("jetTwoTag", &jetTwoTag);

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
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",20,0.5,20.5);

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

    if (!passTrigger)
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
    
    //cout << "New event" << endl;
    vector<TGenParticle*> genParticles;
    std::bitset<4> genEventType; 
    if (!isData) {
        unsigned count = 0;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
            //cout << i << ", " << particle->status << ", " << particle->pdgId  << ", " << particle->parent;
            //cout << "\t" << particle->pt << ", " << particle->eta;
            //cout << endl;
            genParticles.push_back(particle);

            if (
                    particle->status == 23 
                    && (abs(particle->pdgId) < 6 || particle->pdgId == 21) 
                    && particle->parent != -2
               ) {
                ++count;
            }

            // This will save Z and J/psi in 4l events
            //if (
            //        params->selection == "4mu"
            //        && (abs(particle->pdgId) == 443 || abs(particle->pdgId) == 23) 
            //        //&& particle->status == 62 
            //   ) {
            //    genParticles.push_back(particle);
            //}

            //// Let's save some quarks too
            //if (
            //        params->selection == "4mu"
            //        && (0 < abs(particle->pdgId) && abs(particle->pdgId) < 7)
            //   ) { 
            //    genParticles.push_back(particle);
            //}

            // This will save leptons (mostly for use with dilepton analyses)
            //if (abs(particle->pdgId) == 11 || abs(particle->pdgId) == 13) {
            //    if (particle->parent != -2) {
            //        TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
            //        //if (fabs(mother->pdgId) == 24 || abs(mother->pdgId) == 15) {
            //            //cout << i << ", " << particle->status << ", " << particle->pdgId  << ", " << particle->parent << ", " << mother->pdgId;
            //            //cout << "\t" << particle->pt << ", " << particle->eta;
            //            //cout << endl;
            //            genParticles.push_back(particle);

            //            //// categorize the event (only works for e+mu events)
            //            //if (abs(particle->pdgId) == 11 && abs(mother->pdgId) == 24) {
            //            //    genEventType.set(0);
            //            //} else  if (abs(particle->pdgId) == 11 && abs(mother->pdgId) == 15) {
            //            //    genEventType.set(1);
            //            //} else  if (abs(particle->pdgId) == 13 && abs(mother->pdgId) == 24) {
            //            //    genEventType.set(2);
            //            //} else  if (abs(particle->pdgId) == 13 && abs(mother->pdgId) == 15) {
            //            //    genEventType.set(3);
            //            //} 
            //        //}
            //    }
            //}

            nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY
        }

        //if (genEventType.to_string() == "0101") { // W->e + W->mu
        //    hTotalEvents->Fill(11);
        //} else if (genEventType.to_string() == "0110") { // tau->e + W->mu
        //    hTotalEvents->Fill(12);
        //} else if (genEventType.to_string() == "1001") { // W->e + tau->mu
        //    hTotalEvents->Fill(13);
        //} else if (genEventType.to_string() == "1010") { // tau->e + tau->mu
        //    hTotalEvents->Fill(14);
        //} else { // everything else
        //    hTotalEvents->Fill(15);
        //}
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
    vector<TLorentzVector> veto_muons;

    for (int i=0; i < fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);
        
    /* Apply Rochester Muon Corrections */
        TLorentzVector muonP4;
        copy_p4(muon, MUON_MASS, muonP4);
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
        muonP4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);
        
        // Remove muons with very small deltaR
        float minDeltaR = 1e6;
        for (unsigned j=0; j < muons.size(); ++j) {
            TLorentzVector tmpMuonP4;
            tmpMuonP4.SetPtEtaPhiM(muons[j]->pt, muons[j]->eta, muons[j]->phi, 0.1051);
            float dr = muonP4.DeltaR(tmpMuonP4);
            if (dr < minDeltaR) {
                minDeltaR = dr;
            }
        }

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
        
        // muons for jet veto
        if (
                muonP4.Pt() > 10
                // tight muon ID and ISO
                && (muon->typeBits & baconhep::kPOGTightMuon)
                && pfiso/muonP4.Pt() < 0.15
           ) {
            veto_muons.push_back(muonP4);
        }
    }
    std::sort(muons_P4.begin(), muons_P4.end(), P4SortCondition);
   

    /* ELECTRONS */
    vector<TElectron*> electrons;
    vector<TLorentzVector> veto_electrons;
    vector<bool> electrons_trigger;
    vector<unsigned int> electrons_index;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        TLorentzVector electronP4;
        electronP4.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, ELE_MASS);

        if (
                electron->pt > 10
                && fabs(electron->eta) < 2.5
                && particleSelector->PassElectronID(electron, cuts->tightElID)
                && particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
                && GetElectronIsolation(electron, fInfo->rhoJet)/electronP4.Pt() < 0.35
                && fabs(electron->d0) < 0.5
                && fabs(electron->dz) < 1.0
                //&& electron->sip3d < 4.0 
           ) {
            electrons.push_back(electron);
            electrons_index.push_back(electron->eleIndex);
            veto_electrons.push_back(electronP4);

            // trigger matching
            bool triggered = false;
            for (unsigned i = 0; i < triggerNames.size(); ++i) {
                triggered |= trigger->passObj(triggerNames[i], 1, electron->hltMatchBits);
            }
            electrons_trigger.push_back(triggered);
        }
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);

    
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
    TLorentzVector hadronicP4;

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
        TLorentzVector jetP4; 
        jetP4.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
        bool muOverlap = false;
        for (const auto& mu: veto_muons) {
            if (jetP4.DeltaR(mu) < 0.4) {
                muOverlap = true;
                break;
            }
        }
        bool elOverlap = false;
        for (const auto& el: veto_electrons) {
            if (jetP4.DeltaR(el) < 0.4) {
                elOverlap = true;
                break;
            }
        }

        if (
                jet->pt > 30 
                && fabs(jet->eta) < 4.7
                && particleSelector->PassJetID(jet, cuts->looseJetID)
                && !muOverlap 
                && !elOverlap
           ) {

            if (fabs(jet->eta) <= 2.4) { 
                hadronicP4 += jetP4;
                jets.push_back(jet);

                if (isData) {
                    if (jet->bmva > 0.9432) { 
                        ++nBJets;
                    } else {
                        ++nJets;
                    }
                } else {
                    if (jet->bmva > 0.9432) {//particleSelector->BTagModifier(jet, "MVAT")) { 
                        ++nBJets;
                    } else {
                        ++nJets;
                    }
                }
            } else {
                if (fabs(jet->eta) > 2.5) {
                    hadronicP4 += jetP4;
                    ++nFwdJets;
                }
            }
        }
    }
    std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);


    /* MET */
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;
    
    /* HT */
    ht    = hadronicP4.Pt();
    htPhi = hadronicP4.Phi();


    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    nMuons     = muons.size();
    nElectrons = electrons.size();
    nPhotons = photons.size();

    if (params->selection == "mumu") {
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

        //float thisZMass = (muonOneP4 + muonTwoP4).M();
        //if (
        //        thisZMass < 2.0 ||
        //        (thisZMass > 4.0 && thisZMass < 75.) ||
        //        thisZMass > 107.
        //   )
        //    return kTRUE;
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
        
        //if (!hasValidZVertex)
        //    return kTRUE; 
        //
        //
        rZValid = hasValidZVertex;
        hTotalEvents->Fill(9);

        if (photons.size() > 0)
            photonOneP4 = photons[0];

        if (!isData) {
            eventWeight *= weights->GetMuonIDEff(muonOneP4);
            eventWeight *= weights->GetMuonISOEff(muonOneP4);
            eventWeight *= weights->GetMuonIDEff(muonTwoP4);
            eventWeight *= weights->GetMuonISOEff(muonTwoP4);
            leptonOneMother = GetGenMotherId(genParticles, leptonOneP4);
            leptonTwoMother = GetGenMotherId(genParticles, leptonTwoP4);

            // trigger weight
            //pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
            //pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
            //eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);
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

        //float thisZMass = (muonOneP4 + muonTwoP4).M();
        //if (
        //        thisZMass < 2.0 ||
        //        (thisZMass > 4.0 && thisZMass < 75.) ||
        //        thisZMass > 107.
        //   )
        //    return kTRUE;
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
        
        //if (!hasValidZVertex)
        //    return kTRUE; 
        //
        //
        rZValid = hasValidZVertex;
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
        
        //float thisJpsiMass = (muonThreeP4 + muonFourP4).M();
        //if (
        //        thisJpsiMass < 2.0 ||
        //        (thisJpsiMass > 4.0 && thisJpsiMass < 75.) ||
        //        thisJpsiMass > 107.
        //   )
        //    return kTRUE;
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
        
        //if (!hasValidJpsiVertex)
        //    return kTRUE;
 
        rJpsiValid = hasValidJpsiVertex;
        hTotalEvents->Fill(13);

        //if (!hasGoodJpsiVertex)
        //    return kTRUE;

        hTotalEvents->Fill(14);

        //if (fabs(rJpsiCand.Z() - rZCand.Z()) > 1.)
        //    return kTRUE;

        //hTotalEvents->Fill(16);
        
        if (photons.size() > 0)
            photonOneP4 = photons[0];

        if (!isData) { 
            leptonOneMother = GetGenMotherId(genParticles, leptonOneP4);
            leptonTwoMother = GetGenMotherId(genParticles, leptonTwoP4);
            leptonThreeMother = GetGenMotherId(genParticles, leptonThreeP4);
            leptonFourMother = GetGenMotherId(genParticles, leptonFourP4);

            cout << "leptonOneMother = " << leptonOneMother << endl;
            cout << "leptonTwoMother = " << leptonTwoMother << endl;
            cout << "leptonThreeMother = " << leptonThreeMother << endl;
            cout << "leptonFourMother = " << leptonFourMother << endl;
        
            eventWeight *= weights->GetMuonIDEff(leptonOneP4);
            eventWeight *= weights->GetMuonISOEff(leptonOneP4);
            eventWeight *= weights->GetMuonIDEff(leptonTwoP4);
            eventWeight *= weights->GetMuonISOEff(leptonTwoP4);
            //eventWeight *= weights->GetMuonIDEff(leptonThreeP4);
            //eventWeight *= weights->GetMuonISOEff(muonP4[2]);
            //eventWeight *= weights->GetMuonIDEff(leptonFourP4);
            //eventWeight *= weights->GetMuonISOEff(muonP4[3]);

            // trigger weight
            //pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
            //pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
            //eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);

        }

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
        
        TLorentzVector electronOneP4;
        unsigned electronOneIndex = 0;
        bool validElectron1 = false;
        TLorentzVector electronTwoP4;
        unsigned electronTwoIndex = electronOneIndex + 1;
        bool validElectron2 = false;

        unsigned int pairing_algo = 0;

        if (pairing_algo == 0) {
        
            // Select the leading electron
            for (unsigned i = 0; i < (electrons.size() - 1); ++i) {
                TLorentzVector tempElectronOne;
                tempElectronOne.SetPtEtaPhiM(electrons[i]->pt, electrons[i]->eta, electrons[i]->phi, ELE_MASS);
                if (tempElectronOne.Pt() > 25.0)
                {
                    electronOneP4 = tempElectronOne;
                    electronOneIndex = i;
                    validElectron1 = true;
                    break;
                }
            }

            // Select the oppositely charged second electron
            for (unsigned i = electronOneIndex+1; i < electrons.size(); ++i) {
                TLorentzVector tempElectronTwo;
                tempElectronTwo.SetPtEtaPhiM(electrons[i]->pt, electrons[i]->eta, electrons[i]->phi, ELE_MASS);
                if (
                        electrons[electronOneIndex]->q != electrons[i]->q
                        && tempElectronTwo.Pt() > 15.0
                   )
                {
                    electronTwoP4 = tempElectronTwo;
                    electronTwoIndex = i;
                    validElectron2 = true;
                    break;
                }
            }
        }

        else if (pairing_algo == 1) {
            float zMass_pdg = 91.188;
            float zMassDiff = 100.;
            for (unsigned int i = 0; i < electrons.size(); ++i) {
                TLorentzVector tempElectronOne;
                tempElectronOne.SetPtEtaPhiM(electrons[i]->pt, electrons[i]->eta, electrons[i]->phi, ELE_MASS);
                for (unsigned int j = i + 1; j < electrons.size(); ++j) {
                    TLorentzVector tempElectronTwo;
                    tempElectronTwo.SetPtEtaPhiM(electrons[j]->pt, electrons[j]->eta, electrons[j]->phi, ELE_MASS);
                    float thisMass = (tempElectronOne + tempElectronTwo).M();
                    if (
                            fabs(thisMass - zMass_pdg) < zMassDiff
                            && tempElectronOne.Pt() > 25.0
                            && tempElectronTwo.Pt() > 15.0
                            && electrons[i]->q != electrons[j]->q
                       ) 
                    {
                        electronOneIndex = i;
                        electronTwoIndex = j;
                        electronOneP4 = tempElectronOne;
                        electronTwoP4 = tempElectronTwo;
                        validElectron1 = true;
                        validElectron2 = true;
                        zMassDiff = fabs(thisMass - zMass_pdg);
                    }
                }
            }
        }

        if (!validElectron1)
            return kTRUE;

        hTotalEvents->Fill(6);

        if (!validElectron2)
            return kTRUE;

        hTotalEvents->Fill(7);

        leptonOneP4      = electronOneP4;
        leptonOneIso    = GetElectronIsolation(electrons[electronOneIndex], fInfo->rhoJet);
        leptonOneFlavor = electrons[electronOneIndex]->q*11;
        leptonOneDZ     = electrons[electronOneIndex]->dz;
        leptonOneD0     = electrons[electronOneIndex]->d0;
        leptonOneTrigger = electrons_trigger[electronOneIndex];

        leptonTwoP4      = electronTwoP4;            
        leptonTwoIso    = GetElectronIsolation(electrons[electronTwoIndex], fInfo->rhoJet);
        leptonTwoFlavor = electrons[electronTwoIndex]->q*11;
        leptonTwoDZ     = electrons[electronTwoIndex]->dz;
        leptonTwoD0     = electrons[electronTwoIndex]->d0;
        leptonTwoTrigger = electrons_trigger[electronTwoIndex];

        //float thisZMass = (electronOneP4 + electronTwoP4).M();
        //if (
        //        thisZMass < 2.0 ||
        //        (thisZMass > 4.0 && thisZMass < 75.) ||
        //        thisZMass > 107.
        //   )
        //    return kTRUE;
        hTotalEvents->Fill(8);

        // Check for dielectron vertex 
        /*unsigned int leptonOneIndex = electrons_index[electronOneIndex];
        unsigned int leptonTwoIndex = electrons_index[electronTwoIndex];

        bool hasValidZVertex = false;
        //bool hasGoodZVertex = false;
         
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
                    rZCandProb = dielectronVert->prob;
                    rZCandRxy = dielectronVert->rxy;
                    rZCandRxyErr = dielectronVert->rxy_err;
                    hasValidZVertex = true; 
                    //if (rZCandProb > 0.005)
                    //    hasGoodZVertex = true;
                    break;
            }
        }
 
        //if (!hasValidZVertex)
        //    return kTRUE; 
       
        rZValid = hasValidZVertex; */
        hTotalEvents->Fill(9);

        //if (!hasGoodZVertex)
        //    return kTRUE;

        //hTotalEvents->Fill(10);
        //

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
                        (muons_iso[i]/muons_P4[i].Pt()) < 0.25
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
                if (muons_q[muonThreeIndex] != muons_q[i])
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
        
        //float thisJpsiMass = (muonThreeP4 + muonFourP4).M();
        //if (
        //        thisJpsiMass < 2.0 ||
        //        (thisJpsiMass > 4.0 && thisJpsiMass < 75.) ||
        //        thisJpsiMass > 107.
        //   )
        //    return kTRUE;
        hTotalEvents->Fill(12);
       
        // Check for dimuon vertex 
        /* unsigned int leptonThreeIndex = muons_index[muonThreeIndex];
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
        
        //if (!hasValidJpsiVertex)
        //    return kTRUE;
    
        rJpsiValid = hasValidJpsiVertex; */
        hTotalEvents->Fill(13);

        //if (!hasGoodJpsiVertex)
        //    return kTRUE;

        hTotalEvents->Fill(14);

        //if (fabs(rJpsiCand.Z() - rZCand.Z()) > 1.)
        //    return kTRUE;

        //hTotalEvents->Fill(16);
        
        if (photons.size() > 0)
            photonOneP4 = photons[0];
        
        if (!isData) {
        leptonOneMother = GetGenMotherId(genParticles, leptonOneP4);
        leptonTwoMother = GetGenMotherId(genParticles, leptonTwoP4);
        leptonThreeMother = GetGenMotherId(genParticles, leptonThreeP4);
        leptonFourMother = GetGenMotherId(genParticles, leptonFourP4);

        eventWeight *= weights->GetElectronRecoIdEff(leptonOneP4);
        eventWeight *= weights->GetElectronRecoIdEff(leptonTwoP4);
        genWeight = fGenEvtInfo->weight;
        //eventWeight *= weights->GetMuonIDEff(leptonThreeP4);
        //eventWeight *= weights->GetMuonIDEff(leptonFourP4);
        }
        
    }
    
//    else if (params->selection == "ee") {
//
//        if (electrons.size() != 2)
//            return kTRUE;
//        hTotalEvents->Fill(5);
//
//        TLorentzVector dielectron;
//        dielectron = electrons[0] + electrons[1];
//        if (dielectron.M() < 12. || dielectron.M() > 70.)
//            return kTRUE;
//        hTotalEvents->Fill(6);
//
//        leptonOneP4      = electrons[0];
//        //leptonOneIso     = electrons_iso[0];
//        leptonOneQ       = electrons_q[0];
//        leptonOneTrigger = electrons_trigger[0];
//        leptonOneFlavor  = 0;
//
//        leptonTwoP4      = electrons[1];
//        //leptonTwoIso     = electrons_iso[1];
//        leptonTwoQ       = electrons_q[1];
//        leptonTwoTrigger = electrons_trigger[1];
//        leptonTwoFlavor  = 0;
//    }

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

    if (jets.size() > 0) {
        jetOneP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
        jetOneTag    = jets[0]->csv;
    } else {
        jetOneP4.SetPtEtaPhiM(0., 0., 0., 0.);
        jetOneTag    = 0.;
    }

    if (jets.size() > 1) {
        jetTwoP4.SetPtEtaPhiM(jets[1]->pt, jets[1]->eta, jets[1]->phi, jets[1]->mass);
        jetTwoTag    = jets[1]->csv;
    } else {
        jetTwoP4.SetPtEtaPhiM(0., 0., 0., 0.);
        jetTwoTag    = 0.;
    } 

    //if (!isData && genParticles.size() == 2) {
    //    genOneId = genParticles[0]->pdgId;
    //    genOneP4.SetPtEtaPhiM(genParticles[0]->pt, genParticles[0]->eta, genParticles[0]->phi, genParticles[0]->mass); 
    //    genTwoId = genParticles[1]->pdgId;
    //    genTwoP4.SetPtEtaPhiM(genParticles[1]->pt, genParticles[1]->eta, genParticles[1]->phi, genParticles[1]->mass); 
    //} else {
    //    genOneP4.SetPtEtaPhiM(0., 0., 0., 0.); 
    //    genTwoP4.SetPtEtaPhiM(0., 0., 0., 0.); 
    //}

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

int zjpsiAnalyzerV2::GetGenMotherId(vector<TGenParticle*> particles, TLorentzVector p4)
{
    int motherId = 0;
    for (unsigned i = 0; i < particles.size(); ++i) {
        TLorentzVector genP4;
        genP4.SetPtEtaPhiM(particles[i]->pt, particles[i]->eta, particles[i]->phi, particles[i]->mass); 
        if (genP4.DeltaR(p4) < 0.3) {
            TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particles[i]->parent);
            motherId = mother->pdgId;
            if (abs(motherId) == 13) {
                bool foundGranny = false;
                while (!foundGranny) {
                    if (mother->parent > 0) {
                        mother = (TGenParticle*) fGenParticleArr->At(mother->parent);
                        motherId = mother->pdgId;
                        if (abs(motherId) != 13)
                            foundGranny = true;
                    }
                    else {
                        motherId = 0;
                        break;
                    }
                }
            }
                
        }
    }
    return motherId;
}

float zjpsiAnalyzerV2::GetMuonIsolation(const baconhep::TMuon* mu)
{
    float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    return combIso;
}

float zjpsiAnalyzerV2::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
{
    int iEta = 0;
    float etaBins[8] = {0., 1., 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};
    float effArea[8] = {0.1703, 0.1715, 0.1213, 0.1230, 0.1635, 0.1937, 0.2393};
    for (unsigned i = 0; i < 8; ++i) {
        if (fabs(el->scEta) > etaBins[i] && fabs(el->scEta) < etaBins[i+1]) {
            iEta = i;
            break;
        }
    }

    float combIso = el->chHadIso + std::max(0., (double)el->neuHadIso + el->gammaIso - rho*effArea[iEta]);

    return combIso;
}
