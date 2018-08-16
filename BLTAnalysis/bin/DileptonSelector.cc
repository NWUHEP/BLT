#include "DileptonSelector.h"
#include <map>

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

const float ZMASS = 91.19;
const float PSIMASS = 3.096;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 
bool sort_by_btag(const baconhep::TJet* lhs, const baconhep::TJet* rhs) 
{
    return lhs->bmva > rhs->bmva;
}


DileptonSelector::DileptonSelector() : BLTSelector()
{

}

DileptonSelector::~DileptonSelector()
{

}

void DileptonSelector::Begin(TTree *tree)
{
    rng = new TRandom3();

    // Parse command line option
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1), std::sregex_token_iterator(), std::back_inserter(options));

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

    if (params->selection == "single_lepton") {
        triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");
        //triggerNames.push_back("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v*");

    } else if (params->selection == "single_muon") {
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");

    } else if (params->selection == "single_electron") {
        triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");

    } else if (params->selection == "double_muon") {
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");

    } else if (params->selection == "double_electron") {
        triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");

    } else if (params->selection == "mueg") {
        triggerNames.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
        triggerNames.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
        triggerNames.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*");
        triggerNames.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*");
    }

    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
    lumiMask.AddJSONFile(jsonFileName);

    // muon momentum corrections
    muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/rcdata.2016.v3");

    // electron scale corrections
    electronScaler = new EnergyScaleCorrection(cmssw_base + "/src/BLT/BLTAnalysis/data");

    // Prepare the output tree
    string outFileName = params->get_output_filename("output");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();

    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);
    outHistName = params->get_output_treename("GenCategory");
    hGenCat = new TH1D(outHistName.c_str(), "WW decay modes",30,0.5,30.5);

    vector<std::string> channelNames = {"mumu", "ee"};
    for (unsigned i = 0; i < channelNames.size(); ++i) {
        string channel = channelNames[i];
        outFile->mkdir(channel.c_str());
        outFile->cd(channel.c_str());
        string treeName = "bltTree_" + params->datasetgroup;
        tree = new TTree(treeName.c_str(), treeName.c_str());

        tree->Branch("runNumber", &runNumber);
        tree->Branch("evtNumber", &evtNumber, "eventNumber/l");
        tree->Branch("lumiSection", &lumiSection);
        tree->Branch("triggerStatus", &triggerStatus);
        tree->Branch("nPV", &nPV);
        tree->Branch("nPU", &nPU);
        tree->Branch("nPartons", &nPartons);
        tree->Branch("rPV", &rPV);

        // weights and their uncertainties
        tree->Branch("eventWeight", &eventWeight);
        tree->Branch("leptonOneRecoWeight", &leptonOneRecoWeight);
        tree->Branch("leptonTwoRecoWeight", &leptonTwoRecoWeight);
        tree->Branch("topPtWeight", &topPtWeight);
        tree->Branch("puWeight", &puWeight);
        tree->Branch("triggerWeight", &triggerWeight);
        tree->Branch("genWeight", &genWeight);

        // met and ht
        tree->Branch("met", &met);
        tree->Branch("metPhi", &metPhi);
        tree->Branch("ht", &ht);
        tree->Branch("htPhi", &htPhi);
        tree->Branch("htSum", &htSum);

        // leptons
        tree->Branch("leptonOneP4", &leptonOneP4);
        tree->Branch("leptonOneIso", &leptonOneIso);
        tree->Branch("leptonOneFlavor", &leptonOneFlavor);
        tree->Branch("leptonOneMother", &leptonOneMother);
        tree->Branch("leptonOneD0", &leptonOneD0);
        tree->Branch("leptonOneDZ", &leptonOneDZ);

        tree->Branch("leptonTwoP4", &leptonTwoP4);
        tree->Branch("leptonTwoIso", &leptonTwoIso);
        tree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
        tree->Branch("leptonTwoMother", &leptonTwoMother);
        tree->Branch("leptonTwoD0", &leptonTwoD0);
        tree->Branch("leptonTwoDZ", &leptonTwoDZ);

        // jets
        tree->Branch("jetOneP4", &jetOneP4);
        tree->Branch("jetOneTag", &jetOneTag);
        tree->Branch("jetOneFlavor", &jetOneFlavor);
        tree->Branch("jetTwoP4", &jetTwoP4);
        tree->Branch("jetTwoTag", &jetTwoTag);
        tree->Branch("jetTwoFlavor", &jetTwoFlavor);

        // gen level objects
        tree->Branch("genCategory", &genCategory);
        tree->Branch("genOneP4", &genOneP4);
        tree->Branch("genOneId", &genOneId);
        //tree->Branch("genOneMother", &genOneMother);
        tree->Branch("genTwoP4", &genTwoP4);
        tree->Branch("genTwoId", &genTwoId);
        //tree->Branch("genTwoMother", &genTwoMother);

        // object counters
        tree->Branch("nMuons", &nMuons);
        tree->Branch("nElectrons", &nElectrons);
        tree->Branch("nTaus", &nTaus);
        tree->Branch("nPhotons", &nPhotons);
        tree->Branch("nJets", &nJets);
        tree->Branch("nFwdJets", &nFwdJets);
        tree->Branch("nBJets", &nBJets);

        outTrees[channel] = tree;

        // event counter
        string outHistName = params->get_output_treename("TotalEvents_" + channel);
        eventCounts[channel] = new TH1D(outHistName.c_str(),"ChannelCounts",10,0.5,10.5);
    }

    ReportPostBegin();
}

Bool_t DileptonSelector::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
    outFile->cd();
    eventWeight = 1.;
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

    ///////////////////////
    // Generator objects //
    ///////////////////////

    vector<TGenParticle*> genParticles;
    vector<int> genMotherId;

    float topSF = 1.;
    topPtWeight = 1.;
    if (!isData) {

        // Set data period for 2016 MC scale factors
        if (rng->Rndm() < 0.468) {
            weights->SetDataPeriod("2016BtoF");    
        } else {
            weights->SetDataPeriod("2016GH");
        }

        // save gen weight for amc@nlo Drell-Yan sample
        genWeight = fGenEvtInfo->weight > 0 ? 1 : -1;
        if (genWeight < 0) {
            hTotalEvents->Fill(10);
        }

        // loop over gen particle collection:
        //   * parton counting for combining inclusive and jet-binned Drell-Yan samples
        //   * categorization of W and tau decays in ttbar/tW events
        //   * calculate scale factors for top pt reweighting
        unsigned count = 0;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);


            //if (abs(particle->pdgId) == 6) {
            //    cout << i  << ", " << particle->parent << ", " << particle->pdgId << ", " << particle->status;
            //    cout << "\t" << particle->pt << ", " << particle->eta;
            //    cout << endl;
            //}

            // parton counting for jet-binned Drell-Yan samples
            if (
                    particle->status == 23 
                    && (abs(particle->pdgId) < 6 || particle->pdgId == 21) 
                    && particle->parent != -2
               ) {
                ++count;
            }

            // top pt reweighting: get the scale factor based on the top quark pt
            if (abs(particle->pdgId) == 6 && particle->status == 62) {
                topSF *= exp(0.0615 - 0.0005*particle->pt);
            }

            // This will find final state leptons (fix this to discriminate between decays)
            if (
                    (abs(particle->pdgId) == 11 || abs(particle->pdgId) == 13) 
                    && particle->status == 1 
                    && particle->parent != -2
               ) {

                // Find if the lepton comes from a top quark
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                //cout << mother->pdgId << ", ";
                int origin = abs(mother->pdgId);
                int intermediary = origin;
                while (origin != 6 && mother->parent != -2) {
                    mother = (TGenParticle*) fGenParticleArr->At(mother->parent);
                    origin = abs(mother->pdgId);
                    if (origin <= 5) { // remove leptons that have been radiated from light quarks
                        intermediary = -1; 
                        break;
                    } else if (origin == 15) {
                        intermediary = 15;
                    } else if (intermediary != 15 && origin == 24) {
                        intermediary = 24;
                    }
                }

                if (origin == 6 && (intermediary == 24 || intermediary == 15)) {
                    genParticles.push_back(particle);
                    genMotherId.push_back(intermediary);
                }
            }

        }
        nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY

        // Account for the top pt weights
        if (params->datasetgroup.substr(0, 5) == "ttbar") {
            topPtWeight *= sqrt(topSF);
        }

        // pileup reweighting
        nPU          = fInfo->nPUmean;
        puWeight     = weights->GetPUWeight(fInfo->nPUmean); 
        eventWeight *= puWeight;

        // top pt reweighting
        eventWeight *= topPtWeight;

        // keep track of total number of generated events with weights
        hTotalEvents->Fill(9, topPtWeight);

    } else {
        nPU = 0;
        nPartons = 0;
    }

    /* Apply lumi mask */
    if (isData) {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(2);

    /* Trigger selection */
    bool passTrigger = false;
    vector<string> passTriggerNames;
    for (unsigned i = 0; i < triggerNames.size(); ++i) {
        bool triggered = false;
        triggered = trigger->pass(triggerNames[i], fInfo->triggerBits);
        passTrigger |= triggered;

        if (triggered) {
            passTriggerNames.push_back(triggerNames[i]);

            // remove overlap between single electron and single muon datastreams
            if (//isData
                    params->selection == "single_lepton" 
                    && params->datasetgroup.substr(0, 8)  == "electron" 
                    && (triggerNames[i] == "HLT_IsoMu24_v*" || triggerNames[i] == "HLT_IsoTkMu24_v*")
               ) {
                passTrigger = false;
                break;
            }
        }
    }

    if (!passTrigger)
        return kTRUE;
    hTotalEvents->Fill(3);


    /////////////////////
    // Fill event info //
    /////////////////////

    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
    nPV           = fPVArr->GetEntries();
    triggerStatus = passTrigger;

    ///////////////////
    // Select objects//
    ///////////////////

    /* Vertices */
    if (fInfo->hasGoodPV) {
        assert(fPVArr->GetEntries() != 0);
        TVector3 pv;
        copy_xyz((TVertex*) fPVArr->At(0), pv);
        rPV = pv; 
        particleSelector->SetPV(pv);
    } else {
        return kTRUE;
    }
    hTotalEvents->Fill(4);
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);

    /* MUONS */
    vector<TMuon*> muons, fail_muons;
    vector<TLorentzVector> veto_muons;
    for (int i=0; i < fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);

        // Apply rochester muon momentum corrections
        TLorentzVector muonP4;
        copy_p4(muon, MUON_MASS, muonP4);
        double muonSF = 1.;
        if (isData) {
            muonSF = muonCorr->kScaleDT(muon->q, muon->pt, muon->eta, muon->phi, 0, 0);
        } else {
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
                muonP4.Pt() > 10.
                && fabs(muonP4.Eta()) < 2.4

                // tight muon ID and ISO
                && (muon->typeBits & baconhep::kPFMuon) 
                && (muon->typeBits & baconhep::kGlobal) 
                && muon->muNchi2    < 10.
                && muon->nMatchStn  > 1
                && muon->nPixHits   > 0
                && muon->nTkLayers  > 5 
                && muon->nValidHits > 0
                && fabs(muon->d0)   < 0.2
                && fabs(muon->dz)   < 0.5
           ) {
            if (GetMuonIsolation(muon)/muonP4.Pt() < 0.15) {
                muons.push_back(muon);
                veto_muons.push_back(muonP4);
            } else {
                fail_muons.push_back(muon);
            }
        }

        // muons for jet veto
        //if (
        //        muonP4.Pt() > 10
        //        && fabs(muonP4.Eta()) < 2.4
        //        // tight muon ID and ISO
        //        && (muon->pogIDBits & baconhep::kPOGTightMuon)
        //        && GetMuonIsolation(muon)/muonP4.Pt() < 0.15
        //   ) {
        //    veto_muons.push_back(muonP4);
        //}
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);
    sort(fail_muons.begin(), fail_muons.end(), sort_by_higher_pt<TMuon>);

    /* ELECTRONS */
    vector<TElectron*> electrons, fail_electrons;
    vector<TLorentzVector> veto_electrons;
    //float eScale = 1.;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        if (isData) {
            scaleData sdata = electronScaler->GetScaleData(electron, runNumber);
            electron->pt *= sdata.scale;
        } else {
            float sFactor = electronScaler->GetSmearingFactor(electron, 0, 0);
            electron->pt *= rng->Gaus(1, sFactor);
        }

        TLorentzVector electronP4;
        electronP4.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, 511e-6);

        if (
                electron->pt > 10
                && fabs(electron->scEta) < 2.5
                && particleSelector->PassElectronID(electron, cuts->tightElID)
           ) {
            if (particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)) {
                electrons.push_back(electron);
                veto_electrons.push_back(electronP4);
            } else {
                fail_electrons.push_back(electron);
            }
        }
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);
    sort(fail_electrons.begin(), fail_electrons.end(), sort_by_higher_pt<TElectron>);


    /* JETS */
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets;
    TLorentzVector hadronicP4;
    float sumJetPt = 0;
    nJets    = 0;
    nBJets   = 0;
    nFwdJets = 0;
    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

        // apply JEC offline and get scale uncertainties
        double jec = particleSelector->JetCorrector(jet, "NONE");
        jet->pt = jet->ptRaw*jec;


        float gRand = 1.;
        if (!isData) { // apply jet energy resolution corrections to simulation
            pair<float, float> resPair = particleSelector->JetResolutionAndSF(jet, 0);
            gRand = rng->Gaus(0, resPair.first);
            float jerc = 1 + gRand*sqrt(std::max((double)resPair.second*resPair.second - 1, 0.));
            jet->pt = jet->pt*jerc;
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
                fabs(jet->eta) < 4.7
                && particleSelector->PassJetID(jet, cuts->looseJetID)
                && !muOverlap 
                && !elOverlap
           ) {

            if (fabs(jet->eta) <= 2.4) { 

                if (jet->pt > 30) {
                    ++nJets;
                    if (jet->bmva > 0.9432) { 
                        ++nBJets;
                    } 
                }

                if (jet->pt > 30) {
                    hadronicP4 += jetP4;
                    sumJetPt += jetP4.Pt();
                    jets.push_back(jet);
                }
 
            } else {
                if (fabs(jet->eta) > 2.4 && jet->pt > 30) {
                    hadronicP4 += jetP4;
                    sumJetPt += jetP4.Pt();
                    ++nFwdJets;
                }
            }
        }
    }
    std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);

    /* MET */
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;

    TVector2 metP2;
    metP2.SetMagPhi(met, metPhi);

    /* HT */
    htSum = sumJetPt;
    ht    = hadronicP4.Pt();
    htPhi = hadronicP4.Phi();

    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    nMuons     = muons.size();
    nElectrons = electrons.size();

    // trigger selections
    bool muonTriggered = find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_IsoMu24_v*") != passTriggerNames.end() 
        || find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_IsoTkMu24_v*") != passTriggerNames.end();
    bool electronTriggered = find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_Ele27_WPTight_Gsf_v*") != passTriggerNames.end();

    string channel = "";
    if (muons.size() >= 2 && electrons.size() == 0) { // mu+mu selection
        channel = "mumu";
        eventCounts[channel]->Fill(1);

        // convert to TLorentzVectors
        TLorentzVector muonOneP4, muonTwoP4, dimuonP4;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, 0.1052);
        dimuonP4 = muonOneP4 + muonTwoP4;

        float muonOneIso = GetMuonIsolation(muons[0]);
        float muonTwoIso = GetMuonIsolation(muons[1]);
        if (
                muonOneIso/muonOneP4.Pt() > 0.15 
                || muonTwoIso/muonTwoP4.Pt() > 0.15 
                || !muonTriggered
           )
            return kTRUE;
        eventCounts[channel]->Fill(2);

        if (muons[0]->pt < 25. || muons[1]->pt < 10.)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        leptonOneP4     = muonOneP4;
        leptonOneIso    = muonOneIso;
        leptonOneFlavor = muons[0]->q*13;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;

        leptonTwoP4     = muonTwoP4;
        leptonTwoIso    = muonTwoIso;
        leptonTwoFlavor = muons[1]->q*13;
        leptonTwoDZ     = muons[1]->dz;
        leptonTwoD0     = muons[1]->d0;

        if (!isData) {
            leptonOneMother = GetGenMotherId(genParticles, muonOneP4);
            leptonTwoMother = GetGenMotherId(genParticles, muonTwoP4);

            // reconstruction weights
            EfficiencyContainer effCont1, effCont2;
            effCont1 = weights->GetMuonRecoEff(muonOneP4);
            effCont2 = weights->GetMuonRecoEff(muonTwoP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;

            effs = effCont2.GetEff();
            errs = effCont2.GetErr();
            leptonTwoRecoWeight = effs.first/effs.second;

            // trigger weights with trigger matching:
            //
            // check if lepton could pass the trigger threshold and is matched
            // to a trigger object.  When both muons pass the trigger, use the
            // efficiency for detecting either
            bitset<2> triggered;
            for (const auto& name: passTriggerNames) {
                if (trigger->passObj(name, 1, muons[0]->hltMatchBits) && muonOneP4.Pt() > 25)
                    triggered.set(0);
                if (trigger->passObj(name, 1, muons[1]->hltMatchBits) && muonTwoP4.Pt() > 25)
                    triggered.set(1);
            }

            if (triggered.all()) {
                effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
                effCont2      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
                triggerWeight = GetTriggerSF(effCont1, effCont2);
            } else if (triggered.test(0)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
            } else if (triggered.test(1)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
            } else {
                return kTRUE;
            }

            // update the event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;
        }
    } else if (electrons.size() >= 2 && muons.size() == 0) { // e+e selection
        channel = "ee";
        eventCounts[channel]->Fill(1);

        if (electrons[0]->pt < 30 || electrons[1]->pt < 10 || !electronTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        TLorentzVector electronOneP4, electronTwoP4, dielectronP4;
        electronOneP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        electronTwoP4.SetPtEtaPhiM(electrons[1]->pt, electrons[1]->eta, electrons[1]->phi, 511e-6);
        dielectronP4 = electronOneP4 + electronTwoP4;

        if (dielectronP4.M() < 12.)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        leptonOneP4     = electronOneP4;
        leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonOneFlavor = 11*electrons[0]->q;
        leptonOneDZ     = electrons[0]->dz;
        leptonOneD0     = electrons[0]->d0;

        leptonTwoP4     = electronTwoP4;
        leptonTwoIso    = GetElectronIsolation(electrons[1], fInfo->rhoJet);
        leptonTwoFlavor = 11*electrons[1]->q;
        leptonTwoDZ     = electrons[1]->dz;
        leptonTwoD0     = electrons[1]->d0;

        if (!isData) {
            leptonOneMother = GetGenMotherId(genParticles, electronOneP4);
            leptonTwoMother = GetGenMotherId(genParticles, electronTwoP4);

            EfficiencyContainer effCont1, effCont2;
            effCont1 = weights->GetElectronRecoEff(electronOneP4);
            effCont2 = weights->GetElectronRecoEff(electronTwoP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;

            effs = effCont2.GetEff();
            errs = effCont2.GetErr();
            leptonTwoRecoWeight = effs.first/effs.second;

            // trigger weights with trigger matching
            bitset<2> triggered;
            for (const auto& name: passTriggerNames) {
                if (trigger->passObj(name, 1, electrons[0]->hltMatchBits) && electronOneP4.Pt() > 30)
                    triggered.set(0);
                if (trigger->passObj(name, 1, electrons[1]->hltMatchBits) && electronTwoP4.Pt() > 30)
                    triggered.set(1);
            }

            if (triggered.all()) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronOneP4);
                effCont2      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronTwoP4);
                triggerWeight = GetTriggerSF(effCont1, effCont2);
            } else if (triggered.test(0)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronOneP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
            } else if (triggered.test(1)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronTwoP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
            } else {
                return kTRUE;
            }

            // update event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;
        }
    } else {
        return kTRUE;
    }

    ///////////////////
    // Fill jet info //
    ///////////////////

    if (jets.size() > 0) {
        jetOneP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
        jetOneTag    = jets[0]->csv;
        jetOneFlavor = jets[0]->hadronFlavor;
    } else {
        jetOneP4.SetPtEtaPhiM(0., 0., 0., 0.);
        jetOneTag    = 0.;
        jetOneFlavor = 0;
    }

    if (jets.size() > 1) {
        jetTwoP4.SetPtEtaPhiM(jets[1]->pt, jets[1]->eta, jets[1]->phi, jets[1]->mass);
        jetTwoTag    = jets[1]->csv;
        jetTwoFlavor = jets[1]->hadronFlavor;
    } else {
        jetTwoP4.SetPtEtaPhiM(0., 0., 0., 0.);
        jetTwoTag    = 0.;
        jetTwoFlavor = 0;
    } 

    if (!isData && genParticles.size() == 2) {
        genOneId = genParticles[0]->pdgId;
        genOneP4.SetPtEtaPhiM(genParticles[0]->pt, genParticles[0]->eta, genParticles[0]->phi, genParticles[0]->mass); 
        genTwoId = genParticles[1]->pdgId;
        genTwoP4.SetPtEtaPhiM(genParticles[1]->pt, genParticles[1]->eta, genParticles[1]->phi, genParticles[1]->mass); 
    } else {
        genOneP4.SetPtEtaPhiM(0., 0., 0., 0.); 
        genTwoP4.SetPtEtaPhiM(0., 0., 0., 0.); 
    }

    outFile->cd(channel.c_str());
    outTrees[channel]->Fill();
    this->passedEvents++;
    return kTRUE;
}

void DileptonSelector::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void DileptonSelector::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void DileptonSelector::ReportPostTerminate()
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
    std::unique_ptr<DileptonSelector> selector(new DileptonSelector());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

float DileptonSelector::GetMuonIsolation(const baconhep::TMuon* mu)
{
    float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    return combIso;
}

float DileptonSelector::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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

int DileptonSelector::GetGenMotherId(vector<TGenParticle*> particles, TLorentzVector p4)
{
    int motherId = 0;
    for (unsigned i = 0; i < particles.size(); ++i) {
        TLorentzVector genP4;
        genP4.SetPtEtaPhiM(particles[i]->pt, particles[i]->eta, particles[i]->phi, particles[i]->mass); 
        if (genP4.DeltaR(p4) < 0.3) {
            TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particles[i]->parent);
            motherId = mother->pdgId;
        }
    }
    return motherId;
}

float DileptonSelector::GetTriggerSF(EfficiencyContainer eff1, EfficiencyContainer eff2)
{
    pair<double, double> trigEff1, trigEff2;
    trigEff1 = eff1.GetEff();
    trigEff2 = eff2.GetEff();

    float sf = 1.;
    if (trigEff1.second > 0 || trigEff2.second > 0) {
        sf = (1 - (1 - trigEff1.first)*(1 - trigEff2.first))/(1 - (1 - trigEff1.second)*(1 - trigEff2.second));
    }

    return sf;
}

float DileptonSelector::GetTriggerSFError(EfficiencyContainer eff1, EfficiencyContainer eff2)
{
    pair<double, double> trigEff1, trigEff2, trigErr1, trigErr2;
    trigEff1 = eff1.GetEff();
    trigEff2 = eff2.GetEff();
    trigErr1 = eff1.GetErr();
    trigErr2 = eff2.GetErr();

    //cout << trigEff1.first << " :: " << trigEff1.second << " :: "
    //     << trigEff2.first << " :: " << trigEff2.second << endl;

    //cout << trigErr1.first << " :: " << trigErr1.second << " :: "
    //     << trigErr2.first << " :: " << trigErr2.second << endl;

    float denom = 1 - (1 - trigEff1.second)*(1 - trigEff2.second);
    float dEffData1 = (1 - trigEff2.first)/denom;
    float dEffData2 = (1 - trigEff1.first)/denom;
    float dEffMC1   = -(1 - trigEff2.second)/pow(denom, 2);
    float dEffMC2   = -(1 - trigEff1.second)/pow(denom, 2);

    float sfVar = pow(trigErr1.first*dEffData1, 2)
        + pow(trigErr2.first*dEffData2, 2)
        + pow(trigErr1.second*dEffMC1, 2)
        + pow(trigErr2.second*dEffMC2, 2);
    return sfVar;
}
