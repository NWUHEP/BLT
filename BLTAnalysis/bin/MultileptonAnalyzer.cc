#include "MultileptonAnalyzer.h"
#include <map>

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 
bool sort_by_btag(const baconhep::TJet* lhs, const baconhep::TJet* rhs) 
{
    return lhs->bmva > rhs->bmva;
}


MultileptonAnalyzer::MultileptonAnalyzer() : BLTSelector()
{

}

MultileptonAnalyzer::~MultileptonAnalyzer()
{

}

void MultileptonAnalyzer::Begin(TTree *tree)
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

    if (params->selection == "mumu") {
        //triggerNames.push_back("HLT_IsoMu24_v*");
        //triggerNames.push_back("HLT_IsoTkMu24_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");

    } else if (params->selection == "ee") {
        triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
    } else if (params->selection == "emu") {
        triggerNames.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
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
    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);
    outTree->Branch("ht", &ht);
    outTree->Branch("htPhi", &htPhi);

    // leptons
    outTree->Branch("leptonOneP4", &leptonOneP4);
    outTree->Branch("leptonOneIso", &leptonOneIso);
    outTree->Branch("leptonOneFlavor", &leptonOneFlavor);
    outTree->Branch("leptonOneD0", &leptonOneD0);
    outTree->Branch("leptonOneDZ", &leptonOneDZ);

    outTree->Branch("leptonTwoP4", &leptonTwoP4);
    outTree->Branch("leptonTwoIso", &leptonTwoIso);
    outTree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
    outTree->Branch("leptonTwoD0", &leptonTwoD0);
    outTree->Branch("leptonTwoDZ", &leptonTwoDZ);

    // lepton vertices
    outTree->Branch("dileptonVertexOne", &dileptonVertexOne);
    outTree->Branch("dileptonVertexErrOne", &dileptonVertexErrOne);
    outTree->Branch("dileptonVertexChi2One", &dileptonVertexChi2One);
    outTree->Branch("dileptonVertexDOFOne", &dileptonVertexDOFOne);

    if (params->selection == "4l") {
        outTree->Branch("leptonThreeP4", &leptonThreeP4);
        outTree->Branch("leptonThreeIso", &leptonThreeIso);
        outTree->Branch("leptonThreeFlavor", &leptonThreeFlavor);
        outTree->Branch("leptonThreeD0", &leptonThreeD0);
        outTree->Branch("leptonThreeDZ", &leptonThreeDZ);

        outTree->Branch("leptonFourP4", &leptonFourP4);
        outTree->Branch("leptonFourIso", &leptonFourIso);
        outTree->Branch("leptonFourFlavor", &leptonFourFlavor);
        outTree->Branch("leptonThreeD0", &leptonThreeD0);
        outTree->Branch("leptonThreeDZ", &leptonThreeDZ);

        outTree->Branch("dileptonVertexTwo", &dileptonVertexTwo);
        outTree->Branch("dileptonVertexErrTwo", &dileptonVertexErrTwo);
        outTree->Branch("dileptonVertexChi2Two", &dileptonVertexChi2Two);
        outTree->Branch("dileptonVertexDOFTwo", &dileptonVertexDOFTwo);
    }

    // jets
    outTree->Branch("jetOneP4", &jetOneP4);
    outTree->Branch("jetOneTag", &jetOneTag);

    outTree->Branch("jetTwoP4", &jetTwoP4);
    outTree->Branch("jetTwoTag", &jetTwoTag);

    // gen level objects
    outTree->Branch("gnuOneP4", &gnuOneP4);
    outTree->Branch("gnuTwoP4", &gnuTwoP4);

    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nFwdJets", &nFwdJets);
    outTree->Branch("nBJets", &nBJets);

    // event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);

    ReportPostBegin();
}

Bool_t MultileptonAnalyzer::Process(Long64_t entry)
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

    /* Apply lumi mask */
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

    if (!passTrigger)// && isData)
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

    vector<TLorentzVector> gNeutrinos;
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

            if (
                    //particle->status == 23 
                    //&& (abs(particle->pdgId) < 12 || particle->pdgId == 14 || particle->pdgId == 16) 
                    (abs(particle->pdgId) == 12 || abs(particle->pdgId) == 14 || abs(particle->pdgId) == 16) 
                    //&& particle->parent != -2
               ) {
                //cout << particle->status << ", "
                //    << particle->pdgId  << ", "
                //    << particle->parent;

                if (particle->parent != -2) {
                    TGenParticle* mom = (TGenParticle*) fGenParticleArr->At(particle->parent);
                    if (abs(mom->pdgId) == 23 || abs(mom->pdgId) == 24) {
                        TLorentzVector nu;
                        nu.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, 0.);
                        gNeutrinos.push_back(nu);

                        //cout << ", " << mom->pdgId;
                        //cout << "\t" << particle->pt << ", " << particle->eta;
                    }
                }
                //cout << endl;

            }
        }
        nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY
        //cout << nPartons << "\n" << endl;
    } else {
        nPartons = 0;
    }
    //cout << endl;

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
    vector<TMuon*> muons;
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
        muonP4.SetPtEtaPhiM(muonSF*muon->pt, muon->eta, muon->phi, MUON_MASS);

        if (
                muon->pt > 3 
                && fabs(muon->eta) < 2.4
                // Loose muon ID
                //&& (muon->typeBits & baconhep::kPFMuon)
                //&& ((muon->typeBits & baconhep::kGlobal) || (muon->typeBits & baconhep::kTracker))
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
           ) {
            muons.push_back(muon);

            // muons for jet veto
            if (muonP4.Pt() > 10) {
                veto_muons.push_back(muonP4);

                // trigger matching
                bool triggered = false;
                for (unsigned i = 0; i < triggerNames.size(); ++i) {
                    triggered |= trigger->passObj(triggerNames[i], 1, muon->hltMatchBits);
                }
            }
        }
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);

    /* ELECTRONS */
    std::vector<TElectron*> electrons;
    vector<TLorentzVector> veto_electrons;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        TLorentzVector electronP4;
        electronP4.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, 511e-6);

        if (
                electron->pt > 10 
                && fabs(electron->eta) < 2.5
                && particleSelector->PassElectronID(electron, cuts->tightElID)
                && particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
           ) {
            electrons.push_back(electron);
            veto_electrons.push_back(electronP4);
        }
    }

    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);

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
            if (jetP4.DeltaR(mu) < 0.5) {
                muOverlap = true;
                break;
            }
        }
        bool elOverlap = false;
        for (const auto& el: veto_electrons) {
            if (jetP4.DeltaR(el) < 0.5) {
                elOverlap = true;
                break;
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
                        && !elOverlap
                   ) { 
                    hadronicP4 += jetP4;
                    jets.push_back(jet);

                    if (isData) {
                        if (jet->bmva > 0.9432) { 
                            ++nBJets;
                        } else {
                            ++nJets;
                        }
                    } else {
                        if (particleSelector->BTagModifier(jet, "MVAT")) { 
                            ++nBJets;
                        } else {
                            ++nJets;
                        }
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

    std::sort(jets.begin(), jets.end(), sort_by_btag);

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

    if (params->selection == "mumu") {
        if (muons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        if (muons[0]->pt < 20. || muons[1]->pt < 10.)
            return kTRUE;

        hTotalEvents->Fill(6);

        // Find leading positive and negatively charged muons and convert to TLorentzVectors
        TLorentzVector muonOneP4, muonTwoP4, dimuonP4;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, 0.1052);
        dimuonP4 = muonOneP4 + muonTwoP4;
        //cout << muonOneP4.Pt() << ", " << muonTwoP4.Pt() << ", " <<  dimuonP4.M() << endl;

        if (dimuonP4.M() < 12.)
            return kTRUE;
        hTotalEvents->Fill(7);

        if (nJets + nBJets < 2)
            return kTRUE;
        hTotalEvents->Fill(8);

        leptonOneP4     = muonOneP4;
        leptonOneIso    = GetMuonIsolation(muons[0]);
        leptonOneFlavor = muons[0]->q*13;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;

        leptonTwoP4     = muonTwoP4;
        leptonTwoIso    = GetMuonIsolation(muons[1]);
        leptonTwoFlavor = muons[1]->q*13;
        leptonTwoDZ     = muons[1]->dz;
        leptonTwoD0     = muons[1]->d0;

        // Check for dimuon vertex (only works on a limited set of bacon data)
        /*unsigned leptonOneIndex = muons[0]->muIndex;
          unsigned leptonTwoIndex = muons[1]->muIndex;
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
          dileptonVertexOne     = vertPos;
          dileptonVertexErrOne  = vertErr;
          dileptonVertexChi2One = dimuonVert->chi2;
          dileptonVertexDOFOne  = dimuonVert->ndof;
          hasValidVertex        = true;
          break;
          }
          }

          if (!hasValidVertex)
          return kTRUE;
          hTotalEvents->Fill(8);*/

        if (!isData) {
            eventWeight *= weights->GetMuonIDEff(muonOneP4);
            eventWeight *= weights->GetMuonISOEff(muonOneP4);
            eventWeight *= weights->GetMuonIDEff(muonTwoP4);
            eventWeight *= weights->GetMuonISOEff(muonTwoP4);

            // trigger weight
            //pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
            //pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
            //eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);

        }


    } else if (params->selection == "ee") {

        if (electrons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        if (electrons[0]->pt < 25 || electrons[1]->pt < 15)
            return kTRUE;
        hTotalEvents->Fill(6);

        TLorentzVector electronOneP4, electronTwoP4, dielectronP4;
        electronOneP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        electronTwoP4.SetPtEtaPhiM(electrons[1]->pt, electrons[1]->eta, electrons[1]->phi, 511e-6);
        dielectronP4 = electronOneP4 + electronTwoP4;
        if (dielectronP4.M() < 12.)
            return kTRUE;
        hTotalEvents->Fill(7);

        if (nJets + nBJets < 2)
            return kTRUE;
        hTotalEvents->Fill(8);

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
            eventWeight *= weights->GetElectronRecoIdEff(electronOneP4);
            eventWeight *= weights->GetElectronRecoIdEff(electronTwoP4);

            // trigger weight
            //pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
            //pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
            //eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);
        }

    } else if (params->selection == "emu") {

        if (muons.size() != 1 || electrons.size() != 1)
            return kTRUE;
        hTotalEvents->Fill(5);

        if (muons[0]->pt < 25 || electrons[0]->pt)
           return kTRUE;
        hTotalEvents->Fill(6);

        TLorentzVector muonP4, electronP4, dilepton;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        electronP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        dilepton = muonP4 + electronP4;
        if (dilepton.M() < 12)
            return kTRUE;
        hTotalEvents->Fill(7);

        leptonOneP4     = muonP4;
        leptonOneIso    = GetMuonIsolation(muons[0]);
        leptonOneFlavor = muons[0]->q*13;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;

        leptonTwoP4     = electronP4;
        leptonTwoIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonTwoFlavor = 11*electrons[0]->q;
        leptonTwoDZ     = electrons[0]->dz;
        leptonTwoD0     = electrons[0]->d0;

        if (!isData) {
            eventWeight *= weights->GetMuonIDEff(muonP4);
            eventWeight *= weights->GetElectronRecoIdEff(electronP4);

            // trigger efficiency
            //pair<float, float> trigEff = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muons[0]);
            //eventWeight *= trigEff.first;
        }
    } else if (params->selection == "4l") {

        if (muons.size() >= 4) { 

            if (muons[0]->pt < 25)
                return kTRUE;
            hTotalEvents->Fill(6);

            /*leptonOneP4        = muons[0];
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
              leptonFourFlavor   = 1;
              leptonFourTrigger  = electrons_trigger[1];*/
        } else {
            hTotalEvents->Fill(5);
            return kTRUE;
        }
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

    if (!isData && gNeutrinos.size() == 2) {
        gnuOneP4.SetPtEtaPhiM(gNeutrinos[0].Pt(), gNeutrinos[0].Eta(), gNeutrinos[0].Phi(), 0.); 
        gnuTwoP4.SetPtEtaPhiM(gNeutrinos[1].Pt(), gNeutrinos[1].Eta(), gNeutrinos[1].Phi(), 0.); 
    } else {
        gnuOneP4.SetPtEtaPhiM(0., 0., 0., 0.); 
        gnuTwoP4.SetPtEtaPhiM(0., 0., 0., 0.); 
    }

    outTree->Fill();
    this->passedEvents++;
    return kTRUE;
}

void MultileptonAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void MultileptonAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void MultileptonAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<MultileptonAnalyzer> selector(new MultileptonAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

float MultileptonAnalyzer::MetKluge(float met)
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

float MultileptonAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    return combIso;
}

float MultileptonAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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


