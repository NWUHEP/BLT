#include "upsilonGammaAnalyzer.h"
#include <map>
#include <fstream>
#include <math.h>

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

upsilonGammaAnalyzer::upsilonGammaAnalyzer() : BLTSelector()
{

}

upsilonGammaAnalyzer::~upsilonGammaAnalyzer()
{

}

void upsilonGammaAnalyzer::Begin(TTree *tree)
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

    if (params->selection == "mumu" || params->selection == "mumug") {
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
    }
    else if (params->selection == "ee" || params->selection == "elelg") {
        triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
        triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
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
    outTree->Branch("puWeight", &puWeight);
    outTree->Branch("genWeight", &genWeight);
    outTree->Branch("triggerWeight", &triggerWeight);
    outTree->Branch("nPV", &nPV);
    outTree->Branch("nPU", &nPU);
    outTree->Branch("nPartons", &nPartons);
    outTree->Branch("rPV", &rPV);

    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);
    outTree->Branch("ht", &ht);
    outTree->Branch("htPhi", &htPhi);
    outTree->Branch("htSum", &htSum);

    // leptons
    outTree->Branch("leptonOneP4", &leptonOneP4);
    outTree->Branch("leptonOneIso", &leptonOneIso);
    outTree->Branch("leptonOneFlavor", &leptonOneFlavor);
    outTree->Branch("leptonOneMother", &leptonOneMother);
    outTree->Branch("leptonOneD0", &leptonOneD0);
    outTree->Branch("leptonOneDZ", &leptonOneDZ);
    outTree->Branch("leptonOneRecoWeight", &leptonOneRecoWeight);

    outTree->Branch("leptonTwoP4", &leptonTwoP4);
    outTree->Branch("leptonTwoIso", &leptonTwoIso);
    outTree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
    outTree->Branch("leptonTwoMother", &leptonTwoMother);
    outTree->Branch("leptonTwoD0", &leptonTwoD0);
    outTree->Branch("leptonTwoDZ", &leptonTwoDZ);
    outTree->Branch("leptonTwoRecoWeight", &leptonTwoRecoWeight);

    // photons
    outTree->Branch("photonOneP4", &photonOneP4);
    outTree->Branch("photonOneR9", &photonOneR9);
    outTree->Branch("photonOneMVA", &photonOneMVA);
    outTree->Branch("passElectronVeto", &passElectronVeto);

    // lepton vertices
    outTree->Branch("dileptonVertex", &dileptonVertex);
    outTree->Branch("dileptonVertexErr", &dileptonVertexErr);
    outTree->Branch("dileptonVertexChi2", &dileptonVertexChi2);
    outTree->Branch("dileptonVertexDOF", &dileptonVertexDOF);
    outTree->Branch("dileptonVertexProb", &dileptonVertexProb);

    // jets
    outTree->Branch("jetOneP4", &jetOneP4);
    outTree->Branch("jetOneTag", &jetOneTag);
    outTree->Branch("jetTwoP4", &jetTwoP4);
    outTree->Branch("jetTwoTag", &jetTwoTag);

    // gen level objects 
    outTree->Branch("genLeptonOneP4", &genLeptonOneP4);
    outTree->Branch("genLeptonOneId", &genLeptonOneId);
    outTree->Branch("genLeptonTwoP4", &genLeptonTwoP4);
    outTree->Branch("genLeptonTwoId", &genLeptonTwoId);
    outTree->Branch("genPhotonP4", &genPhotonP4);
    outTree->Branch("genPhotonFHPFS", &genPhotonFHPFS);
    outTree->Branch("genPhotonIPFS", &genPhotonIPFS);

    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nTaus", &nTaus);
    outTree->Branch("nPhotons", &nPhotons);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nFwdJets", &nFwdJets);
    outTree->Branch("nCentralJets", &nCentralJets);
    outTree->Branch("nBJets", &nBJets);

    // event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",30,0.5,30.5);

    ReportPostBegin();
}

Bool_t upsilonGammaAnalyzer::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;
    hTotalEvents->Fill(1);
    
    const bool isData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isData);
    
    genWeight = 1;
    if (!isData) {
        if (fGenEvtInfo->weight < 0) {
            genWeight = -1;
            int maxBin = hTotalEvents->GetSize() - 2;
            hTotalEvents->Fill(maxBin);
        }
    }

    if (entry%10000==0)  
        std::cout << "... Processing event " << entry 
            << " Run: " << fInfo->runNum 
            << " Lumi: " << fInfo->lumiSec 
            << " Event: " << fInfo->evtNum 
            << std::endl; 
          
    ///////////////////////
    // Generator objects // 
    ///////////////////////

    vector<TGenParticle*> genLeptons;
    vector<TGenParticle*> genPhotons;
    if (!isData) {
        unsigned count = 0;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            //cout << i  << ", " << particle->parent << ", " << particle->pdgId << ", " << particle->status;
            //cout << "\t" << particle->pt << ", " << particle->eta;
            //cout << endl;

            if (
                    particle->status == 23 
                    && (fabs(particle->pdgId) < 6 || particle->pdgId == 21) 
                    && particle->parent != -2
               ) {
                ++count;
            }

            // saving leptons 
            if ((abs(particle->pdgId) == 11 || abs(particle->pdgId) == 13) && 
                (particle->isPromptFinalState || particle->fromHardProcessFinalState)) 
                genLeptons.push_back(particle);    

            // saving photons
            if (abs(particle->pdgId) == 22 && (particle->isPromptFinalState || particle->fromHardProcessFinalState))
                genPhotons.push_back(particle);         
            }

            nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY

    } else {
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
        }
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
        puWeight = weights->GetPUWeight(nPU); // pileup reweighting
        eventWeight *= puWeight;
    } else {
        nPU = 0;

    }    

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

        muon->pt = muonSF*muon->pt; 
        muonP4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);

        if (     
                   // tight muon ID and ISO
               //    muonP4.Pt() > 5.
               //    && fabs(muonP4.Eta()) < 2.4
               //    && (muon->typeBits & baconhep::kPFMuon) 
               //    && (muon->typeBits & baconhep::kGlobal) 
               //    && muon->muNchi2    < 10.
               //    && muon->nMatchStn  > 1
               //    && muon->nPixHits   > 0
               //    && fabs(muon->d0)   < 0.2
               //    && fabs(muon->dz)   < 0.5
               //    && muon->nTkLayers  > 5 
               //    && muon->nValidHits > 0
               //    //&& GetMuonIsolation(muon)/muonP4.Pt() < 0.25
               //){
                //muons.push_back(muon);
            //}

                // h->ZZ->4l "tight" ID and pf_isorel < 0.35
                muonP4.Pt() > 5.
                && fabs(muonP4.Eta()) < 2.4
                && fabs(muon->d0) < 0.5
                && fabs(muon->dz) < 1.0
                && (((muon->typeBits & baconhep::kGlobal) || 
                   ((muon->typeBits & baconhep::kTracker) && muon->nMatchStn > 0)) &&
                   (muon->btt != 2)) // Global muon or (arbitrated) tracker muon
                && fabs(muon->sip3d) < 4.0                 
                && GetMuonIsolation(muon)/muonP4.Pt() < 0.35

                ) {
                    // We now have h->ZZ->4l "loose" muons
                    if (muonP4.Pt() < 200.0) {
                        if (muon->pogIDBits & baconhep::kPOGLooseMuon)
                            muons.push_back(muon);
                    }
                    else {
                        // need to pass the tracker high-pt ID
                        if (
                                (muon->pogIDBits & baconhep::kPOGLooseMuon) ||
                                (muon->nMatchStn > 1
                                && (muon->ptErr/muon->pt) < 0.3 
                                && fabs(muon->d0) < 0.2
                                && fabs(muon->dz) < 0.5
                                && muon->nPixHits > 0
                                && muon->nTkLayers > 5)
                           )
                            muons.push_back(muon);
                    }
                }
                    
        // muons for jet veto
        if (
                muonP4.Pt() > 10
                // tight muon ID and ISO
                && fabs(muonP4.Eta()) < 2.4
                && (muon->typeBits & baconhep::kPFMuon) 
                && (muon->typeBits & baconhep::kGlobal) 
                && muon->muNchi2    < 10.
                && muon->nMatchStn  > 1
                && muon->nPixHits   > 0
                && fabs(muon->d0)   < 0.2
                && fabs(muon->dz)   < 0.5
                && muon->nTkLayers  > 5 
                && muon->nValidHits > 0
                && GetMuonIsolation(muon)/muonP4.Pt() < 0.15
           ) {
            veto_muons.push_back(muonP4);
        }
    }
    
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);

    /* ELECTRONS */
    vector<TElectron*> electrons;
    vector<TLorentzVector> veto_electrons;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        TLorentzVector electronP4;
        electronP4.SetPtEtaPhiM(electron->calibPt, electron->eta, electron->phi, ELE_MASS);

        if (
                electron->calibPt > 7
                && fabs(electron->scEta) < 2.5
                && particleSelector->PassElectronMVA(electron, cuts->hzzMVAID)
                && GetElectronIsolation(electron, fInfo->rhoJet)/electronP4.Pt() < 0.35
                && fabs(electron->d0) < 0.5
                && fabs(electron->dz) < 1.0
                && fabs(electron->sip3d) < 4.0 
           ) {
            electrons.push_back(electron);
            veto_electrons.push_back(electronP4);
        }
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);

    /* TAUS */
    vector<TTau*> taus;
    for (int i=0; i < fTauArr->GetEntries(); i++) {
         TTau *tau = (TTau*) fTauArr->At(i);
        assert(tau);

        TLorentzVector tauP4; 
        tauP4.SetPtEtaPhiM(tau->pt, tau->eta, tau->phi, tau->m);

        // Prevent overlap of muons and jets
        bool muOverlap = false;
        for (const auto& mu: veto_muons) {
            if (tauP4.DeltaR(mu) < 0.4) {
                muOverlap = true;
                break;
            }
        }
        bool elOverlap = false;
        for (const auto& el: veto_electrons) {
            if (tauP4.DeltaR(el) < 0.4) {
                elOverlap = true;
                break;
            }
        }

        if( // Selection adopted from Stany's tau counting selection for ttDM
                tau->pt > 18  
                && abs(tau->eta) < 2.3 
                && !muOverlap
                && !elOverlap
                && (tau->hpsDisc & baconhep::kByDecayModeFinding)
                && (tau->hpsDisc & baconhep::kByLooseCombinedIsolationDBSumPtCorr3Hits)
          ) {
            taus.push_back(tau);
        }
    }
    sort(taus.begin(), taus.end(), sort_by_higher_pt<TTau>);

    /* PHOTONS */
    vector <TPhoton*> photons;
    vector<TLorentzVector> veto_photons;
    for (int i=0; i<fPhotonArr->GetEntries(); i++) {
        TPhoton* photon = (TPhoton*) fPhotonArr->At(i);
        assert(photon);
        
        TLorentzVector photonP4;
        photonP4.SetPtEtaPhiM(photon->calibPt, photon->eta, photon->phi, 0.);
    
    if (
            // ID conditions
            photon->pt > 10
            && fabs(photon->scEta) < 2.5 
            && (fabs(photon->scEta) <= 1.4442 || fabs(photon->scEta) >= 1.566)
            && particleSelector->PassPhotonMVA(photon, cuts->looseMVAPhID)
            && photon->passElectronVeto
        ) {
        photons.push_back(photon);
        veto_photons.push_back(photonP4);
    }

    } 
    sort(photons.begin(), photons.end(), sort_by_higher_pt<TPhoton>);

    /* JETS */
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets;
    TLorentzVector hadronicP4;
    float sumJetPt = 0;

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

        bool phoOverlap = false;
        for (const auto& pho: veto_photons) {
            if (jetP4.DeltaR(pho) < 0.4) {
                phoOverlap = true;
                break;
            }
        }

        if (
                jet->pt > 30 
                && fabs(jet->eta) < 4.7
                && particleSelector->PassJetID(jet, cuts->looseJetID)
                && !muOverlap 
                && !elOverlap
                && !phoOverlap
           ) {
            
            jets.push_back(jet);
            ++nJets;

            if (fabs(jet->eta) <= 2.4) { 
                hadronicP4 += jetP4;
                sumJetPt += jetP4.Pt();

                if (isData) {
                    if (jet->bmva > 0.9432) { 
                        ++nBJets;
                    } else {
                        ++nCentralJets;
                    }
                } else {
                    if (particleSelector->BTagModifier(jet, "MVAT", 0, 0, rng->Uniform(1.))) { 
                        ++nBJets;
                    } else {
                        ++nCentralJets;
                    }
                }
            } else {
                if (fabs(jet->eta) > 2.5) {
                    hadronicP4 += jetP4;
                    sumJetPt += jetP4.Pt();
                    ++nFwdJets;
                }
            }
        }
    }

    //std::sort(jets.begin(), jets.end(), sort_by_btag);
    std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);

    /* MET */
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;

    /* HT */
    htSum = sumJetPt;
    ht    = hadronicP4.Pt();
    htPhi = hadronicP4.Phi();

    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    nMuons     = muons.size();
    nElectrons = electrons.size();
    nTaus      = taus.size();
    nPhotons   = photons.size();

    if (params->selection == "mumu") {
        
        if (muons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        TLorentzVector muonOneP4, muonTwoP4;
        leptonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, MUON_MASS);
        leptonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, MUON_MASS);

        if (leptonOneP4.Pt() < 20.0) 
            return kTRUE;
        hTotalEvents->Fill(6);

        if (leptonTwoP4.Pt() < 10.0)
            return kTRUE;
        hTotalEvents->Fill(7);

        if (muons[0]->q == muons[1]->q)
            return kTRUE;
        hTotalEvents->Fill(8);

        leptonOneIso    = GetMuonIsolation(muons[0]);
        leptonOneFlavor = muons[0]->q*13;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;
            
        leptonTwoP4     = muonTwoP4;
        leptonTwoIso    = GetMuonIsolation(muons[1]);
        leptonTwoFlavor = muons[1]->q*13;
        leptonTwoDZ     = muons[1]->dz;
        leptonTwoD0     = muons[1]->d0;

        // fill dilepton vertex information
        dileptonVertex.SetXYZ(0., 0., 0.);
        dileptonVertexErr.SetXYZ(-1., -1., -1.);
        dileptonVertexChi2 = -1.;
        dileptonVertexDOF = -1.;
        dileptonVertexProb = -1.;
       
        for (int i = 0; i < fDimuonVertexArr->GetEntries(); ++i) {
            TVertex* dimuonVert = (TVertex*) fDimuonVertexArr->At(i);
            if (
                    ((dimuonVert->index1 == muons[0]->muIndex && dimuonVert->index2 == muons[1]->muIndex) ||
                    (dimuonVert->index1 == muons[1]->muIndex && dimuonVert->index2 == muons[0]->muIndex)) &&
                    (dimuonVert->isValid)
               ) {
                dileptonVertex.SetXYZ(dimuonVert->x, dimuonVert->y, dimuonVert->z);
                dileptonVertexErr.SetXYZ(dimuonVert->xerr, dimuonVert->yerr, dimuonVert->zerr);
                dileptonVertexChi2 = dimuonVert->chi2;
                dileptonVertexDOF = dimuonVert->ndof;
                dileptonVertexProb = dimuonVert->prob;
                break;
            }
        } 

        if (!isData) {
            eventWeight *= weights->GetHZZMuonIDEff(*muons[0]); 
            eventWeight *= weights->GetHZZMuonIDEff(*muons[1]);
            eventWeight *= weights->GetMuonISOEff(leptonOneP4);
            eventWeight *= weights->GetMuonISOEff(leptonTwoP4);

            pair<float, float> eff11 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg1", *muons[0]);
            pair<float, float> eff12 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg1", *muons[1]);
            pair<float, float> eff21 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg2", *muons[0]);
            pair<float, float> eff22 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg2", *muons[1]);
            float eff_data = eff11.first*eff22.first + eff12.first*eff21.first - eff11.first*eff12.first;
            float eff_mc = eff11.second*eff22.second + eff12.second*eff21.second - eff11.second*eff12.second;
            triggerWeight = eff_data/eff_mc;
            eventWeight *= triggerWeight;
        }

    } // end mumu selector

    else if (params->selection == "mumug") {
        if (muons.size() < 2) 
            return kTRUE;
        hTotalEvents->Fill(5);

        if (photons.size() < 1)
            return kTRUE;
        hTotalEvents->Fill(6);
        
        leptonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, MUON_MASS);
        leptonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, MUON_MASS);

        if (leptonOneP4.Pt() < 20.0) 
            return kTRUE;
        hTotalEvents->Fill(6);

        if (leptonTwoP4.Pt() < 10.0)
            return kTRUE;
        hTotalEvents->Fill(7);

        if (muons[0]->q == muons[1]->q)
            return kTRUE;
        hTotalEvents->Fill(8);
        
        // L1EMTF cut 
        if (
            fabs(leptonOneP4.DeltaPhi(leptonTwoP4)) < 70.0*(M_PI/180.0)
            && fabs(leptonOneP4.Eta()) > 1.2 
            && fabs(leptonTwoP4.Eta()) > 1.2
            && leptonOneP4.Eta()*leptonTwoP4.Eta() > 0
           )
            return kTRUE;
        hTotalEvents->Fill(9); 
        
        photonOneP4.SetPtEtaPhiM(photons[0]->calibPt, photons[0]->eta, photons[0]->phi, 0.);
        if (photonOneP4.Pt() < 15.0)
            return kTRUE;
        hTotalEvents->Fill(10);
 
        leptonOneIso    = GetMuonIsolation(muons[0]);
        leptonOneFlavor = muons[0]->q*13;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;
            
        leptonTwoIso    = GetMuonIsolation(muons[1]);
        leptonTwoFlavor = muons[1]->q*13;
        leptonTwoDZ     = muons[1]->dz;
        leptonTwoD0     = muons[1]->d0;
        
        // fill dilepton vertex information
        dileptonVertex.SetXYZ(0., 0., 0.);
        dileptonVertexErr.SetXYZ(-1., -1., -1.);
        dileptonVertexChi2 = -1.;
        dileptonVertexDOF = -1.;
        dileptonVertexProb = -1.;
       
        for (int i = 0; i < fDimuonVertexArr->GetEntries(); ++i) {
            TVertex* dimuonVert = (TVertex*) fDimuonVertexArr->At(i);
            if (
                    ((dimuonVert->index1 == muons[0]->muIndex && dimuonVert->index2 == muons[1]->muIndex) ||
                    (dimuonVert->index1 == muons[1]->muIndex && dimuonVert->index2 == muons[0]->muIndex)) &&
                    (dimuonVert->isValid)
               ) {
                dileptonVertex.SetXYZ(dimuonVert->x, dimuonVert->y, dimuonVert->z);
                dileptonVertexErr.SetXYZ(dimuonVert->xerr, dimuonVert->yerr, dimuonVert->zerr);
                dileptonVertexChi2 = dimuonVert->chi2;
                dileptonVertexDOF = dimuonVert->ndof;
                dileptonVertexProb = dimuonVert->prob;
                break;
            }
        } 

        photonOneMVA = photons[0]->mva;
        passElectronVeto = photons[0]->passElectronVeto;  

        if (!isData)
            photonOneR9 = weights->GetCorrectedPhotonR9(*photons[0]);
        else 
            photonOneR9 = photons[0]->r9;


        if (!isData) {

            eventWeight *= weights->GetHZZMuonIDEff(*muons[0]); 
            eventWeight *= weights->GetHZZMuonIDEff(*muons[1]);
            eventWeight *= weights->GetMuonISOEff(leptonOneP4);
            eventWeight *= weights->GetMuonISOEff(leptonTwoP4);

            pair<float, float> eff11 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg1", *muons[0]);
            pair<float, float> eff12 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg1", *muons[1]);
            pair<float, float> eff21 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg2", *muons[0]);
            pair<float, float> eff22 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg2", *muons[1]);
            float eff_data = eff11.first*eff22.first + eff12.first*eff21.first - eff11.first*eff12.first;
            float eff_mc = eff11.second*eff22.second + eff12.second*eff21.second - eff11.second*eff12.second;
            triggerWeight = eff_data/eff_mc;
            eventWeight *= triggerWeight;

            eventWeight *= weights->GetPhotonMVAIdEff(*photons[0]);
        }

    } // end mumug selection
    
    else if (params->selection == "ee") {
        if (electrons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        leptonOneP4.SetPtEtaPhiM(electrons[0]->calibPt, electrons[0]->eta, electrons[0]->phi, ELE_MASS);
        leptonTwoP4.SetPtEtaPhiM(electrons[1]->calibPt, electrons[1]->eta, electrons[1]->phi, ELE_MASS);

        if (leptonOneP4.Pt() < 25.0) 
            return kTRUE;
        hTotalEvents->Fill(6);

        if (leptonTwoP4.Pt() < 15.0)
            return kTRUE;
        hTotalEvents->Fill(7);

        if (electrons[0]->q == electrons[1]->q)
            return kTRUE;
        hTotalEvents->Fill(8);
        
        leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonOneFlavor = electrons[0]->q*11;
        leptonOneDZ     = electrons[0]->dz;
        leptonOneD0     = electrons[0]->d0;
            
        leptonTwoIso    = GetElectronIsolation(electrons[1], fInfo->rhoJet);
        leptonTwoFlavor = electrons[1]->q*11;
        leptonTwoDZ     = electrons[1]->dz;
        leptonTwoD0     = electrons[1]->d0;
        
        // fill dilepton vertex information
        dileptonVertex.SetXYZ(0., 0., 0.);
        dileptonVertexErr.SetXYZ(-1., -1., -1.);
        dileptonVertexChi2 = -1.;
        dileptonVertexDOF = -1.;
        dileptonVertexProb = -1.;
       
        for (int i = 0; i < fDielectronVertexArr->GetEntries(); ++i) {
            TVertex* dielectronVert = (TVertex*) fDielectronVertexArr->At(i);
            if (
                    ((dielectronVert->index1 == electrons[0]->eleIndex && dielectronVert->index2 == electrons[1]->eleIndex) ||
                    (dielectronVert->index1 == electrons[1]->eleIndex && dielectronVert->index2 == electrons[0]->eleIndex)) &&
                    (dielectronVert->isValid)
               ) {
                dileptonVertex.SetXYZ(dielectronVert->x, dielectronVert->y, dielectronVert->z);
                dileptonVertexErr.SetXYZ(dielectronVert->xerr, dielectronVert->yerr, dielectronVert->zerr);
                dileptonVertexChi2 = dielectronVert->chi2;
                dileptonVertexDOF = dielectronVert->ndof;
                dileptonVertexProb = dielectronVert->prob;
                break;
            }
        } 
           
        if (!isData) {

            eventWeight *= weights->GetHZZElectronRecoIdEff(*electrons[0]); 
            eventWeight *= weights->GetHZZElectronRecoIdEff(*electrons[1]); 

            pair<float, float> eff11 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg1", *electrons[0]);
            pair<float, float> eff12 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg1", *electrons[1]);
            pair<float, float> eff21 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg2", *electrons[0]);
            pair<float, float> eff22 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg2", *electrons[1]);
            float eff_data = eff11.first*eff22.first + eff12.first*eff21.first - eff11.first*eff12.first;
            float eff_mc = eff11.second*eff22.second + eff12.second*eff21.second - eff11.second*eff12.second;
            triggerWeight = eff_data/eff_mc;
            eventWeight *= triggerWeight;

        }

    } // end ee selection
    
    else if (params->selection == "elelg") {
        if (electrons.size() < 2) 
            return kTRUE;
        hTotalEvents->Fill(5);

        if (photons.size() < 1)
            return kTRUE;
        hTotalEvents->Fill(6);
        
        leptonOneP4.SetPtEtaPhiM(electrons[0]->calibPt, electrons[0]->eta, electrons[0]->phi, ELE_MASS);
        leptonTwoP4.SetPtEtaPhiM(electrons[1]->calibPt, electrons[1]->eta, electrons[1]->phi, ELE_MASS);

        if (leptonOneP4.Pt() < 25.0) 
            return kTRUE;
        hTotalEvents->Fill(6);

        if (leptonTwoP4.Pt() < 15.0)
            return kTRUE;
        hTotalEvents->Fill(7);

        if (electrons[0]->q == electrons[1]->q)
            return kTRUE;
        hTotalEvents->Fill(8);
       
        photonOneP4.SetPtEtaPhiM(photons[0]->calibPt, photons[0]->eta, photons[0]->phi, 0.);
        if (photonOneP4.Pt() < 15.0)
            return kTRUE;
        hTotalEvents->Fill(9);
         
        leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonOneFlavor = electrons[0]->q*11;
        leptonOneDZ     = electrons[0]->dz;
        leptonOneD0     = electrons[0]->d0;
            
        leptonTwoIso    = GetElectronIsolation(electrons[1], fInfo->rhoJet);
        leptonTwoFlavor = electrons[1]->q*11;
        leptonTwoDZ     = electrons[1]->dz;
        leptonTwoD0     = electrons[1]->d0;
        
        // fill dilepton vertex information
        dileptonVertex.SetXYZ(0., 0., 0.);
        dileptonVertexErr.SetXYZ(-1., -1., -1.);
        dileptonVertexChi2 = -1.;
        dileptonVertexDOF = -1.;
        dileptonVertexProb = -1.;
       
        for (int i = 0; i < fDielectronVertexArr->GetEntries(); ++i) {
            TVertex* dielectronVert = (TVertex*) fDielectronVertexArr->At(i);
            if (
                    ((dielectronVert->index1 == electrons[0]->eleIndex && dielectronVert->index2 == electrons[1]->eleIndex) ||
                    (dielectronVert->index1 == electrons[1]->eleIndex && dielectronVert->index2 == electrons[0]->eleIndex)) &&
                    (dielectronVert->isValid)
               ) {
                dileptonVertex.SetXYZ(dielectronVert->x, dielectronVert->y, dielectronVert->z);
                dileptonVertexErr.SetXYZ(dielectronVert->xerr, dielectronVert->yerr, dielectronVert->zerr);
                dileptonVertexChi2 = dielectronVert->chi2;
                dileptonVertexDOF = dielectronVert->ndof;
                dileptonVertexProb = dielectronVert->prob;
                break;
            }
        } 

        photonOneMVA = photons[0]->mva;
        passElectronVeto = photons[0]->passElectronVeto;  

        if (!isData)
            photonOneR9 = weights->GetCorrectedPhotonR9(*photons[0]);
        else 
            photonOneR9 = photons[0]->r9;

        if (!isData) {

            eventWeight *= weights->GetHZZElectronRecoIdEff(*electrons[0]); 
            eventWeight *= weights->GetHZZElectronRecoIdEff(*electrons[1]); 

            pair<float, float> eff11 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg1", *electrons[0]);
            pair<float, float> eff12 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg1", *electrons[1]);
            pair<float, float> eff21 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg2", *electrons[0]);
            pair<float, float> eff22 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg2", *electrons[1]);
            float eff_data = eff11.first*eff22.first + eff12.first*eff21.first - eff11.first*eff12.first;
            float eff_mc = eff11.second*eff22.second + eff12.second*eff21.second - eff11.second*eff12.second;
            triggerWeight = eff_data/eff_mc;
            eventWeight *= triggerWeight;
            
            eventWeight *= weights->GetPhotonMVAIdEff(*photons[0]);
        } 

    } // end elelg selection

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
 
    if (!isData && genLeptons.size() == 2) {
        genLeptonOneId = genLeptons[0]->pdgId;
        genLeptonOneP4.SetPtEtaPhiM(genLeptons[0]->pt, genLeptons[0]->eta, genLeptons[0]->phi, genLeptons[0]->mass);
        genLeptonTwoId = genLeptons[1]->pdgId;
        genLeptonTwoP4.SetPtEtaPhiM(genLeptons[1]->pt, genLeptons[1]->eta, genLeptons[1]->phi, genLeptons[1]->mass);
    }
    else {
        genLeptonOneP4.SetPtEtaPhiM(0., 0., 0., 0.);
        genLeptonTwoP4.SetPtEtaPhiM(0., 0., 0., 0.);
    }

    if (!isData && genPhotons.size() == 1) 
        genPhotonP4.SetPtEtaPhiM(genPhotons[0]->pt, genPhotons[0]->eta, genPhotons[0]->phi, 0.);
    else
        genPhotonP4.SetPtEtaPhiM(0., 0., 0., 0.);

    outTree->Fill();
    this->passedEvents++;

    return kTRUE;
}

void upsilonGammaAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void upsilonGammaAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void upsilonGammaAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<upsilonGammaAnalyzer> selector(new upsilonGammaAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

float upsilonGammaAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    //float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    float combIso = (mu->chHadIso03 + std::max(0.,(double)mu->neuHadIso03 + mu->gammaIso03 - 0.5*mu->puIso03));
    return combIso;
}

float upsilonGammaAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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

float upsilonGammaAnalyzer::GetPhotonIsolation(const baconhep::TPhoton* pho, const float rho)
{
    int iEta = 0;
    float etaBins[8] = {0., 1., 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};
    float effArea[8] = {0.1703, 0.1715, 0.1213, 0.1230, 0.1635, 0.1937, 0.2393};
    for (unsigned i = 0; i < 8; ++i) {
        if (fabs(pho->scEta) > etaBins[i] && fabs(pho->scEta) < etaBins[i+1]) {
            iEta = i;
            break;
        }
    }

    float combIso = pho->chHadIso + std::max(0., (double)pho->neuHadIso + pho->gammaIso - rho*effArea[iEta]);

    return combIso;
}
