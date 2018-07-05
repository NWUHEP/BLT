#include "hzgAnalyzer.h"
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

hzgAnalyzer::hzgAnalyzer() : BLTSelector()
{

}

hzgAnalyzer::~hzgAnalyzer()
{

}

void hzgAnalyzer::Begin(TTree *tree)
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
        triggerNames.push_back("HLT_IsoMu27_v*"); // unprescaled 2017 single muon trigger
    }
    else if (params->selection == "ee" || params->selection == "elelg") {
        triggerNames.push_back("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v*"); // unprescaled 2017 single electron trigger
    }
        
    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask

    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt"; // should verify this is the right one!
    lumiMask.AddJSONFile(jsonFileName);

    // muon momentum corrections
    muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/RoccoR2017v0.txt");
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
    outTree->Branch("metNC", &metNC);
    outTree->Branch("metPhiNC", &metPhiNC);
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

    outTree->Branch("isLeptonTag", &isLeptonTag);

    // photons
    outTree->Branch("photonOneP4", &photonOneP4);
    outTree->Branch("photonOneR9", &photonOneR9);
    outTree->Branch("photonOneMVA", &photonOneMVA);
    outTree->Branch("passElectronVeto", &passElectronVeto);

    // fsr photons
    /*outTree->Branch("leptonOneFSRPhotons", &leptonOneFSRPhotons);
    outTree->Branch("leptonOneFSRSum", &leptonOneFSRSum);
    outTree->Branch("leptonOneFSRMatchP4", &leptonOneFSRMatchP4);
    outTree->Branch("leptonOneFSRIsoSum", &leptonOneFSRIsoSum);
    outTree->Branch("fsrPhotonOneMVA", &fsrPhotonOneMVA);
    outTree->Branch("fsrPhotonOneR9", &fsrPhotonOneR9);
    outTree->Branch("fsrPhotonOneIso", &fsrPhotonOneIso);
    outTree->Branch("fsrPassElectronVetoOne", &fsrPassElectronVetoOne);
    outTree->Branch("genFSRPhotonOneP4", &genFSRPhotonOneP4);
    outTree->Branch("genTrueFSRPhoOneP4", &genTrueFSRPhoOneP4);
    outTree->Branch("genFSRPhotonOneFHPFS", &genFSRPhotonOneFHPFS);
    outTree->Branch("genFSRPhotonOneIPFS", &genFSRPhotonOneIPFS);
    outTree->Branch("leptonOneHasFSRPhoton", &leptonOneHasFSRPhoton);
    outTree->Branch("leptonOneHasRecoveredFSRPhoton", &leptonOneHasRecoveredFSRPhoton);
    outTree->Branch("leptonOneHasFakeFSRPhoton", &leptonOneHasFakeFSRPhoton);
    
    outTree->Branch("leptonTwoFSRPhotons", &leptonTwoFSRPhotons);
    outTree->Branch("leptonTwoFSRSum", &leptonTwoFSRSum);
    outTree->Branch("leptonTwoFSRMatchP4", &leptonTwoFSRMatchP4);
    outTree->Branch("leptonTwoFSRIsoSum", &leptonTwoFSRIsoSum);
    outTree->Branch("fsrPhotonTwoMVA", &fsrPhotonTwoMVA);
    outTree->Branch("fsrPhotonTwoR9", &fsrPhotonTwoR9);
    outTree->Branch("fsrPhotonTwoIso", &fsrPhotonTwoIso);
    outTree->Branch("fsrPassElectronVetoTwo", &fsrPassElectronVetoTwo);
    outTree->Branch("genFSRPhotonTwoP4", &genFSRPhotonTwoP4);
    outTree->Branch("genTrueFSRPhoTwoP4", &genTrueFSRPhoTwoP4);
    outTree->Branch("genFSRPhotonTwoFHPFS", &genFSRPhotonTwoFHPFS);
    outTree->Branch("genFSRPhotonTwoIPFS", &genFSRPhotonTwoIPFS);
    outTree->Branch("leptonTwoHasFSRPhoton", &leptonTwoHasFSRPhoton);
    outTree->Branch("leptonTwoHasRecoveredFSRPhoton", &leptonTwoHasRecoveredFSRPhoton);
    outTree->Branch("leptonTwoHasFakeFSRPhoton", &leptonTwoHasFakeFSRPhoton); */

    // lepton vertices
    outTree->Branch("dileptonVertexOne", &dileptonVertexOne);
    outTree->Branch("dileptonVertexErrOne", &dileptonVertexErrOne);
    outTree->Branch("dileptonVertexChi2One", &dileptonVertexChi2One);
    outTree->Branch("dileptonVertexDOFOne", &dileptonVertexDOFOne);

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
    outTree->Branch("vetoDY", &vetoDY);
    //outTree->Branch("brianVetoDY", &brianVetoDY);

    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nTaus", &nTaus);
    outTree->Branch("nPhotons", &nPhotons);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nFwdJets", &nFwdJets);
    outTree->Branch("nBJets", &nBJets);

    // event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",30,0.5,30.5);

    ReportPostBegin();
}

Bool_t hzgAnalyzer::Process(Long64_t entry)
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
   
    bool sync_print = false;
    bool sync_print_precut = false;

    if (sync_print_precut) {          
        if (!(
                (fInfo->runNum == 275657 && fInfo->lumiSec == 45 && fInfo->evtNum == 89361789) ||
                (fInfo->runNum == 275657 && fInfo->lumiSec == 45 && fInfo->evtNum == 88720008) ||
                (fInfo->runNum == 275657 && fInfo->lumiSec == 45 && fInfo->evtNum == 87753250) ||
                (fInfo->runNum == 275657 && fInfo->lumiSec == 45 && fInfo->evtNum == 89277803) ||
                (fInfo->runNum == 275657 && fInfo->lumiSec == 45 && fInfo->evtNum == 88609823) ||
                (fInfo->runNum == 275657 && fInfo->lumiSec == 45 && fInfo->evtNum == 87863789) ||
                (fInfo->runNum == 275657 && fInfo->lumiSec == 45 && fInfo->evtNum == 88896784) ||
                (fInfo->runNum == 275657 && fInfo->lumiSec == 45 && fInfo->evtNum == 88885728) ||
                (fInfo->runNum == 275657 && fInfo->lumiSec == 45 && fInfo->evtNum == 87755395) ||
                (fInfo->runNum == 275657 && fInfo->lumiSec == 45 && fInfo->evtNum == 87761828))
            ) return kTRUE;

        cout << "run, lumi, event" << endl;
        cout << fInfo->runNum << ", " << fInfo->lumiSec << ", " << fInfo->evtNum << endl;
    }
          
    ///////////////////////
    // Generator objects // 
    ///////////////////////

    //brianVetoDY = false; 
    //vector<TGenParticle*> genParticles;
    vector<TGenParticle*> genLeptons;
    vector<TGenParticle*> genPhotons;
    //vector<TGenParticle*> genFSRPhotons;
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
            //if ((abs(particle->pdgId) == 11 || abs(particle->pdgId) == 13) && particle->isPromptFinalState) 
            //    genLeptons.push_back(particle);
            //
            if ((fabs(particle->pdgId) == 11 || fabs(particle->pdgId) == 13) and particle->parent > 0) {
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                if (fabs(mother->pdgId) == 23) 
                    genLeptons.push_back(particle);
            }
                

            // saving photons
            //if (abs(particle->pdgId) == 22 && (particle->isPromptFinalState || particle->fromHardProcessFinalState))
            //    genPhotons.push_back(particle);
            //
            
            //if (abs(particle->pdgId) == 22 && particle->parent > 0) {
            //    TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
            //    while (abs(mother->pdgId) == 22) {
            //        particle = mother;
            //        if (particle->parent > 0)
            //            mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
            //        else
            //            break;
            //    }
            //    genPhotons.push_back(particle);
            //}
            
            // Brian's DY photon veto
            //if (fabs(particle->pdgId) == 22 && particle->parent > 0) {
            //    TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
            //    if (fabs(mother->pdgId) < 22 && particle->pt > 5) {
            //        TLorentzVector particleP4;
            //        TLorentzVector motherP4;
            //        particleP4.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, particle->mass); 
            //        motherP4.SetPtEtaPhiM(mother->pt, mother->eta, mother->phi, mother->mass); 
            //        if (particleP4.DeltaR(motherP4) > 0.3) 
            //            brianVetoDY = true;
            //    }                    
            //}
           
           
           if (fabs(particle->pdgId) == 22) {
               genPhotons.push_back(particle);
               //TGenParticle* phoMom = (TGenParticle*) fGenParticleArr->At(particle->parent);
               //if (fabs(phoMom->pdgId) == 13)
               //    genFSRPhotons.push_back(particle);
           }
 
                

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
    //if (!passTrigger && isData)
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
    TVertex* thePV;
    if (fInfo->hasGoodPV) {
        assert(fPVArr->GetEntries() != 0);
        TVector3 pv;
        thePV = (TVertex *)fPVArr->At(0);
        copy_xyz((TVertex*) fPVArr->At(0), pv);
        rPV = pv; 
        particleSelector->SetPV(pv);
    } else {
        return kTRUE;
    }
    hTotalEvents->Fill(4);
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);

    if (sync_print_precut) {
        cout << "pvx, pvy, pvz, ndof" << endl;
        cout << thePV->x << ", " << thePV->y << ", " << thePV->z << ", " << thePV->ndof << endl;
    }

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

        //if (sync_print_precut) {
        //    cout << "pt_before_roccor, pt_after_roccor" << endl;
        //    cout << muon->pt << ", ";
        //}

        muon->pt = muonSF*muon->pt; 
        //if (sync_print_precut)
        //    cout << muon->pt << endl;
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

        if (sync_print_precut) {
           cout << "pt , eta, phi, isGlobal, isTracker, nMatchStn, nValidHits, d0, dz, sip3d, nPixHits, nTkLayers, pt_track, ptErr, iso, trkIso" << endl;
           cout << muonP4.Pt() << "," << muonP4.Eta() << "," << 
                   muonP4.Phi() << "," <<
                   (muon->typeBits & baconhep::kGlobal) << "," << 
                   (muon->typeBits & baconhep::kTracker) << "," <<  
                   muon->nMatchStn << "," << muon->nValidHits << ", " << muon->d0 << "," << 
                   muon->dz << "," << muon->sip3d << "," << 
                   muon->nPixHits << "," << muon->nTkLayers << "," << 
                   muon->pt << "," << muon->ptErr << "," << 
                   GetMuonIsolation(muon) << ", " << muon->trkIso << endl;
        }

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
        electronP4.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, ELE_MASS);
        
        if (sync_print_precut) {
           cout << "pt , eta, sc_eta, phi, d0, dz, sip3d, iso, pass_mva" << endl;
           cout << electronP4.Pt() << "," << electronP4.Eta() << "," << 
                   electron->scEta << ", " << electronP4.Phi() << "," << electron->d0 << "," << 
                   electron->dz << "," << electron->sip3d << "," << 
                   GetElectronIsolation(electron, fInfo->rhoJet) << ", " << 
                   particleSelector->PassElectronMVA(electron, cuts->HZZMVAElIDNoIso) <<endl;
        }

        if (
                electron->pt > 7
                && fabs(electron->scEta) < 2.5
                //&& particleSelector->PassElectronID(electron, cuts->tightElID)
                //&& particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
                && particleSelector->PassElectronMVA(electron, cuts->HZZMVAElIDNoIso)
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
    //vector <TPhoton*> fsr_photons;
    vector<TLorentzVector> veto_photons;
    for (int i=0; i<fPhotonArr->GetEntries(); i++) {
        TPhoton* photon = (TPhoton*) fPhotonArr->At(i);
        assert(photon);
        
        TLorentzVector photonP4;
        photonP4.SetPtEtaPhiM(photon->pt, photon->eta, photon->phi, 0.);
    
        if (sync_print_precut) {
            cout << "photon_pt, photon_eta, photon_sc_eta, photon_phi, photon_mva, pass_electron_veto" << endl;
            cout << photon->pt << ", " << photon->eta << ", " << photon->scEta << ", " << photon->phi << ", " << photon->mva 
                 << ", " << photon->passElectronVeto << endl;
        }

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

    // FSR photons
    /*if (
            photon->pt > 2 
            && fabs(photon->scEta) < 2.4
            && GetPhotonIsolation(photon, fInfo->rhoJet) < 1.8
       ) {
        fsr_photons.push_back(photon);
    }*/

    } 
    sort(photons.begin(), photons.end(), sort_by_higher_pt<TPhoton>);
    //sort(fsr_photons.begin(), fsr_photons.end(), sort_by_higher_pt<TPhoton>);

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

        /*if (isData) { // fix for broken bacon JEC
            double jec = particleSelector->JetCorrector(jet, "NONE");
            jet->pt = jet->ptRaw*jec;
        }*/

        // Prevent overlap of muons and jets
        TLorentzVector jetP4; 
        jetP4.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
        if (sync_print_precut) {
            std::cout << "jet pt, jet eta, jet phi" << std::endl;
            std::cout << jetP4.Pt() << ", " << jetP4.Eta() << ", " << jetP4.Phi() << std::endl;
        }
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
                && particleSelector->PassJetID(jet, cuts->tightJetID)
                && !muOverlap 
                && !elOverlap
                && !phoOverlap
           ) {
            
            jets.push_back(jet);

            if (fabs(jet->eta) <= 2.4) { 
                hadronicP4 += jetP4;
                sumJetPt += jetP4.Pt();

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
    metNC  = fInfo->pfMET;
    metPhiNC = fInfo->pfMETphi;

    if (sync_print_precut) {
        std::cout << "met, metPhi, metNC, metPhiNC" << std::endl;
        std::cout << met << ", " << metPhi << ", " << metNC << ", " << metPhiNC << std::endl;
    }

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
        //cout << "passed muon number cut" << endl;
        hTotalEvents->Fill(5);
        TLorentzVector muonOneP4, muonTwoP4;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, MUON_MASS);
        muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, MUON_MASS);
        if (muonOneP4.Pt() < 20.0) 
            return kTRUE;
        //cout << "passed pt1 cut" << endl;
        hTotalEvents->Fill(6);
        if (muonTwoP4.Pt() < 10.0)
            return kTRUE;
        //cout << "passed pt2 cut" << endl;
        hTotalEvents->Fill(7);
        TLorentzVector dimuon = muonOneP4 + muonTwoP4;
        if (muons[0]->q == muons[1]->q)
            return kTRUE;
        //cout << "passed charge requirement" << endl;
        hTotalEvents->Fill(8);
        if (dimuon.M() < 80.0 || dimuon.M() > 100.0)
            return kTRUE;
        //cout << "passed Z window cut" << endl;
        hTotalEvents->Fill(9);

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

        if (!isData) {
            eventWeight *= weights->GetHZZMuonIDEff(*muons[0]); 
            eventWeight *= weights->GetHZZMuonIDEff(*muons[1]);
        }


        int isGlobalOne = (muons[0]->typeBits & baconhep::kGlobal) > 0;
        int isGlobalTwo = (muons[1]->typeBits & baconhep::kGlobal) > 0;
        int isTrackerOne = (muons[0]->typeBits & baconhep::kTracker) > 0;
        int isTrackerTwo = (muons[1]->typeBits & baconhep::kTracker) > 0;

        if (sync_print)
            cout << runNumber << "," << lumiSection << "," << evtNumber << ","  << 
                    thePV->x << ", " << thePV->y << ", " << thePV->z << ", " << 
                    muonOneP4.Pt() << "," << muonOneP4.Eta() << "," << 
                    muonOneP4.Phi() << "," <<
                    isGlobalOne << "," << isTrackerOne << "," <<  
                    muons[0]->nMatchStn << "," << muons[0]->d0 << "," << 
                    muons[0]->dz << "," << muons[0]->sip3d << "," << 
                    muons[0]->nPixHits << "," << muons[0]->nTkLayers << "," << 
                    muons[0]->pt << "," << muons[0]->ptErr << "," << 
                    GetMuonIsolation(muons[0]) << "," <<
                    muonTwoP4.Pt() << "," << muonTwoP4.Eta() << "," << 
                    muonTwoP4.Phi() << "," << 
                    isGlobalTwo << "," << isTrackerTwo << "," <<
                    muons[1]->nMatchStn << "," << muons[1]->d0 << "," << 
                    muons[1]->dz << "," << muons[1]->sip3d << "," << 
                    muons[1]->nPixHits << "," << muons[1]->nTkLayers << "," << 
                    muons[1]->pt << "," << muons[1]->ptErr << "," << 
                    GetMuonIsolation(muons[1]) << endl;
    }

    else if (params->selection == "mumug") {
        if (muons.size() < 2) 
            return kTRUE;
        hTotalEvents->Fill(5);
        if (photons.size() < 1)
            return kTRUE;
        hTotalEvents->Fill(6);
        
        unsigned int muonOneIndex = 0;
        unsigned int muonTwoIndex = 1;
        bool hasValidPair = false;
        float zMassDiff = 100.;
        for (unsigned int i = 0; i < muons.size(); ++i) {
            for (unsigned int j = i+1; j < muons.size(); ++j) {
                TLorentzVector tempMuonOne, tempMuonTwo;
                tempMuonOne.SetPtEtaPhiM(muons[i]->pt, muons[i]->eta, muons[i]->phi, MUON_MASS);
                tempMuonTwo.SetPtEtaPhiM(muons[j]->pt, muons[j]->eta, muons[j]->phi, MUON_MASS);
                float thisMass = (tempMuonOne + tempMuonTwo).M();
                if (
                        muons[i]->q != muons[j]->q
                        //muons[i]->q == muons[j]->q // fake selection
                        && muons[i]->pt > 20.0 
                        && muons[j]->pt > 10.0
                        && thisMass > 50.0
                   ) {
                    if (hasValidPair) {
                        if (fabs(thisMass - ZMASS) < zMassDiff) {
                            zMassDiff = fabs(thisMass - ZMASS);
                            leptonOneP4 = tempMuonOne;
                            leptonTwoP4 = tempMuonTwo;
                            muonOneIndex = i;
                            muonTwoIndex = j;
                        }
                    }
                    else {
                        zMassDiff = fabs(thisMass - ZMASS);
                        leptonOneP4 = tempMuonOne;
                        leptonTwoP4 = tempMuonTwo;
                        muonOneIndex = i;
                        muonTwoIndex = j;
                        hasValidPair = true;
                    }
                }
            }
        }


        if (!hasValidPair)
            return kTRUE;
        hTotalEvents->Fill(7);

        // trigger matching:
        bool mu1_fired_leg1 = false;
        bool mu1_fired_leg2 = false;
        bool mu2_fired_leg1 = false;
        bool mu2_fired_leg2 = false;
        for (unsigned int iT = 0; iT < triggerNames.size(); ++iT) {
            mu1_fired_leg1 |= trigger->passObj(triggerNames.at(iT), 1, muons[muonOneIndex]->hltMatchBits);
            mu1_fired_leg2 |= trigger->passObj(triggerNames.at(iT), 2, muons[muonOneIndex]->hltMatchBits);
            mu2_fired_leg1 |= trigger->passObj(triggerNames.at(iT), 1, muons[muonTwoIndex]->hltMatchBits);
            mu2_fired_leg2 |= trigger->passObj(triggerNames.at(iT), 2, muons[muonTwoIndex]->hltMatchBits);
        }
        if (sync_print_precut) {
            cout << "mu1_fired_leg1, mu1_fired_leg2, mu2_fired_leg1, mu2_fired_leg2" << endl;
            cout << mu1_fired_leg1 << ", " << mu1_fired_leg2 << ", " << mu2_fired_leg1 << ", " << mu2_fired_leg2 << endl;
        }

        // L1EMTF cut 
        if (
            fabs(leptonOneP4.DeltaPhi(leptonTwoP4)) < 70.0*(M_PI/180.0)
            && fabs(leptonOneP4.Eta()) > 1.2 
            && fabs(leptonTwoP4.Eta()) > 1.2
            && leptonOneP4.Eta()*leptonTwoP4.Eta() > 0
           )
            return kTRUE;
        hTotalEvents->Fill(8); 

        bool hasValidPhoton = false;
        unsigned int photonIndex = 0;

        for (unsigned int i = 0; i < photons.size(); ++i) {
            TLorentzVector tempPhoton;
            TLorentzVector tempDilepton;
            TLorentzVector tempLLG;
            tempPhoton.SetPtEtaPhiM(photons[i]->pt, photons[i]->eta, photons[i]->phi, 0.);
            //tempPhoton.SetPtEtaPhiM(photons[i]->pt, photons[i]->eta, photons[i]->phi, 0.);
            tempDilepton = leptonOneP4 + leptonTwoP4;
            tempLLG = leptonOneP4 + leptonTwoP4 + tempPhoton;
            float this_dr1 = leptonOneP4.DeltaR(tempPhoton);
            float this_dr2 = leptonTwoP4.DeltaR(tempPhoton);
            if (
                tempPhoton.Pt() > 15.0 &&
                tempPhoton.Et()/tempLLG.M() > (15.0/110.0) &&
                tempDilepton.M() + tempLLG.M() > 185.0 &&
                tempLLG.M() > 100. && tempLLG.M() < 180. &&
                this_dr1 > 0.4 && this_dr2 > 0.4
                ) {
                hasValidPhoton = true;
                photonIndex = i;
                break;
            }
            //cout << "photon " << i << " valid: " << hasValidPhoton << endl;
        }

        if (!hasValidPhoton)
            return kTRUE;
        hTotalEvents->Fill(9);

        photonOneP4.SetPtEtaPhiM(photons[photonIndex]->pt, photons[photonIndex]->eta, photons[photonIndex]->phi, 0.);
        //photonOneP4.SetPtEtaPhiM(photons[photonIndex]->pt, photons[photonIndex]->eta, photons[photonIndex]->phi, 0.);
        if (photonOneP4.Pt() < 15.0)
            return kTRUE;
        hTotalEvents->Fill(10);

        // DY photon overlap removal
        vetoDY = false;
        for (unsigned int i = 0; i < genPhotons.size(); ++i) {
            TGenParticle *pho = genPhotons.at(i);
            if (pho->fromHardProcessFinalState || pho->isPromptFinalState) {
                TLorentzVector thisGenPhotonP4;
                thisGenPhotonP4.SetPtEtaPhiM(pho->pt, pho->eta, pho->phi, 0.);
                if (thisGenPhotonP4.DeltaR(photonOneP4) < 0.1) {
                    vetoDY = true;
                    break;
                }
            }
        }
        
        // checking for lepton tag
        isLeptonTag = false;
        for (unsigned int i = 0; i < muons.size(); ++i) {
            TLorentzVector tempMuon;
            tempMuon.SetPtEtaPhiM(muons[i]->pt, muons[i]->eta, muons[i]->phi, MUON_MASS);
            if (leptonOneP4.DeltaR(tempMuon) < 0.4 ||
                leptonTwoP4.DeltaR(tempMuon) < 0.4 || 
                photonOneP4.DeltaR(tempMuon) < 0.4)
                continue;
            else
                isLeptonTag = true;
        }

        if (!isLeptonTag) {
            for (unsigned int i = 0; i < electrons.size(); ++i) {
                TLorentzVector tempElectron;
                tempElectron.SetPtEtaPhiM(electrons[i]->pt, electrons[i]->eta, electrons[i]->phi, ELE_MASS);
                if (leptonOneP4.DeltaR(tempElectron) < 0.4 ||
                    leptonTwoP4.DeltaR(tempElectron) < 0.4 || 
                    photonOneP4.DeltaR(tempElectron) < 0.4)
                    continue;
                else
                    isLeptonTag = true;
            }
        }

        //cout << "total number of leptons in the event: " << nMuons + nElectrons << endl;
        //cout << "isLeptonTag: " << isLeptonTag << endl;

        // FSR photon recovery
        /*float minDROne = 100.;
        float minDRTwo = 100.;
        int fsrPhoOneIndex = -1;
        int fsrPhoTwoIndex = -1;
        
        leptonOneFSRPhotons = 0;
        leptonTwoFSRPhotons = 0;
        leptonOneFSRSum.SetPtEtaPhiM(0., 0., 0., 0.);
        leptonTwoFSRSum.SetPtEtaPhiM(0., 0., 0., 0.);
        leptonOneFSRMatchP4.SetPtEtaPhiM(0., 0., 0., 0.);
        leptonTwoFSRMatchP4.SetPtEtaPhiM(0., 0., 0., 0.);
        leptonOneFSRIsoSum = 0.;
        leptonTwoFSRIsoSum = 0.;

        // check truth information; does either lepton have a true FSR photon?
        leptonOneHasFSRPhoton = false;
        leptonTwoHasFSRPhoton = false;
        leptonOneHasRecoveredFSRPhoton = false;
        leptonTwoHasRecoveredFSRPhoton = false;
        leptonOneHasFakeFSRPhoton = false;
        leptonTwoHasFakeFSRPhoton = false;
        int genFSRPhoOneIndex = -1;
        int genFSRPhoTwoIndex = -1;
        genTrueFSRPhoOneP4.SetPtEtaPhiM(0., 0., 0., 0.);
        genTrueFSRPhoTwoP4.SetPtEtaPhiM(0., 0., 0., 0.);
        if (!isData) {
            for (unsigned int i = 0; i < genFSRPhotons.size(); i++) {
                TGenParticle* thisFSRPho = genFSRPhotons.at(i);
                TGenParticle* thisFSRPhoMom = (TGenParticle*)fGenParticleArr->At(thisFSRPho->parent);
                TLorentzVector genMuonP4;
                genMuonP4.SetPtEtaPhiM(thisFSRPhoMom->pt, thisFSRPhoMom->eta, thisFSRPhoMom->phi, thisFSRPhoMom->mass);
                if (genMuonP4.DeltaR(leptonOneP4) < 0.3) {
                    leptonOneHasFSRPhoton = true;
                    genFSRPhoOneIndex = i;
                    genTrueFSRPhoOneP4.SetPtEtaPhiM(thisFSRPho->pt, thisFSRPho->eta, thisFSRPho->phi, thisFSRPho->mass);
                }
                else {
                    if (genMuonP4.DeltaR(leptonTwoP4) < 0.3) {
                    leptonTwoHasFSRPhoton = true;
                    genFSRPhoTwoIndex = i;
                    genTrueFSRPhoTwoP4.SetPtEtaPhiM(thisFSRPho->pt, thisFSRPho->eta, thisFSRPho->phi, thisFSRPho->mass);
                    }
                }
            }
        }

        cout << "lepton 1 hasFSRPhoton = " << leptonOneHasFSRPhoton << endl;
        cout << "lepton 2 hasFSRPhoton = " << leptonTwoHasFSRPhoton << endl;
        cout << "number of reco fsr candidates = " << fsr_photons.size() << endl;
        // apply FSR photon recovery algorithm at reco level
        for (unsigned int i = 0; i < fsr_photons.size(); ++i) {
            TPhoton *pho = fsr_photons.at(i);
            TLorentzVector thisFSRPhoP4;
            thisFSRPhoP4.SetPtEtaPhiM(pho->pt, pho->eta, pho->phi, 0.);
            if (thisFSRPhoP4.DeltaR(photonOneP4) < 0.01) // reject the selected prompt photon
                continue;

            float phoLepOneDR = thisFSRPhoP4.DeltaR(leptonOneP4);
            float phoLepTwoDR = thisFSRPhoP4.DeltaR(leptonTwoP4);

            cout << "phoLepOneDR, phoLepTwoDR = " << phoLepOneDR << ", " << phoLepTwoDR << ", " << endl;
            cout << "thisFSRPhoP4.Et() squared = " << pow(thisFSRPhoP4.Et(), 2) << endl;

            //lepton 1
            if ((phoLepOneDR/pow(thisFSRPhoP4.Et(), 2)) < 0.012 && phoLepOneDR < 0.5) {
                ++leptonOneFSRPhotons;
                leptonOneFSRSum += thisFSRPhoP4;
                if (phoLepOneDR > 0.01 && phoLepOneDR < 0.4)
                    leptonOneFSRIsoSum += GetPhotonIsolation(pho, fInfo->rhoJet);
                if (phoLepOneDR/pow(thisFSRPhoP4.Et(), 2) < minDROne) {
                    leptonOneFSRMatchP4 = thisFSRPhoP4;
                    fsrPhoOneIndex = i;
                    minDROne = phoLepOneDR/pow(thisFSRPhoP4.Et(), 2);
                }
            }
            
            //lepton 2
            if ((phoLepTwoDR/pow(thisFSRPhoP4.Et(), 2)) < 0.012 && phoLepTwoDR < 0.5) {
                ++leptonTwoFSRPhotons;
                leptonTwoFSRSum += thisFSRPhoP4;
                if (phoLepTwoDR > 0.01 && phoLepTwoDR < 0.4)
                    leptonTwoFSRIsoSum += GetPhotonIsolation(pho, fInfo->rhoJet);
                if (phoLepTwoDR/pow(thisFSRPhoP4.Et(), 2) < minDRTwo) {
                    leptonTwoFSRMatchP4 = thisFSRPhoP4;
                    fsrPhoTwoIndex = i;
                    minDRTwo = phoLepTwoDR/pow(thisFSRPhoP4.Et(), 2);
                }
            }
        } 
     
        genFSRPhotonOneP4.SetPtEtaPhiM(0., 0., 0., 0.);
        genFSRPhotonTwoP4.SetPtEtaPhiM(0., 0., 0., 0.);
        genFSRPhotonOneFHPFS = 0;
        genFSRPhotonTwoFHPFS = 0;
        genFSRPhotonOneIPFS = 0;
        genFSRPhotonTwoIPFS = 0;
        if (fsrPhoOneIndex > 0) {
            fsrPhotonOneMVA = fsr_photons[fsrPhoOneIndex]->mva;
            fsrPassElectronVetoOne = fsr_photons[fsrPhoOneIndex]->passElectronVeto;  
            fsrPhotonOneIso = GetPhotonIsolation(fsr_photons[fsrPhoOneIndex], fInfo->rhoJet);
            if (!isData)
                fsrPhotonOneR9 = weights->GetCorrectedPhotonR9(*fsr_photons[fsrPhoOneIndex]);
            else 
                fsrPhotonOneR9 = fsr_photons[fsrPhoOneIndex]->r9;

            // find the matching gen photon
            float min_fsr_dr_one = 1000.;
            for (unsigned int i = 0; i < genPhotons.size(); i++) {
                TLorentzVector tmpGenFSRPhoOne;
                tmpGenFSRPhoOne.SetPtEtaPhiM(genPhotons[i]->pt, genPhotons[i]->eta, genPhotons[i]->phi, genPhotons[i]->mass);
                float this_fsr_dr_one = tmpGenFSRPhoOne.DeltaR(leptonOneFSRMatchP4);
                if (this_fsr_dr_one < min_fsr_dr_one) {
                    genFSRPhotonOneP4.SetPtEtaPhiM(genPhotons[i]->pt, genPhotons[i]->eta, genPhotons[i]->phi, genPhotons[i]->mass);
                    genFSRPhotonOneFHPFS = genPhotons[i]->fromHardProcessFinalState;
                    genFSRPhotonOneIPFS = genPhotons[i]->isPromptFinalState;
                    min_fsr_dr_one = this_fsr_dr_one;
                    if (min_fsr_dr_one < 0.3) { // we have a match; check if consistent with FSR
                        if (genFSRPhoOneIndex > 0) {                        
                            TGenParticle *genPho = genFSRPhotons.at(genFSRPhoOneIndex);
                            TLorentzVector trueFSRPhoton;
                            trueFSRPhoton.SetPtEtaPhiM(genPho->pt, genPho->eta, genPho->phi, genPho->mass);
                            if (trueFSRPhoton.DeltaR(genFSRPhotonOneP4) < 0.01)
                                leptonOneHasRecoveredFSRPhoton = true;
                            else
                                leptonOneHasFakeFSRPhoton = true;
                            break;
                        }
                    }        
                }
            } 
        }

       
        if (fsrPhoTwoIndex > 0) {
            fsrPhotonTwoMVA = fsr_photons[fsrPhoTwoIndex]->mva;
            fsrPassElectronVetoTwo = fsr_photons[fsrPhoTwoIndex]->passElectronVeto;  
            fsrPhotonTwoIso = GetPhotonIsolation(fsr_photons[fsrPhoTwoIndex], fInfo->rhoJet);
            if (!isData)
                fsrPhotonTwoR9 = weights->GetCorrectedPhotonR9(*fsr_photons[fsrPhoTwoIndex]);
            else 
                fsrPhotonTwoR9 = fsr_photons[fsrPhoTwoIndex]->r9;
            
            // find the matching gen photon
            float min_fsr_dr_two = 1000.;
            for (unsigned int i = 0; i < genPhotons.size(); i++) {
                TLorentzVector tmpGenFSRPhoTwo;
                tmpGenFSRPhoTwo.SetPtEtaPhiM(genPhotons[i]->pt, genPhotons[i]->eta, genPhotons[i]->phi, genPhotons[i]->mass);
                float this_fsr_dr_two = tmpGenFSRPhoTwo.DeltaR(leptonTwoFSRMatchP4);
                if (this_fsr_dr_two < min_fsr_dr_two) {
                    genFSRPhotonTwoP4.SetPtEtaPhiM(genPhotons[i]->pt, genPhotons[i]->eta, genPhotons[i]->phi, genPhotons[i]->mass);
                    genFSRPhotonTwoFHPFS = genPhotons[i]->fromHardProcessFinalState;
                    genFSRPhotonTwoIPFS = genPhotons[i]->isPromptFinalState;
                    min_fsr_dr_two = this_fsr_dr_two;
                    if (min_fsr_dr_two < 0.3) { // we have a match; check if consistent with FSR
                        if (genFSRPhoTwoIndex > 0) {                        
                            TGenParticle *genPho = genFSRPhotons.at(genFSRPhoTwoIndex);
                            TLorentzVector trueFSRPhoton;
                            trueFSRPhoton.SetPtEtaPhiM(genPho->pt, genPho->eta, genPho->phi, genPho->mass);
                            if (trueFSRPhoton.DeltaR(genFSRPhotonTwoP4) < 0.01)
                                leptonTwoHasRecoveredFSRPhoton = true;
                            else 
                                leptonTwoHasFakeFSRPhoton = true;
                            break;
                        }
                    }        
                }
            } 
        }

        //if (fsrPhoOneIndex > 0) 
        //    cout << "lepton 1 hasFSRPhoton, hasRecoveredFSRPhoton, hasFakeFSRPhoton: " << 
        //            leptonOneHasFSRPhoton << ", " << leptonOneHasRecoveredFSRPhoton << ", " << leptonOneHasFakeFSRPhoton << endl;
        //if (fsrPhoTwoIndex > 0)
        //    cout << "lepton 2 hasFSRPhoton, hasRecoveredFSRPhoton, hasFakeFSRPhoton: " << 
        //            leptonTwoHasFSRPhoton << ", " << leptonTwoHasRecoveredFSRPhoton << ", " << leptonTwoHasFakeFSRPhoton << endl; */
        
        if (sync_print_precut) {
            cout << "lepton1_pt, lepton1_eta, lepton1_phi, lepton2_pt, lepton2_eta, lepton2_phi, dilepton_mass, dr1, dr2" << endl;
            cout << leptonOneP4.Pt() << ", " << leptonOneP4.Eta() << ", " << leptonOneP4.Phi() << ", " 
                 << leptonTwoP4.Pt() << ", " << leptonTwoP4.Eta() << ", " << leptonTwoP4.Phi() << ", "
                 << (leptonOneP4 + leptonTwoP4).M() << ", " << leptonOneP4.DeltaR(photonOneP4) << ", " << leptonTwoP4.DeltaR(photonOneP4) << endl;
        }

        //float m_ll = (leptonOneP4 + leptonTwoP4).M();
        //float m_llg = (leptonOneP4 + leptonTwoP4 + photonOneP4).M();
        //if (photonOneP4.Et() < 0.14*m_llg)
        //    return kTRUE;
        //hTotalEvents->Fill(9);
        //if (leptonOneP4.DeltaR(photonOneP4) < 0.4)
        //    return kTRUE;
        //hTotalEvents->Fill(10);
        //if (leptonTwoP4.DeltaR(photonOneP4) < 0.4)
        //    return kTRUE;
        //hTotalEvents->Fill(11);
        //if (m_ll + m_llg < 185.0)
        //    return kTRUE;
        //hTotalEvents->Fill(12);
        //if (m_ll < 100.0 || m_ll > 180.0)
        //    return kTRUE;
        //hTotalEvents->Fill(13);
 
        leptonOneIso    = GetMuonIsolation(muons[muonOneIndex]);
        leptonOneFlavor = muons[muonOneIndex]->q*13;
        leptonOneDZ     = muons[muonOneIndex]->dz;
        leptonOneD0     = muons[muonOneIndex]->d0;
            
        leptonTwoIso    = GetMuonIsolation(muons[muonTwoIndex]);
        leptonTwoFlavor = muons[muonTwoIndex]->q*13;
        leptonTwoDZ     = muons[muonTwoIndex]->dz;
        leptonTwoD0     = muons[muonTwoIndex]->d0;

        photonOneMVA = photons[photonIndex]->mva;
        passElectronVeto = photons[photonIndex]->passElectronVeto;  
        if (!isData)
            photonOneR9 = weights->GetCorrectedPhotonR9(*photons[photonIndex]);
        else 
            photonOneR9 = photons[photonIndex]->r9;

        if (sync_print_precut) {
            cout << "event is still alive" << endl;
        }

        if (!isData) {

            eventWeight *= weights->GetHZZMuonIDEff(*muons[muonOneIndex]); 
            eventWeight *= weights->GetHZZMuonIDEff(*muons[muonTwoIndex]);
            eventWeight *= weights->GetMuonISOEff(leptonOneP4);
            eventWeight *= weights->GetMuonISOEff(leptonTwoP4);

            pair<float, float> eff11 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg1", *muons[muonOneIndex]);
            pair<float, float> eff12 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg1", *muons[muonTwoIndex]);
            pair<float, float> eff21 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg2", *muons[muonOneIndex]);
            pair<float, float> eff22 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg2", *muons[muonTwoIndex]);
            float eff_data = eff11.first*eff22.first + eff12.first*eff21.first - eff11.first*eff12.first;
            float eff_mc = eff11.second*eff22.second + eff12.second*eff21.second - eff11.second*eff12.second;
            triggerWeight = eff_data/eff_mc;
            //cout << "triggerWeight: " << triggerWeight << endl;
            //cout << "event weight before trigger: " << eventWeight << endl;
            eventWeight *= triggerWeight;
            //cout << "event weight after trigger: " << eventWeight << endl;
            
            eventWeight *= weights->GetPhotonMVAIdEff(*photons[photonIndex]);
        
            // photon r9 reweighting
            //if (fabs(photons[photonIndex]->scEta) < 1.444)
            //    photonOneR9 = 1.0045*photonOneR9 + 0.0010;
            //else
            //    photonOneR9 = 1.0086*photonOneR9 - 0.0007; 
        }
        if (sync_print_precut) {
            cout << "event still alive after event weights" << endl;
        }
    } // end mumug selection
    
    else if (params->selection == "ee") {
        if (electrons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);
        TLorentzVector electronOneP4, electronTwoP4;
        electronOneP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, ELE_MASS);
        electronTwoP4.SetPtEtaPhiM(electrons[1]->pt, electrons[1]->eta, electrons[1]->phi, ELE_MASS);
        if (electronOneP4.Pt() < 25.0) 
            return kTRUE;
        hTotalEvents->Fill(6);
        if (electronTwoP4.Pt() < 15.0)
            return kTRUE;
        hTotalEvents->Fill(7);
        TLorentzVector dielectron = electronOneP4 + electronTwoP4;
        if (electrons[0]->q == electrons[1]->q)
            return kTRUE;
        hTotalEvents->Fill(8);
        if (dielectron.M() < 80.0 || dielectron.M() > 100.0)
            return kTRUE;
        hTotalEvents->Fill(9);   
        
        leptonOneP4     = electronOneP4;
        leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonOneFlavor = electrons[0]->q*11;
        leptonOneDZ     = electrons[0]->dz;
        leptonOneD0     = electrons[0]->d0;
            
        leptonTwoP4     = electronTwoP4;
        leptonTwoIso    = GetElectronIsolation(electrons[1], fInfo->rhoJet);
        leptonTwoFlavor = electrons[1]->q*11;
        leptonTwoDZ     = electrons[1]->dz;
        leptonTwoD0     = electrons[1]->d0;
           
        if (!isData) {

            eventWeight *= weights->GetHZZElectronRecoIdEff(*electrons[0]); 
            eventWeight *= weights->GetHZZElectronRecoIdEff(*electrons[1]); 

            //pair<float, float> eff11 = weights->GetTriggerEffWeight("HLT_DoubleEG_leg1", leptonOneP4);
            //pair<float, float> eff12 = weights->GetTriggerEffWeight("HLT_DoubleEG_leg1", leptonTwoP4);
            //pair<float, float> eff21 = weights->GetTriggerEffWeight("HLT_DoubleEG_leg2", leptonOneP4);
            //pair<float, float> eff22 = weights->GetTriggerEffWeight("HLT_DoubleEG_leg2", leptonTwoP4);
            pair<float, float> eff11 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg1", *electrons[0]);
            pair<float, float> eff12 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg1", *electrons[1]);
            pair<float, float> eff21 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg2", *electrons[0]);
            pair<float, float> eff22 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg2", *electrons[1]);
            float eff_data = eff11.first*eff22.first + eff12.first*eff21.first - eff11.first*eff12.first;
            float eff_mc = eff11.second*eff22.second + eff12.second*eff21.second - eff11.second*eff12.second;
            triggerWeight = eff_data/eff_mc;
            //eventWeight *= triggerWeight;

        }
    }
    
    else if (params->selection == "elelg") {
        if (electrons.size() < 2) 
            return kTRUE;
        hTotalEvents->Fill(5);
        if (sync_print_precut)
            cout << "at least two electrons present" << endl;
        if (photons.size() < 1)
            return kTRUE;
        hTotalEvents->Fill(6);
        if (sync_print_precut)
        cout << "at least one photon present" << endl;
        
        unsigned int electronOneIndex = 0;
        unsigned int electronTwoIndex = 1;
        bool hasValidPair = false;
        float zMassDiff = 100.;
        for (unsigned int i = 0; i < electrons.size(); ++i) {
            for (unsigned int j = i+1; j < electrons.size(); ++j) {
                TLorentzVector tempElectronOne, tempElectronTwo;
                tempElectronOne.SetPtEtaPhiM(electrons[i]->pt, electrons[i]->eta, electrons[i]->phi, ELE_MASS);
                tempElectronTwo.SetPtEtaPhiM(electrons[j]->pt, electrons[j]->eta, electrons[j]->phi, ELE_MASS);
                float thisMass = (tempElectronOne + tempElectronTwo).M();
                if (
                        electrons[i]->q != electrons[j]->q
                        && electrons[i]->pt > 25.0 
                        //&& electrons[i]->pt > 30.0 
                        && electrons[j]->pt > 15.0
                        && thisMass > 50.0
                   ) {
                    if (hasValidPair) {
                        if (fabs(thisMass - ZMASS) < zMassDiff) {
                            zMassDiff = fabs(thisMass - ZMASS);
                            leptonOneP4 = tempElectronOne;
                            leptonTwoP4 = tempElectronTwo;
                            electronOneIndex = i;
                            electronTwoIndex = j;
                        }
                    }
                    else {
                        zMassDiff = fabs(thisMass - ZMASS);
                        leptonOneP4 = tempElectronOne;
                        leptonTwoP4 = tempElectronTwo;
                        electronOneIndex = i;
                        electronTwoIndex = j;
                        hasValidPair = true;
                    }
                }
            }
        }

        if (!hasValidPair)
            return kTRUE;
        hTotalEvents->Fill(7);
        if (sync_print_precut)
            cout << "valid pair" << endl;
        
        bool hasValidPhoton = false;
        unsigned int photonIndex = 0;

        for (unsigned int i = 0; i < photons.size(); ++i) {
            TLorentzVector tempPhoton;
            TLorentzVector tempDilepton;
            TLorentzVector tempLLG;
            tempPhoton.SetPtEtaPhiM(photons[i]->pt, photons[i]->eta, photons[i]->phi, 0.);
            tempDilepton = leptonOneP4 + leptonTwoP4;
            tempLLG = leptonOneP4 + leptonTwoP4 + tempPhoton;
            float this_dr1 = leptonOneP4.DeltaR(tempPhoton);
            float this_dr2 = leptonTwoP4.DeltaR(tempPhoton);
            if (
                tempPhoton.Pt() > 15.0 &&
                tempPhoton.Et()/tempLLG.M() > (15.0/110.0) &&
                tempDilepton.M() + tempLLG.M() > 185.0 &&
                tempLLG.M() > 100. && tempLLG.M() < 180. &&
                this_dr1 > 0.4 && this_dr2 > 0.4
                ) {
                hasValidPhoton = true;
                photonIndex = i;
                break;
            }
            else {
                if (sync_print_precut) {
                    cout << "photon variables in the loop" << endl;
                    cout << "pt, et/m_llg, mass_sum, llg_mass, dr1, dr2" << endl;
                    cout << tempPhoton.Pt() << ", " << tempPhoton.Et()/tempLLG.M() << ", " << 
                            tempDilepton.M() + tempLLG.M() << ", " << tempLLG.M() << ", " << 
                            this_dr1 << ", " << this_dr2 << endl;
                }
            }
            //cout << "photon " << i << " valid: " << hasValidPhoton << endl;
        }

        if (!hasValidPhoton)
            return kTRUE; 
        hTotalEvents->Fill(8);
        if (sync_print_precut)
            cout << "valid photon" << endl;

        //photonOneP4.SetPtEtaPhiM(photons[0]->pt, photons[0]->eta, photons[0]->phi, 0.);
        photonOneP4.SetPtEtaPhiM(photons[photonIndex]->pt, photons[photonIndex]->eta, photons[photonIndex]->phi, 0.);
        if (photonOneP4.Pt() < 15.0)
            return kTRUE;
        hTotalEvents->Fill(9);
        
        // DY photon overlap removal
        vetoDY = false;
        for (unsigned int i = 0; i < genPhotons.size(); ++i) {
            TGenParticle *pho = genPhotons.at(i);
            if (pho->fromHardProcessFinalState || pho->isPromptFinalState) {
                TLorentzVector thisGenPhotonP4;
                thisGenPhotonP4.SetPtEtaPhiM(pho->pt, pho->eta, pho->phi, 0.);
                if (thisGenPhotonP4.DeltaR(photonOneP4) < 0.1) {
                    vetoDY = true;
                    break;
                }
            }
        }  
        
        // checking for lepton tag
        isLeptonTag = false;
        for (unsigned int i = 0; i < muons.size(); ++i) {
            TLorentzVector tempMuon;
            tempMuon.SetPtEtaPhiM(muons[i]->pt, muons[i]->eta, muons[i]->phi, MUON_MASS);
            if (leptonOneP4.DeltaR(tempMuon) < 0.4 ||
                leptonTwoP4.DeltaR(tempMuon) < 0.4 || 
                photonOneP4.DeltaR(tempMuon) < 0.4)
                continue;
            else
                isLeptonTag = true;
        }

        if (!isLeptonTag) {
            for (unsigned int i = 0; i < electrons.size(); ++i) {
                TLorentzVector tempElectron;
                tempElectron.SetPtEtaPhiM(electrons[i]->pt, electrons[i]->eta, electrons[i]->phi, ELE_MASS);
                if (leptonOneP4.DeltaR(tempElectron) < 0.4 ||
                    leptonTwoP4.DeltaR(tempElectron) < 0.4 || 
                    photonOneP4.DeltaR(tempElectron) < 0.4)
                    continue;
                else
                    isLeptonTag = true;
            }
        }
        
        /*// FSR photon recovery
        float minDROne = 100.;
        float minDRTwo = 100.;
        int fsrPhoOneIndex = -1;
        int fsrPhoTwoIndex = -1;

        leptonOneFSRPhotons = 0;
        leptonTwoFSRPhotons = 0;
        leptonOneFSRSum.SetPtEtaPhiM(0., 0., 0., 0.);
        leptonTwoFSRSum.SetPtEtaPhiM(0., 0., 0., 0.);
        leptonOneFSRMatchP4.SetPtEtaPhiM(0., 0., 0., 0.);
        leptonTwoFSRMatchP4.SetPtEtaPhiM(0., 0., 0., 0.);
        leptonOneFSRIsoSum = 0.;
        leptonTwoFSRIsoSum = 0.;
    
        for (unsigned int i = 0; i < fsr_photons.size(); ++i) {
            TPhoton *pho = fsr_photons.at(i);
            TLorentzVector thisFSRPhoP4;
            thisFSRPhoP4.SetPtEtaPhiM(pho->pt, pho->eta, pho->phi, 0.);
            if (thisFSRPhoP4.DeltaR(photonOneP4) < 0.01)
                continue;

            float phoLepOneDR = thisFSRPhoP4.DeltaR(leptonOneP4);
            float phoLepTwoDR = thisFSRPhoP4.DeltaR(leptonTwoP4);
            //lepton 1
            if ((phoLepOneDR/pow(thisFSRPhoP4.Et(), 2)) < 0.012 && phoLepOneDR < 0.5) {
                ++leptonOneFSRPhotons;
                leptonOneFSRSum += thisFSRPhoP4;
                if ((phoLepOneDR > 0.08 || electrons[electronOneIndex]->scEta < 1.479) && phoLepOneDR < 0.4)
                    leptonOneFSRIsoSum += GetPhotonIsolation(pho, fInfo->rhoJet);
                if (phoLepOneDR/pow(thisFSRPhoP4.Et(), 2) < minDROne) {
                    leptonOneFSRMatchP4 = thisFSRPhoP4;
                    fsrPhoOneIndex = i;
                    minDROne = phoLepOneDR/pow(thisFSRPhoP4.Et(), 2);
                }
            }
            
            //lepton 2
            if ((phoLepTwoDR/pow(thisFSRPhoP4.Et(), 2)) < 0.012 && phoLepTwoDR < 0.5) {
                ++leptonTwoFSRPhotons;
                leptonTwoFSRSum += thisFSRPhoP4;
                if ((phoLepTwoDR > 0.08 || electrons[electronTwoIndex]->scEta < 1.479) && phoLepTwoDR < 0.4)
                    leptonTwoFSRIsoSum += GetPhotonIsolation(pho, fInfo->rhoJet);
                if (phoLepTwoDR/pow(thisFSRPhoP4.Et(), 2) < minDRTwo) {
                    leptonTwoFSRMatchP4 = thisFSRPhoP4;
                    fsrPhoTwoIndex = i;
                    minDRTwo = phoLepTwoDR/pow(thisFSRPhoP4.Et(), 2);
                }
            }
        }
       
        genFSRPhotonOneP4.SetPtEtaPhiM(0., 0., 0., 0.);
        genFSRPhotonTwoP4.SetPtEtaPhiM(0., 0., 0., 0.);
        genFSRPhotonOneFHPFS = 0;
        genFSRPhotonTwoFHPFS = 0;
        genFSRPhotonOneIPFS = 0;
        genFSRPhotonTwoIPFS = 0;
        if (fsrPhoOneIndex > 0) {
            fsrPhotonOneMVA = fsr_photons[fsrPhoOneIndex]->mva;
            fsrPassElectronVetoOne = fsr_photons[fsrPhoOneIndex]->passElectronVeto;  
            fsrPhotonOneIso = GetPhotonIsolation(fsr_photons[fsrPhoOneIndex], fInfo->rhoJet);
            if (!isData)
                fsrPhotonOneR9 = weights->GetCorrectedPhotonR9(*fsr_photons[fsrPhoOneIndex]);
            else 
                fsrPhotonOneR9 = fsr_photons[fsrPhoOneIndex]->r9;
            
            // find the matching gen photon
            float min_fsr_dr_one = 1000.;
            for (unsigned int i = 0; i < genPhotons.size(); i++) {
                TLorentzVector tmpGenFSRPhoOne;
                tmpGenFSRPhoOne.SetPtEtaPhiM(genPhotons[i]->pt, genPhotons[i]->eta, genPhotons[i]->phi, genPhotons[i]->mass);
                float this_fsr_dr_one = tmpGenFSRPhoOne.DeltaR(leptonOneFSRMatchP4);
                if (this_fsr_dr_one < min_fsr_dr_one) {
                    genFSRPhotonOneP4.SetPtEtaPhiM(genPhotons[i]->pt, genPhotons[i]->eta, genPhotons[i]->phi, genPhotons[i]->mass);
                    genFSRPhotonOneFHPFS = genPhotons[i]->fromHardProcessFinalState;
                    genFSRPhotonOneIPFS = genPhotons[i]->isPromptFinalState;
                    min_fsr_dr_one = this_fsr_dr_one;
                }
            } 
        }
       
        if (fsrPhoTwoIndex > 0) {
            fsrPhotonTwoMVA = fsr_photons[fsrPhoTwoIndex]->mva;
            fsrPassElectronVetoTwo = fsr_photons[fsrPhoTwoIndex]->passElectronVeto;  
            fsrPhotonTwoIso = GetPhotonIsolation(fsr_photons[fsrPhoTwoIndex], fInfo->rhoJet);
            if (!isData)
                fsrPhotonTwoR9 = weights->GetCorrectedPhotonR9(*fsr_photons[fsrPhoTwoIndex]);
            else 
                fsrPhotonTwoR9 = fsr_photons[fsrPhoTwoIndex]->r9;
            
            // find the matching gen photon
            float min_fsr_dr_two = 1000.;
            for (unsigned int i = 0; i < genPhotons.size(); i++) {
                TLorentzVector tmpGenFSRPhoTwo;
                tmpGenFSRPhoTwo.SetPtEtaPhiM(genPhotons[i]->pt, genPhotons[i]->eta, genPhotons[i]->phi, genPhotons[i]->mass);
                float this_fsr_dr_two = tmpGenFSRPhoTwo.DeltaR(leptonTwoFSRMatchP4);
                if (this_fsr_dr_two < min_fsr_dr_two) {
                    genFSRPhotonTwoP4.SetPtEtaPhiM(genPhotons[i]->pt, genPhotons[i]->eta, genPhotons[i]->phi, genPhotons[i]->mass);
                    genFSRPhotonTwoFHPFS = genPhotons[i]->fromHardProcessFinalState;
                    genFSRPhotonTwoIPFS = genPhotons[i]->isPromptFinalState;
                    min_fsr_dr_two = this_fsr_dr_two;
                }
            } 
        } */
        
        if (sync_print_precut) {
            cout << "lepton1_pt, lepton1_eta, lepton1_phi, lepton2_pt, lepton2_eta, lepton2_phi, dilepton_mass, dr1, dr2" << endl;
            cout << leptonOneP4.Pt() << ", " << leptonOneP4.Eta() << ", " << leptonOneP4.Phi() << ", " 
                 << leptonTwoP4.Pt() << ", " << leptonTwoP4.Eta() << ", " << leptonTwoP4.Phi() << ", "
                 << (leptonOneP4 + leptonTwoP4).M() << ", " << leptonOneP4.DeltaR(photonOneP4) << ", " << leptonTwoP4.DeltaR(photonOneP4) << endl;
        }

        //float m_ll = (leptonOneP4 + leptonTwoP4).M();
        //float m_llg = (leptonOneP4 + leptonTwoP4 + photonOneP4).M();
        //if (photonOneP4.Et() < 0.14*m_llg)
        //    return kTRUE;
        //hTotalEvents->Fill(9);
        //if (leptonOneP4.DeltaR(photonOneP4) < 0.4)
        //    return kTRUE;
        //hTotalEvents->Fill(10);
        //if (leptonTwoP4.DeltaR(photonOneP4) < 0.4)
        //    return kTRUE;
        //hTotalEvents->Fill(11);
        //if (m_ll + m_llg < 185.0)
        //    return kTRUE;
        //hTotalEvents->Fill(12);
        //if (m_ll < 100.0 || m_ll > 180.0)
        //    return kTRUE;
        //hTotalEvents->Fill(13);
 
        leptonOneIso    = GetElectronIsolation(electrons[electronOneIndex], fInfo->rhoJet);
        leptonOneFlavor = electrons[electronOneIndex]->q*11;
        leptonOneDZ     = electrons[electronOneIndex]->dz;
        leptonOneD0     = electrons[electronOneIndex]->d0;
            
        leptonTwoIso    = GetElectronIsolation(electrons[electronTwoIndex], fInfo->rhoJet);
        leptonTwoFlavor = electrons[electronTwoIndex]->q*11;
        leptonTwoDZ     = electrons[electronTwoIndex]->dz;
        leptonTwoD0     = electrons[electronTwoIndex]->d0;

        photonOneMVA = photons[photonIndex]->mva;
        passElectronVeto = photons[photonIndex]->passElectronVeto;  
        if (!isData)
            photonOneR9 = weights->GetCorrectedPhotonR9(*photons[photonIndex]);
        else 
            photonOneR9 = photons[photonIndex]->r9;

        if (!isData) {

            eventWeight *= weights->GetHZZElectronRecoIdEff(*electrons[electronOneIndex]); 
            eventWeight *= weights->GetHZZElectronRecoIdEff(*electrons[electronTwoIndex]); 

            //pair<float, float> eff11 = weights->GetTriggerEffWeight("HLT_DoubleEG_leg1", leptonOneP4);
            //pair<float, float> eff12 = weights->GetTriggerEffWeight("HLT_DoubleEG_leg1", leptonTwoP4);
            //pair<float, float> eff21 = weights->GetTriggerEffWeight("HLT_DoubleEG_leg2", leptonOneP4);
            //pair<float, float> eff22 = weights->GetTriggerEffWeight("HLT_DoubleEG_leg2", leptonTwoP4);
            pair<float, float> eff11 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg1", *electrons[electronOneIndex]);
            pair<float, float> eff12 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg1", *electrons[electronTwoIndex]);
            pair<float, float> eff21 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg2", *electrons[electronOneIndex]);
            pair<float, float> eff22 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg2", *electrons[electronTwoIndex]);
            float eff_data = eff11.first*eff22.first + eff12.first*eff21.first - eff11.first*eff12.first;
            float eff_mc = eff11.second*eff22.second + eff12.second*eff21.second - eff11.second*eff12.second;
            triggerWeight = eff_data/eff_mc;
            //eventWeight *= triggerWeight;
            
            eventWeight *= weights->GetPhotonMVAIdEff(*photons[photonIndex]);
        
            // photon r9 reweighting
            //if (fabs(photons[photonIndex]->scEta) < 1.444)
            //    photonOneR9 = 1.0045*photonOneR9 + 0.0010;
            //else
            //    photonOneR9 = 1.0086*photonOneR9 - 0.0007;
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

    //cout << "new event with " << genPhotons.size() << " gen photons" << endl;
    if (!isData && genPhotons.size() > 0) {
        float min_phot_dr = 1000.;
        for (unsigned int i = 0; i < genPhotons.size(); i++) {
            //cout << "genPhotonFHPFS = " << genPhotons[i]->fromHardProcessFinalState << endl;
            //cout << "genPhotonIPFS = " << genPhotons[i]->isPromptFinalState << endl;
            TLorentzVector tmpGenPhot;
            tmpGenPhot.SetPtEtaPhiM(genPhotons[i]->pt, genPhotons[i]->eta, genPhotons[i]->phi, genPhotons[i]->mass);
            float this_dr = tmpGenPhot.DeltaR(photonOneP4);
            if (this_dr < min_phot_dr) {
                genPhotonP4.SetPtEtaPhiM(genPhotons[i]->pt, genPhotons[i]->eta, genPhotons[i]->phi, genPhotons[i]->mass);
                genPhotonFHPFS = genPhotons[i]->fromHardProcessFinalState;
                genPhotonIPFS = genPhotons[i]->isPromptFinalState;
                min_phot_dr = this_dr;
            }
        }
    }

    outTree->Fill();
    this->passedEvents++;

    if (sync_print_precut) {
        cout << "event should have been filled" << endl;
        }
    return kTRUE;
}

void hzgAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void hzgAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void hzgAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<hzgAnalyzer> selector(new hzgAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

float hzgAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    //float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    float combIso = (mu->chHadIso03 + std::max(0.,(double)mu->neuHadIso03 + mu->gammaIso03 - 0.5*mu->puIso03));
    return combIso;
}

float hzgAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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

float hzgAnalyzer::GetPhotonIsolation(const baconhep::TPhoton* pho, const float rho)
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
