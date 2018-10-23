#include "gbrTrainAnalyzer.h"
#include <map>
#include <fstream>
#include <math.h>

#include <TSystem.h>
#include <TF2.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

// prints a message and exits gracefully
#define FATAL(msg) do { fprintf(stderr, "FATAL: %s\n", msg); gSystem->Exit(1); } while (0)

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

gbrTrainAnalyzer::gbrTrainAnalyzer() : BLTSelector()
{

}

gbrTrainAnalyzer::~gbrTrainAnalyzer()
{

}

void gbrTrainAnalyzer::Begin(TTree *tree)
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

    if (params->selection == "dilepton") {
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
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

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();

    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",30,0.5,30.5);

    vector<std::string> channelNames = {"mm", "ee"};

    for (unsigned int i = 0; i < channelNames.size(); ++i) {
        string channel = channelNames[i];
        outFile->mkdir(channel.c_str());
        outFile->cd(channel.c_str());
        string treeName = "bltTree_" + params->datasetgroup;
        tree = new TTree(treeName.c_str(), treeName.c_str());

        // event data
        tree->Branch("runNumber", &runNumber);
        tree->Branch("evtNumber", &evtNumber, "eventNumber/l");
        tree->Branch("lumiSection", &lumiSection);
        tree->Branch("triggerStatus", &triggerStatus);
        tree->Branch("nPV", &nPV);
        tree->Branch("nPU", &nPU);
        tree->Branch("nPartons", &nPartons);
        tree->Branch("rPV", &rPV);
        tree->Branch("rho", &rho);
        
        // weights and their uncertainties
        tree->Branch("eventWeight", &eventWeight);
        tree->Branch("leptonOneRecoWeight", &leptonOneRecoWeight);
        tree->Branch("leptonTwoRecoWeight", &leptonTwoRecoWeight);
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
        tree->Branch("leptonOneP4KinFit", &leptonOneP4KinFit);
        tree->Branch("leptonOneIso", &leptonOneIso);
        tree->Branch("leptonOneFlavor", &leptonOneFlavor);
        tree->Branch("leptonOneMother", &leptonOneMother);
        tree->Branch("leptonOneD0", &leptonOneD0);
        tree->Branch("leptonOneDZ", &leptonOneDZ);
        tree->Branch("leptonOnePt", &leptonOnePt);
        tree->Branch("leptonOneEta", &leptonOneEta);
        tree->Branch("leptonOneEnergy", &leptonOneEnergy);
        tree->Branch("leptonOneERes", &leptonOneERes);
        tree->Branch("leptonOnePFIsoCH", &leptonOnePFIsoCH);
        tree->Branch("leptonOnePFIsoNH", &leptonOnePFIsoNH);
        tree->Branch("leptonOnePFIsoPho", &leptonOnePFIsoPho);
        tree->Branch("leptonOnePFIsoPU", &leptonOnePFIsoPU);
        tree->Branch("leptonOneMean", &leptonOneMean);
        tree->Branch("leptonOneSigma", &leptonOneSigma);
        tree->Branch("leptonOneAlphaL", &leptonOneAlphaL);
        tree->Branch("leptonOneAlphaR", &leptonOneAlphaR);
        tree->Branch("leptonOnePowerL", &leptonOnePowerL);
        tree->Branch("leptonOnePowerR", &leptonOnePowerR);

        tree->Branch("leptonTwoP4", &leptonTwoP4);
        tree->Branch("leptonTwoP4KinFit", &leptonTwoP4KinFit);
        tree->Branch("leptonTwoIso", &leptonTwoIso);
        tree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
        tree->Branch("leptonTwoMother", &leptonTwoMother);
        tree->Branch("leptonTwoD0", &leptonTwoD0);
        tree->Branch("leptonTwoDZ", &leptonTwoDZ);
        tree->Branch("leptonTwoPt", &leptonTwoPt);
        tree->Branch("leptonTwoEta", &leptonTwoEta);
        tree->Branch("leptonTwoEnergy", &leptonTwoEnergy);
        tree->Branch("leptonTwoERes", &leptonTwoERes);
        tree->Branch("leptonTwoPFIsoCH", &leptonTwoPFIsoCH);
        tree->Branch("leptonTwoPFIsoNH", &leptonTwoPFIsoNH);
        tree->Branch("leptonTwoPFIsoPho", &leptonTwoPFIsoPho);
        tree->Branch("leptonTwoPFIsoPU", &leptonTwoPFIsoPU);
        tree->Branch("leptonTwoMean", &leptonTwoMean);
        tree->Branch("leptonTwoSigma", &leptonTwoSigma);
        tree->Branch("leptonTwoAlphaL", &leptonTwoAlphaL);
        tree->Branch("leptonTwoAlphaR", &leptonTwoAlphaR);
        tree->Branch("leptonTwoPowerL", &leptonTwoPowerL);
        tree->Branch("leptonTwoPowerR", &leptonTwoPowerR);
 
        // lepton vertices
        tree->Branch("dileptonVertexOne", &dileptonVertexOne);
        tree->Branch("dileptonVertexErrOne", &dileptonVertexErrOne);
        tree->Branch("dileptonVertexChi2One", &dileptonVertexChi2One);
        tree->Branch("dileptonVertexDOFOne", &dileptonVertexDOFOne);

        // gen level objects 
        tree->Branch("genLeptonOneP4", &genLeptonOneP4);
        tree->Branch("genLeptonOneId", &genLeptonOneId);
        tree->Branch("genLeptonTwoP4", &genLeptonTwoP4);
        tree->Branch("genLeptonTwoId", &genLeptonTwoId);
        tree->Branch("genPhotonP4", &genPhotonP4);
        tree->Branch("genPhotonFHPFS", &genPhotonFHPFS);
        tree->Branch("genPhotonIPFS", &genPhotonIPFS);
        tree->Branch("vetoDY", &vetoDY);

        // object counters
        tree->Branch("nMuons", &nMuons);
        tree->Branch("nElectrons", &nElectrons);
        tree->Branch("nTaus", &nTaus);
        tree->Branch("nPhotons", &nPhotons);
        tree->Branch("nJets", &nJets);
        tree->Branch("nFwdJets", &nFwdJets);
        tree->Branch("nBJets", &nBJets);

        if (channel == "mm") {
            tree->Branch("leptonOneTkNChi2", &leptonOneTkNChi2);
            tree->Branch("leptonOneMuNChi2", &leptonOneMuNChi2);
            tree->Branch("leptonOneNTkLayers", &leptonOneNTkLayers);
            tree->Branch("leptonOneNPixHits", &leptonOneNPixHits);
            tree->Branch("leptonOneNValidHits", &leptonOneNValidHits);
            tree->Branch("leptonOneNMatchStn", &leptonOneNMatchStn);

            tree->Branch("leptonTwoTkNChi2", &leptonTwoTkNChi2);
            tree->Branch("leptonTwoMuNChi2", &leptonTwoMuNChi2);
            tree->Branch("leptonTwoNTkLayers", &leptonTwoNTkLayers);
            tree->Branch("leptonTwoNPixHits", &leptonTwoNPixHits);
            tree->Branch("leptonTwoNValidHits", &leptonTwoNValidHits);
            tree->Branch("leptonTwoNMatchStn", &leptonTwoNMatchStn);
        }

        else if (channel == "ee") {
            tree->Branch("leptonOneScEta", &leptonOneScEta);
            tree->Branch("leptonOneScPhi", &leptonOneScPhi);
            tree->Branch("leptonOneR9", &leptonOneR9);
            tree->Branch("leptonOneE1x5OverE", &leptonOneE1x5OverE);
            tree->Branch("leptonOneE2x5OverE", &leptonOneE2x5OverE);
            tree->Branch("leptonOneE5x5OverE", &leptonOneE5x5OverE);
            tree->Branch("leptonOneFBrem", &leptonOneFBrem);
            tree->Branch("leptonOneHOverE", &leptonOneHOverE);
            tree->Branch("leptonOneSigmaIEtaIEta", &leptonOneSigmaIEtaIEta);
            
            tree->Branch("leptonTwoScEta", &leptonTwoScEta);
            tree->Branch("leptonTwoScPhi", &leptonTwoScPhi);
            tree->Branch("leptonTwoR9", &leptonTwoR9);
            tree->Branch("leptonTwoE1x5OverE", &leptonTwoE1x5OverE);
            tree->Branch("leptonTwoE2x5OverE", &leptonTwoE2x5OverE);
            tree->Branch("leptonTwoE5x5OverE", &leptonTwoE5x5OverE);
            tree->Branch("leptonTwoFBrem", &leptonTwoFBrem);
            tree->Branch("leptonTwoHOverE", &leptonTwoHOverE);
            tree->Branch("leptonTwoSigmaIEtaIEta", &leptonTwoSigmaIEtaIEta);
        }
        
        outTrees[channel] = tree;
        
        // event counter
        string outHistName = params->get_output_treename("TotalEvents_" + channel);
        eventCounts[channel] = new TH1D(outHistName.c_str(),"ChannelCounts",10,0.5,10.5);
    }

    ReportPostBegin();
}

Bool_t gbrTrainAnalyzer::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
    outFile->cd();
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
    
    genWeight = 1;
    if (!isData) {
        if (fGenEvtInfo->weight < 0) {
            genWeight = -1;
            int maxBin = hTotalEvents->GetSize() - 2;
            hTotalEvents->Fill(maxBin);
        }
    }
     
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

            if ((fabs(particle->pdgId) == 11 || fabs(particle->pdgId) == 13) and particle->parent > 0) {
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                if (fabs(mother->pdgId) == 23) 
                    genLeptons.push_back(particle);
            }
                
            // saving photons     
            if (fabs(particle->pdgId) == 22) 
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
    for (unsigned int i = 0; i < triggerNames.size(); ++i) {
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

        // Remove muons with very small deltaR // THIS IS NOT YET USED
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
               /*// tight muon ID and ISO
               muonP4.Pt() > 5.
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
           ){
            muons.push_back(muon);
        } */

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
                muonP4.Pt() > 10.
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
    vector<TLorentzVector> veto_taus;
    for (int i=0; i < fTauArr->GetEntries(); i++) {
        TTau *tau = (TTau*) fTauArr->At(i);
        assert(tau);

        TLorentzVector tauP4; 
        tauP4.SetPtEtaPhiM(tau->pt, tau->eta, tau->phi, tau->m);

        // Prevent overlap of muons and jets
        bool muOverlap = false;
        for (const auto& mu: veto_muons) {
            if (tauP4.DeltaR(mu) < 0.3) {
                muOverlap = true;
                break;
            }
        }
        bool elOverlap = false;
        for (const auto& el: veto_electrons) {
            if (tauP4.DeltaR(el) < 0.3) {
                elOverlap = true;
                break;
            }
        }

        // apply tau energy scale correction (https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Tau_energy_scale)
        if (!isData) {
            if (tau->decaymode == 0) {
                tau->pt *= 0.995;
            } else if (tau->decaymode == 1) {
                tau->pt *= 1.01;
            } else if (tau->decaymode == 10) {
                tau->pt *= 1.006;
            }
        }

        if ( 
                tau->pt > 18
                && abs(tau->eta) < 2.3
                && !muOverlap
                && !elOverlap
                && (tau->hpsDisc & baconhep::kByDecayModeFinding)
                && (tau->hpsDisc & baconhep::kByTightIsolationMVA3newDMwLT)
                && (tau->hpsDisc & baconhep::kByMVA6VTightElectronRejection)
                && (tau->hpsDisc & baconhep::kByTightMuonRejection3)
          ) {
            taus.push_back(tau);
            veto_taus.push_back(tauP4);
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

        /*bool muOverlap = false;
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
        }*/

        if (
                jet->pt > 30 
                && fabs(jet->eta) < 4.7
                && particleSelector->PassJetID(jet, cuts->looseJetID)
                //&& !muOverlap 
                //&& !elOverlap
                //&& !phoOverlap
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
                    if (particleSelector->BTagModifier(jet, "MVAT", 0, 0, rng->Uniform(1.))) { 
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
    
    // trigger selections
    bool muonTriggered = find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*") != passTriggerNames.end() 
                      || find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*") != passTriggerNames.end()
                      || find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*") != passTriggerNames.end()
                      || find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*") != passTriggerNames.end();
    
    bool electronTriggered = find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*") != passTriggerNames.end()
                          || find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*") != passTriggerNames.end();

    string channel = "";
    if (params->selection == "dilepton") {
       
        if (muons.size() > 1 && muonTriggered) { // mm channel
            channel = "mm";
            eventCounts[channel]->Fill(1);
        
            // convert to TLorentzVectors
            TLorentzVector muonOneP4, muonTwoP4, dimuonP4;
            muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, MUON_MASS);
            muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, MUON_MASS);
            dimuonP4 = muonOneP4 + muonTwoP4;

            if (muons[0]->pt < 20.0)
                return kTRUE;
            eventCounts[channel]->Fill(2);
            if (muons[1]->pt < 10.0)
                return kTRUE;
            eventCounts[channel]->Fill(3);
            if (muons[0]->q == muons[1]->q)
                return kTRUE;
            eventCounts[channel]->Fill(4);

            if (dimuonP4.M() < 80.0 || dimuonP4.M() > 100.0)
                return kTRUE;
            eventCounts[channel]->Fill(5);
        
            // L1EMTF cut 
            if (
                fabs(muonOneP4.DeltaPhi(muonTwoP4)) < 70.0*(M_PI/180.0)
                && fabs(muonOneP4.Eta()) > 1.2 
                && fabs(muonTwoP4.Eta()) > 1.2
                && muonOneP4.Eta()*muonTwoP4.Eta() > 0
               )
                return kTRUE;
            eventCounts[channel]->Fill(6); 
        
            /*// trigger matching:
            bool mu1_fired_leg1 = false;
            bool mu1_fired_leg2 = false;
            bool mu2_fired_leg1 = false;
            bool mu2_fired_leg2 = false;
            for (unsigned int iT = 0; iT < triggerNames.size(); ++iT) {
                mu1_fired_leg1 |= trigger->passObj(triggerNames.at(iT), 1, muons[0]->hltMatchBits);
                mu1_fired_leg2 |= trigger->passObj(triggerNames.at(iT), 2, muons[0]->hltMatchBits);
                mu2_fired_leg1 |= trigger->passObj(triggerNames.at(iT), 1, muons[1]->hltMatchBits);
                mu2_fired_leg2 |= trigger->passObj(triggerNames.at(iT), 2, muons[1]->hltMatchBits);
            }
            bool isTriggerMatched = (mu1_fired_leg1 || mu1_fired_leg2) && (mu2_fired_leg1 || mu2_fired_leg2);
            if (!isTriggerMatched)
                return kTRUE;
            eventCounts[channel]->Fill(7); */

            // gen matching
            bool isGenMatchedOne = false;
            bool isGenMatchedTwo = false;
            float EResOne = 0.;
            float EResTwo = 0.;
            for (unsigned int i = 0; i < genLeptons.size(); ++i) {
                if (fabs(genLeptons[i]->pdgId) == 13) {
                    TLorentzVector genMuonP4;
                    genMuonP4.SetPtEtaPhiM(genLeptons[i]->pt, genLeptons[i]->eta, genLeptons[i]->phi, MUON_MASS);
                    float genRecoDROne = genMuonP4.DeltaR(muonOneP4);
                    float genRecoDRTwo = genMuonP4.DeltaR(muonTwoP4);
                    if (genRecoDROne < 0.1) {
                        isGenMatchedOne = true;
                        EResOne = genMuonP4.E()/muonOneP4.E();
                        continue;
                    }
                    if (genRecoDRTwo < 0.1) {
                        isGenMatchedTwo = true;
                        EResTwo = genMuonP4.E()/muonTwoP4.E();
                    }
                }
            }

            if (!isGenMatchedOne)
                return kTRUE;
            eventCounts[channel]->Fill(7);
            if (!isGenMatchedTwo)
                return kTRUE;
            eventCounts[channel]->Fill(8); 

            float muonOneMean, muonOneSigma, muonOneAlphaL, muonOneAlphaR, muonOnePowerL, muonOnePowerR;
            float muonTwoMean, muonTwoSigma, muonTwoAlphaL, muonTwoAlphaR, muonTwoPowerL, muonTwoPowerR;
            std::map<string, float> mva_input_floats;
            std::map<string, int> mva_input_ints;
   
            // BDT eval for muon one
            mva_input_floats["rho"] = fInfo->rhoJet;
            mva_input_floats["muEnergy"] = muonOneP4.E();
            mva_input_floats["muEta"] = muons[0]->eta;
            mva_input_floats["muTkChi2"] = muons[0]->tkNchi2;
            mva_input_ints["muNumberOfValidTrkLayers"] = muons[0]->nTkLayers;
            mva_input_ints["muNumberOfValidPixelHits"] = muons[0]->nPixHits;
            mva_input_ints["muNumberOfValidMuonHits"] = muons[0]->nValidHits;
            mva_input_ints["muStations"] = muons[0]->nMatchStn;
            mva_input_floats["muPFIsoR04_CH"] = muons[0]->chHadIso;
            mva_input_floats["muPFIsoR04_NH"] = muons[0]->neuHadIso;
            mva_input_floats["muPFIsoR04_Pho"] = muons[0]->gammaIso;
            mva_input_floats["muPFIsoR04_PU"] = muons[0]->puIso;
            EvalMuonEnergyResolution(mva_input_floats, mva_input_ints, muonOneMean, muonOneSigma, 
                                     muonOneAlphaL, muonOnePowerL, muonOneAlphaR, muonOnePowerR);
            
            // BDT eval for muon two
            mva_input_floats["muEnergy"] = muonTwoP4.E();
            mva_input_floats["muEta"] = muons[1]->eta;
            mva_input_floats["muTkChi2"] = muons[1]->tkNchi2;
            mva_input_ints["muNumberOfValidTrkLayers"] = muons[1]->nTkLayers;
            mva_input_ints["muNumberOfValidPixelHits"] = muons[1]->nPixHits;
            mva_input_ints["muNumberOfValidMuonHits"] = muons[1]->nValidHits;
            mva_input_ints["muStations"] = muons[1]->nMatchStn;
            mva_input_floats["muPFIsoR04_CH"] = muons[1]->chHadIso;
            mva_input_floats["muPFIsoR04_NH"] = muons[1]->neuHadIso;
            mva_input_floats["muPFIsoR04_Pho"] = muons[1]->gammaIso;
            mva_input_floats["muPFIsoR04_PU"] = muons[1]->puIso;
            EvalMuonEnergyResolution(mva_input_floats, mva_input_ints, muonTwoMean, muonTwoSigma, 
                                     muonTwoAlphaL, muonTwoPowerL, muonTwoAlphaR, muonTwoPowerR);

            leptonOneP4     = muonOneP4;
            leptonOneIso    = GetMuonIsolation(muons[0]);
            leptonOneFlavor = muons[0]->q*13;
            leptonOneDZ     = muons[0]->dz;
            leptonOneD0     = muons[0]->d0;
            leptonOnePt     = muons[0]->pt;
            leptonOneEta    = muons[0]->eta;
            leptonOneEnergy = muonOneP4.E();
            leptonOneERes   = EResOne;
            leptonOneTkNChi2 = muons[0]->tkNchi2;
            leptonOneMuNChi2 = muons[0]->muNchi2;
            leptonOneNTkLayers = muons[0]->nTkLayers;
            leptonOneNPixHits = muons[0]->nPixHits;
            leptonOneNValidHits = muons[0]->nValidHits;
            leptonOneNMatchStn = muons[0]->nMatchStn;
            leptonOnePFIsoCH = muons[0]->chHadIso;
            leptonOnePFIsoNH = muons[0]->neuHadIso;
            leptonOnePFIsoPho = muons[0]->gammaIso;
            leptonOnePFIsoPU = muons[0]->puIso;
            leptonOneMean = muonOneMean;
            leptonOneSigma = muonOneSigma;
            leptonOneAlphaL = muonOneAlphaL;
            leptonOneAlphaR = muonOneAlphaR;
                
            leptonTwoP4     = muonTwoP4;
            leptonTwoIso    = GetMuonIsolation(muons[1]);
            leptonTwoFlavor = muons[1]->q*13;
            leptonTwoDZ     = muons[1]->dz;
            leptonTwoD0     = muons[1]->d0;
            leptonTwoPt     = muons[1]->pt;
            leptonTwoEta    = muons[1]->eta;
            leptonTwoEnergy = muonTwoP4.E();
            leptonTwoERes   = EResTwo;
            leptonTwoTkNChi2 = muons[1]->tkNchi2;
            leptonTwoMuNChi2 = muons[1]->muNchi2;
            leptonTwoNTkLayers = muons[1]->nTkLayers;
            leptonTwoNPixHits = muons[1]->nPixHits;
            leptonTwoNValidHits = muons[1]->nValidHits;
            leptonTwoNMatchStn = muons[1]->nMatchStn;
            leptonTwoPFIsoCH = muons[1]->chHadIso;
            leptonTwoPFIsoNH = muons[1]->neuHadIso;
            leptonTwoPFIsoPho = muons[1]->gammaIso;
            leptonTwoPFIsoPU = muons[1]->puIso;
            leptonTwoMean = muonTwoMean;
            leptonTwoSigma = muonTwoSigma;
            leptonTwoAlphaL = muonTwoAlphaL;
            leptonTwoAlphaR = muonTwoAlphaR;

            rho             = fInfo->rhoJet;


            // Kinematic Fit
            TVector3 v1, v2;
            v1.SetPtEtaPhi(muons[0]->pt, muons[0]->eta, muons[0]->phi);
            v2.SetPtEtaPhi(muons[1]->pt, muons[1]->eta, muons[1]->phi);
            Double_t p[16];
            p[0] = MUON_MASS*MUON_MASS;
            p[1] = cos(v1.Angle(v2));
            p[2] = muons[0]->pt*cosh(muons[0]->eta);
            p[3] = muons[1]->pt*cosh(muons[1]->eta);
            p[4] = muonOneMean;
            p[5] = muonOneSigma;
            p[6] = muonOneAlphaL;
            p[7] = muonOnePowerL;
            p[8] = muonOneAlphaR;
            p[9] = muonOnePowerR;
            p[10] = muonTwoMean;
            p[11] = muonTwoSigma;
            p[12] = muonTwoAlphaL;
            p[13] = muonTwoPowerL;
            p[14] = muonTwoAlphaR;
            p[15] = muonTwoPowerR;

      
            double e1, e2;
            find_optimized(p, e1, e2);

            double muonOnePtKinFit = e1/cosh(muons[0]->eta);
            double muonTwoPtKinFit = e2/cosh(muons[1]->eta);

            TLorentzVector muonOneP4KinFit, muonTwoP4KinFit;
            muonOneP4KinFit.SetPtEtaPhiM(muonOnePtKinFit, muons[0]->eta, muons[0]->phi, MUON_MASS);
            muonTwoP4KinFit.SetPtEtaPhiM(muonTwoPtKinFit, muons[1]->eta, muons[1]->phi, MUON_MASS);

            leptonOneP4KinFit = muonOneP4KinFit;
            leptonTwoP4KinFit = muonTwoP4KinFit;

    
            if (!isData) {
                eventWeight *= weights->GetHZZMuonIDEff(*muons[0]); 
                eventWeight *= weights->GetHZZMuonIDEff(*muons[1]);
            
                pair<float, float> eff11 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg1", *muons[0]);
                pair<float, float> eff12 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg1", *muons[1]);
                pair<float, float> eff21 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg2", *muons[0]);
                pair<float, float> eff22 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg2", *muons[1]);
                float eff_data = eff11.first*eff22.first + eff12.first*eff21.first - eff11.first*eff12.first;
                float eff_mc = eff11.second*eff22.second + eff12.second*eff21.second - eff11.second*eff12.second;
                triggerWeight = eff_data/eff_mc;
                eventWeight *= triggerWeight;
            }

        } // end mm channel

        else if (electrons.size() > 1 && electronTriggered && !muonTriggered) { // ee channel
            channel = "ee";

            eventCounts[channel]->Fill(1);
            
            // convert to TLorentzVectors
            TLorentzVector electronOneP4, electronTwoP4, dielectronP4;
            electronOneP4.SetPtEtaPhiM(electrons[0]->calibPt, electrons[0]->eta, electrons[0]->phi, ELE_MASS);
            electronTwoP4.SetPtEtaPhiM(electrons[1]->calibPt, electrons[1]->eta, electrons[1]->phi, ELE_MASS);
            dielectronP4 = electronOneP4 + electronTwoP4;

            if (electrons[0]->pt < 25.0)
                return kTRUE;
            eventCounts[channel]->Fill(2);
            if (electrons[1]->pt < 15.0)
                return kTRUE;
            eventCounts[channel]->Fill(3);
            if (electrons[0]->q == electrons[1]->q)
                return kTRUE;
            eventCounts[channel]->Fill(4);

            if (dielectronP4.M() < 80.0 || dielectronP4.M() > 100.0)
                return kTRUE;
            eventCounts[channel]->Fill(5);
         
            /*// trigger matching:
            bool ele1_fired_leg1 = false;
            bool ele1_fired_leg2 = false;
            bool ele2_fired_leg1 = false;
            bool ele2_fired_leg2 = false;
            for (unsigned int iT = 0; iT < triggerNames.size(); ++iT) {
                ele1_fired_leg1 |= trigger->passObj(triggerNames.at(iT), 1, electrons[0]->hltMatchBits);
                ele1_fired_leg2 |= trigger->passObj(triggerNames.at(iT), 2, electrons[0]->hltMatchBits);
                ele2_fired_leg1 |= trigger->passObj(triggerNames.at(iT), 1, electrons[1]->hltMatchBits);
                ele2_fired_leg2 |= trigger->passObj(triggerNames.at(iT), 2, electrons[1]->hltMatchBits);
            }
            bool isTriggerMatched = (ele1_fired_leg1 || ele1_fired_leg2) && (ele2_fired_leg1 || ele2_fired_leg2);
            if (!isTriggerMatched)
                return kTRUE;
            eventCounts[channel]->Fill(7); */

            // gen matching
            bool isGenMatchedOne = false;
            bool isGenMatchedTwo = false;
            float EResOne = 0.;
            float EResTwo = 0.;
            for (unsigned int i = 0; i < genLeptons.size(); ++i) {
                if (fabs(genLeptons[i]->pdgId) == 11) {
                    TLorentzVector genElectronP4;
                    genElectronP4.SetPtEtaPhiM(genLeptons[i]->pt, genLeptons[i]->eta, genLeptons[i]->phi, ELE_MASS);
                    float genRecoDROne = genElectronP4.DeltaR(electronOneP4);
                    float genRecoDRTwo = genElectronP4.DeltaR(electronTwoP4);
                    if (genRecoDROne < 0.1) {
                        isGenMatchedOne = true;
                        EResOne = genElectronP4.E()/electronOneP4.E();
                        continue;
                    }
                    if (genRecoDRTwo < 0.1) {
                        isGenMatchedTwo = true;
                        EResTwo = genElectronP4.E()/electronTwoP4.E();
                    }
                }
            }

            if (!isGenMatchedOne)
                return kTRUE;
            eventCounts[channel]->Fill(6);
            if (!isGenMatchedTwo)
                return kTRUE;
            eventCounts[channel]->Fill(7); 
            
            float electronOneMean, electronOneSigma, electronOneAlphaL, electronOneAlphaR, electronOnePowerL, electronOnePowerR;
            float electronTwoMean, electronTwoSigma, electronTwoAlphaL, electronTwoAlphaR, electronTwoPowerL, electronTwoPowerR;
            std::map<string, float> mva_inputs;
   
            // BDT eval for electron one
            mva_inputs["rho"] = fInfo->rhoJet;
            mva_inputs["elEnergy"] = electronOneP4.E();
            mva_inputs["elScEta"] = electrons[0]->scEta;
            mva_inputs["elScPhi"] = electrons[0]->scPhi;
            mva_inputs["elR9"] = electrons[0]->r9;
            mva_inputs["elE1x5OverE"] = electrons[0]->e1x5/electronOneP4.E();
            mva_inputs["elE2x5OverE"] = electrons[0]->e2x5/electronOneP4.E();
            mva_inputs["elE5x5OverE"] = electrons[0]->e5x5/electronOneP4.E();
            mva_inputs["elFBrem"] = electrons[0]->fbrem;
            mva_inputs["elHOverE"] = electrons[0]->hovere;
            mva_inputs["elSigmaIEtaIEta"] = electrons[0]->sieie;
            mva_inputs["elPFIsoR04_CH"] = electrons[0]->chHadIso;
            mva_inputs["elPFIsoR04_NH"] = electrons[0]->neuHadIso;
            mva_inputs["elPFIsoR04_Pho"] = electrons[0]->gammaIso;
            mva_inputs["elPFIsoR04_PU"] = electrons[0]->puIso;
            EvalElectronEnergyResolution(mva_inputs, electronOneMean, electronOneSigma, 
                                         electronOneAlphaL, electronOnePowerL, electronOneAlphaR, electronOnePowerR);
            
            // BDT eval for electron two
            mva_inputs["elEnergy"] = electronTwoP4.E();
            mva_inputs["elScEta"] = electrons[1]->scEta;
            mva_inputs["elScPhi"] = electrons[1]->scPhi;
            mva_inputs["elR9"] = electrons[1]->r9;
            mva_inputs["elE1x5OverE"] = electrons[1]->e1x5/electronTwoP4.E();
            mva_inputs["elE2x5OverE"] = electrons[1]->e2x5/electronTwoP4.E();
            mva_inputs["elE5x5OverE"] = electrons[1]->e5x5/electronTwoP4.E();
            mva_inputs["elFBrem"] = electrons[1]->fbrem;
            mva_inputs["elHOverE"] = electrons[1]->hovere;
            mva_inputs["elSigmaIEtaIEta"] = electrons[1]->sieie;
            mva_inputs["elPFIsoR04_CH"] = electrons[1]->chHadIso;
            mva_inputs["elPFIsoR04_NH"] = electrons[1]->neuHadIso;
            mva_inputs["elPFIsoR04_Pho"] = electrons[1]->gammaIso;
            mva_inputs["elPFIsoR04_PU"] = electrons[1]->puIso;
            EvalElectronEnergyResolution(mva_inputs, electronTwoMean, electronTwoSigma, 
                                         electronTwoAlphaL, electronTwoPowerL, electronTwoAlphaR, electronTwoPowerR);

            leptonOneP4     = electronOneP4;
            leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
            leptonOneFlavor = electrons[0]->q*11;
            leptonOneDZ     = electrons[0]->dz;
            leptonOneD0     = electrons[0]->d0;
            leptonOnePt     = electrons[0]->calibPt;
            leptonOneEnergy = electronOneP4.E();
            leptonOneERes   = EResOne;
            leptonOneScEta    = electrons[0]->scEta;
            leptonOneScPhi    = electrons[0]->scPhi;
            leptonOneR9     = electrons[0]->r9;
            leptonOneE1x5OverE = electrons[0]->e1x5/electronOneP4.E();
            leptonOneE2x5OverE = electrons[0]->e2x5/electronOneP4.E();
            leptonOneE5x5OverE = electrons[0]->e5x5/electronOneP4.E();
            leptonOneFBrem  = electrons[0]->fbrem;
            leptonOneHOverE = electrons[0]->hovere;
            leptonOneSigmaIEtaIEta = electrons[0]->sieie;
            leptonOnePFIsoCH = electrons[0]->chHadIso;
            leptonOnePFIsoNH = electrons[0]->neuHadIso;
            leptonOnePFIsoPho = electrons[0]->gammaIso;
            leptonOnePFIsoPU = electrons[0]->puIso;
            leptonOneMean = electronOneMean;
            leptonOneSigma = electronOneSigma;
            leptonOneAlphaL = electronOneAlphaL;
            leptonOneAlphaR = electronOneAlphaR;
            leptonOnePowerL = electronOnePowerL;
            leptonOnePowerR = electronOnePowerR;
            
            leptonTwoP4     = electronTwoP4;
            leptonTwoIso    = GetElectronIsolation(electrons[1], fInfo->rhoJet);
            leptonTwoFlavor = electrons[1]->q*11;
            leptonTwoDZ     = electrons[1]->dz;
            leptonTwoD0     = electrons[1]->d0;
            leptonTwoPt     = electrons[1]->calibPt;
            leptonTwoEnergy = electronTwoP4.E();
            leptonTwoERes   = EResTwo;
            leptonTwoScEta    = electrons[1]->scEta;
            leptonTwoScPhi    = electrons[1]->scPhi;
            leptonTwoR9     = electrons[1]->r9;
            leptonTwoE1x5OverE = electrons[1]->e1x5/electronTwoP4.E();
            leptonTwoE2x5OverE = electrons[1]->e2x5/electronTwoP4.E();
            leptonTwoE5x5OverE = electrons[1]->e5x5/electronTwoP4.E();
            leptonTwoFBrem  = electrons[1]->fbrem;
            leptonTwoHOverE = electrons[1]->hovere;
            leptonTwoSigmaIEtaIEta = electrons[1]->sieie;
            leptonTwoPFIsoCH = electrons[1]->chHadIso;
            leptonTwoPFIsoNH = electrons[1]->neuHadIso;
            leptonTwoPFIsoPho = electrons[1]->gammaIso;
            leptonTwoPFIsoPU = electrons[1]->puIso;
            leptonTwoMean = electronTwoMean;
            leptonTwoSigma = electronTwoSigma;
            leptonTwoAlphaL = electronTwoAlphaL;
            leptonTwoAlphaR = electronTwoAlphaR;
            leptonTwoPowerL = electronTwoPowerL;
            leptonTwoPowerR = electronTwoPowerR;
                
            rho             = fInfo->rhoJet;

            // Kinematic Fit
            TVector3 v1, v2;
            v1.SetPtEtaPhi(electrons[0]->calibPt, electrons[0]->eta, electrons[0]->phi);
            v2.SetPtEtaPhi(electrons[1]->calibPt, electrons[1]->eta, electrons[1]->phi);
            Double_t p[16];
            p[0] = ELE_MASS*ELE_MASS;
            p[1] = cos(v1.Angle(v2));
            p[2] = electrons[0]->calibPt*cosh(electrons[0]->eta);
            p[3] = electrons[1]->calibPt*cosh(electrons[1]->eta);
            p[4] = electronOneMean;
            p[5] = electronOneSigma;
            p[6] = electronOneAlphaL;
            p[7] = electronOnePowerL;
            p[8] = electronOneAlphaR;
            p[9] = electronOnePowerR;
            p[10] = electronTwoMean;
            p[11] = electronTwoSigma;
            p[12] = electronTwoAlphaL;
            p[13] = electronTwoPowerL;
            p[14] = electronTwoAlphaR;
            p[15] = electronTwoPowerR;
    
            double e1, e2;
            find_optimized(p, e1, e2);

            double electronOnePtKinFit = e1/cosh(electrons[0]->eta);
            double electronTwoPtKinFit = e2/cosh(electrons[1]->eta);

            TLorentzVector electronOneP4KinFit, electronTwoP4KinFit;
            electronOneP4KinFit.SetPtEtaPhiM(electronOnePtKinFit, electrons[0]->eta, electrons[0]->phi, ELE_MASS);
            electronTwoP4KinFit.SetPtEtaPhiM(electronTwoPtKinFit, electrons[1]->eta, electrons[1]->phi, ELE_MASS);

            leptonOneP4KinFit = electronOneP4KinFit;
            leptonTwoP4KinFit = electronTwoP4KinFit; 
    
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
        
        } // end ee channel

        else 
            return kTRUE;

        // DY photon overlap removal
        vetoDY = false;
        if (photons.size() > 0) {
            TLorentzVector photonP4;
            photonP4.SetPtEtaPhiM(photons[0]->calibPt, photons[0]->eta, photons[0]->phi, 0.);
            for (unsigned int i = 0; i < genPhotons.size(); ++i) {
                TGenParticle *pho = genPhotons.at(i);
                if (pho->fromHardProcessFinalState || pho->isPromptFinalState) {
                    TLorentzVector thisGenPhotonP4;
                    thisGenPhotonP4.SetPtEtaPhiM(pho->pt, pho->eta, pho->phi, 0.);
                    if (thisGenPhotonP4.DeltaR(photonP4) < 0.1) {
                        vetoDY = true;
                        break;
                    }
                }
            }
        }

    } // end dilepton selection
     
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

    if (!isData && genPhotons.size() > 0) {
        float min_phot_dr = 1000.;
        for (unsigned int i = 0; i < genPhotons.size(); i++) {
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

    outFile->cd(channel.c_str());
    outTrees[channel]->Fill();
    this->passedEvents++;
    return kTRUE;
}

void gbrTrainAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void gbrTrainAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void gbrTrainAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<gbrTrainAnalyzer> selector(new gbrTrainAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

float gbrTrainAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    //float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    float combIso = (mu->chHadIso03 + std::max(0.,(double)mu->neuHadIso03 + mu->gammaIso03 - 0.5*mu->puIso03));
    return combIso;
}

float gbrTrainAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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

float gbrTrainAnalyzer::GetPhotonIsolation(const baconhep::TPhoton* pho, const float rho)
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

void gbrTrainAnalyzer::EvalMuonEnergyResolution(std::map<string, float> mva_input_floats, std::map<string, int> mva_input_ints, 
                                                float &mean, float &sigma, float &alphaL, float &powerL, float &alphaR, float &powerR) 
{
    // Evaluates and returns the estimate of muon energy resolution function.
    // semi-parametric MVAs' inputs and outputs
    static RooRealVar* invar[99]; // [varnum]
    static RooAbsReal* mvaMean = NULL;
    static RooAbsReal* mvaSigma;
    static RooAbsReal* mvaAlphaL;
    static RooAbsReal* mvaAlphaR;

    // one-time MVA initialization: get trainings
    if (!mvaMean) {
        // current working directory
        TDirectory* wd = gDirectory;
        
        const std::string cmssw_base = getenv("CMSSW_BASE");
        string fileName;
        fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_training.root";
        TFile f(fileName.c_str());
        if (f.IsZombie())
            FATAL("TFile::Open() failed");

        // cd back into previous current working directory
        if (wd) wd->cd();
        else gDirectory = 0;

        RooWorkspace* ws = dynamic_cast<RooWorkspace*>(f.Get("ws_mva_muons"));
        if (!ws) FATAL("TFile::Get() failed");

        invar[0] = ws->var("var01");
        invar[1] = ws->var("var02");
        invar[2] = ws->var("var03");
        invar[3] = ws->var("var04");
        invar[4] = ws->var("var05");
        invar[5] = ws->var("var06");
        invar[6] = ws->var("var07");
        invar[7] = ws->var("var08");
        invar[8] = ws->var("var09");
        invar[9] = ws->var("var10");
        invar[10] = ws->var("var11");
        invar[11] = ws->var("var12");

        mvaMean = ws->function("limMean");
        mvaSigma = ws->function("limSigma");
        mvaAlphaL = ws->function("limAlphaL");
        mvaAlphaR = ws->function("limAlphaR");

    } // one-time initialization

    // load necessary tree branches
    float rho                      = mva_input_floats["rho"];
    float muEnergy                 = mva_input_floats["muEnergy"];
    float muEta                    = mva_input_floats["muEta"];
    float muChi2NDF                = mva_input_floats["muTkChi2"];
    Int_t muNumberOfValidTrkLayers = mva_input_ints["muNumberOfValidTrkLayers"];
    Int_t muNumberOfValidPixelHits = mva_input_ints["muNumberOfValidPixelHits"];
    Int_t muNumberOfValidMuonHits  = mva_input_ints["muNumberOfValidMuonHits"];
    Int_t muStations               = mva_input_ints["muStations"];
    float muPFIsoR04_CH            = mva_input_floats["muPFIsoR04_CH"];
    float muPFIsoR04_NH            = mva_input_floats["muPFIsoR04_NH"];
    float muPFIsoR04_Pho           = mva_input_floats["muPFIsoR04_Pho"];
    float muPFIsoR04_PU            = mva_input_floats["muPFIsoR04_PU"];

    // set input variables associated with the GBRLikelihood trainings
    *invar[0] = rho;
    *invar[1] = muEnergy;
    *invar[2] = muEta;
    *invar[3] = muChi2NDF;
    *invar[4] = muNumberOfValidTrkLayers;
    *invar[5] = muNumberOfValidPixelHits;
    *invar[6] = muNumberOfValidMuonHits;
    *invar[7] = muStations;
    *invar[8] = muPFIsoR04_CH;
    *invar[9] = muPFIsoR04_NH;
    *invar[10] = muPFIsoR04_Pho;
    *invar[11] = muPFIsoR04_PU;

    mean = mvaMean->getVal();
    sigma = mvaSigma->getVal();
    alphaL = mvaAlphaL->getVal();
    alphaR = mvaAlphaR->getVal();

    // NOTE: negative = infinite; powers were fixed at the training level
    powerL = -1;
    powerR = -1;

}

void gbrTrainAnalyzer::EvalElectronEnergyResolution(std::map<string, float> mva_inputs, float &mean, float &sigma, 
                                                    float &alphaL, float &powerL, float &alphaR, float &powerR) 
{
    // Evaluates and returns the estimate of muon energy resolution function.
    // semi-parametric MVAs' inputs and outputs
    static RooRealVar* invar[99]; // [varnum]
    static RooAbsReal* mvaMean = NULL;
    static RooAbsReal* mvaSigma;
    static RooAbsReal* mvaAlphaL;
    static RooAbsReal* mvaAlphaR;
    static RooAbsReal* mvaPowerL;
    static RooAbsReal* mvaPowerR;

    // one-time MVA initialization: get trainings
    if (!mvaMean) {
        // current working directory
        TDirectory* wd = gDirectory;
        
        const std::string cmssw_base = getenv("CMSSW_BASE");
        string fileName;
        fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_training.root";
        TFile f(fileName.c_str());
        if (f.IsZombie())
            FATAL("TFile::Open() failed");

        // cd back into previous current working directory
        if (wd) wd->cd();
        else gDirectory = 0;

        RooWorkspace* ws = dynamic_cast<RooWorkspace*>(f.Get("ws_mva_electrons"));
        if (!ws) FATAL("TFile::Get() failed");

        invar[0] = ws->var("var01");
        invar[1] = ws->var("var02");
        invar[2] = ws->var("var03");
        invar[3] = ws->var("var04");
        invar[4] = ws->var("var05");
        invar[5] = ws->var("var06");
        invar[6] = ws->var("var07");
        invar[7] = ws->var("var08");
        invar[8] = ws->var("var09");
        invar[9] = ws->var("var10");
        invar[10] = ws->var("var11");
        invar[11] = ws->var("var12");
        invar[12] = ws->var("var13");
        invar[13] = ws->var("var14");
        invar[14] = ws->var("var15");

        mvaMean = ws->function("limMean");
        mvaSigma = ws->function("limSigma");
        mvaAlphaL = ws->function("limAlphaL");
        mvaAlphaR = ws->function("limAlphaR");
        mvaPowerL = ws->function("limPowerL");
        mvaPowerR = ws->function("limPowerR");

    } // one-time initialization

    // load necessary tree branches
    float rho = mva_inputs["rho"];
    float elEnergy = mva_inputs["elEnergy"];
    float elScEta = mva_inputs["elScEta"];
    float elScPhi = mva_inputs["elScPhi"];
    float elR9 = mva_inputs["elR9"];
    float elE1x5OverE = mva_inputs["elE1x5OverE"];
    float elE2x5OverE = mva_inputs["elE2x5OverE"];
    float elE5x5OverE = mva_inputs["elE5x5OverE"];
    float elFBrem = mva_inputs["elFBrem"];
    float elHOverE = mva_inputs["elHOverE"];
    float elSigmaIEtaIEta = mva_inputs["elSigmaIEtaIEta"];
    float elPFIsoR04_CH = mva_inputs["elPFIsoR04_CH"];
    float elPFIsoR04_NH = mva_inputs["elPFIsoR04_NH"];
    float elPFIsoR04_Pho = mva_inputs["elPFIsoR04_Pho"];
    float elPFIsoR04_PU = mva_inputs["elPFIsoR04_PU"];

    // set input variables associated with the GBRLikelihood trainings
    *invar[0] = rho;
    *invar[1] = elEnergy;
    *invar[2] = elScEta;
    *invar[3] = elScPhi;
    *invar[4] = elR9;
    *invar[5] = elE1x5OverE;
    *invar[6] = elE2x5OverE;
    *invar[7] = elE5x5OverE;
    *invar[8] = elFBrem;
    *invar[9] = elHOverE;
    *invar[10] = elSigmaIEtaIEta;
    *invar[11] = elPFIsoR04_CH;
    *invar[12] = elPFIsoR04_NH;
    *invar[13] = elPFIsoR04_Pho;
    *invar[14] = elPFIsoR04_PU;

    mean = mvaMean->getVal();
    sigma = mvaSigma->getVal();
    alphaL = mvaAlphaL->getVal();
    alphaR = mvaAlphaR->getVal();
    powerL = mvaPowerL->getVal();
    powerR = mvaPowerR->getVal();
}

void gbrTrainAnalyzer::find_optimized(double* p, double &e1, double& e2)
{
   // Returns best-fitted (most probable) energies.

   static TF2* fun = NULL;

   // one-time initialization
   if (!fun) {
      fun = new TF2("fun", NegativeProbability, 1, 2000, 1, 2000, 16);
      fun->SetNpx(50);
      fun->SetNpy(50);
   }

   for (int i = 0; i < 16; i++)
      fun->SetParameter(i, p[i]);

   // limit the search range by +-3sigma regions
   fun->SetRange(p[2] * (1 - 3 * p[5]), p[3] * (1 - 3 * p[11]),
                 p[2] * (1 + 3 * p[5]), p[3] * (1 + 3 * p[11]));

   fun->GetMinimumXY(e1, e2);
}
