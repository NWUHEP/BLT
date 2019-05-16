#include "hzgAnalyzer.h"
#include <map>
#include <fstream>
#include <math.h>

#include <TSystem.h>
#include <TF2.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include "KinZfitter/KinZfitter/interface/KinZfitter.h"

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
        if (params->period == "2016") {
            triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
            triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
            triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
            triggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
        }
        else if (params->period == "2017") {
            triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*");
            triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*");
        }
        else if (params->period == "2018") {
            triggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*");
        }
    }

    else if (params->selection == "ee" || params->selection == "elelg") {
        if (params->period == "2016") {
            triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
            triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
        }
        else if (params->period == "2017") {
            triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
        }
        else if (params->period == "2018") {
            triggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
        }
    }

    else if (params->selection == "tautaug") { // select one muon plus one hadronic tau (for now)
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");
    }
        
    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask

    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    string jsonFileName;
    if (params->period == "2016") {
        jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/json/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
    }
    else if (params->period == "2017") {
        jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/json/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt";
    }
    else if (params->period == "2018") {
        jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/json/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt";
    }
    lumiMask.AddJSONFile(jsonFileName);

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
    outTree->Branch("nPV", &nPV);
    outTree->Branch("nPU", &nPU);
    outTree->Branch("nPartons", &nPartons);
    outTree->Branch("xPV", &xPV);
    outTree->Branch("yPV", &yPV);
    outTree->Branch("zPV", &zPV);

    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);
    outTree->Branch("metNC", &metNC);
    outTree->Branch("metPhiNC", &metPhiNC);
    outTree->Branch("ht", &ht);
    outTree->Branch("htPhi", &htPhi);
    outTree->Branch("htSum", &htSum);

    // weights
    outTree->Branch("genWeight", &genWeight);
    outTree->Branch("eventWeight", &eventWeight);
    outTree->Branch("puWeight", &puWeight);
    outTree->Branch("triggerWeight", &triggerWeight);
    outTree->Branch("elIDWeightOne", &elIDWeightOne);
    outTree->Branch("elIDWeightTwo", &elIDWeightTwo);
    outTree->Branch("elTrigWeightOne", &elTrigWeightOne);
    outTree->Branch("elTrigWeightTwo", &elTrigWeightTwo);
    outTree->Branch("muonIDWeightOne", &muonIDWeightOne);
    outTree->Branch("muonIDWeightTwo", &muonIDWeightTwo);
    outTree->Branch("muonTrigWeightOne", &muonTrigWeightOne);
    outTree->Branch("muonTrigWeightTwo", &muonTrigWeightTwo);
    outTree->Branch("photonIDWeight", &photonIDWeight);

    // leptons
    outTree->Branch("leptonOnePt", &leptonOnePt);
    outTree->Branch("leptonOneEta", &leptonOneEta);
    outTree->Branch("leptonOnePhi", &leptonOnePhi);
    outTree->Branch("leptonOnePtKin", &leptonOnePtKin);
    outTree->Branch("leptonOneIso", &leptonOneIso);
    outTree->Branch("leptonOneFlavor", &leptonOneFlavor);
    outTree->Branch("leptonOneMother", &leptonOneMother);
    outTree->Branch("leptonOneD0", &leptonOneD0);
    outTree->Branch("leptonOneDZ", &leptonOneDZ);
    outTree->Branch("leptonOneRecoWeight", &leptonOneRecoWeight);
    outTree->Branch("leptonOneECALDriven", &leptonOneECALDriven);

    outTree->Branch("leptonTwoPt", &leptonTwoPt);
    outTree->Branch("leptonTwoEta", &leptonTwoEta);
    outTree->Branch("leptonTwoPhi", &leptonTwoPhi);
    outTree->Branch("leptonTwoPtKin", &leptonTwoPtKin);
    outTree->Branch("leptonTwoIso", &leptonTwoIso);
    outTree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
    outTree->Branch("leptonTwoMother", &leptonTwoMother);
    outTree->Branch("leptonTwoD0", &leptonTwoD0);
    outTree->Branch("leptonTwoDZ", &leptonTwoDZ);
    outTree->Branch("leptonTwoRecoWeight", &leptonTwoRecoWeight);
    outTree->Branch("leptonTwoECALDriven", &leptonTwoECALDriven);
 
    outTree->Branch("tauDecayMode", &tauDecayMode);
    outTree->Branch("tauMVA", &tauMVA);

    outTree->Branch("isLeptonTag", &isLeptonTag);
    outTree->Branch("isDijetTag", &isDijetTag);
    outTree->Branch("isTightDijetTag", &isTightDijetTag);

    // photons
    outTree->Branch("photonOnePt", &photonOnePt);
    outTree->Branch("photonOneEta", &photonOneEta);
    outTree->Branch("photonOnePhi", &photonOnePhi);
    outTree->Branch("photonOneR9", &photonOneR9);
    outTree->Branch("photonOneMVA", &photonOneMVA);
    outTree->Branch("photonOneERes", &photonOneERes);
    outTree->Branch("passElectronVeto", &passElectronVeto);

    // jets
    outTree->Branch("jetOnePt", &jetOnePt);
    outTree->Branch("jetOneEta", &jetOneEta);
    outTree->Branch("jetOnePhi", &jetOnePhi);
    outTree->Branch("jetOneTag", &jetOneTag);
    outTree->Branch("jetTwoPt", &jetTwoPt);
    outTree->Branch("jetTwoEta", &jetTwoEta);
    outTree->Branch("jetTwoPhi", &jetTwoPhi);
    outTree->Branch("jetTwoTag", &jetTwoTag);

    // gen level objects 
    outTree->Branch("genLeptonOnePt", &genLeptonOnePt);
    outTree->Branch("genLeptonOneEta", &genLeptonOneEta);
    outTree->Branch("genLeptonOnePhi", &genLeptonOnePhi);
    outTree->Branch("genLeptonOneId", &genLeptonOneId);
    outTree->Branch("genLeptonTwoPt", &genLeptonTwoPt);
    outTree->Branch("genLeptonTwoEta", &genLeptonTwoEta);
    outTree->Branch("genLeptonTwoPhi", &genLeptonTwoPhi);
    outTree->Branch("genLeptonTwoId", &genLeptonTwoId);
    outTree->Branch("genPhotonPt", &genPhotonPt);
    outTree->Branch("genPhotonEta", &genPhotonEta);
    outTree->Branch("genPhotonPhi", &genPhotonPhi);
    outTree->Branch("genPhotonFHPFS", &genPhotonFHPFS);
    outTree->Branch("genPhotonIPFS", &genPhotonIPFS);
    outTree->Branch("vetoDY", &vetoDY);

    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nTaus", &nTaus);
    outTree->Branch("nPhotons", &nPhotons);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nFwdJets", &nFwdJets);
    outTree->Branch("nCentralJets", &nCentralJets);
    outTree->Branch("nBJets", &nBJets);
    
    // dilepton
    outTree->Branch("dileptonPt", &dileptonPt);
    outTree->Branch("dileptonEta", &dileptonEta);
    outTree->Branch("dileptonPhi", &dileptonPhi);
    outTree->Branch("dileptonM", &dileptonM);
    outTree->Branch("dileptonDEta", &dileptonDEta);
    outTree->Branch("dileptonDPhi", &dileptonDPhi);
    outTree->Branch("dileptonDR", &dileptonDR);
    outTree->Branch("dileptonMKin", &dileptonMKin);
 
    // dijet
    outTree->Branch("dijetPt", &dijetPt);
    outTree->Branch("dijetEta", &dijetEta);
    outTree->Branch("dijetPhi", &dijetPhi);
    outTree->Branch("dijetM", &dijetM);
    outTree->Branch("dijetDEta", &dijetDEta);
    outTree->Branch("dijetDPhi", &dijetDPhi);
    outTree->Branch("dijetDR", &dijetDR);

    // jet, lepton
    outTree->Branch("l1j1DEta", &l1j1DEta);
    outTree->Branch("l1j1DPhi", &l1j1DPhi);
    outTree->Branch("l1j1DR", &l1j1DR);
    outTree->Branch("l1j2DEta", &l1j2DEta);
    outTree->Branch("l1j2DPhi", &l1j2DPhi);
    outTree->Branch("l1j2DR", &l1j2DR);
    outTree->Branch("l2j1DEta", &l2j1DEta);
    outTree->Branch("l2j1DPhi", &l2j1DPhi);
    outTree->Branch("l2j1DR", &l2j1DR);
    outTree->Branch("l2j2DEta", &l2j2DEta);
    outTree->Branch("l2j2DPhi", &l2j2DPhi);
    outTree->Branch("l2j2DR", &l2j2DR);

    // jet, photon
    outTree->Branch("j1PhotonDEta", &j1PhotonDEta);
    outTree->Branch("j1PhotonDPhi", &j1PhotonDPhi);
    outTree->Branch("j1PhotonDR", &j1PhotonDR);
    outTree->Branch("j2PhotonDEta", &j2PhotonDEta);
    outTree->Branch("j2PhotonDPhi", &j2PhotonDPhi);
    outTree->Branch("j2PhotonDR", &j2PhotonDR);
    outTree->Branch("jPhotonDRMax", &jPhotonDRMax);
    outTree->Branch("jPhotonDRMin", &jPhotonDRMin);

    // three body
    outTree->Branch("llgPt", &llgPt);
    outTree->Branch("llgEta", &llgEta);
    outTree->Branch("llgPhi", &llgPhi);
    outTree->Branch("llgM", &llgM);
    outTree->Branch("llgPtOverM", &llgPtOverM);
    outTree->Branch("llgMKin", &llgMKin);
    outTree->Branch("l1PhotonDEta", &l1PhotonDEta);
    outTree->Branch("l1PhotonDPhi", &l1PhotonDPhi);
    outTree->Branch("l1PhotonDR", &l1PhotonDR);
    outTree->Branch("l2PhotonDEta", &l2PhotonDEta);
    outTree->Branch("l2PhotonDPhi", &l2PhotonDPhi);
    outTree->Branch("l2PhotonDR", &l2PhotonDR);
    outTree->Branch("lPhotonDRMax", &lPhotonDRMax);
    outTree->Branch("lPhotonDRMin", &lPhotonDRMin);
    outTree->Branch("dileptonPhotonDEta", &dileptonPhotonDEta);
    outTree->Branch("dileptonPhotonDPhi", &dileptonPhotonDPhi);
    outTree->Branch("dileptonPhotonDR", &dileptonPhotonDR);
    outTree->Branch("ptt", &ptt);
    outTree->Branch("zgBigTheta", &zgBigTheta);
    outTree->Branch("zgLittleTheta", &zgLittleTheta);
    outTree->Branch("zgPhi", &zgPhi);
    outTree->Branch("zgBigThetaMY", &zgBigThetaMY);
    outTree->Branch("zgLittleThetaMY", &zgLittleThetaMY);
    outTree->Branch("zgPhiMY", &zgPhiMY);
    outTree->Branch("zgBigThetaJames", &zgBigThetaJames);
    outTree->Branch("zgLittleThetaJames", &zgLittleThetaJames);
    outTree->Branch("zgPhiJames", &zgPhiJames);
    outTree->Branch("genBigTheta", &genBigTheta);
    outTree->Branch("genLittleTheta", &genLittleTheta);
    outTree->Branch("genPhi", &genPhi);

    // other 
    outTree->Branch("llgJJDEta", &llgJJDEta);
    outTree->Branch("llgJJDPhi", &llgJJDPhi);
    outTree->Branch("llgJJDR", &llgJJDR);
    outTree->Branch("zepp", &zepp);

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
        xPV = pv.X();
        yPV = pv.Y();
        zPV = pv.Z();
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
        particleSelector->ApplyMuonMomentumCorrection(muon, isData);
   
        // muons for analysis
        if (
                muon->pt > 5. 
                && fabs(muon->eta) < 2.4
                && particleSelector->PassMuonID(muon, "HZZ")
                && particleSelector->GetMuonIsolation(muon)/muon->pt < 0.35
           ) {
            muons.push_back(muon);
        }
                    
        // muons for jet veto
        if (
                muon->pt > 10.
                && fabs(muon->eta) < 2.4
                && particleSelector->PassMuonID(muon, "tight")
                && particleSelector->GetMuonIsolation(muon)/muon->pt < 0.15
           ) {
            TLorentzVector muonP4;
            muonP4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);
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
    
        if (
                electron->calibPt > 7.
                && fabs(electron->scEta) < 2.5
                && particleSelector->PassElectronMVA(electron, "HZZ")
                && particleSelector->GetElectronIsolation(electron, fInfo->rhoJet)/electron->calibPt < 0.35
                && fabs(electron->d0) < 0.5
                && fabs(electron->dz) < 1.0
                && fabs(electron->sip3d) < 4.0 
           ) {
            electrons.push_back(electron); // electrons for analysis
            TLorentzVector electronP4;
            electronP4.SetPtEtaPhiM(electron->calibPt, electron->eta, electron->phi, ELE_MASS);
            veto_electrons.push_back(electronP4); // electrons for jet veto
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

        // Prevent overlap of muons and electrons
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

        // apply tau energy scale correction 
        if (!isData) 
            particleSelector->ApplyTauEnergyScaleCorrection(tau);

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
        
        if (
                photon->pt > 10
                && fabs(photon->scEta) < 2.5 
                && (fabs(photon->scEta) <= 1.4442 || fabs(photon->scEta) >= 1.566)
                && particleSelector->PassPhotonMVA(photon, "loose")
                && photon->passElectronVeto
            ) {
            photons.push_back(photon);
            TLorentzVector photonP4;
            photonP4.SetPtEtaPhiM(photon->calibPt, photon->eta, photon->phi, 0.);
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

    nJets           = 0;
    nCentralJets    = 0;
    nFwdJets        = 0;
    nBJets          = 0;
    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

        double jec = particleSelector->JetCorrector(jet, "NONE");
        jet->pt = jet->ptRaw*jec;

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

        string jetIDName;
        if (params->period == "2016") {
            jetIDName == "loose";
        }
        else { // 2017 or 2018
            jetIDName == "tight";
        }

        if (
                jet->pt > 30 
                && fabs(jet->eta) < 4.7
                && particleSelector->PassJetID(jet, jetIDName)
                //&& !muOverlap 
                //&& !elOverlap
                //&& !phoOverlap
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
                    if (particleSelector->BTagModifier(jet, "MVAT", 0, 0)) { 
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

    std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);
    //std::sort(jets.begin(), jets.end(), sort_by_btag);

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

    if (params->selection == "mumu") {
        
        if (muons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        TLorentzVector muonOneP4, muonTwoP4;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, MUON_MASS);
        muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, MUON_MASS);

        if (muonOneP4.Pt() < 25.0) 
            return kTRUE;
        hTotalEvents->Fill(6);

        if (muonTwoP4.Pt() < 10.0)
            return kTRUE;
        hTotalEvents->Fill(7);

        if (muons[0]->q == muons[1]->q)
            return kTRUE;
        hTotalEvents->Fill(8);
        
        TLorentzVector dimuon = muonOneP4 + muonTwoP4;
        if (dimuon.M() < 80.0 || dimuon.M() > 100.0)
            return kTRUE;
        hTotalEvents->Fill(9);

        leptonOnePt     = muonOneP4.Pt();
        leptonOneEta    = muonOneP4.Eta();
        leptonOnePhi    = muonOneP4.Phi();
        leptonOneIso    = particleSelector->GetMuonIsolation(muons[0]);
        leptonOneFlavor = muons[0]->q*13;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;
          
        leptonTwoPt     = muonTwoP4.Pt();
        leptonTwoEta    = muonTwoP4.Eta();
        leptonTwoPhi    = muonTwoP4.Phi();
        leptonTwoIso    = particleSelector->GetMuonIsolation(muons[1]);
        leptonTwoFlavor = muons[1]->q*13;
        leptonTwoDZ     = muons[1]->dz;
        leptonTwoD0     = muons[1]->d0;

        if (!isData) {
            eventWeight *= weights->GetHZZMuonIDEff(*muons[0]); 
            eventWeight *= weights->GetHZZMuonIDEff(*muons[1]);
        }

    } // end mm selection
    
    else if (params->selection == "ee") {

        if (electrons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        TLorentzVector electronOneP4, electronTwoP4;
        electronOneP4.SetPtEtaPhiM(electrons[0]->calibPt, electrons[0]->eta, electrons[0]->phi, ELE_MASS);
        electronTwoP4.SetPtEtaPhiM(electrons[1]->calibPt, electrons[1]->eta, electrons[1]->phi, ELE_MASS);

        if (electronOneP4.Pt() < 25.0) 
            return kTRUE;
        hTotalEvents->Fill(6);

        if (electronTwoP4.Pt() < 15.0)
            return kTRUE;
        hTotalEvents->Fill(7);

        if (electrons[0]->q == electrons[1]->q)
            return kTRUE;
        hTotalEvents->Fill(8);

        TLorentzVector dielectron = electronOneP4 + electronTwoP4;
        if (dielectron.M() < 80.0 || dielectron.M() > 100.0)
            return kTRUE;
        hTotalEvents->Fill(9);   
        
        leptonOnePt     = electronOneP4.Pt();
        leptonOneEta    = electronOneP4.Eta();
        leptonOnePhi    = electronOneP4.Phi();
        leptonOneIso    = particleSelector->GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonOneFlavor = electrons[0]->q*11;
        leptonOneDZ     = electrons[0]->dz;
        leptonOneD0     = electrons[0]->d0;
            
        leptonTwoPt     = electronTwoP4.Pt();
        leptonTwoEta    = electronTwoP4.Eta();
        leptonTwoPhi    = electronTwoP4.Phi();
        leptonTwoIso    = particleSelector->GetElectronIsolation(electrons[1], fInfo->rhoJet);
        leptonTwoFlavor = electrons[1]->q*11;
        leptonTwoDZ     = electrons[1]->dz;
        leptonTwoD0     = electrons[1]->d0;
           
        if (!isData) {
            eventWeight *= weights->GetHZZElectronRecoIdEff(*electrons[0]); 
            eventWeight *= weights->GetHZZElectronRecoIdEff(*electrons[1]); 
        }

    } // end ee selection
    
    else if (params->selection == "tautaug") {

        if (muons.size() != 1) // avoid multi-muon events (DY contamination)
            return kTRUE;
        hTotalEvents->Fill(5);

        if (taus.size() < 1)
            return kTRUE;
        hTotalEvents->Fill(6);

        if (photons.size() < 1)
            return kTRUE;
        hTotalEvents->Fill(7);

        if (muons[0]->pt <= 25.) 
            return kTRUE;
        hTotalEvents->Fill(8);

        TLorentzVector muonP4;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, MUON_MASS);

        unsigned int tau_index = 0;
        for (unsigned int i = 0; i < taus.size(); ++i) {
            if (taus[i]->q != muons[0]->q && taus[i]->pt > 20.) {
                tau_index = i;
                break;
            }
        }

        if (taus[tau_index]->pt <= 20.)
            return kTRUE;
        hTotalEvents->Fill(9);

        TLorentzVector tauP4;
        tauP4.SetPtEtaPhiM(taus[tau_index]->pt, taus[tau_index]->eta, taus[tau_index]->phi, taus[tau_index]->m);

        if (photons[0]->calibPt <= 15.)
            return kTRUE;
        hTotalEvents->Fill(10);

        // event passed the selection; fill output variables
        TLorentzVector photonOneP4;
        photonOneP4.SetPtEtaPhiM(photons[0]->calibPt, photons[0]->eta, photons[0]->phi, 0.);
        photonOnePt = photonOneP4.Pt();
        photonOneEta = photonOneP4.Eta();
        photonOnePhi = photonOneP4.Phi();
        photonOneMVA = photons[0]->mvaSpring16;
        passElectronVeto = photons[0]->passElectronVeto;  
        if (!isData)
            photonOneR9 = weights->GetCorrectedPhotonR9(*photons[0]);
        else 
            photonOneR9 = photons[0]->r9;

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

        tauDecayMode    = taus[tau_index]->decaymode;
        tauMVA          = taus[tau_index]->rawIsoMVA3newDMwLT;

        if (muonP4.Pt() > tauP4.Pt()) {
            leptonOnePt     = muonP4.Pt();
            leptonOneEta     = muonP4.Eta();
            leptonOnePhi     = muonP4.Phi();
            leptonOneIso    = particleSelector->GetMuonIsolation(muons[0]);
            leptonOneFlavor = muons[0]->q*13;
            leptonOneDZ     = muons[0]->dz;
            leptonOneD0     = muons[0]->d0;

            leptonTwoPt     = tauP4.Pt();
            leptonTwoEta     = tauP4.Eta();
            leptonTwoPhi     = tauP4.Phi();
            leptonTwoIso    = 0.;
            leptonTwoFlavor = 15*taus[tau_index]->q;
            leptonTwoDZ     = taus[tau_index]->dzLeadChHad;
            leptonTwoD0     = taus[tau_index]->d0LeadChHad;
        }
        else {
            leptonOnePt     = tauP4.Pt();
            leptonOneEta     = tauP4.Eta();
            leptonOnePhi     = tauP4.Phi();
            leptonOneIso    = 0.;
            leptonOneFlavor = 15*taus[tau_index]->q;
            leptonOneDZ     = taus[tau_index]->dzLeadChHad;
            leptonOneD0     = taus[tau_index]->d0LeadChHad;

            leptonTwoPt     = muonP4.Pt();
            leptonTwoEta     = muonP4.Eta();
            leptonTwoPhi     = muonP4.Phi();
            leptonTwoIso    = particleSelector->GetMuonIsolation(muons[0]);
            leptonTwoFlavor = muons[0]->q*13;
            leptonTwoDZ     = muons[0]->dz;
            leptonTwoD0     = muons[0]->d0;
        }
    
        isDijetTag = false; // to ensure proper jet filling

        // MC event weights
        if (!isData) {

            eventWeight *= weights->GetHZZMuonIDEff(*muons[0]); 
            //eventWeight *= weights->GetMuonISOEff(muonP4);
            eventWeight *= weights->GetPhotonMVAIdEff(*photons[0]);
            eventWeight *= 0.95; // flat tau id scale factor

            float eff_data = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4).first; 
            float eff_mc = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4).second; 
            triggerWeight = eff_data/eff_mc;
            eventWeight *= triggerWeight;
            cout << "eventWeight = " << eventWeight << endl;
            
        }

    } // end tautaug selection

    else { // llg selection (combines mumug and elelg)
        assert(params->selection == "mumug" || params->selection == "elelg");

        unsigned int nLeptons = 0;
        if (params->selection == "mumug") 
            nLeptons = muons.size();
        else if (params->selection == "elelg")
            nLeptons = electrons.size();

        if (nLeptons < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        if (photons.size() < 1)
            return kTRUE;
        hTotalEvents->Fill(6);
           
        TLorentzVector leptonOneP4, leptonTwoP4;
        unsigned int leptonOneIndex = 0;
        unsigned int leptonTwoIndex = 1;
        bool hasValidPair = false;
        float zMassDiff = 100.;
        for (unsigned int i = 0; i < nLeptons; ++i) {
            for (unsigned int j = i+1; j < nLeptons; ++j) {
                TLorentzVector tempLeptonOne, tempLeptonTwo;
                if (params->selection == "mumug") {
                    tempLeptonOne.SetPtEtaPhiM(muons[i]->pt, muons[i]->eta, muons[i]->phi, MUON_MASS);
                    tempLeptonTwo.SetPtEtaPhiM(muons[j]->pt, muons[j]->eta, muons[j]->phi, MUON_MASS);
                }
                else if (params->selection == "elelg") {
                    tempLeptonOne.SetPtEtaPhiM(electrons[i]->calibPt, electrons[i]->eta, electrons[i]->phi, ELE_MASS);
                    tempLeptonTwo.SetPtEtaPhiM(electrons[j]->calibPt, electrons[j]->eta, electrons[j]->phi, ELE_MASS);
                }

                float thisMass = (tempLeptonOne + tempLeptonTwo).M();
                if (thisMass > 50.0) {
                    if (hasValidPair) {
                        if (fabs(thisMass - ZMASS) < zMassDiff) {
                            zMassDiff = fabs(thisMass - ZMASS);
                            leptonOneP4 = tempLeptonOne;
                            leptonTwoP4 = tempLeptonTwo;
                            leptonOneIndex = i;
                            leptonTwoIndex = j;
                        }
                    }
                    else {
                        zMassDiff = fabs(thisMass - ZMASS);
                        leptonOneP4 = tempLeptonOne;
                        leptonTwoP4 = tempLeptonTwo;
                        leptonOneIndex = i;
                        leptonTwoIndex = j;
                        hasValidPair = true;
                    }
                }

            }
        }

        if (!hasValidPair)
            return kTRUE;
        hTotalEvents->Fill(7);

        if (params->selection == "mumug") {
            if (leptonOneP4.Pt() <= 20.0) 
                return kTRUE;
            if (leptonTwoP4.Pt() <= 10.0)
                return kTRUE;
        }
        else if (params->selection == "elelg") {  
            if (leptonOneP4.Pt() <= 25.0)
                return kTRUE;
            if (leptonTwoP4.Pt() <= 15.0)
                return kTRUE;
        }

        TLorentzVector dileptonP4 = leptonOneP4 + leptonTwoP4;

        // L1EMTF cut 
        if (params->selection == "mumug") {
            if (
                fabs(leptonOneP4.DeltaPhi(leptonTwoP4)) < 70.0*(M_PI/180.0)
                && fabs(leptonOneP4.Eta()) > 1.2 
                && fabs(leptonTwoP4.Eta()) > 1.2
                && leptonOneP4.Eta()*leptonTwoP4.Eta() > 0
               )
                return kTRUE;
        }
        hTotalEvents->Fill(8); 

        bool hasValidPhoton = false;
        unsigned int photonIndex = 0;

        for (unsigned int i = 0; i < photons.size(); ++i) {
            TLorentzVector tempPhoton;
            TLorentzVector tempLLG;
            tempPhoton.SetPtEtaPhiM(photons[i]->calibPt, photons[i]->eta, photons[i]->phi, 0.);
            tempLLG = dileptonP4 + tempPhoton;
            float this_dr1 = leptonOneP4.DeltaR(tempPhoton);
            float this_dr2 = leptonTwoP4.DeltaR(tempPhoton);
            if (
                tempPhoton.Pt() > 15.0 &&
                tempPhoton.Et()/tempLLG.M() > (15.0/110.0) &&
                dileptonP4.M() + tempLLG.M() > 185.0 &&
                tempLLG.M() > 100. && tempLLG.M() < 180. &&
                this_dr1 > 0.4 && this_dr2 > 0.4
                ) {
                hasValidPhoton = true;
                photonIndex = i;
                break;
            }
        }

        if (!hasValidPhoton)
            return kTRUE;
        hTotalEvents->Fill(9);

        TLorentzVector photonOneP4;
        photonOneP4.SetPtEtaPhiM(photons[photonIndex]->calibPt, photons[photonIndex]->eta, photons[photonIndex]->phi, 0.);
        if (photonOneP4.Pt() < 15.0)
            return kTRUE;
        hTotalEvents->Fill(10);

        TLorentzVector llgP4 = dileptonP4 + photonOneP4;

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
                tempElectron.SetPtEtaPhiM(electrons[i]->calibPt, electrons[i]->eta, electrons[i]->phi, ELE_MASS);
                if (leptonOneP4.DeltaR(tempElectron) < 0.4 ||
                    leptonTwoP4.DeltaR(tempElectron) < 0.4 || 
                    photonOneP4.DeltaR(tempElectron) < 0.4)
                    continue;
                else
                    isLeptonTag = true;
            }
        }
            
        // checking for dijet tag
        isDijetTag = false;
        isTightDijetTag = false;
        unsigned int jetOneIndex = 0;
        unsigned int jetTwoIndex = 0;
        unsigned int purityLevel = 0;
        TLorentzVector jetOneP4, jetTwoP4;
        if (!isLeptonTag) {
            if (jets.size() > 1)  {
                for (unsigned int i = 0; i < jets.size(); ++i) {
                    for (unsigned int j = i+1; j < jets.size(); ++j) {
                        TLorentzVector tempJetOne;
                        TLorentzVector tempJetTwo;
                        tempJetOne.SetPtEtaPhiM(jets[i]->pt, jets[i]->eta, jets[i]->phi, jets[i]->mass);
                        tempJetTwo.SetPtEtaPhiM(jets[j]->pt, jets[j]->eta, jets[j]->phi, jets[j]->mass);
                        TLorentzVector tempDijet = tempJetOne + tempJetTwo;
                        float zeppen = llgP4.Eta() - (tempJetOne.Eta() + tempJetTwo.Eta())/2.;
                        if (    purityLevel < 5  
                            &&  tempJetOne.DeltaR(leptonOneP4) >= 0.4 && tempJetOne.DeltaR(leptonTwoP4) >= 0.4 
                            &&  tempJetTwo.DeltaR(leptonOneP4) >= 0.4 && tempJetTwo.DeltaR(leptonTwoP4) >= 0.4 
                            &&  tempJetOne.DeltaR(photonOneP4) >= 0.4 && tempJetTwo.DeltaR(photonOneP4) >= 0.4
                            &&  fabs(tempJetOne.Eta() - tempJetTwo.Eta()) >= 3.5 
                            &&  fabs(zeppen) <= 2.5 && tempDijet.M() >= 500.
                            &&  fabs(tempDijet.DeltaPhi(llgP4)) >= 2.4
                           ) {
                            isDijetTag = true;
                            isTightDijetTag = true;
                            jetOneIndex = i;
                            jetTwoIndex = j;
                            purityLevel = 5;
                            break;
                        }
                        else if (   purityLevel < 4 
                                &&  tempJetOne.DeltaR(leptonOneP4) >= 0.4 && tempJetOne.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetTwo.DeltaR(leptonOneP4) >= 0.4 && tempJetTwo.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetOne.DeltaR(photonOneP4) >= 0.4 && tempJetTwo.DeltaR(photonOneP4) >= 0.4
                                &&  fabs(tempJetOne.Eta() - tempJetTwo.Eta()) >= 3.5 
                                &&  fabs(zeppen) <= 2.5
                                &&  fabs(tempDijet.DeltaPhi(llgP4)) >= 2.4
                                ) {
                            isDijetTag = true;
                            jetOneIndex = i;
                            jetTwoIndex = j;
                            purityLevel = 4;
                        }
                        else if (   purityLevel < 3 
                                &&  tempJetOne.DeltaR(leptonOneP4) >= 0.4 && tempJetOne.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetTwo.DeltaR(leptonOneP4) >= 0.4 && tempJetTwo.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetOne.DeltaR(photonOneP4) >= 0.4 && tempJetTwo.DeltaR(photonOneP4) >= 0.4
                                &&  fabs(tempJetOne.Eta() - tempJetTwo.Eta()) >= 3.5 
                                &&  fabs(tempDijet.DeltaPhi(llgP4)) >= 2.4
                                ) {
                            isDijetTag = true;
                            jetOneIndex = i;
                            jetTwoIndex = j;
                            purityLevel = 3;
                        }
                        else if (   purityLevel < 2
                                &&  tempJetOne.DeltaR(leptonOneP4) >= 0.4 && tempJetOne.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetTwo.DeltaR(leptonOneP4) >= 0.4 && tempJetTwo.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetOne.DeltaR(photonOneP4) >= 0.4 && tempJetTwo.DeltaR(photonOneP4) >= 0.4
                                &&  fabs(tempDijet.DeltaPhi(llgP4)) >= 2.4
                                ) {
                            isDijetTag = true;
                            jetOneIndex = i;
                            jetTwoIndex = j;
                            purityLevel = 2;
                        }
                        else if (   purityLevel < 1
                                &&  tempJetOne.DeltaR(leptonOneP4) >= 0.4 && tempJetOne.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetTwo.DeltaR(leptonOneP4) >= 0.4 && tempJetTwo.DeltaR(leptonTwoP4) >= 0.4 
                                &&  tempJetOne.DeltaR(photonOneP4) >= 0.4 && tempJetTwo.DeltaR(photonOneP4) >= 0.4
                                ) {
                            isDijetTag = true;
                            jetOneIndex = i;
                            jetTwoIndex = j;
                            purityLevel = 1;
                        }
                                
                    }
                }
            }
        }

        if (isDijetTag) {
            jetOneP4.SetPtEtaPhiM(jets[jetOneIndex]->pt, jets[jetOneIndex]->eta, jets[jetOneIndex]->phi, jets[jetOneIndex]->mass);
            jetTwoP4.SetPtEtaPhiM(jets[jetTwoIndex]->pt, jets[jetTwoIndex]->eta, jets[jetTwoIndex]->phi, jets[jetTwoIndex]->mass);
            jetOnePt = jetOneP4.Pt();
            jetOneEta = jetOneP4.Eta();
            jetOnePhi = jetOneP4.Phi();
            jetOneM   = jetOneP4.M();
            
            jetTwoPt = jetTwoP4.Pt();
            jetTwoEta = jetTwoP4.Eta();
            jetTwoPhi = jetTwoP4.Phi();
            jetTwoM   = jetTwoP4.M();

            TLorentzVector dijet = jetOneP4 + jetTwoP4;
            dijetPt = dijet.Pt();
            dijetEta = dijet.Eta();
            dijetPhi = dijet.Phi();
            dijetM = dijet.M();
            dijetDEta = fabs(jetOneP4.Eta() - jetTwoP4.Eta());
            dijetDPhi = fabs(jetOneP4.DeltaPhi(jetTwoP4));
            dijetDR = jetOneP4.DeltaR(jetTwoP4);

            l1j1DEta = fabs(leptonOneP4.Eta() - jetOneP4.Eta());
            l1j1DPhi = fabs(leptonOneP4.DeltaPhi(jetOneP4));
            l1j1DR = leptonOneP4.DeltaR(jetOneP4);
            l1j2DEta = fabs(leptonOneP4.Eta() - jetTwoP4.Eta());
            l1j2DPhi = fabs(leptonOneP4.DeltaPhi(jetTwoP4));
            l1j2DR = leptonOneP4.DeltaR(jetTwoP4);
            l2j1DEta = fabs(leptonTwoP4.Eta() - jetOneP4.Eta());
            l2j1DPhi = fabs(leptonTwoP4.DeltaPhi(jetOneP4));
            l2j1DR = leptonTwoP4.DeltaR(jetOneP4);
            l2j2DEta = fabs(leptonTwoP4.Eta() - jetTwoP4.Eta());
            l2j2DPhi = fabs(leptonTwoP4.DeltaPhi(jetTwoP4));
            l2j2DR = leptonTwoP4.DeltaR(jetTwoP4);

            j1PhotonDEta = fabs(jetOneP4.Eta() - photonOneP4.Eta());
            j1PhotonDPhi = fabs(jetOneP4.DeltaPhi(photonOneP4));
            j1PhotonDR = jetOneP4.DeltaR(photonOneP4);
            j2PhotonDEta = fabs(jetTwoP4.Eta() - photonOneP4.Eta());
            j2PhotonDPhi = fabs(jetTwoP4.DeltaPhi(photonOneP4));
            j2PhotonDR = jetTwoP4.DeltaR(photonOneP4);
        
            if (j1PhotonDR > j2PhotonDR) {
                jPhotonDRMax = j1PhotonDR;
                jPhotonDRMin = j2PhotonDR; 
            }
            else {
                jPhotonDRMax = j2PhotonDR;
                jPhotonDRMin = j1PhotonDR;
            }
            
            zepp = llgP4.Eta() - (jetOneP4.Eta() + jetTwoP4.Eta())/2.;
            llgJJDEta = fabs(llgP4.Eta() - dijet.Eta());
            llgJJDPhi = fabs(llgP4.DeltaPhi(dijet));
            llgJJDR = llgP4.DeltaR(dijet);
        }
     
        leptonOnePt     = leptonOneP4.Pt();
        leptonOneEta     = leptonOneP4.Eta();
        leptonOnePhi     = leptonOneP4.Phi();
        leptonTwoPt     = leptonTwoP4.Pt();
        leptonTwoEta     = leptonTwoP4.Eta();
        leptonTwoPhi     = leptonTwoP4.Phi();

        if (params->selection == "mumug") {
            leptonOneIso    = particleSelector->GetMuonIsolation(muons[leptonOneIndex]);
            leptonOneFlavor = muons[leptonOneIndex]->q*13;
            leptonOneDZ     = muons[leptonOneIndex]->dz;
            leptonOneD0     = muons[leptonOneIndex]->d0;     
            leptonTwoIso    = particleSelector->GetMuonIsolation(muons[leptonTwoIndex]);
            leptonTwoFlavor = muons[leptonTwoIndex]->q*13;
            leptonTwoDZ     = muons[leptonTwoIndex]->dz;
            leptonTwoD0     = muons[leptonTwoIndex]->d0;
        }
        else if (params->selection == "elelg") {
            leptonOneIso    = particleSelector->GetElectronIsolation(electrons[leptonOneIndex], fInfo->rhoJet);
            leptonOneFlavor = electrons[leptonOneIndex]->q*11;
            leptonOneDZ     = electrons[leptonOneIndex]->dz;
            leptonOneD0     = electrons[leptonOneIndex]->d0; 
            leptonTwoIso    = particleSelector->GetElectronIsolation(electrons[leptonTwoIndex], fInfo->rhoJet);
            leptonTwoFlavor = electrons[leptonTwoIndex]->q*11;
            leptonTwoDZ     = electrons[leptonTwoIndex]->dz;
            leptonTwoD0     = electrons[leptonTwoIndex]->d0;
        }

        photonOnePt  = photonOneP4.Pt();
        photonOneEta  = photonOneP4.Eta();
        photonOnePhi  = photonOneP4.Phi();
        photonOneMVA = photons[photonIndex]->mvaSpring16;
        photonOneERes = photons[photonIndex]->eRes;
        passElectronVeto = photons[photonIndex]->passElectronVeto;  

        if (!isData)
            photonOneR9 = weights->GetCorrectedPhotonR9(*photons[photonIndex]);
        else 
            photonOneR9 = photons[photonIndex]->r9;

        // kinematic fit
        KinZfitter* kinZfitter = new KinZfitter(isData);
        std::vector<TObject *> selectedLeptons;
        if (params->selection == "mumug") {
            selectedLeptons.push_back(muons[leptonOneIndex]);
            selectedLeptons.push_back(muons[leptonTwoIndex]);
        }
        else if (params->selection == "elelg") {
            selectedLeptons.push_back(electrons[leptonOneIndex]);
            selectedLeptons.push_back(electrons[leptonTwoIndex]);
        }

        std::map<unsigned int, TLorentzVector> selectedFsrMap;
        kinZfitter->Setup(selectedLeptons, selectedFsrMap);
        kinZfitter->KinRefitZ();
        std::vector<TLorentzVector> refitP4s;
        refitP4s = kinZfitter->GetRefitP4s();
        delete kinZfitter;

        TLorentzVector leptonOneP4KinFit = refitP4s[0];
        TLorentzVector leptonTwoP4KinFit = refitP4s[1];

        leptonOnePtKin = leptonOneP4KinFit.Pt();
        leptonTwoPtKin = leptonTwoP4KinFit.Pt();

        dileptonPt = dileptonP4.Pt();
        dileptonEta = dileptonP4.Eta();
        dileptonPhi = dileptonP4.Phi();
        dileptonM = dileptonP4.M();
        dileptonMKin = (leptonOneP4KinFit + leptonTwoP4KinFit).M();
        dileptonDEta = fabs(leptonOneP4.Eta() - leptonTwoP4.Eta());
        dileptonDPhi = fabs(leptonOneP4.DeltaPhi(leptonTwoP4));
        dileptonDR = leptonOneP4.DeltaR(leptonTwoP4);

        llgPt = llgP4.Pt();
        llgEta = llgP4.Eta();
        llgPhi = llgP4.Phi();
        llgM = llgP4.M();
        llgPtOverM = llgP4.Pt()/llgP4.M();
        llgMKin = (leptonOneP4KinFit + leptonTwoP4KinFit + photonOneP4).M();
        
        l1PhotonDEta = fabs(leptonOneP4.Eta() - photonOneP4.Eta());
        l1PhotonDPhi = fabs(leptonOneP4.DeltaPhi(photonOneP4));
        l1PhotonDR = leptonOneP4.DeltaR(photonOneP4);
        l2PhotonDEta = fabs(leptonTwoP4.Eta() - photonOneP4.Eta());
        l2PhotonDPhi = fabs(leptonTwoP4.DeltaPhi(photonOneP4));
        l2PhotonDR = leptonTwoP4.DeltaR(photonOneP4);
        
        if (l1PhotonDR > l2PhotonDR) {
            lPhotonDRMax = l1PhotonDR;
            lPhotonDRMin = l2PhotonDR; 
        }
        else {
            lPhotonDRMax = l2PhotonDR;
            lPhotonDRMin = l1PhotonDR;
        }

        dileptonPhotonDEta = fabs(dileptonP4.Eta() - photonOneP4.Eta());
        dileptonPhotonDPhi = fabs(dileptonP4.DeltaPhi(photonOneP4));
        dileptonPhotonDR = dileptonP4.DeltaR(photonOneP4);
        ptt = 2*fabs(dileptonP4.Px()*photonOneP4.Py() - photonOneP4.Px()*dileptonP4.Py())/llgP4.Pt();

        // calculate angles like Brian
        TVector3 Xframe = llgP4.BoostVector();
        TVector3 Z1frame = dileptonP4.BoostVector();

        // "partons"
        TLorentzVector kq, kqbar, veckq_in_Xframe, veckqbar_in_Xframe;
        kq.SetPxPyPzE(0., 0., (llgP4.E() + llgP4.Pz())/2., (llgP4.E() + llgP4.Pz())/2.);
        kqbar.SetPxPyPzE(0., 0., (llgP4.Pz() - llgP4.E())/2., (llgP4.E() - llgP4.Pz())/2.);
        veckq_in_Xframe = kq;
        veckqbar_in_Xframe = kqbar;
        veckq_in_Xframe.Boost(-1*Xframe);
        veckqbar_in_Xframe.Boost(-1*Xframe);
   
        // Z vectors
        TLorentzVector vecz_in_Xframe = dileptonP4;
        TLorentzVector vecg_in_Xframe = photonOneP4;
        TLorentzVector vecz_in_Z1frame = dileptonP4;
        vecz_in_Xframe.Boost(-1*Xframe);
        vecg_in_Xframe.Boost(-1*Xframe);
        vecz_in_Z1frame.Boost(-1*Z1frame);

        // coord system in the CM frame
        TVector3 uz_in_Xframe = vecz_in_Xframe.Vect().Unit();
        TVector3 uy_in_Xframe = (veckq_in_Xframe.Vect().Unit().Cross(uz_in_Xframe.Unit())).Unit();
        TVector3 ux_in_Xframe = (uy_in_Xframe.Unit().Cross(uz_in_Xframe.Unit())).Unit();
        TRotation rotation;
        rotation = rotation.RotateAxes(ux_in_Xframe, uy_in_Xframe, uz_in_Xframe).Inverse();

        // for going to the Z frames from the CM frame, boost after transform
        TLorentzVector vecz_in_Xframe_newcoords = vecz_in_Xframe;
        vecz_in_Xframe_newcoords.Transform(rotation);
        TVector3 Z1frame_from_Xframe_newcoords = vecz_in_Xframe_newcoords.BoostVector();

        // define the positive and negative leptons
        TLorentzVector l_minus_james, l_plus_james; 
        if (leptonOneFlavor > 0) {
            l_minus_james = leptonOneP4;
            l_plus_james = leptonTwoP4;
        }
        else {
            l_minus_james = leptonTwoP4;
            l_plus_james = leptonOneP4;
        }
       
        // little theta, phi in Z1 frame; first boost to CM, then redefine coords
        TLorentzVector veclm_in_Z1frame = l_minus_james;
        TLorentzVector veclp_in_Z1frame = l_plus_james;
        veclm_in_Z1frame.Boost(-1*Xframe);
        veclm_in_Z1frame.Transform(rotation);
        veclp_in_Z1frame.Boost(-1*Xframe);
        veclp_in_Z1frame.Transform(rotation);

        // then boost to Z1
        veclm_in_Z1frame.Boost(-1*Z1frame_from_Xframe_newcoords);
        veclp_in_Z1frame.Boost(-1*Z1frame_from_Xframe_newcoords);
        
        // now get angles
        zgPhiJames = veclm_in_Z1frame.Phi();
        zgLittleThetaJames = veclm_in_Z1frame.CosTheta();

        //if (zgPhiJames < 0) 
        //    zgPhiJames += 2*M_PI;

        // Big Theta in X frame
        TLorentzVector veczg_in_Xframe = llgP4;
        veczg_in_Xframe.Transform(rotation);

        TLorentzVector veczg_in_Xframe_newcoords = llgP4;
        veczg_in_Xframe_newcoords.Transform(rotation);
        //zgBigThetaJames = (-1*veczg_in_Xframe_newcoords.Vect()).CosTheta();
        zgBigThetaJames = (veczg_in_Xframe_newcoords.Vect()).CosTheta();

        /*cout << "rotation matrix Brian" << endl;
        cout << rotation.XX() << " " << rotation.XY() << " " << rotation.XZ() << endl;
        cout << rotation.YX() << " " << rotation.YY() << " " << rotation.YZ() << endl;
        cout << rotation.ZX() << " " << rotation.ZY() << " " << rotation.ZZ() << endl;*/

        // calculate angles like Ming-Yan but with l_minus and l_plus
        TLorentzVector l_minus, l_plus; 

        if (leptonOneFlavor > 0) {
            l_minus.SetPtEtaPhiM(leptonOneP4.Pt(), leptonOneP4.Eta(), leptonOneP4.Phi(), leptonOneP4.M());
            l_plus.SetPtEtaPhiM(leptonTwoP4.Pt(), leptonTwoP4.Eta(), leptonTwoP4.Phi(), leptonTwoP4.M());
        }
        else {
            l_minus.SetPtEtaPhiM(leptonTwoP4.Pt(), leptonTwoP4.Eta(), leptonTwoP4.Phi(), leptonTwoP4.M());
            l_plus.SetPtEtaPhiM(leptonOneP4.Pt(), leptonOneP4.Eta(), leptonOneP4.Phi(), leptonOneP4.M());
        }
        
        TVector3 llgFrame = -1*llgP4.BoostVector();
        dileptonP4.Boost(llgFrame);
        l_minus.Boost(llgFrame);
        l_minus.Boost(-dileptonP4.BoostVector());
        zgLittleTheta = cos(dileptonP4.Angle(l_minus.Vect()));
        zgBigTheta = cos(dileptonP4.Angle(llgP4.Vect()));
       
        TVector3 ppAxis(0, 0, 1);
        //TVector3 ppAxis = veckq_in_Xframe.Vect().Unit();
        TVector3 zAxis = dileptonP4.Vect().Unit();
        TVector3 yAxis = ppAxis.Cross(zAxis.Unit()).Unit();
        TVector3 xAxis = (yAxis.Unit().Cross(zAxis.Unit())).Unit();

        TRotation rot;
        rot = rot.RotateAxes(xAxis, yAxis, zAxis).Inverse();

        dileptonP4.Transform(rot);
        l_minus.Transform(rot);
        zgPhi = l_minus.Phi();
        
        /*cout << "rotation matrix Ming-Yan flavor" << endl;
        cout << rot.XX() << " " << rot.XY() << " " << rot.XZ() << endl;
        cout << rot.YX() << " " << rot.YY() << " " << rot.YZ() << endl;
        cout << rot.ZX() << " " << rot.ZY() << " " << rot.ZZ() << endl;*/
        
        // the way MY does it (I think wrong)
        TLorentzVector lep1, lep2;
        lep1.SetPtEtaPhiM(leptonOneP4.Pt(), leptonOneP4.Eta(), leptonOneP4.Phi(), leptonOneP4.M());
        lep2.SetPtEtaPhiM(leptonTwoP4.Pt(), leptonTwoP4.Eta(), leptonTwoP4.Phi(), leptonTwoP4.M());
        dileptonP4 = lep1 + lep2;
        
        llgFrame = -1*llgP4.BoostVector();
        dileptonP4.Boost(llgFrame);
        lep1.Boost(llgFrame);
        lep1.Boost(-dileptonP4.BoostVector());
        zgLittleThetaMY = cos(dileptonP4.Angle(lep1.Vect()));
        zgBigThetaMY = cos(dileptonP4.Angle(llgP4.Vect()));
        
        zAxis = dileptonP4.Vect().Unit();
        yAxis = ppAxis.Cross(zAxis.Unit()).Unit();
        xAxis = (yAxis.Unit().Cross(zAxis.Unit())).Unit();

        TRotation rot_my;
        rot_my = rot_my.RotateAxes(xAxis, yAxis, zAxis).Inverse();

        dileptonP4.Transform(rot_my);
        lep1.Transform(rot_my);
        zgPhiMY = lep1.Phi(); 

        /*cout << "values for cos big theta: James = " << zgBigThetaJames << ",  Ming-Yan + flavor = " << zgBigTheta << ", Ming-Yan + leading = " << zgBigThetaMY << endl; 
        cout << "values for cos little theta: James = " << zgLittleThetaJames << ",  Ming-Yan + flavor = " << zgLittleTheta << ", Ming-Yan + leading = " << zgLittleThetaMY << endl; 
        cout << "values for phi: James = " << zgPhiJames << ",  Ming-Yan + flavor = " << zgPhi << ", Ming-Yan + leading = " << zgPhiMY << endl; */
            
        if (!isData) {
            if (params->selection == "mumug") {
                muonIDWeightOne = weights->GetHZZMuonIDEff(*muons[leptonOneIndex]); 
                muonIDWeightTwo = weights->GetHZZMuonIDEff(*muons[leptonTwoIndex]);
                eventWeight *= muonIDWeightOne;
                eventWeight *= muonIDWeightTwo;

                float sf11 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg1", *muons[leptonOneIndex]);
                float sf12 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg1", *muons[leptonTwoIndex]);
                float sf21 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg2", *muons[leptonOneIndex]);
                float sf22 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg2", *muons[leptonTwoIndex]);
                
                if (leptonTwoPt < 20.) {
                    muonTrigWeightOne = sf11;
                    muonTrigWeightTwo = sf22;
                }
                else {
                    float prod1 = sf11*sf22;
                    float prod2 = sf21*sf12;
                    if (prod1 > prod2) {
                        muonTrigWeightOne = sf11;
                        muonTrigWeightTwo = sf22;
                    }
                    else {
                        muonTrigWeightOne = sf21;
                        muonTrigWeightTwo = sf12;
                    }
                }

                triggerWeight = muonTrigWeightOne*muonTrigWeightTwo;
            }
            else if (params->selection == "elelg") {
                elIDWeightOne = weights->GetHZZElectronRecoIdEff(*electrons[leptonOneIndex]); 
                elIDWeightTwo = weights->GetHZZElectronRecoIdEff(*electrons[leptonTwoIndex]); 
                eventWeight *= elIDWeightOne;
                eventWeight *= elIDWeightTwo;
                
                float sf11 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg1", *electrons[leptonOneIndex]);
                float sf22 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg2", *electrons[leptonTwoIndex]);

                elTrigWeightOne = sf11;
                elTrigWeightTwo = sf22;

                triggerWeight = elTrigWeightOne*elTrigWeightTwo;
            }

            eventWeight *= triggerWeight;
           
            photonIDWeight = weights->GetPhotonMVAIdEff(*photons[photonIndex]); 
            eventWeight *= photonIDWeight;
        }

    } // end llg selection
     
    ///////////////////
    // Fill jet info //
    ///////////////////
    
    if (!isDijetTag) {
        if (jets.size() > 0) {
            jetOnePt = jets[0]->pt;
            jetOneEta = jets[0]->eta;
            jetOnePhi = jets[0]->phi;
            jetOneM = jets[0]->mass;
            jetOneTag    = jets[0]->csv;
        } else {
            jetOnePt = 0.;
            jetOneEta = 0.;
            jetOnePhi = 0.;
            jetOneM = 0.;
            jetOneTag    = 0.;
        }

        if (jets.size() > 1) {
            jetTwoPt = jets[1]->pt;
            jetTwoEta = jets[1]->eta;
            jetTwoPhi = jets[1]->phi;
            jetTwoM = jets[1]->mass;
            jetTwoTag    = jets[1]->csv;
        } else {
            jetTwoPt = 0.;
            jetTwoEta = 0.;
            jetTwoPhi = 0.;
            jetTwoM = 0.;
            jetTwoTag    = 0.;
        }
    }
 
    if (!isData && genLeptons.size() == 2) {
        genLeptonOneId = genLeptons[0]->pdgId;
        genLeptonOnePt = genLeptons[0]->pt;
        genLeptonOneEta = genLeptons[0]->eta;
        genLeptonOnePhi = genLeptons[0]->phi;
        genLeptonTwoId = genLeptons[1]->pdgId;
        genLeptonTwoPt = genLeptons[1]->pt;
        genLeptonTwoEta = genLeptons[1]->eta;
        genLeptonTwoPhi = genLeptons[1]->phi;   
    }
    else {
        genLeptonOnePt = 0.;
        genLeptonOneEta = 0.;
        genLeptonOnePhi = 0.;
        genLeptonTwoPt = 0.;
        genLeptonTwoEta = 0.;
        genLeptonTwoPhi = 0.;
    }

    TLorentzVector genPhotonP4;
    if (!isData && genPhotons.size() > 0) {
        TLorentzVector photonOneP4;
        photonOneP4.SetPtEtaPhiM(photonOnePt, photonOneEta, photonOnePhi, 0.);
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
        genPhotonPt = genPhotonP4.Pt();
        genPhotonEta = genPhotonP4.Eta();
        genPhotonPhi = genPhotonP4.Phi();
    }
    else {
        genPhotonPt = 0.;
        genPhotonEta = 0.;
        genPhotonPhi = 0.;
    }
        
    // gen angles
    if (!isData && genLeptons.size() == 2 && genPhotons.size() > 0) {        
        TLorentzVector l_minus, l_plus; 
        if (genLeptonOneId > 0) {
            l_minus.SetPtEtaPhiM(genLeptonOnePt, genLeptonOneEta, genLeptonOnePhi, genLeptons[0]->mass);
            l_plus.SetPtEtaPhiM(genLeptonTwoPt, genLeptonTwoEta, genLeptonTwoPhi, genLeptons[1]->mass);
        }
        else {
            l_minus.SetPtEtaPhiM(genLeptonTwoPt, genLeptonTwoEta, genLeptonTwoPhi, genLeptons[1]->mass);
            l_plus.SetPtEtaPhiM(genLeptonOnePt, genLeptonOneEta, genLeptonOnePhi, genLeptons[0]->mass);
        }
     
        TLorentzVector dileptonP4Gen, llgP4Gen;
        dileptonP4Gen = l_minus + l_plus;
        llgP4Gen = dileptonP4Gen + genPhotonP4;

        TVector3 llgFrame = -1*llgP4Gen.BoostVector();
        dileptonP4Gen.Boost(llgFrame);
        l_minus.Boost(llgFrame);
        l_minus.Boost(dileptonP4Gen.BoostVector());
        genLittleTheta = cos(dileptonP4Gen.Angle(l_minus.Vect()));
        genBigTheta = cos(dileptonP4Gen.Angle(llgP4Gen.Vect()));
        
        TVector3 ppAxis(0, 0, 1);
        TVector3 zAxis = dileptonP4Gen.Vect().Unit();
        TVector3 yAxis = ppAxis.Cross(zAxis.Unit()).Unit();
        TVector3 xAxis = (yAxis.Unit().Cross(zAxis.Unit())).Unit();

        TRotation rotation;
        rotation = rotation.RotateAxes(xAxis, yAxis, zAxis).Inverse();

        dileptonP4Gen.Transform(rotation);
        l_minus.Transform(rotation);
        genPhi = l_minus.Phi();
    }
    else {
        genLittleTheta = -9.;
        genBigTheta = -9.;
        genPhi = -9.;
    }

    outTree->Fill();
    this->passedEvents++;

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
