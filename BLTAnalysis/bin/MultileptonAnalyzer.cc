#include "MultileptonAnalyzer.h"
#include <map>

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 

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

    // leptons
    outTree->Branch("leptonOneP4", &leptonOneP4);
    outTree->Branch("leptonOneIso", &leptonOneIso);
    outTree->Branch("leptonOneQ", &leptonOneQ);
    outTree->Branch("leptonOneFlavor", &leptonOneFlavor);
    outTree->Branch("leptonOneTrigger", &leptonOneTrigger);
    outTree->Branch("leptonOneD0", &leptonOneD0);
    outTree->Branch("leptonOneDZ", &leptonOneDZ);

    outTree->Branch("leptonTwoP4", &leptonTwoP4);
    outTree->Branch("leptonTwoIso", &leptonTwoIso);
    outTree->Branch("leptonTwoQ", &leptonTwoQ);
    outTree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
    outTree->Branch("leptonTwoTrigger", &leptonTwoTrigger);
    outTree->Branch("leptonTwoD0", &leptonTwoD0);
    outTree->Branch("leptonTwoDZ", &leptonTwoDZ);
    
    outTree->Branch("leptonThreeP4", &leptonThreeP4);
    outTree->Branch("leptonThreeIso", &leptonThreeIso);
    outTree->Branch("leptonThreeQ", &leptonThreeQ);
    outTree->Branch("leptonThreeFlavor", &leptonThreeFlavor);
    outTree->Branch("leptonThreeTrigger", &leptonThreeTrigger);

    outTree->Branch("leptonFourP4", &leptonFourP4);
    outTree->Branch("leptonFourIso", &leptonFourIso);
    outTree->Branch("leptonFourQ", &leptonFourQ);
    outTree->Branch("leptonFourFlavor", &leptonFourFlavor);
    outTree->Branch("leptonFourTrigger", &leptonFourTrigger);

    // lepton vertices
    outTree->Branch("dileptonVertexOne", &dileptonVertexOne);
    outTree->Branch("dileptonVertexErrOne", &dileptonVertexErrOne);
    outTree->Branch("dileptonVertexChi2One", &dileptonVertexChi2One);
    outTree->Branch("dileptonVertexDOFOne", &dileptonVertexDOFOne);

    outTree->Branch("dileptonVertexTwo", &dileptonVertexTwo);
    outTree->Branch("dileptonVertexErrTwo", &dileptonVertexErrTwo);
    outTree->Branch("dileptonVertexChi2Two", &dileptonVertexChi2Two);
    outTree->Branch("dileptonVertexDOFTwo", &dileptonVertexDOFTwo);

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
        copy_p4(muons[i], MUON_MASS, muonP4);
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
            if (muonP4.Pt() > 20) {
                veto_muons.push_back(muonP4);
            }
        }
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);
    return kTRUE;

    /* ELECTRONS */
    std::vector<TElectron*> electrons;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        if (
                electron->pt > 20 
                && fabs(electron->eta) < 2.5
                && particleSelector->PassElectronID(electron, cuts->tightElID)
                && particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
           ) {
            electrons.push_back(electron);
        }
    }

    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);

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
        TLorentzVector vJet; 
        vJet.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
        bool muOverlap = false;
        for (const auto& mu: veto_muons) {
            if (vJet.DeltaR(mu) < 0.5) {
                muOverlap = true;
                break;
            }
        }
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


    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    nMuons     = muons.size();
    nElectrons = electrons.size();

    if (params->selection == "mumu") {
        if (muons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        if (muons[0]->pt < 25)
            return kTRUE;
        hTotalEvents->Fill(6);

        // Find leading positive and negatively charged muons and convert to TLorentzVectors
        TLorentzVector muonOneP4, muonTwoP4, dimuonP4;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 105.2);
        muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, 105.2);
        dimuonP4 = muonOneP4 + muonTwoP4;

        if (dimuonP4.M() > 20.)
            return kTRUE;
        hTotalEvents->Fill(7);

        leptonOneP4     = muonOneP4;
        //leptonOneIso    = muons[0];
        leptonOneFlavor = muons[0]->q*13;
        leptonOneD0     = muons[0]->d0;
        leptonOneDZ     = muons[0]->dz;
        
        leptonTwoP4      = muonTwoP4;
        //leptonTwoIso     = muons_iso[1];
        leptonTwoFlavor  = muons[1]->q*13;
        leptonTwoD0      = muons[1]->d0;
        leptonTwoDZ      = muons[1]->dz;

        // Check for dimuon vertex
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
            pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
            pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
            eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);
        }


    } else if (params->selection == "ee") {

        if (electrons.size() != 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        /*TLorentzVector dielectron;
        dielectron = electrons[0] + electrons[1];
        if (dielectron.M() < 12. || dielectron.M() > 70.)
            return kTRUE;
        hTotalEvents->Fill(6);

        leptonOneP4      = electrons[0];
        leptonOneIso     = electrons_iso[0];
        leptonOneQ       = electrons_q[0];
        leptonOneTrigger = electrons_trigger[0];
        leptonOneFlavor  = 0;

        leptonTwoP4      = electrons[1];
        leptonTwoIso     = electrons_iso[1];
        leptonTwoQ       = electrons_q[1];
        leptonTwoTrigger = electrons_trigger[1];
        leptonTwoFlavor  = 0;*/
    } else if (params->selection == "emu") {

        if (muons.size() != 1 || electrons.size() != 1)
            return kTRUE;
        hTotalEvents->Fill(5);

        /*if (muons[0].Pt() < 25 || fabs(muons[0].Eta()) > 2.1)
            return kTRUE;
        hTotalEvents->Fill(6);

        TLorentzVector dilepton;
        dilepton = muons[0] + electrons[0];
        if (dilepton.M() < 12 || dilepton.M() > 70)
            return kTRUE;
        hTotalEvents->Fill(7);

        leptonOneP4      = muons[0];
        leptonOneIso     = muons_iso[0];
        leptonOneQ       = muons_q[0];
        leptonOneTrigger = muons_trigger[0];
        leptonOneFlavor  = 0;

        leptonTwoP4      = electrons[0];
        leptonTwoIso     = electrons_iso[0];
        leptonTwoQ       = electrons_q[0];
        leptonTwoTrigger = electrons_trigger[0];
        leptonTwoFlavor  = 1;

        if (!isData && false) {
            eventWeight *= weights->GetMuonIDEff(muons[0]);

            // trigger efficiency
            pair<float, float> trigEff = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muons[0]);
            eventWeight *= trigEff.first;
        }
    */
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
