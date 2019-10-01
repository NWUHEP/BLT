#include "MultileptonAnalyzer.h"
#include <map>

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

const float ZMASS = 91.19;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 
// bool sort_by_btag(const baconhep::TJet* lhs, const baconhep::TJet* rhs) 
// {
//     return lhs->bmva > rhs->bmva;
// }


MultileptonAnalyzer::MultileptonAnalyzer() : BLTSelector()
{

}

MultileptonAnalyzer::~MultileptonAnalyzer()
{

}

void MultileptonAnalyzer::Begin(TTree *tree)
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

    vector<std::string> channelNames = {"mumu", "ee", "emu", 
                                        "etau", "mutau", 
                                        "e4j", "mu4j", 
                                        "e4j_fakes", "mu4j_fakes",
                                        "mumu_tau","ee_tau","emu_tau"
                                        //"mumu0", "ee0","etau0", "mutau0"
                                        };
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
        tree->Branch("triggerLeptonStatus",&triggerLeptonStatus); // ziheng, 1- only lepton1 fires trigger, 2- only lepton2 fires trigger, 3- both lepton12 fires trigger
        tree->Branch("rPV", &rPV);

        // weights and their uncertainties
        tree->Branch("eventWeight", &eventWeight);
        tree->Branch("leptonOneRecoWeight", &leptonOneRecoWeight);
        tree->Branch("leptonTwoRecoWeight", &leptonTwoRecoWeight);
        tree->Branch("topPtWeight", &topPtWeight);
        tree->Branch("puWeight", &puWeight);
        tree->Branch("triggerWeight", &triggerWeight);
        tree->Branch("genWeight", &genWeight);

        tree->Branch("leptonOneRecoVar", &leptonOneRecoVar);
        tree->Branch("leptonTwoRecoVar", &leptonTwoRecoVar);
        tree->Branch("topPtVar", &topPtVar);
        tree->Branch("puVar", &puVar);
        tree->Branch("triggerVar", &triggerVar);

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

        if (channel != "mu4j" && channel != "e4j" && channel != "mu4j_fakes" && channel != "e4j_fakes") {
            tree->Branch("leptonTwoP4", &leptonTwoP4);
            tree->Branch("leptonTwoIso", &leptonTwoIso);
            tree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
            tree->Branch("leptonTwoMother", &leptonTwoMother);
            tree->Branch("leptonTwoD0", &leptonTwoD0);
            tree->Branch("leptonTwoDZ", &leptonTwoDZ);
        }

        if (channel == "mutau" || channel == "etau" || channel == "mumu_tau"|| channel == "ee_tau"|| channel == "emu_tau") {

            tree->Branch("tauDecayMode",     &tauDecayMode);
            tree->Branch("tauMVA",           &tauMVA);
            tree->Branch("tauPuppiChHadIso", &tauPuppiChHadIso);
            tree->Branch("tauPuppiGammaIso", &tauPuppiGammaIso);
            tree->Branch("tauPuppiNeuHadIso",&tauPuppiNeuHadIso);
            
            tree->Branch("tauGenFlavor",     &tauGenFlavor);
            tree->Branch("tauGenFlavorHad",  &tauGenFlavorHad);
            
            tree->Branch("tauVetoedJetPT",   &tauVetoedJetPt);
            tree->Branch("tauVetoedJetPTUnc",&tauVetoedJetPtUnc);
            

            if (channel == "mumu_tau"|| channel == "ee_tau"|| channel == "emu_tau"){
                tree->Branch("tauP4", &tauP4);
            }
        }




        // jets
        tree->Branch("jetOneP4", &jetOneP4);
        tree->Branch("jetOneTag", &jetOneTag);
        tree->Branch("jetOneFlavor", &jetOneFlavor);
        tree->Branch("jetTwoP4", &jetTwoP4);
        tree->Branch("jetTwoTag", &jetTwoTag);
        tree->Branch("jetTwoFlavor", &jetTwoFlavor);

        if (channel == "mu4j" || channel == "e4j" || channel == "mu4j_fakes" || channel == "e4j_fakes") {
            tree->Branch("jetThreeP4", &jetThreeP4);
            tree->Branch("jetThreeTag", &jetThreeTag);
            tree->Branch("jetThreeFlavor", &jetThreeFlavor);
            tree->Branch("jetFourP4", &jetFourP4);
            tree->Branch("jetFourTag", &jetFourTag);
            tree->Branch("jetFourFlavor", &jetFourFlavor);
        }

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

        // jet counters for systematics
        tree->Branch("nJetsJESUp", &nJetsJESUp);
        tree->Branch("nJetsJESDown", &nJetsJESDown);
        tree->Branch("nJetsJERUp", &nJetsJERUp);
        tree->Branch("nJetsJERDown", &nJetsJERDown);

        tree->Branch("nBJetsJESUp", &nBJetsJESUp);
        tree->Branch("nBJetsJESDown", &nBJetsJESDown);
        tree->Branch("nBJetsJERUp", &nBJetsJERUp);
        tree->Branch("nBJetsJERDown", &nBJetsJERDown);
        tree->Branch("nBJetsBTagUp", &nBJetsBTagUp);
        tree->Branch("nBJetsBTagDown", &nBJetsBTagDown);
        tree->Branch("nBJetsMistagUp", &nBJetsMistagUp);
        tree->Branch("nBJetsMistagDown", &nBJetsMistagDown);

        // theory uncertainties
        tree->Branch("qcdWeights", &qcdWeights);
        tree->Branch("pdfWeight", &pdfWeight);
        tree->Branch("alphaS", &alphaS);

        outTrees[channel] = tree;

        // event counter
        string outHistName = params->get_output_treename("TotalEvents_" + channel);
        eventCounts[channel] = new TH1D(outHistName.c_str(),"ChannelCounts",10,0.5,10.5);
    }

    ReportPostBegin();
}

Bool_t MultileptonAnalyzer::Process(Long64_t entry)
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
    vector<TGenParticle*> genTaus, genTausLep, genTausHad;
    vector<int> genMotherId;

    unsigned wCount = 0;
    unsigned tauCount = 0;
    unsigned topCount = 0;
    bitset<6> wDecay;
    bitset<4> tauDecay;
    float topSF = 1.;
    topPtWeight = 1.;
    topPtVar = 0.;
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

        // get weights for assessing PDF and QCD scale systematics
        pdfWeight = 0.;
        alphaS    = 1.;
        qcdWeights.clear();
        if (fLHEWeightArr != 0) {
            for (int i = 0; i < fLHEWeightArr->GetEntries(); ++i) {
                float lheWeight = ((TLHEWeight*)fLHEWeightArr->At(i))->weight;
                int id = ((TLHEWeight*)fLHEWeightArr->At(i))->id;
                if (id >= 1001 && id <= 1009) {
                    qcdWeights.push_back(lheWeight);
                } else if (id >= 2001 && id <= 2100) {		
                    pdfWeight += pow(qcdWeights[0] - lheWeight, 2.);
                } else if (id == 2101 || id == 2102) {
                    alphaS = lheWeight;
                }
            }
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
                ++topCount;
            }

            // Tag events based on neutrino flavors
            unsigned flavor = abs(particle->pdgId);
            if (
                    (flavor == 12 || flavor == 14 || flavor == 16)
                    && particle->parent != -2
               ) {
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                //cout << i << ", " << particle->pdgId << ", " << mother->pdgId;
                if (abs(mother->pdgId) == 24) {
                    wDecay.set(3*wCount + (flavor - 12)/2);
                    ++wCount;
                } 

                if (abs(mother->pdgId) == 15 && flavor != 16 && mother->parent != -2) {
                    TGenParticle* gmother = (TGenParticle*) fGenParticleArr->At(mother->parent);
                    if (gmother->parent == -2 || tauCount == 2) continue;
                    tauDecay.set(2*tauCount + (flavor - 12)/2);
                    ++tauCount;
                } 
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


            // This saves gen level tau and leptonic taus
            if ( particle->parent != -2 ) {
                if ( abs(particle->pdgId) == 16) {
                    TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                    if (abs(mother->pdgId) == 15) { 
                        genTaus.push_back(mother); 
                    }
                }

                if ( abs(particle->pdgId) == 11 || abs(particle->pdgId) == 13 ) {
                    TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                    if (abs(mother->pdgId) == 15) {
                        genTausLep.push_back(mother);
                    }
                }
                
            }
            
        }


        // genTausHad = genTaus - genTausLep
        for (unsigned i = 0; i < genTaus.size(); ++i) {
            
            bool isHad = true;

            TLorentzVector genTauP4;
            genTauP4.SetPtEtaPhiM(genTaus[i]->pt, genTaus[i]->eta, genTaus[i]->phi, genTaus[i]->mass); 

            for (unsigned j = 0; j < genTausLep.size(); ++j) {
                TLorentzVector genTauLepP4;
                genTauLepP4.SetPtEtaPhiM(genTausLep[j]->pt, genTausLep[j]->eta, genTausLep[j]->phi, genTausLep[j]->mass); 

                if (genTauP4.DeltaR(genTauLepP4) < 0.2) {
                    isHad = false;
                }
            }
            if (isHad) genTausHad.push_back(genTaus[i]);
        }


        // if(genTaus.size() >2){

        //     cout <<  genTaus.size() <<  "," << genTausLep.size()  <<  "," << genTausHad.size() << endl;
        //     for (unsigned i = 0; i < genTaus.size(); ++i){
        //         cout << "tau " << i <<" parent index is "<<genTaus[i]->parent <<endl;
        //     }
        //     cout << endl;
        // }


        nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY

        // categorize events
        genCategory = 0;
        if (wCount == 2) { // fully leptonic
            if (wDecay.test(0) && wDecay.test(3)) {// W->e, W->e
                hGenCat->Fill(1); 
                genCategory = 1;
            } else if (wDecay.test(1) && wDecay.test(4)) {// W->mu, W->mu
                hGenCat->Fill(2); 
                genCategory = 2;
            } else if ((wDecay.test(0) && wDecay.test(4)) || (wDecay.test(1) && wDecay.test(3))) {// W->e, W->mu
                hGenCat->Fill(3); 
                genCategory = 3;
            } else if (wDecay.test(2) && wDecay.test(5)) {// W->tau, W->tau
                if (tauDecay.test(0) && tauDecay.test(2)) {// tau->e, tau->e
                    hGenCat->Fill(4); 
                    genCategory = 4;
                } else if (tauDecay.test(1) && tauDecay.test(3)) {// tau->mu, tau->mu
                    hGenCat->Fill(5); 
                    genCategory = 5;
                } else if ((tauDecay.test(1) && tauDecay.test(2)) || (tauDecay.test(0) && tauDecay.test(3))) {// tau->e, tau->mu
                    hGenCat->Fill(6); 
                    genCategory = 6;
                } else if (tauDecay.test(0) && tauCount == 1) {// tau->e, tau->h
                    hGenCat->Fill(7); 
                    genCategory = 7;
                } else if (tauDecay.test(1) && tauCount == 1) {// tau->mu, tau->h
                    hGenCat->Fill(8); 
                    genCategory = 8;
                } else if (tauCount == 0) {// tau->h, tau->h
                    hGenCat->Fill(9); 
                    genCategory = 9;
                }
            } else if ((wDecay.test(0) && wDecay.test(5)) || (wDecay.test(2) && wDecay.test(3))) {// W->e, W->tau
                if (tauDecay.test(0)) {// tau->e
                    hGenCat->Fill(10);
                    genCategory = 10;
                } else if (tauDecay.test(1)) {// tau->mu
                    hGenCat->Fill(11);
                    genCategory = 11;
                } else if (tauDecay.none()) {// tau->h
                    hGenCat->Fill(12);
                    genCategory = 12;
                }
            } else if ((wDecay.test(1) && wDecay.test(5)) || (wDecay.test(2) && wDecay.test(4))) {// W->mu, W->tau
                if (tauDecay.test(0)) {// tau->e
                    hGenCat->Fill(13);
                    genCategory = 13;
                } else if (tauDecay.test(1)) {// tau->mu
                    hGenCat->Fill(14);
                    genCategory = 14;
                } else if (tauDecay.none()) {// tau->h
                    hGenCat->Fill(15);
                    genCategory = 15;
                }
            } 
        } else if (wCount == 1) { // semileptonic
            if (wDecay.test(0)) {// W->e, W->h
                hGenCat->Fill(16); 
                genCategory = 16;
            } else if (wDecay.test(1)) {// W->mu, W->h
                hGenCat->Fill(17); 
                genCategory = 17;
            } else if (wDecay.test(2)) {// W->tau, W->h
                if (tauDecay.test(0)) {// tau->e
                    hGenCat->Fill(18);
                    genCategory = 18;
                } else if (tauDecay.test(1)) {// tau->mu
                    hGenCat->Fill(19);
                    genCategory = 19;
                } else if (tauDecay.none()) {// tau->h
                    hGenCat->Fill(20);
                    genCategory = 20;
                }
            }
        } else { // hadronic
            hGenCat->Fill(21); 
            genCategory = 21;
        }
        // Account for the top pt weights
        if (params->datasetgroup.substr(0, 5) == "ttbar") {
            topPtWeight *= sqrt(topSF);
            topPtVar    += pow(0.01*topPtWeight, 2);
        }

        // pileup reweighting
        nPU          = fInfo->nPUmean;
        puWeight     = weights->GetPUWeight(fInfo->nPUmean); 
        puVar        = 0.01*puWeight;
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

    

    /* -------- Vertices ---------*/
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



    
    /* -------- MUONS ---------*/
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
            if (GetMuonIsolation(muon)/muonP4.Pt() <  0.15){
                muons.push_back(muon);
                veto_muons.push_back(muonP4);
            } else {
                fail_muons.push_back(muon);
            }
        }
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);
    sort(fail_muons.begin(), fail_muons.end(), sort_by_higher_pt<TMuon>);




    
    /* -------- ELECTRONS ---------*/
    vector<TElectron*> electrons,fail_electrons;
    vector<TLorentzVector> veto_electrons;

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
            if (particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl) ){
                electrons.push_back(electron);
                veto_electrons.push_back(electronP4);
            } else {
                fail_electrons.push_back(electron);
            }
        }
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);
    sort(fail_electrons.begin(), fail_electrons.end(), sort_by_higher_pt<TElectron>);



    
    /* -------- TAUS ---------*/
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

            bool genTauHadOverlap = false;
            // check if can be tagged as hadronic tau
            for (unsigned i = 0; i < genTausHad.size(); ++i) {
                TLorentzVector genP4;
                genP4.SetPtEtaPhiM(genTausHad[i]->pt, genTausHad[i]->eta, genTausHad[i]->phi, genTausHad[i]->mass); 
                if (tauP4.DeltaR(genP4) < 0.3) {
                    genTauHadOverlap = true;
                    break;
                }
            }

            // apply the tau energy scale correction only for true tau ID
            if (genTauHadOverlap) {
                if (tau->decaymode == 0) {
                    tau->pt *= 0.995;
                } else if (tau->decaymode == 1) {
                    tau->pt *= 1.01;
                } else if (tau->decaymode == 10) {
                    tau->pt *= 1.006;
                }
            }

        }



        if( 
                tau->pt > 20  
                && abs(tau->eta) < 2.3 
                && !muOverlap
                && !elOverlap
                && (tau->hpsDisc & baconhep::kByDecayModeFinding)
                
                && (tau->hpsDisc & baconhep::kByVTightIsolationMVA3newDMwLT)
                && (tau->hpsDisc & baconhep::kByMVA6VTightElectronRejection)
                && (tau->hpsDisc & baconhep::kByTightMuonRejection3)
          ) {
            taus.push_back(tau);
            veto_taus.push_back(tauP4);
        }
    }
    sort(taus.begin(), taus.end(), sort_by_higher_pt<TTau>);
    

    
    /* -------- JETS ---------*/
    TClonesArray* jetCollection;
    
    jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets, vetoedJets;
    TLorentzVector hadronicP4;
    float sumJetPt = 0;

    ResetJetCounters();
    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

        // apply JEC offline and get scale uncertainties
        double jec = particleSelector->JetCorrector(jet, "NONE");
        jet->pt = jet->ptRaw*jec;


        float gRand = 1.;
        float jerc = 1.;
        if (!isData) { // apply jet energy resolution corrections to simulation
            pair<float, float> resPair = particleSelector->JetResolutionAndSF(jet, 0);
            gRand = rng->Gaus(0, resPair.first);
            jerc = 1 + gRand*sqrt(std::max((double)resPair.second*resPair.second - 1, 0.));
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

        bool tauOverlap = false;
        for (const auto& tau: veto_taus) {
            if (jetP4.DeltaR(tau) < 0.4) {
                tauOverlap = true;
                vetoedJets.push_back(jet);
                break;
            }
        }

        if (
                fabs(jet->eta) < 4.7
                && particleSelector->PassJetID(jet, cuts->looseJetID)
                && !muOverlap 
                && !elOverlap
                && !tauOverlap
           ) {

            if (fabs(jet->eta) <= 2.4) { 

                if (isData) {
                    if (jet->pt > 30) {
                        ++nJets;
                        if (jet->bmva > 0.9432) { 
                            ++nBJets;
                        } 

                        // if (jet->csv > 0.8484) { 
                        //      ++nBJets;
                        // } 
                    }
                } else {
                    JetCounting(jet,jerc, gRand);
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
    std::sort(jets.begin(), jets.end(), sort_by_btag<TJet>);

    //cout << nBJets << ", " << nBJetsMistagUp << ", " << nBJetsMistagDown << endl;

    // use the highest jet multiplicities given all systematic variations
    unsigned nJetList[] = {nJets, nJetsJESUp, nJetsJESDown, nJetsJERUp, nJetsJERDown};
    nJetsCut  = *std::max_element(nJetList, nJetList+sizeof(nJetList)/sizeof(unsigned));
    nJetsCut2 = *std::min_element(nJetList, nJetList+sizeof(nJetList)/sizeof(unsigned));

    unsigned nBJetList[] = {nBJets, nBJetsJESUp, nBJetsJESDown, 
                            nBJetsJERUp, nBJetsJERDown, 
                            nBJetsBTagUp, nBJetsBTagDown, 
                            nBJetsMistagUp, nBJetsMistagDown
                           };
    nBJetsCut = *std::max_element(nBJetList, nBJetList+sizeof(nBJetList)/sizeof(unsigned));

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
    nTaus      = taus.size();

    // trigger selections
    bool muonTriggered = find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_IsoMu24_v*") != passTriggerNames.end() 
        || find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_IsoTkMu24_v*") != passTriggerNames.end();
    bool electronTriggered = find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_Ele27_WPTight_Gsf_v*") != passTriggerNames.end();

    string channel = "";
    if        (electrons.size() == 0 && muons.size() == 2 && taus.size() == 0 ){ // mu+mu selection
        
    
        // if (nJetsCut >= 2 && nBJetsCut >= 1 ){ //
        //     channel = "mumu";
        // } else if (nJetsCut2 ==0 ){
        //     //channel = "mumu0";
        //     return kTRUE;
        // } else {
        //     return kTRUE;
        // }

        channel = "mumu";
        eventCounts[channel]->Fill(1);

        // convert to TLorentzVectors
        TLorentzVector muonOneP4, muonTwoP4, dimuonP4;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, 0.1052);
        dimuonP4 = muonOneP4 + muonTwoP4;

        float muonOneIso = GetMuonIsolation(muons[0]);
        float muonTwoIso = GetMuonIsolation(muons[1]);
        if (    !muonTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        if (muons[0]->pt < 25. || muons[1]->pt < 10.)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        if (nJetsCut < 2 || nBJetsCut < 1)
            return kTRUE;
        eventCounts[channel]->Fill(4);

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

        // trigger weights with trigger matching
        bitset<2> triggered;
        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, muons[0]->hltMatchBits) && muonOneP4.Pt() > 25)
                triggered.set(0);
            if (trigger->passObj(name, 1, muons[1]->hltMatchBits) && muonTwoP4.Pt() > 25)
                triggered.set(1);
        }

        if (!triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 0;
        if ( triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 1;
        if (!triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 2;
        if ( triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 3;


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
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            effs = effCont2.GetEff();
            errs = effCont2.GetErr();
            leptonTwoRecoWeight = effs.first/effs.second;
            leptonTwoRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            

            // trigger weights with trigger matching:
            //
            // check if lepton could pass the trigger threshold and is matched
            // to a trigger object.  When both muons pass the trigger, use the
            // efficiency for detecting either

            if (triggered.all()) {
                effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
                effCont2      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
                triggerWeight = GetTriggerSF(effCont1, effCont2);
                triggerVar    = GetTriggerSFError(effCont1, effCont2);
            } else if (triggered.test(0)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else if (triggered.test(1)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }

            // update the event weight
            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;
            eventWeight *= triggerWeight;
            
        }
    } else if (electrons.size() == 2 && muons.size() == 0 && taus.size() == 0) { // e+e selection

        // if (nJetsCut >= 2 ){// && nBJetsCut >= 1){
        //     channel = "ee";
        // } else if (nJetsCut2 == 0 ){
        //     //channel = "ee0";
        //     return kTRUE;
        // } else {
        //     return kTRUE;
        // }

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

        if (nJetsCut < 2 || nBJetsCut < 1)
            return kTRUE;
        eventCounts[channel]->Fill(4);

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

        // trigger weights with trigger matching
        bitset<2> triggered;
        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, electrons[0]->hltMatchBits) && electronOneP4.Pt() > 30)
                triggered.set(0);
            if (trigger->passObj(name, 1, electrons[1]->hltMatchBits) && electronTwoP4.Pt() > 30)
                triggered.set(1);
        }

        if (!triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 0;
        if ( triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 1;
        if (!triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 2;
        if ( triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 3;

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
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            effs = effCont2.GetEff();
            errs = effCont2.GetErr();
            leptonTwoRecoWeight = effs.first/effs.second;
            leptonTwoRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));


            if (triggered.all()) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronOneP4);
                effCont2      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronTwoP4);
                triggerWeight = GetTriggerSF(effCont1, effCont2);
                triggerVar    = GetTriggerSFError(effCont1, effCont2);
            } else if (triggered.test(0)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronOneP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else if (triggered.test(1)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronTwoP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }

            // update event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;
        }
    } else if (electrons.size() == 1 && muons.size() == 1 && taus.size() == 0) { // e+mu selection
        channel = "emu";
        eventCounts[channel]->Fill(1);

        // trigger matching for thresholds
        
        float muPtThreshold = 10.;
        float elPtThreshold = 10.;
        if (
                muonTriggered
           ) {
            muPtThreshold = 25.;
        } else if (electronTriggered) {
            elPtThreshold = 30.;
        }

        if (
                muons[0]->pt < muPtThreshold 
                || electrons[0]->pt < elPtThreshold
                || !(muonTriggered || electronTriggered)
           )
            return kTRUE;
        eventCounts[channel]->Fill(2);

        TLorentzVector muonP4, electronP4, dilepton;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        electronP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        dilepton = muonP4 + electronP4;
        if (dilepton.M() < 12)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        if (nJetsCut < 2 || nBJetsCut < 1)
            return kTRUE;
        eventCounts[channel]->Fill(4);

        leptonOneP4     = muonP4;
        leptonOneIso    = GetMuonIsolation(muons[0]);
        leptonOneFlavor = 13*muons[0]->q;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;

        leptonTwoP4     = electronP4;
        leptonTwoIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonTwoFlavor = 11*electrons[0]->q;
        leptonTwoDZ     = electrons[0]->dz;
        leptonTwoD0     = electrons[0]->d0;

        // trigger weights with trigger matching
        bitset<2> triggered;
        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, muons[0]->hltMatchBits))
                triggered.set(0);
            if (trigger->passObj(name, 1, electrons[0]->hltMatchBits))
                triggered.set(1);
        }

        if (!triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 0;
        if ( triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 1;
        if (!triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 2;
        if ( triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 3;



        if (!isData) {
            leptonOneMother = GetGenMotherId(genParticles, muonP4);
            leptonTwoMother = GetGenMotherId(genParticles, electronP4);

            // reconstruction weights
            EfficiencyContainer effCont1, effCont2;
            effCont1 = weights->GetMuonRecoEff(muonP4);
            effCont2 = weights->GetElectronRecoEff(electronP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            effs = effCont2.GetEff();
            errs = effCont2.GetErr();
            leptonTwoRecoWeight = effs.first/effs.second;
            leptonTwoRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));


            if (triggered.all()) {
                effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
                effCont2      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronP4);
                triggerWeight = GetTriggerSF(effCont1, effCont2);
                triggerVar    = GetTriggerSFError(effCont1, effCont2);
            } else if (triggered.test(0)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else if (triggered.test(1)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }

            // update event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;
        }
    } else if (electrons.size() == 1 && muons.size() == 0 && taus.size() == 1) { // e+tau selection
        

        // if (nJetsCut >= 2 && nBJetsCut >= 1){
        //     channel = "etau";
        // } else if (nJetsCut2 == 0 ){
        //     //channel = "etau0";
        //     return kTRUE;
        // } else {
        //     return kTRUE;
        // }

        channel = "etau";
        eventCounts[channel]->Fill(1);

        if (electrons[0]->pt < 30 || !electronTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        TLorentzVector electronP4, tauP4, dilepton;
        electronP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        tauP4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, 1.776);
        dilepton = electronP4 + tauP4;
        if (dilepton.M() < 12)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        if (nJetsCut < 2 || nBJetsCut < 1)
            return kTRUE;
        eventCounts[channel]->Fill(4);

        leptonOneP4     = electronP4;
        leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonOneFlavor = 11*electrons[0]->q;
        leptonOneDZ     = electrons[0]->dz;
        leptonOneD0     = electrons[0]->d0;

        leptonTwoP4     = tauP4;
        leptonTwoIso    = 0.;
        leptonTwoFlavor = 15*taus[0]->q;
        leptonTwoDZ     = taus[0]->dzLeadChHad;
        leptonTwoD0     = taus[0]->d0LeadChHad;

        tauDecayMode      = taus[0]->decaymode;
        tauMVA            = taus[0]->rawIsoMVA3newDMwLT;
        tauPuppiChHadIso  = taus[0]->puppiChHadIso;
        tauPuppiGammaIso  = taus[0]->puppiGammaIso;
        tauPuppiNeuHadIso = taus[0]->puppiNeuHadIso;

        pair <float, float> tauVetoedJetPtPair;
        tauVetoedJetPtPair= GetTauVetoedJetPt(tauP4, vetoedJets);
        tauVetoedJetPt    = tauVetoedJetPtPair.first;
        tauVetoedJetPtUnc = tauVetoedJetPtPair.second;
        // cout<<tauVetoedJetPt<<"+/-"<<tauVetoedJetPtUnc<<endl;



        // trigger weights
        bool triggered = false;
        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, electrons[0]->hltMatchBits)) {
                triggered = true;
                break;
            }
        }

        if (!triggered) triggerLeptonStatus = 0;
        if ( triggered) triggerLeptonStatus = 1;


        if (!isData) {
            leptonOneMother = GetGenMotherId(genParticles, electronP4);
            tauGenFlavor    = GetTauGenFlavor(tauP4,genTausHad,vetoedJets, false); // useHadronFlavor = false
            tauGenFlavorHad = GetTauGenFlavor(tauP4,genTausHad,vetoedJets, true);  // useHadronFlavor = true

            // reconstruction weights
            EfficiencyContainer effCont1;
            effCont1 = weights->GetElectronRecoEff(electronP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            leptonTwoRecoWeight = 0.95;
            leptonTwoRecoVar    = 0.05;


            if (triggered) {
                effCont1 = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronP4);
                effs = effCont1.GetEff();
                errs = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(triggerWeight, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }

            // update event weight
            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;
            eventWeight *= triggerWeight;

        }
    } else if (electrons.size() == 0 && muons.size() == 1 && taus.size() == 1) { // mu+tau selection
        
        // if (nJetsCut >= 2 && nBJetsCut >= 1){
        //     channel = "mutau";
        // } else if (nJetsCut2 == 0 ){
        //     //channel = "mutau0";
        //     return kTRUE;
        // } else {
        //     return kTRUE;
        // }

        channel = "mutau";
        eventCounts[channel]->Fill(1);

        if (muons[0]->pt < 25 || !muonTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        TLorentzVector muonP4, tauP4, dilepton;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        tauP4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, 1.776);
        dilepton = muonP4 + tauP4;
        if (dilepton.M() < 12)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        if (nJetsCut < 2 || nBJetsCut < 1)
            return kTRUE;
        eventCounts[channel]->Fill(4);

        leptonOneP4     = muonP4;
        leptonOneIso    = GetMuonIsolation(muons[0]);
        leptonOneFlavor = 13*muons[0]->q;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;

        leptonTwoP4     = tauP4;
        leptonTwoIso    = 0.;
        leptonTwoFlavor = 15*taus[0]->q;
        leptonTwoDZ     = taus[0]->dzLeadChHad;
        leptonTwoD0     = taus[0]->d0LeadChHad;


        tauDecayMode      = taus[0]->decaymode;
        tauMVA            = taus[0]->rawIsoMVA3newDMwLT;
        tauPuppiChHadIso  = taus[0]->puppiChHadIso;
        tauPuppiGammaIso  = taus[0]->puppiGammaIso;
        tauPuppiNeuHadIso = taus[0]->puppiNeuHadIso;

        pair <float, float> tauVetoedJetPtPair;
        tauVetoedJetPtPair= GetTauVetoedJetPt(tauP4, vetoedJets);
        tauVetoedJetPt    = tauVetoedJetPtPair.first;
        tauVetoedJetPtUnc = tauVetoedJetPtPair.second;

        // cout<<tauVetoedJetPt<<"+/-"<<tauVetoedJetPtUnc<<endl;



        // trigger weights
        bool triggered = false;
        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, muons[0]->hltMatchBits)) {
                triggered = true;
                break;
            }
        }

        if (!triggered) triggerLeptonStatus = 0;
        if ( triggered) triggerLeptonStatus = 1;


        if (!isData) {
            leptonOneMother = GetGenMotherId(genParticles, muonP4);
            tauGenFlavor    = GetTauGenFlavor(tauP4,genTausHad,vetoedJets, false); // useHadronFlavor = false
            tauGenFlavorHad = GetTauGenFlavor(tauP4,genTausHad,vetoedJets, true);  // useHadronFlavor = true

            // reconstruction weights
            EfficiencyContainer effCont1, effCont2;
            effCont1 = weights->GetMuonRecoEff(muonP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            leptonTwoRecoWeight = 0.95;
            leptonTwoRecoVar    = 0.05;


            if (triggered) {
                effCont1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
                effs = effCont1.GetEff();
                errs = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(triggerWeight, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }


            // update event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;

        }
    } else if (electrons.size() == 1 && muons.size() == 0 && taus.size() == 0) { // e+h selection
        channel = "e4j";
        eventCounts[channel]->Fill(1);

        // convert to TLorentzVectors
        TLorentzVector electronP4;
        electronP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        float electronIso = GetElectronIsolation(electrons[0], fInfo->rhoJet);

        if (electrons[0]->pt < 30. || !electronTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        if (nJetsCut < 4 || nBJetsCut < 1)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        leptonOneP4     = electronP4;
        leptonOneIso    = electronIso;
        leptonOneFlavor = electrons[0]->q*11;
        leptonOneDZ     = electrons[0]->dz;
        leptonOneD0     = electrons[0]->d0;

        // Collect the highest pt jets in the event
        std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);
        //jets = KinematicTopTag(jets, metP2, electronP4);
        if (nJets >= 4) {
            jetOneP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
            jetOneTag      = jets[0]->csv;
            jetOneFlavor   = jets[0]->hadronFlavor;
            jetTwoP4.SetPtEtaPhiM(jets[1]->pt, jets[1]->eta, jets[1]->phi, jets[1]->mass);
            jetTwoTag      = jets[1]->csv;
            jetTwoFlavor   = jets[1]->hadronFlavor;
            jetThreeP4.SetPtEtaPhiM(jets[2]->pt, jets[2]->eta, jets[2]->phi, jets[2]->mass);
            jetThreeTag    = jets[2]->csv;
            jetThreeFlavor = jets[2]->hadronFlavor;
            jetFourP4.SetPtEtaPhiM(jets[3]->pt, jets[3]->eta, jets[3]->phi, jets[3]->mass);
            jetFourTag     = jets[3]->csv;
            jetFourFlavor  = jets[3]->hadronFlavor;
        }


        // trigger weights
        bool triggered = false;
        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, electrons[0]->hltMatchBits)) {
                triggered = true;
                break;
            }
        }
        if (!triggered) triggerLeptonStatus = 0;
        if ( triggered) triggerLeptonStatus = 1;


        if (!isData) {
            leptonOneMother = GetGenMotherId(genParticles, electronP4);

            // reconstruction weights
            EfficiencyContainer effCont1;
            effCont1 = weights->GetElectronRecoEff(electronP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            leptonTwoRecoWeight = 0.;



            if (triggered) {
                effCont1 = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronP4);
                effs = effCont1.GetEff();
                errs = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(triggerWeight, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }



            // update event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneRecoWeight;
        }
    } else if (electrons.size() == 0 && muons.size() == 1 && taus.size() == 0) { // mu+h selection
        channel = "mu4j";
        eventCounts[channel]->Fill(1);

        // convert to TLorentzVectors
        TLorentzVector muonP4;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        float muonIso = GetMuonIsolation(muons[0]);
        if (muons[0]->pt < 25. || !muonTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(2);


        if (nJetsCut < 4 || nBJetsCut < 1)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        leptonOneP4     = muonP4;
        leptonOneIso    = muonIso;
        leptonOneFlavor = muons[0]->q*13;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;

        // Collect the highest pt jets in the event
        std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);
        //jets = KinematicTopTag(jets, metP2, muonP4);
        if (nJets >= 4) {
            jetOneP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
            jetOneTag      = jets[0]->csv;
            jetOneFlavor   = jets[0]->hadronFlavor;
            jetTwoP4.SetPtEtaPhiM(jets[1]->pt, jets[1]->eta, jets[1]->phi, jets[1]->mass);
            jetTwoTag      = jets[1]->csv;
            jetTwoFlavor   = jets[1]->hadronFlavor;
            jetThreeP4.SetPtEtaPhiM(jets[2]->pt, jets[2]->eta, jets[2]->phi, jets[2]->mass);
            jetThreeTag    = jets[2]->csv;
            jetThreeFlavor = jets[2]->hadronFlavor;
            jetFourP4.SetPtEtaPhiM(jets[3]->pt, jets[3]->eta, jets[3]->phi, jets[3]->mass);
            jetFourTag     = jets[3]->csv;
            jetFourFlavor  = jets[3]->hadronFlavor;
        }

        // trigger weights
        bool triggered = false;
        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, muons[0]->hltMatchBits)) {
                triggered = true;
                break;
            }
        }

        if (!triggered) triggerLeptonStatus = 0;
        if ( triggered) triggerLeptonStatus = 1;
        
        if (!isData) {
            leptonOneMother = GetGenMotherId(genParticles, muonP4);

            // reconstruction weights
            EfficiencyContainer effCont1, effCont2;
            effCont1 = weights->GetMuonRecoEff(muonP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            leptonTwoRecoWeight = 0.;



            if (triggered) {
                effCont1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
                effs = effCont1.GetEff();
                errs = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(triggerWeight, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }

            // update event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneRecoWeight;
        }
    } else if (electrons.size() == 0 && muons.size() == 0 && taus.size() == 0 && fail_muons.size() == 1) { // mu_fake+h selection
        channel = "mu4j_fakes";
        eventCounts[channel]->Fill(1);
        
        // convert to TLorentzVectors
        TLorentzVector muonP4;
        muonP4.SetPtEtaPhiM(fail_muons[0]->pt, fail_muons[0]->eta, fail_muons[0]->phi, 0.1052);
        float muonIso = GetMuonIsolation(fail_muons[0]);
        if (fail_muons[0]->pt < 25. || !muonTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        // need 4 jets with one btag
        if (nJetsCut < 4 || nBJetsCut < 1)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        leptonOneP4     = muonP4;
        leptonOneIso    = muonIso;
        leptonOneFlavor = fail_muons[0]->q*13;
        leptonOneDZ     = fail_muons[0]->dz;
        leptonOneD0     = fail_muons[0]->d0;

        // Collect the highest pt jets in the event
        std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);
        //jets = KinematicTopTag(jets, metP2, muonP4);
        if (nJets >= 4) {
            jetOneP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
            jetOneTag      = jets[0]->csv;
            jetOneFlavor   = jets[0]->hadronFlavor;
            jetTwoP4.SetPtEtaPhiM(jets[1]->pt, jets[1]->eta, jets[1]->phi, jets[1]->mass);
            jetTwoTag      = jets[1]->csv;
            jetTwoFlavor   = jets[1]->hadronFlavor;
            jetThreeP4.SetPtEtaPhiM(jets[2]->pt, jets[2]->eta, jets[2]->phi, jets[2]->mass);
            jetThreeTag    = jets[2]->csv;
            jetThreeFlavor = jets[2]->hadronFlavor;
            jetFourP4.SetPtEtaPhiM(jets[3]->pt, jets[3]->eta, jets[3]->phi, jets[3]->mass);
            jetFourTag     = jets[3]->csv;
            jetFourFlavor  = jets[3]->hadronFlavor;
        }

        // trigger weights
        bool triggered = false;
        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, fail_muons[0]->hltMatchBits)) {
                triggered = true;
                break;
            }
        }

        if (!triggered) triggerLeptonStatus = 0;
        if ( triggered) triggerLeptonStatus = 1;
        
        if (!isData) {
            // leptonOneMother = GetGenMotherId(genParticles, muonP4);

            // reconstruction weights
            EfficiencyContainer effCont1, effCont2;
            effCont1 = weights->GetMuonRecoEff(muonP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            leptonTwoRecoWeight = 0.;



            if (triggered) {
                effCont1 = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
                effs = effCont1.GetEff();
                errs = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(triggerWeight, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }

            // update event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneRecoWeight;
        }

    } else if (electrons.size() == 0 && muons.size() == 0 && taus.size() == 0 && fail_electrons.size() == 1) { // e_fake+h selection
        channel = "e4j_fakes";
        eventCounts[channel]->Fill(1);

        

        // convert to TLorentzVectors
        TLorentzVector electronP4;
        electronP4.SetPtEtaPhiM(fail_electrons[0]->pt, fail_electrons[0]->eta, fail_electrons[0]->phi, 511e-6);
        float electronIso = GetElectronIsolation(fail_electrons[0], fInfo->rhoJet);
        

        if (fail_electrons[0]->pt < 30. || !electronTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(2);
        

        // need 4 jets with one btag
        if (nJetsCut < 4 || nBJetsCut < 1)
            return kTRUE;
        eventCounts[channel]->Fill(3);
        
        leptonOneP4     = electronP4;
        leptonOneIso    = electronIso;
        leptonOneFlavor = fail_electrons[0]->q*11;
        leptonOneDZ     = fail_electrons[0]->dz;
        leptonOneD0     = fail_electrons[0]->d0;
        
        
        // Collect the highest pt jets in the event
        std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);
        //jets = KinematicTopTag(jets, metP2, electronP4);
        if (nJets >= 4) {
            jetOneP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
            jetOneTag      = jets[0]->csv;
            jetOneFlavor   = jets[0]->hadronFlavor;
            jetTwoP4.SetPtEtaPhiM(jets[1]->pt, jets[1]->eta, jets[1]->phi, jets[1]->mass);
            jetTwoTag      = jets[1]->csv;
            jetTwoFlavor   = jets[1]->hadronFlavor;
            jetThreeP4.SetPtEtaPhiM(jets[2]->pt, jets[2]->eta, jets[2]->phi, jets[2]->mass);
            jetThreeTag    = jets[2]->csv;
            jetThreeFlavor = jets[2]->hadronFlavor;
            jetFourP4.SetPtEtaPhiM(jets[3]->pt, jets[3]->eta, jets[3]->phi, jets[3]->mass);
            jetFourTag     = jets[3]->csv;
            jetFourFlavor  = jets[3]->hadronFlavor;
        }

        
        // trigger weights
        bool triggered = false;
        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, fail_electrons[0]->hltMatchBits)) {
                triggered = true;
                break;
            }
        }
        if (!triggered) triggerLeptonStatus = 0;
        if ( triggered) triggerLeptonStatus = 1;

        
        if (!isData) {
            // leptonOneMother = GetGenMotherId(genParticles, electronP4);

            // reconstruction weights
            EfficiencyContainer effCont1;
            effCont1 = weights->GetElectronRecoEff(electronP4);

            pair<float, float> effs, errs;
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            leptonTwoRecoWeight = 0.;

            if (triggered) {
                effCont1 = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronP4);
                effs = effCont1.GetEff();
                errs = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(triggerWeight, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }

            // update event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneRecoWeight;
        }
    // } else if (electrons.size() == 0 && muons.size() == 2 && taus.size() == 1) { // mu+mu+tau selection
        
    
    //     // if (nJetsCut >= 2 && nBJetsCut >= 1 ){ //
    //     //     channel = "mumu";
    //     // } else if (nJetsCut2 ==0 ){
    //     //     //channel = "mumu0";
    //     //     return kTRUE;
    //     // } else {
    //     //     return kTRUE;
    //     // }

    //     channel = "mumu_tau";
    //     eventCounts[channel]->Fill(1);

    //     // convert to TLorentzVectors
    //     TLorentzVector muonOneP4, muonTwoP4, dimuonP4;
    //     muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
    //     muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, 0.1052);
    //     dimuonP4 = muonOneP4 + muonTwoP4;

    //     float muonOneIso = GetMuonIsolation(muons[0]);
    //     float muonTwoIso = GetMuonIsolation(muons[1]);
    //     if (    !muonTriggered)
    //         return kTRUE;
    //     eventCounts[channel]->Fill(2);

    //     if (muons[0]->pt < 25. || muons[1]->pt < 10.)
    //         return kTRUE;
    //     eventCounts[channel]->Fill(3);

    //     // if (nJetsCut < 2 || nBJetsCut < 1)
    //     //     return kTRUE;
    //     // eventCounts[channel]->Fill(4);

    //     leptonOneP4     = muonOneP4;
    //     leptonOneIso    = muonOneIso;
    //     leptonOneFlavor = muons[0]->q*13;
    //     leptonOneDZ     = muons[0]->dz;
    //     leptonOneD0     = muons[0]->d0;

    //     leptonTwoP4     = muonTwoP4;
    //     leptonTwoIso    = muonTwoIso;
    //     leptonTwoFlavor = muons[1]->q*13;
    //     leptonTwoDZ     = muons[1]->dz;
    //     leptonTwoD0     = muons[1]->d0;

    //     // trigger weights with trigger matching
    //     bitset<2> triggered;
    //     for (const auto& name: passTriggerNames) {
    //         if (trigger->passObj(name, 1, muons[0]->hltMatchBits) && muonOneP4.Pt() > 25)
    //             triggered.set(0);
    //         if (trigger->passObj(name, 1, muons[1]->hltMatchBits) && muonTwoP4.Pt() > 25)
    //             triggered.set(1);
    //     }

    //     if (!triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 0;
    //     if ( triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 1;
    //     if (!triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 2;
    //     if ( triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 3;

    //     ////////////////////////////////////////
    //     // tau reconstruct info
    //     TLorentzVector tau0P4;
    //     tau0P4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, 1.77682);

    //     tauP4             = tau0P4;

    //     tauDecayMode      = taus[0]->decaymode;
    //     tauMVA            = taus[0]->rawIsoMVA3newDMwLT;
    //     tauPuppiChHadIso  = taus[0]->puppiChHadIso;
    //     tauPuppiGammaIso  = taus[0]->puppiGammaIso;
    //     tauPuppiNeuHadIso = taus[0]->puppiNeuHadIso;
        
    //     pair <float, float> tauVetoedJetPtPair;
    //     tauVetoedJetPtPair= GetTauVetoedJetPt(tauP4, vetoedJets);
    //     tauVetoedJetPt    = tauVetoedJetPtPair.first;
    //     tauVetoedJetPtUnc = tauVetoedJetPtPair.second;

    //     ////////////////////////////////////////

    //     if (!isData) {
    //         leptonOneMother = GetGenMotherId(genParticles, muonOneP4);
    //         leptonTwoMother = GetGenMotherId(genParticles, muonTwoP4);
    //         tauGenFlavor    = GetTauGenFlavor(tau0P4, genTausHad, vetoedJets, false); // useHadronFlavor = false
    //         tauGenFlavorHad = GetTauGenFlavor(tau0P4, genTausHad, vetoedJets, true);  // useHadronFlavor = true

    //         // reconstruction weights
    //         EfficiencyContainer effCont1, effCont2;
    //         effCont1 = weights->GetMuonRecoEff(muonOneP4);
    //         effCont2 = weights->GetMuonRecoEff(muonTwoP4);

    //         pair<float, float> effs, errs;
    //         effs = effCont1.GetEff();
    //         errs = effCont1.GetErr();
    //         leptonOneRecoWeight = effs.first/effs.second;
    //         leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

    //         effs = effCont2.GetEff();
    //         errs = effCont2.GetErr();
    //         leptonTwoRecoWeight = effs.first/effs.second;
    //         leptonTwoRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

    //         // trigger weights with trigger matching:
    //         //
    //         // check if lepton could pass the trigger threshold and is matched
    //         // to a trigger object.  When both muons pass the trigger, use the
    //         // efficiency for detecting either

    //         if (triggered.all()) {
    //             effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
    //             effCont2      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
    //             triggerWeight = GetTriggerSF(effCont1, effCont2);
    //             triggerVar    = GetTriggerSFError(effCont1, effCont2);
    //         } else if (triggered.test(0)) {
    //             effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
    //             effs          = effCont1.GetEff();
    //             errs          = effCont1.GetErr();
    //             triggerWeight = effs.first/effs.second;
    //             triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
    //         } else if (triggered.test(1)) {
    //             effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
    //             effs          = effCont1.GetEff();
    //             errs          = effCont1.GetErr();
    //             triggerWeight = effs.first/effs.second;
    //             triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
    //         } else {
    //             return kTRUE;
    //         }

    //         // update the event weight
    //         eventWeight *= triggerWeight;
    //         eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;
    //     }
    // } else if (electrons.size() == 2 && muons.size() == 0 && taus.size() == 1) { // e+e+tau selection

    //     // if (nJetsCut >= 2 ){// && nBJetsCut >= 1){
    //     //     channel = "ee";
    //     // } else if (nJetsCut2 == 0 ){
    //     //     //channel = "ee0";
    //     //     return kTRUE;
    //     // } else {
    //     //     return kTRUE;
    //     // }

    //     channel = "ee_tau";
    //     eventCounts[channel]->Fill(1);


    //     if (electrons[0]->pt < 30 || electrons[1]->pt < 10 || !electronTriggered)
    //         return kTRUE;
    //     eventCounts[channel]->Fill(2);

    //     TLorentzVector electronOneP4, electronTwoP4, dielectronP4;
    //     electronOneP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
    //     electronTwoP4.SetPtEtaPhiM(electrons[1]->pt, electrons[1]->eta, electrons[1]->phi, 511e-6);
    //     dielectronP4 = electronOneP4 + electronTwoP4;
    //     if (dielectronP4.M() < 12.)
    //         return kTRUE;
    //     eventCounts[channel]->Fill(3);

    //     // if (nJetsCut < 2 || nBJetsCut < 1)
    //     //     return kTRUE;
    //     // eventCounts[channel]->Fill(4);

    //     leptonOneP4     = electronOneP4;
    //     leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
    //     leptonOneFlavor = 11*electrons[0]->q;
    //     leptonOneDZ     = electrons[0]->dz;
    //     leptonOneD0     = electrons[0]->d0;

    //     leptonTwoP4     = electronTwoP4;
    //     leptonTwoIso    = GetElectronIsolation(electrons[1], fInfo->rhoJet);
    //     leptonTwoFlavor = 11*electrons[1]->q;
    //     leptonTwoDZ     = electrons[1]->dz;
    //     leptonTwoD0     = electrons[1]->d0;

    //     // trigger weights with trigger matching
    //     bitset<2> triggered;
    //     for (const auto& name: passTriggerNames) {
    //         if (trigger->passObj(name, 1, electrons[0]->hltMatchBits) && electronOneP4.Pt() > 30)
    //             triggered.set(0);
    //         if (trigger->passObj(name, 1, electrons[1]->hltMatchBits) && electronTwoP4.Pt() > 30)
    //             triggered.set(1);
    //     }

    //     if (!triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 0;
    //     if ( triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 1;
    //     if (!triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 2;
    //     if ( triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 3;


    //     ////////////////////////////////////////
    //     // tau reconstruct info
    //     TLorentzVector tau0P4;
    //     tau0P4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, 1.77682);

    //     tauP4             = tau0P4;

    //     tauDecayMode      = taus[0]->decaymode;
    //     tauMVA            = taus[0]->rawIsoMVA3newDMwLT;
    //     tauPuppiChHadIso  = taus[0]->puppiChHadIso;
    //     tauPuppiGammaIso  = taus[0]->puppiGammaIso;
    //     tauPuppiNeuHadIso = taus[0]->puppiNeuHadIso;

    //     pair <float, float> tauVetoedJetPtPair;
    //     tauVetoedJetPtPair= GetTauVetoedJetPt(tauP4, vetoedJets);
    //     tauVetoedJetPt    = tauVetoedJetPtPair.first;
    //     tauVetoedJetPtUnc = tauVetoedJetPtPair.second;

    //     ////////////////////////////////////////

    //     if (!isData) {
    //         leptonOneMother = GetGenMotherId(genParticles, electronOneP4);
    //         leptonTwoMother = GetGenMotherId(genParticles, electronTwoP4);
    //         tauGenFlavor    = GetTauGenFlavor(tau0P4, genTausHad, vetoedJets, false); // useHadronFlavor = false
    //         tauGenFlavorHad = GetTauGenFlavor(tau0P4, genTausHad, vetoedJets, true);  // useHadronFlavor = true

    //         EfficiencyContainer effCont1, effCont2;
    //         effCont1 = weights->GetElectronRecoEff(electronOneP4);
    //         effCont2 = weights->GetElectronRecoEff(electronTwoP4);

    //         pair<float, float> effs, errs;
    //         effs = effCont1.GetEff();
    //         errs = effCont1.GetErr();
    //         leptonOneRecoWeight = effs.first/effs.second;
    //         leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

    //         effs = effCont2.GetEff();
    //         errs = effCont2.GetErr();
    //         leptonTwoRecoWeight = effs.first/effs.second;
    //         leptonTwoRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));


    //         if (triggered.all()) {
    //             effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronOneP4);
    //             effCont2      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronTwoP4);
    //             triggerWeight = GetTriggerSF(effCont1, effCont2);
    //             triggerVar    = GetTriggerSFError(effCont1, effCont2);
    //         } else if (triggered.test(0)) {
    //             effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronOneP4);
    //             effs          = effCont1.GetEff();
    //             errs          = effCont1.GetErr();
    //             triggerWeight = effs.first/effs.second;
    //             triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
    //         } else if (triggered.test(1)) {
    //             effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronTwoP4);
    //             effs          = effCont1.GetEff();
    //             errs          = effCont1.GetErr();
    //             triggerWeight = effs.first/effs.second;
    //             triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
    //         } else {
    //             return kTRUE;
    //         }

    //         // update event weight
    //         eventWeight *= triggerWeight;
    //         eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;
    //     }
    // } else if (electrons.size() == 1 && muons.size() == 1 && taus.size() == 1) { // e+mu+tau
    //     channel = "emu_tau";
    //     eventCounts[channel]->Fill(1);

    //     // trigger matching for thresholds
        
    //     float muPtThreshold = 10.;
    //     float elPtThreshold = 10.;
    //     if ( muonTriggered ) {
    //         muPtThreshold = 25.;
    //     } else if (electronTriggered) {
    //         elPtThreshold = 30.;
    //     }

    //     if (
    //             muons[0]->pt < muPtThreshold 
    //             || electrons[0]->pt < elPtThreshold
    //             || !(muonTriggered || electronTriggered)
    //        )
    //         return kTRUE;
    //     eventCounts[channel]->Fill(2);

    //     TLorentzVector muonP4, electronP4, dilepton;
    //     muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
    //     electronP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
    //     dilepton = muonP4 + electronP4;
    //     if (dilepton.M() < 12)
    //         return kTRUE;
    //     eventCounts[channel]->Fill(3);

    //     // if (nJetsCut < 2) // || nBJetsCut < 1)
    //     //     return kTRUE;
    //     // eventCounts[channel]->Fill(4);

    //     leptonOneP4     = muonP4;
    //     leptonOneIso    = GetMuonIsolation(muons[0]);
    //     leptonOneFlavor = 13*muons[0]->q;
    //     leptonOneDZ     = muons[0]->dz;
    //     leptonOneD0     = muons[0]->d0;

    //     leptonTwoP4     = electronP4;
    //     leptonTwoIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
    //     leptonTwoFlavor = 11*electrons[0]->q;
    //     leptonTwoDZ     = electrons[0]->dz;
    //     leptonTwoD0     = electrons[0]->d0;

    //     // trigger weights with trigger matching
    //     bitset<2> triggered;
    //     for (const auto& name: passTriggerNames) {
    //         if (trigger->passObj(name, 1, muons[0]->hltMatchBits))
    //             triggered.set(0);
    //         if (trigger->passObj(name, 1, electrons[0]->hltMatchBits))
    //             triggered.set(1);
    //     }

    //     if (!triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 0;
    //     if ( triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 1;
    //     if (!triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 2;
    //     if ( triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 3;

    //     ////////////////////////////////////////
    //     // tau reconstruct info
    //     TLorentzVector tau0P4;
    //     tau0P4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, 1.77682);

    //     tauP4             = tau0P4;

    //     tauDecayMode      = taus[0]->decaymode;
    //     tauMVA            = taus[0]->rawIsoMVA3newDMwLT;
    //     tauPuppiChHadIso  = taus[0]->puppiChHadIso;
    //     tauPuppiGammaIso  = taus[0]->puppiGammaIso;
    //     tauPuppiNeuHadIso = taus[0]->puppiNeuHadIso;

    //     pair <float, float> tauVetoedJetPtPair;
    //     tauVetoedJetPtPair= GetTauVetoedJetPt(tauP4, vetoedJets);
    //     tauVetoedJetPt    = tauVetoedJetPtPair.first;
    //     tauVetoedJetPtUnc = tauVetoedJetPtPair.second;
    //     ////////////////////////////////////////

    //     if (!isData) {
    //         leptonOneMother = GetGenMotherId(genParticles, muonP4);
    //         leptonTwoMother = GetGenMotherId(genParticles, electronP4);
    //         tauGenFlavor    = GetTauGenFlavor(tau0P4, genTausHad, vetoedJets, false); // useHadronFlavor = false
    //         tauGenFlavorHad = GetTauGenFlavor(tau0P4, genTausHad, vetoedJets, true);  // useHadronFlavor = true
            

    //         // reconstruction weights
    //         EfficiencyContainer effCont1, effCont2;
    //         effCont1 = weights->GetMuonRecoEff(muonP4);
    //         effCont2 = weights->GetElectronRecoEff(electronP4);

    //         pair<float, float> effs, errs;
    //         effs = effCont1.GetEff();
    //         errs = effCont1.GetErr();
    //         leptonOneRecoWeight = effs.first/effs.second;
    //         leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

    //         effs = effCont2.GetEff();
    //         errs = effCont2.GetErr();
    //         leptonTwoRecoWeight = effs.first/effs.second;
    //         leptonTwoRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));


    //         if (triggered.all()) {
    //             effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
    //             effCont2      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronP4);
    //             triggerWeight = GetTriggerSF(effCont1, effCont2);
    //             triggerVar    = GetTriggerSFError(effCont1, effCont2);
    //         } else if (triggered.test(0)) {
    //             effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
    //             effs          = effCont1.GetEff();
    //             errs          = effCont1.GetErr();
    //             triggerWeight = effs.first/effs.second;
    //             triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
    //         } else if (triggered.test(1)) {
    //             effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronP4);
    //             effs          = effCont1.GetEff();
    //             errs          = effCont1.GetErr();
    //             triggerWeight = effs.first/effs.second;
    //             triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
    //         } else {
    //             return kTRUE;
    //         }

    //         // update event weight
    //         eventWeight *= triggerWeight;
    //         eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;
    //     }
    
    } else {
        return kTRUE;
    }


    ///////////////////
    // Fill jet info //
    ///////////////////

    if (channel != "mu4j" && channel != "e4j" && channel != "mu4j_fakes" && channel != "e4j_fakes") { // jets are handled differently for the lepton + jet selection
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

vector<unsigned> MultileptonAnalyzer::PairDileptonToZ(vector<TLorentzVector> leptons)
{
    unsigned pairings[3][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}};
    float minMassDiff = 1e6;
    unsigned bestPairing[4] = {0, 1, 2, 3};
    for (unsigned i = 0; i < 3; ++i) {
        unsigned *p = pairings[i];
        TLorentzVector dilepton1 = leptons[p[0]] + leptons[p[1]];
        TLorentzVector dilepton2 = leptons[p[2]] + leptons[p[3]];

        if (fabs(dilepton1.M() - 91.2) < minMassDiff) {
            minMassDiff = fabs(dilepton1.M() - 91.2);

            // This makes me loathe c++
            bestPairing[0] = p[0];
            bestPairing[1] = p[1];
            bestPairing[2] = p[2];
            bestPairing[3] = p[3];
        }
        if (fabs(dilepton2.M() - 91.2) < minMassDiff) {
            minMassDiff = fabs(dilepton2.M() - 91.2);
            bestPairing[0] = p[2];
            bestPairing[1] = p[3];
            bestPairing[2] = p[0];
            bestPairing[3] = p[1];
        }
    }

    // this too
    vector<unsigned> outPairing(bestPairing, bestPairing + sizeof(bestPairing)/sizeof(unsigned));
    return outPairing;
}

int MultileptonAnalyzer::GetGenMotherId(vector<TGenParticle*> particles, TLorentzVector p4)
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

int MultileptonAnalyzer::GetTauGenFlavor(TLorentzVector p4, vector<TGenParticle*> genTausHad, vector<TJet*> vetoedJets, bool useHadronFlavor )
{
    int flavor = 26;
    
    
    // check if can be tagged as hadronic tau
    for (unsigned i = 0; i < genTausHad.size(); ++i) {
        TLorentzVector genP4;
        genP4.SetPtEtaPhiM(genTausHad[i]->pt, genTausHad[i]->eta, genTausHad[i]->phi, genTausHad[i]->mass); 
        if (genP4.DeltaR(p4) < 0.3) {
            flavor = 15;
        }
    }

    // check if can be tagged by jet flavor
    if (flavor!= 15){
        float jetPtMax = - 1.0;
        for (unsigned i = 0; i < vetoedJets.size(); ++i) {


            TLorentzVector jetP4; 
            jetP4.SetPtEtaPhiM(vetoedJets[i]->pt, vetoedJets[i]->eta, vetoedJets[i]->phi, vetoedJets[i]->mass);
            if (jetP4.DeltaR(p4) < 0.4 && jetP4.Pt()>jetPtMax) {
                // cout<<"flavor Jet: "<< jetPtMax<<','<< jetP4.Pt()<< endl;
                if(useHadronFlavor) {
                    flavor = vetoedJets[i]->hadronFlavor;
                } else {
                    flavor = abs(vetoedJets[i]->partonFlavor);
                }
                jetPtMax = jetP4.Pt();
            }
        }
    }
    return flavor;
}


pair<float, float> MultileptonAnalyzer::GetTauVetoedJetPt(TLorentzVector p4, vector<TJet*> vetoedJets)
{

    float jetPt = -1;
    float jetPtUnc = -1;
    

    // check if can be tagged by jet flavor

    float jetPtMax = - 1.0;
    for (unsigned i = 0; i < vetoedJets.size(); ++i) {

        TLorentzVector jetP4; 
        jetP4.SetPtEtaPhiM(vetoedJets[i]->pt, vetoedJets[i]->eta, vetoedJets[i]->phi, vetoedJets[i]->mass);

        if (jetP4.DeltaR(p4) < 0.4 && jetP4.Pt()>jetPtMax) {
            // cout<<"Energy Jet: "<<jetPtMax<<','<< jetP4.Pt() << endl;

            // save 
            jetPt       = jetP4.Pt();
            jetPtUnc    = jetPt * float(particleSelector->JetUncertainty(vetoedJets[i]));
            jetPtUnc    = abs(jetPtUnc);

            jetPtMax = jetP4.Pt();
        }    
    }
    
    return make_pair(jetPt,jetPtUnc);
}


 
vector<TJet*> MultileptonAnalyzer::KinematicTopTag(vector<TJet*> jets, TVector2 metP2, TLorentzVector leptonP4)
{
    float topMass = 172.4;
    float wMass   = 80.385;
    float minChi2 = 1e9;

    vector<int> ix(jets.size());
    std::iota(std::begin(ix), std::end(ix), 0);
    vector<int> optimalOrder = ix;
    do {
        if (ix[0] > ix[1]) continue;

        vector<TLorentzVector> jetP4;
        for (unsigned i = 0; i < 4; ++i) {
            TLorentzVector p;
            p.SetPtEtaPhiM(jets[ix[i]]->pt, jets[ix[i]]->eta, jets[ix[i]]->phi, jets[ix[i]]->mass);
            jetP4.push_back(p);
            //cout << ix[i];
        }

        TLorentzVector wCandidate    = jetP4[0] + jetP4[1];
        TLorentzVector hTopCandidate = wCandidate + jetP4[2];
        TLorentzVector lTopCandidate = leptonP4 + jetP4[3];

        // hadronic top mass tag
        float chi2 = pow(wCandidate.M() - wMass, 2) + pow(hTopCandidate.M() - topMass, 2);

        // leptonic top mass tag

        // b jet assignment

        // angular correlations

        if (chi2 < minChi2) {
            minChi2 = chi2;
            optimalOrder = ix;
        }
        //cout << " " << chi2 << endl;
        std::reverse(ix.begin() + 4, ix.end());
    } while (next_permutation(ix.begin(), ix.end()));
    //cout << minChi2 << "\n\n";

    vector<TJet*> newJets;
    for (unsigned i = 0; i < 4; ++i) {
        newJets.push_back(jets[optimalOrder[i]]);
    }
    return newJets;
}

float MultileptonAnalyzer::GetTriggerSF(EfficiencyContainer eff1, EfficiencyContainer eff2)
{
    pair<double, double> trigEff1, trigEff2;
    trigEff1 = eff1.GetEff();
    trigEff2 = eff2.GetEff();

    float sf = 1.;
    if (trigEff1.second > 0 && trigEff2.second <= 0){
        sf *= trigEff1.first/trigEff1.second;
    }
    
    if (trigEff1.second <= 0 && trigEff2.second > 0){
        sf *= trigEff2.first/trigEff2.second;
    }
    if (trigEff1.second > 0 && trigEff2.second > 0){
        sf *= trigEff1.first/trigEff1.second;
        sf *= trigEff2.first/trigEff2.second;
    }

    //if (trigEff1.second > 0 || trigEff2.second > 0) {
    //  sf = (1 - (1 - trigEff1.first)*(1 - trigEff2.first))/(1 - (1 - trigEff1.second)*(1 - trigEff2.second));
    //}


    return sf;
}

float MultileptonAnalyzer::GetTriggerSFError(EfficiencyContainer eff1, EfficiencyContainer eff2)
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

void MultileptonAnalyzer::ResetJetCounters()
{
    nJets = nBJets = nFwdJets = 0;
    nJetsCut       = nBJetsCut        = 0;
    nJetsJESUp     = nJetsJESDown     = 0;
    nJetsJERUp     = nJetsJERDown     = 0;
    nBJetsJESUp    = nBJetsJESDown    = 0;
    nBJetsJERUp    = nBJetsJERDown    = 0;
    nBJetsBTagUp   = nBJetsBTagDown   = 0;
    nBJetsMistagUp = nBJetsMistagDown = 0;

    nJetsCut2 = 0;
}

void MultileptonAnalyzer::JetCounting(TJet* jet, float jerc_nominal, float resRand)
{
    float jetPt = jet->pt;

    float rNumber = rng->Uniform(1.);
    if (jet->pt > 30) {

        // nominal
        ++nJets;
        if (particleSelector->BTagModifier(jet, "MVAT", 0, 0, rNumber)) ++nBJets;

        // b tag up
        if (particleSelector->BTagModifier(jet, "MVAT", 1, 0, rNumber)) ++nBJetsBTagUp;
             
        // b tag down
        if (particleSelector->BTagModifier(jet, "MVAT", -1, 0, rNumber))  ++nBJetsBTagDown;

        // misttag up
        if (particleSelector->BTagModifier(jet, "MVAT", 0, 1, rNumber)) ++nBJetsMistagUp;
        
        // mistag down
        if (particleSelector->BTagModifier(jet, "MVAT", 0, -1, rNumber)) ++nBJetsMistagDown;
    }

    // JES up
    double jec = particleSelector->JetCorrector(jet, "NONE");
    float jecUnc = particleSelector->JetUncertainty(jet);
    jet->pt = jet->ptRaw*jec*(1 + jecUnc)*jerc_nominal;
    if (jet->pt > 30) {
        ++nJetsJESUp;
        if (particleSelector->BTagModifier(jet, "MVAT", 0, 0, rNumber)) { 
            ++nBJetsJESUp;
        } 
    }

    // JES down
    jet->pt = jet->ptRaw*jec*(1 - jecUnc)*jerc_nominal;
    if (jet->pt > 30) {
        ++nJetsJESDown;
        if (particleSelector->BTagModifier(jet, "MVAT", 0, 0, rNumber)) { 
            ++nBJetsJESDown;
        } 
    }

    // JER up
    pair<float, float> resPair = particleSelector->JetResolutionAndSF(jet, 1);
    float jerc = 1 + resRand*sqrt(std::max((double)resPair.second*resPair.second - 1, 0.));
    jet->pt = jet->ptRaw*jec*jerc;
    if (jet->pt > 30) {
        ++nJetsJERUp;
        if (particleSelector->BTagModifier(jet, "MVAT", 0, 0, rNumber)) { 
            ++nBJetsJERUp;
        }     
    }

    // JER down
    resPair     = particleSelector->JetResolutionAndSF(jet, -1);
    jerc        = 1 + resRand*sqrt(std::max((double)resPair.second*resPair.second - 1, 0.));
    jet->pt     = jet->ptRaw*jec*jerc;
    if (jet->pt > 30) {
        ++nJetsJERDown;
        if (particleSelector->BTagModifier(jet, "MVAT", 0, 0, rNumber)) { 
            ++nBJetsJERDown;
        } 
    }
    jet->pt = jetPt;

}