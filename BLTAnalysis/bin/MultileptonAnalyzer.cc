#include "MultileptonAnalyzer.h"
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
    return lhs->csv > rhs->csv;
}


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

    // counts W decay modes for acceptance calculation
    outHistName = params->get_output_treename("GenCategory");
    hGenCat = new TH1D(outHistName.c_str(), "WW decay modes",30,0.5,30.5);

    // record PDF variations to generated number of events based on MC replicas
    outHistName = params->get_output_treename("var_PDF");
    pdfCountsInit = new TH1D(outHistName.c_str(), "pdf variations (initial)", 101, 0.5, 101.5);
    outHistName = params->get_output_treename("var_QCD");
    qcdCountsInit = new TH1D(outHistName.c_str(), "qcd variations (initial)", 9, 0.5, 9.5);

    // saving pdf variations
    outHistName = params->get_output_treename("var_PDF_jets_init");
    pdfCountsJetsInit    = new TH2D(outHistName.c_str(),"pdf variations", 101, 0.5, 101.5, 4, -0.5, 3.5);
    outHistName = params->get_output_treename("var_PDF_partons_init");
    pdfCountsPartonsInit = new TH2D(outHistName.c_str(),"pdf variations", 101, 0.5, 101.5, 4, -0.5, 3.5);

    outHistName = params->get_output_treename("var_QCD_jets_init");
    qcdCountsJetsInit    = new TH2D(outHistName.c_str(),"qcd variations", 9, 0.5, 9.5, 4, -0.5, 3.5);
    outHistName = params->get_output_treename("var_QCD_partons_init");
    qcdCountsPartonsInit = new TH2D(outHistName.c_str(),"qcd variations", 9, 0.5, 9.5, 4, -0.5, 3.5);

    vector<std::string> channelNames = {"mumu", "ee", "emu", 
                                        "etau", "mutau", 
                                        "e4j", "mu4j", "tau4j",
                                        "mutau_fakes", "mu4j_fakes",
                                        "etau_fakes", "e4j_fakes"
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
        tree->Branch("rPV", &rPV);

        // weights and their uncertainties
        tree->Branch("eventWeight", &eventWeight);
        tree->Branch("leptonOneRecoWeight", &leptonOneRecoWeight);
        tree->Branch("leptonTwoRecoWeight", &leptonTwoRecoWeight);
        tree->Branch("leptonOneIDWeight", &leptonOneIDWeight);
        tree->Branch("leptonTwoIDWeight", &leptonTwoIDWeight);
        tree->Branch("topPtWeight", &topPtWeight);
        tree->Branch("puWeight", &puWeight);
        tree->Branch("triggerWeight", &triggerWeight);
        tree->Branch("genWeight", &genWeight);

        tree->Branch("leptonOneRecoVar", &leptonOneRecoVar);
        tree->Branch("leptonTwoRecoVar", &leptonTwoRecoVar);
        tree->Branch("leptonOneIDVar", &leptonOneIDVar);
        tree->Branch("leptonTwoIDVar", &leptonTwoIDVar);
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
        tree->Branch("leptonOnePtCorr", &leptonOnePtCorr);
        tree->Branch("leptonOneIso", &leptonOneIso);
        tree->Branch("leptonOneFlavor", &leptonOneFlavor);
        tree->Branch("leptonOneMother", &leptonOneMother);
        tree->Branch("leptonOneGenId", &leptonOneGenId);
        tree->Branch("leptonOneD0", &leptonOneD0);
        tree->Branch("leptonOneDZ", &leptonOneDZ);

        tree->Branch("leptonTwoP4", &leptonTwoP4);
        tree->Branch("leptonTwoPtCorr", &leptonTwoPtCorr);
        tree->Branch("leptonTwoIso", &leptonTwoIso);
        tree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
        tree->Branch("leptonTwoMother", &leptonTwoMother);
        tree->Branch("leptonTwoGenId", &leptonTwoGenId);
        tree->Branch("leptonTwoD0", &leptonTwoD0);
        tree->Branch("leptonTwoDZ", &leptonTwoDZ);

        if (channel == "mutau" || channel == "etau" 
                || channel == "mutau_fakes" || channel == "etau_fakes"
                || channel == "tau4j"
                ) {
            //tree->Branch("tauChHadMult",  &tauChHadMult);
            //tree->Branch("tauPhotonMult", &tauPhotonMult);
            tree->Branch("tauDecayMode",  &tauDecayMode);
            tree->Branch("tauMVA",        &tauMVA);
        }

        // jets
        tree->Branch("jetOneP4", &jetOneP4);
        tree->Branch("jetOneTag", &jetOneTag);
        tree->Branch("jetOneFlavor", &jetOneFlavor);
        tree->Branch("jetTwoP4", &jetTwoP4);
        tree->Branch("jetTwoTag", &jetTwoTag);
        tree->Branch("jetTwoFlavor", &jetTwoFlavor);

        if (channel == "mu4j" || channel == "mu4j_fakes" 
                || channel == "e4j" || channel == "e4j_fakes"
                || channel == "tau4j") {
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
        tree->Branch("nBJetsRaw", &nBJetsRaw);

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
        tree->Branch("nBJetsCTagUp", &nBJetsCTagUp);
        tree->Branch("nBJetsCTagDown", &nBJetsCTagDown);
        tree->Branch("nBJetsMistagUp", &nBJetsMistagUp);
        tree->Branch("nBJetsMistagDown", &nBJetsMistagDown);

        // theory uncertainties
        tree->Branch("qcdWeights", &qcdWeights);
        tree->Branch("pdfWeight", &pdfWeight);
        tree->Branch("alphaS", &alphaS);

        outTrees[channel] = tree;

        // event counter
        outHistName = params->get_output_treename("TotalEvents_" + channel);
        eventCounts[channel] = new TH1D(outHistName.c_str(),"ChannelCounts",10,0.5,10.5);

        // saving pdf variations
        outHistName = params->get_output_treename("var_PDF_jets_" + channel);
        pdfCountsJets[channel]    = new TH2D(outHistName.c_str(),"pdf variations", 101, 0.5, 101.5, 5, -0.5, 4.5);
        outHistName = params->get_output_treename("var_PDF_partons_" + channel);
        pdfCountsPartons[channel] = new TH2D(outHistName.c_str(),"pdf variations", 101, 0.5, 101.5, 4, -0.5, 3.5);

        outHistName = params->get_output_treename("var_QCD_jets_" + channel);
        qcdCountsJets[channel]    = new TH2D(outHistName.c_str(),"qcd variations", 9, 0.5, 9.5, 5, -0.5, 4.5);
        outHistName = params->get_output_treename("var_QCD_partons_" + channel);
        qcdCountsPartons[channel] = new TH2D(outHistName.c_str(),"qcd variations", 9, 0.5, 9.5, 4, -0.5, 3.5);
    }

    // initialize jet counters
    nJets = nBJets = nBJetsRaw = 0;
    nJetsCut     = nBJetsCut = 0;
    nJetsJERUp   = nJetsJERDown   = nBJetsJERUp    = nBJetsJERDown    = 0;
    nBJetsCTagUp = nBJetsCTagDown = nBJetsMistagUp = nBJetsMistagDown = 0;

    vector<unsigned> counters(particleSelector->GetJECSourceNames().size(), 0);
    nJetsJESUp    = counters;
    nJetsJESDown  = counters;
    nBJetsJESUp   = counters;
    nBJetsJESDown = counters;

    vector<unsigned> bcounters(particleSelector->GetBTagSourceNames().size(), 0);
    nBJetsBTagUp   = bcounters;
    nBJetsBTagDown = bcounters;

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
    vector<int> genMotherId;
    vector<float> pdfVariations;

    unsigned wCount = 0;
    unsigned tauCount = 0;
    unsigned topCount = 0;
    bitset<6> wDecay;
    bitset<4> tauDecay;
    float topSF = 1.;
    topPtWeight = 1.;
    topPtVar = 0.;
    nPartons = 0;
    if (!isData) {

        // loop over gen particle collection:
        //   * parton counting for combining inclusive and jet-binned Drell-Yan samples
        //   * categorization of W and tau decays in ttbar/tW events
        //   * calculate scale factors for top pt reweighting
        unsigned count = 0;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
            unsigned pid = abs(particle->pdgId);


            //if (abs(particle->pdgId) == 6) {
                //cout << i  << ", " << particle->parent << ", " << particle->pdgId << ", " << particle->status;
                //cout << "\t" << particle->pt << ", " << particle->eta;
                //cout << endl;
            //}

            // parton counting for jet-binned Drell-Yan samples
            if (
                    particle->status == 23 
                    && (pid < 6 || pid == 21) 
                    && particle->parent != -2
               ) {
                ++count;
            }

            // top pt reweighting: get the scale factor based on the top quark pt
            if (pid == 6 && particle->status == 62) {
                topSF *= exp(0.0615 - 0.0005*particle->pt);
                ++topCount;
            }

            // Tag events based on neutrino flavors
            unsigned flavor = pid;
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
            //if (abs(particle->pdgId) == 15) {
            //    cout << particle->status << " " << particle->parent << endl;
            //}
            if (
                    (pid == 11 || pid == 13 || pid == 15)
                    && particle->parent != -2
               ) {

                // Find promptly produced leptons (decaying from a W or Z) or
                // e/mu originating from taus
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                //cout << particle->pdgId << ", " << mother->pdgId << endl; //", ";
                int origin = abs(mother->pdgId);
                bool isPrompt = false;
                while (mother->parent != -2) {
                    int momID = abs(mother->pdgId);
                    if (momID == 23 || momID == 24) { 
                        isPrompt = true;
                        if (((pid == 11 || pid == 13) && origin != 15) || pid == 15) {
                            origin = momID;
                        }                         
                        break;
                    } else if (momID == 15 && (pid == 11 || pid == 13)) {
                        origin = 15;
                    } else if (momID < 6 || momID == 22) {
                        break;
                    }
                    mother = (TGenParticle*) fGenParticleArr->At(mother->parent);
                }

                if (isPrompt) {
                    //cout << abs(particle->pdgId) << " " << origin << endl;
                    genParticles.push_back(particle);
                    genMotherId.push_back(origin);
                }
            }
        }
        nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY

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
                //cout << id << " " << lheWeight << endl;
                if (id >= 1001 && id <= 1009) {
                    qcdWeights.push_back(lheWeight);
                    qcdCountsInit->Fill(id - 1000, lheWeight);
                } else if (id >= 2001 && id <= 2100) {      
                    pdfWeight += pow(qcdWeights[0] - lheWeight, 2.);
                    pdfVariations.push_back(lheWeight);
                    pdfCountsInit->Fill(id - 2000, lheWeight);
                } else if (id == 2101 || id == 2102) {
                    alphaS = lheWeight;
                    pdfVariations.push_back(lheWeight);
                    pdfCountsInit->Fill(101, lheWeight);
                }
            }
        }


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
            topPtWeight = sqrt(topSF);
            topPtVar    = topPtWeight; // the official prescription is to use no topPt weight as the down variation and twice the weight as the up variation
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

        // record theory uncertainties before cuts in differenct parton bins
        for (unsigned i = 0; i < qcdWeights.size(); ++i) {
            qcdCountsPartonsInit->Fill(i+1, nPartons, qcdWeights[i]);
        }

        for (unsigned i = 0; i < pdfVariations.size(); ++i) {
            pdfCountsPartonsInit->Fill(i+1, nPartons, pdfVariations[i]);
        }   

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
            if (
                params->selection == "single_lepton" 
                && params->datasetgroup.substr(0, 8)  == "electron" 
                && (triggerNames[i] == "HLT_IsoMu24_v*" || triggerNames[i] == "HLT_IsoTkMu24_v*")
               ) {
                passTrigger = false;
                break;
            }
        }
    }

    if (isData && !passTrigger)
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
        muon->ptErr = muonSF; // I'm not using ptErr so put the Roch Corr there 
        muonP4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);

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
    map<int, float> electronScaleCorr;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        if (isData) {
            scaleData sdata = electronScaler->GetScaleData(electron, runNumber);
            electron->pt    *= sdata.scale;
            electron->pfPt   = sdata.scale; // since it's not used anywhere else...
        } else {
            float sFactor = electronScaler->GetSmearingFactor(electron, 0, 0);
            float eScale = rng->Gaus(1, sFactor);
            electron->pt    *= eScale; 
            electron->pfPt   = eScale; // since it's not used anywhere else...
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

    /* TAUS */
    vector<TTau*> taus;
    vector<TLorentzVector> veto_taus;
    for (int i=0; i < fTauArr->GetEntries(); i++) {
        TTau *tau = (TTau*) fTauArr->At(i);
        assert(tau);

        TLorentzVector tauP4; 
        tauP4.SetPtEtaPhiM(tau->pt, tau->eta, tau->phi, tau->m);

        // Prevent overlap of muons or electrons and taus
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
                && (tau->hpsDisc & baconhep::kByTightIsolationMVA3oldDMwLT)
                && (tau->hpsDisc & baconhep::kByMVA6VTightElectronRejection)
                && (tau->hpsDisc & baconhep::kByTightMuonRejection3)
          ) {
            taus.push_back(tau);
            veto_taus.push_back(tauP4);
            //cout << tau->m << ", " << tau->decaymode << endl;
        }
    }
    sort(taus.begin(), taus.end(), sort_by_higher_pt<TTau>);

    /* JETS */
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets;
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
                        if (jet->csv > 0.8484) { 
                            ++nBJets;
                        } 
                    }
                } else {
                    JetCounting(jet, jerc, gRand);
                    if (jet->pt > 30 && jet->csv > 0.8484) { 
                        ++nBJetsRaw;
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

    // use the highest jet multiplicities given all systematic variations
    if (!isData) {
        vector<unsigned> nJetList {nJets, nJetsJERUp, nJetsJERDown};
        nJetList.insert(nJetList.end(), nJetsJESUp.begin(), nJetsJESUp.end());
        nJetList.insert(nJetList.end(), nJetsJESDown.begin(), nJetsJESDown.end());
        nJetsCut = *max_element(begin(nJetList), end(nJetList));

        vector<unsigned> nBJetList {nBJets, nBJetsRaw,
                                    nBJetsJERUp, nBJetsJERDown, 
                                    nBJetsMistagUp, nBJetsMistagDown,
                                    nBJetsCTagUp, nBJetsCTagDown
                                   };

        nBJetList.insert(nBJetList.end(), nBJetsJESUp.begin(), nBJetsJESUp.end());
        nBJetList.insert(nBJetList.end(), nBJetsJESDown.begin(), nBJetsJESDown.end());
        nBJetList.insert(nBJetList.end(), nBJetsBTagUp.begin(), nBJetsBTagUp.end());
        nBJetList.insert(nBJetList.end(), nBJetsBTagDown.begin(), nBJetsBTagDown.end());
        nBJetsCut = *max_element(begin(nBJetList), end(nBJetList));
    } else {
        nJetsCut = nJets;
        nBJetsCut = nBJets;
    }

    //if ((nJetsJESUp > nJets && nJetsJESDown > nJets) || (nJetsJESUp > nJets && nJetsJESDown > nJets)) {
    //    cout << nJets << ", " << nJetsJESUp << ", " << nJetsJESDown << endl;
    //}

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
    if (muons.size() == 2 && electrons.size() == 0 && taus.size() == 0) { // mu+mu selection
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

        // save PDF MC replica variation by jet multiplicity (do this before jet selection)
        if (!isData && fLHEWeightArr != 0) {
            FillPDFHist(pdfVariations, channel);
        }

        if (nJetsCut < 2) //|| nBJetsCut < 1)
            return kTRUE;
        eventCounts[channel]->Fill(4);

        leptonOneP4     = muonOneP4;
        leptonOneIso    = muonOneIso;
        leptonOnePtCorr = muons[0]->ptErr;
        leptonOneFlavor = muons[0]->q*13;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;

        leptonTwoP4     = muonTwoP4;
        leptonTwoIso    = muonTwoIso;
        leptonTwoPtCorr = muons[1]->ptErr;
        leptonTwoFlavor = muons[1]->q*13;
        leptonTwoDZ     = muons[1]->dz;
        leptonTwoD0     = muons[1]->d0;

        if (!isData) {
            pair<int, int> pdgId;
            pdgId           = GetGenId(genParticles, genMotherId, muonOneP4);
            leptonOneGenId  = pdgId.first;
            leptonOneMother = pdgId.second;
            pdgId           = GetGenId(genParticles, genMotherId, muonTwoP4);
            leptonTwoGenId  = pdgId.first;
            leptonTwoMother = pdgId.second;

            // ID and ISO weights
            EfficiencyContainer effCont;
            pair<float, float> effs, errs;

            effCont = weights->GetMuonIDEff(muonOneP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonOneIDWeight = effs.first/effs.second;
            leptonOneIDVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            effCont = weights->GetMuonIDEff(muonTwoP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonTwoIDWeight = effs.first/effs.second;
            leptonTwoIDVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            effCont = weights->GetMuonISOEff(muonOneP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            effCont = weights->GetMuonISOEff(muonTwoP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonTwoRecoWeight = effs.first/effs.second;
            leptonTwoRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

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
                EfficiencyContainer effCont1, effCont2;
                effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
                effCont2      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
                triggerWeight = GetTriggerSF(effCont1, effCont2);
                triggerVar    = GetTriggerSFError(effCont1, effCont2);
            } else if (triggered.test(0)) {
                effCont       = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonOneP4);
                effs          = effCont.GetEff();
                errs          = effCont.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else if (triggered.test(1)) {
                effCont       = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonTwoP4);
                effs          = effCont.GetEff();
                errs          = effCont.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }

            //cout << leptonOneRecoVar << ", " << leptonTwoRecoVar << ", " << triggerVar << endl;

            // update the event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneIDWeight*leptonTwoIDWeight*leptonOneRecoWeight*leptonTwoRecoWeight;
        }
    }  else if (electrons.size() == 2 && muons.size() == 0 && taus.size() == 0) { // e+e selection
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

        // save PDF MC replica variation by jet multiplicity (do this before jet selection)
        if (!isData && fLHEWeightArr != 0) {
            FillPDFHist(pdfVariations, channel);
        }

        if (nJetsCut < 2) //|| nBJetsCut < 1)
            return kTRUE;
        eventCounts[channel]->Fill(4);


        leptonOneP4     = electronOneP4;
        leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonOnePtCorr = electrons[0]->pfPt;
        leptonOneFlavor = 11*electrons[0]->q;
        leptonOneDZ     = electrons[0]->dz;
        leptonOneD0     = electrons[0]->d0;

        leptonTwoP4     = electronTwoP4;
        leptonTwoIso    = GetElectronIsolation(electrons[1], fInfo->rhoJet);
        leptonTwoPtCorr = electrons[1]->pfPt;
        leptonTwoFlavor = 11*electrons[1]->q;
        leptonTwoDZ     = electrons[1]->dz;
        leptonTwoD0     = electrons[1]->d0;

        if (!isData) {
            pair<int, int> pdgId;           
            pdgId           = GetGenId(genParticles, genMotherId, electronOneP4);
            leptonOneGenId  = pdgId.first;
            leptonOneMother = pdgId.second;
            pdgId           = GetGenId(genParticles, genMotherId, electronTwoP4);
            leptonTwoGenId  = pdgId.first;
            leptonTwoMother = pdgId.second;

            EfficiencyContainer effCont;
            pair<float, float> effs, errs;

            // reco scale factors
            effCont = weights->GetElectronRecoEff(electronOneP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            effCont = weights->GetElectronRecoEff(electronTwoP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonTwoRecoWeight = effs.first/effs.second;
            leptonTwoRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            // id/iso scale factors
            effCont = weights->GetElectronIDEff(electronOneP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonOneIDWeight = effs.first/effs.second;
            leptonOneIDVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            effCont = weights->GetElectronIDEff(electronTwoP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonTwoIDWeight = effs.first/effs.second;
            leptonTwoIDVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            // trigger weights with trigger matching
            bitset<2> triggered;
            for (const auto& name: passTriggerNames) {
                if (trigger->passObj(name, 1, electrons[0]->hltMatchBits) && electronOneP4.Pt() > 30)
                    triggered.set(0);
                if (trigger->passObj(name, 1, electrons[1]->hltMatchBits) && electronTwoP4.Pt() > 30)
                    triggered.set(1);
            }

            if (triggered.all()) {
                EfficiencyContainer effCont1, effCont2;
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronOneP4);
                effCont2      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronTwoP4);
                triggerWeight = GetTriggerSF(effCont1, effCont2);
                triggerVar    = GetTriggerSFError(effCont1, effCont2);
            } else if (triggered.test(0)) {
                effCont       = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronOneP4);
                effs          = effCont.GetEff();
                errs          = effCont.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else if (triggered.test(1)) {
                effCont       = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronTwoP4);
                effs          = effCont.GetEff();
                errs          = effCont.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }

            // update event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneIDWeight*leptonTwoIDWeight*leptonOneRecoWeight*leptonTwoRecoWeight;
        }
    } else if (muons.size() == 1 && electrons.size() == 1 && taus.size() == 0) { // e+mu selection
        channel = "emu";
        eventCounts[channel]->Fill(1);

        // trigger matching for thresholds
        float muPtThreshold = 10.;
        float elPtThreshold = 10.;
        if (
                find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_IsoMu24_v*") != passTriggerNames.end() 
                || find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_IsoTkMu24_v*") != passTriggerNames.end()
           ) {
            muPtThreshold = 25.;
        } else if (find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_Ele27_WPTight_Gsf_v*") != passTriggerNames.end()) {
            elPtThreshold = 30.;
        }

        //for (const auto& name: passTriggerNames) {
        //    cout << name << ", ";
        //}
        //cout << muPtThreshold << ", " << elPtThreshold << endl;

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

        // save PDF MC replica variation by jet multiplicity (do this before jet selection)
        if (!isData && fLHEWeightArr != 0) {
            FillPDFHist(pdfVariations, channel);
        }

        //if (nJetsCut < 2 || nBJetsCut < 1)
        //    return kTRUE;
        eventCounts[channel]->Fill(4);

        leptonOneP4     = muonP4;
        leptonOneIso    = GetMuonIsolation(muons[0]);
        leptonOnePtCorr = muons[0]->ptErr;
        leptonOneFlavor = 13*muons[0]->q;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;

        leptonTwoP4     = electronP4;
        leptonTwoIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonTwoPtCorr = electrons[0]->pfPt;
        leptonTwoFlavor = 11*electrons[0]->q;
        leptonTwoDZ     = electrons[0]->dz;
        leptonTwoD0     = electrons[0]->d0;

        if (!isData) {
            pair<int, int> pdgId;           
            pdgId           = GetGenId(genParticles, genMotherId, muonP4);
            leptonOneGenId  = pdgId.first;
            leptonOneMother = pdgId.second;
            pdgId           = GetGenId(genParticles, genMotherId, electronP4);
            leptonTwoGenId  = pdgId.first;
            leptonTwoMother = pdgId.second;

            // reconstruction weights
            EfficiencyContainer effCont;
            pair<float, float> effs, errs;

            effCont = weights->GetMuonISOEff(muonP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            effCont = weights->GetElectronRecoEff(electronP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonTwoRecoWeight = effs.first/effs.second;
            leptonTwoRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            // id weights
            effCont = weights->GetMuonIDEff(muonP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonOneIDWeight = effs.first/effs.second;
            leptonOneIDVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            effCont = weights->GetElectronIDEff(electronP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonTwoIDWeight = effs.first/effs.second;
            leptonTwoIDVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            // trigger weights with trigger matching
            bitset<2> triggered;
            for (const auto& name: passTriggerNames) {
                if (trigger->passObj(name, 1, muons[0]->hltMatchBits))
                    triggered.set(0);
                if (trigger->passObj(name, 1, electrons[0]->hltMatchBits))
                    triggered.set(1);
            }

            if (triggered.all()) {
                EfficiencyContainer effCont1, effCont2;
                effCont1      = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
                effCont2      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronP4);
                triggerWeight = GetTriggerSF(effCont1, effCont2);
                triggerVar    = GetTriggerSFError(effCont1, effCont2);
            } else if (triggered.test(0)) {
                effCont       = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
                effs          = effCont.GetEff();
                errs          = effCont.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else if (triggered.test(1)) {
                effCont       = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronP4);
                effs          = effCont.GetEff();
                errs          = effCont.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }

            // update event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneIDWeight*leptonTwoIDWeight*leptonOneRecoWeight*leptonTwoRecoWeight;
        }
    } else if (electrons.size() == 1 && taus.size() == 1 && muons.size() == 0) { // e+tau selection
        channel = "etau";
        eventCounts[channel]->Fill(1);

        if (electrons[0]->pt < 30 || !electronTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        TLorentzVector electronP4, tauP4, dilepton;
        electronP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        tauP4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, taus[0]->m);
        dilepton = electronP4 + tauP4;
        if (dilepton.M() < 12)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        // save PDF MC replica variation by jet multiplicity (do this before jet selection)
        if (!isData && fLHEWeightArr != 0) {
            FillPDFHist(pdfVariations, channel);
        }

        //if (nJetsCut < 2 || nBJetsCut < 1)
        //    return kTRUE;
        //eventCounts[channel]->Fill(4);

        leptonOneP4     = electronP4;
        leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonOnePtCorr = electrons[0]->pfPt;
        leptonOneFlavor = 11*electrons[0]->q;
        leptonOneDZ     = electrons[0]->dz;
        leptonOneD0     = electrons[0]->d0;

        leptonTwoP4     = tauP4;
        leptonTwoIso    = 0.;
        leptonTwoPtCorr = 0.;
        leptonTwoFlavor = 15*taus[0]->q;
        leptonTwoDZ     = taus[0]->dzLeadChHad;
        leptonTwoD0     = taus[0]->d0LeadChHad;

        //tauChHadMult  = taus[0]->nSignalChHad;
        //tauPhotonMult = taus[0]->nSignalGamma;
        tauDecayMode  = taus[0]->decaymode;
        tauMVA        = taus[0]->rawIsoMVA3newDMwLT;

        if (!isData) {
            pair<int, int> pdgId;           
            pdgId           = GetGenId(genParticles, genMotherId, electronP4);
            leptonOneGenId  = pdgId.first;
            leptonOneMother = pdgId.second;
            pdgId           = GetGenId(genParticles, genMotherId, tauP4);
            leptonTwoGenId  = pdgId.first;
            leptonTwoMother = pdgId.second;

            //cout << "etau: " << leptonOneMother << " " << leptonTwoMother << endl;

            // reconstruction weights
            EfficiencyContainer effCont;
            pair<float, float> effs, errs;

            effCont = weights->GetElectronRecoEff(electronP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            leptonTwoRecoWeight = 1.;
            leptonTwoRecoVar    = 0.;

            // id/iso weights
            effCont = weights->GetElectronIDEff(electronP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonOneIDWeight = effs.first/effs.second;
            leptonOneIDVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            leptonTwoIDWeight = 0.95;
            leptonTwoIDVar    = 0.05;

            // trigger weights
            bool triggered = false;
            for (const auto& name: passTriggerNames) {
                if (trigger->passObj(name, 1, electrons[0]->hltMatchBits)) {
                    triggered = true;
                    break;
                }
            }

            if (triggered) {
                effCont = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronP4);
                effs = effCont.GetEff();
                errs = effCont.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(triggerWeight, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }

            // update event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneIDWeight*leptonTwoIDWeight*leptonOneRecoWeight*leptonTwoRecoWeight;

        }
    } else if (muons.size() == 1 && electrons.size() == 0 && taus.size() == 1) { // mu+tau selection
        channel = "mutau";
        eventCounts[channel]->Fill(1);

        if (muons[0]->pt < 25 || !muonTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        TLorentzVector muonP4, tauP4, dilepton;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        tauP4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, taus[0]->m);
        dilepton = muonP4 + tauP4;
        if (dilepton.M() < 12)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        // save PDF MC replica variation by jet multiplicity (do this before jet selection)
        if (!isData && fLHEWeightArr != 0) {
            FillPDFHist(pdfVariations, channel);
        }

        //if (nJetsCut < 2 || nBJetsCut < 1)
        //    return kTRUE;
        //eventCounts[channel]->Fill(4);

        leptonOneP4     = muonP4;
        leptonOneIso    = GetMuonIsolation(muons[0]);
        leptonOnePtCorr = muons[0]->ptErr;
        leptonOneFlavor = 13*muons[0]->q;
        leptonOneDZ     = muons[0]->dz;
        leptonOneD0     = muons[0]->d0;

        leptonTwoP4     = tauP4;
        leptonTwoIso    = 0.;
        leptonOnePtCorr = 0;
        leptonTwoFlavor = 15*taus[0]->q;
        leptonTwoDZ     = taus[0]->dzLeadChHad;
        leptonTwoD0     = taus[0]->d0LeadChHad;

        //tauChHadMult  = taus[0]->nSignalChHad;
        //tauPhotonMult = taus[0]->nSignalGamma;
        tauDecayMode  = taus[0]->decaymode;
        tauMVA        = taus[0]->rawIsoMVA3newDMwLT;

        if (!isData) {
            pair<int, int> pdgId;           
            pdgId           = GetGenId(genParticles, genMotherId, muonP4);
            leptonOneGenId  = pdgId.first;
            leptonOneMother = pdgId.second;
            pdgId           = GetGenId(genParticles, genMotherId, tauP4);
            leptonTwoGenId  = pdgId.first;
            leptonTwoMother = pdgId.second;

            //cout << "mutau: " << leptonOneMother << " " << leptonTwoMother << endl;

            // reconstruction weights
            EfficiencyContainer effCont;
            pair<float, float> effs, errs;

            effCont = weights->GetMuonIDEff(muonP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonOneIDWeight = effs.first/effs.second;
            leptonOneIDVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            effCont = weights->GetMuonISOEff(muonP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

            leptonTwoRecoWeight = 1.;
            leptonTwoRecoVar    = 0.;
            leptonTwoIDWeight = 0.95;
            leptonTwoIDVar    = 0.05;

            // trigger weights
            bool triggered = false;
            for (const auto& name: passTriggerNames) {
                if (trigger->passObj(name, 1, muons[0]->hltMatchBits)) {
                    triggered = true;
                    break;
                }
            }

            if (triggered) {
                effCont       = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
                effs          = effCont.GetEff();
                errs          = effCont.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(triggerWeight, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }

            // update event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneIDWeight*leptonTwoIDWeight*leptonOneRecoWeight*leptonTwoRecoWeight;
        }
    } else if (muons.size() == 0 && electrons.size() == 1 && taus.size() == 0) { // e+h selection
        channel = "e4j";
        eventCounts[channel]->Fill(1);

        // convert to TLorentzVectors
        TLorentzVector electronP4;
        electronP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        float electronIso = GetElectronIsolation(electrons[0], fInfo->rhoJet);

        if (electrons[0]->pt < 30. || !electronTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        // save PDF MC replica variation by jet multiplicity (do this before jet selection)
        if (!isData && fLHEWeightArr != 0) {
            FillPDFHist(pdfVariations, channel);
        }

        if (nJetsCut < 4 || nBJetsCut < 1)
            return kTRUE;
        eventCounts[channel]->Fill(3);


        leptonOneP4     = electronP4;
        leptonOneIso    = electronIso;
        leptonOnePtCorr = electrons[0]->pfPt;
        leptonOneFlavor = electrons[0]->q*13;
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

        if (!isData) {
            pair<int, int> pdgId;           
            pdgId           = GetGenId(genParticles, genMotherId, electronP4);
            leptonOneGenId  = pdgId.first;
            leptonOneMother = pdgId.second;
            leptonTwoGenId  = 0;
            leptonTwoMother = 0;

            // reconstruction weights
            EfficiencyContainer effCont;
            pair<float, float> effs, errs;

            effCont = weights->GetElectronRecoEff(electronP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            leptonTwoRecoWeight = 1.;
            leptonTwoRecoVar    = 0.;

            // id/iso weights
            effCont = weights->GetElectronIDEff(electronP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonOneIDWeight = effs.first/effs.second;
            leptonOneIDVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            leptonTwoIDWeight = 1.;
            leptonTwoIDVar    = 0.;

            // trigger weights
            bool triggered = false;
            for (const auto& name: passTriggerNames) {
                if (trigger->passObj(name, 1, electrons[0]->hltMatchBits)) {
                    triggered = true;
                    break;
                }
            }

            if (triggered) {
                effCont = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronP4);
                effs = effCont.GetEff();
                errs = effCont.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(triggerWeight, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }

            // update event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneRecoWeight*leptonOneIDWeight;
        }
    }  else if (muons.size() == 1 && electrons.size() == 0 && taus.size() == 0) { // mu+h selection
        channel = "mu4j";
        eventCounts[channel]->Fill(1);

        // convert to TLorentzVectors
        TLorentzVector muonP4;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.1052);
        float muonIso = GetMuonIsolation(muons[0]);
        if (muons[0]->pt < 25.)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        if (muonIso/muonP4.Pt() > 0.15 )
            return kTRUE;
        eventCounts[channel]->Fill(3);

        // save PDF MC replica variation by jet multiplicity (do this before jet selection)
        if (!isData && fLHEWeightArr != 0) {
            FillPDFHist(pdfVariations, channel);
        }

        if (nJetsCut < 4 || nBJetsCut < 1)
            return kTRUE;
        eventCounts[channel]->Fill(4);

        leptonOneP4     = muonP4;
        leptonOneIso    = muonIso;
        leptonOnePtCorr = muons[0]->ptErr;
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

        if (!isData) {
            pair<int, int> pdgId;           
            pdgId           = GetGenId(genParticles, genMotherId, muonP4);
            leptonOneGenId  = pdgId.first;
            leptonOneMother = pdgId.second;
            leptonTwoGenId  = 0;
            leptonTwoMother = 0;

            // reconstruction weights
            EfficiencyContainer effCont;
            pair<float, float> effs, errs;

            effCont = weights->GetMuonISOEff(muonP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            leptonTwoRecoWeight = 1.;
            leptonTwoRecoVar    = 0.;

            // id/iso weights 
            effCont = weights->GetMuonIDEff(muonP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonOneIDWeight = effs.first/effs.second;
            leptonOneIDVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            leptonTwoIDWeight = 1.;
            leptonTwoIDVar    = 0.;

            // trigger weights
            bool triggered = false;
            for (const auto& name: passTriggerNames) {
                if (trigger->passObj(name, 1, muons[0]->hltMatchBits)) {
                    triggered = true;
                    break;
                }
            }

            if (triggered) {
                effCont = weights->GetTriggerEffWeight("HLT_IsoMu24_v*", muonP4);
                effs = effCont.GetEff();
                errs = effCont.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(triggerWeight, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            } else {
                return kTRUE;
            }

            // update event weight
            eventWeight *= triggerWeight;
            eventWeight *= leptonOneRecoWeight*leptonOneIDWeight;
        }
    } else if (fail_muons.size() >= 1 && muons.size() == 0 && electrons.size() == 0) {

        // remove fake muon candidate from the jet collection
        unsigned ix = 0;
        for (const auto& jet: jets) {
            TLorentzVector jetP4, muonP4; 
            jetP4.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
            muonP4.SetPtEtaPhiM(fail_muons[0]->pt, fail_muons[0]->eta, fail_muons[0]->phi, 0.1052);
            if (jetP4.DeltaR(muonP4) < 0.4) {
                --nJets;
                if (jet->csv > 0.8484) { 
                    --nBJets;
                } 
                jets.erase(jets.begin() + ix);
                break;
            }
            ++ix;  
        }

        if (taus.size() >= 1) {
            channel = "mutau_fakes";
            eventCounts[channel]->Fill(1);
            nMuons = 1;

            if (fail_muons[0]->pt < 25)
                return kTRUE;
            eventCounts[channel]->Fill(2);

            TLorentzVector muonP4, tauP4, dilepton;
            muonP4.SetPtEtaPhiM(fail_muons[0]->pt, fail_muons[0]->eta, fail_muons[0]->phi, 0.1052);
            tauP4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, taus[0]->m);
            dilepton = muonP4 + tauP4;
            if (dilepton.M() < 12)
                return kTRUE;
            eventCounts[channel]->Fill(3);

            //if (nJets < 2 || nBJets < 1)
            //    return kTRUE;
            //eventCounts[channel]->Fill(4);

            leptonOneP4     = muonP4;
            leptonOneIso    = GetMuonIsolation(fail_muons[0]);
            leptonOneFlavor = 13*fail_muons[0]->q;
            leptonOneDZ     = fail_muons[0]->dz;
            leptonOneD0     = fail_muons[0]->d0;

            leptonTwoP4     = tauP4;
            leptonTwoIso    = 0.;
            leptonTwoFlavor = 15*taus[0]->q;
            leptonTwoDZ     = taus[0]->dzLeadChHad;
            leptonTwoD0     = taus[0]->d0LeadChHad;

        } else if (isData && nJets >= 4 && nBJets >= 1) {
            channel = "mu4j_fakes";
            eventCounts[channel]->Fill(1);
            nMuons = 1;

            // convert to TLorentzVectors
            TLorentzVector muonP4;
            muonP4.SetPtEtaPhiM(fail_muons[0]->pt, fail_muons[0]->eta, fail_muons[0]->phi, 0.1052);
            float muonIso = GetMuonIsolation(fail_muons[0]);
            if (fail_muons[0]->pt < 25.)
                return kTRUE;
            eventCounts[channel]->Fill(2);

            leptonOneP4     = muonP4;
            leptonOneIso    = muonIso;
            leptonOneFlavor = fail_muons[0]->q*13;
            leptonOneDZ     = fail_muons[0]->dz;
            leptonOneD0     = fail_muons[0]->d0;

            // Collect the highest pt jets in the event
            std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);
            //jets = KinematicTopTag(jets, metP2, muonP4);
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
        } else {
            return kTRUE;
        }
    } else if (fail_electrons.size() >= 1 && muons.size() == 0 && electrons.size() == 0) {

        // remove fake electron candidate from the jet collection
        unsigned ix = 0;
        for (const auto& jet: jets) {
            TLorentzVector jetP4, electronP4; 
            jetP4.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
            electronP4.SetPtEtaPhiM(fail_electrons[0]->pt, fail_electrons[0]->eta, fail_electrons[0]->phi, 0.1052);
            if (jetP4.DeltaR(electronP4) < 0.4) {
                --nJets;
                if (jet->csv > 0.8484) { 
                    --nBJets;
                } 
                jets.erase(jets.begin() + ix);
                break;
            }
            ++ix;  
        }

        if (taus.size() >= 1) {
            channel = "etau_fakes";
            eventCounts[channel]->Fill(1);
            nElectrons = 1;

            if (fail_electrons[0]->pt < 25)
                return kTRUE;
            eventCounts[channel]->Fill(2);

            TLorentzVector electronP4, tauP4, dilepton;
            electronP4.SetPtEtaPhiM(fail_electrons[0]->pt, fail_electrons[0]->eta, fail_electrons[0]->phi, 0.1052);
            tauP4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, taus[0]->m);
            dilepton = electronP4 + tauP4;
            if (dilepton.M() < 12)
                return kTRUE;
            eventCounts[channel]->Fill(3);

            //if (nJets < 2 || nBJets < 1)
            //    return kTRUE;
            //eventCounts[channel]->Fill(4);

            leptonOneP4     = electronP4;
            leptonOneIso    = GetElectronIsolation(fail_electrons[0], fInfo->rhoJet);
            leptonOneFlavor = 13*fail_electrons[0]->q;
            leptonOneDZ     = fail_electrons[0]->dz;
            leptonOneD0     = fail_electrons[0]->d0;

            leptonTwoP4     = tauP4;
            leptonTwoIso    = 0.;
            leptonTwoFlavor = 15*taus[0]->q;
            leptonTwoDZ     = taus[0]->dzLeadChHad;
            leptonTwoD0     = taus[0]->d0LeadChHad;

        } else if (isData && nJets >= 4 && nBJets >= 1) {
            channel = "e4j_fakes";
            eventCounts[channel]->Fill(1);
            nElectrons = 1;

            // convert to TLorentzVectors
            TLorentzVector electronP4;
            electronP4.SetPtEtaPhiM(fail_electrons[0]->pt, fail_electrons[0]->eta, fail_electrons[0]->phi, 0.1052);
            float electronIso = GetElectronIsolation(fail_electrons[0], fInfo->rhoJet);
            if (fail_electrons[0]->pt < 25.)
                return kTRUE;
            eventCounts[channel]->Fill(2);

            leptonOneP4     = electronP4;
            leptonOneIso    = electronIso;
            leptonOneFlavor = fail_electrons[0]->q*13;
            leptonOneDZ     = fail_electrons[0]->dz;
            leptonOneD0     = fail_electrons[0]->d0;

            // Collect the highest pt jets in the event
            std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);
            //jets = KinematicTopTag(jets, metP2, electronP4);
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
        } else {
            return kTRUE;
        }
    } else {
        return kTRUE;
    }

    ///////////////////
    // Fill jet info //
    ///////////////////

    if (channel != " tau4j" && channel != "mu4j" && channel != "e4j" && channel != "mu4j_fakes" && channel != "e4j_fakes") { // jets are handled differently for the lepton + jet selection
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

    if (!isData) {
        if (genParticles.size() >= 1) {
            genOneId = genParticles[0]->pdgId;
            genOneP4.SetPtEtaPhiM(genParticles[0]->pt, genParticles[0]->eta, genParticles[0]->phi, genParticles[0]->mass); 
        } else if (genParticles.size() >= 2) {
            genTwoId = genParticles[1]->pdgId;
            genTwoP4.SetPtEtaPhiM(genParticles[1]->pt, genParticles[1]->eta, genParticles[1]->phi, genParticles[1]->mass); 
        }

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

pair<int, int> MultileptonAnalyzer::GetGenId(vector<TGenParticle*> particles, vector<int> motherIDs, TLorentzVector p4)
{
    int motherId = 0;
    int pid = 0;
    for (unsigned i = 0; i < particles.size(); ++i) {
        TLorentzVector genP4;
        genP4.SetPtEtaPhiM(particles[i]->pt, particles[i]->eta, particles[i]->phi, particles[i]->mass); 
        if (genP4.DeltaR(p4) < 0.3) {
            pid = particles[i]->pdgId;
            motherId = motherIDs[i];
        }
    }
    return std::make_pair(pid, motherId);
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
    if (trigEff1.second > 0 || trigEff2.second > 0) {
        sf = (1 - (1 - trigEff1.first)*(1 - trigEff2.first))/(1 - (1 - trigEff1.second)*(1 - trigEff2.second));
    }

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
    nJetsJERUp     = nJetsJERDown     = 0;
    nBJetsJERUp    = nBJetsJERDown    = 0;
    nBJetsCTagUp   = nBJetsCTagDown   = 0;
    nBJetsMistagUp = nBJetsMistagDown = 0;

    std::fill(nJetsJESUp.begin(), nJetsJESUp.end(), 0);
    std::fill(nJetsJESDown.begin(), nJetsJESDown.end(), 0);
    std::fill(nBJetsJESUp.begin(), nBJetsJESUp.end(), 0);
    std::fill(nBJetsJESDown.begin(), nBJetsJESDown.end(), 0);
    std::fill(nBJetsBTagUp.begin(), nBJetsBTagUp.end(), 0);
    std::fill(nBJetsBTagDown.begin(), nBJetsBTagDown.end(), 0);
}

void MultileptonAnalyzer::JetCounting(TJet* jet, float jerc_nominal, float resRand)
{
    float jetPt = jet->pt;
    //cout << jetPt << " ";

    float rNumber = rng->Uniform(1.);
    if (jet->pt > 30) {

        // nominal
        ++nJets;
        if (particleSelector->BTagModifier(jet, "CSVM", "central", rNumber)) {
            ++nBJets;

            if (jet->hadronFlavor == 5) {
                ++nBJetsCTagUp;
                ++nBJetsCTagDown;
                ++nBJetsMistagUp;
                ++nBJetsMistagDown;
            } else if (jet->hadronFlavor == 4) {
                ++nBJetsMistagUp;
                ++nBJetsMistagDown;

                unsigned bcount = 0;
                for (const auto& name: particleSelector->GetBTagSourceNames()) {
                    ++nBJetsBTagUp[bcount];
                    ++nBJetsBTagDown[bcount];
                    bcount++;
                }
            } else {
                ++nBJetsCTagUp;
                ++nBJetsCTagDown;

                unsigned bcount = 0;
                for (const auto& name: particleSelector->GetBTagSourceNames()) {
                    ++nBJetsBTagUp[bcount];
                    ++nBJetsBTagDown[bcount];
                    bcount++;
                }
            }
        }

        if (jet->hadronFlavor == 5) {
            unsigned bcount = 0;
            for (const auto& name: particleSelector->GetBTagSourceNames()) {
                // b tag up
                if (particleSelector->BTagModifier(jet, "CSVM", "up_" + name, rNumber)) 
                    ++nBJetsBTagUp[bcount];
                     
                // b tag down
                if (particleSelector->BTagModifier(jet, "CSVM", "down_" + name, rNumber)) 
                    ++nBJetsBTagDown[bcount];
                bcount++;
            }

        } else if (jet->hadronFlavor == 4){
            // c tag up
            if (particleSelector->BTagModifier(jet, "CSVM", "up", rNumber)) 
                ++nBJetsCTagUp;
                 
            // c tag down
            if (particleSelector->BTagModifier(jet, "CSVM", "down", rNumber)) 
                ++nBJetsCTagDown;

        } else { 
            // mistag up
            if (particleSelector->BTagModifier(jet, "CSVM", "up", rNumber)) 
                ++nBJetsMistagUp;
                 
            // mistag down
            if (particleSelector->BTagModifier(jet, "CSVM", "down", rNumber)) 
                ++nBJetsMistagDown;

        }
    }

    double jec = particleSelector->JetCorrector(jet, "NONE");
    unsigned count = 0;
    for (const auto& name: particleSelector->GetJECSourceNames()) {

        // JES up
        float jecUnc = particleSelector->JetUncertainty(jet, name);
        jet->pt = jet->ptRaw*jec*(1 + jecUnc)*jerc_nominal;
        if (jet->pt > 30) {
            ++nJetsJESUp[count];
            if (particleSelector->BTagModifier(jet, "CSVM", "central", rNumber)) { 
                ++nBJetsJESUp[count];
            } 
        }

        // JES down
        jet->pt = jet->ptRaw*jec*(1 - jecUnc)*jerc_nominal;
        //cout << jet->pt << endl;
        if (jet->pt > 30) {
            ++nJetsJESDown[count];
            if (particleSelector->BTagModifier(jet, "CSVM", "central", rNumber)) { 
                ++nBJetsJESDown[count];
            }
        }
        //cout << name << ": " << jecUnc << endl;
        count++;
    }

    //cout << jet->ptRaw << " " << jerc_nominal << " " << jec << " " << jecUnc << endl;

    // JER up
    pair<float, float> resPair = particleSelector->JetResolutionAndSF(jet, 1);
    float jerc = 1 + resRand*sqrt(std::max((double)resPair.second*resPair.second - 1, 0.));
    jet->pt = jet->ptRaw*jec*jerc;
    if (jet->pt > 30) {
        ++nJetsJERUp;
        if (particleSelector->BTagModifier(jet, "CSVM", "central", rNumber)) { 
            ++nBJetsJERUp;
        }     
    }

    // JER down
    resPair     = particleSelector->JetResolutionAndSF(jet, -1);
    jerc        = 1 + resRand*sqrt(std::max((double)resPair.second*resPair.second - 1, 0.));
    jet->pt     = jet->ptRaw*jec*jerc;
    if (jet->pt > 30) {
        ++nJetsJERDown;
        if (particleSelector->BTagModifier(jet, "CSVM", "central", rNumber)) { 
            ++nBJetsJERDown;
        } 
    }
    jet->pt = jetPt;

}

void MultileptonAnalyzer::FillPDFHist(vector<float> pdfVariations, string channel)
{
    // N.B. nPartons and nJets are both defined at global scope
    for (unsigned i = 0; i < qcdWeights.size(); ++i) {
        if (nJets <= 4) {
            qcdCountsJets[channel]->Fill(i+1, nJets, qcdWeights[i]);
        } else {
            qcdCountsJets[channel]->Fill(i+1, 4, qcdWeights[i]);
        }
        qcdCountsPartons[channel]->Fill(i+1, nPartons, qcdWeights[i]);
    }

    for (unsigned i = 0; i < pdfVariations.size(); ++i) {
        if (nJets <= 4) {
            pdfCountsJets[channel]->Fill(i+1, nJets, pdfVariations[i]);
        } else {
            pdfCountsJets[channel]->Fill(i+1, 4, pdfVariations[i]);
        }
        pdfCountsPartons[channel]->Fill(i+1, nPartons, pdfVariations[i]);
    }
}
