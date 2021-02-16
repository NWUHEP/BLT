#include "LLTauAnalyzer.h"
#include <map>


//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;



LLTauAnalyzer::LLTauAnalyzer() : BLTSelector()
{

}

LLTauAnalyzer::~LLTauAnalyzer()
{

}

void LLTauAnalyzer::Begin(TTree *tree)
{   
    cout<< "begin"<<endl;
    rng = new TRandom3();


    // Parse command line option
    cout<< "Parse command line option"<<endl;
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1), std::sregex_token_iterator(), std::back_inserter(options));

    // Set the parameters
    cout<< "Set the parameters"<<endl;
    params.reset(new Parameters());
    params->setup(options);

    // Set the cuts
    cout<< "Set the cuts"<<endl;
    cuts.reset(new Cuts());
    particleSelector.reset(new ParticleSelector(*params, *cuts));

    // Trigger bits mapping file
    cout<< "Trigger bits mapping file"<<endl;
    const std::string cmssw_base = getenv("CMSSW_BASE");
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns";
    trigger.reset(new baconhep::TTrigger(trigfilename));


    singleElectronTriggerName = "HLT_Ele27_WPTight_Gsf_v*";
    triggerNames.push_back(singleElectronTriggerName);
    triggerNames.push_back("HLT_IsoMu24_v*");
    triggerNames.push_back("HLT_IsoTkMu24_v*");
    

    // Weight utility class
    cout<< "Weight utility class"<<endl;
    cout<< "-- period " << params->period << ", selection " << params->selection <<endl;
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    string jsonFileName = "";
    jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/json/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"; // 2016 mask
    
   

    cout<< "-- use lumi mask: " << jsonFileName <<endl;
    lumiMask.AddJSONFile(jsonFileName);
    


    // muon momentum corrections
    cout<< "muon momentum corrections"<<endl;
    string muonCorrFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_es/rcdata.2016.v3";
    cout<< "-- use muon Roc Correction: " << muonCorrFileName <<endl;
    muonCorr = new RoccoR(muonCorrFileName);
  
    // electron scale corrections
    cout<< "electron scale corrections"<<endl;
    electronScaler = new EnergyScaleCorrection(cmssw_base + "/src/BLT/BLTAnalysis/data/electron_es");

    // Prepare the output tree
    cout<< "Prepare the output tree"<<endl;
    string outFileName = params->get_output_filename("output");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();

    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);
    outHistName = params->get_output_treename("GenCategory");
    hGenCat = new TH1D(outHistName.c_str(), "WW decay modes",30,0.5,30.5);

    vector<std::string> channelNames = {"mumutau", "eetau", "emutau"};

    for (unsigned i = 0; i < channelNames.size(); ++i) {
        string channel = channelNames[i];
        outFile->mkdir(channel.c_str());
        outFile->cd(channel.c_str());
        string treeName = "bltTree_" + params->datasetgroup;
        tree = new TTree(treeName.c_str(), treeName.c_str());

        // event data
        tree->Branch("runNumber", &runNumber);
        tree->Branch("evtNumber", &evtNumber, "eventNumber/l");
        tree->Branch("lumiSection", &lumiSection);
        tree->Branch("nPV", &nPV);
        tree->Branch("nPU", &nPU);
        tree->Branch("triggerLeptonStatus",&triggerLeptonStatus);
        tree->Branch("mcEra",&mcEra);

        tree->Branch("genTauOneDaughters", &genTauOneDaughters);
        tree->Branch("genTauTwoDaughters", &genTauTwoDaughters);


        // gen level objects
        tree->Branch("nPartons", &nPartons);
        tree->Branch("genWeight", &genWeight);
        tree->Branch("genCategory", &genCategory);
        
        // weights and their uncertainties
        tree->Branch("eventWeight", &eventWeight);
        tree->Branch("topPtWeight", &topPtWeight);
        tree->Branch("puWeight", &puWeight);
        tree->Branch("triggerWeight", &triggerWeight);
        tree->Branch("topPtVar", &topPtVar);
        if ( channel == "eetau" || channel == "emutau"){
            tree->Branch("eleTriggerVarTagSyst", &eleTriggerVarTagSyst);
            tree->Branch("eleTriggerVarProbeSyst", &eleTriggerVarProbeSyst);
            tree->Branch("prefiringWeight", &prefiringWeight);
            tree->Branch("prefiringVar", &prefiringVar);
        }

        tree->Branch("leptonOneIDWeight", &leptonOneIDWeight);
        tree->Branch("leptonTwoIDWeight", &leptonTwoIDWeight);
        tree->Branch("leptonOneRecoWeight", &leptonOneRecoWeight);
        tree->Branch("leptonTwoRecoWeight", &leptonTwoRecoWeight);
        
        tree->Branch("leptonOneIDVar", &leptonOneIDVar);
        tree->Branch("leptonTwoIDVar", &leptonTwoIDVar);
        tree->Branch("leptonOneRecoVar", &leptonOneRecoVar);
        tree->Branch("leptonTwoRecoVar", &leptonTwoRecoVar);


        // leptons
        tree->Branch("leptonOneP4" , &leptonOneP4);
        tree->Branch("leptonOneIso", &leptonOneIso);
        tree->Branch("leptonOneFlavor", &leptonOneFlavor);
        tree->Branch("leptonOneD0" , &leptonOneD0);
        tree->Branch("leptonOneDz" , &leptonOneDz);
        tree->Branch("leptonOneSip3d" , &leptonOneSip3d);

        tree->Branch("leptonTwoP4" , &leptonTwoP4);
        tree->Branch("leptonTwoIso", &leptonTwoIso);
        tree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
        tree->Branch("leptonTwoD0" , &leptonTwoD0);
        tree->Branch("leptonTwoDz" , &leptonTwoDz);
        tree->Branch("leptonTwoSip3d" , &leptonTwoSip3d);


        // taus
        tree->Branch("leptonThreeP4",    &leptonThreeP4);
        tree->Branch("leptonThreeIso",   &leptonThreeIso); 
        tree->Branch("tauDecayMode",     &tauDecayMode);
        tree->Branch("tauGenFlavor",     &tauGenFlavor);
        tree->Branch("tauGenFlavorHad",  &tauGenFlavorHad);
        tree->Branch("tauVetoedJetPt",   &tauVetoedJetPt);
        tree->Branch("tauVetoedJetPtUnc",&tauVetoedJetPtUnc);


        // ht
        tree->Branch("htSum", &htSum);
        tree->Branch("ht", &ht);
        tree->Branch("htPhi", &htPhi);

        // met
        tree->Branch("met", &met);
        tree->Branch("metPhi", &metPhi);

        // object counters
        tree->Branch("nMuons", &nMuons);
        tree->Branch("nElectrons", &nElectrons);
        tree->Branch("nTaus", &nTaus);
        tree->Branch("nJets", &nJets);
        tree->Branch("nBJets", &nBJets);

        // jet systematics
        tree->Branch("nJetsJESUp",      &nJetsJESUp);
        tree->Branch("nJetsJESDown",    &nJetsJESDown);
        tree->Branch("nJetsJERUp",      &nJetsJERUp);
        tree->Branch("nJetsJERDown",    &nJetsJERDown);
        tree->Branch("nBJetsJESUp",     &nBJetsJESUp);
        tree->Branch("nBJetsJESDown",   &nBJetsJESDown);
        tree->Branch("nBJetsJERUp",     &nBJetsJERUp);
        tree->Branch("nBJetsJERDown",   &nBJetsJERDown);
        tree->Branch("nBJetsBTagUp",    &nBJetsBTagUp);
        tree->Branch("nBJetsBTagDown",  &nBJetsBTagDown);
        tree->Branch("nBJetsMistagUp",  &nBJetsMistagUp);
        tree->Branch("nBJetsMistagDown",&nBJetsMistagDown);



        outTrees[channel] = tree;

        // event counter
        string outHistName = params->get_output_treename("TotalEvents_" + channel);
        eventCounts[channel] = new TH1D(outHistName.c_str(),"ChannelCounts",10,0.5,10.5);
    }
    cout<< "Done with initialization"<<endl;
    ReportPostBegin();

}

Bool_t LLTauAnalyzer::Process(Long64_t entry)
{
    /* beginning of event loop */
    GetEntry(entry, 1);  // load all branches
    outFile->cd();
    eventWeight = 1.;
    this->totalEvents++;
    hTotalEvents->Fill(1);
    
    if (entry%10000==0) {
        std::cout << "... Processing event " << entry 
            << " Run: " << fInfo->runNum 
            << " Lumi: " << fInfo->lumiSec 
            << " Event: " << fInfo->evtNum 
            << std::endl;
    }


    const bool isData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isData);

    /* Apply lumi mask */
    if (isData) {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(2);

    // Set data period for 2016 MC scale factors
    mcEra = -1;
    if (!isData) {
        if (rng->Rndm() < 0.548) {  //0.468) {
            weights->SetDataPeriod("2016BtoF"); 
            mcEra = 0;
        } else {
            weights->SetDataPeriod("2016GH");
            mcEra = 1;
        }
    }

    /////////////////////
    // Fill event info //
    /////////////////////

    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
    nPV           = fPVArr->GetEntries();
    
    
    ///////////////////////
    // Generator objects //
    ///////////////////////

    vector<TGenParticle*> genTausHad, genElectrons, genMuons;

    if (!isData) {


        // save gen weight for amc@nlo Drell-Yan sample
        genWeight = fGenEvtInfo->weight > 0 ? 1 : -1;
        if (genWeight < 0) {
            hTotalEvents->Fill(10);
        }


        // pileup reweighting
        nPU          = fInfo->nPUmean;
        puWeight     = weights->GetPUWeight(fInfo->nPUmean); 
        eventWeight *= puWeight;


        // loop over gen particle collection:

        // 1.initialization
        //   1.a) counting parton for Drell-Yan samples
        unsigned partonCount= 0;

        //   1.b) counting top for top reweighting
        float topSF = 1.;
        unsigned topCount   = 0;
        topPtWeight = 1.;
        topPtVar = 0.;

        float zPt = -1.;
        TLorentzVector zP4(0., 0., 0., 0);
        zPtWeight = 1.;
        zPtVar    = 0.;

        //   1.c) counting W,tau for WW categorization
        unsigned wCount     = 0;
        unsigned tauCount   = 0;
        bitset<6> wDecay;
        bitset<4> tauDecay;


        //   1.d) save gen level tau and leptonic taus, and e.mu
        vector<TGenParticle*> genTaus, genTausLep, genTausDaughters;
        
        
        // 2. loop over gen particles
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
            
            // // [test] print out gen particles
            // if  (particle->parent >=0 ){
            //     TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
            //     if (abs(mother->pdgId) != abs(particle->pdgId) ) std::cout<< " (" << abs(mother->pdgId) << ")" << abs(particle->pdgId); 
            // }
            

            // 2.a) parton counting: update partonCount
            if (
                    particle->status == 23 
                    && (abs(particle->pdgId) < 6 || particle->pdgId == 21) 
                    && particle->parent >=0
               ) {
                ++partonCount;
            }

            // 2.b) top pt reweighting: update topSF,topCount. get the scale factor based on the top quark pt
            if (abs(particle->pdgId) == 6 && particle->status == 62) {
                topSF *= exp(0.0615 - 0.0005*particle->pt);
                ++topCount;
            }
            // identify Z boson
            if (abs(particle->pdgId) == 23 && particle->status == 62) {
                zPt = particle->pt;
            }
            // Z not included in low mass DY sample :(
            if ((abs(particle->pdgId) == 11 || abs(particle->pdgId) == 13 || abs(particle->pdgId) == 15) && particle->status == 23 && particle->parent == 0) {
                TLorentzVector ellP4(0., 0., 0., 0.);
                ellP4.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, particle->mass);
                zP4 += ellP4;
            }


            // 2.c) Tag WW decay based on neutrino flavors: update wCount,tauCount,wDecay,tauDecay
            unsigned flavor = abs(particle->pdgId);
            if (
                    (flavor == 12 || flavor == 14 || flavor == 16)
                    && particle->parent >=0
               ) {
                
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);
                
                if (abs(mother->pdgId) == 24) {
                    wDecay.set(3*wCount + (flavor - 12)/2);
                    ++wCount;
                } 

                if (abs(mother->pdgId) == 15 && flavor != 16 && mother->parent >=0) {
                    TGenParticle* gmother = (TGenParticle*) fGenParticleArr->At(mother->parent); // gmother should not from initial quark
                    
                    bool isHard = gmother->status != 2 && gmother->parent != -2;
                    if (isHard && tauCount<2){
                        tauDecay.set(2*tauCount + (flavor - 12)/2);
                        ++tauCount;
                        // cout<<abs(gmother->pdgId)<<"-"<<abs(mother->pdgId)<<"-"<<abs(particle->pdgId)<<endl;
                    }                      
                      
                } 
            }


            // 2.d) This saves gen level tau, leptonic taus, electrons and muons
            if ( particle->parent >=0 ) {
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);

                if (abs(particle->pdgId) != 22 && abs(particle->pdgId) != 15 && abs(mother->pdgId) == 15 && mother->parent >=0) {
                    TGenParticle* gmother = (TGenParticle*) fGenParticleArr->At(mother->parent);
                    bool isHard = gmother->status != 2 && gmother->parent != -2;

                    if (isHard && abs(particle->pdgId) == 16) genTaus.push_back(mother); 
                    if (isHard && abs(particle->pdgId) != 16) genTausDaughters.push_back(particle);
                }

                if ( abs(particle->pdgId) == 11 || abs(particle->pdgId) == 13 ) {
                    if (abs(mother->pdgId) == 15) {
                        genTausLep.push_back(mother);
                        if (abs(particle->pdgId) == 11) genElectrons.push_back(particle);
                        if (abs(particle->pdgId) == 13) genMuons.push_back(particle);
                    } else if (abs(mother->pdgId) == 23 || abs(mother->pdgId) == 24) {
                        if (abs(particle->pdgId) == 11) genElectrons.push_back(particle);
                        if (abs(particle->pdgId) == 13) genMuons.push_back(particle);
                    }
                }
                
                
            }
        }
        // cout << endl;
        // cout << "nGenElectron=" << genElectrons.size() << ", nGenMuon= " << genMuons.size() << ", nGenTau= " << genTaus.size() << endl;
        // 3. finish loop over gen particles. now conclude gen level infomation

        // 3.a) counting partons: save to nPartons
        nPartons = partonCount; 

        // 3.b) Account for the top pt weights
        if (params->datasetgroup.substr(0, 5) == "ttbar") {
            topPtWeight *= sqrt(topSF);
            topPtVar    += pow(0.01*topPtWeight, 2);
            eventWeight *= topPtWeight;
        }

        // 3.b) Z pt weights
        if (
                (params->datasetgroup.substr(0, 5) == "zjets" 
                || params->datasetgroup.substr(0, 6) == "z0jets" 
                || params->datasetgroup.substr(0, 6) == "z1jets" 
                || params->datasetgroup.substr(0, 6) == "z2jets") 
                && (zPt > 0. || zP4.Pt() > 0.)) {
            if (zP4.Pt() > 0) 
                zPt = zP4.Pt();

            if (zPt < 140) {
                zPtWeight = (0.876979 + 4.11598e-3*zPt - 2.3552e-5*zPt*zPt);
                zPtWeight *= 1.10211*(0.958512 - 0.131835*erf((zPt - 14.1972)/10.1525));
            } else {
                zPtWeight = 0.891188;
            }
            zPtVar    = 0.01*zPtWeight;

            // apply zPt reweighting
            eventWeight *= zPtWeight;
        }

        
        // 3.c) categorize events: save to genCategory
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

        
        // 3.d) genTausHad = genTaus - genTausLep
        for (unsigned i = 0; i < genTaus.size(); ++i) {
            
            bool isHad = true;

            TLorentzVector genTauP4;
            genTauP4.SetPtEtaPhiM(genTaus[i]->pt, genTaus[i]->eta, genTaus[i]->phi, genTaus[i]->mass); 

            for (unsigned j = 0; j < genTausLep.size(); ++j) {
                TLorentzVector genTauLepP4;
                genTauLepP4.SetPtEtaPhiM(genTausLep[j]->pt, genTausLep[j]->eta, genTausLep[j]->phi, genTausLep[j]->mass); 

                if (genTauP4.DeltaR(genTauLepP4) < 0.1) {
                    isHad = false;
                }
            }
            if (isHad) genTausHad.push_back(genTaus[i]);
        }

        // 3.e) get code for tauDecayMode;
        std::string genTau1Daughters("0000000");
        std::string genTau2Daughters("0000000");
        for (unsigned i = 0; i < genTaus.size(); ++i) {
          for (unsigned j = 0; j < genTausDaughters.size(); ++j) {
            TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(genTausDaughters[j]->parent);
            if (genTaus[i]->pt == mother->pt && genTaus[i]->eta== mother->eta && genTaus[i]->phi==mother->phi){
              int pdgid = abs(genTausDaughters[j]->pdgId);
              int daughterDigitPos = -1;

              if (pdgid == 11) { daughterDigitPos = 0; }
              else if (pdgid == 13) {daughterDigitPos = 1;}
              else if (pdgid ==111) {daughterDigitPos = 2;}
              else if (pdgid ==211) {daughterDigitPos = 3;}
              else if (pdgid ==311) {daughterDigitPos = 4;}
              else if (pdgid ==321) {daughterDigitPos = 5;}
              else if (pdgid >=100) {daughterDigitPos = 6;}

              if (i == 0 && daughterDigitPos>=0) genTau1Daughters[daughterDigitPos] = char(int(genTau1Daughters[daughterDigitPos])+1);
              if (i == 1 && daughterDigitPos>=0) genTau2Daughters[daughterDigitPos] = char(int(genTau2Daughters[daughterDigitPos])+1);
            }

          }
        }

        genTauOneDaughters = std::atoi(genTau1Daughters.c_str());
        genTauTwoDaughters = std::atoi(genTau2Daughters.c_str());
        if (genTaus.size()==2) {
          if (genTaus[1]->pt>genTaus[0]->pt){
            genTauTwoDaughters = std::atoi(genTau1Daughters.c_str());
            genTauOneDaughters = std::atoi(genTau2Daughters.c_str());
          }
        }
        
        // cout << "tau Decay modes are [e,mu,pi0,pi+,k0,k+,h]: " <<genTau1Daughters << "," <<  genTau2Daughters << endl;
        // cout << "------------" << endl;
    } else {
        genWeight = 1.0;
        nPU = 0;
        nPartons = 0;
        genCategory = 0;

    }

    // cout << nPU << endl;
    ///////////////////
    // Select objects//
    ///////////////////

    /* -------- Vertices ---------*/

    if (fInfo->hasGoodPV) {
        assert(fPVArr->GetEntries() != 0);
        TVector3 pv;
        copy_xyz((TVertex*) fPVArr->At(0), pv);
        particleSelector->SetPV(pv);
    } else {
        return kTRUE;
    }
    hTotalEvents->Fill(3);
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);



    /* -------- Trigger --------- */
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

    bool muonTriggered = find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_IsoMu24_v*") != passTriggerNames.end();
    muonTriggered = muonTriggered || find(passTriggerNames.begin(), passTriggerNames.end(), "HLT_IsoTkMu24_v*") != passTriggerNames.end();
    bool electronTriggered = find(passTriggerNames.begin(), passTriggerNames.end(), singleElectronTriggerName) != passTriggerNames.end();

    if (!passTrigger)
        return kTRUE;
    hTotalEvents->Fill(4);




    /* -------- MUONS ---------*/
    vector<TMuon*> muons, fail_muons;
    vector<TLorentzVector> veto_muons;
    for (int i=0; i < fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);

        // Apply rochester muon momentum corrections
        double muonSF = 1.;
        if (isData) {
            muonSF = muonCorr->kScaleDT(muon->q, muon->pt, muon->eta, muon->phi, 0, 0);
        } else {
            muonSF = muonCorr->kScaleAndSmearMC(muon->q, muon->pt, muon->eta, muon->phi,
                    muon->nTkLayers, rng->Rndm(), rng->Rndm(), 
                    0, 0);
        }
        muon->pt = muonSF*muon->pt; 
        
        TLorentzVector muonP4;
        muonP4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, 0.10566);

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
                && GetMuonIsolation(muon)/muon->pt < 0.25 // loose muon
           ) {
            

            if (GetMuonIsolation(muon)/muonP4.Pt() <  0.15){
                muons.push_back(muon);
            } else {
                fail_muons.push_back(muon);
            }
            veto_muons.push_back(muonP4);
        }
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);
    sort(fail_muons.begin(), fail_muons.end(), sort_by_higher_pt<TMuon>);
    nMuons = muons.size();
    nFailMuons = fail_muons.size();

    // fill hist if muonTriggered and leading muon is good.
    if (muonTriggered && nMuons>0){
        if (muons[0]->pt>25){
            hTotalEvents->Fill(5);
        }
    }
    


    /* -------- ELECTRONS ---------*/
    vector<TElectron*> electrons,fail_electrons;
    vector<TLorentzVector> veto_electrons;

    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);
        
        // apply electron scale and smear
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
                electron->pt > 20
                && fabs(electron->scEta) < 2.5
                && particleSelector->PassElectronID(electron, cuts->tightElID)
                && particleSelector->PassElectronIso(electron, cuts->looseElIso, cuts->EAEl)

           ) {
            
            if (particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl) ){
                electrons.push_back(electron);
            } else {
                fail_electrons.push_back(electron);
            }
            veto_electrons.push_back(electronP4);
        }
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);
    sort(fail_electrons.begin(), fail_electrons.end(), sort_by_higher_pt<TElectron>);

    nElectrons = electrons.size();
    nFailElectrons = fail_electrons.size();

    // fill hist if muonTriggered and leading muon is good.
    if (electronTriggered && nElectrons>0){
        if (electrons[0]->pt>30){
            hTotalEvents->Fill(6);
        }
    }

    
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

       
        // check apply correction for tau_h->tau_h and e->tau_h
        if (!isData) {
            bool genOverlap = false;

            // 1. MC events with tau_h -> tau_h
            if (genOverlap == false){
                // 1.a) check overlap with hadronic tau
                bool genTauHadOverlap = false;
                for (unsigned i = 0; i < genTausHad.size(); ++i) {
                    TLorentzVector genP4;
                    genP4.SetPtEtaPhiM(genTausHad[i]->pt, genTausHad[i]->eta, genTausHad[i]->phi, genTausHad[i]->mass); 
                    if (tauP4.DeltaR(genP4) < 0.3) {
                        genTauHadOverlap = true;
                        genOverlap = true;
                        break;
                    }
                }
                // 1.b) apply correction for tau_h -> tau_h 
                // (https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV#Tau_energy_scale)
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
            
            // // 2. MC events with electron -> tau_h
            // if (genOverlap == false){
            //     // 2.a) check overlap with  electron
            //     bool genElectronOverlap = false;
            //     for (unsigned i = 0; i < genElectrons.size(); ++i) {
            //         TLorentzVector genP4;
            //         genP4.SetPtEtaPhiM(genElectrons[i]->pt, genElectrons[i]->eta, genElectrons[i]->phi, genElectrons[i]->mass); 
            //         if (tauP4.DeltaR(genP4) < 0.3) {
            //             genElectronOverlap = true;
            //             genOverlap = true;
            //             break;
            //         }
            //     }
            //     // 2.b) apply correction for electron -> tau_h
            //     if (genElectronOverlap) {
            //         tau->pt *= 1.1;
            //     }
            // }
        }

        if( 
                tau->pt > 20  
                && abs(tau->eta) < 2.3 
                && !muOverlap
                && !elOverlap
                && (tau->hpsDisc & baconhep::kByDecayModeFinding)
                && (tau->hpsDisc & baconhep::kByTightIsolationMVA3newDMwLT) //<----------------------------- tau ID
                && (tau->hpsDisc & baconhep::kByMVA6VTightElectronRejection)
                && (tau->hpsDisc & baconhep::kByTightMuonRejection3)
            ) {
            taus.push_back(tau);
            veto_taus.push_back(tauP4);
        }
    }
    sort(taus.begin(), taus.end(), sort_by_higher_pt<TTau>);
    nTaus = taus.size();
    

    /* -------- PHOTONS --------  */
    vector<TLorentzVector> prefiring_photons;
    for (int i=0; i<fPhotonArr->GetEntries(); i++) {
        TPhoton* photon = (TPhoton*) fPhotonArr->At(i);
        TLorentzVector photonP4;
        photonP4.SetPtEtaPhiM(photon->pt, photon->eta, photon->phi, 0.);
        // prefiring_photons
        if( photon->pt>20 && fabs(photon->eta)<3.0 && fabs(photon->eta)>2.0 ){
            prefiring_photons.push_back(photonP4);
        }
    }


    /* -------- JETS ---------*/
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    //std::vector<TJet*> jets;
    std::vector<TJet*> vetoed_jets;
    vector<TLorentzVector> prefiring_jets;

    TLorentzVector hadronicP4;
    float sumJetPt = 0;

    // set nJet and nBJet counters to 0
    ResetJetCounters();
    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);
        
        // prefiring_jets
        if (jet->pt>20 && fabs(jet->eta) < 3.0 && fabs(jet->eta)>2.0 ){
            TLorentzVector jetP4_raw; 
            jetP4_raw.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
            prefiring_jets.push_back(jetP4_raw);
        }
        
        
        // apply JEC offline and get scale uncertainties
        double jec = particleSelector->JetCorrector(jet, "DummyName");
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
                vetoed_jets.push_back(jet);
                break;
            }
        }

        // save good jets
        if (
                fabs(jet->eta) <= 2.4
                && particleSelector->PassJetID(jet, cuts->looseJetID)
                && !muOverlap 
                && !elOverlap
                && !tauOverlap
            ) {
            if (isData){
                // std::cout<<jet->pt<<std::endl;
                // for data pass pt threshold
                if (jet->pt > 30){
                    ++nJets;
                    if (jet->csv > 0.8484) ++nBJets;
                }

            } else {
                // for MC apply corrections and variate systematics
                JetCounting(jet,jerc, gRand);
            }

            if (jet->pt > 30) {
                hadronicP4 += jetP4;
                sumJetPt += jetP4.Pt();
            }

        }
    }

    // sort jet by higher btag
    // std::sort(jets .begin(), jets .end(), sort_by_btag<TJet>);
    

    // use the highest jet multiplicities given all systematic variations
    unsigned nJetList[] = {nJets, nJetsJESUp, nJetsJESDown, nJetsJERUp, nJetsJERDown};
    unsigned nBJetList[] = {nBJets, nBJetsJESUp, nBJetsJESDown, nBJetsJERUp, nBJetsJERDown, 
                            nBJetsBTagUp, nBJetsBTagDown, nBJetsMistagUp, nBJetsMistagDown};
    // cut are maximum of counters which systematic variation
    nJetsCut  = *std::max_element(nJetList, nJetList+sizeof(nJetList)/sizeof(unsigned));
    nBJetsCut = *std::max_element(nBJetList, nBJetList+sizeof(nBJetList)/sizeof(unsigned));
    // std::cout<<nJetsCut<<", "<<nBJetsCut<<std::endl;
    // save ht
    htSum = sumJetPt;
    ht    = hadronicP4.Pt();
    htPhi = hadronicP4.Phi();

    /* -------- MET ---------*/
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;
 

    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    hTotalEvents->Fill(7);
    string channel = "";
    if        (nElectrons == 0 && nMuons == 2 && nTaus == 1 ){ // mu+mu selection

        channel = "mumutau";
        eventCounts[channel]->Fill(1);

        if (muons[0]->pt < 25. || muons[1]->pt < 10.)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        if (!muonTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        // if (nJetsCut < 2 || nBJetsCut < 1)
        //     return kTRUE;
        // eventCounts[channel]->Fill(4);




        // convert to TLorentzVectors
        TLorentzVector muonOneP4, muonTwoP4, dimuonP4;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.10566);
        muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, 0.10566);
        dimuonP4 = muonOneP4 + muonTwoP4;

        if  ( (muons[0]->q * muons[1]->q) >0)
            return kTRUE;
        eventCounts[channel]->Fill(5);

        if (dimuonP4.M()<12)
            return kTRUE;
        eventCounts[channel]->Fill(6);


        leptonOneP4     = muonOneP4;
        leptonOneIso    = GetMuonIsolation(muons[0]);
        leptonOneFlavor = muons[0]->q*13;
        leptonOneD0     = muons[0]->d0;
        leptonOneDz     = muons[0]->dz;
        leptonOneSip3d  = muons[0]->sip3d;
        leptonTwoP4     = muonTwoP4;
        leptonTwoIso    = GetMuonIsolation(muons[1]);
        leptonTwoFlavor = muons[1]->q*13;
        leptonTwoD0     = muons[1]->d0;
        leptonTwoDz     = muons[1]->dz;
        leptonTwoSip3d  = muons[1]->sip3d;

        // test if selected lepton fire the trigger.
        bitset<2> triggered;
        
        for (const auto& name: passTriggerNames) {
        if (trigger->passObj(name, 1, muons[0]->hltMatchBits) && muonOneP4.Pt() > 25)
            triggered.set(0);
        if (trigger->passObj(name, 1, muons[1]->hltMatchBits) && muonTwoP4.Pt() > 25)
            triggered.set(1);
        }
        if (!triggered.test(0) && !triggered.test(1)) {triggerLeptonStatus = 0;}
        if ( triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 1;
        if (!triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 2;
        if ( triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 3;
    
        eventCounts[channel]->Fill(5);


        ///////////tau info///////////////////
        TLorentzVector tauP4;
        tauP4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, taus[0]->m);
        leptonThreeP4     = tauP4;
        tauDecayMode      = taus[0]->decaymode;
        leptonThreeIso    = taus[0]->rawIsoMVA3newDMwLT;
        pair <float, float> tauVetoedJetPtPair;
        tauVetoedJetPtPair= GetTauVetoedJetPt(tauP4, vetoed_jets);
        tauVetoedJetPt    = tauVetoedJetPtPair.first;
        tauVetoedJetPtUnc = tauVetoedJetPtPair.second;
        //////////////////////////////////////

        
        // correct for MC, including reconstruction and trigger
        if (!isData) {

            tauGenFlavor    = GetTauGenFlavor(tauP4,genTausHad,genElectrons,genMuons,vetoed_jets, false); // useHadronFlavor = false
            tauGenFlavorHad = GetTauGenFlavor(tauP4,genTausHad,genElectrons,genMuons,vetoed_jets, true);  // useHadronFlavor = true

            // correct for id,reco weights
            EfficiencyContainer effCont;
            pair<float, float> effs, errs;

            // id weight
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

            eventWeight *= leptonOneIDWeight*leptonTwoIDWeight;

            // reconstruction weight
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

            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;

            // correct for trigger.

            // check if lepton could pass the trigger threshold and is matched
            // to a trigger object.  When both muons pass the trigger, use the
            // efficiency for detecting either
            
            EfficiencyContainer effCont1, effCont2;
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
            eventWeight *= triggerWeight;
        }
    } else if (nElectrons == 2 && nMuons == 0 && nTaus == 1) { // e+e selection
        channel = "eetau";
        eventCounts[channel]->Fill(1);
        // coutchannel


        if (electrons[0]->pt < 30 || electrons[1]->pt < 10 )
            return kTRUE;
        eventCounts[channel]->Fill(2);

        if (!electronTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        // if (nJetsCut < 2 || nBJetsCut < 1)
        //     return kTRUE;
        // eventCounts[channel]->Fill(4);


        TLorentzVector electronOneP4, electronTwoP4, dielectronP4;
        electronOneP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        electronTwoP4.SetPtEtaPhiM(electrons[1]->pt, electrons[1]->eta, electrons[1]->phi, 511e-6);
        dielectronP4 = electronOneP4 + electronTwoP4;

        if  ( (electrons[0]->q * electrons[1]->q) >0)
            return kTRUE;
        eventCounts[channel]->Fill(5);

        if (dielectronP4.M()<12)
            return kTRUE;
        eventCounts[channel]->Fill(6);

        leptonOneP4     = electronOneP4;
        leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonOneFlavor = 11*electrons[0]->q;
        leptonOneD0     = electrons[0]->d0;
        leptonOneDz     = electrons[0]->dz;
        leptonOneSip3d  = electrons[0]->sip3d;

        leptonTwoP4     = electronTwoP4;
        leptonTwoIso    = GetElectronIsolation(electrons[1], fInfo->rhoJet);
        leptonTwoFlavor = 11*electrons[1]->q;
        leptonTwoD0     = electrons[1]->d0;
        leptonTwoDz     = electrons[1]->dz;
        leptonTwoSip3d  = electrons[1]->sip3d;



        // test if selected lepton fire the trigger.
        bitset<2> triggered;


        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, electrons[0]->hltMatchBits) && electronOneP4.Pt() > 30)
                triggered.set(0);
            if (trigger->passObj(name, 1, electrons[1]->hltMatchBits) && electronTwoP4.Pt() > 30)
                triggered.set(1);
        }
        if (!triggered.test(0) && !triggered.test(1)) {triggerLeptonStatus = 0;}
        if ( triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 1;
        if (!triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 2;
        if ( triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 3;
    
        eventCounts[channel]->Fill(5);


        ///////////tau info///////////////////
        TLorentzVector tauP4;
        tauP4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, taus[0]->m);
        leptonThreeP4     = tauP4;
        tauDecayMode      = taus[0]->decaymode;
        leptonThreeIso    = taus[0]->rawIsoMVA3newDMwLT;
        pair <float, float> tauVetoedJetPtPair;
        tauVetoedJetPtPair= GetTauVetoedJetPt(tauP4, vetoed_jets);
        tauVetoedJetPt    = tauVetoedJetPtPair.first;
        tauVetoedJetPtUnc = tauVetoedJetPtPair.second;
        //////////////////////////////////////

        // correct for MC, including reconstruction and trigger
        if (!isData) {
            tauGenFlavor    = GetTauGenFlavor(tauP4,genTausHad,genElectrons,genMuons,vetoed_jets, false); // useHadronFlavor = false
            tauGenFlavorHad = GetTauGenFlavor(tauP4,genTausHad,genElectrons,genMuons,vetoed_jets, true);  // useHadronFlavor = true

            // correct for id,reco weights
            EfficiencyContainer effCont;
            pair<float, float> effs, errs;

            // id weight
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

            eventWeight *= leptonOneIDWeight*leptonTwoIDWeight;

            // reconstruction weight
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

            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;


            // correct for trigger.

            // check if lepton could pass the trigger threshold and is matched
            // to a trigger object.  When both muons pass the trigger, use the
            // efficiency for detecting either
            EfficiencyContainer effCont1, effCont2;
            if (triggered.all()) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronOneP4);
                effCont2      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronTwoP4);
                triggerWeight = GetTriggerSF(effCont1, effCont2);
                triggerVar    = GetTriggerSFError(effCont1, effCont2);

                eleTriggerVarTagSyst   = weights->GetEleTriggerSyst("tag", electronOneP4);
                eleTriggerVarProbeSyst = weights->GetEleTriggerSyst("probe", electronOneP4);
            } else if (triggered.test(0)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronOneP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

                eleTriggerVarTagSyst   = weights->GetEleTriggerSyst("tag", electronOneP4);
                eleTriggerVarProbeSyst = weights->GetEleTriggerSyst("probe", electronOneP4);

            } else if (triggered.test(1)) {
                effCont1      = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electronTwoP4);
                effs          = effCont1.GetEff();
                errs          = effCont1.GetErr();
                triggerWeight = effs.first/effs.second;
                triggerVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));

                eleTriggerVarTagSyst   = weights->GetEleTriggerSyst("tag", electronTwoP4);
                eleTriggerVarProbeSyst = weights->GetEleTriggerSyst("probe", electronTwoP4);
            } else {
                return kTRUE;
            }
            eventWeight *= triggerWeight;

            // l1 prefiring
            effCont = weights->GetElectronPrefiringWeight(prefiring_photons, prefiring_jets);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            prefiringWeight = effs.first/effs.second;
            prefiringVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            eventWeight *= prefiringWeight;
            

        }
    } else if (nElectrons == 1 && nMuons == 1 && nTaus == 1) { // e+mu selection
        channel = "emutau";
        eventCounts[channel]->Fill(1);
        // coutchannel


        // trigger matching for thresholds
        
        float muPtThreshold = 10.;
        float elPtThreshold = 10.;

        if ( muonTriggered ) { muPtThreshold = 25.; } 
        else if (electronTriggered) {elPtThreshold = 30.;}
        // if (electronTriggered) {elPtThreshold = 30.;}

        if ( muons[0]->pt < muPtThreshold || electrons[0]->pt < elPtThreshold)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        if (!(muonTriggered || electronTriggered))
            return kTRUE;
        eventCounts[channel]->Fill(3);

        // if (nJetsCut < 2 || nBJetsCut < 1)
        //     return kTRUE;
        // eventCounts[channel]->Fill(4);


        TLorentzVector muonP4, electronP4, dielectronP4;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.10566);
        electronP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        dielectronP4 = muonP4 + electronP4;


        leptonOneP4     = muonP4;
        leptonOneIso    = GetMuonIsolation(muons[0]);
        leptonOneFlavor = 13*muons[0]->q;
        leptonOneD0     = muons[0]->d0;
        leptonOneDz     = muons[0]->dz;
        leptonOneSip3d  = muons[0]->sip3d;

        leptonTwoP4     = electronP4;
        leptonTwoIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonTwoFlavor = 11*electrons[0]->q;
        leptonTwoD0     = electrons[0]->d0;
        leptonTwoDz     = electrons[0]->dz;
        leptonTwoSip3d  = electrons[0]->sip3d;

        // test if selected lepton fire the trigger.
        bitset<2> triggered;        


        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, muons[0]->hltMatchBits))
                triggered.set(0);
            if (trigger->passObj(name, 1, electrons[0]->hltMatchBits))
                triggered.set(1);
        }
        if (!triggered.test(0) && !triggered.test(1)) {triggerLeptonStatus = 0;}
        if ( triggered.test(0) && !triggered.test(1)) triggerLeptonStatus = 1;
        if (!triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 2;
        if ( triggered.test(0) &&  triggered.test(1)) triggerLeptonStatus = 3;
    

        eventCounts[channel]->Fill(5);


        ///////////tau info///////////////////
        TLorentzVector tauP4;
        tauP4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, taus[0]->m);
        leptonThreeP4     = tauP4;
        tauDecayMode      = taus[0]->decaymode;
        leptonThreeIso    = taus[0]->rawIsoMVA3newDMwLT;
        pair <float, float> tauVetoedJetPtPair;
        tauVetoedJetPtPair= GetTauVetoedJetPt(tauP4, vetoed_jets);
        tauVetoedJetPt    = tauVetoedJetPtPair.first;
        tauVetoedJetPtUnc = tauVetoedJetPtPair.second;
        //////////////////////////////////////

        // correct for MC, including reconstruction and trigger
        if (!isData) {

            tauGenFlavor    = GetTauGenFlavor(tauP4,genTausHad,genElectrons,genMuons,vetoed_jets, false); // useHadronFlavor = false
            tauGenFlavorHad = GetTauGenFlavor(tauP4,genTausHad,genElectrons,genMuons,vetoed_jets, true);  // useHadronFlavor = true


            // correct for id,reco weights
            EfficiencyContainer effCont;
            pair<float, float> effs, errs;

            // id weight
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

            eventWeight *= leptonOneIDWeight*leptonTwoIDWeight;

            // reconstruction weight
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

            eventWeight *= leptonOneRecoWeight*leptonTwoRecoWeight;

            // correct for trigger.

            // check if lepton could pass the trigger threshold and is matched
            // to a trigger object.  When both muons pass the trigger, use the
            // efficiency for detecting either
            EfficiencyContainer effCont1, effCont2;
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

                eleTriggerVarTagSyst   = weights->GetEleTriggerSyst("tag", electronP4);
                eleTriggerVarProbeSyst = weights->GetEleTriggerSyst("probe", electronP4);

                
            } else {
                return kTRUE;
            }
            eventWeight *= triggerWeight;

            // l1 prefiring
            if(!muonTriggered) { // no prefiring effect if triggerred on muon HLT
                effCont = weights->GetElectronPrefiringWeight(prefiring_photons, prefiring_jets);
                effs = effCont.GetEff();
                errs = effCont.GetErr();
                prefiringWeight = effs.first/effs.second;
                prefiringVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
                eventWeight *= prefiringWeight;
            }
            

        }
    } else {
        return kTRUE;
    }


    ///////////////////
    // Fill jet info //
    ///////////////////

    outFile->cd(channel.c_str());
    outTrees[channel]->Fill();
    this->passedEvents++;
    return kTRUE;
}

void LLTauAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void LLTauAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void LLTauAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<LLTauAnalyzer> selector(new LLTauAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


float LLTauAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    return combIso;
}

float LLTauAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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





int LLTauAnalyzer::GetTauGenFlavor(  TLorentzVector p4, 
                                        vector<TGenParticle*> genTausHad, 
                                        vector<TGenParticle*> genElectrons, 
                                        vector<TGenParticle*> genMuons,
                                        vector<TJet*> vetoed_jets, 
                                        bool useHadronFlavor )
{
    int flavor = 26;
    
    // check if can be tagged as hadronic tau
    if (flavor==26){
        for (unsigned i = 0; i < genTausHad.size(); ++i) {
            TLorentzVector genP4;
            genP4.SetPtEtaPhiM(genTausHad[i]->pt, genTausHad[i]->eta, genTausHad[i]->phi, genTausHad[i]->mass); 
            if (genP4.DeltaR(p4) < 0.3) {
                flavor = 15;
            }
        }
    }

    
    // check if can be tagged as electron
    if (flavor==26){
        for (unsigned i = 0; i < genElectrons.size(); ++i) {
            TLorentzVector genP4;
            genP4.SetPtEtaPhiM(genElectrons[i]->pt, genElectrons[i]->eta, genElectrons[i]->phi, genElectrons[i]->mass); 
            if (genP4.DeltaR(p4) < 0.3) {
                flavor = 11;
            }
        }
    }

    // check if can be tagged as muon
    if (flavor==26){
        for (unsigned i = 0; i < genMuons.size(); ++i) {
            TLorentzVector genP4;
            genP4.SetPtEtaPhiM(genMuons[i]->pt, genMuons[i]->eta, genMuons[i]->phi, genMuons[i]->mass); 
            if (genP4.DeltaR(p4) < 0.3) {
                flavor = 13;
            }
        }
    }


    // check if can be tagged by jet flavor
    if (flavor==26){
        float jetPtMax = - 1.0;
        for (unsigned i = 0; i < vetoed_jets.size(); ++i) {


            TLorentzVector jetP4; 
            jetP4.SetPtEtaPhiM(vetoed_jets[i]->pt, vetoed_jets[i]->eta, vetoed_jets[i]->phi, vetoed_jets[i]->mass);
            if (jetP4.DeltaR(p4) < 0.4 && jetP4.Pt()>jetPtMax) {
                if(useHadronFlavor) {
                    flavor = vetoed_jets[i]->hadronFlavor;
                } else {
                    flavor = abs(vetoed_jets[i]->partonFlavor);
                }
                jetPtMax = jetP4.Pt();
            }
        }
    }



    return flavor;
}


pair<float, float> LLTauAnalyzer::GetTauVetoedJetPt(TLorentzVector p4, vector<TJet*> vetoed_jets)
{

    float jetPt = -1;
    float jetPtUnc = -1;
    

    // check if can be tagged by jet flavor

    float jetPtMax = - 1.0;
    for (unsigned i = 0; i < vetoed_jets.size(); ++i) {

        TLorentzVector jetP4; 
        jetP4.SetPtEtaPhiM(vetoed_jets[i]->pt, vetoed_jets[i]->eta, vetoed_jets[i]->phi, vetoed_jets[i]->mass);

        if (jetP4.DeltaR(p4) < 0.4 && jetP4.Pt()>jetPtMax) {
            // cout<<"Energy Jet: "<<jetPtMax<<","<< jetP4.Pt() << endl;

            // save 
            jetPt       = jetP4.Pt();
            jetPtUnc    = jetPt * float(particleSelector->JetUncertainty(vetoed_jets[i], "Total"));
            jetPtUnc    = abs(jetPtUnc);

            jetPtMax = jetP4.Pt();
        }    
    }
    
    return make_pair(jetPt,jetPtUnc);
}




float LLTauAnalyzer::GetTriggerSF(EfficiencyContainer eff1, EfficiencyContainer eff2)
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

float LLTauAnalyzer::GetTriggerSFError(EfficiencyContainer eff1, EfficiencyContainer eff2)
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

void LLTauAnalyzer::ResetJetCounters()
{
    nJets = nBJets = 0;
    nJetsCut       = nBJetsCut        = 0;
    nJetsJESUp     = nJetsJESDown     = 0;
    nJetsJERUp     = nJetsJERDown     = 0;
    nBJetsJESUp    = nBJetsJESDown    = 0;
    nBJetsJERUp    = nBJetsJERDown    = 0;
    nBJetsBTagUp   = nBJetsBTagDown   = 0;
    nBJetsMistagUp = nBJetsMistagDown = 0;
}

void LLTauAnalyzer::JetCounting(TJet* jet, float jerc_nominal, float resRand)
{
    float jetPt = jet->pt;
    //std::string bTagMethod = "MVAT";
    std::string bTagMethod = "CSVM";

    float rNumber = rng->Uniform(1.);
    if (jet->pt > 30) {

        // nominal
        ++nJets;
        if (particleSelector->BTagModifier(jet, bTagMethod, "central", rNumber)) ++nBJets;

        // b tag up
        if (particleSelector->BTagModifier(jet, bTagMethod, "up", rNumber)) ++nBJetsBTagUp;
             
        // b tag down
        if (particleSelector->BTagModifier(jet, bTagMethod, "down", rNumber))  ++nBJetsBTagDown;

        // misttag up
        if (particleSelector->BTagModifier(jet, bTagMethod, "upMistag", rNumber)) ++nBJetsMistagUp;
        
        // mistag down
        if (particleSelector->BTagModifier(jet, bTagMethod, "downMistag", rNumber)) ++nBJetsMistagDown;
    }

    // jet energy corrections
    double jec = particleSelector->JetCorrector(jet, "DummyName");
    float jecUnc = particleSelector->JetUncertainty(jet,"Total");
    
    // JES up
    jet->pt = jet->ptRaw*jec*(1 + jecUnc)*jerc_nominal;
    if (jet->pt > 30) {
        ++nJetsJESUp;
        if (particleSelector->BTagModifier(jet, bTagMethod, "central", rNumber)) { 
            ++nBJetsJESUp;
        } 
    }

    // JES down
    jet->pt = jet->ptRaw*jec*(1 - jecUnc)*jerc_nominal;
    if (jet->pt > 30) {
        ++nJetsJESDown;
        if (particleSelector->BTagModifier(jet, bTagMethod, "central", rNumber)) { 
            ++nBJetsJESDown;
        } 
    }

    // Jet resolution corrections
    pair<float, float> resPair = particleSelector->JetResolutionAndSF(jet, 1);
    float jerc = 1 + resRand*sqrt(std::max((double)resPair.second*resPair.second - 1, 0.));

    // JER up
    jet->pt = jet->ptRaw*jec*jerc;
    if (jet->pt > 30) {
        ++nJetsJERUp;
        if (particleSelector->BTagModifier(jet, bTagMethod,"central", rNumber)) { 
            ++nBJetsJERUp;
        }     
    }

    // JER down
    resPair     = particleSelector->JetResolutionAndSF(jet, -1);
    jerc        = 1 + resRand*sqrt(std::max((double)resPair.second*resPair.second - 1, 0.));
    jet->pt     = jet->ptRaw*jec*jerc;
    if (jet->pt > 30) {
        ++nJetsJERDown;
        if (particleSelector->BTagModifier(jet, bTagMethod, "central", rNumber)) { 
            ++nBJetsJERDown;
        } 
    }
    jet->pt = jetPt;
}
