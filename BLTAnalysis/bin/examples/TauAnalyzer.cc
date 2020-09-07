#include "TauAnalyzer.h"
#include <map>


//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;



TauAnalyzer::TauAnalyzer() : BLTSelector()
{

}

TauAnalyzer::~TauAnalyzer()
{

}

void TauAnalyzer::Begin(TTree *tree)
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

    if (params->selection == "single_lepton") {
        triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");

    } else if (params->selection == "single_muon") {
        triggerNames.push_back("HLT_IsoMu24_v*");
        triggerNames.push_back("HLT_IsoTkMu24_v*");

    } else if (params->selection == "single_electron") {
        triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");
    }

    // Weight utility class
    cout<< "Weight utility class"<<endl;
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"; // 2016 mask
    lumiMask.AddJSONFile(jsonFileName);
    

    // muon momentum corrections
    cout<< "muon momentum corrections"<<endl;
    muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/rcdata.2016.v3");

    // electron scale corrections
    cout<< "electron scale corrections"<<endl;
    electronScaler = new EnergyScaleCorrection(cmssw_base + "/src/BLT/BLTAnalysis/data");

    // Prepare the output tree
    cout<< "Prepare the output tree"<<endl;
    string outFileName = params->get_output_filename("output");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();

    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);

    vector<std::string> channelNames = {"tau"};

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
        tree->Branch("mcEra",&mcEra);
        
        // weights and their uncertainties
        tree->Branch("eventWeight", &eventWeight);
        tree->Branch("genWeight", &genWeight);
        tree->Branch("puWeight", &puWeight);
        tree->Branch("topPtWeight", &topPtWeight);

        // tau staff
        tree->Branch("tauPt",     &tauPt);
        tree->Branch("tauEta",     &tauEta);
        tree->Branch("tauDecayMode",     &tauDecayMode);
        tree->Branch("tauMVA",           &tauMVA); 
        tree->Branch("tauGenFlavor",     &tauGenFlavor);
        tree->Branch("tauGenFlavorHad",  &tauGenFlavorHad);
        tree->Branch("tauVetoedJetPt",   &tauVetoedJetPt);
        tree->Branch("tauVetoedJetPtUnc",&tauVetoedJetPtUnc);
      
        // object counters
        tree->Branch("nMuons", &nMuons);
        tree->Branch("nElectrons", &nElectrons);
        tree->Branch("nTaus", &nTaus);
        tree->Branch("nJets", &nJets);
        tree->Branch("nBJets", &nBJets);
        // ht,met
        tree->Branch("htSum", &htSum);
        tree->Branch("ht", &ht);
        tree->Branch("htPhi", &htPhi);
        tree->Branch("met", &met);
        tree->Branch("metPhi", &metPhi);


        outTrees[channel] = tree;
    }
    cout<< "Done with initialization"<<endl;
    ReportPostBegin();

}

Bool_t TauAnalyzer::Process(Long64_t entry)
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

        //   1.d) save gen level tau and leptonic taus, and e.mu
        vector<TGenParticle*> genTaus, genTausLep;
        
        // 2. loop over gen particles
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
            cout<<i << particle->pdgId<<endl;

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

            // 2.d) This saves gen level tau, leptonic taus, electrons and muons
            if ( particle->parent >=0 ) {
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);

                if (abs(particle->pdgId) != 22 && abs(particle->pdgId) != 15 && abs(mother->pdgId) == 15 && mother->parent >=0) {
                    TGenParticle* gmother = (TGenParticle*) fGenParticleArr->At(mother->parent);
                    bool isHard = gmother->status != 2 && gmother->parent != -2;

                    if (isHard && abs(particle->pdgId) == 16) genTaus.push_back(mother); 
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
        } cout<<"finish loop"<<endl;

        // 3. finish loop over gen particles. now conclude gen level infomation

        // 3.a) counting partons: save to nPartons
        nPartons = partonCount; 

        // 3.b) Account for the top pt weights
        if (params->datasetgroup.substr(0, 5) == "ttbar") {
            topPtWeight *= sqrt(topSF);
            topPtVar    += pow(0.01*topPtWeight, 2);
            eventWeight *= topPtWeight;
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
    } else {
        genWeight = 1.0;
        nPU = 0;
        nPartons = 0;
        genCategory = 0;

    }


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
    nMuons = muons.size();
    nFailMuons = fail_muons.size();
    cout<<"finish muon loop"<<endl;
    




    
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

    nElectrons = electrons.size();
    nFailElectrons = fail_electrons.size();
    cout<<"finish electron loop"<<endl;




    
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
    nTaus = taus.size();
    cout<<"finish tau loop"<<endl;
    

    
    /* -------- JETS ---------*/
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets, vetoedJets;
    TLorentzVector hadronicP4;
    float sumJetPt = 0;

    // set nJet and nBJet counters to 0
    nJets = 0;
    nBJets = 0;
    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

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
                vetoedJets.push_back(jet);
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

            jets.push_back(jet);
            ++nJets;
            hadronicP4 += jetP4;
            sumJetPt += jetP4.Pt();

            if (isData){ 
                if (jet->csv > 0.8484) ++nBJets;
            } else {
                float rNumber = rng->Uniform(1.);
                if (particleSelector->BTagModifier(jet, "CSVM", "central", rNumber)) ++nBJets;
            }

        }
    }

    cout<<"finish jet loop"<<endl;
    // sort jet by higher btag
    // std::sort(jets .begin(), jets .end(), sort_by_btag<TJet>);


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


    hTotalEvents->Fill(5);
    string channel = "";
    if        (nTaus > 0 ){
        
      channel = "tau";

    //   TLorentzVector tauP4;
    //   tauP4.SetPtEtaPhiM(taus[0]->pt, taus[0]->eta, taus[0]->phi, taus[0]->m);
    //   tauPt  = taus[0]->pt;
    //   tauEta = taus[0]->eta;
    //   tauDecayMode      = taus[0]->decaymode;
    //   tauMVA            = taus[0]->rawIsoMVA3newDMwLT;

    //   pair <float, float> tauVetoedJetPtPair;
    //   tauVetoedJetPtPair= GetTauVetoedJetPt(tauP4, vetoedJets);
    //   tauVetoedJetPt    = tauVetoedJetPtPair.first;
    //   tauVetoedJetPtUnc = tauVetoedJetPtPair.second;

    //   if (!isData) {
    //     tauGenFlavor    = GetTauGenFlavor(tauP4,genTausHad,genElectrons,genMuons,vetoedJets, false); // useHadronFlavor = false
    //     tauGenFlavorHad = GetTauGenFlavor(tauP4,genTausHad,genElectrons,genMuons,vetoedJets, true);  // useHadronFlavor = true
    //   }

      hTotalEvents->Fill(6);
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

void TauAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void TauAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void TauAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<TauAnalyzer> selector(new TauAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


float TauAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    return combIso;
}

float TauAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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





int TauAnalyzer::GetTauGenFlavor(  TLorentzVector p4, 
                                        vector<TGenParticle*> genTausHad, 
                                        vector<TGenParticle*> genElectrons, 
                                        vector<TGenParticle*> genMuons,
                                        vector<TJet*> vetoedJets, 
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
        for (unsigned i = 0; i < vetoedJets.size(); ++i) {


            TLorentzVector jetP4; 
            jetP4.SetPtEtaPhiM(vetoedJets[i]->pt, vetoedJets[i]->eta, vetoedJets[i]->phi, vetoedJets[i]->mass);
            if (jetP4.DeltaR(p4) < 0.4 && jetP4.Pt()>jetPtMax) {
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


pair<float, float> TauAnalyzer::GetTauVetoedJetPt(TLorentzVector p4, vector<TJet*> vetoedJets)
{

    float jetPt = -1;
    float jetPtUnc = -1;
    

    // check if can be tagged by jet flavor

    float jetPtMax = - 1.0;
    for (unsigned i = 0; i < vetoedJets.size(); ++i) {

        TLorentzVector jetP4; 
        jetP4.SetPtEtaPhiM(vetoedJets[i]->pt, vetoedJets[i]->eta, vetoedJets[i]->phi, vetoedJets[i]->mass);

        if (jetP4.DeltaR(p4) < 0.4 && jetP4.Pt()>jetPtMax) {
            // cout<<"Energy Jet: "<<jetPtMax<<","<< jetP4.Pt() << endl;

            // save 
            jetPt       = jetP4.Pt();
            jetPtUnc    = jetPt * float(particleSelector->JetUncertainty(vetoedJets[i], "Total"));
            jetPtUnc    = abs(jetPtUnc);

            jetPtMax = jetP4.Pt();
        }    
    }
    
    return make_pair(jetPt,jetPtUnc);
}


