#include "GenTauTauAnalyzer.h"
#include <map>

// GenTauTauAnalyzer /eos/uscms/store/group/lpcbacon/12/Summer16_DYJetsToLL_M-50_amcatnlo/Summer16_DYJetsToLL_M-50_amcatnlo_bacon_00.root 100000 dy dy single_lepton 2016 1

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;



GenTauTauAnalyzer::GenTauTauAnalyzer() : BLTSelector()
{

}

GenTauTauAnalyzer::~GenTauTauAnalyzer()
{

}

void GenTauTauAnalyzer::Begin(TTree *tree)
{   
    cout << __func__ << endl;
    // rng = new TRandom3();


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
    // cuts.reset(new Cuts());
    // particleSelector.reset(new ParticleSelector(*params, *cuts));

    // Trigger bits mapping file
    cout<< "Trigger bits mapping file"<<endl;
    // const std::string cmssw_base = getenv("CMSSW_BASE");
    // std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns";
    // trigger.reset(new baconhep::TTrigger(trigfilename));

    // if (params->selection == "single_lepton") {
    //     triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");
    //     triggerNames.push_back("HLT_IsoMu24_v*");
    //     triggerNames.push_back("HLT_IsoTkMu24_v*");
    //     // triggerNames.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*");
    //     // triggerNames.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");

    // } else if (params->selection == "single_muon") {
    //     triggerNames.push_back("HLT_IsoMu24_v*");
    //     triggerNames.push_back("HLT_IsoTkMu24_v*");

    // } else if (params->selection == "single_electron") {
    //     triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");
    // }

    // Weight utility class
    cout<< "Weight utility class"<<endl;
    // weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask
    // // Set up object to handle good run-lumi filtering if necessary
    // lumiMask = RunLumiRangeMap();
    // string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"; // 2016 mask
    // lumiMask.AddJSONFile(jsonFileName);
    

    // muon momentum corrections
    cout<< "Muon momentum corrections"<<endl;
    // muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/rcdata.2016.v3");

    // electron scale corrections
    cout<< "Electron scale corrections"<<endl;
    // electronScaler = new EnergyScaleCorrection(cmssw_base + "/src/BLT/BLTAnalysis/data");

    // Prepare the output tree
    cout<< "Prepare the output tree"<<endl;
    string outFileName = params->get_output_filename("output");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();

    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);

    vector<std::string> channelNames = {"etau","mutau","emu"};

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
        // tree->Branch("mcEra",&mcEra);
        // tree->Branch("triggerLeptonStatus",&triggerLeptonStatus);
        
        // weights and their uncertainties
        tree->Branch("eventWeight", &eventWeight);
        tree->Branch("genWeight", &genWeight);
        // tree->Branch("puWeight", &puWeight);
        // tree->Branch("topPtWeight", &topPtWeight);

        // tau staff
        if (channel == "etau" || channel == "mutau"){
          tree->Branch("tauDecayMode",     &tauDecayMode);
          // tree->Branch("tauMVA",           &tauMVA); 
          tree->Branch("tauGenFlavor",     &tauGenFlavor);
          tree->Branch("tauGenFlavorHad",  &tauGenFlavorHad);
          // tree->Branch("tauVetoedJetPt",   &tauVetoedJetPt);
          // tree->Branch("tauVetoedJetPtUnc",&tauVetoedJetPtUnc);
        }

        // vectors
        tree->Branch("leptonOneP4",     &leptonOneP4);
        tree->Branch("leptonTwoP4",     &leptonTwoP4);
        tree->Branch("leptonOneFlavor",     &leptonOneFlavor);
        tree->Branch("leptonTwoFlavor",     &leptonTwoFlavor);
        // tree->Branch("photonP4",        &photonP4);

      
        // object counters
        tree->Branch("nMuons", &nMuons);
        tree->Branch("nElectrons", &nElectrons);
        tree->Branch("nTaus", &nTaus);
        // tree->Branch("nPhotons", &nPhotons);
        // tree->Branch("nJets", &nJets);
        // tree->Branch("nBJets", &nBJets);
        tree->Branch("nGenTausHad", &nGenTausHad);
        tree->Branch("nGenTausLep", &nGenTausLep);
        tree->Branch("nGenElectrons", &nGenElectrons);
        tree->Branch("nGenMuons", &nGenMuons);

        // ht,met
        // tree->Branch("htSum", &htSum);
        // tree->Branch("ht", &ht);
        // tree->Branch("htPhi", &htPhi);
        // tree->Branch("met", &met);
        // tree->Branch("metPhi", &metPhi);
        // tree->Branch("covMet00", &covMet00);
        // tree->Branch("covMet01", &covMet01);
        // tree->Branch("covMet11", &covMet11);

	//SVFit variables
	// tree->Branch("massSVFit", &massSVFit);
	// tree->Branch("massErrSVFit", &massErrSVFit);
	// tree->Branch("svFitStatus", &svFitStatus);

        outTrees[channel] = tree;
        // event counter
        string outHistName = params->get_output_treename("TotalEvents_" + channel);
        eventCounts[channel] = new TH1D(outHistName.c_str(),"ChannelCounts",10,0.5,10.5);

    }
    cout<< "Done with initialization"<<endl;
    ReportPostBegin();

}

Bool_t GenTauTauAnalyzer::Process(Long64_t entry)
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

    //Whether or not to filter out events where selected
    //lepton didn't fire the trigger
    const bool isData = (fInfo->runNum != 1);
    // particleSelector->SetRealData(isData);

    /* Data not defined for gen info */
    if (isData) {
      return kTRUE;
    }
    hTotalEvents->Fill(2);

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

    // save gen weight for amc@nlo Drell-Yan sample
    genWeight = fGenEvtInfo->weight > 0 ? 1 : -1;
    if (genWeight < 0) {
      hTotalEvents->Fill(10);
    }


    // pileup reweighting
    nPU          = fInfo->nPUmean;
    // puWeight     = weights->GetPUWeight(fInfo->nPUmean); 
    // eventWeight *= puWeight;


    // loop over gen particle collection:

    // 1.initialization
    //   1.a) counting parton for Drell-Yan samples
    unsigned partonCount= 0;

    //   1.b) counting top for top reweighting
    float topSF = 1.;
    unsigned topCount   = 0;
    // topPtWeight = 1.;
    // topPtVar = 0.;

    //   1.d) save gen level tau and leptonic taus, and e.mu
    vector<TGenParticle*> genTaus, genTausLep;
    vector<TGenParticle*> genTausHad, genElectrons, genMuons;

        
        
    // 2. loop over gen particles
    for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
      TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

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
    }

    // 3. finish loop over gen particles. now conclude gen level infomation

    // 3.a) counting partons: save to nPartons
    nPartons = partonCount; 

    // 3.b) Account for the top pt weights
    // if (params->datasetgroup.substr(0, 5) == "ttbar") {
    //   topPtWeight *= sqrt(topSF);
    //   topPtVar    += pow(0.01*topPtWeight, 2);
    //   eventWeight *= topPtWeight;
    // } 
        
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
    nGenTausHad = genTausHad.size();
    nGenTausLep = genTausLep.size();
    nGenElectrons = genElectrons.size();
    nGenMuons = genMuons.size();


    ///////////////////
    // Select objects//
    ///////////////////

    //vertex check, skipped
    hTotalEvents->Fill(3);

    //trigger check, skipped
    hTotalEvents->Fill(4);

    vector<TGenParticle*> good_muons;
    for(unsigned int i = 0; i < genMuons.size(); ++i) {
      if(
	 genMuons[i]->pt > 10.
	 && fabs(genMuons[i]->eta) < 2.4
	 ) {
	good_muons.push_back(genMuons[i]);
      }
    }

    vector<TGenParticle*> good_electrons;
    for(unsigned int i = 0; i < genElectrons.size(); ++i) {
      if(
	 genElectrons[i]->pt > 15.
	 && fabs(genElectrons[i]->eta) < 2.5
	 ) {
	good_electrons.push_back(genElectrons[i]);
      }
    }

    vector<TGenParticle*> good_taus;
    for(unsigned int i = 0; i < genTaus.size(); ++i) {
      if(
	 genTaus[i]->pt > 20.
	 && fabs(genTaus[i]->eta) < 2.3
	 ) {
	good_taus.push_back(genTaus[i]);
      }
    }

    
    sort(good_muons.begin(), good_muons.end(), sort_by_higher_pt<TGenParticle>);
    sort(good_electrons.begin(), good_electrons.end(), sort_by_higher_pt<TGenParticle>);
    sort(good_taus.begin(), good_taus.end(), sort_by_higher_pt<TGenParticle>);    

    nElectrons = good_electrons.size();
    nMuons = good_muons.size();
    nTaus = good_taus.size();
    
    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////


    hTotalEvents->Fill(5);
    string channel = "";
    if (nElectrons == 1 && nMuons == 0 && nTaus == 1 ) { // e+tau selection

        channel = "etau";
        eventCounts[channel]->Fill(1);

        if (good_electrons[0]->pt < 30. )
            return kTRUE;
        eventCounts[channel]->Fill(2);

	//trigger check, skip
        eventCounts[channel]->Fill(3);



        TLorentzVector electronP4, tauP4;
        electronP4.SetPtEtaPhiM(good_electrons[0]->pt, good_electrons[0]->eta, good_electrons[0]->phi, 511e-6);
        tauP4.SetPtEtaPhiM(good_taus[0]->pt, good_taus[0]->eta, good_taus[0]->phi, good_taus[0]->mass);

        leptonOneP4     = electronP4;
        leptonOneFlavor = good_electrons[0]->pdgId;
        leptonTwoP4     = tauP4;
        leptonTwoFlavor = good_taus[0]->pdgId;


        ///////////tau info///////////////////
        // tauDecayMode      = taus[0]->decaymode;
        // tauMVA            = taus[0]->rawIsoMVA3newDMwLT;
        // pair <float, float> tauVetoedJetPtPair;
        // tauVetoedJetPtPair= GetTauVetoedJetPt(tauP4, vetoedJets);
        // tauVetoedJetPt    = tauVetoedJetPtPair.first;
        // tauVetoedJetPtUnc = tauVetoedJetPtPair.second;
        //////////////////////////////////////


        // test if selected lepton fire the trigger.
        eventCounts[channel]->Fill(4);

        // correct for MC, including reconstruction and trigger

	eventCounts[channel]->Fill(5);

    } else if (nElectrons == 0 && nMuons == 1 && nTaus == 1 ) { // mu+tau selection
        channel = "mutau";
        eventCounts[channel]->Fill(1);

        if (good_muons[0]->pt < 25. )
            return kTRUE;
        eventCounts[channel]->Fill(2);

	//trigger check, skip
        eventCounts[channel]->Fill(3);

        TLorentzVector muonP4, tauP4;
        muonP4.SetPtEtaPhiM(good_muons[0]->pt, good_muons[0]->eta, good_muons[0]->phi, 0.10566);
        tauP4.SetPtEtaPhiM(good_taus[0]->pt, good_taus[0]->eta, good_taus[0]->phi, good_taus[0]->mass);

        leptonOneP4     = muonP4;
        leptonOneFlavor = good_muons[0]->pdgId;
        leptonTwoP4     = tauP4;
        leptonTwoFlavor = good_taus[0]->pdgId;


        ///////////tau info///////////////////
        // tauDecayMode      = taus[0]->decaymode;
        // tauMVA            = taus[0]->rawIsoMVA3newDMwLT;
        // pair <float, float> tauVetoedJetPtPair;
        // tauVetoedJetPtPair= GetTauVetoedJetPt(tauP4, vetoedJets);
        // tauVetoedJetPt    = tauVetoedJetPtPair.first;
        // tauVetoedJetPtUnc = tauVetoedJetPtPair.second;
        //////////////////////////////////////


        // test if selected lepton fire the trigger.
        eventCounts[channel]->Fill(4);

        // correct for MC, including reconstruction and trigger
	eventCounts[channel]->Fill(5);

    } else if (nElectrons == 1 && nMuons == 1  ) { // e+mu selection

        channel = "emu";
        eventCounts[channel]->Fill(1);

	//trigger check, skip
        eventCounts[channel]->Fill(2);


        TLorentzVector electronP4, muonP4;
        electronP4.SetPtEtaPhiM(good_electrons[0]->pt, good_electrons[0]->eta, good_electrons[0]->phi, 511e-6);
        muonP4.SetPtEtaPhiM(good_muons[0]->pt, good_muons[0]->eta, good_muons[0]->phi, 0.10566);

        leptonOneP4     = electronP4;
        leptonOneFlavor = good_electrons[0]->pdgId;
        leptonTwoP4     = muonP4;
        leptonTwoFlavor = good_muons[0]->pdgId;


        // test if selected lepton fire the trigger.
        eventCounts[channel]->Fill(3);

        // correct for MC, including reconstruction and trigger
        eventCounts[channel]->Fill(4);

    } else {
        return kTRUE;
    }


    outFile->cd(channel.c_str());
    outTrees[channel]->Fill();
    this->passedEvents++;
    return kTRUE;
}

void GenTauTauAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void GenTauTauAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void GenTauTauAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<GenTauTauAnalyzer> selector(new GenTauTauAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


