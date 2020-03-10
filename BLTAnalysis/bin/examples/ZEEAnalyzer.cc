#include "ZEEAnalyzer.h"
#include <map>

// ZEEAnalyzer /eos/uscms/store/group/lpcbacon/12/Summer16_DYJetsToLL_M-50_amcatnlo/Summer16_DYJetsToLL_M-50_amcatnlo_bacon_00.root 10000 dy dy single_electron 2016 1

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;



ZEEAnalyzer::ZEEAnalyzer() : BLTSelector()
{

}

ZEEAnalyzer::~ZEEAnalyzer()
{

}

void ZEEAnalyzer::Begin(TTree *tree)
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

    if (params->selection == "single_electron") {
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

    vector<std::string> channelNames = {"ee"};

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
        
        // weights and their uncertainties
        tree->Branch("eventWeight", &eventWeight);
        tree->Branch("genWeight", &genWeight);


        // vectors
        tree->Branch("leptonOnePt",     &leptonOnePt);
        tree->Branch("leptonTwoPt",     &leptonTwoPt);
        tree->Branch("leptonOneEta",    &leptonOneEta);
        tree->Branch("leptonTwoEta",    &leptonTwoEta);
        tree->Branch("leptonOnePhi",    &leptonOnePhi);
        tree->Branch("leptonTwoPhi",    &leptonTwoPhi);

        tree->Branch("leptonOnePassTriggerTest", &leptonOnePassTriggerTest);
        tree->Branch("leptonTwoPassTriggerTest", &leptonTwoPassTriggerTest);
        tree->Branch("leptonOneNotInGap", &leptonOneNotInGap);
        tree->Branch("leptonTwoNotInGap", &leptonTwoNotInGap);
        

        tree->Branch("dileptonPt",      &dileptonPt);
        tree->Branch("dileptonEta",     &dileptonEta);
        tree->Branch("dileptonPhi",     &dileptonPhi);
        tree->Branch("dileptonMass",    &dileptonMass);

      
        // object counters
        tree->Branch("nElectrons", &nElectrons);

        // met
        tree->Branch("met", &met);
        tree->Branch("metPhi", &metPhi);


        outTrees[channel] = tree;
        // event counter
        string outHistName = params->get_output_treename("TotalEvents_" + channel);
        eventCounts[channel] = new TH1D(outHistName.c_str(),"ChannelCounts",10,0.5,10.5);

    }
    cout<< "Done with initialization"<<endl;
    ReportPostBegin();

}

Bool_t ZEEAnalyzer::Process(Long64_t entry)
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
        unsigned partonCount = 0;
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
        }
        // 3. finish loop over gen particles. now conclude gen level infomation
        //    3.a) counting partons: save to nPartons
        nPartons = partonCount; 
    } else {
        genWeight = 1.0;
        nPU = 0;
        nPartons = 0;
    }
    // cout<< "Gen"<<endl;


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
    // cout<< "Vertices"<<endl;


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
    // cout<< "Trigger"<<endl;

    



    /* -------- ELECTRONS ---------*/
    vector<TElectron*> electrons;

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
                && particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
           ) {
            electrons.push_back(electron);
        }
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);
    nElectrons = electrons.size();
    // cout<< "ELECTRONS"<<endl;





    /* -------- MET ---------*/
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;


    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    

    hTotalEvents->Fill(5);
    string channel = "";
    if (nElectrons == 2 ) { // ee selection

        channel = "ee";
        eventCounts[channel]->Fill(1);

        // test if selected lepton fire the trigger.
        leptonOnePassTriggerTest = false;
        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, electrons[0]->hltMatchBits)) {
                leptonOnePassTriggerTest = true;
                break;
            }
        }
        leptonTwoPassTriggerTest = false;
        for (const auto& name: passTriggerNames) {
            if (trigger->passObj(name, 1, electrons[1]->hltMatchBits)) {
                leptonTwoPassTriggerTest = true;
                break;
            }
        }

        // must pass trigger
        if (!passTrigger)
            return kTRUE;
        eventCounts[channel]->Fill(2);


        // must have opposite signs
        if (electrons[0]->q * electrons[1]->q > 0 )
            return kTRUE;
        eventCounts[channel]->Fill(3);


        TLorentzVector electronOneP4, electronTwoP4, dielectronP4;
        electronOneP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
        electronTwoP4.SetPtEtaPhiM(electrons[1]->pt, electrons[1]->eta, electrons[1]->phi, 511e-6);
        dielectronP4 = electronOneP4 + electronTwoP4;

        leptonOnePt     = electrons[0]->pt;
        leptonOneEta    = electrons[0]->eta;
        leptonOnePhi    = electrons[0]->phi;
        leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet);
        leptonOneFlavor = 11*electrons[0]->q;
        leptonOneNotInGap = abs(leptonOneEta)<1.4442 || abs(leptonOneEta)>1.5660;
        


        leptonTwoPt     = electrons[1]->pt;
        leptonTwoEta    = electrons[1]->eta;
        leptonTwoPhi    = electrons[1]->phi;
        leptonTwoIso    = GetElectronIsolation(electrons[1], fInfo->rhoJet);
        leptonTwoFlavor = 11*electrons[1]->q;
        leptonTwoNotInGap = abs(leptonTwoEta)<1.4442 || abs(leptonTwoEta)>1.5660;


        dileptonPt = dielectronP4.Pt();
        dileptonEta = dielectronP4.Eta();
        dileptonPhi = dielectronP4.Phi();
        dileptonMass = dielectronP4.M();

        
        // must in Z window
        if (dileptonMass < 60 || dileptonMass > 120)
            return kTRUE;
        eventCounts[channel]->Fill(4);


        // must have >=1 tag
        // leptonOneIsTag  = int(leptonOnePt>35 && leptonOneNotInGap && leptonOnePassTriggerTest);
        // leptonTwoIsTag  = int(leptonTwoPt>35 && leptonTwoNotInGap && leptonTwoPassTriggerTest);
        // if ( (!leptonOneIsTag) && (!leptonTwoIsTag) )
        //     return kTRUE;
        // eventCounts[channel]->Fill(5);


        // correct for MC, including reconstruction and trigger
        if (!isData) {

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

        }
    } else {
        return kTRUE;
    }

    outFile->cd(channel.c_str());
    outTrees[channel]->Fill();
    this->passedEvents++;
    return kTRUE;
}

void ZEEAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void ZEEAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void ZEEAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<ZEEAnalyzer> selector(new ZEEAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


float ZEEAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    return combIso;
}

float ZEEAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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

