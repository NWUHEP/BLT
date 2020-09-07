#include "SinglelepAnalyzer.h"
#include <map>


//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;



SinglelepAnalyzer::SinglelepAnalyzer() : BLTSelector()
{

}

SinglelepAnalyzer::~SinglelepAnalyzer()
{

}

void SinglelepAnalyzer::Begin(TTree *tree)
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
    singleMuonTriggerName = "HLT_IsoMu24_v*";
    singleMuonTriggerName2 = "HLT_IsoTkMu24_v*";

    triggerNames.push_back(singleElectronTriggerName);
    triggerNames.push_back(singleMuonTriggerName);
    triggerNames.push_back(singleMuonTriggerName2);
    

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

    vector<std::string> channelNames = {"mu", "e"};

    for (unsigned i = 0; i < channelNames.size(); ++i) {
        string channel = channelNames[i];
        outFile->mkdir(channel.c_str());
        outFile->cd(channel.c_str());
        string treeName = "bltTree_" + params->datasetgroup;
        tree = new TTree(treeName.c_str(), treeName.c_str());

        // event data
        // tree->Branch("runNumber", &runNumber);
        // tree->Branch("evtNumber", &evtNumber, "eventNumber/l");
        // tree->Branch("lumiSection", &lumiSection);
        // tree->Branch("nPV", &nPV);

        // gen level objects
        // // tree->Branch("genWeight", &genWeight);
        // tree->Branch("nPU", &nPU);
        // tree->Branch("nPartons", &nPartons);

        
        // weights and their uncertainties
        tree->Branch("eventWeight", &eventWeight);

        // leptons
        tree->Branch("leptonOneMetMt" , &leptonOneMetMt);
        tree->Branch("leptonOnePt" , &leptonOnePt);
        tree->Branch("leptonOneEta" , &leptonOneEta);
        tree->Branch("leptonOnePhi" , &leptonOnePhi);
        tree->Branch("leptonOneIso", &leptonOneIso);
        tree->Branch("leptonOneIsoPass", &leptonOneIsoPass);

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
        tree->Branch("nJets", &nJets);
        tree->Branch("nBJets", &nBJets);


        outTrees[channel] = tree;
        // event counter
        string outHistName = params->get_output_treename("TotalEvents_" + channel);
        eventCounts[channel] = new TH1D(outHistName.c_str(),"ChannelCounts",10,0.5,10.5);
    }
    cout<< "Done with initialization"<<endl;
    ReportPostBegin();

}

Bool_t SinglelepAnalyzer::Process(Long64_t entry)
{
    /* beginning of event loop */
    GetEntry(entry, 1);  // load all branches
    outFile->cd();
    eventWeight = 1.;
    this->totalEvents++;
    hTotalEvents->Fill(1);

    const bool isData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isData);
    if (!isData) {
        genWeight = fGenEvtInfo->weight > 0 ? 1 : -1;
        if (genWeight < 0) {
            hTotalEvents->Fill(10);
        }
        eventWeight *= genWeight;
    }
    
    if (entry%10000==0) {
        std::cout << "... Processing event " << entry 
            << " Run: " << fInfo->runNum 
            << " Lumi: " << fInfo->lumiSec 
            << " Event: " << fInfo->evtNum 
            << std::endl;
    }





     /* -------- Apply lumi mask  ---------*/
    if (isData) {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(2);

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
    bool muonTriggered = find(passTriggerNames.begin(), passTriggerNames.end(), singleMuonTriggerName) != passTriggerNames.end();
    muonTriggered = muonTriggered || find(passTriggerNames.begin(), passTriggerNames.end(), singleMuonTriggerName2) != passTriggerNames.end();
    bool electronTriggered = find(passTriggerNames.begin(), passTriggerNames.end(), singleElectronTriggerName) != passTriggerNames.end();

    if (!passTrigger)
        return kTRUE;
    hTotalEvents->Fill(4);



    // Set data period for 2016 MC scale factors
    if (!isData) {
        if (rng->Rndm() < 0.548) {  //0.468) {
            weights->SetDataPeriod("2016BtoF"); 
        } else {
            weights->SetDataPeriod("2016GH");
        }
    }

    /////////////////////
    // Fill event info //
    /////////////////////

    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
    // nPV           = fPVArr->GetEntries();
    
    
    ///////////////////////
    // Generator objects //
    ///////////////////////

    vector<TGenParticle*> genElectrons, genMuons;

    if (!isData) {


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

        }
        // cout << endl;

        

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
                && (zPt > 0. || zP4.Pt() > 0.)
                ) {
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
 

    } else {
        genWeight = 1.0;
        nPU = 0;
        nPartons = 0;
    }

    // cout << nPU << endl;
    ///////////////////
    // Select objects//
    ///////////////////




    /* -------- MUONS ---------*/
    vector<TMuon*> muons;
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
                && (GetMuonIsolation(muon)/muonP4.Pt() <  0.25)
           ) {
                muons.push_back(muon);
                veto_muons.push_back(muonP4);
        } 
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);
    nMuons = muons.size();

    // fill hist if muonTriggered and leading muon is good.
    if (muonTriggered && nMuons>0){
        if (muons[0]->pt>25){
            hTotalEvents->Fill(5);
        }
    }
    


    /* -------- ELECTRONS ---------*/
    vector<TElectron*> electrons;
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

                electrons.push_back(electron);
                veto_electrons.push_back(electronP4);
        }
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);

    nElectrons = electrons.size();

    // fill hist if electronTriggered and leading nElectrons is good.
    if (electronTriggered && nElectrons>0){
        if (electrons[0]->pt>30){
            hTotalEvents->Fill(6);
        }
    }

    
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

    std::vector<TJet*> jets;
    vector<TLorentzVector> prefiring_jets;
    TLorentzVector hadronicP4;
    float sumJetPt = 0;

    // set nJet and nBJet counters to 0
    nJets = nBJets = 0;
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


        // save good jets
        if (
                fabs(jet->eta) <= 2.4
                && (jet->pt > 30)
                && particleSelector->PassJetID(jet, cuts->looseJetID)
                && !muOverlap 
                && !elOverlap
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

    hTotalEvents->Fill(7);
    string channel = "";
    if (nElectrons == 1 && nMuons == 0) { // e selection

        channel = "e";
        eventCounts[channel]->Fill(1);

        // pass pt threshold
        if (electrons[0]->pt < 30)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        // triggerred
        if (!electronTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        // nJets
        if (nJets<1)
            return kTRUE;
        eventCounts[channel]->Fill(4);




        TLorentzVector electronP4;
        electronP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, 511e-6);
	    // mt
        leptonOneMetMt = sqrt(2*electronP4.Pt()*met*(1-cos(electronP4.Phi() - metPhi)));
 
        if (nJets<=3){
           if (leptonOneMetMt>40)
                return kTRUE;
            eventCounts[channel]->Fill(6);
        }



        // save other lepton variables
        leptonOnePt     = electronP4.Pt();
        leptonOneEta    = electronP4.Eta();
        leptonOnePhi    = electronP4.Phi();
        leptonOneIso    = GetElectronIsolation(electrons[0], fInfo->rhoJet)/electronP4.Pt();
        leptonOneIsoPass= particleSelector->PassElectronIso(electrons[0], cuts->tightElIso, cuts->EAEl);

        // correct for MC, including reconstruction and trigger
        if (!isData) {


            // correct for id,reco weights
            EfficiencyContainer effCont;
            pair<float, float> effs, errs;

            // id weight
            effCont = weights->GetElectronIDEff(electronP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonOneIDWeight = effs.first/effs.second;
            leptonOneIDVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            eventWeight *= leptonOneIDWeight;

            // reconstruction weight
            effCont = weights->GetElectronRecoEff(electronP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonOneRecoWeight = effs.first/effs.second;
            leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            eventWeight *= leptonOneRecoWeight;

            // correct for trigger.
            EfficiencyContainer effCont1;
            effCont1 = weights->GetTriggerEffWeight(singleElectronTriggerName, electronP4);
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            triggerWeight = effs.first/effs.second;
            triggerVar    = pow(triggerWeight, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            eventWeight *= triggerWeight;

            // prefiring
            effCont = weights->GetElectronPrefiringWeight(prefiring_photons, prefiring_jets);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            prefiringWeight = effs.first/effs.second;
            prefiringVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            eventWeight *= prefiringWeight;
            
        }

    } else if (nMuons == 1 && nElectrons == 0) { // mu selection

        channel = "mu";
        eventCounts[channel]->Fill(1);
        // coutchannel

        // pass pt threshold
        if (muons[0]->pt < 25)
            return kTRUE;
        eventCounts[channel]->Fill(2);

        // triggerred
        if (!muonTriggered)
            return kTRUE;
        eventCounts[channel]->Fill(3);

        // nJets
        if (nJets<1)
            return kTRUE;
        eventCounts[channel]->Fill(4);


        TLorentzVector muonP4;
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, 0.10566);
        // mt cut
        leptonOneMetMt = sqrt(2*muonP4.Pt()*met*(1-cos(muonP4.Phi() - metPhi)));

        if (nJets<=3){
           if (leptonOneMetMt>40)
                return kTRUE;
            eventCounts[channel]->Fill(6);
        }

        // save other lepton variables
        leptonOnePt     = muonP4.Pt();
        leptonOneEta    = muonP4.Eta();
        leptonOnePhi    = muonP4.Phi();
        leptonOneIso    = GetMuonIsolation(muons[0])/muonP4.Pt();
        leptonOneIsoPass=(GetMuonIsolation(muons[0])/muonP4.Pt()<0.15);


        // correct for MC, including reconstruction and trigger
        if (!isData) {

            // correct for id,reco weights
            EfficiencyContainer effCont;
            pair<float, float> effs, errs;

            // id weight
            effCont = weights->GetMuonIDEff(muonP4);
            effs = effCont.GetEff();
            errs = effCont.GetErr();
            leptonOneIDWeight = effs.first/effs.second;
            leptonOneIDVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            eventWeight *= leptonOneIDWeight;

            // iso weight
            if (leptonOneIsoPass){
                effCont = weights->GetMuonISOEff(muonP4);
                effs = effCont.GetEff();
                errs = effCont.GetErr();
                leptonOneRecoWeight = effs.first/effs.second;
                leptonOneRecoVar    = pow(effs.first/effs.second, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
                eventWeight *= leptonOneRecoWeight;
            }



            // correct for trigger.
            EfficiencyContainer effCont1;
            effCont1 = weights->GetTriggerEffWeight(singleMuonTriggerName, muonP4);
            effs = effCont1.GetEff();
            errs = effCont1.GetErr();
            triggerWeight = effs.first/effs.second;
            triggerVar    = pow(triggerWeight, 2)*(pow(errs.first/effs.first, 2) + pow(errs.second/effs.second, 2));
            eventWeight *= triggerWeight;
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

void SinglelepAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void SinglelepAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void SinglelepAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<SinglelepAnalyzer> selector(new SinglelepAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


float SinglelepAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    float combIso = (mu->chHadIso + std::max(0.,(double)mu->neuHadIso + mu->gammaIso - 0.5*mu->puIso));
    return combIso;
}

float SinglelepAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, const float rho)
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
