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

hzgAnalyzer::hzgAnalyzer() : BLTSelector()
{

}

hzgAnalyzer::~hzgAnalyzer()
{

}

void hzgAnalyzer::Begin(TTree *tree)
{
    rng = new TRandom3();

    // Parse command line option
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1),
            std::sregex_token_iterator(), std::back_inserter(options));

    // Set the parameters
    params.reset(new Parameters());
    params->setup(options);

    particleSelector.reset(new ParticleSelector(*params));

    // Trigger bits mapping file
    const std::string cmssw_base = getenv("CMSSW_BASE");
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns";
    trigger.reset(new baconhep::TTrigger(trigfilename));

    if (params->selection == "mm" || params->selection == "mmg") {
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

    else if (params->selection == "ee" || params->selection == "eeg") {
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
    weights.reset(new WeightUtils(params->period, params->selection)); // Lumi mask

    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    string jsonFileName;
    if (params->period == "2016") {
        jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/json/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt";
    }
    else if (params->period == "2017") {
        jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/json/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt";
    }
    else if (params->period == "2018") {
        jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/json/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt";
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
    outTree->Branch("nPV", &nPV);
    outTree->Branch("nPU", &nPU);
    outTree->Branch("nPartons", &nPartons);
    outTree->Branch("xPV", &xPV);
    outTree->Branch("yPV", &yPV);
    outTree->Branch("zPV", &zPV);
    outTree->Branch("rho", &rho);
    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);
    outTree->Branch("ht", &ht);
    outTree->Branch("htPhi", &htPhi);
    outTree->Branch("htSum", &htSum);
   
    // category tags
    outTree->Branch("isLeptonTag", &isLeptonTag);
    outTree->Branch("isDijetTag", &isDijetTag);
    outTree->Branch("isTightDijetTag", &isTightDijetTag);

    // weights
    outTree->Branch("genWeight", &genWeight);
    outTree->Branch("eventWeight", &eventWeight);
    outTree->Branch("puWeight", &puWeight);
    outTree->Branch("puWeightUp", &puWeightUp);
    outTree->Branch("puWeightDown", &puWeightDown); 
    outTree->Branch("leptonOneIDWeight", &leptonOneIDWeight);
    outTree->Branch("leptonTwoIDWeight", &leptonTwoIDWeight);
    outTree->Branch("photonIDWeight", &photonIDWeight);
    outTree->Branch("trigOneWeight", &trigOneWeight);
    outTree->Branch("trigTwoWeight", &trigTwoWeight);
    outTree->Branch("triggerWeight", &triggerWeight);
    outTree->Branch("prefWeight", &prefWeight);
    outTree->Branch("prefWeightUp", &prefWeightUp);
    outTree->Branch("prefWeightDown", &prefWeightDown);
    outTree->Branch("photonR9Weight", &photonR9Weight);

    // leptons
    outTree->Branch("leptonOnePt", &leptonOnePt);
    outTree->Branch("leptonTwoPt", &leptonTwoPt);
    outTree->Branch("leptonOneEta", &leptonOneEta);
    outTree->Branch("leptonTwoEta", &leptonTwoEta);
    outTree->Branch("leptonOnePhi", &leptonOnePhi);
    outTree->Branch("leptonTwoPhi", &leptonTwoPhi);
    outTree->Branch("leptonOnePtKin", &leptonOnePtKin);
    outTree->Branch("leptonTwoPtKin", &leptonTwoPtKin);
    outTree->Branch("leptonOnePtKinErr", &leptonOnePtKinErr);
    outTree->Branch("leptonTwoPtKinErr", &leptonTwoPtKinErr);
    outTree->Branch("leptonOneIso", &leptonOneIso);
    outTree->Branch("leptonTwoIso", &leptonTwoIso);
    outTree->Branch("leptonOneFlavor", &leptonOneFlavor);
    outTree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
    outTree->Branch("leptonOneD0", &leptonOneD0);
    outTree->Branch("leptonTwoD0", &leptonTwoD0);
    outTree->Branch("leptonOneDZ", &leptonOneDZ);
    outTree->Branch("leptonTwoDZ", &leptonTwoDZ);

    if (params->selection == "tautaug") {        
        outTree->Branch("tauDecayMode", &tauDecayMode);
        outTree->Branch("tauMVA", &tauMVA);
    }

    // photons
    outTree->Branch("photonPt", &photonPt);
    outTree->Branch("photonEta", &photonEta);
    outTree->Branch("photonPhi", &photonPhi);
    outTree->Branch("photonRawPt", &photonRawPt);
    outTree->Branch("photonSCPt", &photonSCPt);
    outTree->Branch("photonSCEta", &photonSCEta);
    outTree->Branch("photonSCPhi", &photonSCPhi);
    outTree->Branch("photonR9", &photonR9);
    outTree->Branch("photonRawR9", &photonRawR9);
    outTree->Branch("photonMVA", &photonMVA);
    outTree->Branch("photonERes", &photonERes);
    outTree->Branch("photonE", &photonE);
    outTree->Branch("photonErrE", &photonErrE);
    //outTree->Branch("photonWorstChIso", &photonWorstChIso); // on hold
    outTree->Branch("passElectronVeto", &passElectronVeto);

    // jets
    outTree->Branch("jetOnePt", &jetOnePt);
    outTree->Branch("jetTwoPt", &jetTwoPt);
    outTree->Branch("jetOneRawPt", &jetOneRawPt);
    outTree->Branch("jetTwoRawPt", &jetTwoRawPt);
    outTree->Branch("jetOneEta", &jetOneEta);
    outTree->Branch("jetTwoEta", &jetTwoEta);
    outTree->Branch("jetOnePhi", &jetOnePhi);
    outTree->Branch("jetTwoPhi", &jetTwoPhi);
    outTree->Branch("jetOneTag", &jetOneTag);
    outTree->Branch("jetTwoTag", &jetTwoTag);
    outTree->Branch("jetOneArea", &jetOneArea);
    outTree->Branch("jetTwoArea", &jetTwoArea);
    outTree->Branch("jetOneL1Corr", &jetOneL1Corr);
    outTree->Branch("jetTwoL1Corr", &jetTwoL1Corr);
    outTree->Branch("jetOneL2Corr", &jetOneL2Corr);
    outTree->Branch("jetTwoL2Corr", &jetTwoL2Corr);
    outTree->Branch("jetOneL3Corr", &jetOneL3Corr);
    outTree->Branch("jetTwoL3Corr", &jetTwoL3Corr);
    outTree->Branch("jetOneL4Corr", &jetOneL4Corr);
    outTree->Branch("jetTwoL4Corr", &jetTwoL4Corr);

    // gen level objects 
    outTree->Branch("genLeptonOneId", &genLeptonOneId);
    outTree->Branch("genLeptonTwoId", &genLeptonTwoId);
    outTree->Branch("genLeptonOnePt", &genLeptonOnePt);
    outTree->Branch("genLeptonTwoPt", &genLeptonTwoPt);
    outTree->Branch("genLeptonOneEta", &genLeptonOneEta);
    outTree->Branch("genLeptonTwoEta", &genLeptonTwoEta);
    outTree->Branch("genLeptonOnePhi", &genLeptonOnePhi);
    outTree->Branch("genLeptonTwoPhi", &genLeptonTwoPhi);
    outTree->Branch("genPhotonPt", &genPhotonPt);
    outTree->Branch("genPhotonEta", &genPhotonEta);
    outTree->Branch("genPhotonPhi", &genPhotonPhi);
    outTree->Branch("vetoDY", &vetoDY);

    // multiplicities
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nGenMuons", &nGenMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nGenElectrons", &nGenElectrons);
    outTree->Branch("nTaus", &nTaus);
    outTree->Branch("nGenTaus", &nGenTaus);
    outTree->Branch("nPhotons", &nPhotons);
    outTree->Branch("nGenPhotons", &nGenPhotons);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nFwdJets", &nFwdJets);
    outTree->Branch("nCentralJets", &nCentralJets);
    outTree->Branch("nBJets", &nBJets);
    
    // dilepton
    outTree->Branch("dileptonPt", &dileptonPt);
    outTree->Branch("dileptonEta", &dileptonEta);
    outTree->Branch("dileptonPhi", &dileptonPhi);
    outTree->Branch("dileptonM", &dileptonM);
    outTree->Branch("dileptonMKin", &dileptonMKin);
    outTree->Branch("dileptonDEta", &dileptonDEta);
    outTree->Branch("dileptonDPhi", &dileptonDPhi);
    outTree->Branch("dileptonDR", &dileptonDR);
 
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
    outTree->Branch("llgMKin", &llgMKin);
    outTree->Branch("llgPtOverM", &llgPtOverM);
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

    // angles
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
    outTree->Branch("photonZepp", &photonZepp);
    outTree->Branch("vbfPtBalance", &vbfPtBalance);
    outTree->Branch("jetOneMatched", &jetOneMatched);
    outTree->Branch("jetTwoMatched", &jetTwoMatched);
    outTree->Branch("leptonOneMatched", &leptonOneMatched);
    outTree->Branch("leptonTwoMatched", &leptonTwoMatched);
    outTree->Branch("photonMatched", &photonMatched);

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

    bool sync = false;
    if (sync) {
        if (!(
                (fInfo->runNum == 1 && fInfo->lumiSec == 137 && fInfo->evtNum == 26380) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 142 && fInfo->evtNum == 27518) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 107 && fInfo->evtNum == 20615) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 109 && fInfo->evtNum == 21031) || 
                (fInfo->runNum == 1 && fInfo->lumiSec == 106 && fInfo->evtNum == 20535) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 148 && fInfo->evtNum == 28598) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 148 && fInfo->evtNum == 28660) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 151 && fInfo->evtNum == 29133) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 130 && fInfo->evtNum == 25077) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 124 && fInfo->evtNum == 23877) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 120 && fInfo->evtNum == 23152) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 126 && fInfo->evtNum == 24285) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 127 && fInfo->evtNum == 24597) || 
                (fInfo->runNum == 1 && fInfo->lumiSec == 128 && fInfo->evtNum == 24683) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 128 && fInfo->evtNum == 24845) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 102 && fInfo->evtNum == 19797) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 80  && fInfo->evtNum == 15364) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 72  && fInfo->evtNum == 13886) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 74  && fInfo->evtNum == 14352) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 69  && fInfo->evtNum == 13193) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 69  && fInfo->evtNum == 13209) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 76  && fInfo->evtNum == 14769) || 
                (fInfo->runNum == 1 && fInfo->lumiSec == 91  && fInfo->evtNum == 17510) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 91  && fInfo->evtNum == 17513) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 91  && fInfo->evtNum == 17544) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 91  && fInfo->evtNum == 17563) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 91  && fInfo->evtNum == 17568) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 67  && fInfo->evtNum == 12820) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 67  && fInfo->evtNum == 12945) ||
                (fInfo->runNum == 1 && fInfo->lumiSec == 68  && fInfo->evtNum == 13127)
            )) 
        {
            return kTRUE;
        }
    }

    if (sync) {
        std::cout << "run, lumi, evt = " << fInfo->runNum << ", " << fInfo->lumiSec << ", " << fInfo->evtNum << std::endl;
    }

          
    ///////////////////////
    // Generator objects // 
    ///////////////////////

    vector<TGenParticle*> genLeptons;
    vector<TGenParticle*> genPhotons;
    vector<TGenParticle*> genFSPartons;
    if (!isData) {
        unsigned count = 0;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            if (
                    particle->status == 23 
                    && (fabs(particle->pdgId) < 6 || particle->pdgId == 21) 
                    && particle->parent != -2
               ) {
                ++count;
            }

            if (    (fabs(particle->pdgId) == 11 || fabs(particle->pdgId) == 13) 
                    && particle->fromHardProcessFinalState 
                ) {
                    genLeptons.push_back(particle);
            }
                
            // saving photons     
            if (fabs(particle->pdgId) == 22) {
                genPhotons.push_back(particle);    
            }
            
            if(     particle->isHardProcess 
                    && (fabs(particle->pdgId) <= 6 || particle->pdgId == 21)
               ) {
                genFSPartons.push_back(particle);
            }

            }

            nPartons = count; // This is saved for reweighting inclusive DY and combining it with parton binned DY


    } else {
        nPartons = 0;
    }

    nGenPhotons = genPhotons.size();

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

            // remove overlap between electron and muon channels
            if (params->selection == "eeg") {
                if (params->period == "2016") {
                    if (
                            triggerNames[i] == "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*" || 
                            triggerNames[i] == "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*" ||
                            triggerNames[i] == "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*" || 
                            triggerNames[i] == "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*"
                        ) {
                        passTrigger = false;
                        break;
                    }
                }
                else if (params->period == "2017") {
                    if (
                            triggerNames[i] == "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*" || 
                            triggerNames[i] == "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*"
                        ) {
                        passTrigger = false;
                        break;
                    }
                }
                else if (params->period == "2018") {
                    if (triggerNames[i] == "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*" ) {
                        passTrigger = false;
                        break;
                    }
                }
            }

        }
    }

    if (!passTrigger) // && isData)
        return kTRUE;
    hTotalEvents->Fill(3);

    /////////////////////
    // Fill event info //
    /////////////////////

    eventWeight   = 1;
    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
    nPV           = fPVArr->GetEntries();
    rho           = fInfo->rhoJet;

    if (!isData) {
        nPU = fInfo->nPUmean;
        std::map<std::string, float> puWeights = weights->GetPUWeight(nPU);
        puWeight = puWeights["nom"]; // pileup reweighting
        puWeightUp = puWeights["up"]; // pileup reweighting
        puWeightDown = puWeights["down"]; // pileup reweighting
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
    } else {
        return kTRUE;
    }
    hTotalEvents->Fill(4);

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

        if (sync) {
            std::cout << "electron calibPt, scEta, reliso, dz, d0" << std::endl;
            std::cout << electron->calibPt << ", " << electron->scEta << ", "
                      << particleSelector->GetElectronIsolation(electron, fInfo->rhoJet)/electron->calibPt << ", " 
                      << electron->d0 << ", " << electron->dz << std::endl;
        }
    
        if (
                electron->calibPt > 7.
                && fabs(electron->scEta) < 2.5
                && particleSelector->PassElectronMVA(electron, "looseFall17V2")
                && fabs(electron->d0) < 0.5
                && fabs(electron->dz) < 1.0
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
        if (!isData) particleSelector->ApplyTauEnergyScaleCorrection(tau);

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

        // 2017 EE noise veto
        if (params->period == "2017") {
            if (
                jet->ptRaw < 50. 
                && fabs(jet->eta) > 2.65 
                && fabs(jet->eta) < 3.139
                ) {
                continue;
            }
        }

        double jec = particleSelector->JetCorrector(jet, "NONE");
        jet->pt = jet->ptRaw*jec;

        float gRand = 1.;
        float jerc = 1.;
        bool hasGenJetMatch = false;
        if (!isData) { // apply jet energy resolution corrections to simulation
            pair<float, float> resPair = particleSelector->JetResolutionAndSF(jet, 0);

            if (!(jet->genpt == 0 && jet->geneta == 0 && jet->genphi == 0 && jet->genm == 0)) {
                TLorentzVector tmpJetP4, tmpGenJetP4;
                tmpJetP4.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
                tmpGenJetP4.SetPtEtaPhiM(jet->genpt, jet->geneta, jet->genphi, jet->genm);
                if (tmpJetP4.DeltaR(tmpGenJetP4) < 0.2 && fabs(jet->pt - jet->genpt) < 3*resPair.first*jet->pt) {
                    hasGenJetMatch = true;
                }
            }
                


            if (hasGenJetMatch) {
                jerc = 1 + (resPair.second - 1)*((jet->pt - jet->genpt)/jet->pt);
            }
            else {
                gRand = rng->Gaus(0, resPair.first);
                jerc = 1 + gRand*sqrt(std::max((double)resPair.second*resPair.second - 1, 0.));
            }

            //jerc = std::max((double)jerc, 0.); // truncate jerc at 0
            jet->pt = jet->pt*jerc;
        }

        TLorentzVector jetP4; 
        jetP4.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
        
        // Prevent overlap of other objects and jets

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
            jetIDName = "loose";
        }
        else { // 2017 or 2018
            jetIDName = "tight";
        }

        if (sync) {
            std::cout << "jet passes id = " << particleSelector->PassJetID(jet, jetIDName) << std::endl;
            std::cout << "jec, jerc = " << jec << ", " << jerc << std::endl;
            std::cout << "raw pt, corr pt, eta, phi, area, rho " << 
                          jet->ptRaw << ", " << jet->pt << ", " 
                          << jet->eta << ", " << jet->phi << ", " 
                          << jet->area << ", " << fInfo->rhoJet << std::endl;
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
        
    if (sync) {
        std::cout << "only preselected jets below" << std::endl;
        for (unsigned int i = 0; i < jets.size(); ++i) {
            std::cout << "raw pt, corr pt, eta, phi, area, rho " << 
                          jets.at(i)->ptRaw << ", " << jets.at(i)->pt << ", " 
                          << jets.at(i)->eta << ", " << jets.at(i)->phi << ", " 
                          << jets.at(i)->area << ", " << fInfo->rhoJet << std::endl;
        }
    }

    //std::sort(jets.begin(), jets.end(), sort_by_btag);

    /* MET */
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;

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

    if (params->selection == "mm") {
        
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
            eventWeight *= weights->GetElectronMVARecoIdEff(*electrons[0]); 
            eventWeight *= weights->GetElectronMVARecoIdEff(*electrons[1]); 
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
        TLorentzVector photonP4;
        photonP4.SetPtEtaPhiM(photons[0]->calibPt, photons[0]->eta, photons[0]->phi, 0.);
        photonPt = photonP4.Pt();
        photonEta = photonP4.Eta();
        photonPhi = photonP4.Phi();
        photonMVA = photons[0]->mvaFall17V2;
        passElectronVeto = photons[0]->passElectronVeto;  
        if (!isData)
            photonR9 = weights->GetCorrectedPhotonR9(*photons[0]);
        else 
            photonR9 = photons[0]->r9_full5x5;

        // DY photon overlap removal
        vetoDY = false;
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
            
        }

    } // end tautaug selection

    else { // llg selection (combines mmg and eeg)
        assert(params->selection == "mmg" || params->selection == "eeg");

        unsigned int nLeptons = 0;
        if (params->selection == "mmg") 
            nLeptons = muons.size();
        else if (params->selection == "eeg")
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
                if (params->selection == "mmg") {
                    tempLeptonOne.SetPtEtaPhiM(muons[i]->pt, muons[i]->eta, muons[i]->phi, MUON_MASS);
                    tempLeptonTwo.SetPtEtaPhiM(muons[j]->pt, muons[j]->eta, muons[j]->phi, MUON_MASS);
                }
                else if (params->selection == "eeg") {
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

        if (params->selection == "mmg") {
            if (leptonOneP4.Pt() <= 20.0) 
                return kTRUE;
            if (leptonTwoP4.Pt() <= 10.0)
                return kTRUE;
        }
        else if (params->selection == "eeg") {  
            if (leptonOneP4.Pt() <= 25.0)
                return kTRUE;
            if (leptonTwoP4.Pt() <= 15.0)
                return kTRUE;
        }

        if (sync) {
            std::cout << "lepton one pt, eta, phi: " << leptonOneP4.Pt() << ", " << leptonOneP4.Eta() << ", " << leptonOneP4.Phi() << std::endl;
            std::cout << "lepton two pt, eta, phi: " << leptonTwoP4.Pt() << ", " << leptonTwoP4.Eta() << ", " << leptonTwoP4.Phi() << std::endl;
        }

        TLorentzVector dileptonP4 = leptonOneP4 + leptonTwoP4;

        // L1EMTF cut 
        if (params->selection == "mmg") {
            if (
                fabs(leptonOneP4.DeltaPhi(leptonTwoP4)) < 70.0*(M_PI/180.0)
                && fabs(leptonOneP4.Eta()) > 1.2 
                && fabs(leptonTwoP4.Eta()) > 1.2
                && leptonOneP4.Eta()*leptonTwoP4.Eta() > 0
               )
                return kTRUE;
        }
        hTotalEvents->Fill(8); 
        
        // checking lepton matching
        leptonOneMatched = false;
        leptonTwoMatched = false;
        if (!isData && genLeptons.size() > 0) {
            for (unsigned int i = 0; i < genLeptons.size(); ++i) {
                TLorentzVector tmpGenLepton;
                tmpGenLepton.SetPtEtaPhiM(genLeptons[i]->pt, genLeptons[i]->eta, genLeptons[i]->phi, genLeptons[i]->mass);
                if (tmpGenLepton.DeltaR(leptonOneP4) < 0.1) leptonOneMatched = true;
                if (tmpGenLepton.DeltaR(leptonTwoP4) < 0.1) leptonTwoMatched = true;
            }
        }

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

        TLorentzVector photonP4;
        photonP4.SetPtEtaPhiM(photons[photonIndex]->calibPt, photons[photonIndex]->eta, photons[photonIndex]->phi, 0.);
        if (photonP4.Pt() < 15.0)
            return kTRUE;
        hTotalEvents->Fill(10);
        
        if (sync) {
            std::cout << "photon pt, eta, phi: " << photonP4.Pt() << ", " << photonP4.Eta() << ", " << photonP4.Phi() << std::endl;
        }

        TLorentzVector llgP4 = dileptonP4 + photonP4;

        // DY photon overlap removal
        vetoDY = false;
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

        // checking photon matching
        photonMatched = false;
        if (!isData && genPhotons.size() > 0) {
            for (unsigned int i = 0; i < genPhotons.size(); ++i) {
                TLorentzVector tmpGenPhoton;
                tmpGenPhoton.SetPtEtaPhiM(genPhotons[i]->pt, genPhotons[i]->eta, genPhotons[i]->phi, genPhotons[i]->mass);
                if (tmpGenPhoton.DeltaR(photonP4) < 0.1) photonMatched = true;
            }
        }
            
        // checking for lepton tag
        isLeptonTag = false;
        for (unsigned int i = 0; i < muons.size(); ++i) {
            TLorentzVector tempMuon;
            tempMuon.SetPtEtaPhiM(muons[i]->pt, muons[i]->eta, muons[i]->phi, MUON_MASS);
            if (leptonOneP4.DeltaR(tempMuon) < 0.4 ||
                leptonTwoP4.DeltaR(tempMuon) < 0.4 || 
                photonP4.DeltaR(tempMuon) < 0.4)
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
                    photonP4.DeltaR(tempElectron) < 0.4)
                    continue;
                else
                    isLeptonTag = true;
            }
        }

        //std::cout << "isLeptonTag = " << isLeptonTag << std::endl;
        
        // checking for dijet tag
        isDijetTag = false;
        isTightDijetTag = false;
        unsigned int jetOneIndex = 0;
        unsigned int jetTwoIndex = 0;
        TLorentzVector jetOneP4, jetTwoP4;
        jetOneMatched = false;
        jetTwoMatched = false;  
        if (!isLeptonTag) {
            if (jets.size() > 1)  {
                for (unsigned int i = 0; i < jets.size(); ++i) {
                    for (unsigned int j = i+1; j < jets.size(); ++j) {
                        TLorentzVector tempJetOne;
                        TLorentzVector tempJetTwo;
                        tempJetOne.SetPtEtaPhiM(jets[i]->pt, jets[i]->eta, jets[i]->phi, jets[i]->mass);
                        tempJetTwo.SetPtEtaPhiM(jets[j]->pt, jets[j]->eta, jets[j]->phi, jets[j]->mass);
                        TLorentzVector tempDijet = tempJetOne + tempJetTwo;
                        if (sync) {
                            std::cout << "jet indices = " << i << ", " << j << std::endl;
                            std::cout << "j1L1DR, j1L2DR, j1GamDR: " << 
                                tempJetOne.DeltaR(leptonOneP4) << ", " << 
                                tempJetOne.DeltaR(leptonTwoP4) << ", " << 
                                tempJetOne.DeltaR(photonP4) << std::endl;
                            std::cout << "j2L1DR, j2L2DR, j2GamDR: " << 
                                tempJetTwo.DeltaR(leptonOneP4) << ", " << 
                                tempJetTwo.DeltaR(leptonTwoP4) << ", " << 
                                tempJetTwo.DeltaR(photonP4) << std::endl;
                        }
                        if  (   tempJetOne.DeltaR(leptonOneP4) >= 0.4 && tempJetTwo.DeltaR(leptonOneP4) >= 0.4 
                            &&  tempJetOne.DeltaR(leptonTwoP4) >= 0.4 && tempJetTwo.DeltaR(leptonTwoP4) >= 0.4 
                            &&  tempJetOne.DeltaR(photonP4) >= 0.4 && tempJetTwo.DeltaR(photonP4) >= 0.4
                            ) {
                            isDijetTag = true;
                            jetOneIndex = i;
                            jetTwoIndex = j;
                            float zeppen = llgP4.Eta() - (tempJetOne.Eta() + tempJetTwo.Eta())/2.;
                            if  (   fabs(tempJetOne.Eta() - tempJetTwo.Eta()) >= 3.5 
                                &&  fabs(zeppen) <= 2.5 && tempDijet.M() >= 500.
                                &&  fabs(tempDijet.DeltaPhi(llgP4)) >= 2.4
                                ) {
                                isTightDijetTag = true;
                            }
                            break;
                        }
                    }
                    if (isDijetTag) break;
                }
            }
        }

        if (sync) {
            std::cout << "isDijetTag = " << isDijetTag << std::endl;
        }

        if (isDijetTag) {
            jetOneP4.SetPtEtaPhiM(jets[jetOneIndex]->pt, jets[jetOneIndex]->eta, jets[jetOneIndex]->phi, jets[jetOneIndex]->mass);
            jetTwoP4.SetPtEtaPhiM(jets[jetTwoIndex]->pt, jets[jetTwoIndex]->eta, jets[jetTwoIndex]->phi, jets[jetTwoIndex]->mass);
            jetOnePt = jetOneP4.Pt();
            jetOneEta = jetOneP4.Eta();
            jetOnePhi = jetOneP4.Phi();
            jetOneM   = jetOneP4.M();
            jetOneRawPt = jets[jetOneIndex]->ptRaw;
            jetOneArea = jets[jetOneIndex]->area;

            std::vector<float> jetOneSubcorrections = particleSelector->GetJetSubcorrections(jets[jetOneIndex], "NONE");
            jetOneL1Corr = jetOneSubcorrections.at(0);
            jetOneL2Corr = jetOneSubcorrections.at(1);
            jetOneL3Corr = jetOneSubcorrections.at(2);
            jetOneL4Corr = jetOneSubcorrections.at(3);
            
            jetTwoPt = jetTwoP4.Pt();
            jetTwoEta = jetTwoP4.Eta();
            jetTwoPhi = jetTwoP4.Phi();
            jetTwoM   = jetTwoP4.M();
            jetTwoRawPt = jets[jetTwoIndex]->ptRaw;
            jetTwoArea = jets[jetTwoIndex]->area;
            
            std::vector<float> jetTwoSubcorrections = particleSelector->GetJetSubcorrections(jets[jetTwoIndex], "NONE");
            jetTwoL1Corr = jetTwoSubcorrections.at(0);
            jetTwoL2Corr = jetTwoSubcorrections.at(1);
            jetTwoL3Corr = jetTwoSubcorrections.at(2);
            jetTwoL4Corr = jetTwoSubcorrections.at(3);

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

            j1PhotonDEta = fabs(jetOneP4.Eta() - photonP4.Eta());
            j1PhotonDPhi = fabs(jetOneP4.DeltaPhi(photonP4));
            j1PhotonDR = jetOneP4.DeltaR(photonP4);
            j2PhotonDEta = fabs(jetTwoP4.Eta() - photonP4.Eta());
            j2PhotonDPhi = fabs(jetTwoP4.DeltaPhi(photonP4));
            j2PhotonDR = jetTwoP4.DeltaR(photonP4);
        
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
            photonZepp = photonP4.Eta() - (jetOneP4.Eta() + jetTwoP4.Eta())/2.;
            vbfPtBalance = (photonP4.Px() + photonP4.Py() + dileptonP4.Px() + dileptonP4.Py() + jetOneP4.Px() + jetOneP4.Py() + jetTwoP4.Px() + jetTwoP4.Py()) / (photonP4.Pt() + dileptonP4.Pt() + jetOneP4.Pt() + jetTwoP4.Pt());

            // jet truth information
            if (!isData && genFSPartons.size() > 0) {
                TLorentzVector genJetOne, genJetTwo;
                genJetOne.SetPtEtaPhiM(jets[jetOneIndex]->genpt, jets[jetOneIndex]->geneta, jets[jetOneIndex]->genphi, jets[jetOneIndex]->genm);
                genJetTwo.SetPtEtaPhiM(jets[jetTwoIndex]->genpt, jets[jetTwoIndex]->geneta, jets[jetTwoIndex]->genphi, jets[jetTwoIndex]->genm);
                for (unsigned int i = 0; i < genFSPartons.size(); i++) {
                    TLorentzVector tmpGenParton;
                    tmpGenParton.SetPtEtaPhiM(genFSPartons[i]->pt, genFSPartons[i]->eta, genFSPartons[i]->phi, genFSPartons[i]->mass);
                    if (tmpGenParton.DeltaR(genJetOne) < 0.1) jetOneMatched = true;
                    if (tmpGenParton.DeltaR(genJetTwo) < 0.1) jetTwoMatched = true;
                }
            }
                
            

        }
     
        leptonOnePt     = leptonOneP4.Pt();
        leptonOneEta     = leptonOneP4.Eta();
        leptonOnePhi     = leptonOneP4.Phi();
        leptonTwoPt     = leptonTwoP4.Pt();
        leptonTwoEta     = leptonTwoP4.Eta();
        leptonTwoPhi     = leptonTwoP4.Phi();

        if (params->selection == "mmg") {
            leptonOneIso    = particleSelector->GetMuonIsolation(muons[leptonOneIndex]);
            leptonOneFlavor = muons[leptonOneIndex]->q*13;
            leptonOneDZ     = muons[leptonOneIndex]->dz;
            leptonOneD0     = muons[leptonOneIndex]->d0;     
            leptonTwoIso    = particleSelector->GetMuonIsolation(muons[leptonTwoIndex]);
            leptonTwoFlavor = muons[leptonTwoIndex]->q*13;
            leptonTwoDZ     = muons[leptonTwoIndex]->dz;
            leptonTwoD0     = muons[leptonTwoIndex]->d0;
        }
        else if (params->selection == "eeg") {
            leptonOneIso    = particleSelector->GetElectronIsolation(electrons[leptonOneIndex], fInfo->rhoJet);
            leptonOneFlavor = electrons[leptonOneIndex]->q*11;
            leptonOneDZ     = electrons[leptonOneIndex]->dz;
            leptonOneD0     = electrons[leptonOneIndex]->d0; 
            leptonTwoIso    = particleSelector->GetElectronIsolation(electrons[leptonTwoIndex], fInfo->rhoJet);
            leptonTwoFlavor = electrons[leptonTwoIndex]->q*11;
            leptonTwoDZ     = electrons[leptonTwoIndex]->dz;
            leptonTwoD0     = electrons[leptonTwoIndex]->d0;
        }

        photonPt  = photonP4.Pt();
        photonEta  = photonP4.Eta();
        photonPhi  = photonP4.Phi();
        photonRawPt = photons[photonIndex]->pt;
        photonSCPt = photons[photonIndex]->scEt;
        photonSCEta = photons[photonIndex]->scEta;
        photonSCPhi = photons[photonIndex]->scPhi;
        photonMVA = photons[photonIndex]->mvaFall17V2;
        
        photonERes = photons[photonIndex]->ecalEnergyErrPostCorr / photonP4.E();
        photonE = photonP4.E();
        photonErrE = photons[photonIndex]->ecalEnergyErrPostCorr;
        if (sync) {
            std::cout << "photon resolution = " << photonERes << std::endl;
        }
        
        passElectronVeto = photons[photonIndex]->passElectronVeto;  

        photonRawR9 = photons[photonIndex]->r9_full5x5;
        if (!isData) {
            photonR9 = weights->GetCorrectedPhotonR9(*photons[photonIndex]);
            if (sync) {
                std::cout << "event, uncorrected r9, corrected r9: " << 
                             evtNumber << ", " <<   
                             photons[photonIndex]->r9_full5x5 << ", " << 
                             weights->GetCorrectedPhotonR9(*photons[photonIndex]) << std::endl;
            }
        }
        else photonR9 = photons[photonIndex]->r9_full5x5;
        if (photonRawR9 != 0.) {
            photonR9Weight = photonR9 / photonRawR9;
        }
        else {
            photonR9Weight = 1.;
        }

        // kinematic fit
        KinZfitter* kinZfitter = new KinZfitter(isData);
        std::vector<TObject *> selectedLeptons;
        std::vector<double> pTerr;
        if (params->selection == "mmg") {
            selectedLeptons.push_back(muons[leptonOneIndex]);
            selectedLeptons.push_back(muons[leptonTwoIndex]);
        }
        else if (params->selection == "eeg") {
            selectedLeptons.push_back(electrons[leptonOneIndex]);
            selectedLeptons.push_back(electrons[leptonTwoIndex]);
        }

        std::map<unsigned int, TLorentzVector> selectedFsrMap;
        kinZfitter->Setup(selectedLeptons, selectedFsrMap, pTerr);
        kinZfitter->KinRefitZ();
        std::vector<TLorentzVector> refitP4s;
        refitP4s = kinZfitter->GetRefitP4s();
        delete kinZfitter;

        TLorentzVector leptonOneP4KinFit = refitP4s[0];
        TLorentzVector leptonTwoP4KinFit = refitP4s[1];

        leptonOnePtKin = leptonOneP4KinFit.Pt();
        leptonTwoPtKin = leptonTwoP4KinFit.Pt();
        leptonOnePtKinErr = pTerr[0];
        leptonTwoPtKinErr = pTerr[1];

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
        llgMKin = (leptonOneP4KinFit + leptonTwoP4KinFit + photonP4).M();
        
        l1PhotonDEta = fabs(leptonOneP4.Eta() - photonP4.Eta());
        l1PhotonDPhi = fabs(leptonOneP4.DeltaPhi(photonP4));
        l1PhotonDR = leptonOneP4.DeltaR(photonP4);
        l2PhotonDEta = fabs(leptonTwoP4.Eta() - photonP4.Eta());
        l2PhotonDPhi = fabs(leptonTwoP4.DeltaPhi(photonP4));
        l2PhotonDR = leptonTwoP4.DeltaR(photonP4);
        
        if (l1PhotonDR > l2PhotonDR) {
            lPhotonDRMax = l1PhotonDR;
            lPhotonDRMin = l2PhotonDR; 
        }
        else {
            lPhotonDRMax = l2PhotonDR;
            lPhotonDRMin = l1PhotonDR;
        }

        dileptonPhotonDEta = fabs(dileptonP4.Eta() - photonP4.Eta());
        dileptonPhotonDPhi = fabs(dileptonP4.DeltaPhi(photonP4));
        dileptonPhotonDR = dileptonP4.DeltaR(photonP4);
        ptt = 2*fabs(dileptonP4.Px()*photonP4.Py() - photonP4.Px()*dileptonP4.Py())/llgP4.Pt();

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
        TLorentzVector vecg_in_Xframe = photonP4;
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
    
        if (!isData) {
            if (params->selection == "mmg") {
                leptonOneIDWeight = weights->GetHZZMuonIDEff(*muons[leptonOneIndex]); 
                leptonTwoIDWeight = weights->GetHZZMuonIDEff(*muons[leptonTwoIndex]);

                float sf11 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg1", *muons[leptonOneIndex]);
                float sf12 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg1", *muons[leptonTwoIndex]);
                float sf21 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg2", *muons[leptonOneIndex]);
                float sf22 = weights->GetDoubleMuonTriggerEffWeight("HLT_DoubleMuon_leg2", *muons[leptonTwoIndex]);
                
                if (leptonTwoPt < 20.) {
                    trigOneWeight = sf11;
                    trigTwoWeight = sf22;
                }
                else {
                    float prod1 = sf11*sf22;
                    float prod2 = sf21*sf12;
                    if (prod1 > prod2) {
                        trigOneWeight = sf11;
                        trigTwoWeight = sf22;
                    }
                    else {
                        trigOneWeight = sf21;
                        trigTwoWeight = sf12;
                    }
                }
            }

            else if (params->selection == "eeg") {
                leptonOneIDWeight = weights->GetElectronMVARecoIdEff(*electrons[leptonOneIndex]); 
                leptonTwoIDWeight = weights->GetElectronMVARecoIdEff(*electrons[leptonTwoIndex]); 
                
                float sf11 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg1", *electrons[leptonOneIndex]);
                float sf22 = weights->GetDoubleEGTriggerEffWeight("HLT_DoubleEG_leg2", *electrons[leptonTwoIndex]);

                trigOneWeight = sf11;
                trigTwoWeight = sf22;
            }

            eventWeight *= leptonOneIDWeight;
            eventWeight *= leptonTwoIDWeight;
            
            triggerWeight = trigOneWeight*trigTwoWeight;
            eventWeight *= triggerWeight;
           
            photonIDWeight = weights->GetPhotonMVAIdEff(*photons[photonIndex]); 
            eventWeight *= photonIDWeight;

            prefWeight = fInfo->prefweight;
            prefWeightUp = fInfo->prefweightUp;
            prefWeightDown = fInfo->prefweightDown;

            eventWeight *= prefWeight;
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
        TLorentzVector photonP4;
        photonP4.SetPtEtaPhiM(photonPt, photonEta, photonPhi, 0.);
        float min_phot_dr = 1000.;
        for (unsigned int i = 0; i < genPhotons.size(); i++) {
            TLorentzVector tmpGenPhot;
            tmpGenPhot.SetPtEtaPhiM(genPhotons[i]->pt, genPhotons[i]->eta, genPhotons[i]->phi, genPhotons[i]->mass);
            float this_dr = tmpGenPhot.DeltaR(photonP4);
            if (this_dr < min_phot_dr) {
                genPhotonP4.SetPtEtaPhiM(genPhotons[i]->pt, genPhotons[i]->eta, genPhotons[i]->phi, genPhotons[i]->mass);
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

    if (sync) {
        std::cout << "ABOUT TO FILL THE TREE" << std::endl;
        std::cout << "isDijetTag = " << isDijetTag << std::endl;
    }

    outTree->Fill();
    this->passedEvents++;

    //std::cout << "run, event, PU weight: " << runNumber << ", " << evtNumber << ", " << puWeight << std::endl;
    //std::cout << "run, event, ID weight 1, ID weight 2: " << 
    //              runNumber << ", " << evtNumber << ", " <<  elIDWeightOne << ", " << elIDWeightTwo << std::endl;
    //std::cout << "run, event, photon ID weight" << runNumber << ", " << evtNumber << ", " <<  photonIDWeight << std::endl;
    //std::cout << "run, event, trig weight 1, trig weight 2: " << 
    //              runNumber << ", " << evtNumber << ", " <<  elTrigWeightOne << ", " << elTrigWeightTwo << std::endl;
    //std::cout << "run, event, trig weight 1, trig weight 2: " << 
    //              runNumber << ", " << evtNumber << ", " <<  muonTrigWeightOne << ", " << muonTrigWeightTwo << std::endl;

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
