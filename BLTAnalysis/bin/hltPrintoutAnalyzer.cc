#include "hltPrintoutAnalyzer.h"
#include <map>
#include <fstream>
#include <math.h>

#include <TSystem.h>
#include <TF2.h>

// prints a message and exits gracefully
#define FATAL(msg) do { fprintf(stderr, "FATAL: %s\n", msg); gSystem->Exit(1); } while (0)

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

hltPrintoutAnalyzer::hltPrintoutAnalyzer() : BLTSelector()
{

}

hltPrintoutAnalyzer::~hltPrintoutAnalyzer()
{

}

void hltPrintoutAnalyzer::Begin(TTree *tree)
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

    // Prepare the output file
    string outFileName = params->get_output_filename("output");
    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",30,0.5,30.5);
    
    // Prepare the output tree
    string treeName = params->get_output_treename("tree");
    outTree = new TTree(treeName.c_str(), "bltTree");

    // event data
    outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber, "eventNumber/l");
    outTree->Branch("lumiSection", &lumiSection);
       
    ReportPostBegin();
}

Bool_t hltPrintoutAnalyzer::Process(Long64_t entry)
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

    /* Apply lumi mask */
    RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
    if(!lumiMask.HasRunLumi(rl)) 
        return kTRUE;
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

    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;

    outTree->Fill();
    this->passedEvents++;

    return kTRUE;
}

void hltPrintoutAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void hltPrintoutAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void hltPrintoutAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<hltPrintoutAnalyzer> selector(new hltPrintoutAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
