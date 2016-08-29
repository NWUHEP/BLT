#include "DimuonAnalyzerGen.hh"
#include <cmath>

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 

DimuonAnalyzer::DimuonAnalyzer() : BLTSelector()
{

}

DimuonAnalyzer::~DimuonAnalyzer()
{

}

void DimuonAnalyzer::Begin(TTree *tree)
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
    std::string trigfilename = cmssw_base + std::string("/src/BaconAna/DataFormats/data/HLTFile_25ns");
    trigger.reset(new baconhep::TTrigger(trigfilename));

    // Prepare the output tree
    outFileName = params->get_output_filename("output");
    outTreeName = params->get_output_treename("data");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "demoTree");

    outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber);
    outTree->Branch("lumiSection", &lumiSection);

    // Generator branches
    outTree->Branch("GenZd", &GenZd);
    outTree->Branch("GenMuonLead", &GenMuonLead);
    outTree->Branch("GenMuonTrail", &GenMuonTrail);
    outTree->Branch("GenBQuark", &GenBQuark);
    outTree->Branch("GenQPrime", &GenQPrime);


    // Histograms for muon difference variables
    //TH1F *muonDeltaPhi = new TH1F("muonDeltaPhi", "muonDeltaPhi", 20, 0., 5.);
    //TH1F *muonDeltaEta = new TH1F("muonDeltaEta", "muonDeltaEta", 20, 0., 5.);
    //TH1F *muonDeltaR = new TH1F("muonDeltaR", "muonDeltaR", 20, 0., 5.);

    outTree->Branch("muonDeltaPhi", &muonDeltaPhi);
    outTree->Branch("muonDeltaEta", &muonDeltaEta);
    outTree->Branch("muonDeltaR", &muonDeltaR);

    outTree->Branch("costheta", &costheta);
        

    ReportPostBegin();
}

Bool_t DimuonAnalyzer::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;

    //if (entry%1==0)  std::cout << "... Processing event: " << entry << "." << std::endl;
    if (entry%10000==0)  std::cout << "... Processing event: " << entry << " Run: " << fInfo->runNum << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum << "." << std::endl;

    const bool isRealData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isRealData);

    /* Trigger selection */
    bool passTrigger;
    passTrigger= true;//trigger->pass("HLT_IsoMu24_eta2p1_v*", fInfo->triggerBits);
    if (!passTrigger)
        return kTRUE;

    ///////////////////
    // Select objects//
    ///////////////////

    bool printEvent = false;

    //////////////////
    // GenParticles //
    //////////////////

    std::vector<TGenParticle*> GenMuons;
    std::vector<TGenParticle*> GenJets;
    TLorentzVector GenMuonMinus;

    if (!isRealData) {
        for (int i=0; i<fGenParticleArr->GetEntries(); i++) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
            assert(particle);

            if (printEvent)
                std::cout << "GenParticle " << i << ": " << particle << std::endl;
            
            /*if (particle->status == 23)
            {*/
                if (abs(particle->pdgId) == 13)
                {    
                    int iParent = particle->parent;
                    if (iParent > 0)
                    {
                        TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(iParent);
                        if (abs(mother->pdgId) == 4900023)
                        {
                            copy_p4(mother, mother->mass, GenZd);
                            GenMuons.push_back(particle); 
                            //std::cout << "mother status " << mother->status << std::endl;
                            //std::cout << "muon status " << particle->status << std::endl;
                        }
                    }
                }

                else
                {
                    if (particle->status == 23)
                    {
                        if (abs(particle->pdgId) == 5)
                        {
                            copy_p4(particle, particle->mass, GenBQuark);
                            GenJets.push_back(particle);
                           // std::cout << "b jet status " << particle->status << std::endl;
                        }
                        else if (!(abs(particle->pdgId) == 5) && abs(particle->pdgId) < 6)
                        {    
                            copy_p4(particle, particle->mass, GenQPrime);
                            GenJets.push_back(particle);
                           // std::cout << "other jet status " << particle->status << std::endl;
                        }
                    }
                }
            //}
        }
    }

    std::sort(GenMuons.begin(), GenMuons.end(), sort_by_higher_pt<TGenParticle>);

    TGenParticle* GenMuon1;
    TGenParticle* GenMuon2;
    if (GenMuons.size() == 2)
    {
        GenMuon1 = GenMuons[0];
        GenMuon2 = GenMuons[1];
        if (GenMuon1->pt > GenMuon2->pt)
        {
            copy_p4(GenMuon1, MUON_MASS, GenMuonLead);
            copy_p4(GenMuon2, MUON_MASS, GenMuonTrail);
        }
        else
        {
            copy_p4(GenMuon1, MUON_MASS, GenMuonTrail);
            copy_p4(GenMuon2, MUON_MASS, GenMuonLead);
        }

        if (GenMuon1->pdgId == 13)
            copy_p4(GenMuon1, MUON_MASS, GenMuonMinus);
        else
            copy_p4(GenMuon2, MUON_MASS, GenMuonMinus);


        /* Basic Gen Muon Selection */
        /*if ((GenMuon1->pt < 5.0) || (GenMuon2->pt < 5.0))
            return kTRUE;*/

        if ((abs(GenMuon1->eta) > 2.4) || (abs(GenMuon2->eta) > 2.4))
            return kTRUE;

        /* Muon Difference Variables */
        muonDeltaPhi = abs(GenMuon1->phi - GenMuon2->phi);
        muonDeltaEta = abs(GenMuon1->eta - GenMuon2->eta);
        muonDeltaR = sqrt(pow(muonDeltaPhi, 2.0) + pow(muonDeltaEta, 2.0));


    }
    else
        return kTRUE;

    /* Cos(theta) variable */
    TLorentzVector Zd_clone = GenZd;
    TLorentzVector B_clone = GenBQuark;
    TLorentzVector mu_minus_clone = GenMuonMinus; // Define angle with respect to the negative muon

    TVector3 boost = -1*Zd_clone.BoostVector(); // The boost to move to the Zd rest frame
    B_clone.Boost(boost);  // Boost the B to Zd rest frame
    Zd_clone.Boost(boost); // And the Zd itself

    TVector3 B_threevect = B_clone.Vect();
    TLorentzVector beamAxis(0,0,1,1);
    beamAxis.Boost(boost); 

    // Rotate the axes of the frame to complete the Lorentz transformation
    TVector3 axis_z_CM = (-1*B_threevect).Unit(); // Define z axis (angular reference axis) as opposite the direction of the mother particle in the Zd frame
    TVector3 axis_y_CM = beamAxis.Vect().Cross(B_threevect).Unit();
    TVector3 axis_x_CM = axis_y_CM.Cross(axis_z_CM).Unit();
    TRotation rotation;
    rotation = rotation.RotateAxes(axis_x_CM, axis_y_CM, axis_z_CM).Inverse();

    mu_minus_clone.Boost(boost);
    mu_minus_clone.Transform(rotation); // rotate and boost the muon

    costheta = mu_minus_clone.CosTheta(); // Feed costheta to the root tree branch


    /* Basic Gen Jet Selection */

    std::sort(GenJets.begin(), GenJets.end(), sort_by_higher_pt<TGenParticle>);
    TGenParticle* GenJet1;
    TGenParticle* GenJet2;

    if (GenJets.size()==2)
    {
        GenJet1 = GenJets[0];
        GenJet2 = GenJets[1];
        if ((abs(GenJet1->eta) > 4.7) || (abs(GenJet2->eta) > 4.7))
            return kTRUE;
    }
    else
    {
        return kTRUE;
    }

    if (abs(GenJet1->pdgId) == 5)
    {
        if (abs(GenJet1->eta) > 2.4)
            return kTRUE;
    }

    else if (abs(GenJet2->pdgId) == 5)
    {
        if (abs(GenJet2->eta) > 2.4)
            return kTRUE;
    }

    
    /* RECO Acceptance */
    std::vector<TMuon*> muons;

    for (int i=0; i<fMuonArr->GetEntries(); i++) 
    {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);
        muons.push_back(muon);
    }

    std::sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);

    TMuon* Muon1;
    TMuon* Muon2;
    if (muons.size() >= 2)
    {
        Muon1 = muons[0];
        Muon2 = muons[1];
        if ((abs(Muon1->eta) > 2.4) || (abs(Muon2->eta) > 2.4))
            return kTRUE;    
        if ((Muon1->pt < 5.0) || (Muon2->pt < 5.0))
            return kTRUE;    
    }
    
    else
        return kTRUE;

    /* Jets */        
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    std::vector<TJet*> bjets;
    std::vector<TJet*> otherjets;
    std::vector<TJet*> alljets;

    for (int i=0; i < jetCollection->GetEntries(); i++) 
    {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

        if (jet->csv > 0.898)
        {
            bjets.push_back(jet);
            alljets.push_back(jet);
        }
        else
        {
            otherjets.push_back(jet);
            alljets.push_back(jet);
        }

    }

    std::sort(bjets.begin(), bjets.end(), sort_by_higher_pt<TJet>);
    std::sort(otherjets.begin(), otherjets.end(), sort_by_higher_pt<TJet>);
    std::sort(alljets.begin(), alljets.end(), sort_by_higher_pt<TJet>);

    if (alljets.size() < 2)
        return kTRUE;

    TJet* recoBJet;
    if (bjets.size() > 0)
    {
        recoBJet = bjets[0];
        if (abs(recoBJet->eta) > 2.4)
            return kTRUE;
    }
    else 
        return kTRUE;

    TJet* recoOtherJet;
    if (otherjets.size() > 0)
    {
        recoOtherJet = otherjets[0];
        if (abs(recoOtherJet->eta) > 4.7)
            return kTRUE;
    }
    else
        return kTRUE;

    /* Gen-RECO Matching */
    TGenParticle* GenBJet;
    TGenParticle* GenOtherJet;
    if (abs(GenJet1->pdgId) == 5)
    {
        GenBJet = GenJet1;
        GenOtherJet = GenJet2;
    }
    else
    {
        GenBJet = GenJet2;
        GenOtherJet = GenJet1;
    }
   
    float DeltaRLeading = sqrt(pow(abs(GenMuon1->phi - Muon1->phi), 2.0) + pow(abs(GenMuon1->eta - Muon1->eta), 2.0));
    float DeltaRTrailing = sqrt(pow(abs(GenMuon2->phi - Muon2->phi), 2.0) + pow(abs(GenMuon2->eta - Muon2->eta), 2.0));


    std::vector<float> BDeltaR_vec; 
    for (unsigned i = 0; i < alljets.size(); i++)
    {
        float myDeltaRval = sqrt(pow(abs(GenBJet->phi - alljets[i]->phi), 2.0) + pow(abs(GenBJet->eta - alljets[i]->eta), 2.0));
        BDeltaR_vec.push_back(myDeltaRval);
    }
    std::sort(BDeltaR_vec.begin(), BDeltaR_vec.end(), std::less<float>());
    float DeltaRBJet = BDeltaR_vec[0];

    std::vector<float> otherDeltaR_vec; 
    for (unsigned i = 0; i < alljets.size(); i++)
    {
        float myDeltaRval = sqrt(pow(abs(GenOtherJet->phi - alljets[i]->phi), 2.0) + pow(abs(GenOtherJet->eta - alljets[i]->eta), 2.0));
        otherDeltaR_vec.push_back(myDeltaRval);
    }
    std::sort(otherDeltaR_vec.begin(), otherDeltaR_vec.end(), std::less<float>());
    float DeltaROtherJet = otherDeltaR_vec[0];
    
    if ((DeltaRLeading > 0.5) || (DeltaRTrailing > 0.5) || (DeltaRBJet > 0.5) || (DeltaROtherJet > 0.5))
        return kTRUE;  

    runNumber = fInfo->runNum;
    evtNumber = fInfo->evtNum;
    lumiSection = fInfo->lumiSec;
    outTree->Fill();
    this->passedEvents++;
    return kTRUE;
}

void DimuonAnalyzer::Terminate()
{
    outFile->Write();   
    outFile->Close();

    ReportPostTerminate();
}

void DimuonAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void DimuonAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<DimuonAnalyzer> selector(new DimuonAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
