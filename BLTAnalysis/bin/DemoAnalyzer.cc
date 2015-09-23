#include "DemoAnalyzer.hh"

//
// See header file for class documentation
//

// _____________________________________________________________________________
// Useful functions

namespace baconhep {
    class TMET : public TObject
    {
    public:
        TMET():  pt(0), phi(0) {}
        ~TMET() {}
        float pt, phi;
    };
}

std::ostream& operator<<(std::ostream& os, const baconhep::TGenParticle* p) {
    return os << "pdgId: " << p->pdgId << " status: " << p->status << " parent: " << p->parent << " pt: " << p->pt << " eta: " << p->eta << " phi: " << p->phi << " mass: " << p->mass << " y: " << p->y;
}

std::ostream& operator<<(std::ostream& os, const baconhep::TVertex* p) {
    return os << "z: " << p->z << " perp: " << quad_sum(p->x, p->y) << " ndof: " << p->ndof;
}

std::ostream& operator<<(std::ostream& os, const baconhep::TMuon* p) {
    return os << "pt: " << p->pt << " eta: " << p->eta << " phi: " << p->phi << " q: " << p->q;
}

std::ostream& operator<<(std::ostream& os, const baconhep::TElectron* p) {
    return os << "pt: " << p->pt << " eta: " << p->eta << " phi: " << p->phi << " q: " << p->q;
}

std::ostream& operator<<(std::ostream& os, const baconhep::TPhoton* p) {
    return os << "pt: " << p->pt << " eta: " << p->eta << " phi: " << p->phi;
}

std::ostream& operator<<(std::ostream& os, const baconhep::TTau* p) {
    return os << "pt: " << p->pt << " eta: " << p->eta << " phi: " << p->phi << " q: " << p->q;
}

std::ostream& operator<<(std::ostream& os, const baconhep::TJet* p) {
    return os << "pt: " << p->pt << " eta: " << p->eta << " phi: " << p->phi << " mass: " << p->mass;
}

std::ostream& operator<<(std::ostream& os, const baconhep::TMET* p) {
    return os << "pt: " << p->pt << " phi: " << p->phi;
}

template<class T>
void copy_p4(const T* lhs, float mass, TLorentzVector& rhs) {
    rhs.SetPtEtaPhiM(lhs->pt, lhs->eta, lhs->phi, mass);
}

//template<class T>
//bool sort_by_higher_pt(const T& lhs, const T& rhs) {
//    return lhs.pt > rhs.pt;
//}

//template<class T>
//bool sort_by_higher_pt(const T* lhs, const T* rhs) {
//    return lhs->pt > rhs->pt;
//}

template<class T>
bool sort_by_higher_pt(const T* lhs, const T* rhs) {
    return lhs->pt > rhs->pt;
}

//bool sort_by_higher_pt(const baconhep::TMuon* lhs, const baconhep::TMuon* rhs) {
//    return lhs->pt > rhs->pt;
//}

static const double ELE_MASS  = 0.000511;
static const double MUON_MASS = 0.105658369;

static const int ELE_PDGID  = 11;  // e-
static const int MUON_PDGID = 13;  // mu-
static const int Z_PDGID = 23;


// _____________________________________________________________________________
// DemoAnalyzer implementations

using namespace baconhep;

DemoAnalyzer::DemoAnalyzer() : BLTSelector()
{

}

DemoAnalyzer::~DemoAnalyzer()
{

}

void DemoAnalyzer::Begin(TTree *tree)
{
    // Parse command line option
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1),
              std::sregex_token_iterator(), std::back_inserter(options));
    assert (options.size() == 7);

    // Set the parameters
    params.reset(new Parameters());
    params->suffix    = options[0];
    params->abcd      = options[1];
    params->selection = options[2];
    params->period    = options[3];
    params->dataname  = options[4];
    params->jobcount  = options[5];
    params->pileup    = options[6];

    // Set the cuts
    cuts.reset(new Cuts());
    triggerSelector.reset(new TriggerSelector());
    particleSelector.reset(new ParticleSelector(*params, *cuts));

    // Prepare the output tree
    outFileName = params->get_output_filename("demoFile");
    outTreeName = params->get_output_treename("demoTree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "demoTree");

    outTree->Branch("muonOne", &muonOne);
    outTree->Branch("muonTwo", &muonTwo);
    outTree->Branch("dimuon", &dimuon);
    outTree->Branch("genMuonOne", &genMuonOne);
    outTree->Branch("genMuonTwo", &genMuonTwo);
    outTree->Branch("genZ", &genZ);
    outTree->Branch("genMuonOneCharge", &genMuonOneCharge);
    outTree->Branch("genMuonTwoCharge", &genMuonTwoCharge);

    ReportPostBegin();
}

Bool_t DemoAnalyzer::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;

    //if (entry%1==0)  std::cout << "... Processing event: " << entry << "." << std::endl;
    //if (entry%1==0)  std::cout << "... Processing event: " << entry << " Run: " << fInfo->runNum << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum << "." << std::endl;

    const bool isRealData = !(fInfo->runNum == 1);

    bool printEvent = false;

    //////////////////
    // GenParticles //
    //////////////////

    if (printEvent) {
        if (!isRealData) {
            for (int i=0; i<fGenParticleArr->GetEntries(); i++) {
                const TGenParticle* thisGenParticle = (TGenParticle*) fGenParticleArr->At(i);
                assert(thisGenParticle);
                std::cout << "GenParticle " << i << ": " << thisGenParticle << std::endl;
            }
        }
    }

    //////////////
    // Vertices //
    //////////////

    if (printEvent) {
        for (int i=0; i<fPVArr->GetEntries(); i++) {
            const TVertex* thisVertex = (TVertex*) fPVArr->At(i);
            assert(thisVertex);
            std::cout << "Vertex " << i << ": " << thisVertex << std::endl;
        }
    }

    ///////////
    // Muons //
    ///////////

    if (printEvent) {
        for (int i=0; i<fMuonArr->GetEntries(); i++) {
            const TMuon* thisMuon = (TMuon*) fMuonArr->At(i);
            assert(thisMuon);
            std::cout << "Muon " << i << ": " << thisMuon << std::endl;
        }
    }

    ///////////////
    // Electrons //
    ///////////////

    if (printEvent) {
        for (int i=0; i<fElectronArr->GetEntries(); i++) {
            const TElectron* thisElectron = (TElectron*) fElectronArr->At(i);
            assert(thisElectron);
            std::cout << "Electron " << i << ": " << thisElectron << std::endl;
        }
    }

    /////////////
    // Photons //
    /////////////

    if (printEvent) {
        for (int i=0; i<fPhotonArr->GetEntries(); i++) {
            const TPhoton* thisPhoton = (TPhoton*) fPhotonArr->At(i);
            assert(thisPhoton);
            std::cout << "Photon " << i << ": " << thisPhoton << std::endl;
        }
    }

    //////////
    // Taus //
    //////////

    if (printEvent) {
        for (int i=0; i<fTauArr->GetEntries(); i++) {
            const TTau* thisTau = (TTau*) fTauArr->At(i);
            assert(thisTau);
            std::cout << "Tau " << i << ": " << thisTau << std::endl;
        }
    }

    //////////
    // Jets //
    //////////

    if (printEvent) {
        for (int i=0; i<fAK4CHSArr->GetEntries(); i++) {
            const TJet* thisJet = (TJet*) fAK4CHSArr->At(i);
            assert(thisJet);
            std::cout << "Jet " << i << ": " << thisJet << std::endl;
        }
    }

    /////////
    // MET //
    /////////

    TMET* pfMET = new TMET();
    pfMET->pt = fInfo->pfMET;
    pfMET->phi = fInfo->pfMETphi;

    TMET* caloMET = new TMET();
    caloMET->pt = fInfo->caloMET;
    caloMET->phi = fInfo->caloMETphi;

    if (printEvent) {
        std::cout << "MET " << "(PF)" << ": " << pfMET << std::endl;
        std::cout << "MET " << "(C) " << ": " << caloMET << std::endl;
    }

    ////////////
    // Select //
    ////////////

    std::vector<const TMuon*> muons;

    for (int i=0; i<fMuonArr->GetEntries(); i++) {
        const TMuon* thisMuon = (TMuon*) fMuonArr->At(i);
        assert(thisMuon);

        if (thisMuon->pt > 20 && std::abs(thisMuon->eta) < 2.4)
            muons.push_back(thisMuon);
    }

    std::sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);

    bool found_dimuon = false;
    TLorentzVector tmp_muonOne;
    TLorentzVector tmp_muonTwo;
    TLorentzVector tmp_dimuon;

    for (unsigned i=0; i<muons.size(); ++i) {
        const TMuon* imuon = muons.at(i);
        copy_p4(imuon, MUON_MASS, tmp_muonOne);

        for (unsigned j=i+1; j<muons.size(); ++j) {
            const TMuon* jmuon = muons.at(j);
            copy_p4(jmuon, MUON_MASS, tmp_muonTwo);

            tmp_dimuon = tmp_muonOne + tmp_muonTwo;

            if (imuon->q != jmuon->q) {
                if (70 < tmp_dimuon.M() && tmp_dimuon.M() < 120) {
                    found_dimuon = true;
                    break;
                }
            }
        }

        if (found_dimuon)
            break;
    }

    if (!found_dimuon)
        return kFALSE;

    TLorentzVector tmp_genMuonOne;
    TLorentzVector tmp_genMuonTwo;
    TLorentzVector tmp_genZ;
    int tmp_genMuonOneCharge = 0;
    int tmp_genMuonTwoCharge = 0;

/*
    int idx_genZ = -1;
    for (int i=0; i<fGenParticleArr->GetEntries(); i++) {
        const TGenParticle* thisGenParticle = (TGenParticle*) fGenParticleArr->At(i);
        assert(thisGenParticle);

        if (std::abs(thisGenParticle->pdgId) == Z_PDGID && thisGenParticle->status == 22) {
            idx_genZ = i;
            copy_p4(thisGenParticle, thisGenParticle->mass, tmp_genZ);
        }
    }

    assert(idx_genZ != -1);
*/

    int idx_genMuonOne = -1;
    int idx_genMuonTwo = -1;
    for (int i=0; i<fGenParticleArr->GetEntries(); i++) {
        const TGenParticle* thisGenParticle = (TGenParticle*) fGenParticleArr->At(i);
        assert(thisGenParticle);

        if (std::abs(thisGenParticle->pdgId) == 13 && thisGenParticle->status == 1) {
            // Check parent
            const TGenParticle* thisGenParticleParent = 0;
            if (thisGenParticle->parent >= 0) {
                thisGenParticleParent = (TGenParticle*) fGenParticleArr->At(thisGenParticle->parent);
                while (thisGenParticleParent->pdgId == thisGenParticle->pdgId && thisGenParticleParent->parent >= 0) {
                    thisGenParticleParent = (TGenParticle*) fGenParticleArr->At(thisGenParticleParent->parent);
                }
            }

            if (!thisGenParticleParent || thisGenParticleParent->pdgId != Z_PDGID)
                continue;

            if (idx_genMuonOne == -1) {
                idx_genMuonOne = i;
                copy_p4(thisGenParticle, thisGenParticle->mass, tmp_genMuonOne);
                tmp_genMuonOneCharge = -1 * thisGenParticle->pdgId / std::abs(thisGenParticle->pdgId);

            } else if (idx_genMuonTwo == -1) {
                idx_genMuonTwo = i;
                copy_p4(thisGenParticle, thisGenParticle->mass, tmp_genMuonTwo);
                tmp_genMuonTwoCharge = -1 * thisGenParticle->pdgId / std::abs(thisGenParticle->pdgId);

            } else {
                assert(false);
            }
        }
    }

    assert((idx_genMuonOne == -1 && idx_genMuonTwo == -1) || (idx_genMuonOne != -1 && idx_genMuonTwo != -1));

    if (idx_genMuonOne == -1 || idx_genMuonTwo == -1)  // FIXME
        return kFALSE;


    //////////
    // Fill //
    //////////

    muonOne = tmp_muonOne;
    muonTwo = tmp_muonTwo;
    dimuon  = tmp_dimuon;

    if (tmp_genMuonOne.Pt() > tmp_genMuonTwo.Pt()) {
        genMuonOne = tmp_genMuonOne;
        genMuonTwo = tmp_genMuonTwo;
        genMuonOneCharge = tmp_genMuonOneCharge;
        genMuonTwoCharge = tmp_genMuonTwoCharge;
    } else {
        genMuonOne = tmp_genMuonTwo;
        genMuonTwo = tmp_genMuonOne;
        genMuonOneCharge = tmp_genMuonTwoCharge;
        genMuonTwoCharge = tmp_genMuonOneCharge;
    }
    genZ = genMuonOne + genMuonTwo;

    outTree->Fill();
    this->passedEvents++;

    delete pfMET;
    delete caloMET;

    return kTRUE;
}

void DemoAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void DemoAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job ==========" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  =========================" << std::endl;
}

void DemoAnalyzer::ReportPostTerminate()
{
    std::cout << "  ==== Terminate Job ======" << std::endl;
    std::cout << "  output   : " << params->get_output_filename("demoFile") << std::endl;
    std::cout << "           : Processed " << this->fileCount << " files with " << this->unskimmedEventCount << " unskimmed events." << std::endl;
    std::cout << "           : Selected " << this->passedEvents << " / " << this->totalEvents << " events." << std::endl;
    std::cout << "  =========================" << std::endl;
}


// _____________________________________________________________________________
// Main function

int main(int argc, char **argv)
{
    std::unique_ptr<DemoAnalyzer> selector(new DemoAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
