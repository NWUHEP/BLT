#ifndef BLTHELPER_HH
#define BLTHELPER_HH

// Bacon header files
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TTau.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TAddJet.hh"

#include <TLorentzVector.h>
#include <TVector3.h>

#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>


// Constants
static const double ELE_MASS  = 0.000511;
static const double MUON_MASS = 0.105658369;

static const int ELE_PDGID  = 11;  // e-
static const int MUON_PDGID = 13;  // mu-
static const int TAU_PDGID = 15;   // tau-
static const int Z_PDGID = 23;

// Math
double quad_sum(double x, double y);

// String manipulation
std::string get_file_extension(const std::string& filename);
bool starts_with(const std::string& str1, const std::string& str2);
bool ends_with(const std::string& str1, const std::string& str2);

// Copy functions
template<class T>
void copy_xyz(const T* lhs, TVector3& rhs) {
    rhs.SetXYZ(lhs->x, lhs->y, lhs->z);
}

template<class T>
void copy_p4(const T* lhs, float mass, TLorentzVector& rhs) {
    rhs.SetPtEtaPhiM(lhs->pt, lhs->eta, lhs->phi, mass);
}

// Sort functions
template<class T>
bool sort_by_higher_pt(const T* lhs, const T* rhs) {
    return lhs->pt > rhs->pt;
}

template<class T>
bool sort_by_btag(const T* lhs, const T* rhs) {
    return lhs->bmva > rhs->bmva;
}

// Basic classes
namespace baconhep {
    class TMET : public TObject
    {
    public:
        TMET():  pt(0), phi(0) {}
        ~TMET() {}
        float pt, phi;
    };
}

// Stuff for ostream
std::ostream& operator<<(std::ostream& os, const baconhep::TGenParticle* p);
std::ostream& operator<<(std::ostream& os, const baconhep::TVertex* p);
std::ostream& operator<<(std::ostream& os, const baconhep::TMuon* p);
std::ostream& operator<<(std::ostream& os, const baconhep::TElectron* p);
std::ostream& operator<<(std::ostream& os, const baconhep::TPhoton* p);
std::ostream& operator<<(std::ostream& os, const baconhep::TTau* p);
std::ostream& operator<<(std::ostream& os, const baconhep::TJet* p);
std::ostream& operator<<(std::ostream& os, const baconhep::TMET* p);

// Stuff for cout
std::string error();
std::string warning();
std::string info();
std::string debug();

#endif  // BLTHELPER_HH

