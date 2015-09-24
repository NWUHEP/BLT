#include "BLT/BLTAnalysis/interface/BLTHelper.hh"

#include <iostream>

double quad_sum(double x, double y) {
    return std::sqrt(x*x + y*y);
}

std::string get_file_extension(const std::string& filename) {
    if (filename.find_last_of(".") != std::string::npos)
        return filename.substr(filename.find_last_of(".")+1);
    return "";
}

bool starts_with(const std::string& str1, const std::string& str2) {
    return str1.size() >= str2.size() && str1.compare(0, str2.size(), str2) == 0;
}

bool ends_with(const std::string& str1, const std::string& str2) {
    return str1.size() >= str2.size() && str1.compare(str1.size()-str2.size(), str2.size(), str2) == 0;
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

std::string error() {
    return "[ERROR  ] ";
}

std::string warning() {
    return "[WARNING] ";
}

std::string info() {
    return "[INFO   ] ";
}

std::string debug() {
    return "[DEBUG  ] ";
}

