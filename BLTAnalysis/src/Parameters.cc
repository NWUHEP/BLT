#include "BLT/BLTAnalysis/interface/Parameters.hh"

#include "BLT/BLTAnalysis/interface/BLTHelper.hh"

#include <iostream>
#include <stdexcept>

void Parameters::setup(const std::vector<std::string>& options) {
    if (options.size() < 5) {
        std::cout << error() << "Expect 5 command line arguments, received " << options.size() << "." << std::endl;
        throw std::runtime_error("bad arguments");
    }

    this->dataset      = options.at(0);
    this->datasetgroup = options.at(1);
    this->selection    = options.at(2);
    this->period       = options.at(3);
    this->jobid        = options.at(4);
}


std::ostream& operator<<(std::ostream& os, const Parameters& params) {
    return os << "  dataset      : " << params.dataset      << "\n"
              << "  datasetgroup : " << params.datasetgroup << "\n"
              << "  selection    : " << params.selection    << "\n"
              << "  period       : " << params.period       << "\n"
              << "  jobid        : " << params.jobid        ;
}
