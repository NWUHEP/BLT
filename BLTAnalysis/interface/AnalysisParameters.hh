#ifndef ANALYSISPARAMETERS_HH
#define ANALYSISPARAMETERS_HH

#include <string>

class AnalysisParameters {
public:
    AnalysisParameters() {}
    ~AnalysisParameters() {}

    std::string get_output_filename(const std::string& name) {
        return name + "_" + this->dataname + "_" + this->jobcount + ".root";
    }

    std::string get_output_treename(const std::string& name) {
        return name + "_" + this->suffix;
    }

    std::string suffix;
    std::string abcd;
    std::string selection;
    std::string period;
    std::string dataname;
    std::string jobcount;
    std::string pileup;
};

// Output streams
std::ostream& operator<<(std::ostream& os, const AnalysisParameters& params);

#endif  // ANALYSISPARAMETERS_HH 
