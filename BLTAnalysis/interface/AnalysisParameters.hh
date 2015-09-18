#ifndef ANALYSISPARAMETERS_HH
#define ANALYSISPARAMETERS_HH

#include <string>

class AnalysisParameters {
public:
    AnalysisParameters() {}
    ~AnalysisParameters() {}

    std::string selection;
    std::string period;
    std::string abcd;
    std::string suffix;
    std::string dataname;
    std::string jobCount;
};

#endif  // ANALYSISPARAMETERS_HH 
