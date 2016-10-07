#ifndef PARAMETERS_HH
#define PARAMETERS_HH

#include <vector>
#include <string>

class Parameters {
public:
    Parameters() {}
    ~Parameters() {}

    std::string get_output_filename(const std::string& name) {
        return name + "_" + this->dataset + "_" + this->jobid + ".root";
    }

    std::string get_output_treename(const std::string& name) {
        return name + "_" + this->datasetgroup;
    }

    void setup(const std::vector<std::string>& options);

    std::string dataset;
    std::string datasetgroup;
    std::string selection;
    std::string period;
    std::string jobid;
};

// Output streams
std::ostream& operator<<(std::ostream& os, const Parameters& params);

#endif  // PARAMETERS_HH 
