#ifndef PARAMETERS_HH
#define PARAMETERS_HH

#include <string>

class Parameters {
public:
    Parameters() {}
    ~Parameters() {}

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
std::ostream& operator<<(std::ostream& os, const Parameters& params);

#endif  // PARAMETERS_HH 
