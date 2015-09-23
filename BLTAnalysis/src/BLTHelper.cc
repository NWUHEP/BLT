#include "BLT/BLTAnalysis/interface/BLTHelper.hh"

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

double quad_sum(double x, double y) {
    return std::sqrt(x*x + y*y);
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

