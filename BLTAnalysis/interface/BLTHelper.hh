#ifndef BLTHELPER_HH
#define BLTHELPER_HH

#include <string>

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

std::string Error() {
    return "ERROR  : ";
}

std::string Warning() {
    return "WARNING: ";
}

std::string Info() {
    return "INFO   : ";
}

std::string Debug() {
    return "DEBUG  : ";
}

#endif  // BLTHELPER_HH
