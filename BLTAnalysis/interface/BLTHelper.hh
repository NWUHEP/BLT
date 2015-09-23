#ifndef BLTHELPER_HH
#define BLTHELPER_HH

#include <string>
#include <cmath>

std::string get_file_extension(const std::string& filename);
bool starts_with(const std::string& str1, const std::string& str2);
bool ends_with(const std::string& str1, const std::string& str2);

double quad_sum(double x, double y);

std::string error();
std::string warning();
std::string info();
std::string debug();

#endif  // BLTHELPER_HH

