#pragma once
#include <vector>
#include <string>
#include <cmath>

class Utils {
public:
    static std::vector<std::string> parseLine3(std::string iniStr);
    static bool renorm(std::vector<double>& s);
    static std::vector<double> vecProd(const std::vector<double>& A, const std::vector<double>& B);
    static double scalProd(const std::vector<double>& A, const std::vector<double>& B);
    static std::vector<double> vecProd(double* A, double* B);
    static double scalProd(double* A, double* B);
};