#include "Utils.hpp"

std::vector<std::string> Utils::parseLine3(std::string iniStr) {
    std::vector<std::string> result;
    int startPos(0), endPos(0), status(0);
    std::string tStr;
    for (auto& x : iniStr) {
        if ((x == ' ') || (x == '\t') || (x == ',') || (x == ';')) {
            if (status == 1) {
                status = 0;
                result.push_back(iniStr.substr(startPos, endPos - startPos));
            }
        } else {
            if (status == 0) {
                status = 1;
                startPos = endPos;
            }
        }
        endPos++;
    }
    if (status == 1) {
        result.push_back(iniStr.substr(startPos, endPos - startPos));
    }
    return result;
}

bool Utils::renorm(std::vector<double>& s) {
    double length(0);
    for (size_t i = 0; i < 3; ++i) length += s[i] * s[i];
    length = sqrt(length);
    if (length < 1e-10) return false;
    for (size_t i = 0; i < 3; ++i) s[i] /= length;
    return true;
}

std::vector<double> Utils::vecProd(const std::vector<double>& A, const std::vector<double>& B) {
    std::vector<double> res(3);
    res[0] = A[1] * B[2] - A[2] * B[1];
    res[1] = A[2] * B[0] - A[0] * B[2];
    res[2] = A[0] * B[1] - A[1] * B[0];
    return res;
}

double Utils::scalProd(const std::vector<double>& A, const std::vector<double>& B) {
    double res(0);
    res += A[0] * B[0];
    res += A[1] * B[1];
    res += A[2] * B[2];
    return res;
}

std::vector<double> Utils::vecProd(double* A, double* B) {
    std::vector<double> res(3);
    res[0] = A[1] * B[2] - A[2] * B[1];
    res[1] = A[2] * B[0] - A[0] * B[2];
    res[2] = A[0] * B[1] - A[1] * B[0];
    return res;
}

double Utils::scalProd(double* A, double* B) {
    double res(0);
    res += A[0] * B[0];
    res += A[1] * B[1];
    res += A[2] * B[2];
    return res;
}