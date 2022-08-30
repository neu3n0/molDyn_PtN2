#pragma once
#include <fstream>
#include <iostream>
#include <string>

class Outer {
public:
    Outer(const std::string& out, const std::string& kinen, const std::string& poten,
        const std::string& fullen, const std::string& rmsVel, const std::string& avel,
            const std::string& temp);
    Outer(const std::string& Out);
    ~Outer();
public:
    std::ofstream fout{};
    std::ofstream Kout{};
    std::ofstream Pout{};
    std::ofstream Fout{};
    std::ofstream Rout{};
    std::ofstream Aout{};
    std::ofstream Tout{};
public:
    std::ofstream Out{};
};