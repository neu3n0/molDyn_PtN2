#pragma once

#include <string>

class Parser {
public:
    Parser(const int argc, char* argv[]);
    Parser(const int argc, char* argv[], bool flag);
    const std::string& getInpFile() const;
    const int& getCountConfs() const;
    const std::string& getCfgFile() const;
    const std::string& getOutFile() const;
    const std::string& getCfgInpFile() const;
private:
    std::string inpFile{""};
    std::string cfgFile{""};
    std::string outFile{""};
    std::string cfgInpFile{""};
    bool isLoaded{false};
    int N{0};
};
