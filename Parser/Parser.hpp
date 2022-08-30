#pragma once

#include <string>

class Parser {
public:
    Parser(const int argc, char* argv[]);
    const std::string& getInpFile() const;
    const std::string& getCfgFile() const;
    const std::string& getOutFile() const;
    const std::string& getCfgInpFile() const;
private:
    std::string inpFile{""};
    std::string cfgFile{""};
    std::string outFile{""};
    std::string cfgInpFile{""};
    bool isLoaded{false};
};
