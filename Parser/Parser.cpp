#include "Parser.h"
#include <stdexcept>

Parser::Parser(const int argc, char* argv[]) {
    for (int i = 1; i < argc; i++) {
        auto key = std::string(argv[i]);
        if (key == "--inp") {
            if (++i == argc) {
                throw std::runtime_error(
                    "Expected value after key -inp is empty");
            }
            inpFile = std::string(argv[i]);
            size_t len = inpFile.length();
            if ((len < 5) || (inpFile.substr(len - 4, 4) != ".inp"))
                throw std::runtime_error(
                    "Value after key -inp have unexpected format");
            continue;
        }
        if (key == "--cfg") {
            if (++i == argc) {
                throw std::runtime_error(
                    "Expected value after key -cfg is empty");
            }
            cfgFile = std::string(argv[i]);
            size_t len = cfgFile.length();
            if ((len < 5) || (cfgFile.substr(len - 4, 4) != ".txt"))
                throw std::runtime_error(
                    "Value after key -cfg have unexpected format");
            continue;
        }
        if (key == "--out") {
            if (++i == argc) {
                throw std::runtime_error(
                    "Expected value after key -inp is empty");
            }
            outFile = std::string(argv[i]);
            continue;
        }
        if (key == "--cfginp") {
            if (++i == argc) {
                throw std::runtime_error(
                    "Expected value after key -inp is empty");
            }
            cfgInpFile = std::string(argv[i]);
            continue;
        }
        throw std::runtime_error("Unknown key \"" + key + "\"");
    }
    if (inpFile.empty()) throw std::runtime_error("INP file not specified");
    if (cfgFile.empty()) throw std::runtime_error("CFG file not specified");
    if (cfgInpFile.empty()) throw std::runtime_error("CFGINP file not specified");
    isLoaded = true;
}

const std::string& Parser::getInpFile() const {
    if (not isLoaded) throw std::runtime_error("CommLine not yet parsed");
    return inpFile;
}

const std::string& Parser::getCfgFile() const {
    if (not isLoaded) throw std::runtime_error("CommLine not yet parsed");
    return cfgFile;
}

const std::string& Parser::getOutFile() const {
    if (not isLoaded) throw std::runtime_error("CommLine not yet parsed");
    return outFile;
}

const std::string& Parser::getCfgInpFile() const {
    if (not isLoaded) throw std::runtime_error("CommLine not yet parsed");
    return cfgInpFile;
}
