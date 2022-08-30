#include "Utils.h"

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