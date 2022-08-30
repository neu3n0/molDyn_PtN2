#include "Out.hpp"

Outer::Outer(const std::string& out, const std::string& kinen, const std::string& poten,
        const std::string& fullen, const std::string& rmsVel, const std::string& avel,
            const std::string& temp) {
    fout.open("output/macro/" + out);
    Kout.open("output/macro/" + kinen);
    Pout.open("output/macro/" + poten);
    Fout.open("output/macro/" + fullen);
    Rout.open("output/macro/" + rmsVel);
    Aout.open("output/macro/" + avel);
    Tout.open("output/macro/" + temp);
}

Outer::Outer(const std::string& out) {
    Out.open("output/calcs" + out);
}

Outer::~Outer() {
    if (fout.is_open()) fout.close();
    if (Kout.is_open()) Kout.close();
    if (Pout.is_open()) Pout.close();
    if (Fout.is_open()) Fout.close();
    if (Rout.is_open()) Rout.close();
    if (Aout.is_open()) Aout.close();
    if (Tout.is_open()) Tout.close();
    if (Out.is_open()) Out.close();
}