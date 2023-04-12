#include <iostream>
#include "Space.hpp"
#include "Parser/Parser.hpp"

int main(const int argc, char* argv[]) {
    Parser pars(argc, argv, true);
    size_t N = pars.getCountConfs();
    std::ifstream fin(pars.getCfgInpFile());
    std::vector<double> E_tr;
    std::vector<double> E_rot;
    std::vector<double> E_vib;
    std::vector<double> alpha;
    std::string tStr;
    
    while (std::getline(fin, tStr)) {
        std::vector<std::string> params = Utils::parseLine3(tStr);
        if (params.size() != 4) {
            std::cerr << "Bad input for molecule" << std::endl;
            return 1;
        }
        E_tr.push_back(stod(params[0]));
        E_rot.push_back(stod(params[1]));
        E_vib.push_back(stod(params[2]));
        alpha.push_back(stod(params[3]));
    }
    fin.close();

    Outer out(pars.getOutFile(), true);

    for (size_t i = 0; i < E_tr.size(); ++i) {
        for (size_t j = 0; j < N; ++j) {
            Space* space = new Space();
            space->setConfig(pars.getCfgFile());

            if (!space->prepareMolecule(E_tr[i], E_rot[i], E_vib[i], alpha[i])) {
                std::cerr << "bad_input" << std::endl;
                delete space;
                continue;
            }

            std::vector<double> startInf = {
                space->molsN2[0].atom[0]->coord[0], space->molsN2[0].atom[0]->coord[1], space->molsN2[0].atom[0]->coord[2],
                space->molsN2[0].atom[1]->coord[0], space->molsN2[0].atom[1]->coord[1], space->molsN2[0].atom[1]->coord[2],
                (space->molsN2[0].atom[0]->vel[0] + space->molsN2[0].atom[0]->vel2[0]) / 2,
                (space->molsN2[0].atom[0]->vel[1] + space->molsN2[0].atom[0]->vel2[1]) / 2,
                (space->molsN2[0].atom[0]->vel[2] + space->molsN2[0].atom[0]->vel2[2]) / 2,
                (space->molsN2[0].atom[1]->vel[0] + space->molsN2[0].atom[1]->vel2[0]) / 2,
                (space->molsN2[0].atom[1]->vel[1] + space->molsN2[0].atom[1]->vel2[1]) / 2,
                (space->molsN2[0].atom[1]->vel[2] + space->molsN2[0].atom[1]->vel2[2]) / 2
            };

            for (auto x : startInf) out.Out << x << " ";
            out.Out << std::endl;

            delete space;
        }
    }

    return 0;
}