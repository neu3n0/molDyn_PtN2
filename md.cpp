#include "Space.hpp"
#include "Parser/Parser.hpp"
#include <chrono>
#include <omp.h>

const size_t MAXSTEPS = 300000;

int main(const int argc, char* argv[]) {
    Parser pars(argc, argv);

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

    Outer out("N.out");
    
    auto wholeStart = std::chrono::steady_clock::now();

    size_t N = 1;
    for (size_t i = 0; i < E_tr.size(); ++i) {
        size_t countBadCalcs = 0;
        // std::cout << E_tr[i] << "  " << E_rot[i] << "  " << E_vib[i] << " " << alpha[i] << std::endl;
        std::vector<double> vel(3, 0);
        double eRot = 0, eVib = 0, eTr = 0;
        size_t steps = 0;
        auto confStart = std::chrono::steady_clock::now();
        Outer outT("out_" + std::to_string(i));
        Outer outS("start_inf_" + std::to_string(i));
        double badTime = 0;
        #pragma omp parallel for shared(badTime)
        for (size_t j = 0; j < N; ++j) {
            // std::cout << "----- N = " << j + 1 << " -----" << std::endl;
            Space* space = new Space();
            srand(time(nullptr));
            space->setConfig(pars.getCfgFile());
            space->init(pars.getInpFile());
            if (!space->initFromParams(E_tr[i], E_rot[i], E_vib[i], alpha[i])) {
                delete space;
                continue;
            }

            // need for debug
            std::vector<double> startInf = {
                space->molsN2[0].atom[0]->coord[0], space->molsN2[0].atom[0]->coord[1], space->molsN2[0].atom[0]->coord[2],
                space->molsN2[0].atom[1]->coord[0], space->molsN2[0].atom[1]->coord[1], space->molsN2[0].atom[1]->coord[2],
                space->molsN2[0].atom[0]->vel[0], space->molsN2[0].atom[0]->vel[1], space->molsN2[0].atom[0]->vel[2], 
                space->molsN2[0].atom[1]->vel[0], space->molsN2[0].atom[1]->vel[1], space->molsN2[0].atom[1]->vel[2]
            };

            size_t step = 0;
            auto start = std::chrono::steady_clock::now();
            space->vtkNum = 0;
            space->saveAvg = false;

            while (!space->MDStep() && step < 1000) {
                // if (step % 100 == 0) {
                //     auto tmp = std::chrono::steady_clock::now();
                //     std::cout << "step: " << step << "  |  " << static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(tmp - start).count()) / 1000 << " sec" << std::endl;
                // }
                if (step % 100 == 0) space->writeVTK(pars.getOutFile() + "_" + std::to_string(i) + "_" + std::to_string(j));
                ++step;
            }
            space->writeVTK(pars.getOutFile() + "_" + std::to_string(i) + "_" + std::to_string(j));
            double eV = 0, eR = 0, eT = 0;
            std::vector<double> vTmp(3, 0);
            if (!space->molsN2.empty()) {
                for (size_t in = 0; in < 3; ++in)
                    vTmp[in] = (space->molsN2[0].atom[0]->vel[in] + space->molsN2[0].atom[1]->vel[in]) / 2.0;
                for (size_t in = 0; in < 3; ++in) eT += vTmp[in] * vTmp[in];
                eT *= MASS_FOR_N / KB; 
            }
            space->saveAvg = true;
            for (size_t gr = 0; gr < 150; ++gr) space->MDStep();
            eV = space->avgVibEn / 144;
            eR = space->avgRotEn / 144;
            auto tmp = std::chrono::steady_clock::now();
            double timeT = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(tmp - start).count()) / 1000;
            // std::cout << "duration of calculation: " << timeT << " sec" << std::endl;
            
            if (step != MAXSTEPS) {
                for (size_t in = 0; in < 3; ++in) vel[in] += vTmp[in];
                eVib += eV;
                eRot += eR;
                eTr += eT;
                steps += step;
                saveInfo(outT, E_tr[i], E_rot[i], E_vib[i], alpha[i], vTmp, eT, eV, eR, 1, timeT, step);
                for (auto x : startInf) outS.Out << x << " ";
                    outS.Out << std::endl;
            }
            else { 
                ++countBadCalcs;
                if (countBadCalcs < 8) 
                    --j;
                badTime += timeT;
                std::cout << "BAD_CALC: " << E_tr[i] << " " << E_rot[i] << " " << E_vib[i] << " " << alpha[i] << std::endl;
                for (auto x : startInf) std::cout << x << " ";
                std::cout << std::endl;
            }
            // std::cout << vTmp[0] << " " << vTmp[1] << " " << vTmp[2] 
            //     << "  |  E_Tr = " << eT << " E_rot = " << eR << " E_vib = " << eV << "  |  " << step << " " << timeT << std::endl;
            delete space;
        }
        
        auto tmp = std::chrono::steady_clock::now();
        double confTime = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(tmp - confStart).count()) / 1000;
        saveInfo(out, E_tr[i], E_rot[i], E_vib[i], alpha[i], vel, eTr, eVib, eRot, N, confTime - badTime, steps);
        // std::cout << vel[0] / N << " " << vel[1] / N << " " << vel[2] / N 
        //         << "  |  E_tr = " << eTr / N << " E_rot = " << eRot / N << " E_vib = " << eVib / N 
        //             << "  |  " << steps << " " << confTime << std::endl;

    }

    auto tmp = std::chrono::steady_clock::now();
    std::cout << "full duration of calculation: " << static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(tmp - wholeStart).count()) / 1000 << " sec" << std::endl;
    
    return 0;
}   