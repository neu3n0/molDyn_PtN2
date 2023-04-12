#include "Space.hpp"
#include "Parser/Parser.hpp"
#include <chrono>
#include <omp.h>

const size_t MAXSTEPS = 200000;

int main(const int argc, char* argv[]) {
    Parser pars(argc, argv);
    std::vector<std::vector<double>> coord1, coord2, vel1, vel2;
    if (!Space::getMoleculesParams(pars.getCfgInpFile(), coord1, coord2, vel1, vel2)) return 1;
    Outer outT("out");
    Outer outS("start_inf");
    auto wholeStart = std::chrono::steady_clock::now();

    std::vector<double> config = Space::getCfg(pars.getCfgFile());
    std::vector<std::vector<double>> PlatinumSurf = Space::getPlatinumSurf(pars.getInpFile());

#ifndef DEBUG_INFO
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < coord1.size(); ++i) {
        Space* space = new Space();
        space->setConfig(config);
        space->init(PlatinumSurf);
        std::vector<double> macro = space->initFromCoordsAndVel(coord1[i], coord2[i], vel1[i], vel2[i]);
        std::vector<double> startInf = space->getStartMol();

        size_t step = 0;
        auto start = std::chrono::steady_clock::now();
#ifdef DEBUG_INFO
        std::cout << "conf #" << i << std::endl;
#endif
        while (!space->MDStep() && step < MAXSTEPS) {
            // if (pars.getOutFile() != "" && step % 100 == 0) space->writeVTK(pars.getOutFile() + "_" + std::to_string(i));
            ++step;
        }
        double timeT = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()) / 1000;
        space->saveAndStepsForAvg();
        Space::saveInfo(outT, outS, i, macro[0], macro[1], macro[2], macro[3], space->vMol, space->eTr, space->avgVibEn, space->avgRotEn, 1, timeT, step, startInf, step >= MAXSTEPS ? true : false);

        delete space;
    }

    auto tmp = std::chrono::steady_clock::now();
    std::cout << "full duration of calculation: " << static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(tmp - wholeStart).count()) / 1000 << " sec" << std::endl;
    
    return 0;
}