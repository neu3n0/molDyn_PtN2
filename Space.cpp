#include "Space.hpp"

void Space::setConfig(const std::string& fname) {
    std::ifstream fin;
    fin.open(fname);
    setConfig(fin);
    fin.close();
}

void Space::setConfig(std::istream& inpCfg) {
    inpCfg >> countMolN2 >> countAtPt >> temp >> p >> dt >> lengthCell >> numberCellsX >> numberCellsY >> numberCellsZ;
    spaceLength[0] = lengthCell * numberCellsX;
    spaceLength[1] = lengthCell * numberCellsY;
    spaceLength[2] = lengthCell * numberCellsZ;
}

void Space::init(const std::string& fname) {
    std::ifstream fin;
    fin.open(fname);
    init(fin);
    fin.close();
}

void Space::init(std::istream& inp) {
    double coord[3]{0, 0, 0};
    double vel[3]{0, 0, 0};
    for (size_t c = 0; c < countAtPt; ++c) {
        inp >> coord[0] >> coord[1] >> coord[2] >> vel[0] >> vel[1] >> vel[2];
        int i = static_cast<int>(coord[0] / lengthCell);
        int j = static_cast<int>(coord[1] / lengthCell);
        int k = static_cast<int>(coord[2] / lengthCell);
        if (i >= 0 && j >= 0 && k >= 0 && i < static_cast<int>(numberCellsX) && j < static_cast<int>(numberCellsY) && k < static_cast<int>(numberCellsZ)) {
            Atom* atom = new Atom(coord, vel, MASS_FOR_PT, 1);
            for (size_t ii = 0; ii < 3; ++ii)
                atom->vel2[ii] = atom->vel[ii];
            cells[i][j][k].atoms.push_back(atom);
        }
        else {
            throw std::runtime_error("Incorrect init data for Pt");
        }
    }
}

bool Space::initFromEnergy(const std::string& fname) {
    std::ifstream fin;
    fin.open(fname);
    initFromEnergy(fin);
    fin.close();
    return true;
}

void createVector(std::vector<double>& v, const double length) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    // double fi = M_PI * (2.0 * static_cast<double>(rand()) / RAND_MAX - 1.0);
    double fi = M_PI * (2 * dis(gen) - 1);
    double theta;
    do {
        // theta = M_PI * rand() / RAND_MAX;
        theta = M_PI * dis(gen);
    } while (sin(theta) > 1.0 * static_cast<double>(rand()) / RAND_MAX );
    v[0] = length * sin(theta) * cos(fi);
    v[1] = length * sin(theta) * sin(fi);
    v[2] = length * cos(theta);
}

void createVelAbs(std::vector<double>& v, const double length, const double theta) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    double fi = M_PI * (2 * dis(gen) - 1);
    v[0] = length * sin(theta) * cos(fi);
    v[1] = length * sin(theta) * sin(fi);
    v[2] = length * cos(theta);
}

bool renorm(std::vector<double>& s) {
    double length(0);
    for (size_t i = 0; i < 3; ++i) length += s[i] * s[i];
    length = sqrt(length);
    if (length < 1e-10) return false;
    for (size_t i = 0; i < 3; ++i) s[i] /= length;
    return true;
}

bool Space::initFromParams(const double E_tr, const double E_rot, const double E_vib, const double alpha) {
    //coord
    std::vector<double> e1(3);
    std::vector<double> rC(3);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.1, 34.18);    // need change
    rC[0] = dis(gen);
    rC[1] = dis(gen);
    rC[2] = startHeight;
    std::vector<double> r1(3);
    std::vector<double> r2(3);
    double r_k = std::abs(Constants::bondR0 - sqrt(E_vib * KB * 2 / Constants::bondK));
    createVector(r1, r_k / 2);
    for (size_t i = 0; i < 3; ++i) {
        e1[i] = r1[i];
        r2[i] = -r1[i];
        r1[i] += rC[i];
        r2[i] += rC[i];
    }
    // r1 = {2, -Constants::bondR0 / 2 - Constants::bondR0 / 40  , 2};
    // r2 = {2, Constants::bondR0 / 2 + Constants::bondR0 / 40, 2};
    // r1 = {30.4352, 15.2025, 34.0064};
    // r2 = {30.5754, 14.8615, 32.9936};
    // std::cout << "rC: " <<  (r1[0] + r2[0]) / 2 << ", " << (r1[1] + r2[1]) / 2 << ", " << (r1[2] + r2[2]) / 2 << std::endl; 
    // std::cout << "r1: " << r1[0] << " " << r1[1] << " " << r1[2] << std::endl;
    // std::cout << "r2: " << r2[0] << " " << r2[1] << " " << r2[2] << std::endl;
    // std::cout << "e1: " << e1[0] << " " << e1[1] << " " << e1[2] << std::endl;
    
    double modV = sqrt(E_tr * KB / MASS_FOR_N);     // sqrt(K * Дж/K * 10^(20) / кг * 10^(20))  = sqrt(м^2 / c^2) = (A / t_)  
    double modW = 2.0 / r_k * sqrt(E_rot * KB / MASS_FOR_N);   // 1 / A * sqrt(К * Дж / K  / кг * 10^(20))
    // std::cout << "modules: " << modW << ' ' << modV << std::endl; 

    // vel 
    std::vector<double> vC(3);
    createVelAbs(vC, modV, (180 - alpha) * M_PI / 180);
    // std::cout << "vC: " << vC[0] <<  " " << vC[1] << " " << vC[2] << std::endl;
    // renorm(e1);
    std::vector<double> w_dir(3);
    do {
        std::vector<double> randTmp(3);
        createVector(randTmp, 1);
        w_dir = Utils::vecProd(e1, randTmp);
    } while(!renorm(w_dir));
    
    for (auto& x : w_dir) x *= modW;
    // std::cout << "wdir: " << w_dir[0] << " " << w_dir[1] << " " << w_dir[2] << std::endl;
    // std::cout << "wdirNorm: " << sqrt(w_dir[0] * w_dir[0] + w_dir[1] * w_dir[1] + w_dir[2] * w_dir[2]) << std::endl;
    std::vector<double> v1(3);
    v1 = Utils::vecProd(w_dir, e1);
    // for (auto& x : v1) x *= scalProd(e1, e1);
    std::vector<double> v2(v1);
    for (auto& x : v2) x = -x;
    // std::cout << "v1: " << v1[0] << " " << v1[1] << " " << v1[2] << std::endl;
    std::vector<double> vAbs1(3);
    std::vector<double> vAbs2(3);
    for (size_t i = 0; i < 3; ++i) {
        vAbs1[i] = vC[i] + v1[i];
        vAbs2[i] = vC[i] + v2[i];
    }
    // vAbs1 = {-296.151, -1021.18, -1780.72};
    // vAbs2 = {-233.204, 109.764, -2152.83};
    // std::cout << "vAbs1: " << vAbs1[0] << " " << vAbs1[1] << " " << vAbs1[2] << std::endl;
    // std::cout << "vAbs2: " << vAbs2[0] << " " << vAbs2[1] << " " << vAbs2[2] << std::endl;
    // std::cout << "vC: " << (vAbs1[0] + vAbs2[0]) / 2 << " " << (vAbs1[1] + vAbs2[1]) / 2 << " " << (vAbs1[2] + vAbs2[2]) / 2 << std::endl; 
    // vAbs1 = {0, 0, 0};
    // vAbs2 = {0, 0, 0}; 
    // std::vector<double> wW(3);
    // wW = Utils::vecProd(e1, v1);
    // wW[0] /= Utils::scalProd(e1, e1);
    // wW[1] /= Utils::scalProd(e1, e1);
    // wW[2] /= Utils::scalProd(e1, e1);
    // double I = 2 * MASS_FOR_N * Utils::scalProd(e1, e1);
    // double W2 = Utils::scalProd(wW, wW);
    // std::cout << "EROT = " << I * W2 / 2 / KB << std::endl;
    // std::cout << "alpha: " << vC[2] / sqrt(vC[0] * vC[0] + vC[1] * vC[1]) << std::endl;
    // std::cout << "Etr: " << MASS_FOR_N * (vC[0] * vC[0] + vC[1] * vC[1] + vC[2] * vC[2]) / KB << std::endl;


    double rrr = 0;
    for (size_t i = 0; i < 3; ++i) rrr += (r1[i] - r2[i]) * (r1[i] - r2[i]);
    rrr = sqrt(rrr);
    double force = KX_F(rrr);
    double power[3] = {0, 0, 0};
    for (size_t i = 0; i < 3; ++i)
        power[i] += (r1[i] - r2[i]) / rrr * force;

    MoleculeN2 mol;
    int i = static_cast<int>(r1[0] / lengthCell);
    int j = static_cast<int>(r1[1] / lengthCell);
    int k = static_cast<int>(r1[2] / lengthCell);
    // std::cout << "ind: " << i << "  " << j << "  " << k << std::endl;
    if (i >= 0 && j >= 0 && k >= 0 && i < static_cast<int>(numberCellsX) && j < static_cast<int>(numberCellsY) && k < static_cast<int>(numberCellsZ)) {
        Atom* atom = new Atom(r1, vAbs1, MASS_FOR_N, 0);
        for (size_t ii = 0; ii < 3; ++ii) {
            // atom->vel2[ii] = atom->vel[ii];
            atom->vel2[ii] = atom->vel[ii] + power[ii] * (-dt) / 2 / MASS_FOR_N;
            atom->vel[ii] += power[ii] * dt / 2 / MASS_FOR_N;
        }
        cells[i][j][k].atoms.push_back(atom);
        mol.atom[0] = atom;
    }
    else {
        std::cout << "ind: " << i << " " << j << " " << k << std::endl;
        throw std::runtime_error("Incorrect init data for N2");
    }
    i = static_cast<int>(r2[0] / lengthCell);
    j = static_cast<int>(r2[1] / lengthCell);
    k = static_cast<int>(r2[2] / lengthCell);
    // std::cout << "ind: " << i << "  " << j << "  " << k << std::endl;
    if (i >= 0 && j >= 0 && k >= 0 && i < static_cast<int>(numberCellsX) && j < static_cast<int>(numberCellsY) && k < static_cast<int>(numberCellsZ)) {
        Atom* atom = new Atom(r2, vAbs2, MASS_FOR_N, 0);
        for (size_t ii = 0; ii < 3; ++ii) {
            // atom->vel2[ii] = atom->vel[ii];
            atom->vel2[ii] = atom->vel[ii] - power[ii] * (-dt) / 2 / MASS_FOR_N;
            atom->vel[ii] -= power[ii] * dt / 2 / MASS_FOR_N;
        }
        cells[i][j][k].atoms.push_back(atom);
        mol.atom[1] = atom;
    }
    else {
        std::cout << "ind: " << i << " " << j << " " << k << std::endl;
        throw std::runtime_error("Incorrect init data for N2");
    }
    mol.atom[0]->atMolN2 = mol.atom[1];
    mol.atom[1]->atMolN2 = mol.atom[0];
    molsN2.push_back(mol);

    double eTrr = 0;
    for (size_t in = 0; in < 3; ++in)
        eTrr += pow((molsN2[0].atom[0]->vel[in] + molsN2[0].atom[1]->vel[in]) / 2, 2);
    eTrr *= MASS_FOR_N;
    double ang = std::abs(atan2((molsN2[0].atom[0]->vel[2] + molsN2[0].atom[1]->vel[2]) / 2, sqrt(pow((molsN2[0].atom[0]->vel[1] + molsN2[0].atom[1]->vel[1]) / 2, 2) + pow((molsN2[0].atom[0]->vel[0] + molsN2[0].atom[1]->vel[0]) / 2, 2))) / M_PI * 180);
    ang = 90 - ang;
    if (std::abs(calcVibEn() / KB - E_vib) > 1e-1 || std::abs(calcRotEn() / KB - E_rot) > 1e-1 || std::abs(eTrr / KB - E_tr) > 1e-1 || ang - alpha > 1e-1) {
        std::cerr << "bad initFromParams\n";
        std::cerr << "eVib = " << calcVibEn() / KB << std::endl;
        std::cerr << "eRot = " << calcRotEn() / KB << std::endl;
        std::cerr << "eTr = " << eTrr / KB << std::endl;
        std::cerr << "ang = " << ang << std::endl;
        return false;
    }

    return true; 
}

std::vector<double> Space::initFromCoordsAndVel(
        const std::vector<double>& r1, 
            const std::vector<double>& r2, 
                const std::vector<double>& vAbs1, 
                    const std::vector<double>& vAbs2) {
    double rrr = 0;
    for (size_t i = 0; i < 3; ++i) rrr += (r1[i] - r2[i]) * (r1[i] - r2[i]);
    rrr = sqrt(rrr);
    double force = KX_F(rrr);
    double power[3] = {0, 0, 0};
    for (size_t i = 0; i < 3; ++i)
        power[i] += (r1[i] - r2[i]) / rrr * force;

    MoleculeN2 mol;
    int i = static_cast<int>(r1[0] / lengthCell);
    int j = static_cast<int>(r1[1] / lengthCell);
    int k = static_cast<int>(r1[2] / lengthCell);
    if (i >= 0 && j >= 0 && k >= 0 && i < static_cast<int>(numberCellsX) && j < static_cast<int>(numberCellsY) && k < static_cast<int>(numberCellsZ)) {
        Atom* atom = new Atom(r1, vAbs1, MASS_FOR_N, 0);
        for (size_t ii = 0; ii < 3; ++ii) {
            atom->vel2[ii] = atom->vel[ii] + power[ii] * (-dt) / 2 / MASS_FOR_N;
            atom->vel[ii] += power[ii] * dt / 2 / MASS_FOR_N;
        }
        cells[i][j][k].atoms.push_back(atom);
        mol.atom[0] = atom;
    }
    else {
        std::cout << "ind: " << i << " " << j << " " << k << std::endl;
        throw std::runtime_error("Incorrect init data for N2");
    }
    i = static_cast<int>(r2[0] / lengthCell);
    j = static_cast<int>(r2[1] / lengthCell);
    k = static_cast<int>(r2[2] / lengthCell);
    if (i >= 0 && j >= 0 && k >= 0 && i < static_cast<int>(numberCellsX) && j < static_cast<int>(numberCellsY) && k < static_cast<int>(numberCellsZ)) {
        Atom* atom = new Atom(r2, vAbs2, MASS_FOR_N, 0);
        for (size_t ii = 0; ii < 3; ++ii) {
            atom->vel2[ii] = atom->vel[ii] - power[ii] * (-dt) / 2 / MASS_FOR_N;
            atom->vel[ii] -= power[ii] * dt / 2 / MASS_FOR_N;
        }
        cells[i][j][k].atoms.push_back(atom);
        mol.atom[1] = atom;
    }
    else {
        std::cout << "ind: " << i << " " << j << " " << k << std::endl;
        throw std::runtime_error("Incorrect init data for N2");
    }
    mol.atom[0]->atMolN2 = mol.atom[1];
    mol.atom[1]->atMolN2 = mol.atom[0];
    molsN2.push_back(mol);

    double eTrr = 0;
    for (size_t in = 0; in < 3; ++in)
        eTrr += pow((molsN2[0].atom[0]->vel[in] + molsN2[0].atom[1]->vel[in]) / 2, 2);
    eTrr *= MASS_FOR_N;
    double ang = std::abs(atan2((molsN2[0].atom[0]->vel[2] + molsN2[0].atom[1]->vel[2]) / 2, sqrt(pow((molsN2[0].atom[0]->vel[1] + molsN2[0].atom[1]->vel[1]) / 2, 2) + pow((molsN2[0].atom[0]->vel[0] + molsN2[0].atom[1]->vel[0]) / 2, 2))) / M_PI * 180);
    ang = 90 - ang;

    return {eTrr / KB, calcRotEn() / KB, calcVibEn() / KB, ang}; 
}


bool Space::initFromEnergy(std::istream& inp) {
    double E_tr = 0, E_rot = 0, E_vib = 0, alpha = 0;
    inp >> E_tr >> E_rot >> E_vib >> alpha;
    return initFromParams(E_tr, E_rot, E_vib, alpha);
}

int Space::MDStep() {
    resetChecker();
    if (saveAvg) saveAvgEn();

    int turnOff = 0;
    for (size_t i = 0; i < numberCellsX; ++i)
        for (size_t j = 0; j < numberCellsY; ++j)
            for (size_t k = 0; k < numberCellsZ; ++k)
                for (size_t indAt = 0; indAt < cells[i][j][k].atoms.size(); ++indAt) {
                    turnOff += cells[i][j][k].atoms[indAt]->coordShift(dt, spaceLength, isZ_periodic, hMax);
                    size_t iNew = 0, jNew = 0, kNew = 0;
                    if (cells[i][j][k].atoms[indAt]->checkCell(i, j, k, iNew, jNew, kNew, lengthCell)) {
                        changeCell(i, j, k, iNew, jNew, kNew, indAt);
                        --indAt;
                    }
                }

    SetNullMacro();

    int numberCellsIntX = static_cast<int>(numberCellsX);
    int numberCellsIntY = static_cast<int>(numberCellsY);
    int numberCellsIntZ = static_cast<int>(numberCellsZ);
    for (int i = 0; i < numberCellsIntX; ++i)
        for (int j = 0; j < numberCellsIntY; ++j)
            for (int k = 0; k < numberCellsIntZ; ++k)
                for (size_t indAt = 0; indAt < cells[i][j][k].atoms.size(); ++indAt)
                    for (int i2 = -1; i2 < 2; ++i2)
                        for (int j2 = -1; j2 < 2; ++j2)
                            for (int k2 = 0; k2 < 2; ++k2) {
                                if (k2 == 0 && ((i2 == -1 && j2 == -1) || (i2 == -1 && j2 == 0) || (i2 == 0 && j2 == -1) || (i2 == 1 && j2 == -1))) continue;
                                int i3 = i + i2, j3 = j + j2, k3 = k + k2;
                                double shift[3] = { 0, 0, 0 };
                                if (i3 >= numberCellsIntX) { i3 = 0;    shift[0] = spaceLength[0]; }
                                if (j3 >= numberCellsIntY) { j3 = 0;    shift[1] = spaceLength[1]; }
                                if (k3 >= numberCellsIntZ) { k3 = 0;    shift[2] = spaceLength[2]; }
                                if (i3 < 0) { i3 = numberCellsX - 1;     shift[0] = -spaceLength[0]; }
                                if (j3 < 0) { j3 = numberCellsY - 1;     shift[1] = -spaceLength[1]; }
                                if (k3 < 0) { k3 = numberCellsZ - 1;     shift[2] = -spaceLength[2]; }
                                for (size_t indAt2 = 0; indAt2 < cells[i3][j3][k3].atoms.size(); ++indAt2) {
                                    if (i3 == i && j3 == j && k3 == k && indAt >= indAt2) continue;
                                    if (cells[i][j][k].atoms[indAt]->atMolN2 != cells[i3][j3][k3].atoms[indAt2])
                                        cells[i][j][k].atoms[indAt]->powerLJ(cells[i3][j3][k3].atoms[indAt2], shift);
                                    else
                                        cells[i][j][k].atoms[indAt]->powerKX(cells[i3][j3][k3].atoms[indAt2], shift);
                                }
                            }


    for (size_t i = 0; i < numberCellsX; ++i) 
        for (size_t j = 0; j < numberCellsY; ++j)
            for (size_t k = 0; k < numberCellsZ; ++k) 
                for (size_t indAt = 0; indAt < cells[i][j][k].atoms.size(); ++indAt)
                    cells[i][j][k].atoms[indAt]->velShift(dt);
    
    
    for (auto& mol : molsN2) {
        double kin = mol.atom[0]->kinVib();
        mol.atom[0]->testVib2 += kin;
        mol.atom[1]->testVib2 += kin;
        mol.atom[0]->eVib += kin;
        mol.atom[1]->eVib += kin;
    }

    if (turnOff) return 1;
    return 0;
}

void Space::changeCell(size_t i, size_t j, size_t k, size_t iNew, size_t jNew, size_t kNew, size_t indAt) {
    cells[iNew][jNew][kNew].atoms.push_back(cells[i][j][k].atoms[indAt]);
    std::swap(cells[i][j][k].atoms[indAt], cells[i][j][k].atoms[cells[i][j][k].atoms.size() - 1]);
    cells[i][j][k].atoms.pop_back();
}

void Space::SetNullMacro() {
    kinEn = 0;
    potEn = 0;
    vibEn = 0;
    energy = 0;
    for (size_t i = 0; i < numberCellsX; ++i) 
        for (size_t j = 0; j < numberCellsY; ++j)
            for (size_t k = 0; k < numberCellsZ; ++k) 
                for (size_t indAt = 0; indAt < cells[i][j][k].atoms.size(); ++indAt) {
                    for (size_t ind = 0; ind < 3; ++ind) 
                        cells[i][j][k].atoms[indAt]->power[ind] = 0;

                    cells[i][j][k].atoms[indAt]->kinEn = 0;
                    cells[i][j][k].atoms[indAt]->kxEn = 0;
                    cells[i][j][k].atoms[indAt]->ljEn = 0;
                    cells[i][j][k].atoms[indAt]->eRot = 0;
                    cells[i][j][k].atoms[indAt]->eVib = 0;
                    cells[i][j][k].atoms[indAt]->testVib1 = 0;
                    cells[i][j][k].atoms[indAt]->testVib2 = 0;                    
                }
}

void Space::resetChecker() {
    for (size_t i = 0; i < numberCellsX; ++i)
        for (size_t j = 0; j < numberCellsY; ++j)
            for (size_t k = 0; k < numberCellsZ; ++k)
                for (size_t indAt = 0; indAt < cells[i][j][k].atoms.size(); ++indAt)
                    cells[i][j][k].atoms[indAt]->was = false;    
}


void Space::saveAvgEn() {
    if (indAvg < 14) {
        avgVibEn += molsN2[0].atom[0]->eVib / KB;
        avgRotEn += molsN2[0].atom[0]->eRot / KB;
    }
    ++indAvg;
}

void Space::printConfig() const {
    std::cout << countMolN2 << std::endl <<
    countAtPt << std::endl << temp << std::endl <<
    p << std::endl << dt << std::endl << 
    spaceLength[0] << " " << spaceLength[1] << " " << spaceLength[2] << std::endl  << lengthCell << std::endl  <<
    numberCellsX << std::endl << numberCellsY << std::endl <<
    numberCellsZ << std::endl;
}

void Space::printInit() const {
    for (size_t i = 0; i < numberCellsX; ++i) 
        for (size_t j = 0; j < numberCellsY; ++j)
            for (size_t k = 0; k < numberCellsZ; ++k) 
                for (size_t indAt = 0; indAt < cells[i][j][k].atoms.size(); ++indAt) {
                    std::cout << i << " " << j << " " << k << " " << indAt << std::endl
                        << cells[i][j][k].atoms[indAt]->coord[0] << " " << cells[i][j][k].atoms[indAt]->coord[1] << " " 
                        << cells[i][j][k].atoms[indAt]->coord[2] << " " << cells[i][j][k].atoms[indAt]->vel[0] << " "
                        << cells[i][j][k].atoms[indAt]->vel[1] << " " << cells[i][j][k].atoms[indAt]->vel[2] << " " 
                        << cells[i][j][k].atoms[indAt]->m << " " << cells[i][j][k].atoms[indAt]->atMolN2 << std::endl; 
                }
}

void Space::printInfo() const {
    for (size_t i = 0; i < numberCellsX; ++i) 
        for (size_t j = 0; j < numberCellsY; ++j)
            for (size_t k = 0; k < numberCellsZ; ++k) 
                for (size_t indAt = 0; indAt < cells[i][j][k].atoms.size(); ++indAt) {
                    std::cout << i << " " << j << " " << k << " " << indAt << std::endl
                        << "coord: " << cells[i][j][k].atoms[indAt]->coord[0] << " " << cells[i][j][k].atoms[indAt]->coord[1] << " " << cells[i][j][k].atoms[indAt]->coord[2] << "; " 
                        << "vel: " << cells[i][j][k].atoms[indAt]->vel[0] << " " << cells[i][j][k].atoms[indAt]->vel[1] << " " << cells[i][j][k].atoms[indAt]->vel[2] << "; " 
                        << "power: " << cells[i][j][k].atoms[indAt]->power[0] << " " << cells[i][j][k].atoms[indAt]->power[1] << " " << cells[i][j][k].atoms[indAt]->power[2] << "; "
                        << "m: " << cells[i][j][k].atoms[indAt]->m << "; pointer: " << cells[i][j][k].atoms[indAt]->atMolN2 << "; type: " << cells[i][j][k].atoms[indAt]->type << std::endl; 
                }
}

void Space::printMol() const {
    size_t numMol = 1;
    for (const auto& m : molsN2) {
        std::cout << "mol #" << numMol << std::endl;
        std::cout << m.atom[0] << "; coord: " << m.atom[0]->coord[0] << " " << m.atom[0]->coord[1] << " " << m.atom[0]->coord[2] << "; " 
            << "vel: " << m.atom[0]->vel[0] << " " << m.atom[0]->vel[1] << " " << m.atom[0]->vel[2] << "; " 
            << "power: " << m.atom[0]->power[0] << " " << m.atom[0]->power[1] << " " << m.atom[0]->power[2] << "; "
            << "m: " << m.atom[0]->m << "; pointer: " << m.atom[0]->atMolN2 << "; type: " << m.atom[0]->type << std::endl; 
        std::cout << m.atom[1] << "; coord: " << m.atom[1]->coord[0] << " " << m.atom[1]->coord[1] << " " << m.atom[1]->coord[2] << "; " 
            << "vel: " << m.atom[1]->vel[0] << " " << m.atom[1]->vel[1] << " " << m.atom[1]->vel[2] << "; " 
            << "power: " << m.atom[1]->power[0] << " " << m.atom[1]->power[1] << " " << m.atom[1]->power[2] << "; "
            << "m: " << m.atom[1]->m << "; pointer: " << m.atom[1]->atMolN2 << "; type: " << m.atom[1]->type << std::endl; 
        ++numMol;
    }
}

int Space::writeVTK(const std::string& name) {
    size_t t = 0;
    for (size_t i = 0; i < numberCellsX; ++i) {
        for (size_t j = 0; j < numberCellsY; ++j) 
            for (size_t k = 0; k < numberCellsZ; ++k) {
                for (size_t indAt = 0; indAt < cells[i][j][k].atoms.size(); ++indAt) {
                    if (!cells[i][j][k].atoms[indAt]->atMolN2)
                        t += 1;
                }
                // t += cells[i][j][k].atoms.size();
            }
    }
    t += 2 * countMolN2;
	char fname[100];
	sprintf(fname, "vtk/%s_%010d.vtk", name.c_str(), vtkNum);
	std::ofstream fout(fname);
    fout << "# vtk DataFile Version 2.0\nMolecules states\nASCII\nDATASET POLYDATA\nPOINTS " << t << " float" << std::endl;
	for (size_t i = 0; i < numberCellsX; i++)
		for (size_t j = 0; j < numberCellsY; j++)
			for (size_t k = 0; k < numberCellsZ; k++)
				for (size_t n = 0; n < cells[i][j][k].atoms.size(); n++) {
                    if (cells[i][j][k].atoms[n]->atMolN2) continue;
                    fout << cells[i][j][k].atoms[n]->coord[0] << " " << cells[i][j][k].atoms[n]->coord[1] <<
                        " " << cells[i][j][k].atoms[n]->coord[2] << std::endl;
                }
                    

    for (const auto& m : molsN2) {
        fout << m.atom[0]->coord[0] << " " << m.atom[0]->coord[1] <<
                        " " << m.atom[0]->coord[2] << std::endl;
        fout << m.atom[1]->coord[0] << " " << m.atom[1]->coord[1] <<
                        " " << m.atom[1]->coord[2] << std::endl;
    }
                        
	fout << "POINT_DATA " << t << std::endl << "SCALARS MoleculeType float 1" << 
        std::endl << "LOOKUP_TABLE default" << std::endl;

	for (size_t i = 0; i < numberCellsX; i++)
		for (size_t j = 0; j < numberCellsY; j++)
			for (size_t k = 0; k < numberCellsZ; k++)
				for (size_t n = 0; n < cells[i][j][k].atoms.size(); n++) {
					if (cells[i][j][k].atoms[n]->atMolN2) continue;
                    fout << cells[i][j][k].atoms[n]->type + 1 << " ";
				}

    for (const auto& m : molsN2) fout << m.atom[0]->type + 1 << " " << m.atom[1]->type + 1 << " ";
	++vtkNum;
	fout.close();
	return 0;
}

Cell::~Cell() {
    while(!atoms.empty()) {
        delete atoms.back();
        atoms.pop_back();
    }
}


void Space::saveStruct() {
    std::ofstream fout("output/struct.txt");
    for (size_t i = 0; i < numberCellsX; ++i)
        for (size_t j = 0; j < numberCellsY; ++j)
            for (size_t k = 0; k < numberCellsZ; ++k) 
                for (size_t indAt = 0; indAt < cells[i][j][k].atoms.size(); ++indAt) {
                    for (size_t n = 0; n < 3; ++n) {
                        fout << cells[i][j][k].atoms[indAt]->coord[n] << " ";
                    }
                    for (size_t n = 0; n < 3; ++n) {
                        fout << cells[i][j][k].atoms[indAt]->vel[n] << " ";
                    }
                    fout << std::endl;
                }
    fout.close();
}

double Space::calcRotEn() {
    double vC[3];
    double xC[3];
    double res(0);
    if (!molsN2.empty()) {
        for (size_t i = 0; i < 3; ++i) {
            vC[i] = (molsN2[0].atom[0]->vel[i] + molsN2[0].atom[1]->vel[i]) / 2;
            xC[i] = (molsN2[0].atom[0]->coord[i] + molsN2[0].atom[1]->coord[i]) / 2;
        }
        double e1[3];
        double v[3];
        for (size_t i = 0; i < 3; ++i) {
            e1[i] = molsN2[0].atom[0]->coord[i] - xC[i];
            v[i] = molsN2[0].atom[0]->vel[i] - vC[i];
        }
        double I = 2 * molsN2[0].atom[0]->m * Utils::scalProd(e1, e1);
        std::vector<double> w(3);
        w = Utils::vecProd(e1, v);
        for (auto& x : w)
            x /= Utils::scalProd(e1, e1);

        double W2 = Utils::scalProd(w, w);
        res = I * W2 / 2;
    }
    return res;
}

double Space::calcKinVib() {
    double* line = new double[3];
    for (size_t i = 0; i < 3; ++i) 
        line[i] = molsN2[0].atom[0]->coord[i] - molsN2[0].atom[1]->coord[i];
    double* velRel = new double[3];
    for (size_t i = 0; i < 3; ++i) 
        velRel[i] = (molsN2[0].atom[0]->vel[i] +molsN2[0].atom[0]->vel2[i]) / 2 - ((molsN2[0].atom[0]->vel[i] + molsN2[0].atom[0]->vel2[i]) / 2 + (molsN2[0].atom[1]->vel[i] + molsN2[0].atom[1]->vel2[i]) / 2) / 2;

    double vRel = Utils::scalProd(velRel, line) / sqrt(Utils::scalProd(line, line));
    delete[] line;
    delete[] velRel;
    return molsN2[0].atom[0]->m * vRel * vRel;
}

double Space::calcVibEn() {
    double res{0};
    if (!molsN2.empty()) {
        double r = 0;
        for (size_t ir = 0; ir < 3; ++ir) r += (molsN2[0].atom[0]->coord[ir] - molsN2[0].atom[1]->coord[ir]) * (molsN2[0].atom[0]->coord[ir] - molsN2[0].atom[1]->coord[ir]);
        r = sqrt(r);
        double potential = KX_P(r); 
        double kinetical = calcKinVib();
        res = potential + kinetical;
    }
    return res;
}

void saveInfo(Outer& out, const double E_tr, const double E_rot, const double E_vib,
    const double alpha, std::vector<double>& vel, const double eTr, const double eVib, 
        const double eRot, size_t n, double time, size_t steps) {
    
    out.Out << E_tr << " " << E_rot << " " << E_vib << " " << alpha << " " 
        << vel[0] / n << " " << vel[1] / n << " " << vel[2] / n << " " << eTr / n
            << " " << eRot / n << " " << eVib / n << " " << steps / n << " " << time << std::endl;
}

void Space::getEnergy(Outer& out, int step) {
	out.fout << "Full energy = " << energy << "\t" << "Kinetic energy = " 
		<< kinEn << "\t" << "Potential energy = " << potEn << std::endl;
	out.Kout << step << "\t" <<  kinEn << std::endl;
	out.Pout << step << "\t" << potEn << std::endl;
	out.Fout << step << "\t" << energy << std::endl;
	out.Tout << step << "\t" << tempr() << std::endl;
	out.Rout << step << "\t" << rmsVel() << std::endl;
	std::vector<double> av_vel = averVel();
	out.Aout << step << " " << av_vel[0] << " " << av_vel[1] << " " << av_vel[2] << " " << std::endl;
}

double Space::rmsVel() {
	return sqrt(3 * KB * temp / MASS_FOR_PT);
}

double Space::tempr() {
	double mV2 = 0;
	double a = 0;
	double res = 0;
	for (size_t i = 0; i < numberCellsX; ++i) {
		for (size_t j = 0; j < numberCellsY; ++j) {
			for (size_t k = 0; k < numberCellsZ; ++k) {
				for (size_t l = 0; l < cells[i][j][k].atoms.size(); ++l) {
					Atom* at = cells[i][j][k].atoms[l];
					a = at->m * (pow(at->vel[0], 2) + pow(at->vel[1], 2) + pow(at->vel[2], 2));
					mV2 += a;
				}
			}
		}
	}
	res = mV2 / ((countAtPt + countMolN2) * (3 * KB));
	temp = res;
	return res;
}

std::vector<double> Space::averVel() {
	double v1 = 0;
	double v2 = 0;
	double v3 = 0;
	for (size_t i = 0; i < numberCellsX; ++i) {
		for (size_t j = 0; j < numberCellsY; ++j) {
			for (size_t k = 0; k < numberCellsZ; ++k) {
				for (size_t l = 0; l < cells[i][j][k].atoms.size(); ++l) {
					Atom* at = cells[i][j][k].atoms[l];
					v1 += at->vel[0];
					v2 += at->vel[1];
					v3 += at->vel[2];
				}
			}
		}
	}
	std::vector<double> av_vel;
	av_vel.push_back(v1 / (countAtPt + countMolN2));
	av_vel.push_back(v2 / (countAtPt + countMolN2));
	av_vel.push_back(v3 / (countAtPt + countMolN2));
	return av_vel;
}