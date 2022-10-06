#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <ctime>
#include <fstream>
#include "Utils/Utils.hpp"

const double KB = 0.00138;

const double startHeight = 33.5;

struct Constants {
    inline static double ljSigma12[2][2]{{pow(3.22, 12), pow(3.82, 12)}, {pow(3.82, 12), pow(2.523, 12)}};
    inline static double ljSigma6[2][2]{{pow(3.22, 6), pow(3.82, 6)}, {pow(3.82, 6), pow(2.523, 6)}};
    inline static double ljEps[2][2]{{29.41 * KB, 67.39 * KB}, {67.39 * KB, 3771.5 * KB}};
    inline static double bondK = 2243;
    inline static double bondR0 = 1.097; // Ангстремы
};

class Atom {
public:
    bool was{false};
    Atom();
	Atom(const double* coord0, const double* vel0, const double m0, const int type_);
    Atom(const std::vector<double>& coord0, const std::vector<double>& vel0, const double m0, const int type_);
    int coordShift(const double dt, const double* spaceLength, bool isZ_periodic, const double hMax);
    void velShift(const double dt);
    bool checkCell(const size_t i0, const size_t j0, const size_t k0, size_t& i, size_t& j, size_t& k, const double lengthCell);
    void powerLJ(Atom* atProb, const double* shift, bool);
    void powerKX(Atom* atProb, const double* shift, bool);
    double kinVib(const double*);
    double kinEnergy();
    double testVib1 = 0;
    double testVib2 = 0;
public:
    double coord[3]{0, 0, 0};
    double vel[3]{0, 0, 0};
    double vel2[3]{0, 0, 0};
    double power[3]{0, 0, 0};
    double m{0};
    Atom* atMolN2{nullptr};
    int type{-1};
public:
    double kinEn{0};
    double ljEn{0};
    double kxEn{0};
public:
    double eRot{0};
    double eVib{0};
    double calcEnRot(const double* shift);
    double calcEnVib(const double* shift, const double r);
};

class MoleculeN2 {
public:
	Atom* atom[2]{};
};

double LJ_F(const double r, const int type1, const int type2);
double LJ_P(const double r, const int type1, const int type2);
double KX_F(const double r);
double KX_P(const double r);