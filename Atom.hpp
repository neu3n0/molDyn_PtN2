#pragma once

#include <iostream>
#include <vector>

class Atom {
public:
    Atom();
public:
    double coord[3]{0, 0, 0};
    double vel[3]{0, 0, 0};
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
};

class MoleculeN2 {
public:
	Atom* atom[2]{};
};