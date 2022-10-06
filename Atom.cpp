#include "Atom.hpp"

Atom::Atom() {}

Atom::Atom(const double* coord0, const double* vel0, const double m0, const int type_) {
    for (size_t i = 0; i < 3; ++i) {
        coord[i] = coord0[i];
        vel[i] = vel0[i];
    }
    m = m0;
    type = type_;
}

Atom::Atom(const std::vector<double>& coord0, const std::vector<double>& vel0, const double m0, const int type_) {
    for (size_t i = 0; i < 3; ++i) {
        coord[i] = coord0[i];
        vel[i] = vel0[i];
    }
    m = m0;
    type = type_;
}

int Atom::coordShift(const double dt, const double* spaceLength, bool isZ_periodic, const double hMax) {
    if (was) return 0;
    for (size_t i = 0; i < 3; ++i) {
        coord[i] += vel[i] * dt;
        if (atMolN2 != nullptr && i == 2 && ((coord[i] + atMolN2->coord[i]) / 2) >= hMax + 0.1) {
            was = true;
            return 1;
        }
        if ((!isZ_periodic && i == 2) && (coord[i] < 0 || coord[i] > spaceLength[i])) continue;
        if (!atMolN2) {
            if (coord[i] < 0) coord[i] += spaceLength[i];
            if (coord[i] > spaceLength[i]) coord[i] -= spaceLength[i];
        }
        else {
            if (i == 2) continue;
            if (coord[i] < 0) coord[i] += spaceLength[i];
            if (coord[i] > spaceLength[i]) coord[i] -= spaceLength[i];
        }
        if ((coord[i] < 0 || coord[i] > spaceLength[i]) && !atMolN2) {
            std::cout << i << " " << coord[i] << " " << vel[i] << " " << dt << " " << type << std::endl;
            throw std::runtime_error("Problem with coordShift");
        }
    }
    was = true;
    return 0;
}

void Atom::velShift(const double dt) {
    for (size_t i = 0; i < 3; ++i) {
        vel2[i] = vel[i];
        vel[i] += power[i] * dt / m;
    }
}

double Atom::kinEnergy() {
    return m * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]) / 2;
}

bool Atom::checkCell(const size_t i0, const size_t j0, const size_t k0, size_t& i, size_t& j, size_t& k, const double lengthCell){
    i = static_cast<size_t>(coord[0] / lengthCell);
    j = static_cast<size_t>(coord[1] / lengthCell);
    k = static_cast<size_t>(coord[2] / lengthCell);
    if (i != i0 || j != j0 || k != k0) return true;
    return false;
}

void Atom::powerLJ(Atom* atProb, const double* shift, bool oneCell) {
    double r = 0;
    for (size_t i = 0; i < 3; ++i) 
        r += (coord[i] - atProb->coord[i] - shift[i]) * (coord[i] - atProb->coord[i] - shift[i]);
    r = sqrt(r);
    double force = LJ_F(r, type, atProb->type);
    double potential = LJ_P(r, type, atProb->type);
    ljEn += potential;
    if (!oneCell)
        atProb->ljEn += potential;
    for (size_t i = 0; i < 3; ++i) {
        power[i] += (coord[i] - atProb->coord[i] - shift[i]) / r * force;
        if (!oneCell)
            atProb->power[i] -= (coord[i] - atProb->coord[i] - shift[i]) / r * force;
    }
} 

double Atom::kinVib(const double* shift) {
    double* line = new double[3];
    for (size_t i = 0; i < 3; ++i) 
        line[i] = coord[i] - atMolN2->coord[i] - shift[i];
    double* velRel = new double[3];
    for (size_t i = 0; i < 3; ++i) 
        velRel[i] = (vel[i] +  vel2[i]) / 2 - ((vel[i] +  vel2[i]) / 2 + (atMolN2->vel[i] + atMolN2->vel2[i]) / 2) / 2;
        // velRel[i] = vel[i] - (vel[i] + atMolN2->vel[i]) / 2;


    double vRel = Utils::scalProd(velRel, line) / sqrt(Utils::scalProd(line, line));
    // std::cout << sqrt(Utils::scalProd(velRel, velRel)) << " " << std::abs(vRel) << std::endl;
    delete[] line;
    delete[] velRel;
    return m * vRel * vRel;
}

void Atom::powerKX(Atom* atProb, const double* shift, bool oneCell) {
    double r = 0;
    for (size_t i = 0; i < 3; ++i) r += (coord[i] - atProb->coord[i] - shift[i]) * (coord[i] - atProb->coord[i] - shift[i]);
    r = sqrt(r);
    double force = KX_F(r);
    double potential = KX_P(r); 
    double kinetical = kinVib(shift);
    eVib += potential;
    eVib += kinetical;
    testVib1 += potential;
    testVib2 += kinetical;
    eRot += calcEnRot(shift);
    if (!oneCell) {
        atProb->eVib += potential;
        atProb->eVib += kinetical;
        atProb->testVib1 += potential;
        atProb->testVib2 += kinetical;
        atProb->eRot += calcEnRot(shift);
    }
    for (size_t i = 0; i < 3; ++i) {
        power[i] += (coord[i] - atProb->coord[i] - shift[i]) / r * force; 
        if (!oneCell)
            atProb->power[i] -= (coord[i] - atProb->coord[i] - shift[i]) / r * force;
    }
} 

double Atom::calcEnRot(const double* shift) {
    double vC[3];
    double xC[3];
    double res(0);
    for (size_t i = 0; i < 3; ++i) {
        vC[i] = (vel[i] + atMolN2->vel[i]) / 2;
        xC[i] = (coord[i] + atMolN2->coord[i] - shift[i]) / 2;
    }
    double e1[3];
    double v[3];
    for (size_t i = 0; i < 3; ++i) {
        e1[i] = coord[i] - xC[i];
        if (shift[i] != 0 and std::abs(e1[i]) > 4)
            e1[i] = coord[i] - xC[i] - shift[i];
        v[i] = vel[i] - vC[i];
    }
    double I = 2 * m * Utils::scalProd(e1, e1);
    std::vector<double> w(3);
    w = Utils::vecProd(e1, v);
    for (auto& x : w)
        x /= Utils::scalProd(e1, e1);

    double W2 = Utils::scalProd(w, w);
    res = I * W2 / 2;
    return res;
}

double Atom::calcEnVib(const double* shift) {
    double res(0);
    double r = 0;
    for (size_t ir = 0; ir < 3; ++ir) r += (coord[ir] - atMolN2->coord[ir] - shift[ir]) * (coord[ir] - atMolN2->coord[ir] - shift[ir]);
    r = sqrt(r);
    res = KX_P(r);
    return res;
}

double LJ_F(const double r, const int type1, const int type2) {
    return 24 * Constants::ljEps[type1][type2] * (2 * Constants::ljSigma12[type1][type2] / pow(r, 13) 
        - Constants::ljSigma6[type1][type2] / pow(r, 7));
}

double LJ_P(const double r, const int type1, const int type2) {
    return 4 * Constants::ljEps[type1][type2] * (Constants::ljSigma12[type1][type2] / pow(r, 12) 
        - Constants::ljSigma6[type1][type2] / pow(r, 6));
}

double KX_F(const double r) {
    return Constants::bondK * (Constants::bondR0 - r);
}

double KX_P(const double r) {
	return Constants::bondK * (Constants::bondR0 - r) * (Constants::bondR0 - r) / 2;
}
