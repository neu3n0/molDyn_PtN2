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
        if ((!isZ_periodic && i == 2) && (coord[i] < 0 || coord[i] > spaceLength[i])) continue;
        if (coord[i] < 0) coord[i] += spaceLength[i];
        if (coord[i] > spaceLength[i]) coord[i] -= spaceLength[i];
        if (atMolN2 != nullptr && i == 2 && ((coord[i] + atMolN2->coord[i]) / 2) >= hMax + 0.1) {
            was = true;
            return 1;
        }
        // if ((coord[i] < 0 || coord[i] > spaceLength[i]) && !atMolN2) {
        if ((coord[i] < 0 || coord[i] > spaceLength[i])) {
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

bool Atom::checkCell(const size_t i0, const size_t j0, const size_t k0, size_t& i, size_t& j, size_t& k, const double lengthCell){
    i = static_cast<size_t>(coord[0] / lengthCell);
    j = static_cast<size_t>(coord[1] / lengthCell);
    k = static_cast<size_t>(coord[2] / lengthCell);
    if (i != i0 || j != j0 || k != k0) return true;
    return false;
}

void Atom::powerLJ(Atom* atProb, const double* shift) {
    double dr[3] = {0, 0, 0};
    for (size_t i = 0; i < 3; ++i) 
        dr[i] += coord[i] - atProb->coord[i] - shift[i];
    double r2 = 0;
    for (size_t i = 0; i < 3; ++i) 
        r2 += dr[i] * dr[i];
    
    double force = LJ_F_div_r(r2, type, atProb->type);

#ifdef DEBUG_INFO
    double potential = LJ_P(r2, type, atProb->type);
    ljEn += potential;
    atProb->ljEn += potential;
    if (atMolN2 || atProb->atMolN2) {
        ljEn_Pt_N2 += potential;
        atProb->ljEn_Pt_N2 += potential;
    }
#endif

    for (size_t i = 0; i < 3; ++i) {
        double tmp = dr[i] * force;
        power[i] += tmp;
        atProb->power[i] -= tmp;
    }
}

void Atom:: powerKX(Atom* atProb, const double* shift) {
    double r = 0;
    double dr[3] = {0, 0, 0};
    for (size_t i = 0; i < 3; ++i) 
        dr[i] += coord[i] - atProb->coord[i] - shift[i];
    for (size_t i = 0; i < 3; ++i) r += dr[i] * dr[i];
    r = sqrt(r);
    double force = KX_F(r);

    for (size_t i = 0; i < 3; ++i) {
        double tmp = dr[i] / r * force;
        power[i] += tmp;
        atProb->power[i] -= tmp;
    }
}

double Atom::kinEnergy() {
    return m * ((vel[0] + vel2[0]) / 2 * (vel[0] + vel2[0]) / 2 + (vel[1] + vel2[1]) / 2 * (vel[1] + vel2[1]) / 2 + (vel[2] + vel2[2]) / 2 * (vel[2] + vel2[2]) / 2) / 2;
}

double MoleculeN2::calc_length(const std::vector<double>& shift) {
    double r = 0;
    for (size_t i = 0; i < 3; ++i)
        r += (atom[0]->coord[i] - atom[1]->coord[i] - shift[i]) * (atom[0]->coord[i] - atom[1]->coord[i] - shift[i]); 
    return sqrt(r);
}

void MoleculeN2::calcEnVib(const double* shift) {
    double line[3];
    for (size_t i = 0; i < 3; ++i) 
        line[i] = atom[0]->coord[i] - atom[1]->coord[i] - shift[i];
    double velRel[3];
    for (size_t i = 0; i < 3; ++i)
        velRel[i] = (atom[0]->vel[i] + atom[0]->vel2[i]) / 2 - ((atom[0]->vel[i] +  atom[0]->vel2[i]) / 2 + (atom[1]->vel[i] + atom[1]->vel2[i]) / 2) / 2;

    double vRel = Utils::scalProd(velRel, line) / sqrt(Utils::scalProd(line, line));
    double r = 0;
    for (size_t i = 0; i < 3; ++i)
        r += line[i] * line[i]; 
    r = sqrt(r);
    eVib = atom[0]->m * vRel * vRel + KX_P(r);
}

void MoleculeN2::calcEnRot(const double* shift) {
    double vC[3];
    double xC[3];
    double res(0);
    for (size_t i = 0; i < 3; ++i) {
        vC[i] = ((atom[0]->vel[i] + atom[0]->vel2[i]) / 2 + (atom[1]->vel[i] + atom[1]->vel2[i]) / 2) / 2;
        xC[i] = (atom[0]->coord[i] + atom[1]->coord[i] - shift[i]) / 2;
    }
    double e1[3];
    double v[3];
    for (size_t i = 0; i < 3; ++i) {
        // e1[i] = coord[i] - xC[i];
        // if (shift[i] != 0 and std::abs(e1[i]) > 4)
        e1[i] = atom[0]->coord[i] - xC[i] - shift[i];
        v[i] = (atom[0]->vel[i] + atom[0]->vel2[i]) / 2 - vC[i];
    }
    double I = 2 * atom[0]->m * Utils::scalProd(e1, e1);
    std::vector<double> w(3);
    w = Utils::vecProd(e1, v);
    for (auto& x : w)
        x /= Utils::scalProd(e1, e1);

    double W2 = Utils::scalProd(w, w);
    res = I * W2 / 2;
    eRot = res;
}

double MoleculeN2::calcEnTr() {
    double vSq = 0;
    for (size_t in = 0; in < 3; ++in)
        vSq += ((atom[0]->vel[in] + atom[0]->vel2[in]) / 2 + (atom[1]->vel[in] + atom[1]->vel2[in]) / 2) / 2.0 * 
               ((atom[0]->vel[in] + atom[0]->vel2[in]) / 2 + (atom[1]->vel[in] + atom[1]->vel2[in]) / 2) / 2.0;
    eTr = atom[0]->m * vSq;
    return eTr;
}

double LJ_F(const double r, const double r2, const int type1, const int type2) {
    const double r6 = r2 * r2 * r2;
    const double r7 = r6 * r;
    const double r13 = r6 * r6 * r; 
    return 24 * Constants::ljEps[type1][type2] * (2 * Constants::ljSigma12[type1][type2] / r13
        - Constants::ljSigma6[type1][type2] / r7);
}

double LJ_F_div_r(const double r2, const int type1, const int type2) {
    const double r6 = r2 * r2 * r2;
    const double r12 = r6 * r6;
    return 24 * Constants::ljEps[type1][type2] * (2 * Constants::ljSigma12[type1][type2] / r12
        - Constants::ljSigma6[type1][type2] / r6) / r2;
}

double LJ_P(const double r2, const int type1, const int type2) {
    const double r6 = r2 * r2 * r2;
    const double r12 = r6 * r6;
    return 4 * Constants::ljEps[type1][type2] * (Constants::ljSigma12[type1][type2] / r12
        - Constants::ljSigma6[type1][type2] / r6);
}

double KX_F(const double r) {
    return Constants::bondK * (Constants::bondR0 - r);
}

double KX_P(const double r) {
	return Constants::bondK * (Constants::bondR0 - r) * (Constants::bondR0 - r) / 2;
}
