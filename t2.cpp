#include <iostream>
#include <vector>
#include <math.h>

const double dt = 1e-6;
const double m = 2.33e-6;
const double bondK = 2243;
const double bondR0 = 1.097;
const double KB = 0.00138;

double calcR(std::vector<double>& one, std::vector<double>& two) {
    double r = 0;
    for (size_t i = 0; i < 3; ++i) r += (one[i] - two[i]) * (one[i] - two[i]);
    return sqrt(r);
}

void coordShift(std::vector<double>& coordOne, std::vector<double>& coordTwo,
                std::vector<double>& velOne, std::vector<double>& velTwo) {
    for (size_t i = 0; i < 3; ++i) {
        coordOne[i] += velOne[i] * dt;
        coordTwo[i] += velTwo[i] * dt;
    }
}

void veldShift(std::vector<double>& powerOne, std::vector<double>& powerTwo,
                std::vector<double>& velOne, std::vector<double>& velTwo,
                std::vector<double>& velOne_2, std::vector<double>& velTwo_2) {
    for (size_t i = 0; i < 3; ++i) {
        velOne_2[i] = velOne[i];
        velTwo_2[i] = velTwo[i];
        velOne[i] += powerOne[i] * dt / m;
        velTwo[i] += powerTwo[i] * dt / m;
    }
}

double KX_F(const double r) {
    return bondK * (bondR0 - r);
}

void powerShift(std::vector<double>& powerOne, std::vector<double>& powerTwo,
                std::vector<double>& coordOne, std::vector<double>& coordTwo) {
    double r = calcR(coordOne, coordTwo);
    double force = KX_F(r);
    for (size_t i = 0; i < 3; ++i) {
        powerOne[i] += (coordOne[i] - coordTwo[i]) / r * force;
        powerTwo[i] -= (coordOne[i] - coordTwo[i]) / r * force;
    }
}


void print_vec(std::vector<double>& vec) {
    for (auto x : vec) {
        std::cout << x << " ";
    }
    std::cout << std::endl;
}

void print_vec_avg(std::vector<double>& vec, std::vector<double>& vec2) {
    for (size_t i = 0; i < 3; ++i) {
        std::cout << (vec[i] + vec2[i]) / 2 << " ";
    }
    std::cout << std::endl;
}

void setNullPower(std::vector<double>& powerOne, std::vector<double>& powerTwo) {
    for (auto& x : powerOne) x = 0;
    for (auto& x : powerTwo) x = 0;
}


int main() {
    std::vector<double> mol_coords_one = {2, - bondR0 / 2 - bondR0 / 40, 2};
    std::vector<double> mol_coords_two = {2, bondR0 / 2 + bondR0 / 40, 2};
    std::vector<double> mol_vel_one(3, 0);
    std::vector<double> mol_vel_two(3, 0);
    std::vector<double> mol_vel_2_one(3, 0);
    std::vector<double> mol_vel_2_two(3, 0);
    std::vector<double> mol_power_one(3, 0);
    std::vector<double> mol_power_two(3, 0);

    powerShift(mol_power_one, mol_power_two, mol_coords_one, mol_coords_two);
    mol_vel_one[1] = mol_power_one[1] * dt / 2 / m;
    mol_vel_two[1] = mol_power_two[1] * dt / 2 / m;
    mol_vel_2_one[1] = mol_power_one[1] * (-dt) / 2 / m;
    mol_vel_2_two[1] = mol_power_two[1] * (-dt) / 2 / m;

    for (size_t step = 0; step < 10; ++step) {
        setNullPower(mol_power_one, mol_power_two);
        coordShift(mol_coords_one, mol_coords_two, mol_vel_one, mol_vel_two);
        powerShift(mol_power_one, mol_power_two, mol_coords_one, mol_coords_two);
        veldShift(mol_power_one, mol_power_two, mol_vel_one, mol_vel_two, mol_vel_2_one, mol_vel_2_two);
        print_vec(mol_coords_two);
        print_vec_avg(mol_vel_two, mol_vel_2_two);
        print_vec(mol_power_two);
        std::cout << "------\n";
    }
    return 0;
}