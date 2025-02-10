#pragma once

#include "iostream"
#include "vector"
#include "math.h"
#include "map"
#include "fstream"
#include "random"

// only the header file is included here, the implementation is in the lattice.cpp file
std::vector<std::vector<int>> lattice(int L, int q);

class PottsLattice
{
public:
    PottsLattice(int L, int q, float J);
    ~PottsLattice() = default;
    double Potts_Energy();
    std::pair<float, float> Energy_Limit();

    std::vector<std::vector<int>> GetLattice() { return lattice; }
    void PrintLattice();

    std::vector<std::vector<int>> lattice;

private:
    int L;
    int q;
    double J;
};

class IsingLattice
{
public:
    IsingLattice(int L, float J);
    ~IsingLattice() = default;
    double Ising_Energy();
    std::pair<float, float> Energy_Limit();

    std::vector<std::vector<int>> GetLattice() { return lattice; }
    void PrintLattice();

    std::vector<std::vector<int>> lattice;

private:
    int L;
    double J;

};