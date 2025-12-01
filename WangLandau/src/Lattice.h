#pragma once

#include "iostream"
#include "vector"
#include "math.h"
#include "map"
#include "fstream"
#include "random"

// only the header file is included here, the implementation is in the lattice.cpp file
std::vector<std::vector<int>> lattice(int L, int q);

// As a class

class PottsLattice
{
public:
    PottsLattice(int L, int q, float J);
    ~PottsLattice() = default;
    double Potts_Energy();
    double Energy_Change(int x, int y, int new_spin);
    std::pair<float, float> Energy_Limit();

    std::vector<std::vector<int>> GetLattice() { return lattice; }
    void PrintLattice();

    std::vector<std::vector<int>> lattice;
    
    // Random number generator
    std::mt19937 rng;

private:
    int L;
    int q;
    double J;
};