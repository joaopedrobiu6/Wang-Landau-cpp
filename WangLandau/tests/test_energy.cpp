#include "../src/Lattice.h"
#include <iostream>
#include <cassert>
#include <cmath>

int main() {
    int L = 10;
    int q = 4;
    float J = 1.0;
    
    PottsLattice lat(L, q, J);
    
    // Test 1: Energy Consistency
    double E_full = lat.Potts_Energy();
    
    // Pick a random site
    int x = 5;
    int y = 5;
    int s0 = lat.lattice[x][y];
    int s1 = (s0 % q) + 1; // Different spin
    
    double delta_E = lat.Energy_Change(x, y, s1);
    
    // Manually update lattice
    lat.lattice[x][y] = s1;
    double E_new_full = lat.Potts_Energy();
    
    std::cout << "Old Energy: " << E_full << std::endl;
    std::cout << "Predicted Delta: " << delta_E << std::endl;
    std::cout << "New Energy (Full Calc): " << E_new_full << std::endl;
    std::cout << "Actual Delta: " << E_new_full - E_full << std::endl;
    
    if (std::abs((E_new_full - E_full) - delta_E) < 1e-6) {
        std::cout << "Test Passed: Energy update is consistent." << std::endl;
    } else {
        std::cout << "Test Failed: Energy update mismatch." << std::endl;
        return 1;
    }
    
    return 0;
}
