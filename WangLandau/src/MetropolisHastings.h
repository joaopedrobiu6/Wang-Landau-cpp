#include "iostream"
#include "vector"
#include "math.h"
#include "random"
#include "map"
#include "algorithm"
#include "fstream"
#include "Lattice.h"
#include "utils.h"

void MetropolisStep(double Temperature)
{    
    std::mt19937 rng(std::time(nullptr));
    std::uniform_int_distribution<int> coord_dist(0, L - 1);
    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
}   