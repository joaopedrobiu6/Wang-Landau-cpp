#include "Lattice.h"
#include "WangLandau.h"
#include "utils.h"

int main()
{
    int L = 14;                 // this value needs to be changed each time we want to simulate a system of different size
    float J = 1.0;
    double f_tol = 1.0e-6, h_tol = 0.95;
    bool NoLog = false;
    int sampleInterval = 2000, f_factor = 2;
    std::string filename = "test_ising.txt"; // GetFilename(L, J, q);

    IsingLattice ising(L, J);
    ising.PrintLattice();
    std::cout << "Energy = " << ising.Ising_Energy() << std::endl;

    std::map<int, double> lng = WangLandauIsing(ising, 10000, f_tol, h_tol, sampleInterval, f_factor, filename);
    return 0;
}