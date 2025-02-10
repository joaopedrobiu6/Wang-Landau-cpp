#include "Lattice.h"

// implementation of the PottsLattice class

PottsLattice::PottsLattice(int L, int q, float J) : L(L), q(q), J(J)
{
    std::mt19937 rng(std::time(nullptr));
    std::uniform_int_distribution<int> dist(1, q);

    lattice.reserve(L);

    for (int i = 0; i < L; i++)
    {
        std::vector<int> row;
        row.reserve(L);
        for (int j = 0; j < L; j++)
        {
            row.push_back(dist(rng));
        }
        lattice.push_back(std::move(row));
    }
}

double PottsLattice::Potts_Energy()
{
    int L = lattice.size();
    int E = 0;

    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            int s0 = lattice[i][j];
            
            // Right neighbor (dx=1, dy=0)
            int x = (i + 1) % L;
            int s1 = lattice[x][j];
            E += (s1 == s0);

            // Down neighbor (dx=0, dy=1)
            int y = (j + 1) % L;
            s1 = lattice[i][y];
            E += (s1 == s0);
        }
    }

    return -J * E;
}

std::pair<float, float> PottsLattice::Energy_Limit()
{
    float E_min = -2 * L * L * J;
    float E_max = 0;
    return std::make_pair(E_min, E_max);
};

void PottsLattice::PrintLattice()
{
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            std::cout << lattice[i][j] << ", ";
        }
        std::cout << std::endl;
    }
};

IsingLattice::IsingLattice(int L, float J) : L(L), J(J) {

    std::mt19937 rng(std::time(nullptr));
    std::uniform_int_distribution<int> dist(0, 1);

    lattice.reserve(L);

    for (int i = 0; i < L; i++)
    {
        std::vector<int> row;
        row.reserve(L);
        for (int j = 0; j < L; j++)
        {
            int value = (dist(rng) == 0) ? -1 : 1;
            row.push_back(value);
        }
        lattice.push_back(std::move(row));
    }
};

std::pair<float, float> IsingLattice::Energy_Limit() {
    float E_min = -2 * L * L * J;
    float E_max = 2 * L * L * J;
    return std::make_pair(E_min, E_max);
};

void IsingLattice::PrintLattice()
{
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            std::cout << lattice[i][j] << ", ";
        }
        std::cout << std::endl;
    }
};

double IsingLattice::Ising_Energy() {
    int L = lattice.size();
    int E = 0;

    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            int s0 = lattice[i][j];
            
            // Right neighbor (dx=1, dy=0)
            int x = (i + 1) % L;
            int s1 = lattice[x][j];
            E += s1 * s0;

            // Down neighbor (dx=0, dy=1)
            int y = (j + 1) % L;
            s1 = lattice[i][y];
            E += s1 * s0;
        }
    }

    return -J * E;
};