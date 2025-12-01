#include "Lattice.h"

// implementation of the PottsLattice class

PottsLattice::PottsLattice(int L, int q, float J) : L(L), q(q), J(J)
{
    std::random_device rd;
    rng.seed(rd());
    std::uniform_int_distribution<int> dist(1, q);

    for (int i = 0; i < L; i++)
    {
        std::vector<int> row;
        for (int j = 0; j < L; j++)
        {
            row.push_back(dist(rng));
        }
        lattice.push_back(row);
    }
}

double PottsLattice::Potts_Energy()
{
    int L = lattice.size();
    double E = 0.0;

    // maybe we can iterate over only half of the lattice to avoid double counting
    for (int i = 0; i < L; ++i)
    {
        for (int j = 0; j < L; ++j)
        {
            int s0 = lattice[i][j];
            for (auto [dx, dy] : {std::pair{-1, 0}, std::pair{1, 0}, std::pair{0, -1}, std::pair{0, 1}})
            {
                int x = (i + dx + L) % L;
                int y = (j + dy + L) % L;
                int s1 = lattice[x][y];
                E -= (s1 == s0) ? 1 : 0;
            }
        }
    }
    return J * E * 0.5;
}

double PottsLattice::Energy_Change(int x, int y, int new_spin)
{
    int s0 = lattice[x][y];
    if (s0 == new_spin) return 0.0;

    double E_diff = 0.0;
    // Periodic boundary conditions
    int xm = (x - 1 + L) % L;
    int xp = (x + 1) % L;
    int ym = (y - 1 + L) % L;
    int yp = (y + 1) % L;

    int neighbors[4] = {lattice[xm][y], lattice[xp][y], lattice[x][ym], lattice[x][yp]};

    for (int s_neighbor : neighbors)
    {
        if (new_spin == s_neighbor) E_diff -= J;
        if (s0 == s_neighbor) E_diff += J;
    }
    return E_diff;
}

std::pair<float, float> PottsLattice::Energy_Limit()
{
    float E_min = -2 * L * L * J;
    float E_max = 0;
    return std::make_pair(E_min, E_max);
}

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
}