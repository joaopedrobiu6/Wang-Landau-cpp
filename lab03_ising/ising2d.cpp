// ising2d.cc - 2D Ising model

// Written by Jack Lidmar, 2024.
// This file is made available under the Creative Commons CC0 1.0 Universal Public Domain Dedication.

#include <iostream>
#include <fstream>
#include <cmath>
#include <array>

#include "Rand.h"		// Random number generator

using namespace std;

class Lattice2D
{
public:
  int Lx, Ly, N;
  int8_t* s;  //  Stores the spins +/-1

  Lattice2D() { s = 0; }
  ~Lattice2D() { if (s) delete [] s; }

  void init(const int Lx_, const int Ly_) {
    Lx = Lx_;
    Ly = Ly_;
    N  = Lx*Ly;
    s = new int8_t[N];
  }

  int8_t& operator()(int x, int y) {
    return s[x + Lx*y];
  }

  void set(const int8_t s_) {
    for (int i = 0; i < N; i++)
      s[i] = s_;
  }
};

// The following is an alternative data structure to hold the spins. You can try out which is faster.
class Lattice2Dv
{
public:
  int Lx, Ly, N;
  vector<vector<int8_t>> s; //  Stores the spins +/-1

  void init(const int Lx_, const int Ly_) {
    Lx = Lx_;
    Ly = Ly_;
    N  = Lx*Ly;
    s.resize(Lx);
    for (auto& v : s)
      v.resize(Ly);
  }

  int8_t& operator()(int x, int y) {
    return s[x][y];
  }

  void set(int8_t s_) {
    for (auto& v : s)
      v.assign(Ly, s_);
  }
};

class Lattice3D // (change)
{
  // ... to be filled in ...
};

class Simulation
{
public:
  int L;			  // Size
  double T;			// Temperature.

  double Energy;	// Current energy.
  int M;			    // Current magnetization.

  Lattice2D s;		// Spins.

  int time;			  // Number of sweeps.
  int sweeps;			// Number of sweeps per sample.

  int nsamp;			// Number of samples (= time/sweeps).

  array<double,13> Prob; // Table of exp(-dE/T) to speed up calculation.

  // Averages:
  double e, e2;
  double m, m2, m4, mabs;
  double chi, cv, binder;

  double accratio;		// Acceptance ratio.

  vector<int> plus, minus; // Used to implement periodic boundaries.

  void init(int L_) {
    L = L_;
    s.init(L,L);  // (change)
    init_pbc(L);

    s.set(1);			// All up!
    M = L*L;			// Total magnetization (change)

    Energy = Energy_calc();
  }

  void init_pbc(int L) { // Set up tables for periodic boundary conditions
    plus.resize(L);
    minus.resize(L);
    for (int x = 0; x < L-1; x++) {
      plus[x] = x+1;
      minus[x+1] = x;
    }
    plus[L-1] = 0;
    minus[0] = L-1;
  }

  // Other options to implement periodic boundaries:
  // int plus(int x)  { return (x + 1) % L; } // ... slow due to modulo (%)
  // int minus(int x) { return (x + L - 1) % L; }

  void set_T(double T_) {
    T = T_;
    // Initialize precalculated table of acceptance probabilities:
    for (int k = 0; k < int(Prob.size()); k++)
      Prob[k] = exp(-k/T);
  }

  void mcmc_sweep(int nsweep = 1) {   // (change)
    int acc = 0;
    nsweep *= L*L;

    for (int ii = 0; ii < nsweep; ii++) { // Go through all sites.

      // Choose a spin randomly:
      int x = int(L*rnd());
      int y = int(L*rnd());

      int8_t sxy = s(x,y);

      // There are many ways to implement periodic boundary condtions.
      // Here we use the tables declared in plus and minus:
      int xp = plus[x], yp = plus[y], xm = minus[x], ym = minus[y];

      //   Calculate energy difference for new configuration   //

      int8_t dE = 2*sxy*(s(xp,y) + s(x,yp) + s(xm,y) + s(x,ym));

      // if (dE < 0 or exp(-dE/T) > rnd()) { // Metropolis

      // Avoid calculation of exponential by using a table:
      // (this works when the energy change is an integer)
      if (dE < 0 or Prob[dE] > rnd()) { // Metropolis
        s(x,y) = -sxy;			   // Flip the spin
        Energy += dE;
        M -= 2*sxy;			   // Keep track of magnetization.
        acc++;
      }
    }

    accratio += double(acc)/nsweep;
    time++;

  } // mcmc_sweep

  // (change)
  double Energy_calc() { // Calculate energy for current configuration
    double E = 0;
    for (int x = 0; x < L; x++) {
      for (int y = 0; y < L; y++) {
        E += s(x,y)*(s(plus[x],y) + s(x,plus[y]));
      }
    }
    E = -E;
    return E;
  }

  double Energy_check() {
    double E = Energy_calc();
    if (fabs(Energy-E) > 1e-8)
      cerr << "Energy missmatch: "
           << Energy << ' ' << E << ' ' << Energy - E << endl;
    return E;
  }

  // Equilibrate the system while writing the energy vs time to a file "equil"
  virtual void equil(int t) {
    time = 0;
    accratio = 0;
    mcmc_sweep(2);
    static ofstream out("equil");
    for (int ti = 2; ti < t; ti += ti) {
      double E = 0;
      for (int i = 0; i < ti; i += 1) {
        mcmc_sweep(1);
        E += Energy;
      }
      out << time << '\t' << E/ti << endl;
    }
    out << '&' << endl;
  }

  void progress_bar(int i, int total) {
    if ((i + 1) % 10000 == 0) {
      cerr << "\r  L = " << L << "  T = " << T
           << "    step: " << i+1 << " (" << total << ")"
           << "                              \r";
    }
  }

  void run(long t, int sweeps_) {
    reset();
    time = 0;

    accratio = 0;
    sweeps = sweeps_;

    for (long i = 0; i < t; i++) {
      progress_bar(i,t);
      mcmc_sweep(sweeps);
      samp();
    }

    results();

    Energy = Energy_check();

  }

  void reset() {
    nsamp = 0;

    e = 0;
    e2 = 0;
    m = 0;
    m2 = 0;
    m4 = 0;
    mabs = 0;
  }

  inline void samp() {
    nsamp++;

    e  += Energy;
    e2 += Energy*Energy;

    m  += M;
    mabs += abs(M);		// Average <|M|> rather than <M>.
    double M2 = M*M;  // Avoid integer overflow...
    m2 += M2;
    m4 += M2*M2;
  }

  void results() {
    e  /= nsamp;
    e2 /= nsamp;
    m  /= nsamp;
    m2 /= nsamp;
    m4 /= nsamp;
    mabs /= nsamp;
    accratio /= nsamp;

    chi = (m2 - mabs*mabs)/T;
    cv  = (e2 - e*e)/(T*T);
    binder = 1.0 - m4/(3*m2*m2);

    e /= s.N; // Divide by system size
    e2 /= s.N;
    m  /= s.N;
    m2 /= s.N;
    m4 /= s.N;
    mabs /= s.N;
    chi /= s.N;
    cv  /= s.N;
  }

  void output(ostream & out = cout) {
    constexpr char tab = '\t';
    // Columns: T, m, abs(m), chi, e, cv, binder, L, nsamp, sweeps, accratio
    out << T << tab         // 1
        << m << tab         // 2
        << mabs << tab      // 3
        << chi << tab       // 4
        << e << tab         // 5
        << cv << tab        // 6
        << binder << tab    // 7
        << L << tab         // 8
        << nsamp << tab     // 9
        << sweeps << tab    // 10
        << accratio         // 11
        << endl << flush;
  }
};

int main(int argc, char *argv[]) {
 
  ifstream fin;
  if (argc == 1) {
    std::cerr << "Enter simulation parameters\n"
    "L Tstart Tstep Tstop Nsamp Samples/sweep Nequil" << endl << "> ";
  }
  else if (argc == 2) {
    fin.open(argv[1]);
    if (!fin) { cerr << "Error reading from file " << argv[1] << endl; return(1); }
  }
  else {
    std::cerr << "Usage: " << argv[0] << " [parameter_file]" << endl;
    return(2);
  }

  // Read parameters either from an input file or from stdin:
  istream& in = (argc == 2 && fin.good()) ? fin : std::cin;
  cout.precision(10);

  Simulation s;

  while ( in >> ws ) {

    if (in.peek() == '#') { // Skip lines beginning with #
      in.ignore(1000,'\n');
      continue;
    }

    int L;
    double Tstart, Tstep, Tstop;
    int number, sweeps, equ;

    in >> L >> Tstart >> Tstep >> Tstop >> number >> sweeps >> equ;

    if (!in) break;

    ofstream datafile("data" + to_string(L), ios::trunc);  // Overwrite the file if already present.
    datafile.precision(10);

    s.init(L);

    for (double T = Tstart; T*Tstep < (Tstop+1e-6)*Tstep; T += Tstep) {
      s.set_T(T);
      s.equil(equ);
      s.run(number,sweeps);
      s.output(datafile);
    }

    cerr << endl;
  }

  return 0;
}
