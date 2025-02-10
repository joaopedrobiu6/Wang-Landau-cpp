// Wrapper class to simplify initialization of random number generator.

#pragma once

#include <random>

class Rand
{
public:

  std::mt19937_64 gen64;
  std::uniform_real_distribution<double> U;
  std::normal_distribution<double> Normal;

  Rand() {
    init();
  }

  void init() { // Default initialization using random_device:
    std::random_device rd;
    std::seed_seq seed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    gen64.seed(seed);
  }

  void seed(std::seed_seq& seed) {
    gen64.seed(seed);
  }

  inline double operator()() {	// Return real in U[0,1)
    return U(gen64);
  }
  
  inline double normal() {	// Return a sample ~ N(0,1)
    return Normal(gen64);
  }

};

// Initialize it!
static Rand rnd;
