################# -*-makefile-*-

all:	ising2d # ising3d

ising2d:	ising2d.cpp Rand.h
	$(CCC) $(CCFLAGS) $< -o $@ $(LIBS)

ising3d:	ising3d.cpp Rand.h
	$(CCC) $(CCFLAGS) $< -o $@ $(LIBS)

tar:;	tar -cvf lab03_ising.tar ising2d.cpp Rand.h Makefile params analyze.sh analyze_adv.m exact_solution.data Ising.py

CCC=c++
CCFLAGS=-O3 -march=native -Wall -std=c++17
LIBS=
