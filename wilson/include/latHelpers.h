#ifndef LATHELPERS_H
#define LATHELPERS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>

using namespace std;

//Lattice utilities
//---------------------------------------------------------------------------

// Zero lattice 2D
template<typename T> inline void zeroLat(T v[L][L][2]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int mu=0; mu<2; mu++)
	v[x][y][mu] = 0.0;
}

// Copy lattice 2D
template<typename T> inline void copyLat(T v2[L][L][2], const T v1[L][L][2]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int mu=0; mu<2; mu++)
	v2[x][y][mu] = v1[x][y][mu];
}

// Zero Wilson Loop field 2D.
template<typename T> inline void zeroWL(T psi[L/2][L/2]) {
  for(int x=0; x<L/2; x++)
    for(int y=0; y<L/2; y++)
      psi[x][y] = 0.0;
}

#endif
