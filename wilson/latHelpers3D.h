#ifndef LATHELPERS3D_H
#define LATHELPERS3D_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>

using namespace std;

//Lattice utilities
//---------------------------------------------------------------------------

// Zero lattice
template<typename T> inline void zeroLat(T v[L][L][LZ][3]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<3; mu++)
	  v[x][y][z][mu] = 0.0;
}

// Copy lattice
template<typename T> inline void copyLat(T v2[L][L][LZ][3], T v1[L][L][LZ][3]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<3; mu++)
	  v2[x][y][z][mu] = v1[x][y][z][mu];
}

// Extract 2D lattice slice
template<typename T> inline void extractLatSlice(T gauge[L][L][LZ][3],
						 T gauge2D[L][L][2],
						 int slice) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int mu=0; mu<2; mu++)
	gauge2D[x][y][mu] = gauge[x][y][slice][mu];
}

// Insert 2D slice
template<typename T> inline void insertLatSlice(T gauge[L][L][LZ][3],
						T gauge2D[L][L][2],
						int slice) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int mu=0; mu<2; mu++)
	gauge[x][y][slice][mu] = gauge2D[x][y][mu];
}

// Zero Wilson Loop field.
template<typename T> inline void zeroWL(T v[L/2][L/2][LZ]) {
  for(int x=0; x<L/2; x++)
    for(int y=0; y<L/2; y++)
      for(int z=0; z<LZ; z++)
	v[x][y][z] = 0.0;  
}

void checkGauge(Complex gauge[L][L][LZ][3], Complex gaugeAlt[L][L][LZ][3], int n) {

#ifdef CHECK_GAUGE
  bool fail = false;
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<3; mu++)
	  if(gauge[x][y][z][mu] != gaugeAlt[x][y][z][mu]) {
	    cout << "Gauge Check fail at " << x << " " << y << " " << z << " " << mu << " from call " << n << endl;
	    fail = true;
	  }
  
  if(fail) {
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++)
	for(int z=0; z<LZ; z++)
	  for(int mu=0; mu<3; mu++) {
	    if(gauge[x][y][z][mu] != gaugeAlt[x][y][z][mu]) {
	      cout << "Gauge compare " << x << " " << y << " " << z << " " << mu << gauge[x][y][z][mu] << " " <<  gaugeAlt[x][y][z][mu] << endl;
	    }
	  }
    exit(0);
  }
#endif
}

#endif
