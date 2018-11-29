#ifndef LINALGHELPERS_H
#define LINALGHELPERS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>
#include <omp.h>

using namespace std;

// Zero lattice vectors.
template<typename T> inline void zeroLat(T v[L][L][LZ][D]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<D; mu++)
	  v[x][y][z][mu] = 0.0;
}

// Zero lattice vectors.
template<typename T> inline void zeroLat2D(T v[L][L][D]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int mu=0; mu<D; mu++)
	v[x][y][mu] = 0.0;
}


// Copy lattice vector
template<typename T> inline void copyLat(T v2[L][L][LZ][D],T v1[L][L][LZ][D]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0;mu<D;mu++)
	  v2[x][y][z][mu] =  v1[x][y][z][mu];
}

// Copy lattice vector
template<typename T> inline void copyLat2D(T v2[L][L][2],T v1[L][L][2]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int mu=0; mu<2; mu++)
	v2[x][y][mu] = v1[x][y][mu];
}

// Extract 2D slice
template<typename T> inline void extractLatSlice(T gauge[L][L][LZ][D],
						 T gauge2D[L][L][D],
						 int slice) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int mu=0;mu<D;mu++)
	gauge2D[x][y][mu] =  gauge[x][y][slice][mu];
}


// Insert 2D slice
template<typename T> inline void insertLatSlice(T gauge[L][L][LZ][D],
						T gauge2D[L][L][D],
						int slice) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int mu=0;mu<D;mu++)
	gauge[x][y][slice][mu] = gauge2D[x][y][mu];
}


// Copy lattice field
template<typename T> inline void copyField(T psi2[L][L][LZ], T psi1[L][L][LZ]) {

  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++)
	psi2[x][y][z] = psi1[x][y][z];
}

// Copy lattice field
template<typename T> inline void copyField2D(T psi2[L][L],T psi1[L][L]) {
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      psi2[x][y] =  psi1[x][y];
}


// Zero lattice field.
template<typename T> inline void zeroField(T psi[L][L][LZ]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++) {
	psi[x][y][z] = 0.0;
      }
}

// Zero lattice field.
template<typename T> inline void zeroField2D(T psi[L][L]) {

  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      psi[x][y] = 0.0;
}


// Add Equ conj(v2) dot v1 lattice field
template<typename T> inline T dotField2D(T psi1[L][L], T psi2[L][L]) {
  T scalar = (T) 0.0;
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      scalar += conj(psi1[x][y])*psi2[x][y];
  return scalar;
}

double norm2(Complex psi[L][L]) {

  double norm2 = 0.0;
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      norm2 += psi[x][y].real() * psi[x][y].real() + psi[x][y].imag() * psi[x][y].imag();

  return norm2/(L*L);
}

Complex cDotProduct(Complex psi1[L][L], Complex psi2[L][L]) {

  Complex prod(0.0,0.0);

  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      prod += conj(psi1[x][y])*psi2[x][y];

  return prod;
}


void caxpby(Complex a, Complex X[L][L], Complex b,
            Complex Y[L][L], Complex result[L][L]){

  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      result[x][y] = a*X[x][y] + b*Y[x][y];
}

void printParams(param_t p) {
  cout << endl;
  cout << "Physics:  Size = "<< L << endl;
  cout << "          Beta = "<< p.beta << endl;
  cout << "          BetaZ = "<< p.betaz << endl;
  cout << "          Quenching = " << (p.quenched == true ? "True" : "False") << endl;
  if (!p.quenched) cout << "          Mass = "<< p.m << endl;
  cout << "          Z locked = " << (p.lockedZ == true ? "True" : "False") << endl;
  cout << "HMC:      Therm Sweeps = " << p.therm << endl; 
  cout << "          Data Points = " << p.iterHMC - p.therm << endl;
  cout << "          Time Step = " << p.dt << endl;
  cout << "          Trajectory Steps " << p.nstep << endl;
  cout << "          Trajectory Length = " << p.dt*p.nstep << endl;
  cout << "          Trajectory Length = " << p.dt*p.nstep << endl;
  cout << "Smearing: APE iter = " << p.smearIter << endl;
  cout << "          APE alpha = " << p.alpha << endl;
  #ifdef USE_ARPACK
  cout << "ARPACK:   nkv = " << p.nKv << endl;
  cout << "          nev = " << p.nEv << endl;
  cout << "          tol = " << p.arpackTol << endl;
  cout << "          maxiter = " << p.arpackMaxiter << endl;
#endif  
}

void constructName(string &name, param_t p) {
  name += "_L" + to_string(p.Latsize) + "_LZ" + to_string(LZ) + "_B" + to_string(p.beta);
  if(p.quenched == false) name += "_M"+ to_string(p.m);
  name += "_dt" + to_string(p.dt) + "_nHMCstep" + to_string(p.nstep);
}

void checkGauge(Complex gauge[L][L][LZ][D], Complex gaugeAlt[L][L][LZ][D], int n) {

#ifdef CHECK_GAUGE
  bool fail = false;
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<D; mu++)
	  if(gauge[x][y][z][mu] != gaugeAlt[x][y][z][mu]) {
	    cout << "Gauge Check fail at " << x << " " << y << " " << z << " " << mu << " from call " << n << endl;
	    fail = true;
	  }
  
  if(fail) {
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++)
	for(int z=0; z<LZ; z++)
	  for(int mu=0; mu<D; mu++) {
	    if(gauge[x][y][z][mu] != gaugeAlt[x][y][z][mu]) {
	      cout << "Gauge compare " << x << " " << y << " " << z << " " << mu << gauge[x][y][z][mu] << " " <<  gaugeAlt[x][y][z][mu] << endl;
	    }
	  }
    exit(0);
  }
#endif
}

#endif
