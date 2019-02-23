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
template<typename T> inline void zeroLat(T v[L][L][D]) {
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      for(int mu=0;mu<D;mu++)
	v[x][y][mu] = 0.0;
}

// Copy lattice vector
template<typename T> inline void copyLat(T v2[L][L][D],T v1[L][L][D]) {
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      for(int mu=0;mu<D;mu++)
	v2[x][y][mu] =  v1[x][y][mu];
}

// Add Equ v2 += v1 lattice vector
template<typename T> inline void addEqLat(T v2[L][L][D],T v1[L][L][D]) {
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      for(int mu=0;mu<D;mu++)
	v2[x][y][mu] +=  v1[x][y][mu];
}

// Zero lattice field.
template<typename T> inline void zeroField(T psi[L][L]) {
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      psi[x][y] = 0.0;
}

// Zero half lattice field.
template<typename T> inline void zeroHalfField(T psi[L/2][L/2]) {
  for(int x=0; x<L/2; x++)
    for(int y=0; y<L/2; y++)
      psi[x][y] = 0.0;
}

// Copy lattice field
template<typename T> inline void copyField(T psi2[L][L],T psi1[L][L])
{
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      psi2[x][y] =  psi1[x][y];
}

// Add Equ v2 += v1 lattice field
template<typename T> inline void addEqField(T psi2[L][L],T psi1[L][L]) {
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      psi2[x][y] +=  psi1[x][y];
}

// Add Equ square of real b dot b lattice field
template<typename T> inline T sqrSumField(T b[L][L]) {
  T square = (T) 0.0;
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      square +=  b[x][y]*b[x][y];
  return square;
}

// Add Equ conj(v2) dot v1 lattice field
template<typename T> inline T dotField(T psi1[L][L], T psi2[L][L]) {
  T scalar = (T) 0.0;
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      scalar +=  conj(psi1[x][y])*psi2[x][y];
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

void axpby(double a, Complex X[L][L], double b,
	   Complex Y[L][L], Complex result[L][L]){
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      result[x][y] = a*X[x][y] + b*Y[x][y];
}

void xpaypbz(Complex X[L][L],
	     double a, Complex Y[L][L],
	     double b, Complex Z[L][L]) {
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++) {
      //result[x][y] = X[x][y] + a*Y[x][y] + b*Z[x][y];
      Z[x][y] *= b;
      Z[x][y] += X[x][y] + a*Y[x][y];
    }
}

void ax(double a, Complex X[L][L]){
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      X[x][y] *= a;
}


void printParams(param_t p) {
  cout << endl;
  cout << "Physics:  Size = "<< L << endl;
  cout << "          Beta = "<< p.beta << endl;
  cout << "          Dynamic = " << (p.dynamic == true ? "True" : "False") << endl;
  if (p.dynamic) cout << "          Mass = "<< p.m << endl;
  cout << "HMC:      Therm Sweeps = " << p.therm << endl; 
  cout << "          Data Points = " << p.iterHMC - p.therm << endl;
  cout << "          Time Step = " << p.dt << endl;
  cout << "          Trajectory Steps " << p.nstep << endl;
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
  name += "_L" + to_string(p.Latsize) + "_B" + to_string(p.beta);
  if(p.dynamic == true) name += "_M"+ to_string(p.m);
  name += "_dt" + to_string(p.dt) + "_nHMCstep" + to_string(p.nstep);
}

#endif
