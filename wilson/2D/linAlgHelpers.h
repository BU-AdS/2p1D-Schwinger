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
      for(int mu=0; mu<D; mu++)
	v2[x][y][mu] +=  v1[x][y][mu];
}

// Zero lattice field.
template<typename T> inline void zeroField(T psi[L][L][2]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	psi[x][y][s] = 0.0;
}

// Zero half lattice field.
template<typename T> inline void zeroHalfField(T psi[L/2][L/2]) {
  for(int x=0; x<L/2; x++)
    for(int y=0; y<L/2; y++)
      psi[x][y] = 0.0;
}

// Copy lattice field
template<typename T> inline void copyField(T psi2[L][L][2],T psi1[L][L][2])
{
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	psi2[x][y][s] = psi1[x][y][s];
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
template<typename T> inline T dotField(T psi1[L][L][2], T psi2[L][L][2]) {
  T scalar = (T) 0.0;
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	scalar += conj(psi1[x][y][s])*psi2[x][y][s];
  return scalar;
}

double norm2(Complex psi[L][L][2]) {

  double norm2 = 0.0;
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	norm2 += (psi[x][y][s].real() * psi[x][y][s].real() + psi[x][y][s].imag() * psi[x][y][s].imag());
  
  return norm2/(2*L*L);
}

void caxpby(Complex a, Complex X[L][L][2],
	    Complex b, Complex Y[L][L][2], Complex result[L][L][2]){
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	result[x][y][s] = a*X[x][y][s] + b*Y[x][y][s];
}

void axpby(double a, Complex X[L][L], double b,
	   Complex Y[L][L], Complex result[L][L]){
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      result[x][y] = a*X[x][y] + b*Y[x][y];
}

void axpy(double a, Complex X[L][L][2], Complex Y[L][L][2]){ 
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	Y[x][y][s] = a*X[x][y][s] + Y[x][y][s];
}


void axpy(double a, Complex X[L][L][2], Complex Y[L][L][2], Complex result[L][L][2]){ 
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	result[x][y][s] = a*X[x][y][s] + Y[x][y][s];
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

void ax(double a, Complex X[L][L][2]){
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	X[x][y][s] *= a;
}


void printParams(param_t p) {
  cout << endl;
  cout << "Physics:  Size = "<< L << endl;
  cout << "          Beta = "<< p.beta << endl;
  cout << "          Dynamic = " << (p.dynamic == true ? "True" : "False") << endl;
  if (p.dynamic == true) cout << "          Mass = "<< p.m << endl;
  cout << "HMC:      Therm Sweeps: (" << p.therm << " accept) (" << p.therm << " accept/reject)" << endl; 
  cout << "          Data Points = " << p.iterHMC << endl;
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
