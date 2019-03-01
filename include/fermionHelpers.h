#ifndef WFERMHELPERS_H
#define WFERMHELPERS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>

using namespace std;

//Wilson Fermion Utilities
//---------------------------------------------------------------------------------

// Zero Wilson fermion field
template<typename T> inline void zeroField(T psi[L][L][2]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	psi[x][y][s] = 0.0;
}

// Copy Wilson fermion field
template<typename T> inline void copyField(T psi2[L][L][2],T psi1[L][L][2]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	psi2[x][y][s] = psi1[x][y][s];
}

// Inner product Wilson fermion field
template<typename T> inline T dotField(T psi1[L][L][2], T psi2[L][L][2]) {
  T scalar = (T) 0.0;
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	scalar += conj(psi1[x][y][s])*psi2[x][y][s];
  return scalar;
}

// Norm squared Wilson fermion field
template<typename T> inline double norm2(T psi[L][L][2]) {
  
  double norm2 = 0.0;
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	norm2 += (psi[x][y][s].real() * psi[x][y][s].real() + psi[x][y][s].imag() * psi[x][y][s].imag());
  
  return norm2/(2*L*L);
}


template<typename T> inline void caxpby(T a, T X[L][L][2],
					T b, T Y[L][L][2],
					T result[L][L][2]){
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	result[x][y][s] = a*X[x][y][s] + b*Y[x][y][s];
}

template<typename T> inline void axpby(double a, T X[L][L][2], double b,
				       T Y[L][L][2], T result[L][L][2]){
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	result[x][y][2] = a*X[x][y][2] + b*Y[x][y][2];
}

//Wilson fermion axpy in place 
template<typename T> inline void axpy(double a, T X[L][L][2], T Y[L][L][2]){ 
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	Y[x][y][s] = a*X[x][y][s] + Y[x][y][s];
}

//Wilson fermion axpy in result
template<typename T> inline void axpy(double a, T X[L][L][2], T Y[L][L][2],
				      T result[L][L][2]){ 
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	result[x][y][s] = a*X[x][y][s] + Y[x][y][s];
}

template<typename T> inline void xpaypbz(T X[L][L][2],
					 double a, T Y[L][L][2],
					 double b, T Z[L][L][2]) {
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++) {
	Z[x][y][s] *= b;
	Z[x][y][s] += X[x][y][s] + a*Y[x][y][s];
      }
}

template<typename T> inline void ax(double a, T X[L][L][2]){
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	X[x][y][s] *= a;
}


// Zero Wilson fermion field
template<typename T> inline void zeroField(T psi[L][L]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      psi[x][y] = 0.0;
}

// Copy Wilson fermion field
template<typename T> inline void copyField(T psi2[L][L],T psi1[L][L]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	psi2[x][y] = psi1[x][y];
}

// Inner product Wilson fermion field
template<typename T> inline T dotField(T psi1[L][L], T psi2[L][L]) {
  T scalar = (T) 0.0;
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      scalar += conj(psi1[x][y])*psi2[x][y];
  return scalar;
}

// Norm squared Wilson fermion field
template<typename T> inline double norm2(T psi[L][L]) {
  
  double norm2 = 0.0;
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      norm2 += (psi[x][y].real() * psi[x][y].real() + psi[x][y].imag() * psi[x][y].imag());
  
  return norm2/(2*L*L);
}


template<typename T> inline void caxpby(T a, T X[L][L],
					T b, T Y[L][L],
					T result[L][L]){
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      result[x][y] = a*X[x][y] + b*Y[x][y];
}

template<typename T> inline void axpby(double a, T X[L][L], double b,
				       T Y[L][L], T result[L][L]){
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      result[x][y] = a*X[x][y] + b*Y[x][y];
}

//Staggered fermion axpy in place 
template<typename T> inline void axpy(double a, T X[L][L], T Y[L][L]){ 
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      Y[x][y] = a*X[x][y] + Y[x][y];
}

//Staggered fermion axpy in result
template<typename T> inline void axpy(double a, T X[L][L], T Y[L][L],
				      T result[L][L]){ 
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      result[x][y] = a*X[x][y] + Y[x][y];
}

template<typename T> inline void xpaypbz(T X[L][L],
					 double a, T Y[L][L],
					 double b, T Z[L][L]) {
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++) {
      Z[x][y] *= b;
      Z[x][y] += X[x][y] + a*Y[x][y];
    }
}

template<typename T> inline void ax(double a, T X[L][L]){
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      X[x][y] *= a;
}

#endif
