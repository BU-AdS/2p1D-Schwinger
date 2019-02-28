#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>
#include "latHelpers.h"

using namespace std;

typedef struct{
  
  //HMC
  int nstep = 25;
  double tau = 1.0;
  int iterHMC = 1000;
  int therm = 50;
  int skip = 25;
  int chkpt = 100;
  int checkpointStart = 0;
  int maxIterCG = 1000;
  double eps = 1e-6;
  
  //physics
  int Latsize = L;
  double beta = 3.0;
  double betaz = 0.5;
  double m = -0.06;
  bool dynamic = true;
  bool lockedZ = true;
  
  //Smearing
  double alpha = 0.5;
  int smearIter = 1;
  
  //Arpack params
  int nEv = 16;
  int nKv = 32;
  double arpackTol = 1e-6;
  int arpackMaxiter = 10000;
  int polyACC = 0;
  double amax = 10.0;
  double amin = 0.1;
  int n_poly = 100;
  
} param_t;


void printParams(param_t p) {
  cout << endl;
  cout << "Physics:  Size = "<< L << endl;
  cout << "          Beta = "<< p.beta << endl;
  cout << "          Dynamic = " << (p.dynamic == true ? "True" : "False") << endl;
  cout << "          Z locked (3D only) = " << (p.lockedZ == true ? "True" : "False") << endl;
  if (p.dynamic == true) cout << "          Mass = "<< p.m << endl;
  cout << "HMC:      Therm Sweeps: (" << p.therm << " accept) (" << p.therm << " accept/reject)" << endl; 
  cout << "          Data Points = " << p.iterHMC << endl;
  cout << "          Time Step = " << p.tau/p.nstep << endl;
  cout << "          Trajectory Steps " << p.nstep << endl;
  cout << "          Trajectory Length = " << p.tau << endl;
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
  name += "_tau" + to_string(p.tau) + "_nHMCstep" + to_string(p.nstep);
}

void printLattice(Complex gauge[L][L][2]){
  
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++) {
      cout << " (x, y) = " << x << "," << y << " ";
      cout << gauge[x][y][0]<< " " << gauge[x][y][1] << endl;
    }
  return;
}

double measPlaq(Complex gauge[L][L][2]){
  
  double plaq = 0.0;  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++){
      plaq += real(gauge[x][y][0]*gauge[ (x+1)%L ][y][1]*conj(gauge[x][ (y+1)%L ][0])*conj(gauge[x][y][1]));
    }
  return plaq/(L*L);
}


void writeGaugeLattice(Complex gauge[L][L][2], string name){

  fstream outPutFile;
  outPutFile.open(name,ios::in|ios::out|ios::trunc);  
  outPutFile.setf(ios_base::fixed,ios_base::floatfield); 

  //Plaquette action header
  outPutFile << setprecision(20) <<  setw(20) << measPlaq(gauge) << endl;
  
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      for(int mu=0; mu<2; mu++)
	outPutFile << setprecision(12) <<  setw(20) << arg(gauge[x][y][mu]) << endl;
  
  outPutFile.close();
  return;
  
}
void readGaugeLattice(Complex gauge[L][L][2], string name){

  fstream inPutFile;
  inPutFile.open(name);
  string val;
  if(!inPutFile.is_open()) {
    cout << "Error opening file " << name << endl;
    exit(0);
  }
  
  //Header check
  getline(inPutFile, val);
  double plaq = stod(val);

  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      for(int mu=0; mu<D; mu++){
	getline(inPutFile, val);
	gauge[x][y][mu] = polar(1.0, stod(val));	  
      }
  cout << setprecision(16) << setw(20) << "Plaqette on file  = " << plaq << endl;
  cout << setprecision(16) << setw(20) << "Plaqette measured = " << measPlaq(gauge) << endl;
  double err = fabs(1.0 - plaq/measPlaq(gauge));
  if(err > 1e-12) {
    cout << "Gauge read fail!" <<  endl;
    exit(0);
  }    
  return;
}

/*===============================================================================
  Gaussian numbers with p(theta) = sqrt(beta/ 2 PI) exp( - beta* theta^2/2)
  <Gaussian^2> = 1/beta  
  Perimeter Law:  Wilson Loop = exp[ - 4 sigma L ]   sigma = - Log[ <cos(theta)> ]
  ================================================================================*/ 
void gaussStart(Complex gauge[L][L][2],param_t p){

  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++){
      gauge[x][y][0] = polar(1.0,sqrt(1.0/p.beta)*drand48());
      gauge[x][y][1] = polar(1.0,sqrt(1.0/p.beta)*drand48());
    }
  return;
}  

void coldStart(Complex gauge[L][L][2],param_t p){

  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      for(int mu=0;mu<2;mu++)
	gauge[x][y][mu] = Complex(1.0,0.0);
  return;
}  

void gaussReal_F(double field[L][L][2]) {
  //normalized gaussian exp[ - phi*phi/2]  <phi|phi> = 1
  double r, theta, sum;
  double inv_sqrt2 = 1.0/sqrt(2);
  for(int x=0; x<L;x++)
    for(int y=0; y<L; y++){
      r = sqrt(-2.0*log(drand48()));
      theta = TWO_PI*drand48();
      field[x][y][0] = r*cos(theta);
      
      r = sqrt(-2.0*log(drand48()));
      theta = TWO_PI*drand48();
      field[x][y][1] = r*cos(theta);
      //sum += field[x][y][0]*field[x][y][0] + field[x][y][1]*field[x][y][1];
    }
  
  //cout << "mom check: " << sum/(L*L) << endl;
  
  return;
}

void gaussComplex_F(Complex eta[L][L][2], param_t p) {
  
  //normalized gaussian exp[ - eta*eta/2]  <eta|eta> = 1;
  double r1, theta1, r2, theta2;
  double inv_sqrt2 = 1.0/sqrt(2);
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++) {
      for(int s=0; s<2; s++) {
	r1 = sqrt(-2.0*log(drand48()));
	theta1 = TWO_PI*(drand48());
	r1 = sqrt(-2.0*log(drand48()));
	theta1 = TWO_PI*(drand48());
	
	eta[x][y][s] = Complex(r1*cos(theta1),r2*sin(theta2))*inv_sqrt2;
      }
      return;
    }
}

//staple x is 0th, y is 1st.
//APE smearing: project back on U(1)       
void smearLink(Complex Smeared[L][L][2], const Complex gauge[L][L][2], param_t p){

  double alpha = p.alpha;
  int iter = p.smearIter;
  
  Complex SmearedTmp[L][L][2];
  copyLat(Smeared, gauge);
  copyLat(SmearedTmp, Smeared);
  
  for(int i=0; i<iter; i++) {
    
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++){
	SmearedTmp[x][y][0] += alpha * Smeared[x][y][1] * Smeared[x][(y+1)%L][0] * conj(Smeared[x][y][1]);
	SmearedTmp[x][y][0] += alpha * conj(Smeared[x][(y-1 +L)%L][1]) * Smeared[x][(y-1 +L)%L][0] * Smeared[(x+1)%L][(y-1 +L)%L][1];
	SmearedTmp[x][y][1] += alpha * Smeared[x][y][0]* Smeared[(x+1)%L][y][1] * conj(Smeared[(x+1)%L][(y+1)%L][0]);
	SmearedTmp[x][y][1] += alpha * conj(Smeared[(x-1+L)%L][y][0]) * Smeared[(x-1+L)%L][y][1] * Smeared[(x-1+L)%L][(y+1)%L][0];
      }

    //Project back to U(1)
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++)
	for(int mu=0; mu<2; mu++)
	  Smeared[x][y][mu] = polar(1.0,arg(SmearedTmp[x][y][mu]));
  }
}

#endif
