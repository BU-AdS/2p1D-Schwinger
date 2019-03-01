#ifndef UTILS3D_H
#define UTILS3D_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>

using namespace std;

double measPlaq(const Complex gauge[LX][LY][LZ][3], int z){
  
  double plaq = 0.0;
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      plaq += real(gauge[x][y][z][0]*gauge[ (x+1)%LX ][y][z][1]*conj(gauge[x][ (y+1)%LY ][z][0])*conj(gauge[x][y][z][1]));
  
  return plaq/(LX*LY);
}


void writeGaugeLattice(const Complex gauge[LX][LY][LZ][3], string name){

  fstream outPutFile;
  outPutFile.open(name,ios::in|ios::out|ios::trunc);  
  outPutFile.setf(ios_base::fixed,ios_base::floatfield); 

  //Plaquette action header
  for(int z=0; z<LZ; z++) outPutFile << setprecision(20) <<  setw(20) << measPlaq(gauge, z) << endl;
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<3; mu++)
	  outPutFile << setprecision(20) <<  setw(20) <<  arg(gauge[x][y][z][mu]) << endl;
  
  outPutFile.close();
  return;
}

void readGaugeLattice(Complex gauge[LX][LY][LZ][3], string name){

  fstream inPutFile;
  inPutFile.open(name);
  string val;
  
  //Header check
  double plaq[LZ];
  for(int z=0; z<LZ; z++) {
    getline(inPutFile, val);
    plaq[z] = stod(val);
  }
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<3; mu++) {
	  getline(inPutFile, val);
	  gauge[x][y][z][mu] = polar(1.0, stod(val));	  
	}

  for(int z=0; z<LZ; z++) {
    cout << setprecision(16) << setw(20) << "Plaqette["<<z<<"] on file  = " << plaq[z] << endl;
    cout << setprecision(16) << setw(20) << "Plaqette["<<z<<"] measured = " << measPlaq(gauge, z) << endl;
    double err = fabs(1.0 - plaq[z]/measPlaq(gauge, z));
    if(err > 1e-12) {
      cout << "Gauge read fail at slice " << z << endl;
      exit(0);
    }    
  }
  inPutFile.close();
  return;
}

void printLattice(Complex gauge[LX][LY][LZ][3]){
  
  for(int x=0; x<LX; x++)    
    for(int y=0; y<LY; y++)
      for(int z=0; z<LZ; z++) {
	cout << "["<< x << "," << y << "," << z << "] = ";
	cout << gauge[x][y][z][0] << " "
	     << gauge[x][y][z][1] << " "
	     << gauge[x][y][z][2] << endl;
      }
  return;
}

// Gaussian numbers with p(theta) = sqrt(beta/ 2 PI) exp( - beta* theta^2/2)
// <Gaussian^2> = 1/beta  
// Perimeter Law:  Wilson Loop = exp[ - 4 sigma L ]   sigma = - Log[ <cos(theta)> ]
//=================================================================================

void gaussStart(Complex gauge[LX][LY][LZ][3], param_t p){  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<3; mu++) {
	  gauge[x][y][z][mu] = polar(1.0, TWO_PI*drand48());
	  if(p.lockedZ && mu == 2) gauge[x][y][z][mu] = 1.0;
	}
  return;
}  

void coldStart(Complex gauge[LX][LY][LZ][3], param_t p){

  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int z=0; z<LZ; z++) 
	for(int mu=0; mu<3; mu++)
	  gauge[x][y][z][mu] = Complex(1.0,0.0);
  return;
}  


void gaussReal_F(double field[LX][LY][LZ][3]) {
  //normalized gaussian exp[ - phi*phi/2]  <phi|phi> = 1
  double r, theta, sum;
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int z=0; z<LZ; z++){
	r = sqrt(-2.0*log(drand48()));
	theta = TWO_PI*drand48();
	field[x][y][z][0] = r*cos(theta);
	
	r = sqrt(-2.0*log(drand48()));
	theta = TWO_PI*drand48();
	field[x][y][z][1] = r*cos(theta);
	//sum += field[x][y][0]*field[x][y][0] + field[x][y][1]*field[x][y][1];
      }
  
  //cout << "mom check: " << sum/(L*L) << endl;
  
  return;
}

#endif
