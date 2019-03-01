#ifndef MEASUREMENTHELPERS3D_H
#define MEASUREMENTHELPERS3D_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>
#include "utils.h"
#include "fermionHelpers.h"
#include "dOpHelpers.h"

using namespace std;

//-----------------------------------------------------------------------------------
// 3 Dimensional routines 
//-----------------------------------------------------------------------------------

double measGaugeAction(Complex gauge[L][L][LZ][3], param_t p) {
  
  double beta = p.beta;
  double betaz = p.betaz;
  double Hgauge = 0.0;
  Complex plaq = 0.0;;  
  
  for(int x=0; x<L;x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++) {

	plaq = gauge[x][y][z][0]*gauge[ (x+1)%L ][y][z][1]*conj(gauge[x][ (y+1)%L ][z][0])*conj(gauge[x][y][z][1]);
	Hgauge += beta*real(1.0 - plaq);
	
	//Compute extra dim contribution
	if(z != LZ-1) {
	  //+x, +z, -x, -z
	  plaq = gauge[x][y][z][0] * cUnit * conj(gauge[x][y][ (z+1)%LZ ][0]) * cUnit;
	  Hgauge += betaz*real(1.0 - plaq);
	  
	  //+y, +z, -y, -z
	  plaq = gauge[x][y][z][1] * cUnit * conj(gauge[x][y][ (z+1)%LZ ][1]) * cUnit;
	  Hgauge += betaz*real(1.0 - plaq);	
	}
      }
  
  return Hgauge;
}

double measMomAction(double mom[L][L][LZ][3], param_t p) {
  
  double Hmom = 0.0;
  Complex plaq;
  
  for(int x=0; x<L;x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<3; mu++){
	  Hmom += 0.5*mom[x][y][z][mu] * mom[x][y][z][mu];
	}
  
  
  return Hmom;
}

//Staggered
double measAction(double mom[L][L][LZ][3], Complex gauge[L][L][LZ][3],
		  Complex phi[L][L], param_t p, bool postStep) {
  
  double H = 0.0;
  H += measMomAction(mom, p);
  H += measGaugeAction(gauge, p);
  if (p.dynamic) {
    Complex gauge2D[L][L][2];
    extractLatSlice(gauge, gauge2D, (LZ-1)/2);
    H += measFermAction(gauge2D, phi, p, postStep);
  }
  
  return H;
}

//Wilson
double measAction(double mom[L][L][LZ][3], Complex gauge[L][L][LZ][3],
		  Complex phi[L][L][2], param_t p, bool postStep) {
  
  double H = 0.0;
  H += measMomAction(mom, p);
  H += measGaugeAction(gauge, p);
  if (p.dynamic) {
    Complex gauge2D[L][L][2];
    extractLatSlice(gauge, gauge2D, (LZ-1)/2);
    H += measFermAction(gauge2D, phi, p, postStep);
  }
  return H;
}

#endif
//-----------------------------------------------------------------------------------
