#ifndef MEASUREMENTHELPERS_H
#define MEASUREMENTHELPERS_H

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
// 2 Dimensional routines 
//-----------------------------------------------------------------------------------

//======================================================================
//   Creutz     exp[ -sigma L^2] exp[ -sigma(L-1)(L-1)]
//   ratio:    ---------------------------------------  = exp[ -sigma]
//              exp[ -sigma (L-1)L] exp[-sigma L(L-1)]
//======================================================================

void measWilsonLoops(Complex gauge[LX][LY][2], Complex wLoops[LX/2][LY/2], param_t p){
  
  Complex w;
  int p1, p2;
  Complex smeared[LX][LY][2];
  double inv_Lsq = 1.0/(LX*LY);
  smearLink(smeared, gauge, p);
  
  //Loop over all X side sizes of rectangle 
  for(int Xrect=1; Xrect<LX/2; Xrect++) {
      
    //Loop over all Y side sizes of rectangle
    for(int Yrect=1; Yrect<LY/2; Yrect++) {

      //Loop over all x,y
      for(int x=0; x<LX; x++)
	for(int y=0; y<LY; y++){
	    
	  w = Complex(1.0,0.0);
	    
	  //Move in +x up to p1.
	  for(int dx=0; dx<Xrect; dx++)    w *= smeared[ (x+dx)%LX ][y][0];
	    
	  //Move in +y up to p2 (p1 constant)
	  p1 = (x + Xrect)%LX;
	  for(int dy=0; dy<Yrect; dy++)    w *= smeared[p1][ (y+dy)%LY ][1];
	    
	  //Move in -x from p1 to (p2 constant)
	  p2 = (y + Yrect)%LY;
	  for(int dx=Xrect-1; dx>=0; dx--)  w *= conj(smeared[ (x+dx)%LX ][p2][0]);
	    
	  //Move in -y from p2 to y
	  for(int dy=Yrect-1; dy>=0; dy--)  w *= conj(smeared[x][ (y+dy)%LY ][1]);
	  wLoops[Xrect][Yrect] += w*inv_Lsq;
	}
    }
  }  
  return;
}

//Polyakov loops. x is the spatial dim, y is the temporal dim.
void measPolyakovLoops(Complex gauge[LX][LY][2], Complex pLoops[LX/2]){
  
  Complex w1, w2;
  //Eack polyakov loop correlation is defined by its delta x value.
  //We start at x0, separate to x0 + (x0+L/2-1), and loop over all
  //x0=1 -> x0 = L/2-1.

  //Starting x
  for(int x=0; x<LX/2; x++) {

    //Loop over time
    w1 = Complex(1.0,0.0);
    for(int dy=0; dy<LY; dy++) w1 *= gauge[x][dy][1];
    
    //x separation
    for(int dx=0; dx<LX/2; dx++) {
      
      w2 = Complex(1.0,0.0);
      for(int dy=0; dy<LY; dy++) w2 *= gauge[x+dx][dy][1];
      
      pLoops[dx] += conj(w1)*w2/(1.0*LX/2);
      
    }
  }
    
  return;
}

double measTopCharge(Complex gauge[LX][LY][2], param_t p){
  
  Complex w;
  double top = 0.0;  
  Complex smeared[LX][LY][2];
  smearLink(smeared, gauge, p);
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++){
      w = (smeared[x][y][0] * smeared[ (x+1)%LX ][y][1] *
	   conj(smeared[x][ (y+1)%LY ][0])*conj(smeared[x][y][1]));
      top += arg(w);  // -pi < arg(w) < pi  Geometric value is an integer.
      //print local def here for topology dynamics
      //printf("arg(w) = [ arg(link1) + arg(link2) + c_arg(link3) + c_arg(link4)]
    }
  return top/TWO_PI;
}

double measGaugeAction(Complex gauge[LX][LY][2], param_t p) {

  double beta  = p.beta;
  double Hgauge = 0.0;
  Complex plaq = 0.0;

  
  for(int x=0; x<LX;x++)
    for(int y=0; y<LY; y++){      
      plaq = gauge[x][y][0]*gauge[ (x+1)%LX ][y][1]*conj(gauge[x][ (y+1)%LY ][0])*conj(gauge[x][y][1]);
      Hgauge += beta*real(1.0 - plaq);      
    }  
  return Hgauge;
}

double measMomAction(double mom[LX][LY][2], param_t p) {

  double Hmom = 0.0;
  Complex plaq;
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++){
      for(int mu=0; mu<2; mu++){
	Hmom += 0.5*mom[x][y][mu] * mom[x][y][mu];
      }
    }
  
  return Hmom;
}

//Staggered fermion
double measFermAction(Complex gauge[LX][LY][2], Complex phi[LX][LY],
		      param_t p, bool postStep) {
  
  double Hferm = 0.0;
  Complex phitmp[LX][LY];
  
  // cout << "Before Fermion force H = " << H << endl;
  Complex scalar = Complex(0.0,0.0);
  zeroField(phitmp);
  Ainvpsi(phitmp, phi, phitmp, gauge, p);
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++){
      if((x+y)%2 == 0)
	scalar += conj(phi[x][y])*phitmp[x][y];
    }
  
  Hferm += real(scalar);
  //cout << "After Fermion Force H  = " << H << endl;
  
  return Hferm;
}


//Staggered Action
double measAction(double mom[LX][LY][2], Complex gauge[LX][LY][2],
		  Complex phi[LX][LY], param_t p, bool postStep) {
  
  double H = 0.0;
  H += measMomAction(mom, p);
  H += measGaugeAction(gauge, p);
  if (p.dynamic) H += measFermAction(gauge, phi, p, postStep);
  
  return H;
}

//Wilson fermion
double measFermAction(Complex gauge[LX][LY][2], Complex phi[LX][LY][2],
		      param_t p, bool postStep) {
  
  double Hferm = 0.0;
  Complex phitmp[LX][LY][2];
  
  //cout << "Before Fermion force H = " << H << endl;
  Complex scalar = Complex(0.0,0.0);
  zeroField(phitmp);
  if(postStep) Ainvpsi(phitmp, phi, phitmp, gauge, p);
  else copyField(phitmp, phi);
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++){
      for(int s=0; s<2; s++){
	scalar += conj(phi[x][y][s])*phitmp[x][y][s];
      }
    }    
  
  Hferm += real(scalar);
  
  return Hferm;
}

//Wilson Action
double measAction(double mom[LX][LY][2], Complex gauge[LX][LY][2],
		  Complex phi[LX][LY][2], param_t p, bool postStep) {
  
  double H = 0.0;
  H += measMomAction(mom, p);
  H += measGaugeAction(gauge, p);
  if (p.dynamic) H += measFermAction(gauge, phi, p, postStep);
  
  return H;
}

//-----------------------------------------------------------------------------------








#endif
