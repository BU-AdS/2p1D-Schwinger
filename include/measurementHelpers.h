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

void measWilsonLoops(Complex gauge[L][L][2], Complex wLoops[L/2][L/2], param_t p){
  
  Complex w;
  int p1, p2;
  Complex smeared[L][L][2];
  double inv_Lsq = 1.0/(L*L);
  smearLink(smeared, gauge, p);
  
  //Loop over all X side sizes of rectangle 
  for(int Xrect=1; Xrect<L/2; Xrect++) {
      
    //Loop over all Y side sizes of rectangle
    for(int Yrect=1; Yrect<L/2; Yrect++) {

      //Loop over all x,y
      for(int x=0; x<L; x++)
	for(int y=0; y<L; y++){
	    
	  w = Complex(1.0,0.0);
	    
	  //Move in +x up to p1.
	  for(int dx=0; dx<Xrect; dx++)    w *= smeared[ (x+dx)%L ][y][0];
	    
	  //Move in +y up to p2 (p1 constant)
	  p1 = (x + Xrect)%L;
	  for(int dy=0; dy<Yrect; dy++)    w *= smeared[p1][ (y+dy)%L ][1];
	    
	  //Move in -x from p1 to (p2 constant)
	  p2 = (y + Yrect)%L;
	  for(int dx=Xrect-1; dx>=0; dx--)  w *= conj(smeared[ (x+dx)%L ][p2][0]);
	    
	  //Move in -y from p2 to y
	  for(int dy=Yrect-1; dy>=0; dy--)  w *= conj(smeared[x][ (y+dy)%L ][1]);
	  wLoops[Xrect][Yrect] += w*inv_Lsq;
	}
    }
  }  
  return;
}

//Polyakov loops. x is the spatial dim, y is the temporal dim.
void measPolyakovLoops(Complex gauge[L][L][2], Complex pLoops[L/2]){
  
  Complex w1, w2;
  //Eack polyakov loop correlation is defined by its delta x value.
  //We start at x0, separate to x0 + (x0+L/2-1), and loop over all
  //x0=1 -> x0 = L/2-1.

  //Starting x
  for(int x=0; x<L/2; x++) {

    //Loop over time
    w1 = Complex(1.0,0.0);
    for(int dy=0; dy<L; dy++) w1 *= gauge[x][dy][1];
    
    //x separation
    for(int dx=0; dx<L/2; dx++) {
      
      w2 = Complex(1.0,0.0);
      for(int dy=0; dy<L; dy++) w2 *= gauge[x + dx][dy][1];
      
      pLoops[dx] += conj(w1)*w2/(1.0*L/2);
      
    }
  }
    
  return;
}

double measTopCharge(Complex gauge[L][L][2], param_t p){
  
  Complex w;
  double top = 0.0;  
  Complex smeared[L][L][2];
  smearLink(smeared, gauge, p);
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++){
      w = (smeared[x][y][0] * smeared[ (x+1)%L ][y][1] *
	   conj(smeared[x][ (y+1)%L ][0])*conj(smeared[x][y][1]));
      top += arg(w);  // -pi < arg(w) < pi  Geometric value is an integer.
      //print local def here for topology dynamics
      //printf("arg(w) = [ arg(link1) + arg(link2) + c_arg(link3) + c_arg(link4)]
    }
  return top/TWO_PI;
}

double measGaugeAction(Complex gauge[L][L][2], param_t p) {

  double beta  = p.beta;
  double Hgauge = 0.0;
  Complex plaq = 0.0;

  
  for(int x=0; x<L;x++)
    for(int y=0; y<L; y++){      
      plaq = gauge[x][y][0]*gauge[ (x+1)%L ][y][1]*conj(gauge[x][ (y+1)%L ][0])*conj(gauge[x][y][1]);
      Hgauge += beta*real(1.0 - plaq);      
    }  
  return Hgauge;
}

double measMomAction(double mom[L][L][2], param_t p) {

  double Hmom = 0.0;
  Complex plaq;
  
  for(int x=0; x<L;x++)
    for(int y=0; y<L; y++){
      for(int mu=0; mu<2; mu++){
	Hmom += 0.5*mom[x][y][mu] * mom[x][y][mu];
      }
    }
  
  return Hmom;
}

//Staggered fermion
double measFermAction(Complex gauge[L][L][2], Complex phi[L][L],
		      param_t p, bool postStep) {
  
  double Hferm = 0.0;
  Complex phitmp[L][L];
  
  // cout << "Before Fermion force H = " << H << endl;
  Complex scalar = Complex(0.0,0.0);
  zeroField(phitmp);
  Ainv_psi(phitmp, phi, phitmp, gauge, p);
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++){
      if((x+y)%2 == 0)
	scalar += conj(phi[x][y])*phitmp[x][y];
    }
  
  Hferm += real(scalar);
  //cout << "After Fermion Force H  = " << H << endl;
  
  return Hferm;
}


//Staggered Action
double measAction(double mom[L][L][2], Complex gauge[L][L][2],
		  Complex phi[L][L], param_t p, bool postStep) {
  
  double H = 0.0;
  H += measMomAction(mom, p);
  H += measGaugeAction(gauge, p);
  if (p.dynamic) H += measFermAction(gauge, phi, p, postStep);
  
  return H;
}

//Wilson fermion
double measFermAction(Complex gauge[L][L][2], Complex phi[L][L][2],
		      param_t p, bool postStep) {
  
  double Hferm = 0.0;
  Complex phitmp[L][L][2];
  
  //cout << "Before Fermion force H = " << H << endl;
  Complex scalar = Complex(0.0,0.0);
  zeroField(phitmp);
  if(postStep) Ainv_psi(phitmp, phi, phitmp, gauge, p);
  else copyField(phitmp, phi);
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++){
      for(int s=0; s<2; s++){
	scalar += conj(phi[x][y][s])*phitmp[x][y][s];
      }
    }    
  
  Hferm += real(scalar);
  
  return Hferm;
}

//Wilson Action
double measAction(double mom[L][L][2], Complex gauge[L][L][2],
		  Complex phi[L][L][2], param_t p, bool postStep) {
  
  double H = 0.0;
  H += measMomAction(mom, p);
  H += measGaugeAction(gauge, p);
  if (p.dynamic) H += measFermAction(gauge, phi, p, postStep);
  
  return H;
}

//-----------------------------------------------------------------------------------








#endif