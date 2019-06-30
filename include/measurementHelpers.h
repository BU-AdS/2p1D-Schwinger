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
#include "inverters.h"

#ifdef USE_ARPACK
#include "arpack_interface_wilson.h"
#endif


using namespace std;

//-----------------------------------------------------------------------------------
// 2 Dimensional routines 
//-----------------------------------------------------------------------------------

//   Creutz     exp[ -sigma L^2] exp[ -sigma(L-1)(L-1)]
//   ratio:    ---------------------------------------  = exp[ -sigma]
//              exp[ -sigma (L-1)L] exp[-sigma L(L-1)]
void measWilsonLoops(Complex gauge[LX][LY][2], double plaq, int iter, param_t p){
  
  Complex wLoops[LX/2][LY/2];
  zeroWL(wLoops);

  double sigma[LX/2];  
  for(int i=0; i<LX/2; i++) sigma[i] = 0.0;
  
  Complex w;
  int p1, p2, dx, dy, x, y;
  double inv_Lsq = 1.0/(LX*LY);

  int loopMax = p.loopMax;

  // Polyakov loops and Creutz ratios are measured for 
  // several values of smearing hits to observe the effect
  // of APE smearing.
  for (int hits = 0; hits <= p.smearIter; hits++) {
    
    // CREUTZ RATIOS
    //----------------------------------------------------------------------------
    if(p.measWL) {
      //Smear the gauge field
      Complex smeared[LX][LY][2];
      smearLink(smeared, gauge, p);    
      
      //Loop over all X side sizes of rectangle 
      for(int Xrect=1; Xrect<loopMax; Xrect++) {
	
	//Loop over all Y side sizes of rectangle
	for(int Yrect=1; Yrect<loopMax; Yrect++) {
	  
	  //Loop over all x,y starting points
	  for(x=0; x<LX; x++)
	    for(y=0; y<LY; y++){
	      
	      w = Complex(1.0,0.0);
	      
	      //Move in +x up to p1.
	      for(dx=0; dx<Xrect; dx++)     w *= smeared[ (x+dx)%LX ][y][0];
	      
	      //Move in +y up to p2 (p1 constant)
	      p1 = (x + Xrect)%LX;
	      for(dy=0; dy<Yrect; dy++)     w *= smeared[p1][ (y+dy)%LY ][1];
	      
	      //Move in -x from p1 to (p2 constant)
	      p2 = (y + Yrect)%LY;
	      for(dx=Xrect-1; dx>=0; dx--)  w *= conj(smeared[ (x+dx)%LX ][p2][0]);
	      
	      //Move in -y from p2 to y
	      for(dy=Yrect-1; dy>=0; dy--)  w *= conj(smeared[x][ (y+dy)%LY ][1]);
	      wLoops[Xrect][Yrect] += w*inv_Lsq;
	    }
	}
      }
      
      string name;
      char fname[256];
      FILE *fp;
      
      for(int sizex = 1; sizex<loopMax; sizex++)
	for(int sizey = (sizex == 1 ? 1 : sizex-1); (sizey < loopMax && sizey <= sizex+1); sizey++) {
	  name = "data/rect/rectWL";
	  name += "_hits" + to_string(hits) + "_" + to_string(sizex) + "_" + to_string(sizey);
	  constructName(name, p);
	  name += ".dat";
	  sprintf(fname, "%s", name.c_str());
	  fp = fopen(fname, "a");
	  fprintf(fp, "%d %.16e %.16e\n", iter+1, real(wLoops[sizex][sizey]), imag(wLoops[sizex][sizey]));	    
	  fclose(fp);
	}
    }
    
    // POLYAKOV LOOPS
    //----------------------------------------------------------------------------
    if(p.measPL) {
      Complex pLoops[LX/2];
      Complex loops[LX];
      for(int x=0; x<LX/2; x++) pLoops[x] = 0.0;
      for(int x=0; x<LX; x++) loops[x] = Complex(1.0,0.0);
      
      //Eack polyakov loop correlation is defined by its delta x value.
      //We start at x0, separate to x0 + (x0+L/2-1), and loop over all
      //x0=1 -> x0 = L/2-1.
      
      //We compute all loops at each value of x. Starting x points are 
      //irrelevant for closed temporal loops
      for(int x=0; x<LX; x++) 
	for(int y=0; y<LY; y++) 
	  loops[x] *= gauge[x][y][1];
      
      //Now we have all possible loops, we compute all possible combinations
      //(all possible x) * (all possible dx)
      
      //Loop1
      for(int x1=0; x1<LX; x1++) {
	
	//Loop2
	for(dx=-LX/2; dx<LX/2; dx++) {
	  int x2 = ((x1 + dx) + LX)%LX;
	  
	  if (abs(dx) == LX/2) pLoops[abs(dx)] += conj(loops[x1])*loops[x2];
	  else pLoops[abs(dx)] += conj(loops[x1])*loops[x2]/2.0;
	}
      }

      string name;
      char fname[256];
      FILE *fp;
      
      name= "data/polyakov/polyakov";
      name += "_hits" + to_string(hits);
      constructName(name, p);
      name += ".dat";
      sprintf(fname, "%s", name.c_str());
      
      fp = fopen(fname, "a");
      fprintf(fp, "%d ", iter+1);
      for(dx=0; dx<LX/2; dx++)
	fprintf(fp, "%.16e %.16e ",
		real(pLoops[dx]),
		imag(pLoops[dx]));
      fprintf(fp, "\n");
      fclose(fp);
      
    }
  }
  
  return;
}

//Pion correlation function
//                              |----------------|
//                              |        |-------|---------|
//  < pi(x) | pi(0) > = < ReTr[dn*(x) g3 up(x) | dn(0) g3 up*(0)] >     
//                    = < ReTr(g3 Gd[0,x] g3 Gu[x,0]) >  
//                    = < ReTr(G*[x,0] G[x,0]) >
//
// using g3 G[x,0] g3 = G*[x,0] and Gu[x,0] \eq Gd[x,0]
//
// if H = Hdag, Tr(H * Hdag) = Sum_{n,m} (H_{n,m}) * (H_{n,m})^*,
// i.e., the sum of the modulus squared of each element

void measPionCorrelation(Complex gauge[LX][LY][2], int top, int iter, param_t p) {
  
  //Up type fermion prop
  Complex propUp[LX][LY][2];
  //Down type fermion prop
  Complex propDn[LX][LY][2];
  //fermion prop CG guess
  Complex propGuess[LX][LY][2];
  //Deflation eigenvectors
  Complex defl_evecs[NEV][LX][LY][2];
  //Deflation eigenvalues
  Complex defl_evals[NEV];
      
  double pion_corr[LY];
            
  Complex source[LX][LY][2];
  Complex Dsource[LX][LY][2];

  char fname[256];
  string name;
  FILE *fp;
  
  //Deflate if requested
#ifdef USE_ARPACK
  arpack_solve(gauge, defl_evecs, defl_evals, 0, 0, p);
#endif
  
  //Up type source
  zeroField(source);
  zeroField(Dsource);
  zeroField(propUp);
  zeroField(propGuess);
  source[0][0][0] = cUnit;
  
  // up -> (g3Dg3) * up *****
  // (g3Dg3D)^-1 * (g3Dg3) up = D^-1 * up *****

  g3psi(source);
  g3Dpsi(Dsource, source, gauge, p);  

  if (p.deflate) deflate(propGuess, Dsource, defl_evecs, defl_evals, p);
  Ainvpsi(propUp, Dsource, propGuess, gauge, p);

  //Down type source
  zeroField(source);
  zeroField(Dsource);
  zeroField(propDn);
  zeroField(propGuess);
  source[0][0][1] = cUnit;	    
      
  // dn -> (g3Dg3) * dn *****
  // (g3Dg3D)^-1 * (g3Dg3) dn = D^-1 * dn ***** 

  g3psi(source);
  g3Dpsi(Dsource, source, gauge, p);

  if (p.deflate) deflate(propGuess, Dsource, defl_evecs, defl_evals, p);
  Ainvpsi(propDn, Dsource, propGuess, gauge, p);

  //Get estimate of vacuum trace
  Complex q[2] = {0.0,0.0};
  q[0] = propUp[0][0][0];
  q[1] = propDn[0][0][1];

  double vacuum_trace = (q[0] - q[1]).real();
  
  name = "data/vacuum/estimate_vacuum_Q" + std::to_string(abs(top));
  constructName(name, p);
  name += ".dat";
  sprintf(fname, "%s", name.c_str());
  fp = fopen(fname, "a");
  fprintf(fp, "%d ", iter+1);
  fprintf(fp, "%.16e\n", vacuum_trace);
  fclose(fp);
  
  //Let y be the 'time' dimension
  double corr = 0.0, tmp = 0.0;
  for(int y=0; y<LY; y++) {
    //initialise
    pion_corr[y] = 0.0;
    //Loop over space and spin, fold propagator
    corr = 0.0;
    for(int x=0; x<LX; x++) {
      tmp = abs((conj(propDn[x][y][0]) * propDn[x][y][0] +
      		 conj(propDn[x][y][1]) * propDn[x][y][1] +
      		 conj(propUp[x][y][0]) * propUp[x][y][0] +
      		 conj(propUp[x][y][1]) * propUp[x][y][1]));
      
      corr += tmp;
    }
    
    //Compute folded propagator
    if ( y < ((LY/2)+1) ) pion_corr[y] += corr;
    else {
      pion_corr[LY-y] += corr;
      pion_corr[LY-y] /= 2.0;
    }
  }
  
  //topological sector pion correlation
  name = "data/pion/pion_Q" + std::to_string(abs(top));
  constructName(name, p);
  name += ".dat";  
  sprintf(fname, "%s", name.c_str());
  fp = fopen(fname, "a");
  fprintf(fp, "%d ", iter+1);
  for(int t=0; t<LY/2+1; t++)
    fprintf(fp, "%.16e ", pion_corr[t]);
  fprintf(fp, "\n");
  fclose(fp);

  //Full pion correlation
  name = "data/pion/pion";
  constructName(name, p);
  name += ".dat";  
  sprintf(fname, "%s", name.c_str());
  fp = fopen(fname, "a");
  fprintf(fp, "%d ", iter+1);
  for(int t=0; t<LY/2+1; t++)
    fprintf(fp, "%.16e ", pion_corr[t]);
  fprintf(fp, "\n");
  fclose(fp);

}

void measVacuumTrace(Complex gauge[LX][LY][2], int top, int iter, param_t p) {
  
  //Up type fermion prop
  Complex propUp[LX][LY][2];
  //Down type fermion prop
  Complex propDn[LX][LY][2];
  //fermion prop CG guess
  Complex propGuess[LX][LY][2];
  //Deflation eigenvectors
  Complex defl_evecs[NEV][LX][LY][2];
  //Deflation eigenvalues
  Complex defl_evals[NEV];
  
  double vacuum_trace[2] = {0.0, 0.0};
  
  Complex source[LX][LY][2];
  Complex Dsource[LX][LY][2];
  
  //Disconnected
  //Loop over time slices
  for(int y=0; y<LY; y++) {
    
    //Sum over spatial sites
    for(int x=0; x<LX; x++) {
      
      //Up type source
      zeroField(source);
      zeroField(Dsource);
      zeroField(propUp);
      source[x][y][0] = cUnit;
      
      g3psi(source);
      g3Dpsi(Dsource, source, gauge, p);
      if (p.deflate) deflate(propGuess, Dsource, defl_evecs, defl_evals, p);
      Ainvpsi(propUp, Dsource, propGuess, gauge, p);
      
      //Down type source
      zeroField(source);
      zeroField(Dsource);
      zeroField(propDn);
      source[x][y][1] = cUnit;
      
      g3psi(source);
      g3Dpsi(Dsource, source, gauge, p);
      if (p.deflate) deflate(propGuess, Dsource, defl_evecs, defl_evals, p);
      Ainvpsi(propDn, Dsource, propGuess, gauge, p);
      
      vacuum_trace[0] += (conj(propDn[x][y][0]) * propDn[x][y][0] +
			  conj(propDn[x][y][1]) * propDn[x][y][1] +
			  conj(propUp[x][y][0]) * propUp[x][y][0] +
			  conj(propUp[x][y][1]) * propUp[x][y][1]).real();
      
      vacuum_trace[1] += (conj(propDn[x][y][0]) * propDn[x][y][0] +
			  conj(propDn[x][y][1]) * propDn[x][y][1] +
			  conj(propUp[x][y][0]) * propUp[x][y][0] +
			  conj(propUp[x][y][1]) * propUp[x][y][1]).imag();
      
    }    
  }

  string name = "data/vacuum/vacuum_Q" + std::to_string(abs(top));
  constructName(name, p);
  name += ".dat";
  
  char fname[256];
  sprintf(fname, "%s", name.c_str());
  FILE *fp = fopen(fname, "a");
  fprintf(fp, "%d ", iter+1);
  fprintf(fp, "%.16e %.16e\n", vacuum_trace[0], vacuum_trace[1]);
  fclose(fp);
  
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
