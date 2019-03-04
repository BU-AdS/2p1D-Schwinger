#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>

using namespace std;

#define LX 16
#define LY 32
#define D 2
#define NEV 24
#define NKR 32
#define PI 3.141592653589793
#define TWO_PI 6.283185307179586

typedef complex<double> Complex;
#define I Complex(0,1.0)
#define cUnit Complex(1.0,0)

#include "utils.h"
#include "latHelpers.h"
#include "measurementHelpers.h"
#include "fermionHelpers.h"
#include "dOpHelpers.h"

#ifdef USE_ARPACK
#include "arpack_interface_wilson.h"
#endif

//HMC routines defined by dimension, so kept in main file
//----------------------------------------------------------------------------

//Fermion dependent
void trajectory(double mom[LX][LY][D], Complex gauge[LX][LY][D],
		Complex phi[LX][LY][2], param_t p, int iter);
void forceD(double fU[LX][LY][D], Complex gauge[LX][LY][D], Complex phi[LX][LY][2],
	    Complex guess[LX][LY][2], param_t p);

//Fermion agnostic
int hmc(Complex gauge[LX][LY][D], param_t p, int iter);
void forceU(double fU[LX][LY][D], Complex gauge[LX][LY][D], param_t p);
void update_mom(double fU[LX][LY][D], double fD[LX][LY][D], double mom[LX][LY][D], double dtau);
void update_gauge(Complex gauge[LX][LY][D], double mom[LX][LY][D], double dtau);
void doParityFlip(Complex gauge[LX][LY][D], param_t p);
//----------------------------------------------------------------------------

//Global variables.
int expcount = 0;
double expdHAve = 0.0;
double dHAve = 0.0;

int main(int argc, char **argv) {

  param_t p;
  
  p.beta = atof(argv[1]); 
  p.iterHMC = atoi(argv[2]);
  p.therm = atoi(argv[3]);
  p.skip = atoi(argv[4]);
  p.chkpt = atoi(argv[5]);
  p.checkpointStart = atoi(argv[6]);  
  p.nstep = atoi(argv[7]);
  p.tau = atof(argv[8]);
  
  p.smearIter = atoi(argv[9]);
  p.alpha = atof(argv[10]);  
  long iseed = (long)atoi(argv[11]);
  //Pseudo RNG seed
  srand48(iseed);
  
  if(atoi(argv[12]) == 0) 
    p.dynamic = false;
  else
    p.dynamic = true;

  p.m = atof(argv[13]);
  
  p.maxIterCG = atoi(argv[14]);
  p.eps = atof(argv[15]);
  
  //Arpack params
  p.nKr = NKR;
  p.nEv = NEV;
  p.arpackTol = atof(argv[16]);
  p.arpackMaxiter = atoi(argv[17]);
  p.polyACC = atoi(argv[18]);
  p.amax = atof(argv[19]);
  p.amin = atof(argv[20]);
  p.n_poly = atoi(argv[21]);
  
  //Topology
  double top = 0.0;
  int top_int = 0;
  int top_old = 0;
  int top_stuck = 0;

  int histL = 101;
  int histQ[histL];
  double plaqSum = 0.0;
  int index = 0;
  for(int i = 0; i < histL; i++) histQ[i] = 0;
  
  Complex gauge[LX][LY][D];
  zeroLat(gauge);  
  Complex gaugeFree[LX][LY][D];
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++) {
      gaugeFree[LX][LY][0] = cUnit;
      gaugeFree[LX][LY][1] = cUnit;
    }
    
  double sigma[LX/2];
  Complex pLoops[LX/2];
  Complex wLoops[LX/2][LY/2];  
  for(int x=0; x<LX/2; x++) {
    pLoops[x] = 0.0;
    sigma[x] = 0.0;
    for(int y=0; y<LY/2; y++) {
      wLoops[x][y] = Complex(0.0,0.0);   // for average Creutz ratios
    }
  }

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
  double vacuum_trace[2] = {0.0, 0.0};
  
  int count = 0;
  string name;
  fstream outPutFile;
  
  int accept;
  int accepted = 0;
  char fname[256];
  FILE *fp;

  printParams(p);  
  gaussStart(gauge,p);  // hot start

  //Start simulation
  double time0 = -((double)clock());
  int iter_offset = 0;
  int iter = 0;
  cout << setprecision(16);
  
  if(p.checkpointStart > 0) {

    //Read in gauge field if requested
    //---------------------------------------------------------------------
    name = "gauge/gauge";
    constructName(name, p);
    name += "_traj" + to_string(p.checkpointStart) + ".dat";	
    readGaugeLattice(gauge,name);
    iter_offset = p.checkpointStart;    
  } else {

    //Thermalise from random start
    //---------------------------------------------------------------------
    for(iter=0; iter<p.therm; iter++){  
      //Perform HMC step
      accept = hmc(gauge, p, iter);
      double time = time0 + clock();
      cout << fixed << iter+1 << " "; //Iteration
      cout << time/CLOCKS_PER_SEC << " " << endl;         //Time
    }
    
    for(iter=p.therm; iter<2*p.therm; iter++){  
      //Perform HMC step with accept/reject
      accept = hmc(gauge, p, iter);
      double time = time0 + clock();
      cout << fixed << iter+1 << " "; //Iteration
      cout << time/CLOCKS_PER_SEC << " " << endl;         //Time
    }
    iter_offset = 2*p.therm;    
  }
  
  for(iter=iter_offset; iter<p.iterHMC + iter_offset; iter++){
    
    //Perform HMC step
    accept = hmc(gauge, p, iter);
    
    //HMC acceptance
    accepted += accept;

    //Measure the topological charge
    //---------------------------------------------------------------------
    top = measTopCharge(gauge, p);
    top_int = round(top);
    name = "data/top/top_charge";
    constructName(name, p);
    name += ".dat";
    sprintf(fname, "%s", name.c_str());
    fp = fopen(fname, "a");
    fprintf(fp, "%d %d\n", iter, top_int);
    fclose(fp);

    index = top_int + (histL-1)/2;
    histQ[index]++;
    if(top_old == top_int) top_stuck++;
    top_old = top_int;
    
    //Perform Measurements
    //---------------------------------------------------------------------
    if( (iter+1)%p.skip == 0) {
      
      count++;

      //Checkpoint the gauge field
      if( (iter+1)%p.chkpt == 0) {	  
	name = "gauge/gauge";
	constructName(name, p);
	name += "_traj" + to_string(iter+1) + ".dat";
	writeGaugeLattice(gauge, name);
      }
      
      //Plaquette action
      double plaq = measPlaq(gauge);
      plaqSum += plaq;
	
      //Info dumped to stdout
      double time = time0 + clock();
      cout << fixed << setprecision(16) << iter+1 << " "; //Iteration
      cout << time/CLOCKS_PER_SEC << " ";                 //Time
      cout << plaqSum/count << " ";                       //Action
      cout << (double)top_stuck/(count*p.skip) << " ";    //P(stuck)
      cout << expdHAve/expcount << " ";                   //Average exp(-dH)
      cout << dHAve/expcount << " ";                      //Average dH
      cout << (double)accepted/(count*p.skip) << " ";     //Acceptance
      cout << top_int << endl;                            //T charge
	
      //Dump to file
      name = "data/data/data";//I cannot make bricks without clay!
      constructName(name, p);
      name += ".dat";	
      sprintf(fname, "%s", name.c_str());	
      fp = fopen(fname, "a");
	
      fprintf(fp, "%d %.16e %.16e %.16e %.16e %.16e\n",
	      iter+1,
	      time/CLOCKS_PER_SEC,
	      plaqSum/count,
	      (double)top_stuck/(count*p.skip),
	      expdHAve/expcount,
	      (double)accepted/(count*p.skip));
      fclose(fp);
	
      //Polyakov wLoops      
      for(int x=0; x<LX/2; x++) pLoops[x] = 0.0;
      measPolyakovLoops(gauge, pLoops);
      
      name = "data/polyakov/polyakov";
      constructName(name, p);
      name += ".dat";
      sprintf(fname, "%s", name.c_str());
      fp = fopen(fname, "a");
      fprintf(fp, "%d ", iter + 1);
      for(int size=1; size<LX/2; size++)
	fprintf(fp, "%.16e %.16e ",
		real(pLoops[size]),
		imag(pLoops[size]) );
      fprintf(fp, "\n");
      fclose(fp);

      name = "data/polyakov/polyakov_ratios";
      constructName(name, p);
      name += ".dat";
      sprintf(fname, "%s", name.c_str());
      fp = fopen(fname, "a");
      fprintf(fp, "%d ", iter + 1);
      for(int size=1 ; size < LX/2-1; size++)
	fprintf(fp, "%.16e ",
		real(pLoops[size+1])/real(pLoops[size]));
      fprintf(fp, "\n");
      fclose(fp);
      
      // Creutz Ratios
      zeroWL(wLoops);
      measWilsonLoops(gauge, wLoops, p);
      
      //Compute string tension
      for(int size=1; size<LX/2; size++) {
	sigma[size]  = - log(abs((real(wLoops[size][size])/real(wLoops[size-1][size]))* 
				(real(wLoops[size-1][size-1])/real(wLoops[size][size-1]))));
	
	sigma[size] += - log(abs((real(wLoops[size][size])/real(wLoops[size][size-1]))* 
				 (real(wLoops[size-1][size-1])/real(wLoops[size-1][size]))));
	
	sigma[size] *= 0.5;
	
      }
      
      name = "data/creutz/creutz";
      constructName(name, p);
      name += ".dat";
      sprintf(fname, "%s", name.c_str());
      fp = fopen(fname, "a");
      fprintf(fp, "%d %.16e ", iter+1, -log(abs(plaq)) );
      for(int size=2; size<LX/2; size++)
	fprintf(fp, "%.16e ", sigma[size]);
      fprintf(fp, "\n");
      fclose(fp);
      
      for(int sizex=2; sizex<LX/2; sizex++)
	for(int sizey=sizex-1; (sizey < LY/2 && sizey <= sizex+1); sizey++) {
	  name = "data/rect/rectWL";
	  name += "_" + to_string(sizex) + "_" + to_string(sizey);
	  constructName(name, p);
	  name += ".dat";
	  sprintf(fname, "%s", name.c_str());
	  fp = fopen(fname, "a");
	  fprintf(fp, "%d %.16e %.16e\n", iter+1, real(wLoops[sizex][sizey]), imag(wLoops[sizex][sizey]));	    
	  fclose(fp);
	}

      //Update topoligical charge histogram
      name = "data/top/top_hist";
      constructName(name, p);
      name += ".dat";
      sprintf(fname, "%s", name.c_str());
      fp = fopen(fname, "w");
      for(int i=0; i<histL; i++) fprintf(fp, "%d %d\n", i - (histL-1)/2, histQ[i]);
      fclose(fp);
      
      //Pion correlation function
      //                              |----------------|
      //                              |        |-------|---------|
      //  < pi(x) | pi(0) > = < ReTr[up(x) g5 dn*(x) | up*(0) g5 dn(0)] >     
      //                    = < ReTr(G[x,0] G*[x,0]) >
      //
      // if H = Hdag, Tr(H * Hdag) = Sum_{n,m} (H_{n,m}) * (H_{n,m})^*,
      // i.e., the sum of the modulus squared of each element
      
      Complex source[LX][LY][2];
      Complex Dsource[LX][LY][2];

      //Up type source
      zeroField(source);
      zeroField(Dsource);
      zeroField(propUp);
      zeroField(propGuess);
      source[0][0][0] = cUnit;
      // up -> (g5Dg5) * up
      g5psi(source);
      g5Dpsi(Dsource, source, gauge, p);

#ifdef USE_ARPACK
      arpack_solve(gauge, defl_evecs, defl_evals, 0, 0, p);
#endif
      deflate(propGuess, Dsource, defl_evecs, defl_evals, p);
      // (g5Dg5D)^-1 * (g5Dg5) up = D^-1 * up
      Ainvpsi(propUp, Dsource, propGuess, gauge, p);

      //Down type source
      zeroField(source);
      zeroField(Dsource);
      zeroField(propDn);
      source[0][0][1] = cUnit;	    
      
      // dn -> (g5Dg5) * dn
      g5psi(source);
      g5Dpsi(Dsource, source, gauge, p);

      deflate(propGuess, Dsource, defl_evecs, defl_evals, p);
      // (g5Dg5D)^-1 * (g5Dg5) dn = D^-1 * dn
      Ainvpsi(propDn, Dsource, propGuess, gauge, p);
      
      //Let y be the 'time' dimension
      double corr = 0.0;
      for(int y=0; y<LY; y++) {
	//initialise
	pion_corr[y] = 0.0;
	//Loop over space and spin, fold propagator
	for(int x=0; x<LX; x++)
	  corr = (conj(propDn[x][y][0]) * propDn[x][y][0] +
		  conj(propDn[x][y][1]) * propDn[x][y][1] +
		  conj(propUp[x][y][0]) * propUp[x][y][0] +
		  conj(propUp[x][y][1]) * propUp[x][y][1]).real();
	
	if(y<LY/2+1) pion_corr[y] += corr;
	else pion_corr[LY-y] += corr;	    
      }
      
      
      //pion Correlation
      name = "data/pion/pion";
      constructName(name, p);
      name += ".dat";
      sprintf(fname, "%s", name.c_str());
      fp = fopen(fname, "a");
      fprintf(fp, "%d ", iter+1);
      for(int t=0; t<LY/2; t++)
	fprintf(fp, "%.16e ", pion_corr[t]);
      fprintf(fp, "\n");
      fclose(fp);

      //Disconnected
      //Up type source
      zeroField(source);
      zeroField(Dsource);
      zeroField(propUp);
      for(int x=0; x<LX; x++) {
	for(int y=0; y<LY; y++) {
	  source[x][1][0] = cUnit;

	  g5psi(source);
	  g5Dpsi(Dsource, source, gauge, p);
	  deflate(propGuess, Dsource, defl_evecs, defl_evals, p);
	  Ainvpsi(propUp, Dsource, propGuess, gauge, p);
	  
	  //Down type source
	  zeroField(source);
	  zeroField(Dsource);
	  zeroField(propDn);
	  source[x][y][1] = cUnit;

	  g5psi(source);
	  g5Dpsi(Dsource, source, gauge, p);
	  deflate(propGuess, Dsource, defl_evecs, defl_evals, p);
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
      cout << "Vacuum trace = (" << vacuum_trace[0]/(count*LX*LY) << "," << vacuum_trace[1]/(count*LX*LX) << ")" << endl;
    }
  }
  return 0;
}


// HMC Routines
//---------------------------------------------------------------------

int hmc(Complex gauge[LX][LY][D], param_t p, int iter) {
  
  int accept = 0;
  
  double mom[LX][LY][2];
  Complex gaugeOld[LX][LY][2];
  Complex phi[LX][LY][2], chi[LX][LY][2];
  double H, Hold;

  copyLat(gaugeOld, gauge);
  zeroLat(mom); 
  zeroField(phi);
  zeroField(chi);
  H = 0.0;
  Hold = 0.0;

  // init mom[LX][LY][D]  <mom^2> = 1;
  gaussReal_F(mom); 
  
  if(p.dynamic == true) {    
    //Create gaussian distributed fermion field chi. chi[LX][LY] E exp(-chi^* chi)
    gaussComplex_F(chi, p);
    //Create pseudo fermion field phi = D chi
    g5Dpsi(phi, chi, gauge, p);    
  }
  
  if (iter >= p.therm) Hold = measAction(mom, gauge, chi, p, false);
  trajectory(mom, gauge, phi, p, iter);
  if (iter >= p.therm) H = measAction(mom, gauge, phi, p, true);
  
  if (iter >= 2*p.therm) {      
    expcount++;
    expdHAve += exp(-(H-Hold));
    dHAve += (H-Hold);
  }

  // Metropolis accept/reject step
  if (iter >= p.therm) {    
    if ( drand48() > exp(-(H-Hold)) ) copyLat(gauge, gaugeOld);
    else accept = 1;
  }
  
  return accept;
}

void trajectory(double mom[LX][LY][2], Complex gauge[LX][LY][2],
		Complex phi[LX][LY][2], param_t p, int iter) {  
  
  double dtau = p.tau/p.nstep;
  double H = 0.0;
  Complex guess[LX][LY][2];
#ifdef USE_ARPACK
  //zeroField(guess);
  
  //deflate using phi as source
  //Deflation eigenvectors
  Complex defl_evecs[NEV][LX][LY][2];
  //Deflation eigenvalues
  Complex defl_evals[NEV];
  
  copyField(guess, phi);
  arpack_solve(gauge, defl_evecs, defl_evals, 0, 0, p);
  deflate(guess, phi, defl_evecs, defl_evals, p);
  
#else
  zeroField(guess);
#endif
  
  //gauge force
  double fU[LX][LY][D];
  //fermion fermion
  double fD[LX][LY][D];
  //Both arrays are zeroed in forceU/D function call
  
  //Initial half step.
  //P_{1/2} = P_0 - dtau/2 * (fU - fD)
  forceU(fU, gauge, p);
  forceD(fD, gauge, phi, guess, p);
  update_mom(fU, fD, mom, 0.5*dtau);  
  
  for(int k=1; k<p.nstep; k++) {
    
    //U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
    update_gauge(gauge, mom, dtau);
    
    //P_{k+1/2} = P_{k-1/2} - dtau * (fU - fD)
    forceU(fU, gauge, p);
    forceD(fD, gauge, phi, guess, p);
    update_mom(fU, fD, mom, dtau);
  }
  
  //Final half step.
  //U_{n} = exp(i dtau P_{n-1/2}) * U_{n-1}
  update_gauge(gauge, mom, dtau);
  
  //P_{n} = P_{n-1/2} - dtau/2 * (fU - fD)
  forceU(fU, gauge, p);
  forceD(fD, gauge, phi, guess, p);
  update_mom(fU, fD, mom, 0.5*dtau);
  
  //trajectory complete
}

void forceU(double fU[LX][LY][2], Complex gauge[LX][LY][2], param_t p) {
  
  Complex plaq0;
  Complex plaq;
  zeroLat(fU);
  int xp1, xm1, yp1, ym1;
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++) {

      xp1 = (x+1)%LX;
      xm1 = (x-1+LX)%LX;
      yp1 = (y+1)%LY;
      ym1 = (y-1+LY)%LY;

      
      plaq0 = gauge[x][y][0]*gauge[xp1][y][1]*conj(gauge[x][yp1][0])*conj(gauge[x][y][1]);
      fU[x][y][0] += p.beta*imag(plaq0);
      
      plaq =  gauge[x][ym1][0]*gauge[xp1][ym1][1]*conj(gauge[x][y][0])*conj(gauge[x][ym1][1]);
      fU[x][y][0] -= p.beta*imag(plaq);

      plaq =  gauge[x][y][1]*conj(gauge[xm1][yp1][0])*conj(gauge[xm1][y][1])*gauge[xm1][y][0];
      fU[x][y][1] += p.beta*imag(plaq);

      //This plaquette was aleady computed. We want the conjugate.
      fU[x][y][1] -= p.beta*imag(plaq0);
      
    }
}

// let dD \equiv (d/dtheta D)
//
// d/dtheta (phi^* (DD^dag)^-1 phi) = -((DD^dag)^1 phi)^dag ([dD]*D^dag + D*[dD^dag]) ((DD^dag)^-1 phi)
//
// *****  Should optimize this to operate only on EVEN sites. ****

void forceD(double fD[LX][LY][2], Complex gauge[LX][LY][2], Complex phi[LX][LY][2], Complex guess[LX][LY][2], param_t p){
  
  if(p.dynamic == true) {

    zeroLat(fD);
    
    //phip = (DD^dag)^-1 * phi
    Complex phip[LX][LY][2];
    zeroField(phip);
    
    //Ainvpsi inverts using the DdagD (g5Dg5D) operator, returns
    // phip = (D^-1 * Ddag^-1) phi = (D^-1 * g5 * D^-1 g5) phi.
    Complex guess[LX][LY][2]; //Initial guess to CG
    zeroField(guess);
    Ainvpsi(phip, phi, guess, gauge, p);
    
    //g5Dphi = g5D * phip
    Complex g5Dphi[LX][LY][2];
    zeroField(g5Dphi);
    g5Dpsi(g5Dphi, phip, gauge, p);
    
    int xp1, xm1, yp1, ym1;
    double r = 1.0;
    for(int x=0; x<LX; x++)
      for(int y=0; y<LY; y++) {

	xp1 = (x+1)%LX;
	yp1 = (y+1)%LY;
	xm1 = (x-1+LX)%LX;
	ym1 = (y-1+LY)%LY;	
	
	//mu = 0
	//upper
	// | r  1 | 
	// | 1  r |
	//lower
	// | r -1 |
	// | 1 -r |					
	fD[x][y][0] += real(I*((conj(gauge[x][y][0]) *
				(conj(phip[xp1][y][0]) * (r*g5Dphi[x][y][0] +   g5Dphi[x][y][1]) -
				 conj(phip[xp1][y][1]) * (  g5Dphi[x][y][0] + r*g5Dphi[x][y][1])))
			       -
			       (gauge[x][y][0] *
				(conj(phip[x][y][0]) * (r*g5Dphi[xp1][y][0] -   g5Dphi[xp1][y][1]) +
				 conj(phip[x][y][1]) * (  g5Dphi[xp1][y][0] - r*g5Dphi[xp1][y][1])))
			       )
			    );	
	
	//mu = 1
	//upper
	// | r -i | 
	// | i  r |
	//lower
	// | r  i |
	// | i -r |
	fD[x][y][1] += real(I*((conj(gauge[x][y][1]) *
				(conj(phip[x][yp1][0]) * (r*g5Dphi[x][y][0] - I*g5Dphi[x][y][1]) -
				 conj(phip[x][yp1][1]) * (I*g5Dphi[x][y][0] + r*g5Dphi[x][y][1])))
			       -			       
			       (gauge[x][y][1] *
				(conj(phip[x][y][0]) * (r*g5Dphi[x][yp1][0] + I*g5Dphi[x][yp1][1]) +
				 conj(phip[x][y][1]) * (I*g5Dphi[x][yp1][0] - r*g5Dphi[x][yp1][1])))
			       )
			    );	
      }
  }
}

//P_{k+1/2} = P_{k-1/2} - dtau * (fU + fD)
void update_mom(double fU[LX][LY][D], double fD[LX][LY][D], double mom[LX][LY][D], double dtau){

  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int mu=0; mu<D; mu++)
	mom[x][y][mu] -= (fU[x][y][mu] - fD[x][y][mu])*dtau;
}

//U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
void update_gauge(Complex gauge[LX][LY][D], double mom[LX][LY][D], double dtau){
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int mu=0; mu<D; mu++)
	gauge[x][y][mu] *= polar(1.0, mom[x][y][mu] * dtau);
}


// Do a parity flip based on eqn 4 of arXiv:1203.2560v2
void doParityFlip(Complex gauge[LX][LY][D], param_t p) {
  Complex parity_gauge[LX][LY][D];
  for (int x=0; x<LX; x++) {
    for (int y=0; y<LY; y++) {
      //parity_gauge[(-x-1+2*L)%L][y][0] = conj(gauge[x][y][0]);
      //parity_gauge[(-x+2*L)%L][y][1] = gauge[x][y][1];
    }
  }
  copyLat(gauge, parity_gauge);
}

//-------------------------------------------------------------------------------
