#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>

using namespace std;

#define L 16
#define LX 16
#define LY 16
#define D 2
#define PI 3.141592653589793
#define TWO_PI 6.283185307179586

typedef complex<double> Complex;
#define I Complex(0,1)
#define cUnit Complex(1,0)

#include "utils.h"
#include "latHelpers.h"
#include "measurementHelpers.h"
#include "fermionHelpers.h"
#include "dOpHelpers.h"

#ifdef USE_ARPACK
#include "arpack_interface_staggered.h"
#endif

//HMC routines defined by dimension, so kept in main file
//----------------------------------------------------------------------------
int hmc(Complex gauge[LX][LY][D], param_t p, int iter);

//ACTION DEP
void trajectory(double mom[LX][LY][D], Complex gauge[LX][LY][D],
		Complex chi[LX][LY], param_t p);

//ACTION DEP
void forceD(double fU[LX][LY][D], Complex gauge[LX][LY][D], Complex chi[LX][LY], param_t p);

void forceU(double fU[LX][LY][D], Complex gauge[LX][LY][D], param_t p);
void update_mom(double fU[LX][LY][D], double fD[LX][LY][D], double mom[LX][LY][D], double dtau);
void update_gauge(Complex gauge[LX][LY][D], double mom[LX][LY][D], double dtau);

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
  p.nKv = atoi(argv[16]);
  p.nEv = atoi(argv[17]);
  p.arpackTol = atof(argv[18]);
  p.arpackMaxiter = atoi(argv[19]);
  p.polyACC = atoi(argv[20]);
  p.amax = atof(argv[21]);
  p.amin = atof(argv[22]);
  p.n_poly = atoi(argv[23]);
  
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
      cout << fixed << setprecision(16) << iter+1 << " "; //Iteration
      cout << time/CLOCKS_PER_SEC << " " << endl;         //Time
    }
    
    for(iter=p.therm; iter<2*p.therm; iter++){  
      //Perform HMC step with accept/reject
      accept = hmc(gauge, p, iter);
      double time = time0 + clock();
      cout << fixed << setprecision(16) << iter+1 << " "; //Iteration
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
    //------------------------------
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
    //--------------------
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
      for(int a=0; a<L/2; a++) pLoops[a] = 0.0;
      measPolyakovLoops(gauge, pLoops);
      
      name = "data/polyakov/polyakov";
      constructName(name, p);
      name += ".dat";
      sprintf(fname, "%s", name.c_str());
      fp = fopen(fname, "a");
      fprintf(fp, "%d ", iter + 1);
      for(int size=1 ; size < L/2; size++)
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
      for(int size=1 ; size < L/2-1; size++)
	fprintf(fp, "%.16e ",
		real(pLoops[size+1])/real(pLoops[size]));
      fprintf(fp, "\n");
      fclose(fp);
      
      // Creutz Ratios
      zeroWL(wLoops);
      measWilsonLoops(gauge, wLoops, p);
      
      //Compute string tension
      for(int size=1; size<L/2; size++) {
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
      for(int size=2 ; size < L/2; size++)
	fprintf(fp, "%.16e ", sigma[size]);
      fprintf(fp, "\n");
      fclose(fp);
      
      for(int sizex=2; sizex<L/2; sizex++)
	for(int sizey=sizex-1; (sizey < L/2 && sizey <= sizex+1); sizey++) {
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

    }
  }
  
  return 0;
}

int hmc(Complex gauge[LX][LY][D], param_t p, int iter) {

  int accept = 0;
  
  double mom[LX][LY][D];
  Complex gaugeOld[LX][LY][D];
  Complex phi[LX][LY], chi[LX][LY];
  double H, Hold;

  copyLat(gaugeOld,gauge);
  zeroLat(mom);
  gaussReal_F(mom); 

  H = 0.0;
  Hold = 0.0;
  
  if(p.dynamic) {
    //Create gaussian distributed fermion field chi. chi[LX][LY] E exp(-chi^* chi)
    gaussComplex_F(chi, p);
    //Create pseudo fermion field phi = D chi
    Dpsi(phi, chi, gauge, p);    

    //Mask out odd sites.
    for(int x=0; x<LX; x++)  
      for(int y=0; y<LY; y++)
	if((x+y)%2 == 1) phi[x][y] = 0.0;
  }

  Hold = measAction(mom, gauge, phi, p, false);
  trajectory(mom, gauge, phi, p); // MD trajectory using Verlet
  H = measAction(mom, gauge, phi, p, true);

  if (iter >= 2*p.therm) {      
    expcount++;
    expdHAve += exp(-(H-Hold));
    dHAve += (H-Hold);
  }
  
  // Metropolis accept/reject step
  if (iter >= p.therm) {    
    if ( drand48() > exp(-(H-Hold)) ) copyLat(gauge,gaugeOld);
    else accept = 1;
  }
  
  return accept;
}

void trajectory(double mom[LX][LY][D], Complex gauge[LX][LY][D],
		Complex phi[LX][LY], param_t p) {

  double dtau = p.tau/p.nstep;
  double H = 0.0;
  Complex guess[LX][LY];
#ifdef USE_ARPACK
  gaussComplex_F(guess, p);
  //for(int x=0; x<LX; x++)
  //for(int y=0; y<LY; y++)
  //guess[x][y] = drand48();
#endif
  
  //gauge force
  double fU[LX][LY][D];
  zeroLat(fU);
  //fermion fermion
  double fD[LX][LY][D];
  zeroLat(fD);
  //Both arrays are zeroed in forceU/D function call
  
  //Initial half step.
  //P_{1/2} = P_0 - dtau/2 * (fU + fD)
  forceU(fU, gauge, p);
  forceD(fD, gauge, phi, p);
  update_mom(fU, fD, mom, 0.5*dtau);  
  
  for(int k=1; k<p.nstep; k++) {
    
    //U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
    update_gauge(gauge, mom, dtau);
    
    //P_{k+1/2} = P_{k-1/2} - dtau * (fU + fD)
    forceU(fU, gauge, p);
    forceD(fD, gauge, phi, p);
    update_mom(fU, fD, mom, dtau);  
       
#ifdef USE_ARPACK
    int ARPACK_iter = arpack_solve_double(gauge, p, guess, 1, k, 0);
#endif
    //H = measAction(mom, gauge, phi, p, true);
  }
  
  //Final half step.
  //U_{n} = exp(i dtau P_{n-1/2}) * U_{n-1}
  update_gauge(gauge, mom, dtau);
  
  //P_{n} = P_{n-1/2} - dtau/2 * (fU + fD)
  forceU(fU, gauge, p);
  forceD(fD, gauge, phi, p);
  update_mom(fU, fD, mom, 0.5*dtau);
  
  //trajectory complete
  
}

void forceU(double fU[LX][LY][D],Complex gauge[LX][LY][D], param_t p) {

  Complex plaq0;
  Complex plaq;
  zeroLat(fU);
  int xp1, xm1, yp1, ym1;
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++) {

      xp1 = (x+1)%LX;
      yp1 = (y+1)%LY;
      xm1 = (x-1+LX)%LX;
      ym1 = (y-1+LY)%LY;
      
      plaq0 = gauge[x][y][0]*gauge[xp1][y][1]*conj(gauge[x][yp1][0])*conj(gauge[x][y][1]);
      fU[x][y][0] += p.beta*imag(plaq0);
      
      plaq =  gauge[x][ym1][0]*gauge[xp1][ym1][1]*conj(gauge[x][y][0])*conj(gauge[x][ym1][1]);
      fU[x][y][0] -= p.beta*imag(plaq);

      plaq = gauge[x][y][1]*conj(gauge[xm1][yp1][0])*conj(gauge[xm1][y][1])*gauge[xm1][y][0];
      fU[x][y][1] += p.beta*imag(plaq);

      //This plaquette was aleady computed.
      fU[x][y][1] -= p.beta*imag(plaq0);
  
      //cout << x  << y <<" Forces = "<< fU[x][y][0] << "   " <<fU[x][y][1] << endl;
    }
}




// VD(theta) = chi*_e (1/D D^dag) chi_e
//  
// D chi = phi   and phi gausian.
//  
// fD = - \dd_theta VD(theta) 
//    =   chi*_e (1/D D^dag)_ee [  D^\dag (\dd_theta D) + (\dd_theta D^\dag) D ]_ee
//                                (1/D D^dag)_ee chi_e chi_e 
//  
// Should optimze this to operate only on EVEN sites.
// chi_even, chitmp_even, Dchitmp_odd
  
void forceD(double fD[LX][LY][D], Complex gauge[LX][LY][D], Complex chi[LX][LY], param_t p) {

  if(p.dynamic) {
    Complex chitmp[LX][LY];
    Complex Dchitmp[LX][LY];
    zeroField(chitmp);
    zeroField(Dchitmp);
    zeroLat(fD);
    
    Ainv_psi(chitmp, chi, chitmp, gauge, p); // note chitmp = 0 for ODD
    Dpsi(Dchitmp, chitmp, gauge, p); // restrict to Dslash, m = 0
    
    for(int x=0; x<LX; x++)
      for(int y=0; y<LY; y++){
	if((x+y)%2 == 1)  chitmp[x][y] = Complex(0.0,0.0);
	if((x+y)%2 == 0) Dchitmp[x][y] = Complex(0.0,0.0);
      }
    
    double eta1;
    for(int x=0; x<LX; x++)
      for(int y=0; y<LY; y++) {
	eta1 =(1.0 - 2.0*(x%2));

	//Odd site
	if((x+y+1)%2 == 0){ 
	  fD[x][y][0] += 2.0*imag(conj(Dchitmp[x][y]) * gauge[x][y][0] * chitmp[(x+1)%LX][y]);
	  fD[x][y][1] += 2.0*eta1*imag(conj(Dchitmp[x][y]) * gauge[x][y][1] * chitmp[x][(y+1)%LY]);
	}
	else {
	  //even site
	  fD[x][y][0] += 2.0*imag(conj(Dchitmp[(x+1)%LX][y]) * conj(gauge[x][y][0]) * chitmp[x][y]);
	  fD[x][y][1] += 2.0*eta1*imag(conj(Dchitmp[x][(y+1)%LY]) * conj(gauge[x][y][1]) * chitmp[x][y]);
	}	
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
