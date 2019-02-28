#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>
#include <omp.h>

using namespace std;

#define L 32
#define D 2
#define PI 3.141592653589793
#define TWO_PI 6.283185307179586

typedef complex<double> Complex;
#define I Complex(0,1)
#define cUnit Complex(1,0)

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


//Measurements
//----------------------------------------------------------------------------
void areaLaw(Complex gauge[L][L][D], Complex avW[L/2][L/2], param_t p);
void polyakovLoops(Complex  gauge[L][L][D], Complex polyakov[L/2]);
double getTopCharge(Complex gauge[L][L][D], param_t p);


//HMC
//----------------------------------------------------------------------------
int hmc(Complex gauge[L][L][D], param_t p, int iter);
double calcH(double mom[L][L][D], Complex gauge[L][L][D],
	     Complex phi[L][L][2], param_t p, bool postStep);
void trajectory(double mom[L][L][D], Complex gauge[L][L][D],
		Complex phi[L][L][2], param_t p, int iter);
void forceU(double fU[L][L][D], Complex gauge[L][L][D], param_t p);
void forceD(double fU[L][L][D], Complex gauge[L][L][D], Complex phi[L][L][2], param_t p);
void update_mom(double fU[L][L][D], double fD[L][L][D], double mom[L][L][D], double dtau);
void update_gauge(Complex gauge[L][L][D], double mom[L][L][D], double dtau);
//----------------------------------------------------------------------------

void DoPFlip(Complex gauge[L][L][D], param_t p);

#include "utils.h"
#include "linAlgHelpers.h"
#include "dOpHelpers.h"

#ifdef USE_ARPACK
#include "arpack_interface.h"
#endif

int expcount = 0;
double expdHAve = 0.0;
double dHAve = 0.0;

int main(int argc, char **argv) {

  param_t p;
  p.Latsize = L;
  
  p.beta = atof(argv[1]); 
  p.iterHMC = atoi(argv[2]);
  p.therm = atoi(argv[3]);
  p.skip = atoi(argv[4]);
  p.chkpt = atoi(argv[5]);
  p.checkpointStart = atoi(argv[6]);  
  if(p.checkpointStart > 0) p.iterHMC += p.checkpointStart;
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
  
  Complex gauge[L][L][D];
  Complex gaugeFree[L][L][D];
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++) {
      gaugeFree[L][L][0] = cUnit;
      gaugeFree[L][L][1] = cUnit;
    }
  
  Complex avW[L][L], avWc[L/2][L/2];

#ifdef USE_ARPACK
  Complex guess[L][L][2];
  Complex guessDUM[L][L][2];
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	guess[x][y][s] = I;
  
  int icount1 = 0;
  int icount2 = 0;
  int icount3 = 0;
#endif
  
  int histL = 101;
  int histQ[histL];
  for(int i = 0; i < histL; i++) histQ[i] = 0;
  
  double sigma[L];
  Complex polyakov[L/2];
  //Complex polyakovAveTest(0.0,0.0);
  
  for(int x=0; x<L/2; x++) {
    polyakov[x] = 0.0;
    for(int y=0; y<L/2; y++) {
      avW[x][y] = Complex(0.0,0.0);    // for average wilson loops
      avWc[x][y] = Complex(0.0,0.0);   // for average Creutz ratios
    }
  }
  
  int count = 0;
  double plaqSum = 0.0;
  int index = 0;
  string name;
  fstream outPutFile;
  
  int accept;
  int accepted = 0;
  char fname[256];
  FILE *fp;

  //sran2(iseed);
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
    
    //Measure the plaquette and topological charge
    
    //Perform HMC step
    accept = hmc(gauge, p, iter);
    
    //HMC acceptance
    accepted += accept;
      
    //Topology
    top = getTopCharge(gauge,p);
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
      name = "data/data/data";
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
	
      //Compute gauge observables	
      for(int a=0; a<L/2; a++) polyakov[a] = 0.0;
      polyakovLoops(gauge, polyakov);
	
      zeroWL(avWc);
      areaLaw(gauge, avWc, p);
	
      // Creutz Ratios
      for(int size=1; size<L/2; size++) {
	sigma[size]  = - log(abs((real(avWc[size][size])/real(avWc[size-1][size]))* 
				(real(avWc[size-1][size-1])/real(avWc[size][size-1]))));

	sigma[size] += - log(abs((real(avWc[size][size])/real(avWc[size][size-1]))* 
				(real(avWc[size-1][size-1])/real(avWc[size-1][size]))));
	
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
	  fprintf(fp, "%d %.16e %.16e\n", iter+1, real(avWc[sizex][sizey]), imag(avWc[sizex][sizey]));	    
	  fclose(fp);
	}
	
      name = "data/polyakov/polyakov";
      constructName(name, p);
      name += ".dat";
      sprintf(fname, "%s", name.c_str());
      fp = fopen(fname, "a");
      fprintf(fp, "%d ", iter + 1);
      for(int size=1 ; size < L/2; size++)
	fprintf(fp, "%.16e %.16e ",
		real(polyakov[size]),
		imag(polyakov[size]) );
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
		real(polyakov[size+1])/real(polyakov[size]));
      fprintf(fp, "\n");
      fclose(fp);
	
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


//======================================================================
//   Creutz     exp[ -sigma L^2] exp[ -sigma(L-1)(L-1)]
//   ratio:    ---------------------------------------  = exp[ -sigma]
//              exp[ -sigma (L-1)L] exp[-sigma L(L-1)]
//======================================================================

void areaLaw(Complex gauge[L][L][D], Complex avWc[L/2][L/2], param_t p){

  Complex w;
  int p1, p2;
  Complex smeared[L][L][D];
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
	    avWc[Xrect][Yrect] += w*inv_Lsq;
	  }
      }
    }

    /*
  for(int xL=0; xL<L/2; xL++) 
    for(int yL=0; yL<L/2; yL++){
      for(int x=0; x<L; x++)
	for(int y=0; y<L; y++){
	  w = Complex(1.0,0.0);
	  
	  for(int dx=0; dx<xL; dx++)     w *=      Smeared[ (x+dx)%L ][y][0];
	  for(int dy=0; dy<yL; dy++)     w *=      Smeared[ (x+xL)%L ][ (y+dy)%L ][1];
	  for(int dx=xL-1; dx>-1; dx--)  w *= conj(Smeared[ (x+dx)%L ][ (y+yL)%L ][0]);
	  for(int dy=yL-1; dy>-1; dy--)  w *= conj(Smeared[x][ (y+dy)%L ][1]);
	  
	  avWc[xL][yL] += inv_Lsq*w;
	}
    }
    */
  return;
}

//Polyakov loops. x is the spatial dim, y is the temporal dim.
void polyakovLoops(Complex  gauge[L][L][D], Complex polyakov[L/2]){
  
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

      polyakov[dx] += conj(w1)*w2/(1.0*L/2);
      
    }
  }
    
  return;
}


double getTopCharge(Complex  gauge[L][L][D], param_t p){

  Complex w;
  double top = 0.0;  
  Complex Smeared[L][L][D];
  smearLink(Smeared,gauge,p);
  
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++){
      w = Smeared[x][y][0]*Smeared[(x + 1)%L][y][1]*conj(Smeared[x][(y+ 1)%L][0])*conj(Smeared[x][y][1]);
      top += arg(w);  // -pi < arg(w) < pi  Geometric value is an integer.
      //print local def here
      //print arg(w) - [ arg(link1) + arg(link2) + c_arg(link3) + c_arg(link4)]
    }
  return top/TWO_PI;
}


//======================HMC===========================
// momenta conj to theta.  dU/dtheta = i U = I * U


int hmc(Complex gauge[L][L][D], param_t p, int iter) {

  int accept = 0;
  
  double mom[L][L][D];
  Complex gaugeOld[L][L][D];
  Complex phi[L][L][2], chi[L][L][2];
  double H, Hold, rmsq;

  copyLat(gaugeOld,gauge);
  zeroLat(mom); 
  zeroField(phi);
  zeroField(chi);
  H = 0.0;
  Hold = 0.0;

  // init mom[L][L][D]  <mom^2> = 1;
  gaussReal_F(mom); 
  
  if(p.dynamic == true) {    
    //Create gaussian distributed fermion field chi. chi[L][L] E exp(-chi^* chi)
    gaussComplex_F(chi, p);
    //Create pseudo fermion field phi = D chi
    g5Dpsi(phi, chi, gauge, p);    
  }
  
  if (iter >= p.therm) Hold = calcH(mom, gauge, chi, p, false);
  trajectory(mom, gauge, phi, p, iter);
  if (iter >= p.therm) H = calcH(mom, gauge, phi, p, true);
  
  //cout << "Iter = " << iter << " H = " << H << " Hold = " << Hold << " exp(-(H-Hold)) = " << exp(-(H-Hold)) << endl;
  if (iter >= 2*p.therm) {      
    expcount++;
    expdHAve += exp(-(H-Hold));
    dHAve += (H-Hold);
  }
  //  cout << " H = "<< H << endl;   
  //  if(iter%100 == 0) printf("   Hold = %.15g, Final H = %.15g, dH = H - Hold = %.15g, frac = %.10g\n",
  //			    Hold, H, H-Hold, 1.0 - Hold/H);
  // checkRev(mom, phi, chi, Hold, p);     // Check reversibility

  // Metropolis accept/reject step
  if (iter >= p.therm) {    
    if ( drand48() > exp(-(H-Hold)) ) copyLat(gauge,gaugeOld);
    else accept = 1;
  }
  
  return accept;
}

double sum0 = 0.0;
int count0 = 0;

/*=================
  
  H(mom,theta) = \sum_i mom^2(i)/2 + V(theta_i)  
  = \sum_i mom^2(i)/2 + beta*\sum_P real(1 - U_P)
  
  dtheta/dt = dH/dmom = mom
  dmom/dt = - dH/theta = F
  
  F = dmon/dt = dt^2/dtheta  (Newton's Law)
  
  ===============*/

double calcH(double mom[L][L][D], Complex gauge[L][L][D],
	     Complex phi[L][L][2], param_t p, bool postStep) {
  
  double H = 0.0;
  double Hmom = 0.0, Hgauge = 0.0, Hferm = 0.0;
  Complex plaq;
  Complex phitmp[L][L][2];
 
  for(int x=0; x<L;x++)
    for(int y=0; y<L; y++){
      
      plaq = gauge[x][y][0]*gauge[ (x+1)%L ][y][1]*conj(gauge[x][ (y+1)%L ][0])*conj(gauge[x][y][1]);
      Hgauge += p.beta*real(1.0 - plaq);
      
      for(int mu=0; mu<D; mu++){
	Hmom += 0.5 * mom[x][y][mu] * mom[x][y][mu];
      }
    }
  H = Hmom + Hgauge;
     
  if(p.dynamic == true) {

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
    //cout << "After Fermion Force H  = " << H << endl;
    
  }
  //cout << "Hmom = " << Hmom << " Hgauge = " << Hgauge << " Hferm = " << Hferm << endl;  
  return Hmom + Hgauge + Hferm;
}

void trajectory(double mom[L][L][D], Complex gauge[L][L][D],
		Complex phi[L][L][2], param_t p, int iter) {  
  
  double dtau = p.tau/p.nstep;
  double H = 0.0;
  Complex guess[L][L][2];
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	guess[x][y][s] = 0.0;
    
  //gauge force
  double fU[L][L][D];
  zeroLat(fU);
  //fermion fermion
  double fD[L][L][D];
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
    int ARPACK_iter = arpack_solve_double(gauge, p, guess, 0, k, iter);
#endif
    //H = calcH(mom, gauge, phi, p, true);
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

// x+ = 0th dim
// y+ = 1th dim
//
//   y+1<_______                              x-1/y+1<_____________>_x+1/y+1 
//     |       ^				     |       ^      |
//     |       |				     |       #      |
//     |       |				     |       #      v
//     v======>|				     v----->-#<-----|
//     x/y    x+1                             x-1/y       x/y   x+1/y
//     |       |				   
//     ^       v				   
// y-1 |<------|				
//    x/y       x+1
//
// - p.beta d/dtheat[x][y][mu] (1 -   Real[U_P])  = p.betas (I U  - I U^*)/2  =-  p.beta imag[U] 


void forceU(double fU[L][L][D], Complex gauge[L][L][D], param_t p) {
  
  Complex plaq0;
  Complex plaq;
  zeroLat(fU);
  int xp1, xm1, yp1, ym1;
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++) {

      xp1 = (x+1)%L;
      yp1 = (y+1)%L;
      ym1 = (y-1+L)%L;
      xm1 = (x-1+L)%L;
      
      plaq0 = gauge[x][y][0]*gauge[xp1][y][1]*conj(gauge[x][yp1][0])*conj(gauge[x][y][1]);
      fU[x][y][0] += p.beta*imag(plaq0);
      
      plaq =  gauge[x][ym1][0]*gauge[xp1][ym1][1]*conj(gauge[x][y][0])*conj(gauge[x][ym1][1]);
      fU[x][y][0] -= p.beta*imag(plaq);

      plaq =  gauge[x][y][1]*conj(gauge[xm1][yp1][0])*conj(gauge[xm1][y][1])*gauge[xm1][y][0];
      fU[x][y][1] += p.beta*imag(plaq);

      //This plaquette was aleady computed. We want the conjugate.
      fU[x][y][1] -= p.beta*imag(plaq0);
      
      //cout << "(" << x << "," << y <<") Forces = "<< fU[x][y][0] << " " << fU[x][y][1] << endl;
    }
}

// let dD \equiv (d/dtheta D)
//
// d/dtheta (phi^* (DD^dag)^-1 phi) = -((DD^dag)^1 phi)^dag ([dD]*D^dag + D*[dD^dag]) ((DD^dag)^-1 phi)
//
// *****  Should optimize this to operate only on EVEN sites. ****

void forceD(double fD[L][L][D], Complex gauge[L][L][D], Complex phi[L][L][2], param_t p){

  if(p.dynamic == true) {

    zeroLat(fD);
    
    //phip = (DD^dag)^-1 * phi
    Complex phip[L][L][2];
    zeroField(phip);

    //Ainv_psi inverts using the DdagD (g5Dg5D) operator, returns
    // phip = (D^-1 * Ddag^-1) phi = (D^-1 * g5 * D^-1 g5) phi.
    Complex guess[L][L][2]; //Initial guess to CG
    zeroField(guess);
    Ainv_psi(phip, phi, guess, gauge, p);
    
    //g5Dphi = g5D * phip
    Complex g5Dphi[L][L][2];
    zeroField(g5Dphi);
    g5Dpsi(g5Dphi, phip, gauge, p);
      
    //URBACH calling sequence.
    // 1) Generate momenta in gp array, add to Hold.
    // 2) Generate random fermion field \xi (\xi = g_X in URBACH), add to Hold.
    // 3) Compute (gamma5 D) \xi = \phi (\phi = g_fermion in URBACH)
    // 4) Enter Leapfrog
    //     4.1) Initial half step.
    //          update_momentum:
    //          i) Run CG with \phi as input, place result in \xi, g5Dg5D as operator
    //          ii) Apply g5D to result.
    //          iii) Compute momentum update
    //     4.2) Run through leapfrog sequence.
    //     4.3) Final half step.
    // 5) Compute final energy.
    // 6) Accept/Reject
    
    int xp1, xm1, yp1, ym1;
    double r = 1.0;
    for(int x=0; x<L;x++)
      for(int y=0; y<L;y++) {

	xp1 = (x+1)%L;
	yp1 = (y+1)%L;
	ym1 = (y-1+L)%L;
	xm1 = (x-1+L)%L;
	
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
void update_mom(double fU[L][L][D], double fD[L][L][D], double mom[L][L][D], double dtau){

  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int mu=0; mu<D; mu++)
	mom[x][y][mu] -= (fU[x][y][mu] - fD[x][y][mu])*dtau;
}

//U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
void update_gauge(Complex gauge[L][L][D], double mom[L][L][D], double dtau){
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int mu=0; mu<D; mu++)
	gauge[x][y][mu] *= polar(1.0, mom[x][y][mu] * dtau);
}


// Do a parity flip based on eqn 4 of arXiv:1203.2560v2
void DoPFlip(Complex gauge[L][L][D], param_t p) {
  Complex parity_gauge[L][L][D];
  for (int x = 0; x < L; x++) {
    for (int y = 0; y < L; y++) {
      parity_gauge[(-x-1+2*L)%L][y][0] = conj(gauge[x][y][0]);
      parity_gauge[(-x+2*L)%L][y][1] = gauge[x][y][1];
    }
  }
  copyLat(gauge, parity_gauge);
}



