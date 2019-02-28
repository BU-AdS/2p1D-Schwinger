#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>

using namespace std;

#define L 8
#define LZ 5
#define D 3
#define PI 3.141592653589793
#define TWO_PI 6.283185307179586

typedef complex<double> Complex;
#define I Complex(0,1)
#define cUnit Complex(1,0)

#include "utils.h"
#include "utils3D.h"
#include "latHelpers.h"
#include "latHelpers3D.h"
#include "measurementHelpers.h"
#include "wFermionHelpers.h"
#include "dOpHelpers.h"

#ifdef USE_ARPACK
#include "arpack_interface.h"
#endif

//HMC
//----------------------------------------------------------------------------
int hmc(Complex gauge[L][L][LZ][D], param_t p, int iter);
double calcH(double mom[L][L][LZ][D], Complex gauge[L][L][LZ][D],
	     Complex phi[L][L][2], param_t p, bool postStep);
double calcH2D(double mom[L][L][2], Complex gauge[L][L][2],
	       Complex phi[L][L][2], param_t p, bool postStep);
void trajectory(double mom[L][L][LZ][D], Complex gauge[L][L][LZ][D],
		Complex phi[L][L][2], param_t p);
void forceU(double fU[L][L][LZ][D], const Complex gauge[L][L][LZ][D], param_t p);
void forceD(double fD[L][L][2], const Complex gauge[L][L][2],
	    Complex phi[L][L][2], param_t p);
void update_mom(double fU[L][L][LZ][D], double fD[L][L][2],
		double mom[L][L][LZ][D], param_t p, double dtau);
void update_gauge(Complex gauge[L][L][LZ][D], double mom[L][L][LZ][D],
		  param_t p, double dtau);
//----------------------------------------------------------------------------

//Global variables.
int hmccount = 0;
double expdHAve = 0.0;
double dHAve = 0.0;

int main(int argc, char **argv) {

  param_t p;
  p.Latsize = L;
  
  p.beta = atof(argv[1]);
  p.betaz = atof(argv[2]);
  p.iterHMC = atoi(argv[3]);
  p.therm = atoi(argv[4]);
  p.skip = atoi(argv[5]);
  p.chkpt = atoi(argv[6]);
  p.checkpointStart = atoi(argv[7]);  
  p.nstep = atoi(argv[8]);
  p.tau = atof(argv[9]);
  
  p.smearIter = atoi(argv[10]);
  p.alpha = atof(argv[11]);  
  long iseed = (long)atoi(argv[12]);
  srand48(iseed);

  
  if(atoi(argv[13]) == 0) 
    p.dynamic = false;
  else
    p.dynamic = true;

  if(atoi(argv[14]) != 0)
    p.lockedZ = true;
  else
    p.lockedZ = true;
  
  p.m = atof(argv[15]);
  
  p.maxIterCG = atoi(argv[16]);
  p.eps = atof(argv[17]);
  
  //Arpack params
  p.nKv = atoi(argv[18]);
  p.nEv = atoi(argv[19]);
  p.arpackTol = atof(argv[20]);
  p.arpackMaxiter = atoi(argv[21]);
  p.polyACC = atoi(argv[22]);
  p.amax = atof(argv[23]);
  p.amin = atof(argv[24]);
  p.n_poly = atoi(argv[25]);
  
  //Topology
  double top = 0.0;
  int top_int[LZ];
  int top_old[LZ];
  int top_stuck[LZ];

  int histL = 101;
  int histQ[LZ][histL];
  double plaqSum[LZ];
  double plaq[LZ];
  int index[LZ];
  for(int z=0; z<LZ; z++) {
    plaqSum[z] = 0.0;
    plaq[z] = 0.0;
    index[z] = 0;
    top_stuck[z] = 0;
    top_old[z] = 0;
    top_int[z] = 0;
    for(int i = 0; i < histL; i++) histQ[z][i] = 0;
  }
  
  Complex gauge[L][L][LZ][D];  
  Complex gauge2D[L][L][2];
  
  double sigma[L/2];
  Complex pLoops[L/2];
  Complex wLoops[L/2][L/2];
  
  for(int i=0; i<L/2; i++) {
    sigma[i] = 0.0;
    pLoops[i] = 0.0;
    for(int j=0; j<L/2; j++) {
      wLoops[i][j] = Complex(0.0,0.0);
    }
  }
  
  int count = 0;
  string name;
  fstream outPutFile;

  int accept;
  int accepted = 0;
  char fname[256];
  FILE *fp;

  gaussStart(gauge,p);  // hot start

  //Start simulation
  double time0 = -((double)clock());
  int iter_offset = 0;
  int iter;

  if(p.checkpointStart > 0) {

    //Read in gauge field if requested
    //---------------------------------------------------------------------
    name = "gauge/gauge";
    constructName(name, p);
    name += "_traj" + to_string(p.checkpointStart) + ".dat";	
    readGaugeLattice(gauge, name);
    iter_offset = p.checkpointStart;    
  } else {

    //Thermalise from random start
    //---------------------------------------------------------------------
    printParams(p);  
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

    for(int z=0; z<LZ; z++) {

      //All observables will be measured on the extracted 2D slice 
      extractLatSlice(gauge, gauge2D, z);
      
      //Topology
      top = getTopCharge(gauge2D, p);
      top_int[z] = round(top);
      name = "data/top/top_charge_Lz" + to_string(z);
      constructName(name, p);
      name += ".dat";
      sprintf(fname, "%s", name.c_str());
      fp = fopen(fname, "a");
      fprintf(fp, "%d %d\n", iter, top_int[z]);
      fclose(fp);

      index[z] = top_int[z] + (histL-1)/2;
      histQ[z][index[z]]++;
      if(top_old[z] == top_int[z]) top_stuck[z]++;
      top_old[z] = top_int[z];
      
    }

    if( (iter+1)%p.skip == 0) {
      
      count++; //Number of measurements taken

      //Checkpoint the gauge field
      if( (iter+1)%p.chkpt == 0) {	  
	name = "gauge/gauge";
	constructName(name, p);
	name += "_traj" + to_string(iter+1) + ".dat";
	writeGaugeLattice(gauge, name);
      }
      
      //Plaquette actions
      for(int z=0; z<LZ; z++) {
	plaq[z] = measPlaq(gauge, z);
	plaqSum[z] += plaq[z];
      }
      
      double time = time0 + clock();
      
      //Info dumped to stdout
      cout << fixed << setprecision(16) << iter+1 << " "; //Iteration
      cout << time/(CLOCKS_PER_SEC) << " ";               //Time
      cout << plaqSum[(LZ-1)/2]/count << " ";             //Central Plaquette Action
      cout << (double)top_stuck[(LZ-1)/2]/(count*p.skip) << " "; //P(stuck)
      cout << expdHAve/hmccount << " ";                   //Average exp(-dH)
      cout << dHAve/hmccount << " ";                      //Average dH
      cout << (double)accepted/(count*p.skip) << " ";     //Acceptance
      for(int z=0; z<LZ; z++) cout << top_int[z] << " ";  //Top Charge
      cout << endl;
      
      for(int z=0; z<LZ; z++) {
      
	//Dump to file
	name = "data/data/data_Lz" + to_string(z);//I cannot make bricks without clay!
	constructName(name, p);
	name += ".dat";	
	sprintf(fname, "%s", name.c_str());	
	fp = fopen(fname, "a");	
	fprintf(fp, "%d %.16e %.16e %.16e %.16e %.16e\n",
		iter+1,		
		time/CLOCKS_PER_SEC,
		plaqSum[z]/count,
		(double)top_stuck[z]/(count*p.skip),
		expdHAve/hmccount,
		(double)accepted/(count*p.skip));
	fclose(fp);
      }
      
      //Compute gauge observables
      for(int z=0; z<LZ; z++) {

	//All observables will be measured on the extracted 2D slice 
	extractLatSlice(gauge, gauge2D, z);
	
	//Measure Polyakov loops
	for(int i=0; i<L/2; i++) pLoops[i] = 0.0;
	calcPolyakovLoops(gauge2D, pLoops);	
	
	name = "data/polyakov/polyakov_Lz" + to_string(z);
	constructName(name, p);
	name += ".dat";
	sprintf(fname, "%s", name.c_str());
	fp = fopen(fname, "a");
	fprintf(fp, "%d ", iter + 1);
	for(int size=1; size<L/2; size++)
	  fprintf(fp, "%.16e %.16e ",
		  real(pLoops[size]),
		  imag(pLoops[size]) );
	fprintf(fp, "\n");
	fclose(fp);
	
	name = "data/polyakov/polyakov_ratios_Lz" + to_string(z);
	constructName(name, p);
	name += ".dat";
	sprintf(fname, "%s", name.c_str());
	fp = fopen(fname, "a");
	fprintf(fp, "%d ", iter + 1);
	for(int size=1; size<L/2-1; size++)
	  fprintf(fp, "%.16e ",
		  real(pLoops[size+1])/real(pLoops[size]));
	fprintf(fp, "\n");
	fclose(fp);


	//Creutz Ratios
	zeroWL(wLoops);
	calcWilsonLoops(gauge2D, wLoops, p);
	
	for(int size=1; size<L/2; size++) {
	  sigma[size]  = - log(abs((real(wLoops[size][size])/real(wLoops[size-1][size]))* 
				   (real(wLoops[size-1][size-1])/real(wLoops[size][size-1]))));
	  
	  sigma[size] += - log(abs((real(wLoops[size][size])/real(wLoops[size][size-1]))* 
				   (real(wLoops[size-1][size-1])/real(wLoops[size-1][size]))));
	  
	  sigma[size] *= 0.5;
	}
	
	name = "data/creutz/creutz_Lz" + to_string(z);
	constructName(name, p);
	name += ".dat";
	sprintf(fname, "%s", name.c_str());
	fp = fopen(fname, "a");
	fprintf(fp, "%d %.16e ", iter+1, -log(abs(plaq[z])) );
	for(int size=2 ; size < L/2; size++)
	  fprintf(fp, "%.16e ", sigma[size]);
	fprintf(fp, "\n");
	fclose(fp);
	
	for(int sizex=2; sizex<L/2; sizex++)
	  for(int sizey=sizex-1; (sizey < L/2 && sizey <= sizex+1); sizey++) {
	    name = "data/rect/rectWL_Lz" + to_string(z);
	    name += "_" + to_string(sizex) + "_" + to_string(sizey);
	    constructName(name, p);
	    name += ".dat";
	    sprintf(fname, "%s", name.c_str());
	    fp = fopen(fname, "a");
	    fprintf(fp, "%d %.16e %.16e\n", iter+1, real(wLoops[sizex][sizey]), imag(wLoops[sizex][sizey]));	    
	    fclose(fp);
	  }
	
	
	//Topological Historgram
	name = "data/top/top_hist+Lz" + to_string(z);
	constructName(name, p);
	name += ".dat";
	sprintf(fname, "%s", name.c_str());
	fp = fopen(fname, "w");
	for(int i=0; i<histL; i++) fprintf(fp, "%d %d\n", i - (histL-1)/2, histQ[z][i]);
	fclose(fp);
	
      }
    }
  }
  
  return 0;
}


//======================================================================
//   Creutz     exp[ -sigma L^2] exp[ -sigma(L-1)(L-1)]
//   ratio:    ---------------------------------------  = exp[ -sigma]
//              exp[ -sigma (L-1)L] exp[-sigma L(L-1)]
//======================================================================
/*
void areaLaw(Complex gauge[L][L][LZ][D], Complex avWc[L/2][L/2][LZ], param_t p){
  
  Complex w;
  int p1, p2;
  double inv_Lsq = 1.0/(L*L);
  
  //Loop over z planes
  for(int z=0; z<LZ; z++) {

    Complex smeared[L][L][2];
    Complex gauge2D[L][L][2];
    extractLatSlice(gauge, gauge2D, z);
    smearLink(smeared, gauge2D, p);
    
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
	    avWc[Xrect][Yrect][z] += w*inv_Lsq;
	  }
      }
    }
    
  }
  
  return;
}
*/

//======================HMC===========================
// momenta conj to theta.  dU/dtheta = i U = I * U


int hmc(Complex gauge[L][L][LZ][D], param_t p, int iter) {

  int accept = 0;
  
  double mom[L][L][LZ][D];
  double mom2D[L][L][2];
  Complex gaugeOld[L][L][LZ][D];
  Complex gauge2D[L][L][2];
  Complex phi[L][L][2], chi[L][L][2];
  double H, Hold, rmsq;

  copyLat(gaugeOld, gauge);
  zeroLat(mom); 
  zeroField(phi);
  zeroField(chi);
  H = 0.0;
  Hold = 0.0;

  // init mom[L][L][LZ][D]  <mom^2> = 1;
  gaussReal_F(mom); 
  
  if(p.dynamic == true) {    
    //Create gaussian distributed fermion field chi. chi[L][L] E exp(-chi^* chi)
    gaussComplex_F(chi, p);
    //Create pseudo fermion field phi = D chi
    extractLatSlice(gauge, gauge2D, (LZ-1)/2);
    g5Dpsi(phi, chi, gauge2D, p);    
  }
  
  if (iter >= p.therm) {
    //extractLatSlice(gauge, gauge2D, (LZ-1)/2);
    //extractLatSlice(mom, mom2D, (LZ-1)/2);
    //Hold = calcH2D(mom2D, gauge2D, chi, p, false);
    Hold = calcH(mom, gauge, chi, p, false);    
  }
  
  trajectory(mom, gauge, phi, p);
  
  if (iter >= p.therm) {
    //extractLatSlice(gauge, gauge2D, (LZ-1)/2);
    //extractLatSlice(mom, mom2D, (LZ-1)/2);
    //H = calcH2D(mom2D, gauge2D, phi, p, true);
    H = calcH(mom, gauge, phi, p, true);
  }
  //cout << "Iter = " << iter << " H = " << H << " Hold = " << Hold << " exp(-(H-Hold)) = " << exp(-(H-Hold)) << endl;
  if (iter >= 2*p.therm) {      
    hmccount++;
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

double calcH(double mom[L][L][LZ][D], Complex gauge[L][L][LZ][D],
	     Complex phi[L][L][2], param_t p, bool postStep) {
  
  double H = 0.0;
  double Hmom = 0.0, Hgauge = 0.0, Hferm = 0.0;
  Complex plaq;
  Complex phitmp[L][L][2];
  Complex gauge2D[L][L][2];
  double beta  = p.beta;
  double betaz = p.betaz;
  
  for(int x=0; x<L;x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++) {

	//x, y, -x, -y
	plaq = gauge[x][y][z][0]*gauge[ (x+1)%L ][y][z][1]*conj(gauge[x][ (y+1)%L ][z][0])*conj(gauge[x][y][z][1]);
	Hgauge += beta*real(1.0 - plaq);
	
	if(z != LZ-1) {
	  //x, z, -x, -z
	  plaq = gauge[x][y][z][0] * cUnit * conj(gauge[x][y][ (z+1)%LZ ][0]) * cUnit;
	  Hgauge += betaz*real(1.0 - plaq);
	  
	  //y, z, -y, -z
	  plaq = gauge[x][y][z][1] * cUnit * conj(gauge[x][y][ (z+1)%LZ ][1]) * cUnit;
	  Hgauge += betaz*real(1.0 - plaq);	
	}
	
	//No gauge force in z dim.
	for(int mu=0; mu<2; mu++){	  
	  Hmom += 0.5 * mom[x][y][z][mu]*mom[x][y][z][mu];	  
	  // We do not want the spurious extra dim momentum terms on
	  // the last slice.
	  //if(z == LZ-1 && mu == D-1) {
	  //CHECK
	  //Hmom -= mom[x][y][z][mu]*mom[x][y][z][mu]/2.0;
	  //}
	}
      }
  
  if(p.dynamic == true) {

    //cout << "Before Fermion force H = " << H << endl;
    Complex scalar = Complex(0.0,0.0);
    zeroField(phitmp);    
    if(postStep) {
      extractLatSlice(gauge, gauge2D, (LZ-1)/2);
      Ainv_psi(phitmp, phi, phitmp, gauge2D, p);
    }
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

double calcH2D(double mom[L][L][2], Complex gauge[L][L][2],
	       Complex phi[L][L][2], param_t p, bool postStep) {
  
  double H = 0.0;
  double Hmom = 0.0, Hgauge = 0.0, Hferm = 0.0;
  Complex plaq;
  Complex phitmp[L][L][2];
 
  for(int x=0; x<L;x++)
    for(int y=0; y<L; y++){
      
      plaq = gauge[x][y][0]*gauge[ (x+1)%L ][y][1]*conj(gauge[x][ (y+1)%L ][0])*conj(gauge[x][y][1]);
      Hgauge += p.beta*real(1.0 - plaq);
      
      for(int mu=0; mu<2; mu++){
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

void trajectory(double mom[L][L][LZ][D], Complex gauge[L][L][LZ][D],
		Complex phi[L][L][2], param_t p) {  

  double dtau = p.tau/p.nstep;
  double H = 0.0;
  Complex guess[L][L][2];
#ifdef USE_ARPACK
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int s=0; s<2; s++)
	guess[x][y][s] = drand48();
#endif
  
  //gauge force
  double fU[L][L][LZ][D];
  zeroLat(fU);
  //fermion fermion
  double fD[L][L][2];
  zeroField(fD);
  
  Complex gauge2D[L][L][2];

  //Initial half step.
  //P_{1/2} = P_0 - dtau/2 * (fU + fD)
  forceU(fU, gauge, p);
  extractLatSlice(gauge, gauge2D, (LZ-1)/2);
  forceD(fD, gauge2D, phi, p);
  update_mom(fU, fD, mom, p, 0.5*dtau);  
  
  for(int k=1; k<p.nstep; k++) {
    
    //U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
    update_gauge(gauge, mom, p, dtau);
    
    //P_{k+1/2} = P_{k-1/2} - dtau * (fU - fD)
    forceU(fU, gauge, p);
    extractLatSlice(gauge, gauge2D, (LZ-1)/2);
    forceD(fD, gauge2D, phi, p);
    update_mom(fU, fD, mom, p, dtau);  
       
#ifdef USE_ARPACK
    int ARPACK_iter = arpack_solve_double(gauge2D, p, guess, 1, k, 0);
#endif
    //H = calcH(mom, gauge, phi, p, true);
  }
  
  //Final half step.
  //U_{n} = exp(i dtau P_{n-1/2}) * U_{n-1}
  update_gauge(gauge, mom, p, dtau);
  
  //P_{n} = P_{n-1/2} - dtau/2 * (fU + fD)
  forceU(fU, gauge, p);
  extractLatSlice(gauge, gauge2D, (LZ-1)/2);
  forceD(fD, gauge2D, phi, p);
  update_mom(fU, fD, mom, p, 0.5*dtau);
  
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
// - p.beta d/dtheat[x][y][mu] (1-Real[U_P]) = p.betas (IU-IU^*)/2 = -p.beta*imag[U] 
// All Wilson loops are computed clockwise.
void forceU(double fU[L][L][LZ][D], const Complex gauge[L][L][LZ][D], param_t p) {

  zeroLat(fU);

  Complex plaq, plaq0;
  double beta = p.beta;
  double betaz= p.betaz;
  
  int xp1, xm1, yp1, ym1, zp1, zm1;
  
  for(int x=0; x<L; x++) {
    xp1 = (x+1)%L;
    xm1 = (x-1+L)%L;
    for(int y=0; y<L; y++) {
      yp1 = (y+1)%L;
      ym1 = (y-1+L)%L;
      for(int z=0; z<LZ; z++) {
	zp1 = (z+1)%LZ;
	zm1 = (z-1+LZ)%LZ;
	
	//X dir
	//-------
	// +x, +y, -x, -x
	plaq0 = gauge[x][y][z][0] * gauge[xp1][y][z][1] * conj(gauge[x][yp1][z][0]) * conj(gauge[x][y][z][1]);
	fU[x][y][z][0] += beta*imag(plaq0);

	// -y, +x, +y, -x
       	plaq = conj(gauge[x][ym1][z][1])*gauge[x][ym1][z][0] * gauge[xp1][ym1][z][1]*conj(gauge[x][y][z][0]);
	fU[x][y][z][0] -= beta*imag(plaq);
      
	if(z != LZ-1) {
	  // +x, +z, -x, -z
	  plaq = gauge[x][y][z][0] * cUnit * conj(gauge[x][y][zp1][0]) * cUnit;
	  fU[x][y][z][0] += betaz*imag(plaq);
	}
	  
	if(z != 0) {
	  // -z, +x, +z, -x
	  plaq = cUnit * gauge[x][y][zm1][0] * cUnit * conj(gauge[x][y][z][0]);
	  fU[x][y][z][0] -= betaz*imag(plaq);
	}
	
	//Y dir
	//------
	// +y, -x, -y, +x
	plaq = gauge[x][y][z][1] * conj(gauge[xm1][yp1][z][0]) * conj(gauge[xm1][y][z][1]) * gauge[xm1][y][z][0];
	fU[x][y][z][1] += beta*imag(plaq);
	
	//This plaquette was aleady computed. We want the conjugate.
	fU[x][y][z][1] -= beta*imag(plaq0);
	
	if(z != LZ-1) {
	  // y, z, -y, -z
	  plaq = gauge[x][y][z][1] * cUnit * conj(gauge[x][y][zp1][1]) * cUnit;
	  fU[x][y][z][1] += betaz*imag(plaq);
	}
	
	if(z != 0) {
	  // -z, +y, +z, -y
	  plaq = cUnit * gauge[x][y][zm1][1] * cUnit * conj(gauge[x][y][z][1]);
	  fU[x][y][z][1] -= betaz*imag(plaq);
	}
	
	//Z dir
	//------
	// Only update the z links if not locked
	if(z != LZ-1 && !p.lockedZ) {
	  /*
	  //z, x, -z, -x
	  plaq = gauge[x][y][z][2]*gauge[x][y][zp1][0]*conj(gauge[xp1][y][z][2])*conj(gauge[x][y][z][0]);
	  fU[x][y][z][2] -= betaz*imag(plaq);
	  
	  //z, -x, -z, x
	  plaq = gauge[x][y][z][2]*conj(gauge[xm1][y][zp1][0])*conj(gauge[xm1][y][z][2])*gauge[xm1][y][z][0];
	  fU[x][y][z][2] -= betaz*imag(plaq);
	  
	  //z, y, -z, -y
	  plaq = gauge[x][y][z][2]*gauge[x][y][zp1][1]*conj(gauge[x][yp1][z][2])*conj(gauge[x][y][z][1]);
	  fU[x][y][z][2] -= betaz*imag(plaq);
	  
	  //z, -y, -z, y
	  plaq = gauge[x][y][z][2]*conj(gauge[x][ym1][zp1][1])*conj(gauge[x][ym1][z][2])*gauge[x][ym1][z][1];
	  fU[x][y][z][2] -= betaz*imag(plaq);
	  */
	}
      }
    }
  }
}

// let dD \equiv (d/dtheta D)
//
// d/dtheta (phi^* (DD^dag)^-1 phi) = -((DD^dag)^1 phi)^dag ([dD]*D^dag + D*[dD^dag]) ((DD^dag)^-1 phi)
//
// Should optimise this to operate only on EVEN sites.

void forceD(double fD[L][L][2], const Complex gauge[L][L][2], Complex phi[L][L][2], param_t p){
  
  if(p.dynamic == true) {

    zeroField(fD);
    
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
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++) {

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
void update_mom(double fU[L][L][LZ][D], double fD[L][L][2], double mom[L][L][LZ][D], param_t p, double dtau){
  
  //Always update from the 2D gauge fields
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++) 
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<2; mu++) {
	  mom[x][y][z][mu] -= fU[x][y][z][mu]*dtau;	  
	}
  
  //Update from the fermion field if dynamic
  if(p.dynamic == true) {
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++)
	for(int mu=0; mu<2; mu++) {
	  mom[x][y][(LZ-1)/2][mu] += fD[x][y][mu]*dtau;
	}  
  }
}
      

//U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
void update_gauge(Complex gauge[L][L][LZ][D], double mom[L][L][LZ][D], param_t p, double dtau){

  //Always update from the 2D gauge fields
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<2; mu++) {
	  gauge[x][y][z][mu] *= polar(1.0, mom[x][y][z][mu] * dtau);
	}
  
  //Update from the extra dimension if not z locked.
  if(p.lockedZ == false) {
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++) 
	for(int z=0; z<LZ; z++)
	  for(int mu=2; mu<D; mu++) {
	    gauge[x][y][z][mu] *= polar(1.0, mom[x][y][z][mu] * dtau);
	  }
  }
}
