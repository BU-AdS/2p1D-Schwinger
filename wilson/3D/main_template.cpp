#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>

using namespace std;

#define LX __LX__
#define LY __LY__
#define LZ __LZ__
#define D 3
#define NEV 24
#define NKR 32
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
#include "measurementHelpers3D.h"
#include "fermionHelpers.h"
#include "dOpHelpers.h"
#include "inverters.h"
#include "hmcHelpers3D.h"

#ifdef USE_ARPACK
#include "arpack_interface_wilson.h"
#endif

////Dimension dependent HMC functions defined in main file
//----------------------------------------------------------------------------
void trajectory(double mom[LX][LY][LZ][D], Complex gauge[LX][LY][LZ][D],
		Complex phi[LX][LY][2], param_t p);
int hmc(Complex gauge[LX][LY][LZ][D], param_t p, int iter);
void forceU(double fU[LX][LY][LZ][D], const Complex gauge[LX][LY][LZ][D], param_t p);
void update_mom(double fU[LX][LY][LZ][D], double fD[LX][LY][2],
		double mom[LX][LY][LZ][D], double dtau, param_t p);
void update_gauge(Complex gauge[LX][LY][LZ][D], double mom[LX][LY][LZ][D], double dtau, param_t p);
//----------------------------------------------------------------------------

//Global variables.
int hmccount = 0;
double expdHAve = 0.0;
double dHAve = 0.0;

int main(int argc, char **argv) {

  param_t p;
  
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
  //Pseudo RNG seed
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
  p.nKr = NKR;
  p.nEv = NEV;
  p.arpackTol = atof(argv[18]);
  p.arpackMaxiter = atoi(argv[19]);
  p.polyACC = atoi(argv[20]);
  p.amax = atof(argv[21]);
  p.amin = atof(argv[22]);
  p.n_poly = atoi(argv[23]);

  //Measurements
  if(atoi(argv[24]) == 0) p.measPL = false;
  else p.measPL = true;

  if(atoi(argv[25]) == 0) p.measWL = false;
  else p.measWL = true;
  
  if(atoi(argv[26]) == 0) p.measPC = false;
  else p.measPC = true;
  
  if(atoi(argv[27]) == 0) p.measVT = false;
  else p.measVT = true;  
  
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
  
  Complex gauge[LX][LY][LZ][D];
  zeroLat(gauge);  
  Complex gauge2D[LX][LY][2];
  zeroLat(gauge2D);
  
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
    readGaugeLattice(gauge, name);
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

  //Begin thermalised trajectories
  //---------------------------------------------------------------------
  for(iter=iter_offset; iter<p.iterHMC + iter_offset; iter++){

    //Perform HMC step
    accept = hmc(gauge, p, iter);

    //HMC acceptance rate
    accepted += accept;
          
    //Measure the topological charge if trajectory is accepted
    //---------------------------------------------------------------------
    if(accept == 1) {

      for(int z=0; z<LZ; z++) {
	
	extractLatSlice(gauge, gauge2D, z);
	
	top = measTopCharge(gauge2D, p);
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
    }
    
    //Perform Measurements
    //---------------------------------------------------------------------
    if( (iter+1)%p.skip == 0) {
      
      count++; //Number of measurements taken

      //Checkpoint the gauge field?
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
            
      //Dump simulation data to stdout
      double time = time0 + clock();
      cout << fixed << setprecision(16) << iter+1 << " "; //Iteration
      cout << time/(CLOCKS_PER_SEC) << " ";               //Time
      cout << plaqSum[(LZ-1)/2]/count << " ";             //Central Plaquette Action
      cout << (double)top_stuck[(LZ-1)/2]/(accepted) << " "; //P(stuck)
      cout << expdHAve/hmccount << " ";                   //Average exp(-dH)
      cout << dHAve/hmccount << " ";                      //Average dH
      cout << (double)accepted/(count*p.skip) << " ";     //Acceptance
      for(int z=0; z<LZ; z++) cout << top_int[z] << " ";  //Top Charge
      cout << endl;
      
      for(int z=0; z<LZ; z++) {
      
	//Dump simulation data to file
	name = "data/data/data_Lz" + to_string(z);
	constructName(name, p);
	name += ".dat";	
	sprintf(fname, "%s", name.c_str());	
	fp = fopen(fname, "a");	
	fprintf(fp, "%d %.16e %.16e %.16e %.16e %.16e %.16e %d\n",
		iter+1,		
		time/CLOCKS_PER_SEC,
		plaqSum[z]/count,
		(double)top_stuck[z]/(accepted),
		expdHAve/hmccount,
		dHAve/hmccount,
		(double)accepted/(count*p.skip),
		top_int[z]);
	fclose(fp);
      }	

      //Physical observables
      //-------------------------------------------------------------
      //All observables will be measured on the central 2D slice
      int cz = (LZ-1)/2;
      extractLatSlice(gauge, gauge2D, cz);
      
      //Gauge observables
      if(p.measPL || p.measWL) measWilsonLoops(gauge2D, plaq[cz], iter, p);
            
      //Pion Correlation
      if(p.measPC) measPionCorrelation(gauge2D, top_old[cz], iter, p);
      
      //Vacuum Trace
      if(p.measVT) measVacuumTrace(gauge2D, top_old[cz], iter, p);
      //-------------------------------------------------------------
    }
  }   
  return 0;
}

// HMC Routines
//---------------------------------------------------------------------
int hmc(Complex gauge[LX][LY][LZ][D], param_t p, int iter) {

  int accept = 0;
  
  double mom[LX][LY][LZ][D];
  double mom2D[LX][LY][2];
  Complex gaugeOld[LX][LY][LZ][D];
  Complex gauge2D[LX][LY][2];
  Complex phi[LX][LY][2], chi[LX][LY][2];
  double H, Hold;
  
  copyLat(gaugeOld, gauge);
  zeroLat(mom); 
  zeroField(phi);
  zeroField(chi);
  H = 0.0;
  Hold = 0.0;

  // init mom[LX][LY][LZ][D]  <mom^2> = 1;
  gaussReal_F(mom); 
  
  if(p.dynamic == true) {    
    //Create gaussian distributed fermion field chi. chi[LX][LY] E exp(-chi^* chi)
    gaussComplex_F(chi, p);
    //Create pseudo fermion field, phi = D * chi
    extractLatSlice(gauge, gauge2D, (LZ-1)/2);
    g3Dpsi(phi, chi, gauge2D, p);    
  }
  
  if (iter >= p.therm) Hold = measAction(mom, gauge, chi, p, false);    
  trajectory(mom, gauge, phi, p);
  if (iter >= p.therm) H = measAction(mom, gauge, phi, p, true);
  
  if (iter >= 2*p.therm) {      
    hmccount++;
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

void trajectory(double mom[LX][LY][LZ][D], Complex gauge[LX][LY][LZ][D],
		Complex phi[LX][LY][2], param_t p) {  

  double dtau = p.tau/p.nstep;
  double H = 0.0;
  Complex guess[LX][LY][2];
  Complex gauge2D[LX][LY][2];
  extractLatSlice(gauge, gauge2D, (LZ-1)/2);
#ifdef USE_ARPACK
  //deflate using phi as source
  //Deflation eigenvectors
  //Complex defl_evecs[NEV][LX][LY][2];
  //Deflation eigenvalues
  //Complex defl_evals[NEV];
  //copyField(guess, phi);
  //arpack_solve(gauge2D, defl_evecs, defl_evals, 0, 0, p);
  //deflate(guess, phi, defl_evecs, defl_evals, p);
  zeroField(guess);
#else
  zeroField(guess);
#endif
  
  //gauge force
  double fU[LX][LY][LZ][D];
  zeroLat(fU);
  //fermion fermion
  double fD[LX][LY][2];
  zeroField(fD);
  
  //Initial half step.
  //P_{1/2} = P_0 - dtau/2 * (fU + fD)
  forceU(fU, gauge, p);
  forceD(fD, gauge2D, phi, guess, p);
  update_mom(fU, fD, mom, 0.5*dtau, p);  
  
  for(int k=1; k<p.nstep; k++) {
    
    //U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
    update_gauge(gauge, mom, dtau, p);
    
    //P_{k+1/2} = P_{k-1/2} - dtau * (fU - fD)
    forceU(fU, gauge, p);
    extractLatSlice(gauge, gauge2D, (LZ-1)/2);
    forceD(fD, gauge2D, phi, guess, p);
    update_mom(fU, fD, mom, dtau, p);  
  }
  
  //Final half step.
  //U_{n} = exp(i dtau P_{n-1/2}) * U_{n-1}
  update_gauge(gauge, mom, dtau, p);
  
  //P_{n} = P_{n-1/2} - dtau/2 * (fU + fD)
  forceU(fU, gauge, p);
  extractLatSlice(gauge, gauge2D, (LZ-1)/2);
  forceD(fD, gauge2D, phi, guess, p);
  update_mom(fU, fD, mom, 0.5*dtau, p);
  
  //trajectory complete
}
//---------------------------------------------------------------------------------
