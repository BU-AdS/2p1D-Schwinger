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
#include "measurementHelpers3D.h"
#include "fermionHelpers.h"
#include "dOpHelpers.h"

#ifdef USE_ARPACK
#include "arpack_interface_staggered.h"
#endif

//HMC
//----------------------------------------------------------------------------
int hmc(Complex gauge[LX][LY][LZ][D], param_t p, int iter);
void trajectory(double mom[LX][LY][LZ][D], Complex gauge[LX][LY][LZ][D],
		Complex phi[LX][LY], param_t p);
void forceU(double fU[LX][LY][LZ][D], const Complex gauge[LX][LY][LZ][D], param_t p);
void forceD(double fD[LX][LY][2], const Complex gauge[LX][LY][2],
	    Complex phi[LX][LY], param_t p);
void update_mom(double fU[LX][LY][LZ][D], double fD[LX][LY][2],
		double mom[LX][LY][LZ][D], param_t p, double dtau);
void update_gauge(Complex gauge[LX][LY][LZ][D], double mom[LX][LY][LZ][D],
		  param_t p, double dtau);
//----------------------------------------------------------------------------

//Global variables.
int hmccount = 0;
double expdHAve = 0.0;
double dHAve = 0.0;

int main(int argc, char **argv) {

  param_t p;
  
  p.beta  = atof(argv[1]);
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
  
  Complex gauge[LX][LY][LZ][D];  
  Complex gauge2D[LX][LY][2];
  
  double sigma[LX/2];
  Complex pLoops[LX/2];
  Complex wLoops[LX/2][LY/2];
  
  for(int x=0; x<LX/2; x++) {
    sigma[x] = 0.0;
    pLoops[x] = 0.0;
    for(int y=0; y<LY/2; y++) {
      wLoops[x][y] = Complex(0.0,0.0);
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
	measPolyakovLoops(gauge2D, pLoops);	
	
	name = "data/polyakov/polyakov_Lz" + to_string(z);
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
	
	name = "data/polyakov/polyakov_ratios_Lz" + to_string(z);
	constructName(name, p);
	name += ".dat";
	sprintf(fname, "%s", name.c_str());
	fp = fopen(fname, "a");
	fprintf(fp, "%d ", iter + 1);
	for(int size=1; size<LX/2-1; size++)
	  fprintf(fp, "%.16e ",
		  real(pLoops[size+1])/real(pLoops[size]));
	fprintf(fp, "\n");
	fclose(fp);


	//Creutz Ratios
	zeroWL(wLoops);
	measWilsonLoops(gauge2D, wLoops, p);
	
	for(int size=1; size<LX/2; size++) {
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
	for(int size=2; size<LX/2; size++)
	  fprintf(fp, "%.16e ", sigma[size]);
	fprintf(fp, "\n");
	fclose(fp);
	
	for(int sizex=2; sizex<LX/2; sizex++)
	  for(int sizey=sizex-1; (sizey < LY/2 && sizey <= sizex+1); sizey++) {
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

//======================HMC===========================
// momenta conj to theta.  dU/dtheta = i U = I * U
int hmc(Complex gauge[LX][LY][LZ][D], param_t p, int iter) {

  int accept = 0;
  
  double mom[LX][LY][LZ][D];
  double mom2D[LX][LY][2];
  Complex gaugeOld[LX][LY][LZ][D];
  Complex gauge2D[LX][LY][2];
  Complex phi[LX][LY], chi[LX][LY];
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
    Dpsi(phi, chi, gauge2D, p);
    
    for(int x=0; x<LX; x++)  //Masks out odd sites.
      for(int y=0; y<LY; y++)
	if((x+y)%2 == 1) phi[x][y] = 0.0;
  }
  
  if (iter >= p.therm) Hold = measAction(mom, gauge, phi, p, false);    
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
		Complex phi[LX][LY], param_t p) {

  double dtau = p.tau/p.nstep;
  double H = 0.0;
  Complex guess[LX][LY];
#ifdef USE_ARPACK
  gaussComplex_F(guess, p);
#endif
  
  //gauge force
  double fU[LX][LY][LZ][D];
  zeroLat(fU);
  //fermion fermion
  double fD[LX][LY][2];
  zeroField(fD);
  
  Complex gauge2D[LX][LY][2];

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

// - p.beta d/dtheat[x][y][mu] (1-Real[U_P]) = p.betas (IU-IU^*)/2 = -p.beta*imag[U] 
// All Wilson loops are computed clockwise.
void forceU(double fU[LX][LY][LZ][D], const Complex gauge[LX][LY][LZ][D], param_t p) {

  zeroLat(fU);

  Complex plaq, plaq0;
  double beta = p.beta;
  double betaz= p.betaz;
  
  int xp1, xm1, yp1, ym1, zp1, zm1;
  
  for(int x=0; x<L; x++) {
    xp1 = (x+1)%LX;
    xm1 = (x-1+LX)%LX;
    for(int y=0; y<L; y++) {
      yp1 = (y+1)%LY;
      ym1 = (y-1+LY)%LY;
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

void forceD(double fD[LX][LY][2], const Complex gauge[LX][LY][2], Complex chi[LX][LY],
	    param_t p) {
  
  Complex chitmp[LX][LY];
  Complex Dchitmp[LX][LY];
  zeroField(chitmp);
  zeroField(Dchitmp);
  zeroLat(fD);
  
  Ainv_psi(chitmp, chi, chitmp, gauge, p); // note chitmp = 0 for ODD
  Dpsi(Dchitmp, chitmp, gauge, p); // restrict to Dslash, m = 0

  for(int x=0; x<LX;x++)
    for(int y=0;y<LY; y++){
      if((x+y)%2 == 1)  chitmp[x][y] = Complex(0.0,0.0);
      if((x+y)%2 == 0) Dchitmp[x][y] = Complex(0.0,0.0);
    }
    
  double eta1;
  for(int x=0; x<L;x++)
    for(int y=0; y<L;y++) {
      eta1 =(1.0 - 2.0*(x%2));
	
      if((x+y+1)%2 == 0){ 
	fD[x][y][0] += 2.0*imag(conj(Dchitmp[x][y]) *
				gauge[x][y][0] * chitmp[(x+1)%LX][y]);
      }
      else {
	fD[x][y][0] += 2.0*imag(conj(Dchitmp[(x+1)%LX][y]) *
				conj(gauge[x][y][0]) * chitmp[x][y]);
      };
	
      if((x+y+1)%2 == 0){    
	fD[x][y][1] += 2.0*eta1*imag(conj(Dchitmp[x][y]) *
				     gauge[x][y][1] * chitmp[x][(y+1)%LY]);
      }
      else {
	fD[x][y][1] += 2.0*eta1*imag(conj(Dchitmp[x][(y+1)%LY]) *
				     conj(gauge[x][y][1]) * chitmp[x][y]);
      }
    }	    
}

//P_{k+1/2} = P_{k-1/2} - dtau * (fU + fD)
void update_mom(double fU[LX][LY][LZ][D], double fD[LX][LY][2], double mom[LX][LY][LZ][D],
		param_t p, double dtau){
  
  //Always update from the 2D gauge fields
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++) 
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<2; mu++) {
	  mom[x][y][z][mu] -= fU[x][y][z][mu]*dtau;	  
	}
  
  //Update from the fermion field if dynamic
  if(p.dynamic == true) {
    for(int x=0; x<LX; x++)
      for(int y=0; y<LY; y++)
	for(int mu=0; mu<2; mu++) {
	  mom[x][y][(LZ-1)/2][mu] += fD[x][y][mu]*dtau;
	}  
  }
}
      

//U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
void update_gauge(Complex gauge[LX][LY][LZ][D], double mom[LX][LY][LZ][D],
		  param_t p, double dtau){

  //Always update from the 2D gauge fields
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<2; mu++) {
	  gauge[x][y][z][mu] *= polar(1.0, mom[x][y][z][mu] * dtau);
	}
  
  //Update from the extra dimension if not z locked.
  if(p.lockedZ == false) {
    for(int x=0; x<LX; x++)
      for(int y=0; y<LY; y++) 
	for(int z=0; z<LZ; z++)
	  for(int mu=2; mu<D; mu++) {
	    gauge[x][y][z][mu] *= polar(1.0, mom[x][y][z][mu] * dtau);
	  }
  }
}
