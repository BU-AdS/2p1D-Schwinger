/* ================================================================================= 

   Sat Jan 22 13:43:09 PST 2005 Richard Brower  
   U(1) NON-COMPACT  Gauge Field for 2-d Quenched Lattice

   This does Heat Bath on the phase angles for the non-compact action:

   prob = Prod_{Plaquettes} = sqrt{PI/ beta} exp( - (beta/2)* (theta_Plaquette)^2 ]
                         
   The exact value of the average plaquette is

   Compared with exact value at beta = 16:
   - Log[Sqrt[begta/(2 Pi)] NIntegrate[2 Cos[x]*Exp[- beta*x^2/2],{x,0,Infinity}] = 1/32 =  0.03125

   This is not relevant to the real problem but the code was also debugged by
   running with purely Gussian links at beta = 16: 

   sigma = Log[plaq]/4 = 0.0312485 (Simulation) Gausian Links 
   sigma = Log[plaq] = 0.0312285 (Simulation) Gausian Plaqs (Note factor of 4)


   Wed Aug 17 10:02:02 EDT 2016 Richard Brower
   U(1) COMPACT Full HMC for 2-d U(1)

   Optimization:

   Should introduce an ee CG for DdagD 
   Should introduce an D_oe Dslash.

   Check: Dslash, Inverters, Dforce.

   ===================================================================================*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>
#include <omp.h>

using namespace std;
#include "ran2s.h"
#ifdef USE_ARPACK
#include "arpack_interface.h"
#endif

#define L 48
#define D 2
#define PI 3.141592653589793
#define TWO_PI 6.283185307179586

typedef complex<double> Complex;
#define I Complex(0,1)

typedef struct{
  
  int Latsize;

  //HMC
  int nstep;
  double dt;
  int iterHMC;
  int therm;
  int skip;
  int chkpt;
  int checkpointStart;
  int maxIterCG;
  double eps;
  
  //physics
  double beta;
  double m;
  bool quenched;

  //Smearing
  double alpha;
  int smearIter;

  //Arpack params
  int nEv;
  int nKv;
  double arpackTol;
  int arpackMaxiter;
  
} param_t;

void printParams(param_t p);
void printLattice(Complex gauge[L][L][D]);
void constructName(string &name, param_t p);
double measPlaq(Complex  gauge[L][L][D]);
void getPhase(double phase64[L][L][D], const Complex gauge[L][L][D]);
void getCompact(Complex gauge[L][L][D], const double phase[L][L][D]);
void phaseUpdate(double phase[L][L][D],param_t p);
void gaussStart(Complex gauge[L][L][D], param_t p);
void coldStart(Complex gauge[L][L][D],param_t p);
void areaLaw(Complex gauge[L][L][D], Complex avW[L][L], param_t p);
void areaFast(const Complex gauge[L][L][D],  Complex avW[L][L]);
void polyakovLoops(Complex  gauge[L][L][D], Complex polyakov[L/2]);
void writeGaugeLattice(Complex gauge[L][L][D], string name);
void readGaugeLattice(Complex gauge[L][L][D], string name);
double getTopCharge(Complex gauge[L][L][D], param_t p);
void smearLink(Complex gaugeSmeared[L][L][D], Complex gauge[L][L][D], param_t p);
void gaussReal_F(double field[L][L][D]);
void gaussComplex_F(Complex eta[L][L], double w, param_t p);
void MakeInstanton(Complex gauge[L][L][D], int Q, param_t p);
void MakePointVortex(Complex gauge[L][L][D], int Q, int x0, int y0, param_t p);
int hmc(Complex gauge[L][L][D], param_t p, int iter);
double calcH(double mom[L][L][D], Complex gauge[L][L][D],Complex chi[L][L], param_t p);
void trajectory(double mom[L][L][D], Complex gauge[L][L][D],
		Complex chi[L][L], param_t p);
void forceV(double fV[L][L][D], Complex gauge[L][L][D],param_t p);
void forceD(double fV[L][L][D],Complex gauge[L][L][D], Complex chi[L][L], param_t p);
void DoPFlip(Complex gauge[L][L][D], param_t p);

void arpack_solve_double(Complex gauge[L][L][D], param_t p, Complex guess[L][L]);
void TestCG(Complex gauge[L][L][D], param_t p);
void TestForce(Complex gauge[L][L][D], param_t p);

// Dslash + m
void Dpsi(Complex psi2[L][L], Complex  psi1[L][L], Complex gauge[L][L][D], param_t p );

// Dslash + m
void Ddagpsi(Complex psi2[L][L], Complex  psi1[L][L], Complex gauge[L][L][D], param_t p );

// gamma_5
void gamma_5(Complex psi2[L][L], Complex  psi1[L][L], param_t p );

// gamma_5 (Dslash + m) gamma_5 (Dslash + m)
void DdagDpsi(Complex psi2[L][L], Complex  psi1[L][L],  Complex gauge[L][L][D], param_t p );

// Blas tools
double norm2(Complex psi[L][L]);
Complex cDotProduct(Complex psi1[L][L], Complex psi2[L][L]);

int Ainv_psi(Complex psi[L][L], Complex b[L][L], Complex  psi0[L][L], Complex gauge[L][L][D], param_t p);

// Zero lattice vectors.
template<typename T> inline void zeroLat(T v[L][L][D])
{
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      for(int mu=0;mu<D;mu++)
	v[x][y][mu] = 0.0;
}

// Copy lattice vector
template<typename T> inline void copyLat(T v2[L][L][D],T v1[L][L][D])
{
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      for(int mu=0;mu<D;mu++)
	v2[x][y][mu] =  v1[x][y][mu];
}

// Add Equ v2 += v1 lattice vector
template<typename T> inline void addEqLat(T v2[L][L][D],T v1[L][L][D])
{
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      for(int mu=0;mu<D;mu++)
	v2[x][y][mu] +=  v1[x][y][mu];
}

// Zero lattice field.
template<typename T> inline void zeroField(T psi[L][L])
{
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      psi[x][y] = 0.0;
}

// Copy lattice field
template<typename T> inline void copyField(T psi2[L][L],T psi1[L][L])
{
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      psi2[x][y] =  psi1[x][y];
}

// Add Equ v2 += v1 lattice field
template<typename T> inline void addEqField(T psi2[L][L],T psi1[L][L])
{
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      psi2[x][y] +=  psi1[x][y];
}

// Add Equ square of real b dot b lattice field
template<typename T> inline T sqrSumField(T b[L][L])
{
  T square = (T) 0.0;
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      square +=  b[x][y]*b[x][y];
  return square;
}

// Add Equ conj(v2) dot v1 lattice field
template<typename T> inline T dotField(T psi1[L][L], T psi2[L][L])
{
  T scalar = (T) 0.0;
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      scalar +=  conj(psi1[x][y])*psi2[x][y];
  return scalar;
}

void printParams(param_t p) {
  cout << endl;
  cout << "Physics:  Size = "<< L << endl;
  cout << "          Beta = "<< p.beta << endl;
  cout << "          Quenching = " << (p.quenched == true ? "True" : "False") << endl;
  if (!p.quenched) cout << "         Mass = "<< p.m << endl;
  cout << "HMC:      Therm Sweeps = " << p.therm << endl; 
  cout << "          Data Points = " << p.iterHMC - p.therm << endl;
  cout << "          Time Step = " << p.dt << endl;
  cout << "          Trajectory Steps " << p.nstep << endl;
  cout << "          Trajectory Length = " << p.dt*p.nstep << endl;
  cout << "          Trajectory Length = " << p.dt*p.nstep << endl;
  cout << "Smearing: APE iter = " << p.smearIter << endl;
  cout << "          APE alpha = " << p.alpha << endl;
}

void constructName(string &name, param_t p) {
  name += "_L" + to_string(p.Latsize) + "_B" + to_string(p.beta);
  if(p.quenched == false) name += "_M"+ to_string(p.m);
  name += "_dt" + to_string(p.dt) + "_nHMCstep" + to_string(p.nstep);
}

int expcount = 0;
double expdHAve = 0.0;

int main(int argc, char **argv) {

  param_t p;

  p.beta = atof(argv[1]); 
  p.iterHMC = atoi(argv[2]);
  p.therm = atoi(argv[3]);
  p.skip = atoi(argv[4]);
  p.chkpt = atoi(argv[5]);
  p.checkpointStart = atoi(argv[6]);
  if(p.checkpointStart > 0) p.iterHMC += p.checkpointStart;
  
  p.Latsize = L;
  p.nstep = 100;
  p.dt = 0.01;
  p.quenched = true;

  p.maxIterCG = 1000;
  p.m = 0.032;
  p.eps = 1e-6;

  p.smearIter = atoi(argv[7]);
  p.alpha = atof(argv[8]);

  long iseed = (long)atoi(argv[9]);
    
  //Arpack params
  p.nEv = 16;
  p.nKv = 32;
  p.arpackTol = 1e-6;
  p.arpackMaxiter = 1000000;
  
  int procs, threads;
#ifdef USE_OMP
  procs = omp_get_num_procs();
  threads = omp_get_max_threads();
  
  cout << "Number of processors available = " << omp_get_num_procs ( ) << endl;
  cout << "Number of threads =              " << omp_get_max_threads ( ) << endl;
#else
  procs = 1;
  threads = 1;
#endif

  //Topology
  double top = 0.0;
  int top_int = 0;
  int top_old = 0;
  int top_stuck = 0;
  
  Complex gauge[L][L][D];
  double phase[L][L][D], initPhase[L][L][D];
  Complex avW[L][L], avWc[L][L];
  Complex guess[L][L];
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++) guess[x][y] = I;
  
  int histL = 101;
  int histQ[histL];
  for(int i = 0; i < histL; i++) histQ[i] = 0;
  
  double area[L], sigma[L];
  Complex polyakov[L/2];
  Complex polyakovAveTest(0.0,0.0);
  
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
  double elapsed = 0.0;
  
  int accept;
  int accepted = 0;
  char fname[256];
  FILE *fp;

  sran2(iseed);
  printParams(p);  
  gaussStart(gauge,p);  // hot start

  //Start simulation
  double time0 = -((double)clock());
  
  for(int iter=0; iter<p.iterHMC; iter++){
    
    //Read in gauge field if requested
    if(p.checkpointStart > 0 && iter == 0) {
      name = "gauge";
      constructName(name, p);
      name += "_traj" + to_string(p.checkpointStart) + ".dat";	
      readGaugeLattice(gauge,name);
      iter += p.checkpointStart;
    }
    
    //Perform HMC step
    accept = hmc(gauge, p, iter);

    //Measure the plaquette and topological charge
    if(iter+1 > p.therm){
      
      //HMC acceptance
      accepted += accept;
      
      //Topology
      top = getTopCharge(gauge,p);
      top_int = round(top);
      name = "top_charge";
      constructName(name, p);
      name += ".dat";
      sprintf(fname, "%s", name.c_str());
      fp = fopen(fname, "a");
      fprintf(fp, "%d %d\n", iter, top_int);
      fclose(fp);

      index = top_int + (histL-1)/2;
      histQ[index] += 1;
      if(top_old == top_int) top_stuck++;
      top_old = top_int;
      
      if( (iter+1)%p.skip == 0 && (iter+1) > p.therm) {

	count++;
	
	//Plaquette action
	double plaq = measPlaq(gauge);
	plaqSum += plaq;
	
	double time = time0 + clock();
	
	cout << fixed << setprecision(16) << iter+1 << " ";
	cout << time/(threads*CLOCKS_PER_SEC) << " ";
	cout << plaqSum/count << " ";
	cout << (double)top_stuck/(count*p.skip) << " ";
	cout << expdHAve/expcount << " ";
	cout << (double)accepted/(count*p.skip) << endl;
	
	//Dump to file
	name = "data";
	constructName(name, p);
	name += ".dat";	
	sprintf(fname, "%s", name.c_str());	
	fp = fopen(fname, "a");
	
	fprintf(fp, "%d %.16e %.16e %.16e %.16e %.16e\n",
		iter+1,
		time/(threads*CLOCKS_PER_SEC),
		plaqSum/count,
		(double)top_stuck/(count*p.skip),
		expdHAve/expcount,
		(double)accepted/(count*p.skip));
	fclose(fp);
	
	//Checkpoint the gauge field
	if( (iter+1)%p.chkpt == 0) {	  
	  name = "gauge";
	  constructName(name, p);
	  name += "_traj" + to_string(iter+1) + ".dat";
	  writeGaugeLattice(gauge, name);
	}

#ifdef USE_ARPACK
	arpack_solve_double(gauge, p, guess);
	//exit(0);
#endif
	
	for(int a=0; a<L/2; a++) polyakov[a] = 0.0;
	polyakovLoops(gauge, polyakov);
	
	zeroField(avWc);
	areaLaw(gauge, avWc, p);
	
	// Creutz Ratios
	for(int size=1; size<L/2; size++)
	  sigma[size]   = - log(abs( (real(avWc[size][size])/real(avWc[size-1][size]))* 
				     (real(avWc[size-1][size-1])/real(avWc[size][size-1])))) ;
	name = "creutz";
	constructName(name, p);
	name += ".dat";
	sprintf(fname, "%s", name.c_str());
	fp = fopen(fname, "a");
	fprintf(fp, "%d %.16e ", iter+1, -log(abs(plaq)) );
	for(int size=2 ; size < L/2; size++)
	  fprintf(fp, "%.16e ", sigma[size]);
	fprintf(fp, "\n");
	fclose(fp);

	name = "rectWL";
	constructName(name, p);
	name += ".dat";
	sprintf(fname, "%s", name.c_str());
	fp = fopen(fname, "a");
	for(int sizex=2 ; sizex < L/2; sizex++)
	  for(int sizey=sizex-1 ; (sizey < L/2 && sizey <= sizex+1) ; sizey++) {
	    name = "rectWL";
	    name += "_" + to_string(sizex) + "_" + to_string(sizey);
	    constructName(name, p);
	    name += ".dat";
	    sprintf(fname, "%s", name.c_str());
	    fp = fopen(fname, "a");
	    fprintf(fp, "%d %.16e %.16e\n", iter+1, real(avWc[sizex][sizey]), imag(avWc[sizex][sizey]));	    
	    fclose(fp);
	  }
	
	name = "polyakov";
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

	name = "polyakov_ratios";
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
	
	name = "top_hist";
	constructName(name, p);
	name += ".dat";
	sprintf(fname, "%s", name.c_str());
	fp = fopen(fname, "w");
	for(int i=0; i<histL; i++) fprintf(fp, "%d %d\n", i - (histL-1)/2, histQ[i]);
	fclose(fp);
	
      }
    }
  }
  
  return 0;
}

/*===============================================================================
  Heatbath: sqrt{PI/ beta} exp( - (beta/2)* [(theta + staple1)^2 + (theta + staple2)^2]
  
  =  sqrt{PI/beta} exp( -beta * [(theta + 1/2(staple1+ staple2))^2 + const])
  
  x+ = 0the dim
  y+ = 1th dim
  
  
  y+1<_______                              x-1/y+1<_____________>_x+1/y+1 
  |       ^				     |       ^      |
  |       |				     |       #      |
  |       |				     |       #      v
  v======>|				     v----->-#<-----|
  x/y       x+1                             x-1/y       x/y   x+1/y
  |       |				   
  ^       v				   
  y-1  <-------|				
  x/y       x+1
     

  theta = theta_gaussian -  (1/2)staple  where staple = staple1 + staple2

  <Gaussian^2> = 1/2 beta  

  Area Law:  Wilson Loop = exp[ - sigma L^2 ]   sigma = - Log[ <cos(theta_plaq)> ]
  ================================================================================*/ 

void phaseUpdate(double phase[L][L][D],param_t p){
  Complex w = Complex(1.0,0.0);
  double staple; 

  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++) {
      //x or 0 link update       
      staple  = phase[(x + 1)%L][y][1] - phase[x][(y+ 1)%L][0] - phase[x][y][1];
      staple += - phase[(x+1)%L][(y-1+L)%L][1] - phase[x][(y-1+L)%L][0] + phase[x][(y-1+L)%L][1];
      phase[x][y][0] = sqrt(0.5/p.beta)*rang() - staple/2.0;
    };
  
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++) {
      //y or 1 link update
      staple = - phase[(x-1+L)%L][(y+1)%L][0] -phase[(x-1+L)%L][y][1] + phase[(x-1+L)%L][y][0];
      staple +=  phase[x][(y+1)%L][0] - phase[(x+1)%L][y][1] - phase[x][y][0];
      phase[x][y][1] = sqrt(0.5/p.beta)*rang() - staple/2.0;
    };
  
  return;
}  

//================= UTILILTIES ROUTINES =======================

void getPhase(double phase[L][L][D],const Complex gauge[L][L][D]){

  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      for(int mu=0; mu<D; mu++)
	phase[x][y][mu] = arg(gauge[x][y][mu]);
 
  return;
}

void getCompact(Complex gauge[L][L][D], const double phase[L][L][D]){

  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      for(int mu=0; mu<D; mu++)
	gauge[x][y][mu] = polar(1.0,phase[x][y][mu]);

  return;
}

void printLattice(Complex gauge[L][L][D]){

  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++) {
      cout << " (x, y) = " << x << "," << y << " ";
      cout << gauge[x][y][0]<< " " << gauge[x][y][1] << endl;
    }
  return;
}

/*===============================================================================
  Gaussian numbers with p(theta) = sqrt(beta/ 2 PI) exp( - beta* theta^2/2)
  <Gaussian^2> = 1/beta  
  Perimeter Law:  Wilson Loop = exp[ - 4 sigma L ]   sigma = - Log[ <cos(theta)> ]
  ================================================================================*/ 
void gaussStart(Complex gauge[L][L][D],param_t p){

  //#pragma omp parallel for  
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++){
      gauge[x][y][0] = polar(1.0,sqrt(1.0/p.beta)*rang());
      gauge[x][y][1] = polar(1.0,sqrt(1.0/p.beta)*rang());
    }
  return;
}  

void coldStart(Complex gauge[L][L][D],param_t p){

  //#pragma omp parallel for
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      for(int mu=0;mu<D;mu++)
	gauge[x][y][mu] = Complex(1.0,0.0);
  return;
}  


double measPlaq(Complex  gauge[L][L][D]){

  Complex w(0.0,0.0);
  double plaq = 0.0;
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++){
      w = gauge[x][y][0]*gauge[ (x+1)%L ][y][1]*conj(gauge[x][ (y+1)%L ][0])*conj(gauge[x][y][1]);
      plaq += real(w)/(L*L);
    }
  return plaq;
}


/*=============================================================
  Creutz     exp[ -sigma L^2] exp[ -sigma(L-1)(L-1)]
  ratio:    ---------------------------------------  = exp[ -sigma]
  exp[ -sigma (L-1)L] exp[-sigma L(L-1)]
  ================================================================*/

void areaLaw(Complex gauge[L][L][D], Complex avWc[L][L], param_t p){

  Complex w;
  Complex Smeared[L][L][D];
  smearLink(Smeared,gauge,p);
  
  for(int xL=0; xL<L/2; xL++) 
    for(int yL=0; yL<L/2; yL++){
      for(int x=0; x<L; x++)
	for(int y=0; y<L; y++){
	  w = Complex(1.0,0.0);
	  
	  for(int dx=0; dx<xL; dx++)     w *=      Smeared[ (x+dx)%L ][y][0];
	  for(int dy=0; dy<yL; dy++)     w *=      Smeared[ (x+xL)%L ][ (y+dy)%L ][1];
	  for(int dx=xL-1; dx>-1; dx--)  w *= conj(Smeared[ (x+dx)%L ][ (y+yL)%L ][0]);
	  for(int dy=yL-1; dy>-1; dy--)  w *= conj(Smeared[x][ (y+dy)%L ][1]);
	  
	  avWc[xL][yL] += (1.0/(L*L))*w;
	};
    };
  
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

/* =========================================================================
   Brower et al claim  Wilson = exp[ - sigma size^2] with sigma \simeq 1/(2 beta)
   beta = 2 beta_old  ????  Looks like beta = beta_old
   ============================================================================*/
void areaFast(const Complex  gauge[L][L][D],Complex avW[L][L]){
  
  Complex w;
  for(int xL=0 ; xL < L/2; xL++) {
    for(int x =0;x< L;x++)
      for(int y =0;y< L;y++){
	w = Complex(1.0,0.0);
	for(int dx = 0;dx<xL;dx++)    w = w * gauge[(x + dx)%L][y][0];
	for(int dy = 0;dy<xL;dy++)    w = w * gauge[(x + xL)%L][(y+ dy)%L][1];
	for(int dx =xL-1;dx>-1 ;dx--) w = w * conj(gauge[(x + dx)%L][(y+xL)%L][0]);
	for(int dy =xL-1;dy>-1 ;dy--) w = w * conj(gauge[x][(y+ dy)%L][1]);
	avW[xL][xL] += (1.0/(L*L))*w;
      }
  };
}

void writeGaugeLattice(Complex gauge[L][L][D], string name){

  fstream outPutFile;
  outPutFile.open(name,ios::in|ios::out|ios::trunc);  
  outPutFile.setf(ios_base::fixed,ios_base::floatfield); 

  //Plaquette action header
  outPutFile << setprecision(20) <<  setw(20) << measPlaq(gauge) << endl;
  
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      for(int mu=0; mu<D; mu++)
	outPutFile << setprecision(12) <<  setw(20) << arg(gauge[x][y][mu]) << endl;
  
  outPutFile.close();
  return;
  
}
void readGaugeLattice(Complex gauge[L][L][D], string name){

  fstream inPutFile;
  inPutFile.open(name);
  string val;
  if(!inPutFile.is_open()) {
    cout << "Error opening file " << name << endl;
    exit(0);
  }
  
  //Header check
  getline(inPutFile, val);
  double plaq = stod(val);

  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      for(int mu=0; mu<D; mu++){
	getline(inPutFile, val);
	gauge[x][y][mu] = polar(1.0, stod(val));	  
      }
  cout << setprecision(16) << setw(20) << "Plaqette on file  = " << plaq << endl;
  cout << setprecision(16) << setw(20) << "Plaqette measured = " << measPlaq(gauge) << endl;
  double err = fabs(1.0 - plaq/measPlaq(gauge));
  if(err > 1e-12) {
    cout << "Gauge read fail!" <<  endl;
    exit(0);
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
    }
  return top/TWO_PI;
}


//staple x is 0th, y is 1st.
//APE smearing: project back on U(1)       
void smearLink(Complex Smeared[L][L][D], Complex gauge[L][L][D], param_t p){

  double alpha = p.alpha;
  int iter = p.smearIter;
  
  Complex SmearedTmp[L][L][D];
  copyLat(Smeared,gauge);
  copyLat(SmearedTmp,Smeared);
  
  for(int i=0; i<iter; i++) {
    
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++){
	SmearedTmp[x][y][0] += alpha * Smeared[x][y][1] * Smeared[x][(y+1)%L][0] * conj(Smeared[x][y][1]);
	SmearedTmp[x][y][0] += alpha * conj(Smeared[x][(y-1 +L)%L][1]) * Smeared[x][(y-1 +L)%L][0] * Smeared[(x+1)%L][(y-1 +L)%L][1];
	SmearedTmp[x][y][1] += alpha * Smeared[x][y][0]* Smeared[(x+1)%L][y][1] * conj(Smeared[(x+1)%L][(y+1)%L][0]);
	SmearedTmp[x][y][1] += alpha * conj(Smeared[(x-1+L)%L][y][0]) * Smeared[(x-1+L)%L][y][1] * Smeared[(x-1+L)%L][(y+1)%L][0];
      }

    //Project back to U(1)
    for(int x =0;x< L;x++)
      for(int y =0;y< L;y++)
	for(int mu = 0;mu <D;mu++)
	  Smeared[x][y][mu] = polar(1.0,arg(SmearedTmp[x][y][mu]));
  }
}																   

/* ==================================HMC===========================
   momenta conj to theta.  dU/dtheta = i U = I * U
*/


int hmc(Complex gauge[L][L][D], param_t p, int iter) {

  static int stop = 0;  
  int i;
  int accept = 0;
  
  double mom[L][L][D];
  Complex gaugeOld[L][L][D];
  Complex chi[L][L], eta[L][L];
  double H, Hold, rmsq;

  copyLat(gaugeOld,gauge);
  zeroLat(mom); 
  zeroField(chi);
  zeroField(eta);
  H = 0.0;
  Hold = 0.0;

  // init mom[L][L][D]  <mom^2> = 1;
  gaussReal_F(mom); 
  
  if(!p.quenched) {
    rmsq = 0.5;
    gaussComplex_F(eta,rmsq,p);  //  chi[L][L]
    Dpsi(chi, eta, gauge, p);
    
    for(int x = 0;x< L;x++)  //Masks out odd sites.
      for(int y = 0;y< L;y++)
	if((x + y)%2 == 1) chi[x][y] = 0.0;
  }
  
  Hold = calcH(mom, gauge, chi, p);
  trajectory(mom, gauge, chi, p); // MD trajectory using Verlet
  H = calcH(mom, gauge, chi, p);

  //cout << "Iter = " << iter << " H - Hold = " << H << " - " << Hold << " = " << H - Hold << endl;

  expcount++;
  expdHAve += exp(-(H-Hold));
  
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

void gaussReal_F(double field[L][L][D]) {
  //normalized gaussian exp[ - phi*phi/2]  <eta^2> = 1
  double r, theta;
  double sum = 0.0;
  for(int x=0; x<L;x++)
    for(int y= 0; y<L; y++){
      r = sqrt( -2.0*log(drand48()) );
      theta = TWO_PI*drand48();
      field[x][y][0] = r*cos(theta);
      field[x][y][1] = r*sin(theta);

      //sum0 += (field[x][y][1]*field[x][y][1] + field[x][y][0]*field[x][y][0]);
    }
  //count0++;
  //cout << sum0 /(count0*L*L) << endl;

  return;
}

void gaussComplex_F(Complex eta[L][L], double w, param_t p)
{
  //normalized gaussian exp[ - eta*eta/2]  <eta^2> = 1;
  int i;
  double r, theta;

  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      {
	r = sqrt( - 2.0*w*log(drand48()) );
	//	  r = sqrt( - 2.0*w*log((double)(rand())/RAND_MAX) );
	theta   = TWO_PI*(double)(rand())/RAND_MAX;
	eta[x][y] = Complex(r*cos(theta),r*sin(theta));      //eta[i] = cexp(theta*I);  works too!
      };
  return;
}

/*=================
  
  H(mom,theta) = \sum_i mom^2(i)/2 + V(theta_i)  
  = \sum_i mom^2(i)/2 + beta*\sum_P real(1 - U_P)
  
  dtheta/dt = dH/dmom = mom
  dmom/dt = - dH/theta = F
  
  F = dmon/dt = dt^2/dtheta  (Newton's Law)
  
  ===============*/

double calcH(double mom[L][L][D], Complex gauge[L][L][D], Complex chi[L][L], param_t p)
{
  double H = 0.0;
  Complex plaq;
  Complex chitmp[L][L];
 
  for(int x=0; x<L;x++)
    for(int y=0; y<L; y++){
      
      plaq = gauge[x][y][0]*gauge[ (x+1)%L ][y][1]*conj(gauge[x][ (y+1)%L ][0])*conj(gauge[x][y][1]);
      H += p.beta*real(1.0 - plaq);
      
      for(int mu = 0;mu <D;mu++){
	H += mom[x][y][mu]* mom[x][y][mu]/2.0;
      }
    }
 

  if(!p.quenched)
    {// cout << "Before Fermion force H = " << H << endl;
      Complex scalar = Complex(0.0,0.0);
      zeroField(chitmp);
      Ainv_psi(chitmp,chi,chitmp, gauge,p);  
      for(int x =0;x< L;x++)
	for(int y =0;y< L;y++){
	  if((x+y)%2 ==0)
	    scalar += conj(chi[x][y])*chitmp[x][y];
	}
      //      H +=  real(dotField(chi,chitmp));
      H += real(scalar);
      //       cout << "After Fermion Force H  = " << H << endl;
    }
  return H;
}


void trajectory(double mom[L][L][D], Complex gauge[L][L][D], Complex chi[L][L], param_t p) {

  int i, step;
  double fV[L][L][D];
  double fD[L][L][D];
  
  for(step = 0; step < p.nstep; step++) {
    
    if(step == 0) {
      for(int x =0;x< L;x++)
	for(int y =0;y< L;y++)
	  for(int mu = 0;mu <D;mu++)
	    gauge[x][y][mu] *= polar(1.0,mom[x][y][mu] * p.dt/2.0);
    }
    else {
      for(int x =0;x< L;x++)
	for(int y =0;y< L;y++)
	  for(int mu = 0;mu <D;mu++)
	    gauge[x][y][mu] *= polar(1.0,mom[x][y][mu] * p.dt);
    }
    
    forceV(fV, gauge, p);

    if(!p.quenched) forceD(fD,gauge,chi,p);
    
    if(p.quenched) {
      for(int x=0; x<L; x++)
	for(int y=0; y<L; y++)
	  for(int mu=0; mu<D; mu++)
	    mom[x][y][mu] +=  (fV[x][y][mu])*p.dt;
    } else {
      for(int x=0; x<L; x++)
	for(int y=0; y<L; y++)
	  for(int mu=0; mu<D; mu++)
	    mom[x][y][mu] += (fV[x][y][mu] + fD[x][y][mu])*p.dt;
    } 
  } //end for loop

  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      for(int mu = 0;mu <D;mu++){
	gauge[x][y][mu] *= polar(1.0,mom[x][y][mu] * p.dt/2.0);
      }
}

/*==================== 
  x+ = 0the dim
  y+ = 1th dim
  
  y+1<_______                              x-1/y+1<_____________>_x+1/y+1 
  |       ^				     |       ^      |
  |       |				     |       #      |
  |       |				     |       #      v
  v======>|				     v----->-#<-----|
  x/y       x+1                             x-1/y       x/y   x+1/y
  |       |				   
  ^       v				   
  y-1  <-------|				
  x/y       x+1

  - p.beta d/dtheat[x][y][mu] (1 -   Real[U_P])  = p.betas (I U  - I U^*)/2  =-  p.beta imag[U] 

  ============================*/
  
void forceV(double fV[L][L][D],Complex gauge[L][L][D], param_t p)
{
  Complex plaq ;
  zeroLat<double>(fV);
  //#pragma omp parallel for
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++) {
      
      plaq = gauge[x][y][0]*gauge[(x + 1)%L][y][1]*conj(gauge[x][(y+ 1)%L][0])*conj(gauge[x][y][1]);
      fV[x][y][0] +=  -  p.beta*imag(plaq);
      
      plaq =  gauge[x][y][0]*conj(gauge[(x + 1)%L][(y-1 +L)%L][1])*conj(gauge[x][(y -1 + L)%L][0])*gauge[x][(y -1 + L)%L][1];
      fV[x][y][0] +=  - p.beta*imag(plaq);
      
      plaq =  gauge[x][y][1]*conj(gauge[(x-1 + L)%L][(y+1)%L][0])*conj(gauge[(x-1+L)%L][y][1])*gauge[(x-1 +L)%L][y][0];
      fV[x][y][1] += - p.beta*imag(plaq);
      
      plaq =  gauge[x][y][1]*gauge[x][(y+1)%L][0]*conj(gauge[(x+1)%L][y][1])*conj(gauge[x][y][0]);
      fV[x][y][1] += - p.beta*imag(plaq);
      
      //	cout << x  << y <<" Forces = "<< fV[x][y][0] << "   " <<fV[x][y][1] << endl;
    }
}


/*
  
  VD(theta) = chi*_e (1/D D^dag) chi_e
  
  D chi = phi   and phi gausian.
  
  fD = - \dd_theta VD(theta) 
  =   chi*_e (1/D D^dag)_ee [ D^\dag (\dd_theta D) + (\dd_theta D^\dag) D ]_ee (1/D D^dag)_ee chi_e chi_e 


  *****  Should optimze this to operate only on EVEN sites. ****
  chi_even, chitmp_even, Dchitmp_odd
  
*/
void forceD(double fD[L][L][D],Complex gauge[L][L][D], Complex chi[L][L], param_t p)
{
  Complex chitmp[L][L];
  Complex Dchitmp[L][L];
  zeroField(chitmp);
  zeroField(Dchitmp);
  zeroLat(fD);
    
  Ainv_psi(chitmp,chi,chitmp, gauge,p); // note chitmp = 0 for ODD
  Dpsi(Dchitmp,chitmp,gauge,p); // restrict to Dslash, m = 0

  //#pragma omp parallel for
  for(int x = 0; x< L;x++)
    for(int y = 0;y < L; y++){
      if((x + y)%2 ==1) chitmp[x][y] = Complex(0.0,0.0);
      if((x + y)%2 ==0) Dchitmp[x][y]= Complex(0.0,0.0);
    }
    
  double eta1;
  for(int x=0; x<L;x++)
    for(int y=0; y<L;y++) {
      eta1 =(1.0 - 2.0*(x%2));
	
      if((x + y + 1)%2 == 0){ 
	fD[x][y][0] += 2.0*imag(conj(Dchitmp[x][y]) * gauge[x][y][0] * chitmp[(x+1)%L][y]);
      }
      else {
	fD[x][y][0] += 2.0*imag(conj(Dchitmp[(x+1)%L][y]) * conj(gauge[x][y][0]) * chitmp[x][y]);
      };
	
      if((x + y + 1)%2 == 0){    
	fD[x][y][1] += 2.0*eta1*imag(conj(Dchitmp[x][y]) * gauge[x][y][1] * chitmp[x][(y+1)%L]);
      }
      else {
	fD[x][y][1] += 2.0*eta1*imag(conj(Dchitmp[x][(y+1)%L]) * conj(gauge[x][y][1]) * chitmp[x][y]);
      }
    }	    
}

/*
  For a 2D square lattice, the stencil is:

  psi2[x] = D_xy psi_y = m psi_x delta_xy 
  - eta_mu(x) [U(x,x+mu) psi[x+mu] - U*(x-mu,x) psi[x-mu]]

  The even/odd anti-hermiticity             

  eta_mu(x) [U(x,y)\deta_x+mu,y - U(y,x) delta_y+mu,
 
  1 |  0 -eta1  0 |
  - | +1    0  -1 |  , where eta0 = 1, eta1 = (-)^x = 1 - 2*(x%2)
  2 |  0 +eta1  0 |

*/

void Dpsi(Complex psi2[L][L], Complex  psi1[L][L], Complex gauge[L][L][D], param_t p )
{
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++){
      psi2[x][y] = p.m*psi1[x][y];
      // cout<<"psi1[x,y] = ["<< x <<","<<y <<" ]  "<< psi1[x][y] << endl;
      // cout<<"Dpsi2[x,y] = ["<< x <<","<<y <<" ]  "<< psi2[x][y] << endl <<endl;
    }
  
  double eta1;
  for(int x=0; x<L; x++) {
    eta1 =(1 - 2*(x%2));
    for(int y=0; y<L; y++) {      
      psi2[x][y] += - (gauge[x][y][0] * psi1[(x+1)%L][y] - conj(gauge[ (x-1+L)%L ][y][0]) * psi1[ (x-1+L)%L ][y]);
      psi2[x][y] += - eta1*(gauge[x][y][1]*psi1[x][ (y+1)%L ] - conj(gauge[x][ (y-1+L)%L ][1])*psi1[x][ (y-1+L)%L ]);
      // cout<<"Dpsi2[x,y] = ["<< x <<","<<y <<" ]  "<< psi2[x][y] << endl <<endl;
    }
  }
}

void Ddagpsi(Complex psi2[L][L], Complex  psi1[L][L], Complex gauge[L][L][D], param_t p ) {
  
  // For a 2D square lattice, the stencil is:
  //   1 |  0 -eta1  0 |
  //   - | +1    0  -1 |  , where eta0 = 1, eta1 = (-)^x = 1 - 2(x%L)
  //   2 |  0 +eta1  0 |

  //#pragma omp parallel for
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      psi2[x][y] = p.m*psi1[x][y];
  
  double eta1;
  for(int x=0; x<L; x++) {
    eta1 =(1 - 2*(x%2));
    for(int y=0; y<L; y++) {   
      psi2[x][y] += (gauge[x][y][0] * psi1[(x+1)%L][y] - conj(gauge[(x-1 +L)%L][y][0]) * psi1[(x-1 + L)%L][y]);
      psi2[x][y] += eta1*(gauge[x][y][1]*psi1[x][(y+1)%L] - conj(gauge[x][(y-1+L)%L][1])*psi1[x][(y-1+L)%L]);
    }
  }
}
//=======================//
// Note: Ddag D = D Ddag //
///======================//

void DdagDpsi(Complex psi2[L][L], Complex  psi1[L][L], Complex gauge[L][L][D], param_t p )
{
  Complex psitmp[L][L];
  zeroField(psitmp);
  Dpsi(psitmp, psi1, gauge, p);
  Ddagpsi(psi2, psitmp, gauge, p);
}


/*============

  A_0 = - omega y   with omega = 2 pi Q/(L0*L1)

  A_1 = 0   except 2 pi Q x/L0   at y = L1 -1 : a gauge tranformaion at the boundary.

  F_10  = dd_1 A_0  - dd_0 A_1 = 2pi Q/(L0*L1) 

  4 by 4 lattice 
  ____________      0 omega
  |__|__|__|__|    -3 omega
  |__|__|__|__|   - 2 omega 
  |__|__|__|__|   - 1 omega
  |__|__|__|__|     0 omega  A0 = __

  omega = 2 pi Q/16 

  =============*/

void MakeInstanton(Complex gauge[L][L][D], int Q, param_t p)
{
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++){
      gauge[x][y][0] *= polar(1.0, - Q*TWO_PI*y/((double) L*L));
      if(y == L -1)  gauge[x][y][1] *=  polar(1.0, Q*TWO_PI*x/((double) L));
    }
}

/*====================

  Q = (1/2pi) int d^2x F_01 =  (1/2pi) int d^2x (dd_0 A_1 -  dd_1 A_0)
  =  -(1/2pi)  int (A_0 dx   + A_1 dy )
  = - (1/2pi) Q/r^2 int(y dx -x dy)
  =   Q int dtheta (sin^2  + cos^2)= Q

  This seem to work? Can displace it ot
  rx = x - x0 + 0.5
  ry = y - L/2 + 0.5;

  =====================*/

void MakePointVortex(Complex gauge[L][L][D], int Q, int x0, int y0, param_t p)
{  double rx, ry;
   
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++){
      rx = x - L/2 + 0.5;
      ry = y - L/2 + 0.5;
      gauge[x][y][0] *= polar(1.0, Q*ry/(rx*rx + ry*ry));
      gauge[x][y][1] *= polar(1.0, -Q*rx/(rx*rx + ry*ry));
    }
  // Shift location
  Complex tmp[L][L][D];
  copyLat(tmp, gauge);
  
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      for(int mu = 0; mu < 2; mu++) 
	gauge[x][y][mu] = tmp[(x + x0)%L][(y+ y0)%L][mu];
  
}

// Do a parity flip based on eqn 4 of arXiv:1203.2560v2
void DoPFlip(Complex gauge[L][L][D], param_t p)
{
  Complex parity_gauge[L][L][D];
  for (int x = 0; x < L; x++)
    {
      for (int y = 0; y < L; y++)
	{
	  parity_gauge[(-x-1+2*L)%L][y][0] = conj(gauge[x][y][0]);
	  parity_gauge[(-x+2*L)%L][y][1] = gauge[x][y][1];
	}
    }
  copyLat(gauge, parity_gauge);
}


double norm2(Complex psi[L][L]) {

  double norm2 = 0.0;
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      norm2 += psi[x][y].real() * psi[x][y].real() + psi[x][y].imag() * psi[x][y].imag();

  return norm2/(L*L);
}

Complex cDotProduct(Complex psi1[L][L], Complex psi2[L][L]) {

  Complex prod(0.0,0.0);
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      prod += conj(psi1[x][y])*psi2[x][y];

  return prod;
}


void caxpby(Complex a, Complex X[L][L], Complex b,
	    Complex Y[L][L], Complex result[L][L]){

  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      result[x][y] = a*X[x][y] + b*Y[x][y];
}

/*===================== 
  CG solutions to Apsi = b 
  see http://en.wikipedia.org/wiki/Conjugate_gradient_method
  ======================*/

int Ainv_psi(Complex psi[L][L],Complex b[L][L], Complex psi0[L][L], Complex gauge[L][L][D], param_t p)
{
  Complex res[L][L] , pvec[L][L], Apvec[L][L];
  double alpha, beta, denom ;
  double rsq = 0, rsqNew = 0, truersq=0.0, bsqrt = 0.0;
  
  //Intialize
  
  zeroField(res);
  zeroField(Apvec);
  zeroField(pvec);
  
  // Find norm of rhs.
  bsqrt =  real(dotField(b,b));
  bsqrt = sqrt(bsqrt);
  
  copyField(res,b); // res = b  - A psi0, for now start with phi0 = 0
  copyField(pvec, res);

  rsq =  real(dotField(res,res));
  
  // cout << "# Enter Ainv_p  bsqrt =  " << bsqrt << "  rsq   = " << rsq << endl;
  
  // Compute Ap.
  DdagDpsi(Apvec, pvec, gauge,p);
  
  // iterate till convergence
  int k;
  int success = 0;
  for(k = 0; k< p.maxIterCG; k++)
    {
      denom = real(dotField(pvec,Apvec));
      alpha = rsq/denom;
      
      for(int x =0;x<L;x++)
	for(int y =0;y<L;y++)
	  {
	    psi[x][y] +=  alpha * pvec[x][y];
	    res[x][y] += - alpha*Apvec[x][y];
	  }
      
      // Exit if new residual is small enough
      rsqNew = real(dotField(res,res));
      //   cout << "k = "<< k <<"   rsqNew = " << rsqNew << endl; 
      
      if (sqrt(rsqNew) < p.eps*bsqrt) {
	//    	printf("Final rsq = %g\n", rsqNew);
	break;
      }
      
      // Update vec using new residual
      beta = rsqNew / rsq;
      rsq = rsqNew;
      
      for(int x =0;x<L;x++)
	for(int y =0;y<L;y++)
	  pvec[x][y] = res[x][y] + beta * pvec[x][y];
      
      // Compute the new Ap.
      DdagDpsi(Apvec, pvec, gauge,p);  
    }
  
  
  if( k == p.maxIterCG) {
    printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq); 
    success = 0; /* Failed convergence */
  }
  else
    {
      success = 1; /* Convergence */
      k++;
    }

  DdagDpsi(Apvec,  psi, gauge,p);
  for(int x =0;x<L;x++)
    for(int y =0;y<L;y++)
      res[x][y] =  b[x][y]- Apvec[x][y];
  
  double truesq;
  truersq =  real(dotField(res,res));
  //  printf("CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  return success;
}

void TestCG(Complex gauge[L][L][D], param_t p)
{

  Complex b[L][L];
  double rmsq = 0.5;
  gaussComplex_F(b,rmsq,p);
  Complex psi1[L][L];
  Complex psi2[L][L];
  zeroField(psi1);
  zeroField(psi2);

  Ainv_psi(psi1,b,psi1, gauge,p);
  DdagDpsi(psi2, psi1, gauge,p);
  for(int x = 0;x < L;x++)
    for(int y=0;y<L;y++)
      cout<<" (x,y) = ("<< x <<","<<y<<") psi1 = " << psi1[x][y]
	  << "  b[x][y] =  " <<  b[x][y] << " psi2[x][y] = "<<  psi2[x][y] << endl;
    
  /*
    for(int x = 0;x < L;x++)
    for(int y=0;y<L;y++)
    cout<<" (x,y) = ("<< x <<","<<y<<") " << psi[x][y] << endl;
    DdagDpsi(psi, b, gauge,p);

    for(int x = 0;x < L;x++)
    for(int y=0;y<L;y++)
    cout<<" (x,y) = ("<< x <<","<<y<<") " << psi[x][y] << endl;
   
  
    Ainv_psi(psi,b,psi, gauge,p);
    for(int x = 0;x < L;x++)
    for(int y=0;y<L;y++)
    cout<< " (x,y) = ("<< x <<","<<y<<") "<< psi[x][y] << endl;
  */

}

void TestForce(Complex gauge[L][L][D], param_t p)
{
   
  Complex chi[L][L], eta[L][L];
  double fV[L][L][D], fD[L][L][D], mom[L][L][D];
  Complex Deltagauge[L][L][D];
  zeroField(chi);
  zeroField(eta);
  zeroLat(mom);
  zeroLat(Deltagauge);


  double  rmsq = 0.5;
  gaussComplex_F(eta,rmsq,p);  //  chi[L][L]
  Ddagpsi(chi, eta, gauge, p);  
  for(int x = 0;x< L;x++)  //Masks out odd sites.
    for(int y = 0;y< L;y++)
      if((x + y)%2 == 1) chi[x][y] = 0.0;

  double Hnew, Hold;
  Hnew = 0.0; Hold = 0.0;
  double delta = 0.00001;

  copyLat(Deltagauge,gauge);
  
#if 1
  p.beta = 1.0;
  p.quenched = true;
  forceV(fV, gauge, p);
  Hold = calcH(mom, gauge,chi, p);

  copyLat(Deltagauge,gauge);
 
  for(int x=0;x<L;x++)
    for(int y=0;y<L;y++)
      {

	Deltagauge[x][y][0] =  polar(1.0,delta)*gauge[x][y][0];
	Hnew = calcH(mom,  Deltagauge,chi, p);
	cout<<"fV[x,y,mu] = [ "<< x <<","<< y <<", 0 ] "<< fV[x][y][0] <<" <--> " << - (Hnew - Hold)/delta << endl;
        Deltagauge[x][y][0] =  gauge[x][y][0];
	
	Deltagauge[x][y][1] =  polar(1.0,delta)*gauge[x][y][1];
	Hnew = calcH(mom,  Deltagauge,chi, p);
	cout<<"fV[x,y,mu] = ["<< x <<","<< y <<", 1 ]  "<< fV[x][y][1] <<" <--> " << - (Hnew - Hold)/delta << endl << endl;
	Deltagauge[x][y][1] =  gauge[x][y][1];
      }
#endif

 

   
  p.beta = 0.0;
  p.quenched = false;
  forceD(fD,gauge,chi,p);

  cout << "Hold = " << Hold << endl;
  delta = 0.00001;
  Complex chitmp[L][L];
  zeroField(chitmp);

  // p.m = 0.0;
#if 0 
  for(int x=0;x<L;x++)
    for(int y=0;y<L;y++) 
      cout<<"chi[x,y] = ["<< x <<","<<y <<" ]  "<< chi[x][y] << endl;

  cout << " CALL WITH CHI EVEN "<< endl;
  Dpsi(chitmp,chi, gauge,p);
  for(int x=0;x<L;x++)
    for(int y=0;y<L;y++) 
      cout<<"Dchi[x,y] = ["<< x <<","<<y <<" ]  "<< chitmp[x][y] << endl <<endl;

  Ddagpsi(chitmp,chi, gauge,p);
  for(int x=0;x<L;x++)
    for(int y=0;y<L;y++) 
      cout<<"Ddagchi[x,y] = ["<< x <<","<<y <<" ]  "<< chitmp[x][y] << endl <<endl;
 
  DdagDpsi(chitmp,chi, gauge,p);
  for(int x=0;x<L;x++)
    for(int y=0;y<L;y++) 
      cout<<"DdagDchi[x,y] = ["<< x <<","<<y <<" ]  "<< chitmp[x][y] << endl <<endl;

 
  zeroField(chitmp);
 
  Ainv_psi(chitmp,chi, chitmp, gauge, p);

  for(int x=0;x<L;x++)
    for(int y=0;y<L;y++) 
      cout<<"Ainv_chi[x,y] = ["<< x <<","<<y <<" ]  "<< chitmp[x][y] << endl <<endl;

#endif

  Hold = calcH(mom, gauge, chi, p);
  //   cout << " Original Hold = " << Hold << endl;
  for(int x=0;x<L;x++)
    for(int y=0;y<L;y++)
      {
	Deltagauge[x][y][0] =  polar(1.0,delta)*gauge[x][y][0];
	Hnew = calcH(mom,  Deltagauge,chi, p);
	cout<<"fD[x,y,mu] = ["<< x <<","<<y <<", 0 ]  "<< fD[x][y][0] <<" <--> "  << - (Hnew - Hold)/delta << endl;
	Deltagauge[x][y][0] =  gauge[x][y][0];
       
     
        Deltagauge[x][y][1] =  polar(1.0,delta)*gauge[x][y][1];
	Hnew = calcH(mom,  Deltagauge,chi, p);  
	cout<<"fD[x,y,mu] = ["<< x <<","<<y <<", 1 ]  "<< fD[x][y][1] <<" <-->  " << - (Hnew - Hold)/delta << endl << endl;
	Deltagauge[x][y][1] =  gauge[x][y][1];
      }

}

#ifdef USE_ARPACK
void polyOp(const Dirac &mat,
	    cudaColorSpinorField &out,
	    const cudaColorSpinorField &in,	   
	    QudaArpackParam *arpack_param) {
  
  Float delta,theta;
  Float sigma,sigma1,sigma_old;
  Float d1,d2,d3;
  
  double a = arpack_param->amin;
  double b = arpack_param->amax;
  
  delta = (b-a)/2.0;
  theta = (b+a)/2.0;
  
  sigma1 = -delta/theta;
  
  blas::copy(out,in);
  
  if(arpack_param->polyDeg == 0)
    return;
  
  d1 =  sigma1/delta;
  d2 =  1.0;

  if(arpack_param->useNormOp && arpack_param->useDagger) {
    mat.MMdag(out,in);
  }
  else if(arpack_param->useNormOp && !arpack_param->useDagger) {
    mat.MdagM(out,in);
  }
  else if (!arpack_param->useNormOp && arpack_param->useDagger) {
    mat.Mdag(out,in);
  }
  else {  
    mat.M(out,in);
  }
  
  blas::axpby(d2, const_cast<cudaColorSpinorField&>(in), d1, out);
  
  //C_1(x) = x
  if(arpack_param->polyDeg == 1 )
    return;
  
  // C_0 is the current 'in'  vector.
  // C_1 is the current 'out' vector.
  
  //Clone 'in' to two temporary vectors.
  cudaColorSpinorField *tmp1 = new cudaColorSpinorField(in);
  cudaColorSpinorField *tmp2 = new cudaColorSpinorField(in);
  
  blas::copy(*tmp1,in);
  blas::copy(*tmp2,out);
  
  //Using Chebyshev polynomial recursion relation,
  //C_{m+1}(x) = 2*x*C_{m} - C_{m-1}
  
  sigma_old = sigma1;
  
  //construct C_{m+1}(x)
  for(int i=2; i < arpack_param->polyDeg; i++){
    
    sigma = 1.0/(2.0/sigma1-sigma_old);
    
    d1 = 2.0*sigma/delta;
    d2 = -d1*theta;
    d3 = -sigma*sigma_old;
    
    //mat*C_m
    if(arpack_param->useNormOp && arpack_param->useDagger) {
      mat.MMdag(out, *tmp2);
    }
    else if(arpack_param->useNormOp && !arpack_param->useDagger) {
      mat.MdagM(out, *tmp2);
    }
    else if (!arpack_param->useNormOp && arpack_param->useDagger) {
      mat.Mdag(out, *tmp2);
    }
    else {  
      mat.M(out, *tmp2);
    }
    
    blas::ax(d3,*tmp1);
    std::complex<double> d1c(d1,0.0);
    std::complex<double> d2c(d2,0.0);
    blas::cxpaypbz(*tmp1,d2c,*tmp2,d1c,out);
    
    blas::copy(*tmp1,*tmp2);
    blas::copy(*tmp2,out);
    sigma_old = sigma;
    
  }
  
  delete tmp1;
  delete tmp2;
}


void arpack_solve_double(Complex gauge[L][L][D], param_t p, Complex guess[L][L]) {
  
  //Construct parameters and memory allocation
  //------------------------------------------
  
  // all FORTRAN communication uses underscored 
  int ido_; 
  int info_;
  int *ipntr_ = (int*)malloc(14*sizeof(int));
  int *iparam_ = (int*)malloc(11*sizeof(int));
  int n_    = L*L,
    nev_    = p.nEv,
    nkv_    = p.nKv,
    ldv_    = L*L,
    lworkl_ = (3 * nkv_*nkv_ + 5*nkv_) * 2,
    rvec_   = 1;
  int max_iter = p.arpackMaxiter;

  double tol_ = p.arpackTol;

  Complex Zero(0.0,0.0);
  
  double *mod_evals_sorted  = (double*)malloc(nkv_*sizeof(double));
  int *evals_sorted_idx = (int*)malloc(nkv_*sizeof(int));
  
  //Memory checks
  if((mod_evals_sorted == nullptr) ||
     (evals_sorted_idx == nullptr) ) {
    printf("eigenSolver: not enough memory for host eigenvalue sorting");
    exit(0);
  }
  
  //ARPACK workspace
  std::complex<double> sigma_ = 0.0;
  std::complex<double> *resid_ =
    (std::complex<double> *) malloc(ldv_*sizeof(std::complex<double>));
  std::complex<double> *w_workd_ =
    (std::complex<double> *) malloc(3*ldv_*sizeof(std::complex<double>));
  std::complex<double> *w_workl_ =
    (std::complex<double> *) malloc(lworkl_*sizeof(std::complex<double>)); 
  std::complex<double> *w_workev_=
    (std::complex<double> *) malloc(2*nkv_*sizeof(std::complex<double>));    
  double *w_rwork_  = (double *)malloc(nkv_*sizeof(double));    
  int *select_ = (int*)malloc(nkv_*sizeof(int));
  
  std::complex<double> *evecs = (std::complex<double> *) malloc(nkv_*L*L*sizeof(std::complex<double>));
  std::complex<double> *evals = (std::complex<double> *) malloc(nkv_    *sizeof(std::complex<double>));

  for(int n=0; n<nkv_; n++) {
    evals[n] = Zero;
    for(int x=0; x<L; x++) {
      for(int y=0; y<L; y++) {
	evecs[n*L*L + x*L + y] = Zero;
	resid_[y*L + x] = guess[x][y];
      }
    }
  }
      
  //Alias pointers
  std::complex<double> *evecs_ = nullptr;
  evecs_ = (std::complex<double>*) (double*)(evecs);    
  std::complex<double> *evals_ = nullptr;
  evals_ = (std::complex<double>*) (double*)(evals);
  
  //Memory checks
  if((iparam_ == nullptr) ||
     (ipntr_ == nullptr) || 
     (resid_ == nullptr) ||  
     (w_workd_ == nullptr) || 
     (w_workl_ == nullptr) ||
     (w_workev_ == nullptr) ||
     (w_rwork_ == nullptr) || 
     (select_ == nullptr) ) {
    printf("eigenSolver: not enough memory for ARPACK workspace.\n");
    exit(0);
  }    

  //Assign values to ARPACK params 
  ido_        = 0;
  info_       = 0;
  iparam_[0]  = 1;
  iparam_[1]  = 1;
  iparam_[2]  = max_iter;
  iparam_[3]  = 1;
  iparam_[6]  = 1;
  //iparam_[7]  = 1;
  
  //ARPACK problem type to be solved
  char howmny='P';
  char bmat = 'I';
  char *spectrum;
  spectrum = strdup("SR"); //Initialsed just to stop the compiler warning...
  
  int iter_cnt= 0;

  bool allocate = true;

  //Start ARPACK routines
  //---------------------------------------------------------------------------------
 
  Complex *psi1;
  Complex *psi2;

  Complex psi1_cpy[L][L];
  Complex psi2_cpy[L][L];
  
  for(int x=0; x<L; x++) {
    for(int y=0; y<L; y++) {
      psi1_cpy[x][y] = Zero;
      psi2_cpy[x][y] = Zero;
    }
  }
  
  psi1 = w_workd_;
  psi2 = w_workd_ + L*L;

  double t1;
  double time = 0.0;;
  do {
    
    t1 = -((double)clock());
    
    //Interface to arpack routines
    //----------------------------
    
    ARPACK(znaupd)(&ido_, &bmat, &n_, spectrum, &nev_, &tol_, resid_, &nkv_,
		   evecs_, &n_, iparam_, ipntr_, w_workd_, w_workl_, &lworkl_,
		   w_rwork_, &info_, 1, 2);

    if (info_ != 0) {
      printf("\nError in znaupd info = %d. Exiting...\n",info_);
      arpackErrorHelpNAUPD();
      exit(0);
    }
    
    if (ido_ == 99 || info_ == 1)
      break;
    
    if (ido_ == -1 || ido_ == 1) {

      //Copy from Arpack workspace
      for(int y=0; y<L; y++) {
	for(int x=0; x<L; x++) {
	  psi1_cpy[x][y] = *(psi1 + x + y*L);
	}
      }
      //Apply matrix-vector operation
      DdagDpsi(psi2_cpy, psi1_cpy, gauge, p);
      
      //Copy to Arpack workspace
      for(int y=0; y<L; y++) {
	for(int x=0; x<L; x++) {
	  *(psi2 + x + y*L) = psi2_cpy[x][y];
	}
      }
    }
    
    t1 += clock();
    time += t1;
    if(iter_cnt % 100 == 0) printf("Arpack Iteration: %d (%e secs)\n", iter_cnt, time/CLOCKS_PER_SEC);
    iter_cnt++;
    
  } while (99 != ido_ && iter_cnt < max_iter);
  
  //Subspace calulated sucessfully. Compute nEv eigenvectors and values   
  printf("Finished in %e secs: iter=%04d  info=%d  ido=%d\n", time/CLOCKS_PER_SEC, iter_cnt, info_, ido_);      
  printf("Computing eigenvectors\n");
  
  //Interface to arpack routines
  //----------------------------
  ARPACK(zneupd)(&rvec_, &howmny, select_, evals_, evecs_, &n_, &sigma_,
		 w_workev_, &bmat, &n_, spectrum, &nev_, &tol_,
		 resid_, &nkv_, evecs_, &n_, iparam_, ipntr_, w_workd_,
		 w_workl_, &lworkl_, w_rwork_, &info_, 1, 1, 2);
  if (info_ == -15) {
    printf("\nError in zneupd info = %d. You likely need to\n"
	   "increase the maximum ARPACK iterations. Exiting...\n", info_);
    arpackErrorHelpNEUPD();
    exit(0);
  } else if (info_ != 0) {
    printf("\nError in zneupd info = %d. Exiting...\n", info_);
    arpackErrorHelpNEUPD();
  }
  
  //Print out the computed ritz values, absolute values, and their error estimates
  int nconv = iparam_[4];
  for(int j=0; j<nconv; j++){
  
    printf("RitzValue[%04d]  %+e  %+e  %+e  error= %+e \n",j,
	   real(evals_[j]),
	   imag(evals_[j]),
	   std::abs(evals_[j]),
	   std::abs(*(w_workl_ + ipntr_[10]-1+j)));
  
    evals_sorted_idx[j] = j;
    mod_evals_sorted[j] = std::abs(evals_[j]);
  }
  
  //Sort the eigenvalues in absolute ascending order
  t1 =  -((double)clock());;
  sortAbs(mod_evals_sorted, nconv, true, evals_sorted_idx);
  t1 +=  clock();
  
  //Print sorted evals
  printf("Sorting time: %f sec\n",t1/CLOCKS_PER_SEC);
  printf("Sorted eigenvalues based on their absolute values:\n");
  
  
  for(int j=0; j<nconv; j++){
    printf("RitzValue[%04d]  %+e  %+e  %+e  error = %+e \n",j,
	   real(evals_[evals_sorted_idx[j]]),
	   imag(evals_[evals_sorted_idx[j]]),
	   std::abs(evals_[evals_sorted_idx[j]]),
	   std::abs( *(w_workl_ + ipntr_[10] - 1 + evals_sorted_idx[j])) );
  }
  
  // Print additional convergence information.
  if( (info_) == 1){
    printf("Maximum number of iterations reached.\n");
  }
  else{
    if(info_ == 3){
      printf("Error: No shifts could be applied during implicit\n");
      printf("Error: Arnoldi update, try increasing NkV.\n");
    }
  }

  //Print Evalues
  t1 = -(double)clock();
  Complex psi3[L][L];
  
  for(int i =0 ; i < nev_ ; i++){

    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++) {
	psi3[x][y] = evecs[evals_sorted_idx[i]*L*L + L*y + x];
      }
    
    //apply matrix-vector operation here:
    DdagDpsi(psi2_cpy, psi3, gauge, p);

    //Check norms
    //cout << norm2(psi2_cpy) << " " << norm2(psi3) << endl;
    
    // lambda = v^dag * M*v    
    evals_[i] = cDotProduct(psi3, psi2_cpy);
    
    Complex unit(1.0,0.0);
    Complex m_lambda(-real(evals_[i]),
		     -imag(evals_[i]));
    
    // d_v = ||M*v - lambda*v||
    caxpby(unit, psi2_cpy, m_lambda, psi3, psi1_cpy);

    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++) {
	//cout << psi2_cpy[x][y] << " " << -m_lambda*psi3[x][y] << " " << psi1_cpy[x][y] << endl;
      }
    
    double L2norm = norm2(psi1_cpy);    
    printf("Eval[%04d] = %+e  %+e  Residual: %+e\n",
	   i, real(evals_[i]), imag(evals_[i]), sqrt(L2norm));
    
    
  }
  //exit(0);
  t1 += clock();
  printf("Eigenvalues of Dirac operator computed in: %f sec\n",
	 t1/CLOCKS_PER_SEC);
  
  // cleanup 
  free(ipntr_);
  free(iparam_);
  free(mod_evals_sorted);
  free(evals_sorted_idx);
  free(resid_);
  free(w_workd_);
  free(w_workl_);
  free(w_workev_);
  free(w_rwork_);
  free(select_);
  free(spectrum);
  
  return;  
}  
#endif
