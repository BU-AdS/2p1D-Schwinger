/* ================================================================================= 

   Sat Jan 22 13:43:09 PST 2005 Richard Brower  
   U(1) NON-COMPACT  Gauge Field for 2-d Quenched Lattice

   This does Heat Bath on the phase angles for the non-compact action:

   prob = Prod_{Plaquettes} = sqrt{PI/ beta} exp( - (beta/2)* (theta_Plaquette)^2 ]
                         
   The exact value of the average plaquette is

   Compared with exact value at beta = 16:
   - Log[Sqrt[begta/(2 Pi)] NIntegrate[2 Cos[x]*Exp[- beta*x^2/2],{x,0,Infinity}] = 1/32 =  0.03125

   This is not relevant to the real problem but the code was also debugged by
   running with purely Gaussian links at beta = 16: 

   sigma = Log[plaq]/4 = 0.0312485 (Simulation) Gaussian Links 
   sigma = Log[plaq] = 0.0312285 (Simulation) Gaussian Plaqs (Note factor of 4)


   Wed Aug 17 10:02:02 EDT 2016 Richard Brower
   U(1) COMPACT Full HMC for 2-d U(1)

   Optimization:

   Should introduce an ee CG for DdagD 
   Should introduce an D_oe Dslash.

   Check: Dslash, Inverters, Dforce.   

   ================================================================================= */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>
#include <omp.h>

using namespace std;
#include "ran2s.h"

#define L 16
#define LZ 5
#define D 3
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
  double betaz;
  double m;
  bool quenched;
  bool lockedZ;
  
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
void writeGaugeLattice(Complex gauge[L][L][LZ][D], string name);
void readGaugeLattice(Complex gauge[L][L][LZ][D], string name);
double measPlaq(Complex gauge[L][L][LZ][D], int z);
void gaussStart(Complex (&gauge)[L][L][LZ][D], param_t p);
void areaLaw(const Complex gauge[L][L][LZ][D], Complex (&avWc)[L][L][LZ]);
void polyakovLoops(Complex  gauge[L][L][LZ][D], Complex polyakov[L/2][LZ]);
double getTopCharge(Complex gauge[L][L][LZ][D], param_t p, int z);
void smearLink(Complex (&gaugeSmeared)[L][L][2], Complex gauge[L][L][2], param_t p);
void gaussReal_F(double (&field)[L][L][LZ][D]);

int hmc(Complex (&gauge)[L][L][LZ][D], param_t p, int iter);
double calcH(double mom[L][L][LZ][D], Complex gauge[L][L][LZ][D],
	     Complex chi[L][L], param_t p);
void trajectory(double (&mom)[L][L][LZ][D], Complex (&gauge)[L][L][LZ][D],
		Complex chi[L][L], param_t p);
void forceV(double (&fV)[L][L][LZ][D], double (&mom)[L][L][LZ][D],
	    Complex (&gauge)[L][L][LZ][D], param_t p);
void forceD(double fV[L][L][D], Complex gauge[L][L][D], Complex chi[L][L], param_t p);
void TestCG(Complex gauge[L][L][D], param_t p);

void checkGauge(Complex gauge[L][L][LZ][D], Complex gaugeAlt[L][L][LZ][D], int n);

void gaussComplex_F(Complex eta[L][L], double w, param_t p);

// Blas tools
double norm2(Complex psi[L][L]);
Complex cDotProduct(Complex psi1[L][L], Complex psi2[L][L]);
void caxpby(Complex a, Complex X[L][L], Complex b,
	    Complex Y[L][L], Complex result[L][L]);

// Fermion utils
// Dslash + m
void Dpsi(Complex psi2[L][L], Complex psi1[L][L], Complex gauge[L][L][D], param_t p );

// Dslash + m
void Ddagpsi(Complex psi2[L][L], Complex  psi1[L][L],
	     Complex gauge[L][L][D], param_t p );

// gamma_5
void gamma_5(Complex psi2[L][L], Complex psi1[L][L], param_t p );

// gamma_5 (Dslash + m) gamma_5 (Dslash + m)
void DdagDpsi(Complex psi2[L][L], Complex psi1[L][L],
	      Complex gauge[L][L][D], param_t p );

int Ainv_psi(Complex psi[L][L], Complex b[L][L], Complex psi0[L][L],
	     Complex gauge[L][L][D], param_t p);


// Zero lattice vectors.
template<typename T> inline void zeroLat(T v[L][L][LZ][D]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<D; mu++)
	  v[x][y][z][mu] = 0.0;
}

// Zero lattice vectors.
template<typename T> inline void zeroLat2D(T v[L][L][D]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int mu=0; mu<D; mu++)
	v[x][y][mu] = 0.0;
}


// Copy lattice vector
template<typename T> inline void copyLat(T v2[L][L][LZ][D],T v1[L][L][LZ][D]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0;mu<D;mu++)
	  v2[x][y][z][mu] =  v1[x][y][z][mu];
}

// Copy lattice vector
template<typename T> inline void copyLat2D(T v2[L][L][2],T v1[L][L][2]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int mu=0; mu<2; mu++)
	v2[x][y][mu] = v1[x][y][mu];
}

// Extract 2D slice
template<typename T> inline void extractLatSlice(T gauge[L][L][LZ][D],
						 T gauge2D[L][L][D],
						 int slice) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int mu=0;mu<D;mu++)
	gauge2D[x][y][mu] =  gauge[x][y][slice][mu];
}


// Insert 2D slice
template<typename T> inline void insertLatSlice(T gauge[L][L][LZ][D],
						T gauge2D[L][L][D],
						int slice) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int mu=0;mu<D;mu++)
	gauge[x][y][slice][mu] = gauge2D[x][y][mu];
}


// Copy lattice field
template<typename T> inline void copyField(T psi2[L][L][LZ], T psi1[L][L][LZ]) {

  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++)
	psi2[x][y][z] = psi1[x][y][z];
}

// Copy lattice field
template<typename T> inline void copyField2D(T psi2[L][L],T psi1[L][L]) {
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      psi2[x][y] =  psi1[x][y];
}


// Zero lattice field.
template<typename T> inline void zeroField(T psi[L][L][LZ]) {
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++) {
	psi[x][y][z] = 0.0;
      }
}

// Zero lattice field.
template<typename T> inline void zeroField2D(T psi[L][L]) {

  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      psi[x][y] = 0.0;
}


// Add Equ conj(v2) dot v1 lattice field
template<typename T> inline T dotField2D(T psi1[L][L], T psi2[L][L]) {
  T scalar = (T) 0.0;
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++)
      scalar += conj(psi1[x][y])*psi2[x][y];
  return scalar;
}

void printParams(param_t p) {
  cout << endl;
  cout << "Physics:  Size = "<< L << endl;
  cout << "          Beta = "<< p.beta << endl;
  cout << "          BetaZ = "<< p.betaz << endl;
  cout << "          Quenching = " << (p.quenched == true ? "True" : "False") << endl;
  cout << "          Z locked = " << (p.lockedZ == true ? "True" : "False") << endl;
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
  name += "_L" + to_string(p.Latsize) + "_LZ" + to_string(LZ) + "_B" + to_string(p.beta);
  if(p.quenched == false) name += "_M"+ to_string(p.m);
  name += "_dt" + to_string(p.dt) + "_nHMCstep" + to_string(p.nstep);
}

void checkGauge(Complex gauge[L][L][LZ][D], Complex gaugeAlt[L][L][LZ][D], int n) {

#ifdef CHECK_GAUGE
  bool fail = false;
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<D; mu++)
	  if(gauge[x][y][z][mu] != gaugeAlt[x][y][z][mu]) {
	    cout << "Gauge Check fail at " << x << " " << y << " " << z << " " << mu << " from call " << n << endl;
	    fail = true;
	  }
  
  if(fail) {
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++)
	for(int z=0; z<LZ; z++)
	  for(int mu=0; mu<D; mu++) {
	    if(gauge[x][y][z][mu] != gaugeAlt[x][y][z][mu]) {
	      cout << "Gauge compare " << x << " " << y << " " << z << " " << mu << gauge[x][y][z][mu] << " " <<  gaugeAlt[x][y][z][mu] << endl;
	    }
	  }
    exit(0);
  }
#endif
}

double expdHave = 0.0;
int counthmc = 0;
double sumBM = 0.0;
int countBM = 0;

int main(int argc, char **argv) {

  param_t p;

  p.beta = atof(argv[1]);
  p.betaz= atof(argv[2]);
  p.iterHMC = atoi(argv[3]);
  p.therm = atoi(argv[4]);
  p.skip = atoi(argv[5]);
  p.chkpt = atoi(argv[6]);
  p.checkpointStart = atoi(argv[7]);
  if(p.checkpointStart > 0) p.iterHMC += p.checkpointStart;

  p.Latsize = L;
  p.nstep = 200;
  p.dt = 0.005;
  p.quenched = false;
  p.lockedZ = true;
 
  p.maxIterCG = 1000;
  p.m = 0.032;
  p.eps = 1E-6;

  p.smearIter = atoi(argv[8]);
  p.alpha = atof(argv[9]);

  long iseed = (long)atoi(argv[10]);
  
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
  int histL = 101;
  int top_int[LZ];
  int top_old[LZ];
  int top_stuck[LZ];
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

  double area[L], sigma[L/2][LZ];
  Complex polyakov[L/2][LZ];
  Complex gauge[L][L][LZ][D];  
  Complex avWc[L][L][LZ];
  for(int a=0; a<LZ; a++) 
    for(int i=0; i<L/2; i++) {
      sigma[i][a] = 0.0;
      polyakov[i][a] = 0.0;
      for(int j=0; j<L/2; j++) {
	avWc[i][j][a] = Complex(0.0,0.0);
      }
    }
  
  string name;
  double elapsed = 0.0;

  int count = 0;
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
      
      for(int z=0; z<LZ; z++) {

	//Topology
	top = getTopCharge(gauge, p, z);
	top_int[z] = round(top);
	name = "top_charge_Lz";
	name += to_string(z);
	constructName(name, p);
	name += ".dat";
	sprintf(fname, "%s", name.c_str());
	fp = fopen(fname, "a");
	fprintf(fp, "%d %d\n", iter, top_int[z]);
	fclose(fp);
	
	index[z] = top_int[z] + (histL-1)/2;
	histQ[z][index[z]] += 1;
	if(top_old[z] == top_int[z]) top_stuck[z]++;
	top_old[z] = top_int[z];

      }
      
      if( (iter+1)%p.skip == 0 && (iter+1) > p.therm) {

	count++;
	
	//Plaquette actions
	for(int a=0; a<LZ; a++) {
	  plaq[a] = measPlaq(gauge, a);
	  plaqSum[a] += plaq[a];
	}
      
	double time = time0 + clock();
      
	cout << fixed << setprecision(16) << iter+1 << " ";
	cout << time/(threads*CLOCKS_PER_SEC) << " ";
	cout << plaqSum[(LZ-1)/2]/count << " ";
	cout << (double)top_stuck[(LZ-1)/2]/(count*p.skip) << " ";
	cout << expdHave/counthmc << " ";
	cout << (double)accepted/(count*p.skip) << endl;
	
	//Dump to file
	name = "data";
	constructName(name, p);
	name += ".dat";
	sprintf(fname, "%s", name.c_str());
	fp = fopen(fname, "a");
	fprintf(fp, "%d %.15e %.15e %.15e %.15e %.15e\n",
		iter+1,
		time/(threads*CLOCKS_PER_SEC),
		plaqSum[(LZ-1)/2]/count,
		(double)top_stuck[(LZ-1)/2]/(count*p.skip),
		expdHave/counthmc,
		(double)accepted/(count*p.skip) );
	fclose(fp);

	//Checkpoint the gauge field
	if( (iter+1)%p.chkpt == 0) {	  
	  name = "gauge";
	  constructName(name, p);
	  name += "_traj" + to_string(iter+1) + ".dat";	
	  writeGaugeLattice(gauge,name);
	}
      
	//arpack_solve_double(gauge, p, guess);
	//exit(0);
      
	for(int z=0; z<LZ; z++) 
	  for(int i=0; i<L/2; i++)  polyakov[i][z] = 0.0;	
	polyakovLoops(gauge, polyakov);
	
	zeroField(avWc);
	areaLaw(gauge, avWc);
      
	for(int z=0; z<LZ; z++) {

	  //Plaquette actions
	  name = "plaq_Lz" + to_string(z);
	  constructName(name, p);
	  name += ".dat";
	  sprintf(fname, "%s", name.c_str());
	  fp = fopen(fname, "a");
	  fprintf(fp, "%d %.15e\n", iter+1, plaq[z]);	  
	  fclose(fp);

	  
	  // Creutz Ratios	  
	  for(int size=1; size<L/2; size++) {      
	    sigma[size][z] =
	   -log( abs( ( real(avWc[size][size][z])  *real(avWc[size-1][size-1][z]) )/ 
		      ( real(avWc[size-1][size][z])*real(avWc[size][size-1][z]  ))));
	  }
	  name = "creutz_Lz" + to_string(z);
	  constructName(name, p);
	  name += ".dat";
	  sprintf(fname, "%s", name.c_str());
	  fp = fopen(fname, "a");
	  fprintf(fp, "%d %.16e ", iter+1, -log(abs(plaq[z])) );
	  for(int size=2 ; size < L/2; size++) fprintf(fp, "%.16e ", sigma[size][z] );
	  fprintf(fp, "\n");
	  fclose(fp);

	  for(int sizex=2 ; sizex < L/2; sizex++)
	    for(int sizey=sizex-1 ; (sizey < L/2 && sizey <= sizex+1) ; sizey++) {
	      name = "rectWL_Lz" + to_string(z);
	      name += "_" + to_string(sizex) + "_" + to_string(sizey);
	      constructName(name, p);
	      name += ".dat";
	      sprintf(fname, "%s", name.c_str());
	      fp = fopen(fname, "a");
	      fprintf(fp, "%d %.16e %.16e\n", iter+1, real(avWc[sizex][sizey][z]), imag(avWc[sizex][sizey][z]));	    
	      fclose(fp);
	    }
	  
	  // Polyakov Loops
	  name = "polyakov_Lz" + to_string(z);
	  constructName(name, p);
	  name += ".dat";
	  sprintf(fname, "%s", name.c_str());
	  fp = fopen(fname, "a");
	  fprintf(fp, "%d ", iter + 1);
	  for(int size=1 ; size < L/2; size++)
	    fprintf(fp, "%.16e %.16e ",
		    real(polyakov[size][z]),
		    imag(polyakov[size][z]));
	  fprintf(fp, "\n");
	  fclose(fp);

	  name = "polyakov_ratios";
	  constructName(name, p);
	  name += ".dat";
	  sprintf(fname, "%s", name.c_str());
	  fp = fopen(fname, "a");
	  fprintf(fp, "%d ", iter + 1);
	  for(int size=1 ; size < L/2-1; size++)
	    fprintf(fp, "%.16e ", real(polyakov[size+1][z])/real(polyakov[size][z]));
	  fprintf(fp, "\n");
	  fclose(fp);
	  
	  //Topological Historgram
	  name = "top_hist_Lz" + to_string(z);
	  constructName(name, p);
	  name += ".dat";
	  sprintf(fname, "%s", name.c_str());
	  fp = fopen(fname, "w");
	  for(int i=0; i<histL; i++) fprintf(fp, "%d %d\n", i - (histL-1)/2, histQ[z][i]);
	  fclose(fp);	  
	}	
      }
    }  
  }
  
  return 0;
}

void writeGaugeLattice(Complex gauge[L][L][LZ][D], string name){

  fstream outPutFile;
  outPutFile.open(name,ios::in|ios::out|ios::trunc);  
  outPutFile.setf(ios_base::fixed,ios_base::floatfield); 

  //Plaquette action header
  for(int z=0; z<LZ; z++) outPutFile << setprecision(20) <<  setw(20) << measPlaq(gauge, z) << endl;
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<D; mu++)
	  outPutFile << setprecision(20) <<  setw(20) <<  arg(gauge[x][y][z][mu]) << endl;
  
  outPutFile.close();
  return;
}

void readGaugeLattice(Complex gauge[L][L][LZ][D], string name){

  fstream inPutFile;
  inPutFile.open(name);
  string val;
  
  //Header check
  double plaq[LZ];
  for(int z=0; z<LZ; z++) {
    getline(inPutFile, val);
    plaq[z] = stod(val);
  }
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<D; mu++) {
	  getline(inPutFile, val);
	  gauge[x][y][z][mu] = polar(1.0, stod(val));	  
	}

  for(int z=0; z<LZ; z++) {
    cout << setprecision(16) << setw(20) << "Plaqette["<<z<<"] on file  = " << plaq[z] << endl;
    cout << setprecision(16) << setw(20) << "Plaqette["<<z<<"] measured = " << measPlaq(gauge, z) << endl;
    double err = fabs(1.0 - plaq[z]/measPlaq(gauge, z));
    if(err > 1e-12) {
      cout << "Gauge read fail at slice " << z << endl;
      exit(0);
    }    
  }
  inPutFile.close();
  return;
}

// /*===============================================================================
//   Gaussian numbers with p(theta) = sqrt(beta/ 2 PI) exp( - beta* theta^2/2)
//   <Gaussian^2> = 1/beta  
//   Perimeter Law:  Wilson Loop = exp[ - 4 sigma L ]   sigma = - Log[ <cos(theta)> ]
//   ================================================================================*/ 
void gaussStart(Complex (&gauge)[L][L][LZ][D], param_t p){

  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<D; mu++) {
	  //gauge[x][y][z][mu] = polar(1.0,sqrt(1.0/p.beta)*rang());
	  gauge[x][y][z][mu] = polar(1.0, TWO_PI*drand48() );
	  if(p.lockedZ && mu == 2) gauge[x][y][z][mu] = 1.0;
	}  
  return;
}

double measPlaq(Complex gauge[L][L][LZ][D], int z){
  
  Complex w(0.0,0.0);
  double plaq = 0.0;
  double norm = 1.0/(L*L);
  
  Complex temp[L][L][2];
   for(int x=0; x<L; x++) 
    for(int y=0; y<L; y++){
      temp[x][y][0] = gauge[x][y][z][0];
      temp[x][y][1] = gauge[x][y][z][1];
    }

  for(int x=0; x<L; x++) {
    for(int y=0; y<L; y++){
      w = temp[x][y][0]*temp[ (x+1)%L ][y][1]*conj(temp[x][ (y+1)%L ][0])*conj(temp[x][y][1]);
      //cout << "("<<x<<")("<<y<<") = " << w << " -> " << sqrt(real(w)*real(w) + imag(w)*imag(w)) << endl;
      plaq += real(w)*norm;
    }
  }
  return plaq;
}

double getTopCharge(Complex gauge[L][L][LZ][D], param_t p, int z){
  Complex w;
  double top = 0.0;
  Complex Smeared[L][L][2];
  Complex temp[L][L][2];

  for(int x=0; x<L; x++) 
    for(int y=0; y<L; y++){
      Smeared[x][y][0] = gauge[x][y][z][0];
      Smeared[x][y][1] = gauge[x][y][z][1];
      temp[x][y][0] = gauge[x][y][z][0];
      temp[x][y][1] = gauge[x][y][z][1];
    }
  
  smearLink(Smeared,temp,p);
  
  for(int x=0; x<L; x++){
    for(int y=0; y<L; y++){
      w = Smeared[x][y][0]*Smeared[ (x+1)%L ][y][1]*conj(Smeared[x][ (y+1)%L ][0])*conj(Smeared[x][y][1]);      
      top += arg(w);  // -pi < arg(w) < pi  Geometric value is an integer?
      //if(y == 2 && x == 2 && z == 2) cout << arg(w) << endl;
    }
    //cout<<endl;
  }
  //cout<<endl;
  return top/TWO_PI;
}


// staple x is  0th y is 1th    
// APE smearing: project back on U(1)   

void smearLink(Complex (&Smeared)[L][L][2], Complex gauge[L][L][2], param_t p){

  double alpha = p.alpha;
  int iter = p.smearIter;
  
  Complex SmearedTmp[L][L][2];
  copyLat2D(Smeared,gauge);
  copyLat2D(SmearedTmp, Smeared);
  
  for(int i=0; i<iter; i++) {
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++){
	SmearedTmp[x][y][0] += alpha * Smeared[x][y][1] * Smeared[x][(y+1)%L][0] * conj(Smeared[x][y][1]);
	SmearedTmp[x][y][0] += alpha * conj(Smeared[x][(y-1 +L)%L][1]) * Smeared[x][(y-1 +L)%L][0] * Smeared[(x+1)%L][(y-1 +L)%L][1];
	SmearedTmp[x][y][1] += alpha * Smeared[x][y][0]* Smeared[(x+1)%L][y][1] * conj(Smeared[(x+1)%L][(y+1)%L][0]);
	SmearedTmp[x][y][1] += alpha * conj(Smeared[(x-1+L)%L][y][0]) * Smeared[(x-1+L)%L][y][1] * Smeared[(x-1+L)%L][(y+1)%L][0];
      }
    
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++)
	for(int mu=0; mu<2; mu++){
	  Smeared[x][y][mu] = polar(1.0, arg(SmearedTmp[x][y][mu]));
	}
  }
}

int hmc(Complex (&gauge)[L][L][LZ][D], param_t p, int iter) {

  static int stop = 0;  
  int i;
  int accept = 0;
  
  double mom[L][L][LZ][D];
  Complex gaugeOld[L][L][LZ][D];
  Complex gauge2D[L][L][D];
  Complex chi[L][L], eta[L][L];
  double H, Hold, rmsq;

  copyLat(gaugeOld, gauge);
  zeroLat(mom); 
  zeroField2D(chi);
  zeroField2D(eta);
  H = 0.0;
  Hold = 0.0;

  // init mom[L][L][LZ][D]  <mom^2> = 1;
  gaussReal_F(mom); 

  if(!p.quenched) {
    rmsq = 0.5;
    gaussComplex_F(eta, rmsq, p);  //  chi[L][L]
    extractLatSlice(gauge, gauge2D, 2);
    Dpsi(chi, eta, gauge2D, p);
    
    for(int x = 0;x< L;x++)  //Masks out odd sites.
      for(int y = 0;y< L;y++)
	if((x + y)%2 == 1) chi[x][y] = 0.0;
  }
  
  copyLat(gaugeOld, gauge);
  Hold = calcH(mom, gauge, chi, p);
  checkGauge(gauge, gaugeOld, iter);
  
  trajectory(mom, gauge, chi, p); // MD trajectory using Verlet

  copyLat(gaugeOld, gauge);
  H = calcH(mom, gauge, chi, p);
  checkGauge(gauge, gaugeOld, iter);
  //cout << "Iter = " << iter << " H - Hold = " << H << " - " << Hold << " = " << H - Hold << endl;

  counthmc++;
  expdHave += exp(-(H-Hold));

  //cout << expdHave/counthmc << endl;
  
  // Metropolis accept/reject step
  if (iter >= p.therm) {    
    if ( drand48() > exp(-(H-Hold)) ) copyLat(gauge,gaugeOld);
    else accept = 1;
  }
  
  return accept;
}

//Check... OK!
void gaussReal_F(double (&field)[L][L][LZ][D]) {
  //normalized gaussian exp[ - phi*phi/2]  <eta^2> = 1
  double r, theta;
  double sum = 0.0;

  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++){

	r = sqrt( -2.0*log(drand48()));
	theta = TWO_PI*drand48();
	field[x][y][z][0] = r*cos(theta);
	field[x][y][z][1] = r*sin(theta);

	r = sqrt( -2.0*log(drand48()));
	theta = TWO_PI*drand48();
	field[x][y][z][2] = r*cos(theta);
	
	//Check Box-Muller 
	//sumBM += field[x][y][z][0]*field[x][y][z][0] + field[x][y][z][1]*field[x][y][z][1] + field[x][y][z][2]*field[x][y][z][2];
	
      }
  
  //countBM++;
  //cout << sumBM/(L*L*LZ*countBM) << endl;
  //exit(0);
  return;
}

void gaussComplex_F(Complex eta[L][L], double w, param_t p) {

  //normalized gaussian exp[ - eta*eta/2]  <eta^2> = 1;
  int i;
  double r, theta;
  
  for(int x =0;x< L;x++)
    for(int y =0;y< L;y++) {
      r = sqrt( - 2.0*w*log(drand48()) );
      //r = sqrt( - 2.0*w*log((double)(rand())/RAND_MAX) );
      theta   = TWO_PI*(double)(rand())/RAND_MAX;
      eta[x][y] = Complex(r*cos(theta),r*sin(theta));
      //eta[i] = cexp(theta*I);  works too!
    };
  return;
}


//Check... OK!
double calcH(double mom[L][L][LZ][D], Complex gauge[L][L][LZ][D], Complex chi[L][L], param_t p) {

  double Hgauge = 0.0;
  double Hmom = 0.0;
  double Hferm = 0.0;
  Complex plaq;
  Complex chitmp[L][L];
  double beta = p.beta;
  double betaz= p.betaz;
  Complex gauge2D[L][L][D];
  
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++) {

	//x, y, -x, -y
	plaq = gauge[x][y][z][0]*gauge[ (x+1)%L ][y][z][1]*conj(gauge[x][ (y+1)%L ][z][0])*conj(gauge[x][y][z][1]);
	Hgauge += beta*real(1.0 - plaq);
	
	//x, z, -x, -z
	if(z != LZ-1) {
	  plaq = gauge[x][y][z][0]*gauge[ (x+1)%L ][y][z][2]*conj(gauge[x][y][ (z+1)%LZ ][0])*conj(gauge[x][y][z][2]);
	  Hgauge += betaz*real(1.0 - plaq);
	
	  //y, z, -y, -z
	  plaq = gauge[x][y][z][1]*gauge[x][ (y+1)%L ][z][2]*conj(gauge[x][y][ (z+1)%LZ ][1])*conj(gauge[x][y][z][2]);
	  Hgauge += betaz*real(1.0 - plaq);	
	}
	
	for(int mu=0; mu<D; mu++){
	  Hmom += mom[x][y][z][mu]*mom[x][y][z][mu]/2.0;
	}
	if(LZ == 1) Hmom -= mom[x][y][z][2]*mom[x][y][z][2]/2.0;
      }

  
  if(!p.quenched) {
    extractLatSlice(gauge, gauge2D, 2);
    // cout << "Before Fermion force H = " << H << endl;
    Complex scalar = Complex(0.0,0.0);
    zeroField2D(chitmp);
    Ainv_psi(chitmp, chi, chitmp, gauge2D, p);  
    for(int x =0;x< L;x++)
      for(int y =0;y< L;y++){
	if((x+y)%2 ==0)
	  scalar += conj(chi[x][y])*chitmp[x][y];
      }
    //H +=  real(dotField(chi,chitmp));
    Hferm += real(scalar);
    //cout << "After Fermion Force H  = " << H << endl;
  }
  
  
  return Hgauge + Hmom + Hferm;
}

void trajectory(double (&mom)[L][L][LZ][D], Complex (&gauge)[L][L][LZ][D],
		Complex chi[L][L], param_t p) {
  
  int dim = D;
  if(p.lockedZ) dim -= 1;
  
  int i, step;
  double fV[L][L][LZ][D];
  double fD[L][L][D];
  Complex gauge2D[L][L][D];
  
  //First 1/2 step.
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++) 
	for(int mu=0; mu<dim; mu++) 
	  gauge[x][y][z][mu] *= polar(1.0, mom[x][y][z][mu] * p.dt/2.0);

  for(step = 1; step < p.nstep; step++) {
    
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++)
	for(int z=0; z<LZ; z++)
	  for(int mu=0; mu<dim; mu++) 
	    gauge[x][y][z][mu] *= polar(1.0, mom[x][y][z][mu] * p.dt);
           
    forceV(fV, mom, gauge, p);
    if(p.quenched) {
      for(int x=0; x<L; x++)
	for(int y=0; y<L; y++)
	  for(int z=0; z<LZ; z++)
	    for(int mu=0; mu<dim; mu++)
	      mom[x][y][z][mu] += (fV[x][y][z][mu])*p.dt;
    }
    else {
      extractLatSlice(gauge, gauge2D, 2);
      forceD(fD, gauge2D, chi, p);
      
      for(int z=0; z<LZ; z++) {
	if(z == 2) {	
	  for(int x=0; x<L; x++)
	    for(int y=0; y<L; y++)
	      for(int mu=0; mu<D; mu++)	      
		mom[x][y][z][mu] += (fV[x][y][z][mu] + fD[x][y][mu])*p.dt;
	}
	else {
	  for(int x=0; x<L; x++)
	    for(int y=0; y<L; y++)
	      for(int mu=0; mu<D; mu++)	      
		mom[x][y][z][mu] += (fV[x][y][z][mu])*p.dt;
	}
      }
    }    
  } //end for loop

  //Final 1/2 step.
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      for(int z=0; z<LZ; z++)
	for(int mu=0; mu<dim; mu++){
	  gauge[x][y][z][mu] *= polar(1.0, mom[x][y][z][mu] * p.dt/2.0);	  
	}
}


//Check... OK!
void forceV(double (&fV)[L][L][LZ][D], double (&mom)[L][L][LZ][D],
	    Complex (&gauge)[L][L][LZ][D], param_t p) {

  //Zero out the forces
  zeroLat<double>(fV);
  
  Complex plaq ;
  double beta = p.beta;
  double betaz= p.betaz;
  double dt = p.dt;
  
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

	//if(x == 0 && y == 0 && z == 0) cout << gauge[x][y][z][0] << " " << gauge[x][y][z][1] << " " << gauge[x][y][z][2] << "    ";
	
	//X dir
	//-------	
	//x, y, -x, -y
	plaq = gauge[x][y][z][0] * gauge[xp1][y][z][1] * conj(gauge[x][yp1][z][0]) * conj(gauge[x][y][z][1]);
	fV[x][y][z][0] -= beta*imag(plaq);

	//x, -y, -x, y
       	plaq = gauge[x][y][z][0]*conj(gauge[xp1][ym1][z][1])*conj(gauge[x][ym1][z][0])*gauge[x][ym1][z][1];
	fV[x][y][z][0] -= beta*imag(plaq);

	if(z != LZ-1) {
	  //x, z, -x, -z
	  plaq = gauge[x][y][z][0] * gauge[xp1][y][z][2] * conj(gauge[x][y][zp1][0]) * conj(gauge[x][y][z][2]);
	  fV[x][y][z][0] -= betaz*imag(plaq);
	}
	  
	if(z != 0) {
	  //x, -z, -x, z
	  plaq = gauge[x][y][z][0] * conj(gauge[xp1][y][zm1][2])*conj(gauge[x][y][zm1][0])*gauge[x][y][zm1][2];
	  fV[x][y][z][0] -= betaz*imag(plaq);
	}
	
	//Y dir
	//------
	//y, x, -y, -x
	plaq = gauge[x][y][z][1]*gauge[x][yp1][z][0]*conj(gauge[xp1][y][z][1])*conj(gauge[x][y][z][0]);
	fV[x][y][z][1] -= beta*imag(plaq);
	
	//y, -x, -y, x
	plaq = gauge[x][y][z][1]*conj(gauge[xm1][yp1][z][0])*conj(gauge[xm1][y][z][1])*gauge[xm1][y][z][0];
	fV[x][y][z][1] -= beta*imag(plaq);

	if(z != LZ-1) {
	  //y, z, -y, -z
	  plaq = gauge[x][y][z][1]*gauge[x][yp1][z][2]*conj(gauge[x][y][zp1][1])*conj(gauge[x][y][z][2]);
	  fV[x][y][z][1] -= betaz*imag(plaq);
	}
	
	if(z != 0) {
	  //y, -z, -y, z
	  plaq = gauge[x][y][z][1]*conj(gauge[x][yp1][zm1][2])*conj(gauge[x][y][zm1][1])*gauge[x][y][zm1][2];
	  fV[x][y][z][1] -= betaz*imag(plaq);
	}
	
	//Z dir
	//------
	// Only update the z links if not locked
	if(z != LZ-1 && !p.lockedZ) {
	  //z, x, -z, -x
	  plaq = gauge[x][y][z][2]*gauge[x][y][zp1][0]*conj(gauge[xp1][y][z][2])*conj(gauge[x][y][z][0]);
	  fV[x][y][z][2] -= betaz*imag(plaq);
	  
	  //z, -x, -z, x
	  plaq = gauge[x][y][z][2]*conj(gauge[xm1][y][zp1][0])*conj(gauge[xm1][y][z][2])*gauge[xm1][y][z][0];
	  fV[x][y][z][2] -= betaz*imag(plaq);
	  
	  //z, y, -z, -y
	  plaq = gauge[x][y][z][2]*gauge[x][y][zp1][1]*conj(gauge[x][yp1][z][2])*conj(gauge[x][y][z][1]);
	  fV[x][y][z][2] -= betaz*imag(plaq);
	  
	  //z, -y, -z, y
	  plaq = gauge[x][y][z][2]*conj(gauge[x][ym1][zp1][1])*conj(gauge[x][ym1][z][2])*gauge[x][ym1][z][1];
	  fV[x][y][z][2] -= betaz*imag(plaq);
	}

	//if(x == 0 && y == 0 && z == 0) cout << gauge[x][y][z][0] << " " << gauge[x][y][z][1] << " " << gauge[x][y][z][2] << endl;
	//for(int mu=0; mu<D; mu++) mom[x][y][z][mu] += fV[x][y][z][mu]*dt;	
      }
    }
  }
}


void forceD(double fD[L][L][D],Complex gauge[L][L][D], Complex chi[L][L], param_t p) {

  Complex chitmp[L][L];
  Complex Dchitmp[L][L];
  zeroField2D(chitmp);
  zeroField2D(Dchitmp);
  zeroLat2D(fD);
  
  Ainv_psi(chitmp, chi, chitmp, gauge, p); // note chitmp = 0 for ODD
  Dpsi(Dchitmp, chitmp, gauge, p); // restrict to Dslash, m = 0

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
	fD[x][y][0] += 2.0*imag(conj(Dchitmp[x][y]) *
				gauge[x][y][0] * chitmp[(x+1)%L][y]);
      }
      else {
	fD[x][y][0] += 2.0*imag(conj(Dchitmp[(x+1)%L][y]) *
				conj(gauge[x][y][0]) * chitmp[x][y]);
      };
	
      if((x + y + 1)%2 == 0){    
	fD[x][y][1] += 2.0*eta1*imag(conj(Dchitmp[x][y]) *
				     gauge[x][y][1] * chitmp[x][(y+1)%L]);
      }
      else {
	fD[x][y][1] += 2.0*eta1*imag(conj(Dchitmp[x][(y+1)%L]) *
				     conj(gauge[x][y][1]) * chitmp[x][y]);
      }
    }	    
}




// ===================================================================
//  Creutz     exp[ -sigma L^2] exp[ -sigma(L-1)(L-1)]
//  ratio:    ---------------------------------------  = exp[ -sigma]
//             exp[ -sigma (L-1)L] exp[-sigma L(L-1)]
// ==================================================================

void areaLaw(const Complex gauge[L][L][LZ][D], Complex (&avWc)[L][L][LZ]){
  
  Complex w;
  int p1, p2, p3, p4;
  double denom = 1.0/(L*L);
  int YrectStart = 1;
  int YrectFin = 1;
  
  //Loop over z planes
  for(int z=0; z<LZ; z++) {

    //Loop over all X side sizes of rectangle 
    for(int Xrect=1; Xrect<L/2; Xrect++) {

      //Loop over all Y side sizes of rectangle
      //Xrect > 1 ? YrectStart = Xrect-1 : YrectStart = 1;
      for(int Yrect=1; Yrect<L/2; Yrect++) {

	//Loop over all x,y
	for(int x=0; x<L; x++)
	  for(int y=0; y<L; y++){
	    
	    w = Complex(1.0,0.0);

	    //if(x == 0 && y == 0 && z == 0) cout << gauge[x][y][z][0] << " " << gauge[x][y][z][1] << " " << gauge[x][y][z][2] << endl;
	    
	    //Move in +x up to p1.
	    for(int dx=0; dx<Xrect; dx++)    w *= gauge[ (x+dx)%L ][y][z][0];
	    
	    //Move in +y up to p2 (p1 constant)
	    p1 = (x + Xrect)%L;
	    for(int dy=0; dy<Yrect; dy++)    w *= gauge[p1][ (y+dy)%L ][z][1];
	    
	    //Move in -x from p1 to (p2 constant)
	    p2 = (y + Yrect)%L;
	    for(int dx=Xrect-1; dx>=0; dx--)  w *= conj(gauge[ (x+dx)%L ][p2][z][0]);
	    
	    //Move in -y from p2 to y
	    for(int dy=Yrect-1; dy>=0; dy--)  w *= conj(gauge[x][ (y+dy)%L ][z][1]);
	    avWc[Xrect][Yrect][z] += w*denom;
	  }
      }
    }
    
  }
  
  return;
}

//Polyakov loops. x is the spatial dim, y is the temporal dim.
void polyakovLoops(Complex  gauge[L][L][LZ][D], Complex polyakov[L/2][LZ]){

  double denom = 1.0/(L/2);  
  Complex w1, w2;
  
  //Loop over z planes
  for(int z=0; z<LZ; z++) {
    //Eack polyakov loop correlation is defined by its delta x value.
    //We start at x0, separate to x0 + (x0+L/2-1), and loop over all
    //x0=1 -> x0 = L/2-1.
    
    //Starting x
    for(int x=0; x<L/2; x++) {
      
      //Loop over time
      w1 = Complex(1.0,0.0);
      for(int dy=0; dy<L; dy++) w1 *= gauge[x][dy][z][1];
      
      //x separation
      for(int dx=0; dx<L/2; dx++) {	
	w2 = Complex(1.0,0.0);
	for(int dy=0; dy<L; dy++) w2 *= gauge[x + dx][dy][z][1];	
	polyakov[dx][z] += conj(w1)*w2*denom;	
      }
    }
  }
  return;
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

void Dpsi(Complex psi2[L][L], Complex psi1[L][L], Complex gauge[L][L][D], param_t p ){
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

void DdagDpsi(Complex psi2[L][L], Complex  psi1[L][L], Complex gauge[L][L][D], param_t p ) {
  
  Complex psitmp[L][L];
  zeroField2D(psitmp);
  Dpsi(psitmp, psi1, gauge, p);
  Ddagpsi(psi2, psitmp, gauge, p);
}

/*===================== 
  CG solutions to Apsi = b 
  see http://en.wikipedia.org/wiki/Conjugate_gradient_method
  ======================*/

int Ainv_psi(Complex psi[L][L],Complex b[L][L], Complex psi0[L][L], Complex gauge[L][L][D], param_t p) {
  Complex res[L][L] , pvec[L][L], Apvec[L][L];
  double alpha, beta, denom ;
  double rsq = 0, rsqNew = 0, truersq=0.0, bsqrt = 0.0;
  
  //Intialize  
  zeroField2D(res);
  zeroField2D(Apvec);
  zeroField2D(pvec);
  
  // Find norm of rhs.
  bsqrt =  real(dotField2D(b,b));
  bsqrt = sqrt(bsqrt);
  
  copyField2D(res,b); // res = b  - A psi0, for now start with phi0 = 0
  copyField2D(pvec, res);

  rsq =  real(dotField2D(res,res));
  
  // cout << "# Enter Ainv_p  bsqrt =  " << bsqrt << "  rsq   = " << rsq << endl;
  
  // Compute Ap.
  DdagDpsi(Apvec, pvec, gauge,p);
  
  // iterate till convergence
  int k;
  int success = 0;
  for(k = 0; k< p.maxIterCG; k++) {
    
    denom = real(dotField2D(pvec,Apvec));
    alpha = rsq/denom;
    
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++) {
	psi[x][y] +=  alpha * pvec[x][y];
	res[x][y] += - alpha*Apvec[x][y];
      }
    
    // Exit if new residual is small enough
    rsqNew = real(dotField2D(res,res));
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
  //End loop over k
  
  
  if(k == p.maxIterCG) {
    printf("CG: Failed to converge iter = %d, rsq = %e\n", k,rsq); 
    success = 0; // Failed convergence 
  }
  else {
    success = 1; // Convergence 
    k++;
  }

  DdagDpsi(Apvec,  psi, gauge,p);
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      res[x][y] = b[x][y] - Apvec[x][y];
  
  //truersq =  real(dotField2D(res,res));
  //printf("CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  return success;
}

void TestCG(Complex gauge[L][L][D], param_t p) {
  
  Complex b[L][L];
  double rmsq = 0.5;
  gaussComplex_F(b,rmsq,p);
  Complex psi1[L][L];
  Complex psi2[L][L];
  zeroField2D(psi1);
  zeroField2D(psi2);

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


