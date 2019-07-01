#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>

using namespace std;

void jackknife(double **arrayLsq, double **arrayLm1sq, double **arrayLLm1, double **arrayLm1L, int N, int B, int L);

// Computes the jackknifed Creutz Ratio

//   Creutz     exp[ -sigma L^2] exp[ -sigma(L-1)(L-1)]
//   ratio:    ---------------------------------------  = exp[ -sigma]
//              exp[ -sigma (L-1)L] exp[-sigma L(L-1)]

int main(int argc, char **argv) {

  cout << setprecision(16);
  
  if (argc != 8) {
    cout << "./jackknife <filenameL^2> <filenameLm1^2> <filenameLLm1> <filenameLm1L> <data points> <JK block size> <L value>" << endl;
    exit(0);
  }
  
  string nameLsq(argv[1]);   //Name of the L^2 input file
  string nameLm1sq(argv[2]); //Name of the Lm1^2 input file
  string nameLLm1(argv[3]);  //Name of the LLm1 input file
  string nameLm1L(argv[4]);  //Name of the Lm1L input file
  fstream inputFileLsq;
  fstream inputFileLm1sq;
  fstream inputFileLLm1;
  fstream inputFileLm1L;
  inputFileLsq.open(nameLsq);
  inputFileLm1sq.open(nameLm1sq);
  inputFileLLm1.open(nameLLm1);
  inputFileLm1L.open(nameLm1L);
  if(!inputFileLsq.is_open()) {
    cout << "Error opening file: " << nameLsq << endl;
    exit(0);
  }
  if(!inputFileLm1sq.is_open()) {
    cout << "Error opening file: " << nameLm1sq << endl;
    exit(0);
  }
  if(!inputFileLLm1.is_open()) {
    cout << "Error opening file: " << nameLLm1 << endl;
    exit(0);
  }
  if(!inputFileLm1L.is_open()) {
    cout << "Error opening file: " << nameLm1L << endl;
    exit(0);
  }
  
  int N = atoi(argv[5]); //Number of data points
  int B = atoi(argv[6]); //Number of jk blocks
  int L = atoi(argv[7]); //Value of L

  if( N%B != 0) {
    cout << "Please ensure <data points> " << N << " is divisible by <block size> " << B << endl;
    exit(0);
  }

  int T=3; //3 columns for the data point entry.
  double **arrayLsq, **arrayLm1sq, **arrayLLm1, **arrayLm1L;
  arrayLsq = new double*[N];
  arrayLm1sq = new double*[N];
  arrayLLm1 = new double*[N];
  arrayLm1L = new double*[N];
  for(int n=0; n<N; n++) {
    arrayLsq[n] = new double[T];
    arrayLm1sq[n] = new double[T];
    arrayLLm1[n] = new double[T];
    arrayLm1L[n] = new double[T];
  }

  for(int n=0; n<N; n++) {
    for(int t=0; t<T; t++) {
      inputFileLsq >> arrayLsq[n][t];
      inputFileLm1sq >> arrayLm1sq[n][t];
      inputFileLLm1 >> arrayLLm1[n][t];
      inputFileLm1L >> arrayLm1L[n][t];
      //cout << n << " " << t << " " << array[n][t] << " " << endl;
    }
  }
  
  //Jackknife and dump the data
  jackknife(arrayLsq, arrayLm1sq, arrayLLm1, arrayLm1L, N, B, L);
}

void jackknife(double **arrayLsq, double **arrayLm1sq, double **arrayLLm1, double **arrayLm1L, int N, int B, int L) {

  //First we compute the means from the entire data set
  //---------------------------------------------------------  
  //Get the arithmetic mean for each loop type
  double meanLsq = 0.0;
  double meanLm1sq = 0.0;
  double meanLLm1 = 0.0;
  double meanLm1L = 0.0;
  //Real data is in the t=1 column
  for(int n=0; n<N; n++) {
    meanLsq += arrayLsq[n][1];
    meanLm1sq += arrayLm1sq[n][1];
    meanLLm1 += arrayLLm1[n][1];
    meanLm1L += arrayLm1L[n][1];
  }

  meanLsq /= N;
  meanLm1sq /= N;
  meanLLm1 /= N;
  meanLm1L /= N;
  
  //Get the string tension from creutz ratio.
  double sigma = -log(meanLsq*meanLm1sq/(meanLLm1*meanLm1L));  
  //---------------------------------------------------------
  
  //We now compute the jackknife errors
  //---------------------------------------------------------  
  //Initialise JK
  int Nr = N/B;
  double coeff = (1.0*Nr - 1.0)/(1.0*Nr);
  double resamp_norm = 1.0/(N - B);
  double resampLsq[Nr];
  double resampLm1sq[Nr];
  double resampLLm1[Nr];
  double resampLm1L[Nr];
  double resampSigma[Nr];
  double jkErrSigma = 0.0;
  //all real data in t=1 column
  for(int i=0; i<Nr; i++) {
    resampLsq[i] = 0.0;
    resampLm1sq[i] = 0.0;
    resampLLm1[i] = 0.0;
    resampLm1L[i] = 0.0;
    resampSigma[i] = 0.0;
  }
  
  //Compute resampled means
  for(int i=0; i<Nr; i++) {
    for(int n=0; n<N; n++) {
      if( n<i*B || n>=(i+1)*B ) {
	resampLsq[i] += arrayLsq[n][1];
	resampLm1sq[i] += arrayLm1sq[n][1];
	resampLLm1[i] += arrayLLm1[n][1];
	resampLm1L[i] += arrayLm1L[n][1];
      }
    }
    resampSigma[i] += -log(resampLsq[i]*resampLm1sq[i]/(resampLLm1[i]*resampLm1L[i]));
        
  }

  //Get jk error on sigma 
  for(int i=0; i<Nr; i++) {
    jkErrSigma += pow(sigma - resampSigma[i],2);
  }
  jkErrSigma = sqrt(coeff*jkErrSigma);

  string name("sigma.dat"); //Name of the output file
  fstream outputFile;
  outputFile.open(name, fstream::out);
  if(!outputFile.is_open()) {
    cout << "Error opening file: " << name << endl;
    exit(0);
  }
  outputFile << L << setprecision(16) << " " << sigma << " " << jkErrSigma << endl;
  outputFile.close();
  
}

