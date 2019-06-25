#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#include <complex>

using namespace std;

void jackknife(double **array, int N, int T, int B, int WL);

int main(int argc, char **argv) {

  cout << setprecision(16);
  
  if (argc < 5 || argc > 6) {
    cout << "./jackknife <filename> <data points> <columns> <JK block size> OPTIONAL:<Wilson loop size>" << endl;
    exit(0);
  }
  
  string name(argv[1]);   //Name of the input file
  string val;
  fstream inputFile;
  inputFile.open(name);
  if(!inputFile.is_open()) {
    cout << "Error opening file: " << name << endl;
    exit(0);
  }
  
  int N = atoi(argv[2]);   //Number of data points
  int T = atoi(argv[3]);   //Number of cols in data point
  int B = atoi(argv[4]);   //Number of jk blocks
  int WL = atoi(argv[5]);  //Wilson Loop size
  
  if( N%B != 0) {
    cout << "Please ensure <data points> " << N << " is divisible by <block size> " << B << endl;
    exit(0);
  }

  T++; //Extra column for the data point entry.
  double **array;
  array = new double*[N];
  for(int n=0; n<N; n++) array[n] = new double[T];
  
  for(int n=0; n<N; n++) {
    for(int t=0; t<T; t++) {
      inputFile >> array[n][t];
      //cout << n << " " << t << " " << array[n][t] << " " << endl;
    }
  }

  //Jackknife and dump the data
  jackknife(array, N, T, B, WL);
}

void jackknife(double **array, int N, int T, int B, int WL) {

  //First we compute the mean and the effective mass from the
  //entire data set
  //---------------------------------------------------------  
  //Get the arithmetic mean
  double mean[T];
  for(int t=0; t<T; t++) {
    mean[t] = 0.0;
    for(int n=0; n<N; n++) {
      mean[t] += array[n][t];
    }
    mean[t] /= N;
  }

  //Get the effective masses
  double meff[T];
  meff[0] = 0;
  for(int t=1; t<T; t++) {
    meff[t] = -log(mean[t]/mean[t-1]);
  }
  //---------------------------------------------------------

  
  //We now compute the jackknifed means and effective masses
  //---------------------------------------------------------  
  //Initialise JK
  int Nr = N/B;
  double coeff = (1.0*Nr - 1.0)/(1.0*Nr);
  double resamp_norm = 1.0/(N - B);
  double resamp_mean[T][Nr];
  double resamp_meff[T][Nr];
  double jk_err_mean[T];
  double jk_err_meff[T];
  for(int t=0; t<T; t++) {
    jk_err_mean[t] = 0.0;
    jk_err_meff[t] = 0.0;
    for(int i=0; i<Nr; i++) {
      resamp_mean[t][i] = 0.0;
      resamp_meff[t][i] = 0.0;
    }
  }
    
  //Compute resampled mean
  for(int i=0; i<Nr; i++) {
    for(int t=0; t<T; t++) {
      for(int n=0; n<N; n++) {
	if( n<i*B || n>=(i+1)*B ) {
	  resamp_mean[t][i] += array[n][t];
	}
      }
      resamp_mean[t][i] *= resamp_norm;
    }
  }

  //Compute resampled effective masses
  for(int i=0; i<Nr; i++) {
    resamp_meff[0][i] = 0.0;
    for(int t=1; t<T; t++) {
      resamp_meff[t][i] = -log(resamp_mean[t][i]/resamp_mean[t-1][i]);
    }
  }
  
  //Get jk error on mean and effective mass
  for(int t=0; t<T; t++) {
    for(int i=0; i<Nr; i++) {
      jk_err_meff[t] += pow(meff[t] - resamp_meff[t][i],2);
      jk_err_mean[t] += pow(mean[t] - resamp_mean[t][i],2);
    }
    jk_err_meff[t] = sqrt(coeff*jk_err_meff[t]);
    jk_err_mean[t] = sqrt(coeff*jk_err_mean[t]);
  }
  
  //Dump to disk

  // This is pion data
  if(T>3) {
    string mname("m_eff.dat"); //Name of the output file
    fstream outputFile;
    outputFile.open(mname, fstream::out);
    if(!outputFile.is_open()) {
      cout << "Error opening file: " << mname << endl;
      exit(0);
    }    
    cout << "Computing pion effective mass" << endl;
    //When dumping to disk, forget the 0th column which was just the data point entry.
    for(int t=1; t<T; t++) {
      outputFile << t-1 << setprecision(16) << " " << mean[t] << " " << jk_err_mean[t] << " " << meff[t] << " " << jk_err_meff[t] << endl;
    }  
    outputFile.close();
  }
  // This is wilson loop data
  else {
    string name("NN_loop.dat"); //Name of the output file
    fstream outputFile;
    outputFile.open(name, fstream::out);
    if(!outputFile.is_open()) {
      cout << "Error opening file: " << name << endl;
      exit(0);
    }
    cout << "Computing average NN loop" << endl;
    
    //When dumping to disk, forget the 0th column which was just the data point entry.
    for(int t=1; t<T-1; t++) {
      outputFile << WL << setprecision(16) << " " << mean[t] << " " << jk_err_mean[t] << endl;
    }  
  }  
}
