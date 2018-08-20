/*=====================================================================================

  Thu Aug 02 14:15:21 EDT Dean Howarth
  
  ARPACK interafce for 2D compact U(1) theory.

  ====================================================================================*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <complex>
#include <omp.h>

using namespace std;

#include "ran2s.h"

void arpackErrorHelpNAUPD();
void arpackErrorHelpNEUPD();

static void mergeAbs(double *sort1, int *idx1, int n1, double *sort2,
		     int *idx2, int n2, bool inverse) {
  int i1=0, i2=0;
  int *ord;
  double *result;
  
  ord    = (int *)    malloc(sizeof(int)   *(n1+n2)); 
  result = (double *) malloc(sizeof(double)*(n1+n2)); 
  
  for(int i=0; i<(n1+n2); i++) {
    if((fabs(sort1[i1]) >= fabs(sort2[i2])) != inverse) { //LOGICAL XOR
      result[i] = sort1[i1];
      ord[i] = idx1[i1];
      i1++;
    } else {
      result[i] = sort2[i2];
      ord[i] = idx2[i2];
      i2++;
    }
    
    if(i1 == n1) {
      for(int j=i+1; j<(n1+n2); j++,i2++) {
	result[j] = sort2[i2];
	ord[j] = idx2[i2];
      }
      i = n1+n2;
    } else if (i2 == n2) {
      for(int j=i+1; j<(n1+n2); j++,i1++) {
	result[j] = sort1[i1];
	ord[j] = idx1[i1];
      }
      i = i1+i2;
    }
  }  
  for(int i=0;i<n1;i++) {
    idx1[i] = ord[i];
    sort1[i] = result[i];
  }
  
  for(int i=0;i<n2;i++) {
    idx2[i] = ord[i+n1];
    sort2[i] = result[i+n1];
  }  
  free (ord);
  free (result);
}
  
static void sortAbs(double *unsorted, int n, bool inverse, int *idx) {
  
  if (n <= 1)
    return;
  
  int n1,n2;
  
  n1 = n>>1;
  n2 = n-n1;
  
  double *unsort1 = unsorted;
  double *unsort2 = (double *)((char*)unsorted + n1*sizeof(double));
  int *idx1 = idx;
  int *idx2 = (int *)((char*)idx + n1*sizeof(int));
  
  sortAbs(unsort1, n1, inverse, idx1);
  sortAbs(unsort2, n2, inverse, idx2);
  
  mergeAbs(unsort1, idx1, n1, unsort2, idx2, n2, inverse);
}



#define ARPACK(s) s ## _

#ifdef __cplusplus
extern "C" {
#endif

  /**
   *  Interface functions to the external ARPACK library. These functions utilize 
   *  ARPACK's implemntation of the Implicitly Restarted Arnoldi Method to compute a 
   *  number of eigenvectors/eigenvalues with user specified features, such as those 
   *  with small real part, small magnitude etc. Parallel (OMP/MPI) versions
   *  are also supported.
   */
  
  
  //Serial, single prec complex eigenvectors
  extern int ARPACK(cnaupd) (int *ido, char *bmat, int *n, char *which, int *nev,
			     float *tol, std::complex<float> *resid, int *ncv,
			     std::complex<float> *v, int *ldv, int *iparam, int *ipntr,
			     std::complex<float> *workd, std::complex<float> *workl,
			     int *lworkl, float *rwork, int *info, int bmat_size,
			     int spectrum_size);
  
  //Serial, double prec complex eigenvectors
  extern int ARPACK(znaupd)(int *ido, char *bmat, int *n, char *which, int *nev,
			    double *tol, std::complex<double> *resid, int *ncv,
			    std::complex<double> *v, int *ldv, int *iparam, int *ipntr,
			    std::complex<double> *workd, std::complex<double> *workl, 
			    int *lworkl, double *rwork, int *info, int bmat_size,
			    int spectrum_size);
  
  //Serial, single prec complex eigenvalues
  extern int ARPACK(cneupd) (int *comp_evecs, char *howmany, int *select,
			     std::complex<float> *evals, std::complex<float> *v,
			     int *ldv, std::complex<float> *sigma,
			     std::complex<float> *workev, char *bmat, int *n,
			     char *which, int *nev, float *tol,
			     std::complex<float> *resid, int *ncv,
			     std::complex<float> *v1, int *ldv1, int *iparam,
			     int *ipntr, std::complex<float> *workd,
			     std::complex<float> *workl, int *lworkl,
			     float *rwork, int *info, int howmany_size, int bmat_size,
			     int spectrum_size);			
  
  //Serial, double prec complex eigenvalues
  extern int ARPACK(zneupd) (int *comp_evecs, char *howmany, int *select,
			     std::complex<double> *evals, std::complex<double> *v,
			     int *ldv, std::complex<double> *sigma,
			     std::complex<double> *workev, char *bmat, int *n,
			     char *which, int *nev, double *tol,
			     std::complex<double> *resid, int *ncv,
			     std::complex<double> *v1, int *ldv1, int *iparam,
			     int *ipntr, std::complex<double> *workd,
			     std::complex<double> *workl, int *lworkl,
			     double *rwork, int *info, int howmany_size, int bmat_size,
			     int spectrum_size);
  
  extern int ARPACK(mcinitdebug)(int*,int*,int*,int*,int*,int*,int*,int*);
    
  //ARPACK initlog and finilog routines for printing the ARPACK log  
  extern int ARPACK(initlog) (int*, char*, int);
  extern int ARPACK(finilog) (int*);
  
#ifdef __cplusplus
}
#endif


void arpackErrorHelpNAUPD() {
  printf("\nError help NAUPD\n\n");
  printf("INFO Integer.  (INPUT/OUTPUT)\n");
  printf("     If INFO .EQ. 0, a randomly initial residual vector is used.\n");
  printf("     If INFO .NE. 0, RESID contains the initial residual vector,\n");
  printf("                        possibly from a previous run.\n");
  printf("     Error flag on output.\n");
  printf("     =  0: Normal exit.\n");
  printf("     =  1: Maximum number of iterations taken.\n");
  printf("        All possible eigenvalues of OP has been found. IPARAM(5)\n");
  printf("        returns the number of wanted converged Ritz values.\n");
  printf("     =  2: No longer an informational error. Deprecated starting\n");
  printf("        with release 2 of ARPACK.\n");
  printf("     =  3: No shifts could be applied during a cycle of the\n");
  printf("        Implicitly restarted Arnoldi iteration. One possibility\n");
  printf("        is to increase the size of NCV relative to NEV.\n");
  printf("        See remark 4 below.\n");
  printf("     = -1: N must be positive.\n");
  printf("     = -2: NEV must be positive.\n");
  printf("     = -3: NCV-NEV >= 1 and less than or equal to N.\n");
  printf("     = -4: The maximum number of Arnoldi update iteration\n");
  printf("        must be greater than zero.\n");
  printf("     = -5: WHICH must be 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'\n");
  printf("     = -6: BMAT must be one of 'I' or 'G'.\n");
  printf("     = -7: Length of private work array is not sufficient.\n");
  printf("     = -8: Error return from LAPACK eigenvalue calculation;\n");
  printf("     = -9: Starting vector is zero.\n");
  printf("     = -10: IPARAM(7) must be 1,2,3.\n");
  printf("     = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.\n");
  printf("     = -12: IPARAM(1) must be equal to 0 or 1.\n");
  printf("     = -9999: Could not build an Arnoldi factorization.\n");
  printf("        User input error highly likely.  Please\n");
  printf("        check actual array dimensions and layout.\n");
  printf("        IPARAM(5) returns the size of the current Arnoldi\n");
  printf("        factorization.\n");
}

void arpackErrorHelpNEUPD() {
  printf("\nError help NEUPD\n\n");
  printf("INFO Integer.  (OUTPUT)\n");
  printf("     Error flag on output.\n");
  printf("     =  0: Normal exit.\n");
  printf("     =  1: The Schur form computed by LAPACK routine csheqr\n");
  printf("        could not be reordered by LAPACK routine ztrsen.\n");
  printf("        Re-enter subroutine zneupd with IPARAM(5)=NCV and\n");
  printf("        increase the size of the array D to have\n");
  printf("        dimension at least dimension NCV and allocate at\n");
  printf("        least NCV\n");
  printf("        columns for Z. NOTE: Not necessary if Z and V share\n");
  printf("        the same space. Please notify the authors if this\n");
  printf("        error occurs.\n");
  printf("     = -1: N must be positive.\n");
  printf("     = -2: NEV must be positive.\n");
  printf("     = -3: NCV-NEV >= 1 and less than or equal to N.\n");
  printf("     = -5: WHICH must be 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'\n");
  printf("     = -6: BMAT must be one of 'I' or 'G'.\n");
  printf("     = -7: Length of private work WORKL array is inufficient.\n");
  printf("     = -8: Error return from LAPACK eigenvalue calculation.\n");
  printf("        This should never happened.\n");
  printf("     = -9: Error return from calculation of eigenvectors.\n");
  printf("        Informational error from LAPACK routine ztrevc.\n");
  printf("     = -10: IPARAM(7) must be 1,2,3\n");
  printf("     = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.\n");
  printf("     = -12: HOWMNY = 'S' not yet implemented\n");
  printf("     = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.\n");
  printf("     = -14: ZNAUPD did not find any eigenvalues to sufficient\n");
  printf("        accuracy.\n");
  printf("     = -15: ZNEUPD got a different count of the number of\n");
  printf("        converged Ritz values than ZNAUPD got. This\n");
  printf("        indicates the user probably made an error in\n");
  printf("        passing data from ZNAUPD to ZNEUPD or that the\n");
  printf("        data was modified before entering ZNEUPD\n");
}

