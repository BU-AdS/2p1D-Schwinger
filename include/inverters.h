#ifndef INVERTERS_H
#define INVERTERS_H

#include "dOpHelpers.h"

//===============================================================
// CG solutions to Apsi = b 
// see http://en.wikipedia.org/wiki/Conjugate_gradient_method
//===============================================================
//       x: The solution
//       b: The RHS vector
//      x0: An initial guess
//   gauge: The gauge field defining the operator
//   param: The parameter container

// Wilson g5Dg5D matrix inverter
//---------------------------------------------------------------
int Ainvpsi(Complex x[LX][LY][2], Complex b[LX][LY][2], Complex x0[LX][LY][2],
	    const Complex gauge[LX][LY][2], param_t param) {
  
  int success = 0;
  
  Complex res[LX][LY][2], p[LX][LY][2], Ap[LX][LY][2], tmp[LX][LY][2];
  double alpha, beta, denom;
  double rsq = 0, rsqNew = 0, bsqrt = 0.0, bnorm = 0.0;
  bool deflating = false;
  
  //Intialize
  zeroField(res);
  zeroField(Ap);
  zeroField(p);
  zeroField(x);
  
  // Find norm of rhs.
  bnorm = norm2(b);
  bsqrt = sqrt(bnorm);
  if(bsqrt == 0 || bsqrt != bsqrt) {
    printVector(b);
    cout << "Error in Wilson Ainvpsi: inverting on zero source... or nan!" << endl;
    exit(0);
  }
  copyField(res, b);
  
  // res = b - A*x0
  if (norm2(x0) != 0.0) {
    
    //Solve the deflated system.
    deflating = true;
    DdagDpsi(tmp, x0, gauge, param);    
    axpy(-1.0, tmp, res);
    
    cout << "using initial guess, |x0| = " << sqrt(norm2(x0))
	 << ", |b| = " << bsqrt
	 << ", |res| = " << sqrt(norm2(res)) << endl;
  }
  
  copyField(p, res);
  rsq = norm2(res);
  
  // Iterate until convergence
  int k;
  for (k=0; k<param.maxIterCG; k++) {

    // Compute Ap.
    DdagDpsi(Ap, p, gauge, param);
    
    denom = real(dotField(p, Ap));
    alpha = rsq/denom;
    
    axpy( alpha, p,  x);
    axpy(-alpha, Ap, res);
    
    // Exit if new residual is small enough
    rsqNew = norm2(res);
    //printf("CG iter %d, rsq = %g\n", k+1, rsqNew);
    if (rsqNew < param.eps*bnorm) {
      rsq = rsqNew;
      break;
    }
    
    // Update vec using new residual
    beta = rsqNew/rsq;
    rsq = rsqNew;
    
    axpy(beta, p, res, p);
    
  } // End loop over k
  
  if(k == param.maxIterCG) {
    // Failed convergence 
    printf("CG: Failed to converge iter = %d, rsq = %.16e\n", k+1, rsq); 
    success = 0; 
  } else {
    // Convergence 
    success = 1; 
  }
  
  if(deflating) {
    // x contains the solution to the deflated system b - A*x0.
    // We must add back the exact part
    axpy(1.0, x0, x);
    // x now contains the solution to the RHS b.
  }
  DdagDpsi(tmp, x, gauge, param);
  axpy(-1.0, tmp, b, res);
  
  //double truersq = real(dotField(res, res));
  //printf("CG: Converged iter = %d, rsq = %.16e, truersq = %.16e\n", k+1, rsq, truersq/(bsqrt*bsqrt));
  return success;
  
}

//Staggered
int Ainvpsi(Complex psi[LX][LY], Complex b[LX][LY], Complex psi0[LX][LY], const Complex gauge[LX][LY][2], param_t p) {

  int success = 0;
  
  Complex res[LX][LY] , pvec[LX][LY], Apvec[LX][LY];
  double alpha, beta, denom ;
  double rsq = 0, rsqNew = 0, bsqrt = 0.0;
  
  //Intialize  
  zeroField(res);
  zeroField(Apvec);  
  zeroField(pvec);
  
  // Find norm of rhs.
  bsqrt = real(dotField(b,b));
  bsqrt = sqrt(bsqrt);
  
  copyField(res, b); // res = b  - A psi0, for now start with phi0 = 0
  copyField(pvec, res);

  rsq = real(dotField(res,res));
  
  // Compute Ap.
  DdagDpsi(Apvec, pvec, gauge, p);

  // iterate till convergence
  int k;
  for (k=0; k<p.maxIterCG; k++) {
    
    denom = real(dotField(pvec,Apvec));
    alpha = rsq/denom;

    axpy( alpha, pvec, psi);
    axpy(-alpha, Apvec, res);
    
    // Exit if new residual is small enough
    rsqNew = real(dotField(res,res));
    
    if (sqrt(rsqNew) < p.eps*bsqrt) {
      //printf("Final rsq = %g\n", rsqNew);
      break;
    }
    
    // Update vec using new residual
    beta = rsqNew / rsq;
    rsq = rsqNew;
    
    axpy(beta, pvec, res, pvec);
    
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

  DdagDpsi(Apvec, psi, gauge,p);
  axpy(-1.0, Apvec, b, res);
  
  //double truersq =  real(dotField(res,res));
  //printf("CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  return success;
}

#endif
