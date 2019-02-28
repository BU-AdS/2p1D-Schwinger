#ifndef DOPHELPERS_H
#define DOPHELPERS_H

// Wilson stencil

// D_{W}(n,m) = (m_{0} + 2r)\delta(n,m)
//               - 1/2 Sum [(r-\sigma_{mu}) U_{n,\mu} \delta_{n,m-\hat{\mu}} +
//                          (r+\sigma_{mu}) U^{\dagger}_{m,\mu} \delta_{n,m+\hat{\mu}}]

// sigma_1 = | 0  1 |  sigma_2 = | 0 -i | sigma_3 = i*sigma_1*sigma_2 = | 1  0 |
//           | 1  0 |            | i  0 |                               | 0 -1 |

void Dpsi(Complex psi2[L][L][2], Complex psi1[L][L][2],
	  const Complex gauge[L][L][2], param_t p ){
  
  double m0 = p.m;
  double  r = 1.0;
  double constant = (2*r + m0);
  
  //Sum over 0,1 directions.
  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++) {

      //upper
      psi2[x][y][0] = constant * psi1[x][y][0] -
	
	0.5*(     gauge[x][y][0]          * ( r*psi1[(x+1)%L  ][y][0] - psi1[(x+1)%L  ][y][1]) +
	     conj(gauge[(x-1+L)%L][y][0]) * ( r*psi1[(x-1+L)%L][y][0] + psi1[(x-1+L)%L][y][1]) +
		  
		  gauge[x][y][1]          * ( r*psi1[x][(y+1)%L  ][0] + I*psi1[x][(y+1)%L  ][1]) +
	     conj(gauge[x][(y-1+L)%L][1]) * ( r*psi1[x][(y-1+L)%L][0] - I*psi1[x][(y-1+L)%L][1]));
      
      //lower
      psi2[x][y][1] = constant * psi1[x][y][1] -
	
	0.5*(     gauge[x][y][0]          * (-psi1[(x+1)%L][y  ][0] + r*psi1[(x+1)%L  ][y][1]) -
	     conj(gauge[(x-1+L)%L][y][0]) * (-psi1[(x-1+L)%L][y][0] - r*psi1[(x-1+L)%L][y][1]) +

		  gauge[x][y][1]          * (-I*psi1[x][(y+1)%L  ][0] + r*psi1[x][(y+1)%L  ][1]) -
	     conj(gauge[x][(y-1+L)%L][1]) * (-I*psi1[x][(y-1+L)%L][0] - r*psi1[x][(y-1+L)%L][1]));
      
      
    }
}

void Ddagpsi(Complex psi2[L][L][2], Complex psi1[L][L][2],
	     const Complex gauge[L][L][2], param_t p ){
  
}

void g5psi(Complex psi2[L][L][2], Complex psi1[L][L][2],
	   const Complex gauge[L][L][2], param_t p ){

  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++) {
      psi2[x][y][0] =  psi1[x][y][0];
      psi2[x][y][1] = -psi1[x][y][1];
    }
}

void g5Dpsi(Complex psi2[L][L][2], Complex psi1[L][L][2],
	    const Complex gauge[L][L][2], param_t p ){

  Complex temp[L][L][2];
  Dpsi(temp, psi1, gauge, p);
  g5psi(psi2, temp, gauge, p);
}

void DdagDpsi(Complex psi2[L][L][2], Complex  psi1[L][L][2],
	      const Complex gauge[L][L][2], param_t p) {

  //Hack for now
  Complex temp[L][L][2];
  Dpsi(temp, psi1, gauge, p);
  g5psi(psi2, temp, gauge, p);
  Dpsi(temp, psi2, gauge, p);
  g5psi(psi2, temp, gauge, p);
}



//===============================================================
// CG solutions to Apsi = b 
// see http://en.wikipedia.org/wiki/Conjugate_gradient_method
//===============================================================

int Ainv_psi(Complex psi[L][L][2], Complex b[L][L][2], Complex psi0[L][L][2],
	     const Complex gauge[L][L][2], param_t p) {

  int success = 0;

  Complex res[L][L][2], pvec[L][L][2], Apvec[L][L][2];
  double alpha, beta, denom;
  double rsq = 0, rsqNew = 0, bsqrt = 0.0;
  
  //Intialize  
  zeroField(res);
  zeroField(Apvec);
  zeroField(pvec);
  
  // Find norm of rhs.
  bsqrt = real(dotField(b,b));
  bsqrt = sqrt(bsqrt);
  
  // res = b  - A psi0, for now start with phi0 = 0
  copyField(res,b); 
  copyField(pvec, res);

  rsq = real(dotField(res,res));
  
  // Compute Ap.
  DdagDpsi(Apvec, pvec, gauge,p);
  
  // iterate till convergence
  int k;
  for(k = 0; k< p.maxIterCG; k++) {
    
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

  DdagDpsi(Apvec,  psi, gauge,p);
  axpy(-1.0, Apvec, b, res);
  
  //double truersq =  real(dotField(res,res));
  //printf("CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  return success;
  
}

#endif
