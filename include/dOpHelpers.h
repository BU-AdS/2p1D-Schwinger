#ifndef DOPHELPERS_H
#define DOPHELPERS_H

// Wilson stencil
// D_{W}(n,m) = (m_{0} + 2r)\delta(n,m)
//               - 1/2 Sum [(r-\sigma_{mu}) U_{n,\mu} \delta_{n,m-\hat{\mu}} +
//                          (r+\sigma_{mu}) U^{\dagger}_{m,\mu} \delta_{n,m+\hat{\mu}}]
// sigma_1 = | 0  1 |  sigma_2 = | 0 -i | sigma_3 = i*sigma_1*sigma_2 = | 1  0 |
//           | 1  0 |            | i  0 |                               | 0 -1 |

void Dpsi(Complex psi2[LX][LY][2], Complex psi1[LX][LY][2],
	  const Complex gauge[LX][LY][2], param_t p ){
  
  double m0 = p.m;
  double  r = 1.0;
  double constant = (2*r + m0);
  int xp1, xm1, yp1, ym1;
  
  //Sum over 0,1 directions.
  for(int x=0; x<LX; x++) {
    xp1 = (x+1)%LX;
    xm1 = (x-1+LX)%LX;
    for(int y=0; y<LY; y++) {
      yp1 = (y+1)%LY;
      ym1 = (y-1+LY)%LY;
      
      //upper
      psi2[x][y][0] = constant * psi1[x][y][0] -
	
	0.5*(     gauge[x][y][0]    * (r*psi1[xp1][y][0] - psi1[xp1][y][1]) +
	     conj(gauge[xm1][y][0]) * (r*psi1[xm1][y][0] + psi1[xm1][y][1]) +
		  
		  gauge[x][y][1]    * (r*psi1[x][yp1][0] + I*psi1[x][yp1][1]) +
	     conj(gauge[x][ym1][1]) * (r*psi1[x][ym1][0] - I*psi1[x][ym1][1]));
      
      //lower
      psi2[x][y][1] = constant * psi1[x][y][1] -
	
	0.5*(     gauge[x][y][0]    * (-psi1[xp1][y][0] + r*psi1[xp1][y][1]) -
	     conj(gauge[xm1][y][0]) * (-psi1[xm1][y][0] - r*psi1[xm1][y][1]) +

		  gauge[x][y][1]    * (-I*psi1[x][yp1][0] + r*psi1[x][yp1][1]) -
	     conj(gauge[x][ym1][1]) * (-I*psi1[x][ym1][0] - r*psi1[x][ym1][1]));
      
      
    }
  }
}

void Ddagpsi(Complex psi2[LX][LY][2], Complex psi1[LX][LY][2],
	     const Complex gauge[LX][LY][2], param_t p ){
  
}

void g5psi(Complex psi2[LX][LY][2], Complex psi1[LX][LY][2]){
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++) {
      psi2[x][y][0] =  psi1[x][y][0];
      psi2[x][y][1] = -psi1[x][y][1];
    }
}

void g5psi(Complex psi1[LX][LY][2]){
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++) {
      psi1[x][y][1] *= -1.0;
    }
}


void g5Dpsi(Complex psi2[LX][LY][2], Complex psi1[LX][LY][2],
	    const Complex gauge[LX][LY][2], param_t p ){

  Complex temp[LX][LY][2];
  Dpsi(temp, psi1, gauge, p);
  g5psi(psi2, temp);
}

void DdagDpsi(Complex psi2[LX][LY][2], Complex psi1[LX][LY][2],
	      const Complex gauge[LX][LY][2], param_t p) {
  
  //Hack for now
  Complex temp[LX][LY][2];
  Dpsi(temp, psi1, gauge, p);
  g5psi(temp);
  Dpsi(psi2, temp, gauge, p);
  g5psi(psi2);
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

void Dpsi(Complex psi2[LX][LY], Complex psi1[LX][LY],
	  const Complex gauge[LX][LY][2], param_t p ){

  int xp1, xm1, yp1, ym1;
  double eta1;  
  
  for(int x=0; x<LX; x++) {
    eta1 =(1-2*(x%2));
    xp1 = (x+1)%LX;
    xm1 = (x-1+LX)%LX;    
    for(int y=0; y<LY; y++) {
      yp1 = (y+1)%LY;
      ym1 = (y-1+LY)%LY;

      psi2[x][y] = p.m*psi1[x][y];
      
      psi2[x][y] += - (gauge[x][y][0] * psi1[xp1][y]
		       - conj(gauge[xm1][y][0]) * psi1[xm1][y]);
      
      psi2[x][y] += - eta1*(gauge[x][y][1]*psi1[x][yp1]
			    - conj(gauge[x][ym1][1])*psi1[x][ym1]);
      
    }
  }
}

void Ddagpsi(Complex psi2[LX][LY], Complex  psi1[LX][LY],
	     const Complex gauge[LX][LY][2], param_t p ) {
  
  // For a 2D square lattice, the stencil is:
  //   1 |  0 -eta1  0 |
  //   - | +1    0  -1 |  , where eta0 = 1, eta1 = (-)^x = 1 - 2(x%L)
  //   2 |  0 +eta1  0 |

  int xp1, xm1, yp1, ym1;
  double eta1;

  for(int x=0; x<LX; x++) {
    eta1 =(1-2*(x%2));
    xp1 = (x+1)%LX;
    xm1 = (x-1+LX)%LX;    
    for(int y=0; y<LY; y++) {
      yp1 = (y+1)%LY;
      ym1 = (y-1+LY)%LY;

      psi2[x][y] = p.m*psi1[x][y];
      
      psi2[x][y] += (gauge[x][y][0] * psi1[xp1][y]
		     - conj(gauge[xm1][y][0]) * psi1[xm1][y]);
      
      psi2[x][y] += eta1*(gauge[x][y][1]*psi1[x][yp1]
			  - conj(gauge[x][ym1][1])*psi1[x][ym1]);
    }
  }
}

//=======================//
// Note: Ddag D = D Ddag //
///======================//

void DdagDpsi(Complex psi2[LX][LY], Complex  psi1[LX][LY],
	      const Complex gauge[LX][LY][2], param_t p ) {
  
  Complex psitmp[LX][LY];
  zeroField(psitmp);
  Dpsi(psitmp, psi1, gauge, p);
  Ddagpsi(psi2, psitmp, gauge, p);
}

#include "arpack_interface_wilson.h"

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
  double rsq = 0, rsqNew = 0, bsqrt = 0.0;
  bool deflating = false;
  
  //Intialize
  zeroField(res);
  zeroField(Ap);
  zeroField(p);
  zeroField(x);
  
  // Find norm of rhs.
  bsqrt = sqrt(norm2(b));
  
  // res = b - A*x0
  if (norm2(x0) != 0.0) {

    //Solve the deflated system.
    deflating = true;
    DdagDpsi(tmp, x0, gauge, param);
    axpby(1.0, b, -1.0, tmp, res);  
    //cout << "using initial guess, |x0| = " << sqrt(norm2(x0))
    //<< ", |b| = " << bsqrt
    //<< ", |res| = " << sqrt(norm2(res)) << endl;
  } else {
    //Solve the given system
    copyField(res, b);  
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
    if (sqrt(rsqNew) < sqrt(param.eps)*bsqrt) {
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
    // Add back the deflated guess 
    axpy(1.0, x0, x, res);
    copyField(x, res);
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
