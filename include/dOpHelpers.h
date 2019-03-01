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

int Ainv_psi(Complex psi[L][L][2], Complex b[L][L][2], Complex psi0[L][L][2], const Complex gauge[L][L][2], param_t p) {

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

void Dpsi(Complex psi2[L][L], Complex psi1[L][L],
	  const Complex gauge[L][L][2], param_t p ){
  
  for(int x=0; x<L;x++)
    for(int y=0; y<L; y++){
      psi2[x][y] = p.m*psi1[x][y];
      // cout<<"psi1[x,y] = ["<< x <<","<<y <<" ]  "<< psi1[x][y] << endl;
      // cout<<"Dpsi2[x,y] = ["<< x <<","<<y <<" ]  "<< psi2[x][y] << endl <<endl;
    }
  
  double eta1;
  for(int x=0; x<L; x++) {
    eta1 =(1-2*(x%2));
    for(int y=0; y<L; y++) {      
      psi2[x][y] += - (gauge[x][y][0] * psi1[(x+1)%L][y]
		       - conj(gauge[ (x-1+L)%L ][y][0]) * psi1[ (x-1+L)%L ][y]);
      
      psi2[x][y] += - eta1*(gauge[x][y][1]*psi1[x][ (y+1)%L ]
			    - conj(gauge[x][ (y-1+L)%L ][1])*psi1[x][ (y-1+L)%L ]);
      
      // cout<<"Dpsi2[x,y] = ["<< x <<","<<y <<" ]  "<< psi2[x][y] << endl <<endl;
    }
  }
}

void Ddagpsi(Complex psi2[L][L], Complex  psi1[L][L],
	     const Complex gauge[L][L][2], param_t p ) {
  
  // For a 2D square lattice, the stencil is:
  //   1 |  0 -eta1  0 |
  //   - | +1    0  -1 |  , where eta0 = 1, eta1 = (-)^x = 1 - 2(x%L)
  //   2 |  0 +eta1  0 |

  for(int x=0; x<L; x++)
    for(int y=0; y<L; y++)
      psi2[x][y] = p.m*psi1[x][y];
  
  double eta1;
  for(int x=0; x<L; x++) {
    eta1 =(1-2*(x%2));    
    for(int y=0; y<L; y++) {      
      psi2[x][y] += (gauge[x][y][0] * psi1[(x+1)%L][y]
		     - conj(gauge[(x-1 +L)%L][y][0]) * psi1[(x-1 + L)%L][y]);
      
      psi2[x][y] += eta1*(gauge[x][y][1]*psi1[x][(y+1)%L]
			  - conj(gauge[x][(y-1+L)%L][1])*psi1[x][(y-1+L)%L]);
    }
  }
}

//=======================//
// Note: Ddag D = D Ddag //
///======================//

void DdagDpsi(Complex psi2[L][L], Complex  psi1[L][L],
	      const Complex gauge[L][L][2], param_t p ) {
  
  Complex psitmp[L][L];
  zeroField(psitmp);
  Dpsi(psitmp, psi1, gauge, p);
  Ddagpsi(psi2, psitmp, gauge, p);
}

/*===================== 
  CG solutions to Apsi = b 
  see http://en.wikipedia.org/wiki/Conjugate_gradient_method
  ======================*/

int Ainv_psi(Complex psi[L][L], Complex b[L][L], Complex psi0[L][L], const Complex gauge[L][L][2], param_t p) {

  Complex res[L][L] , pvec[L][L], Apvec[L][L];
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

  rsq =  real(dotField(res,res));
  
  // cout << "# Enter Ainv_p  bsqrt =  " << bsqrt << "  rsq   = " << rsq << endl;
  
  // Compute Ap.
  DdagDpsi(Apvec, pvec, gauge, p);
  
  // iterate till convergence
  int k;
  int success = 0;
  for(k = 0; k< p.maxIterCG; k++) {
    
    denom = real(dotField(pvec,Apvec));
    alpha = rsq/denom;
    
    for(int x=0; x<L; x++)
      for(int y=0; y<L; y++) {
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
  
  //double truersq =  real(dotField(res,res));
  //printf("CG: Converged iter = %d, rsq = %e, truersq = %e\n",k,rsq,truersq);
  return success;
}

void TestCG(Complex gauge[L][L][2], param_t p) {
  
  Complex b[L][L];
  gaussComplex_F(b, p);
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







#endif
