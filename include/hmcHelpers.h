#ifndef HMC_HELPERS
#define HMC_HELPERS

//2D Routines
//-------------------------------------------------------------------------------------
void forceU(double fU[LX][LY][2], Complex gauge[LX][LY][2], param_t p) {
  
  Complex plaq0;
  Complex plaq;
  zeroLat(fU);
  int xp1, xm1, yp1, ym1;
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++) {

      xp1 = (x+1)%LX;
      xm1 = (x-1+LX)%LX;
      yp1 = (y+1)%LY;
      ym1 = (y-1+LY)%LY;
      
      plaq0 = gauge[x][y][0]*gauge[xp1][y][1]*conj(gauge[x][yp1][0])*conj(gauge[x][y][1]);
      fU[x][y][0] += p.beta*imag(plaq0);
      
      plaq =  gauge[x][ym1][0]*gauge[xp1][ym1][1]*conj(gauge[x][y][0])*conj(gauge[x][ym1][1]);
      fU[x][y][0] -= p.beta*imag(plaq);

      plaq =  gauge[x][y][1]*conj(gauge[xm1][yp1][0])*conj(gauge[xm1][y][1])*gauge[xm1][y][0];
      fU[x][y][1] += p.beta*imag(plaq);

      //This plaquette was aleady computed. We want the conjugate.
      fU[x][y][1] -= p.beta*imag(plaq0);      
    }
}

//P_{k+1/2} = P_{k-1/2} - dtau * (fU + fD)
void update_mom(double fU[LX][LY][2], double fD[LX][LY][2], double mom[LX][LY][2], double dtau){

  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int mu=0; mu<2; mu++)
	mom[x][y][mu] -= (fU[x][y][mu] - fD[x][y][mu])*dtau;
}

//U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
void update_gauge(Complex gauge[LX][LY][2], double mom[LX][LY][2], double dtau){
  
  for(int x=0; x<LX; x++)
    for(int y=0; y<LY; y++)
      for(int mu=0; mu<2; mu++)
	gauge[x][y][mu] *= polar(1.0, mom[x][y][mu] * dtau);
}
//----------------------------------------------------------------------------------



//3D Routines
//----------------------------------------------------------------------------------
// All Wilson loops are computed clockwise.
void forceU(double fU[LX][LY][LZ][D], const Complex gauge[LX][LY][LZ][D], param_t p) {

  zeroLat(fU);

  Complex plaq, plaq0;
  double beta = p.beta;
  double betaz= p.betaz;
  
  int xp1, xm1, yp1, ym1, zp1, zm1;
  
  for(int x=0; x<LX; x++) {
    xp1 = (x+1)%LX;
    xm1 = (x-1+LX)%LX;
    for(int y=0; y<LY; y++) {
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

//P_{k+1/2} = P_{k-1/2} - dtau * (fU + fD)
void update_mom(double fU[LX][LY][LZ][D], double fD[LX][LY][2], double mom[LX][LY][LZ][D], double dtau, param_t p){
  
  int x,y,z,mu;
  
  //Always update from the 2D gauge fields
  for(x=0; x<LX; x++)
    for(y=0; y<LY; y++) 
      for(z=0; z<LZ; z++)
	for(mu=0; mu<2; mu++) {
	  mom[x][y][z][mu] -= fU[x][y][z][mu]*dtau;	  
	}
  
  //Update from the fermion field if dynamic
  if(p.dynamic == true) {
    for(x=0; x<LX; x++)
      for(y=0; y<LY; y++)
	for(mu=0; mu<2; mu++) {
	  mom[x][y][(LZ-1)/2][mu] += fD[x][y][mu]*dtau;
	}  
  }
}
      

//U_{k} = exp(i dtau P_{k-1/2}) * U_{k-1}
void update_gauge(Complex gauge[LX][LY][LZ][D], double mom[LX][LY][LZ][D], double dtau, param_t p){
  
  int x,y,z,mu;
  
  //Always update from the 2D gauge fields
  for(x=0; x<LX; x++)
    for(y=0; y<LY; y++)
      for(z=0; z<LZ; z++)
	for(mu=0; mu<2; mu++) {
	  gauge[x][y][z][mu] *= polar(1.0, mom[x][y][z][mu] * dtau);
	}
  
  //Update from the extra dimension if not z locked.
  if(p.lockedZ == false) {
    for(x=0; x<LX; x++)
      for(y=0; y<LY; y++) 
	for(z=0; z<LZ; z++)
	  for(mu=2; mu<D; mu++) {
	    gauge[x][y][z][mu] *= polar(1.0, mom[x][y][z][mu] * dtau);
	  }
  }
}
//------------------------------------------------------------------------------------
#endif
