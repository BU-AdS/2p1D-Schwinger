#include "ran2s.h"

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)

static long idum=987654321, idum2=123456789;
static long iy=0;
static long iv[NTAB];

void
sran2(long iseed)
{
  long k;
  int j;

  //if (*idum <= 0) {              /*Initialize.*/
  if (iseed == 0) iseed = 1;     /*Be sure to prevent idum = 0.*/
  else if (iseed < 0) iseed = -iseed;
  idum = idum2 = iseed;
  for(j=NTAB+7; j>=0; j--) {  /*Load the shuffle table (after 8 warm-ups).*/
    k = idum/IQ1;
    idum = IA1*(idum-k*IQ1) - k*IR1;
    if (idum < 0) idum += IM1;
    if (j < NTAB) iv[j] = idum;
  }
  iy = iv[0];
  //}
}

double
ran2(void)
/*********************
 Long period (> 2x10^18 ) random number generator of L'Ecuyer with
 Bays-Durham shuffle and added safeguards. Returns a uniform random
 deviate between 0.0 and 1.0 (exclusive of the endpoint values). Call
 with idum a negative integer to initialize; thereafter, do not alter
 idum between successive deviates in a sequence.
*********************/
{
  long k;
  int j;

  k = idum/IQ1;                  /* Start here when not initializing.*/
  /* Compute idum=(IA1*idum) % IM1 without overflows by Schrage's method.*/
  idum = IA1*(idum-k*IQ1) - k*IR1;
  if (idum < 0) idum += IM1;

  k = idum2/IQ2;
  /* Compute idum2=(IA2*idum) % IM2 likewise.*/
  idum2 = IA2*(idum2-k*IQ2) - k*IR2;
  if (idum2 < 0) idum2 += IM2;

  j = iy/NDIV;                      /*Will be in the range 0..NTAB-1.*/
  /* Here idum is shuffled, idum and idum2 are combined to generate output.*/
  iy = iv[j]-idum2;
  iv[j] = idum;
  if (iy < 1) iy += IMM1;
  return AM*iy;
}
