#include <math.h>

extern double ran2(void);

double
rang(void)
{
 static long iset=0;
 static double gset;
 double v1,v2,rsq,retval,fac;

 if(iset==0) {
   do {
     v1 = 2.0*ran2() - 1.0;
     v2 = 2.0*ran2() - 1.0;
     rsq = v1*v1 + v2*v2;
   } while((rsq>=1.0)||(rsq==0.0));
   fac = sqrt(-2.0*log(rsq)/rsq);
   gset = v1*fac;
   iset = 1;
   retval = v2*fac;
 } else {
   iset = 0;
   retval = gset;
 }
 return(retval);
}
