/**
 * @file MathUtil.cxx
 * @brief Implementation of MathUtil class.
 * @author Analia Cillis
 * $Header: /nfs/slac/g/glast/ground/cvs/Likelihood/src/likelihood/MathUtil.cxx,v 1.1 2006/09/14 20:39:40 peachey Exp $
*/
#include "MathUtil.h"
#include <math.h>
#include <stdexcept>

#define RAND_ITMAX      10000//before was 1000
#define RAND_EPS        3.0e-7

double  MathUtil::gammln(double x){
   double tmp, sum;
   static double  cof[6] =                                        /* static data! */
         {76.18009173, -86.50532033, 24.01409822,
         -1.231739516, 0.120858003e-2, -0.536382e-5};
   int j;
 
   x -= 1.0;
   tmp = x + 5.5;
   tmp -= (x + 0.5) * log(tmp);
   sum = 1.0;
   for (j = 0; j <= 5; j++) {
      x += 1.0;
      sum += cof[j] / x;
   }
   return (-tmp + log(2.50662827465 * sum));
}

void  MathUtil::gser(double *gamser, double a, double x, double *gln) {
   int n;
   double sum, del, ap; 
   *gln = gammln(a);
   if (x <= 0.0) {
      if (x < 0.0)
          return;
      *gamser = 0.0;
      return;
   } else {
       ap = a;
       del = sum = 1.0 / a;
       for (n = 1; n <= RAND_ITMAX; n++) {
          ap += 1.0;
          del *= x / ap;
          sum += del;
          if (fabs(del) < fabs(sum) * RAND_EPS) {
             *gamser = sum * exp(-x + a * log(x) - (*gln));
             return;
          }
       }
       throw std::runtime_error("a too large, RAND_ITMAX too small in routine GSER");
       return;
   }
 }
 
void  MathUtil::gcf(double *gammcf, double a, double x, double *gln) {
   int n;
   double gold = 0.0, g, fac = 1.0, b1 = 1.0;
   double b0 = 0.0, anf, ana, an, a1, a0 = 1.0;
 
   *gln = gammln(a);
   a1 = x;
   for (n = 1; n <= RAND_ITMAX; n++) {
      an = (double) n;
      ana = an - a;
      a0 = (a1 + a0 * ana) * fac;
      b0 = (b1 + b0 * ana) * fac;
      anf = an * fac;
      a1 = x * a0 + anf * a1;
      b1 = x * b0 + anf * b1;
      if (a1) {
         fac = 1.0 / a1;
         g = b1 * fac;
         if (fabs((g - gold) / g) < RAND_EPS) {
            *gammcf = exp(-x + a * log(x) - (*gln)) * g;
            return;
         }
         gold = g;
      }
   }
 }
 
double  MathUtil::gammq(double a, double x) {
   double gamser, gammcf, gln;
 
   if (x < 0.0 || a <= 0.0)
       return (0.0);
   if (x < (a + 1.0)) {
       gser(&gamser, a, x, &gln);
       return 1.0 - gamser;
   } else {
      gcf(&gammcf, a, x, &gln);
      return (gammcf);
   }
 }
 
double MathUtil::poissonSig(double x, double mean) {
   double significance;
      
   if (mean < x) {
       significance= (1-MathUtil::gammq(x,mean)) + MathUtil::gammq(x,2*x-mean);
   } else {
      if ( (2*x-mean) > 0) {
       significance= MathUtil::gammq(x,mean)+ (1-MathUtil::gammq(x,2*x-mean));
       } else {
          significance= MathUtil::gammq(x,mean);
       }
   }
  return significance;
}
