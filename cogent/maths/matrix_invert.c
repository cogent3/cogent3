/***************************************************
*
* routine for inverting a matrix, extracted from
* Ziheng Yang's source files.
*
****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#define PYCOGENT_VERSION "1.4"

/*****
*
* taking inv of matrix
*
*
******/

int matinv (double x[], int n, int m, double space[])
{
/* x[n*m]  ... m>=n
*/
   register int i,j,k;
   int *irow=(int*) space;
   double ee=1.0e-30, t,t1,xmax;
   double det=1.0;

   for (i = 0; i < n; i++)  {
      xmax = 0.;
      for (j= i; j < n; j++)
         if (xmax < fabs(x[j*m+i])) { xmax = fabs(x[j*m+i]); irow[i]=j; }
      det *= xmax;
      if (xmax < ee)   {
         /* printf("\nDet = %.4e close to zero at %3d!\t\n", xmax,i+1); */
         return(-1);
      }
      if (irow[i] != i) {
         for (j = 0; j < m; j++) {
            t = x[i*m+j];
            x[i*m+j] = x[irow[i]*m+j];
            x[irow[i]*m+j] = t;
         }
      }
      t = 1./x[i*m+i];
      for (j = 0; j < n; j++) {
         if (j == i) continue;
         t1 = t*x[j*m+i];
         for(k = 0; k < m; k++)  x[j*m+k] -= t1*x[i*m+k];
         x[j*m+i] = -t1;
      }
      for(j = 0; j < m; j++)   x[i*m+j] *= t;
      x[i*m+i] = t;
   }                            /* i  */
   for (i=n-1; i>=0; i--) {
      if (irow[i] == i) continue;
      for(j = 0; j < n; j++)  {
         t = x[j*m+i];
         x[j*m+i] = x[j*m + irow[i]];
         x[j*m + irow[i]] = t;
      }
   }
   return (0);
}

