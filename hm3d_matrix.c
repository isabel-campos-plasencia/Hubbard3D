/***********************************************************/
/* hm3d_matrix.c : Simulation of the Hubbard model in d=3  */
/*        Operations with matrix and pointers              */
/*                                                         */
/*  Index convention from NR: [1,...,V]                    */
/*                                                         */
/* I.Campos (December.00)  version 1.0                    */
/***********************************************************/


#include "hm3d.h"


void matrix_mult(int n, double **aa, double **bb, double **cc)
{
  int i,j,k,j0,i0,k0;
  double a00,a01,a10,a11;
  double v1[V + 1],v2[V + 1];   


  for(j=1; j<=n; j+=2)
    {
      for(i=1; i<=n; i++) 
	{
	  v1[i] = bb[i][j];
	  v2[i] = bb[i][j+1];
	}
      for(i=1; i<=n; i+=2)
	{
	  a00 = a01 = a10 = a11 = 0.0;
	
	  for(k=1;k<=n;k++)
	    {
	      a00 += aa[i][k]*v1[k];
	      a01 += aa[i][k]*v2[k];
	      a10 += aa[i+1][k]*v1[k];
	      a11 += aa[i+1][k]*v2[k];
	
	    }

	      
	  cc[i][j] = a00;
	  cc[i][j+1] = a01;
	  cc[i+1][j] = a10;
	  cc[i+1][j+1] = a11;
	

	}
    }
}
  



