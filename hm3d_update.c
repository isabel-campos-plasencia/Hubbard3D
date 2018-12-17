/**********************************************************/
/* hm3d_update.c : Simulation of the Hubbard model in d=3 */
/*          Update related functions                      */
/*                                                        */
/*    I.Campos (November.00)  version 2.0                */
/*                                                        */
/**********************************************************/


#include "hm3d.h"

extern signed char spin[V+1][Ns_p1];               /* Ising model */

extern float t,u,Beta,delta_tau,lambda;

extern float prob0, prob1;

extern int good;

void HeatBath(int slice, int site, double **gup, double **gdown)
{

  int i,j,k,r,s;
  double **a,**b,**cerog_up, **cerog_down;
  float R_up,R_down,R,deltaii_up,deltaii_down,F_up,F_down; 
  
  if (spin[site][slice] == 1)
    {
      deltaii_up = prob0;
      deltaii_down = prob1;
    }
  if (spin[site][slice] == - 1)
    {
      deltaii_up = prob1;
      deltaii_down = prob0;
    }
  
  
  R_up = fabs(1.0 + (1.0 - gup[site][site])*deltaii_up);
  R_down = fabs(1.0 + (1.0 - gdown[site][site])*deltaii_down);
  
  R = (R_up * R_down) /(1.0 + (R_up*R_down));  

#ifdef DEBUG 
  printf("site = %d\n",site);
  printf("G up = %g\n",gup[site][site]);
  printf("G down = %g\n",gdown[site][site]);
  printf("R_up = %g\n",R_up);
  printf("R_down = %g\n",R_down);                      
  printf("Heatbath ratio = %g\n\n",R); 
#endif   
  

  /****************** update of the spin + Green Function  **************/


  
  if(FRANDOM < R)            
    {

      spin[site][slice] = - spin[site][slice];


      good++;
  
      F_up = deltaii_up/R_up;
      F_down = deltaii_down/R_down;
  
      a = dmatrix(1,V,1,V);
      b = dmatrix(1,V,1,V);


      for(i=1;i<=V;i++) 	
	for(j=1;j<=V;j++) 	
	  {
	    a[i][j] = gup[i][j] + (gup[i][site]*gup[site][j])*F_up;
	    b[i][j] = gdown[i][j] + (gdown[i][site]*gdown[site][j])*F_down;
	  }

      for(i=1;i<=V;i++)
	{
	  a[i][site] = gup[i][site] + (gup[i][site]*gup[site][site] - gup[i][site])*F_up;
	  b[i][site] = gdown[i][site] + (gdown[i][site]*gdown[site][site] - gdown[i][site])*F_down;

	}

      for(i=1;i<=V;i++) 	
	for(j=1;j<=V;j++)
	  {
	    gup[i][j] = a[i][j];
	    gdown[i][j] = b[i][j];
	  }

      
      free_dmatrix(a,1,V,1,V);
      free_dmatrix(b,1,V,1,V);

    }      


}






