
/**********************************************************/
/* hm3d_green.c : Simulation of the Hubbard model in d=3  */
/*         Calculation of Green functions                 */
/*                                                        */
/*    I.Campos (Dec.00)  version 1.0                      */
/*                                                        */
/**********************************************************/


#include "hm3d.h"


extern double **K;                          /* Hopping matrix */

extern signed char spin[V+1][Ns_p1];             /* Ising model */

extern int Check_gf;

extern struct s_data data;

extern float t,u,Beta,delta_tau,lambda,mu;

extern int SIGN,times,TIMES;

extern float exp_B_plus, exp_B_minus;

void Compute_B(int slice, double **bup, double **bdown)
{

  int i,j;
  double eNu[V+1];
  float tmp1,tmp2;
                                 /*  Construct the exponential of nu() up  */
 
  for(i=1;i<=V;i++)                 
    {
      if(spin[i][slice] == 1) eNu[i] = exp_B_plus;
      if(spin[i][slice] == -1) eNu[i] = exp_B_minus; 

    }

  for(i=1;i<=V;i++)
    for(j=1;j<=V;j++)
      bup[i][j] = eNu[i]*K[i][j];

                                   /* Construct the exponential of nu() down  */
  for(i=1;i<=V;i++)    
    {
      if(spin[i][slice] == 1) eNu[i] = exp_B_minus;
      if(spin[i][slice] == -1) eNu[i] = exp_B_plus;  
    }
    

  for(i=1;i<=V;i++)
    for(j=1;j<=V;j++)
      bdown[i][j] = eNu[i]*K[i][j];
                                       

  /* #ifdef DEBUG */

  /******************* Check the matrices ************************/
 
  /*  
  printf("Calculation of B at l=%d\n",slice-1);
  for(j=1;j<=V;j++)                   
    for(i=1;i<=V;i++)
      {
	printf("\t\t Bdown(%d,%d) = %g\n",i-1,j-1,bdown[i][j]);    
      }

  for(j=1;j<=V;j++)                   
    for(i=1;i<=V;i++)
      {
	printf("\t\t Bup(%d,%d) = %g\n",i-1,j-1,bup[i][j]);  
      }
  */
  
  /****************************************************************/  
  
    /* #endif */

} 






void Green_function(int slice, double **gup, double **gdown)
{


  double **CheckG, **ident;
  
  int i,j,k,s,*indx,suma;
  int permut[Ns_p1];
  double **y,*col;
  double **tmpBu, **tmpBd;
  double d,Mup,Mdown;  
  double **aa,**bb;
  double **r1,**r2;

  r1 = dmatrix(1,V,1,V);                    /*  Pointer initialization  */
  r2 = dmatrix(1,V,1,V);
  tmpBu = dmatrix(1,V,1,V);
  tmpBd = dmatrix(1,V,1,V);
  aa = dmatrix(1,V,1,V);
  bb = dmatrix(1,V,1,V);

                                                /* Permutation vector */

  for(i=0;i<slice;i++) permut[i+1] = slice - i;

  for(i=0;i < (Ns - slice);i++) permut[slice+i+1] = Ns - i;
    


  for(i=1;i<=V;i++)                             /*  Auxiliar pointers  */
    for(j=1;j<=V;j++)
    {
      r1[i][j] = 0.0;
      r2[i][j] = 0.0;
    }

  for(i=1;i<=V;i++)    
    {
      r1[i][i] = 1.0;
      r2[i][i] = 1.0;
    }



  for(i=Ns;i>0;i--)    
    {

      Compute_B(permut[i],tmpBu,tmpBd);    /* Only one B is stored at a time */
      

      matrix_mult(V,tmpBu,r1,aa);

      
      matrix_mult(V,tmpBd,r2,bb);
      
      
      for(k=1;k<=V;k++)
	for(j=1;j<=V;j++)
	  {
	    r1[k][j] = aa[k][j];
	    r2[k][j] = bb[k][j];
	  }

    }                                    /* End of multiplication loop */

  /*
  for(i=1;i<=V;i++)
    for(j=1;j<=V;j++)
      printf("aa(%d,%d) = %g\n",i,j,aa[i][j]);

  for(i=1;i<=V;i++)
    for(j=1;j<=V;j++)
      printf("bb(%d,%d) = %g\n",i,j,bb[i][j]);
  */
  
  /************  Adding unity to the diagonal of A(l) *********************/
  

  for(i=1;i<=V;i++)
    {
      aa[i][i] = 1.0 + aa[i][i];
      bb[i][i] = 1.0 + bb[i][i];
    }


  /************ Inversion to obtain the  Up Green function ******************/


  indx = ivector(1,V);
  ident =  dmatrix(1,V,1,V);

#ifdef CHECK

  CheckG = dmatrix(1,V,1,V);

  for(i=1;i<=V;i++)
    for(j=1;j<=V;j++)
      CheckG[i][j] = aa[i][j];

#endif

  ludcmp(aa,V,indx,&d);                   


  /***********  Calculation of the sign of the determinant *************/
  
  SIGN = d;

  for(i=1;i<=V;i++)
    {
      if(aa[i][i] < 0.0) SIGN=-SIGN;
      if(aa[i][i] == 0.0) SIGN = 0;         /* exceptional configuration */
    }
      
#ifdef DETERMINANT
  for(i=1;i<=V;i++)  d *= aa[i][i];
  Mup = d;
  printf(" M_up = %g\t",Mup); 
#endif
  /*********************************************************************/



  y=dmatrix(1,V,1,V);
  col = dvector(1,V);

  for(j=1;j<=V;j++)
    {
      for(i=1;i<=V;i++) col[i]=0.0;
      col[j] = 1.0;

      lubksb(aa,V,indx,col);

      for(i=1;i<=V;i++) y[i][j] = col[i];
    }

  for(i=1;i<=V;i++)
    for(j=1;j<=V;j++)
	gup[i][j] = y[i][j];


#ifdef DEBUG

  printf("Inside the function  l =%d\n",slice);

  for(i=1;i<=V;i++)
    for(j=1;j<=V;j++)
      printf("\t Up Inside the function G(%d,%d) = %g\n",i,j,gup[i][j]);

#endif 


#ifdef CHECK

  /************************  Check the inversion ***************************/
  
  matrix_mult(V,CheckG,y,ident);

  for(i=1;i<=V;i++)
    {
      if( fabs(ident[i][i]) - 1.0 > EPSILON) 
	printf("\t Precision error: ident(%d,%d) = %g\n",i,i,ident[i][i]);
      
      for(j=1;j<i;j++)
	{
	  if( (fabs(ident[i][j]) > EPSILON) || (fabs(ident[j][i]) > EPSILON))
	    printf("\t Precision error: ident(%d,%d) = %g\n",i,j,ident[i][j]);
	}
    }
  
		  
  /**************************************************************************/

#endif


  /*********** Inversion to obtains the  Down Green function  **************/



#ifdef CHECK

  for(i=1;i<=V;i++)
    for(j=1;j<=V;j++)
      CheckG[i][j] = bb[i][j];
 
#endif

 
  ludcmp(bb,V,indx,&d);                 


  /*****************  Calculation of the determinant *******************/
  
  SIGN = d;

  for(i=1;i<=V;i++)
    {
      if(bb[i][i] < 0.0) SIGN=-SIGN;
      if(bb[i][i] == 0) SIGN = 0;     /* exceptional configuration */
    }
      

#ifdef DETERMINANT
  for(i=1;i<=V;i++)  d *= bb[i][i];
  Mdown = d;
  printf("M_down = %g\t ratio =%g\n ",Mdown,Mup/Mdown); 
#endif
  /*********************************************************************/

  
  for(j=1;j<=V;j++)
    {
      for(i=1;i<=V;i++) col[i]=0.0;
      col[j] = 1.0;
      
      lubksb(bb,V,indx,col);
      
      for(i=1;i<=V;i++) y[i][j] = col[i];
    }

  
  for(i=1;i<=V;i++)                     /* Asigns result  */
    for(j=1;j<=V;j++)
	gdown[i][j] = y[i][j];




#ifdef DEBUG

  printf("\t\t Iteration =%d\n",data.itmax);

  for(i=1;i<=V;i++)
    for(j=1;j<=V;j++)
      printf("\t Down G(%d,%d) = %g\n",i,j,gdown[i][j]);

#endif


#ifdef CHECK

  /************************  Check ***********************************/
  
  matrix_mult(V,CheckG,y,ident);

  for(i=1;i<=V;i++)
    {
      if( fabs(ident[i][i]) - 1.0 > EPSILON)
	printf("\t Precision error: ident(%d,%d) = %g\n",i,i,ident[i][i]);
      
      for(j=1;j<i;j++)
	{
	  if( (fabs(ident[i][j]) > EPSILON) || (fabs(ident[j][i]) > EPSILON))
	    printf("\t Precision error: ident(%d,%d) = %g\n",i,j,ident[i][j]);
	}
    }
  
		  
  /**************************************************************************/

#endif


  free_dmatrix(tmpBu,1,V,1,V);
  free_dmatrix(tmpBd,1,V,1,V);
  free_dmatrix(aa,1,V,1,V);
  free_dmatrix(bb,1,V,1,V);
  free_dmatrix(r1,1,V,1,V);
  free_dmatrix(r2,1,V,1,V);
  free_dmatrix(ident,1,V,1,V);
#ifdef CHECK
  free_dmatrix(CheckG,1,V,1,V);
#endif
  free_dmatrix(y,1,V,1,V);
  free_dvector(col,1,V);
  free_ivector(indx,1,V);
 
}




void Green_lplusone(int n, double **gup, double **gdown)
{

  double **bup, **bdown, **aa;
  int i,j;
  int *indx;
  double **y,*col;
  double d;

  bup = dmatrix(1,V,1,V);
  bdown = dmatrix(1,V,1,V);
  aa = dmatrix(1,V,1,V);
  

  Compute_B(n+1,bup,bdown);      /* Computation of B(n+1)  */

  matrix_mult(V,bup,gup,aa);     /* B(n+1)*G(n)  */ 


                                 /* Inversion of B(l+1) up */  
  indx = ivector(1,V);  
  ludcmp(bup,V,indx,&d);  
  y=dmatrix(1,V,1,V);   
  col = dvector(1,V);  

  for(j=1;j<=V;j++)     
    {  
      for(i=1;i<=V;i++) col[i]=0.0;       
      col[j] = 1.0;  
      lubksb(bup,V,indx,col);  
      for(i=1;i<=V;i++) 
	y[i][j] = col[i];     
    }
  
  matrix_mult(V,aa,y,gup);      /*  B(n+1)*G(n)*B(n+1)^(-1)  */              



  matrix_mult(V,bdown,gdown,aa);

                                /* Inversion of B(l+1) down */   
  ludcmp(bdown,V,indx,&d);  
  for(j=1;j<=V;j++)     
    {  
      for(i=1;i<=V;i++) col[i]=0.0;       
      col[j] = 1.0;  
      lubksb(bdown,V,indx,col);  
      for(i=1;i<=V;i++) 
	y[i][j] = col[i];     
    }
  
  matrix_mult(V,aa,y,gdown);    /*   B(n+1)*G(n)*B(n+1)^(-1)  */                


  free_dmatrix(bup,1,V,1,V);
  free_dmatrix(bdown,1,V,1,V);
  free_dmatrix(aa,1,V,1,V);
  free_dmatrix(y,1,V,1,V);
  free_dvector(col,1,V);
  free_ivector(indx,1,V);

}

Check_Calc_Green(int n)
{

  double diff_up, diff_down;
  double **gup, **gdown, **aa, **bb;
  int i,j,precision;


  gup = dmatrix(1,V,1,V);
  gdown = dmatrix(1,V,1,V);
  aa = dmatrix(1,V,1,V);
  bb = dmatrix(1,V,1,V);


  stable_Green_function(n,gup,gdown);

  Green_lplusone(n, gup, gdown);

  stable_Green_function(n+1, aa, bb);

  precision = 1;

  for(i=1;i<=V;i++)
    for(j=1;j<=V;j++)
      {

	diff_up = fabs(gup[i][j] - aa[i][j]);
	diff_down = fabs(gdown[i][j] - bb[i][j]);

	if( diff_up > ERROR ) 
	  { 
	    printf("Difference diff_up[%d][%d] = %g\n",i,j,diff_up);
	    precision = 0;
	  }
	else
	  {
	    if( diff_down > ERROR ) 
	      {
		printf("Difference diff_down[%d][%d] = %g\n",i,j,diff_down);
		precision = 0;
	      }
	  }

	if(precision == 0)
	  {
	    TIMES = 1;
	    printf("\t\t salgo de la funcion \n");
	    goto end;
	  }


      }


  TIMES++;

 end:
  printf("precision = %d\t TIMES = %d\n",precision,TIMES);
  free_dmatrix(gup,1,V,1,V);
  free_dmatrix(gdown,1,V,1,V);
  free_dmatrix(aa,1,V,1,V);
  free_dmatrix(bb,1,V,1,V);

}

