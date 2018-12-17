/***********************************************************/
/* hm3d_stable.c : Simulation of the Hubbard model in d=3  */
/*    Stabilization of the Green function computation      */
/*                                                         */
/*    I.Campos (Feb.01)  version 1.0                       */
/*                                                         */
/***********************************************************/

#include "hm3d.h"

extern double **K;                          /* Hopping matrix */

extern signed char spin[V+1][Ns_p1];             /* Ising model */

extern struct s_data data;

extern float t,u,Beta,delta_tau,lambda,mu;

extern int SIGN;

extern float exp_B_plus, exp_B_minus;


void M_GS(int n, double **u, double *s, double **v)
{

  /* Factorizes by the Modified Gram-Schmidt method the */
  /* nth order matrix u = u*w*v, where on output, u is  */
  /* orthonormal, w is diagonal with elements s and v is */
  /* upper triangular, unit diagonal. The original u is overwriten */

  int i,j,k;

  

  for(k=1;k<=n;k++)                    /* for each column */
    {

      s[k] = 0.;

      for(i=1;i<=n;i++)                /* first normalizes it */
	s[k] += (u[i][k]*u[i][k]);
	  
      s[k] = sqrt(s[k]);
      
      for(i=1;i<=n;i++)
	u[i][k] /= s[k]; 

      /* then makes the remaining columns orthogonal to it */

      v[k][k] = 1.0;

      for(j=k+1;j<=n;j++)
	{
	  v[k][j] = 0.;
	  for(i=1;i<=n;i++)
	    v[k][j] += u[i][k]*u[i][j];

	  for(i=1;i<=n;i++)
	    u[i][j] -= v[k][j]*u[i][k];

	  v[k][j] /=s[k];
	}

    }



}



void Inverter_tri(int n, double **inv_tri)   
{                                  
                                    /* on output the original is overwriten */

  int i,j,k,kp1;

  double t;


  for(k=1;k<=n;k++)
    {
      for(i=1;i<k;i++) inv_tri[i][k] = -inv_tri[i][k];

      kp1 = k + 1;
      for(j=kp1;j<=n;j++)
	{
	  t=inv_tri[k][j];
	  inv_tri[k][j]=0.0;

	  for(i=1;i<=k;i++) inv_tri[i][j] += t*inv_tri[i][k];
	}
    }


}

void stable_Green_function(int slice, double **gup, double **gdown)
{


  double **ident, **check_up, **check_down;
  
  int i,j,k,s,suma;
  int permut[Ns_p1],counter;
  double **tmpu, **tmpd;
  double d,Mup,Mdown;  
  double **aa,**bb;
  double **r1,**r2;
  double **U_up, **U_down;
  double **V_up, **V_down;
  double **saveV_up, **saveV_down;
  double *wup, *wdown;
  double **vup, **vdown;

  r1 = dmatrix(1,V,1,V);                    /*  Pointer initialization  */
  r2 = dmatrix(1,V,1,V);
  tmpu = dmatrix(1,V,1,V);
  tmpd = dmatrix(1,V,1,V);
  aa = dmatrix(1,V,1,V);
  bb = dmatrix(1,V,1,V);
  U_up = dmatrix(1,V,1,V);
  U_down = dmatrix(1,V,1,V);
  V_up = dmatrix(1,V,1,V);
  V_down = dmatrix(1,V,1,V);
  saveV_up = dmatrix(1,V,1,V);
  saveV_down = dmatrix(1,V,1,V);
  vup = dmatrix(1,V,1,V);
  vdown = dmatrix(1,V,1,V);
  wup = dvector(1,V);
  wdown = dvector(1,V);
  ident = dmatrix(1,V,1,V);

#ifdef CHECK

    check_up = dmatrix(1,V,1,V);
    check_down = dmatrix(1,V,1,V);

#endif
                                                /* Permutation vector */

  for(i=0;i<slice;i++) permut[i+1] = slice - i;

  for(i=0;i < (Ns - slice);i++) permut[slice+i+1] = Ns - i;
    


  for(i=1;i<=V;i++)                             /*  Auxiliar pointers  */
    for(j=1;j<=V;j++)
    {
      r1[i][j] = 0.0; 
      r2[i][j] = 0.0;
      U_up[i][j] = 0.0; 
      U_down[i][j] = 0.;
      V_up[i][j] = 0.0;
      V_down[i][j] = 0.;
    }

  for(i=1;i<=V;i++)    
    {
      r1[i][i] = 1.0;
      r2[i][i] = 1.0;
      U_up[i][i] = 1.0;
      U_down[i][i] = 1.0;
      V_up[i][i] = 1.0;
      V_down[i][i] = 1.0;
      wup[i] = 1.0;
      wdown[i] = 1.0;
    }

  counter = 0;

  for(i=Ns;i>0;i--)    
    {

      counter++;

      /* printf("counter = %d\n",counter); */


      Compute_B(permut[i],tmpu,tmpd);   /* Only one B is stored at a time */
      

      matrix_mult(V,tmpu,r1,aa);
      
      
      matrix_mult(V,tmpd,r2,bb);
      
      
      for(k=1;k<=V;k++)
	for(j=1;j<=V;j++)
	  {
	    r1[k][j] = aa[k][j];
	    r2[k][j] = bb[k][j];
	  }

      if(counter==MAX_MULT)
	{

	  counter = 0;

	  matrix_mult(V,aa,U_up,r1);
	  matrix_mult(V,bb,U_down,r2);



	  for(j=1;j<=V;j++)
	    for(k=1;k<=V;k++)
	      {
		r1[j][k] = r1[j][k]*wup[k];
		r2[j][k] = r2[j][k]*wdown[k];
	      }



	  M_GS(V,r1,wup,vup);
	  M_GS(V,r2,wdown,vdown);


	  for(j=1;j<=V;j++)
	    for(k=1;k<=V;k++)
	      {
		U_up[j][k] = r1[j][k];
		U_down[j][k] = r2[j][k];
	      }

	  matrix_mult(V,vup,V_up,r1);
	  matrix_mult(V,vdown,V_down,r2);

	  for(j=1;j<=V;j++)
	    for(k=1;k<=V;k++)
	      {
		V_up[j][k] = r1[j][k];
		V_down[j][k] = r2[j][k];
	      }

	  for(k=1;k<=V;k++)        /*  Auxiliar pointers reinitialized  */
	    for(j=1;j<=V;j++)
	      {
		r1[k][j] = 0.0;
		r2[k][j] = 0.0;
	      }

	  for(k=1;k<=V;k++) 
	    {   
	      r1[k][k] = 1.0; 
	      r2[k][k] = 1.0;
	    }

	
	} /******** End stabilization process ***********/



    }  /*********** End multiplication loop **************/


  /**** After the loop the matrix are set such that  *****/
  /**** A_up = U_up * D_up * V_up                      ***/
  /**** A_down = U_down * D_down * V_down              ***/

  /********** Adition of the identity to A(l)  ***********/


  for(i=1;i<=V;i++)              /* U are orthonormal  */
    for(j=1;j<=V;j++)
      {
	r1[i][j] = U_up[j][i];
	aa[i][j] = U_down[j][i];
	saveV_up[i][j] =V_up[i][j];
	saveV_down[i][j] =V_down[i][j];
      }



  Inverter_tri(V,V_up);               /* V are triangular   */

  Inverter_tri(V,V_down);


  matrix_mult(V,r1,V_up,tmpu);
  matrix_mult(V,aa,V_down,tmpd);
  

  for(i=1;i<=V;i++)
    {
      tmpu[i][i] += wup[i];
      tmpd[i][i] += wdown[i];
    }


  M_GS(V,tmpu,wup,vup);
  M_GS(V,tmpd,wdown,vdown);


  matrix_mult(V,U_up,tmpu,r1);      /* These are the finals U's   */
  matrix_mult(V,U_down,tmpd,r2);


  for(i=1;i<=V;i++)
    for(j=1;j<=V;j++)
      {
	U_up[i][j] = r1[i][j];
	U_down[i][j] = r2[i][j];
      }

  matrix_mult(V,vup,saveV_up,r1);      /* These are the final V's  */
  matrix_mult(V,vdown,saveV_down,r2);



#ifdef DETERMINANT
    
  Mup = Mdown = 1.0;
  

  for(i=1;i<=V;i++) Mup *=wup[i];
  for(i=1;i<=V;i++) Mdown *=wdown[i];

  printf("Mup = %g\t",Mup);
  printf("Mdown = %g\n",Mdown);

#endif


#ifdef CHECK                   /* Check the inversion of the matrix */

  /*  printf("Entro al check \n"); */


  for(i=1;i<=V;i++)
    for(j=1;j<=V;j++)
      {
	tmpu[i][j] = U_up[i][j]*wup[j];
	tmpd[i][j] = U_down[i][j]*wdown[j];
      }

  matrix_mult(V,tmpu,r1,check_up);     /* This is [1 + A(l)] up & down  */
  matrix_mult(V,tmpd,r2,check_down);

  /*
  for(i=1;i<=V;i++)
    for(j=1;j<=V;j++)
      printf("check_up[%d][%d] = %g\n",i,j,check_up[i][j]);
  */
#endif


  /******************* [1 + A(l)]_up = U_up * wup * r1 ******************/
  /******************* [1 + A(l)]_dw = U_dw * wdw * r2 ******************/
  /*                                                                    */
  /*           Inversion to obtain the Green function                   */
  /*                                                                    */
  /**********************************************************************/


                            /* r1 & r2 are triangular matrix */

  Inverter_tri(V,r1);
  Inverter_tri(V,r2);

  for(i=1;i<=V;i++)         /* transposed of U_up & U_down */
    for(j=1;j<=V;j++)
      {
	tmpu[j][i] = U_up[i][j];
	tmpd[j][i] = U_down[i][j];
      }

 
  for(i=1;i<=V;i++)         /* Multiplication V^(-1)*(1/w)  */
    for(j=1;j<=V;j++)
      {
	r1[i][j] = (1/wup[j])*r1[i][j];
	r2[i][j] = (1/wdown[j])*r2[i][j];
      }

  matrix_mult(V,r1,tmpu,gup);      /* Mult.  R^(-1) * (1/w) * U^t  */
  matrix_mult(V,r2,tmpd,gdown);


  /*
    for(i=1;i<=V;i++)
    for(j=1;j<=V;j++)
    printf("\t Up Inside the function Gup(%d,%d) = %g\n",i-1,j-1,gup[j][i]);
    
    for(i=1;i<=V;i++)
    for(j=1;j<=V;j++)
    printf("\t Down Inside the function Gdw(%d,%d) = %g\n",i,j,gdown[i][j]);
  */
  

  free_dmatrix(r1,1,V,1,V);                    /*  Pointer reset  */
  free_dmatrix(r2,1,V,1,V);
  free_dmatrix(tmpu,1,V,1,V);
  free_dmatrix(tmpd,1,V,1,V);
  free_dmatrix(aa,1,V,1,V);
  free_dmatrix(bb,1,V,1,V);
  free_dmatrix(U_up,1,V,1,V);
  free_dmatrix(U_down,1,V,1,V);
  free_dmatrix(V_up,1,V,1,V);
  free_dmatrix(V_down,1,V,1,V);
  free_dmatrix(saveV_up,1,V,1,V);
  free_dmatrix(saveV_down,1,V,1,V);
  free_dmatrix(vup,1,V,1,V);
  free_dmatrix(vdown,1,V,1,V);
  free_dvector(wup,1,V);
  free_dvector(wdown,1,V);
  free_dmatrix(ident,1,V,1,V);

#ifdef CHECK

    free_dmatrix(check_up,1,V,1,V);
    free_dmatrix(check_down,1,V,1,V);

#endif              


}


