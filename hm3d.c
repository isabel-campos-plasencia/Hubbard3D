

/**************************************************************************/
/*                                                                        */
/* hm3d.c : Simulation of the Hubbard model in d=3                        */
/*                                                                        */
/* I. Campos : December. 2000                                             */
/*************************************************************************/


                         /* MAIN PROGRAM  */

#define MAIN
#include "hm3d.h"
                             /* Declaracions of functions externs to MAIN */
extern void tiempo(void);
extern void Init_Rand(int); 
extern void Boundary_periodic(void);
extern void Read_Input(void);
extern void Initial(int);
extern void Generate_Ising(int);
extern void Generate_Diagonal_K();
extern void Write_Conf();
extern void Read_Conf();
extern void Read_Input();
extern void Write_Results();

/*
extern void Compute_B(int,double,double);
extern void Green_function(int,double,double);
extern void matrix_mult(int,double,double,double);
*/

extern void Measurements(void);
extern void Energy();

/* extern void HeatBath(int, int, double, double); */

                                           /* Global variables   */

int neigh_px,neigh_py,neigh_pz,neigh_pt,
    neigh_mx,neigh_my,neigh_mz,neigh_mt;

int x_p[L+1],y_p[L+1],z_p[L+1],t_p[L+1],
    x_m[L+1],y_m[L+1],z_m[L+1],t_m[L+1];

int times,TIMES;

double **Bup, **Bdown, **Gup, **Gdown;

double **K;

float exp_B_plus, exp_B_minus;

signed char spin[V+1][Ns_p1];           

int SIGN;

float Beta,t,u,mu,delta_tau,lambda,tz;

float prob0,prob1;                      /* This is for Metropolis */

float Ek,Eu;                            /*  Kinetic and Coulomb energy  */

float mag_pla[3],mag_stag,mag_ferro,S2;
int sigmaij,sigmaiijj;
float correlation[L+1];
float Mag_Ferro,Mag_Stag,Sigma_1rst,Sigma_2nd;

int it;

int neigh_px,neigh_mx,neigh_py,neigh_my,neigh_pz,neigh_mz;

int seed,flag,good;     

struct s_data data;

float obs_dat[n_obs][maxit];

float results[NOBS_HISTER];

int main()
{
                              /* Variables to be used only in 'MAIN' */

   int x,y,z,n,site,mfr,mesfr,ibin,i,j,l;     /* counters */
   int dummy,suma;
   double arg;
   float m_pl;
   float stag_s2, ferro_s2, average,diff, diff_min, diff_max;
   double **Knew;

   Bup = dmatrix(1,V,1,V);      /* Global Pointers Initialization  */
   Bdown = dmatrix(1,V,1,V);
   Gup = dmatrix(1,V,1,V);
   Gdown = dmatrix(1,V,1,V);
   K = dmatrix(1,V,1,V);
   Knew = dmatrix(1,V,1,V);

   Boundary_periodic();       /* Impose p.b.c.   */
 
   Read_Input();              /* Read the 'input' file  */
                   
   Init_Rand(data.seed);     /* Initialize random number generator */

   Initial(data.flag);      /* Sets starting conditions            */

   t=data.t;
   u=data.u;
   Beta = data.Beta;
   mesfr = data.mesfr;

#ifdef ANISOTROPY
   tz=data.tz;
#endif


#ifdef HALF_FILLING            /*  Chemical Potential  */
   mu = data.u/2;
#else
   mu = data.mu;
#endif
   
#ifndef HISTER

   Init_Prob();

   Jacobi_Generate_Diagonal_K(K);    /* Calculation of e^(-K) */
                                     /*  using Jacobi diagonalization */


   /* Generate_K(K); */   /* Calculation of e^(-K) using Trotter approx. */

   
   /***********************  Checking both approaches  ***********************/
   /*
   average = 0.;
   diff_min = 40.;
   diff_max = 0.;

   for(i=1;i<=V;i++)
     for(j=1;j<=V;j++)
       {
	 diff = (float) fabs(K[i][j] - Knew[i][j]);
	 average += diff;
	 diff_min=(diff<diff_min)?diff:diff_min;
	 diff_max=(diff>diff_max)?diff:diff_max;  
       }


   printf("Average difference = %g\n",average/(V*V));
   printf("Min diff = %g\t Max diff = %g\n",diff_min,diff_max);
   */
   /**************************************************************************/
   

#endif

   

   for(ibin=data.itcut;ibin<data.nbin;ibin++)    /* block number */
   {

     Init_Rand(data.seed);     /* Initialize random number generator */

    
#ifdef HISTER

   if(ibin<data.nbin/2)
     data.Beta +=data.dBeta;
   else
     data.Beta -=data.dBeta;
   
   Init_Prob();


   Generate_K(K);  /*needs recompute at each histeresis step */ 


#endif
                           
   good = 0;

   for(it=0;it<data.itmax;it++)       /* Loop in # of iterations   */
     {

       for(mfr=0;mfr<mesfr;mfr++)     /* Loop without measurements */
         { 

	   if(mfr == 0)
	     {
	       Ek = Eu = 0.;
	       S2 = 0.;
	       Sigma_1rst = Sigma_2nd = 0.;
	       Mag_Ferro = Mag_Stag = 0.;
	       ferro_s2 = stag_s2 = 0.;
	       for(l=0;l<3;l++) mag_pla[l] = 0.;
	       for(l=0;l< (L/2+1);l++) correlation[l] = 0.;
	     }

	   

	   for(l=1;l<Ns_p1;l++)       /* Loop in # of time slices  */
	     {

	       if(l == 1) 
		 {
		   times = 0;
		   TIMES = 4;
		   stable_Green_function(1,Gup,Gdown); 

		   if(mfr == mesfr - 1) Measure(1); 
		 }

	       site = 1;

	       for(x=1;x<=L;x++)
		 for(y=1;y<=L;y++)
		   for(z=1;z<=L;z++)
		     {
		       HeatBath(l,site,Gup,Gdown);   
		       site++;
		     }
	     

	       if(l < Ns)
		 {
		  
		   if(times == TIMES)
		     {
		       /* Check_Calc_Green(l);  */
		       times = 0;
		       stable_Green_function(l+1,Gup,Gdown);
		       
		     }
		   else
		     {
		       Green_lplusone(l,Gup,Gdown);
		       times++;		  
		     }
     

		   if(mfr == mesfr - 1) Measure(l+1);   

		 }
	       
	     }                               /* coord Nslices  */

	 }                                      /* end of mesfr */
   

       /*  MEASUREMENTS + WRITING  */


#ifdef DETERMINANT      
       suma =0;
       for(i=1;i<=V;i++)
	 for(j=1;j<Ns_p1;j++)
	   suma += spin[i][j]; 
       
       printf("Sum of ising spins = %d\t",suma);
       printf("factor = %g\n",exp( - lambda*delta_tau*suma));

       for(l=1;l<=Ns;l++)
	 stable_Green_function(l,Gup,Gdown);


#endif 

       obs_dat[0][it] = Ek*Normener; 
       obs_dat[1][it] = Eu*u*Normener; 
       obs_dat[2][it] = Mag_Ferro*Normener;
       obs_dat[3][it] = Mag_Stag*Normener;
       obs_dat[4][it] = Sigma_1rst*Normener/6;
       obs_dat[5][it] = Sigma_2nd*Normener/12;
       obs_dat[6][it] = (float) 0.75*S2* Normener;  


       m_pl=0.;
       for(j=0;j<3;j++)
	 m_pl+=(float) mag_pla[j] * (float) mag_pla[j];

       obs_dat[7][it] = (float)(sqrt(m_pl)*Normener);
       
       for(j=0;j<=L/2;j++)
	 obs_dat[8+j][it] = (float)correlation[j];

#ifndef HISTER  
       printf("\n MC iterations = %d\n",(it +1)*(mfr+1)*(ibin+1));
       printf(" Kinetic Energy =%g\n",Ek*Normener);
       printf(" Coulomb Energy =%g\n",Eu*u*Normener);
       printf("   Ferromagnetic Magnetization = %g\n",obs_dat[2][it]);
       printf("   Staggered Magnetization =%g\n",obs_dat[3][it]);
       printf("   Plane Magnetization =%g\n",obs_dat[7][it]); 
       printf("     Spin correlation 1rst neigh =%g\n",obs_dat[4][it]); 
       printf("     Spin correlation 2nd neigh =%g\n",obs_dat[5][it]); 
       printf("     Local magnetization squared =%g\n",obs_dat[6][it]); 
#endif 
       

     }                       /* End of  'itmax' */
   

#ifdef HISTER

   results[0] = obs_dat[0][data.itmax - 1];
   results[1] = obs_dat[1][data.itmax - 1];
   results[2] = obs_dat[2][data.itmax - 1];
   results[3] = obs_dat[3][data.itmax - 1];
   results[4] = obs_dat[4][data.itmax - 1];
   results[5] = obs_dat[5][data.itmax - 1];
   results[6] = obs_dat[6][data.itmax - 1];
   results[7] = obs_dat[7][data.itmax - 1];

   Write_Hister(ibin);

   printf("Kinetic Energy =%g\n",Ek*Normener);
   printf("Coulomb Energy =%g\n",Eu*u*Normener);
   printf("Ferromagnetic Magnetization = %g\n",Mag_Ferro);
   printf("Staggered Magnetization =%g\n",Mag_Stag);
   printf("Plane Magnetization =%g\n",results[7]); 
   printf("Spin-spin correlation 1rst neigh =%g\n",results[4]); 
   printf("Spin-spin correlation 2nd neigh =%g\n",results[5]); 
   printf("Local magnetization squared =%g\n",obs_dat[6][it]); 
#endif


   
#ifndef HISTER
   Write_Results(ibin);
#endif
  

   
   data.seed=RANDOM;                  /* Generate next random number */
   data.itcut=ibin+1;
   
   
   
#ifndef HISTER
   Write_Conf(0);
#endif
   
                    
}                                     /*  final de NBIN  */   
   

   free_dmatrix(Bup,1,V,1,V);         /* Global Pointers Reset  */
   free_dmatrix(Bdown,1,V,1,V);
   free_dmatrix(Gup,1,V,1,V);
   free_dmatrix(Gdown,1,V,1,V);
   free_dmatrix(K,1,V,1,V);
   free_dmatrix(Knew,1,V,1,V);
   
   return 0;
   
   
}  /*      END   OF   MAIN   PROGRAM    */

