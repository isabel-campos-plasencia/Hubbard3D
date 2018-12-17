
/*******************************************************/
/* hm3d_ini.c : Simulation of the Hubbard model in d=3 */
/*              Initialization routines and I/O        */
/*                                                     */
/* I.Campos (December.00)  version 1.0                 */
/* The Computation of exp(K) using the Trotter aprox.  */
/* has been removed                                    */
/*******************************************************/


#include "hm3d.h"

extern int x_p[L+1],y_p[L+1],x_m[L+1],y_m[L+1],z_p[L+1],z_m[L+1];

extern signed char spin[V+1][Ns_p1];          /* Ising model */

extern int SIGN;

extern struct s_data data;

extern float obs_dat[n_obs][maxit],resultados[NOBS_HISTER];

extern float t,u,Beta,delta_tau,mu,lambda,tz;

extern float exp_B_plus,exp_B_minus,prob0,prob1;

extern float results[NOBS_HISTER];

char dir[LPATH];

void tiempo(void)
{
    static time_t time1=0,
                  time2;
    static clock_t clock1=0,
                   clock2;
    int temp;
    float t_CPU;
 
    time2 = time(NULL);
    clock2 = clock();


    if (time1>0)
    {
        temp=time2-time1;
        t_CPU=(0.5+(float)(clock2-clock1))/CLOCKS_PER_SEC;
        if (temp)
            printf("%4ds (%2.0f%%): ",temp,t_CPU/(temp+0.5)*100.);
        else
            printf("%4ds (??%%): ",temp);
    }        
    time1=time2;
    clock1=clock2;
}           

void Init_Prob()       /* Some constants related to the Update probability */
{
    
  delta_tau = data.Beta/Ns;  
  
  lambda = 2*atan(sqrt(tanh(delta_tau*u/4))) /delta_tau ;    
  
  /* This is to compute B */
  
  exp_B_plus = exp(-delta_tau*lambda - delta_tau*(mu - u/2));
  exp_B_minus = exp(delta_tau*lambda - delta_tau*(mu - u/2));
  
  
  
  prob0 = exp(2*delta_tau*lambda) - 1.0;       /*  Metropolis exponential  */
  prob1 = exp(-2*delta_tau*lambda) - 1.0;
  

}



void Init_Rand(int seed)    
{                              /*  One uses the standard C generator */
  int i;                       /* to initialize a safer one          */ 
                                
  srand((unsigned)seed);

  for (i=0;i<111;i++)         /* We through away the first ones */
    rand();                   /* for further security           */

    ip=128;    
    ip1=ip-24;    
    ip2=ip-55;    
    ip3=ip-61;
    
    for (i=0; i<256; i++)
        ira[i] = (unsigned) rand()+ (unsigned) rand(); 

    for (i=0;i<1111;i++)      
        RANDOM;
}


void Boundary_periodic(void)
{
    int i;

    for (i=1;i<=L;i++)
    {
	x_p[i]= 1;
	x_m[i]=-1;
	y_p[i]= L;
	y_m[i]=-L;
	z_p[i]= L*L;
	z_m[i]= -L*L;
    }
    x_m[1]= L-1;
    y_m[1]=(L-1)*L;
    z_m[1]=(L-1)*L*L;
    x_p[L]=-x_m[1];
    y_p[L]=-y_m[1];
    z_p[L]=-z_m[1];
} 


void Generate_Ising(int flag)
{

  int i,k;
  int x,y,z,sig_stag;
  double rn;

  if(flag==-1)                   /* Staggered field  */
    {
      for(k=1;k<Ns_p1;k++)
	{
	  i=1;

	  for(z=1;z<=L;z++)
	    for(y=1;y<=L;y++)
	      for(x=1;x<=L;x++)
		{
		  sig_stag = (x+y+z)&1?1:-1;   /* Equivalent to (-1)^{x+y+z} */
		  
		  spin[i][k] = sig_stag;
		  i++;
		}
	}
    }
		   

  if (flag >= 0)
    {
      if(flag != 2)
	{ 
	  for(i=1;i<=V;i++)
	    for(k=1;k<Ns_p1;k++)
	      {
		if(flag == 0) 
		  {    
		    rn = RAN();                         /* random start  */
		    spin[i][k] = 2*((int)(rn+0.5))-1;
		    //   printf("spin[%d][%d] =%d\t",i-1,k-1,spin[i][k]);
		    // printf("rn = %g\n",rn);
		  }
		else
		  {
		    spin[i][k] = 1;             /* flag = 1 , ferro start */
		  }
	      }
	}
      
    }

}


void Write_Conf()
{

  char name[LPATH+33],name_dollar[LPATH+33],name_old[LPATH+33];
  int i;

  sprintf(name_dollar,"%s%s",dir,"conf.$$$");
  sprintf(name,"%s%s",dir,"conf");
  sprintf(name_old,"%s%s",dir,"conf.old");


  Fconfig=fopen(name_dollar,"wb");

  fwrite(&data,sizeof(data),1,Fconfig);
  fwrite(&SIGN,sizeof(int),1,Fconfig);

#ifdef DEBUG
    Show_Input(&data);
#endif
    
    for(i=1;i<=V;i++)
      fwrite(&spin[i][1],sizeof(char),Ns,Fconfig);

    fclose(Fconfig);
    remove(name_old);
    rename(name,name_old);
    rename(name_dollar,name);


}


void Read_Conf(int i)
{

  struct s_data datab;
  char name[LPATH+33];
  int j;    
  
  
  sprintf(name,"%s%s",dir,"conf");
  
  Fconfig=fopen(name,"rb");
  if (Fconfig==NULL)
    {
      printf(" The file '%s' does not exist. Aborting program.\a\n",name);
      exit;
    }
  
    fread(&datab,sizeof(datab),1,Fconfig);
    fread(&SIGN,sizeof(int),1,Fconfig);
 
    for(j=1;j<=V;j++)
      fread(&spin[j][1],sizeof(char),Ns,Fconfig);

#ifdef DEBUG
    Show_Input(&datab);
#endif

    if (data.itmax   != datab.itmax ||
        data.mesfr != datab.mesfr ||
        data.Beta   != datab.Beta ||
        data.t   != datab.t ||
	data.mu != datab.mu ||
	data.u != datab.u ||
        data.flag  == 3 )
      {
        printf(" Warning: data in '%s' are not compatible with those of 'input'\n",name);
        
        printf(" \t\t `input` parameters will be used instead. \n");
	
        sprintf(name,"%s%s%03d.DAT",dir,"OUT",data.itcut);
	
        Foutput=fopen(name,"rb");
        if (Foutput!=NULL)
	  {
            fclose(Foutput);
            printf(" %s  already exists.\a\n",name);
            exit;
        }
      }
    else
      {
        data.itcut=datab.itcut;
        data.seed=datab.seed;
      }
    fclose(Fconfig);
}




void Read_Input(void)
{
    int j;
    int * ptdata_int;
    float * ptdata_real;
    char name[33];
       

    Finput=fopen("input","r");
    if (Finput==NULL)
    {
        printf(" The file 'input' does not exist. Aborting program. \n");
        exit(0);
    }
     

    fscanf(Finput,"%s",dir);  

    for (j=0,ptdata_int=&data.itmax;j<NDATINT;j++)
        fscanf(Finput,"%d",ptdata_int++);
    for (j=0,ptdata_real=&data.t;j<NDATFLOAT;j++)
        fscanf(Finput,"%f",ptdata_real++);
    fclose(Finput);


    Show_Input(&data);      /* Print information about run on the screen */
       

    if(data.flag < 2)       /* If this  is a rum continuation  */
    {
       sprintf(name,"%s%s%03d.DAT",dir,"OUT",data.itcut);

        Foutput=fopen(name,"rb");
        if (Foutput!=NULL)
        {
            fclose(Foutput);
            printf("Warning: %s  already exists.\a\n",name);
            exit;
        }

       sprintf(name,"%s%s",dir,"conf");

        Foutput=fopen(name,"rb");
        if (Foutput!=NULL)
        {
            fclose(Foutput);
            printf("Warning:  %s exist already.\a\n",name);
            exit;
        }
    }

}

void Write_Results(int i)
{

  int idat;
  char name[LPATH + 33], nameconf[LPATH + 33];

  sprintf(name,"%s%s%03d.DAT",dir,"OUT",i);
  sprintf(nameconf,"%s%s%03d.DAT",dir,"CONF",i);

  Foutput =fopen(name,"wb");

  fwrite(&data,sizeof(data),1,Foutput);
  fwrite(&SIGN,sizeof(int),1,Foutput);


  for(idat=0;idat<n_obs;idat++)
    fwrite(&obs_dat[idat][0],4*data.itmax ,1,Foutput);
  fclose(Foutput);
    

#ifdef WRITE_CONF
                          /* Here writes the configuration  */

  Foutput =fopen(nameconf,"wb");

  fwrite(&data,sizeof(data),1,Foutput);
  fwrite(&SIGN,sizeof(int),1,Foutput);

  for(i=1;i<=V;i++)
    fwrite(&spin[i][1],sizeof(char),Ns,Foutput);

  fclose(Foutput);

#endif


  /****************************************************************/

}



int Show_Input(struct s_data *dat)   
{
  printf("/*********************************************/ \n");

  printf("/**  Simulation of the Hubbard model in d=3 **/ \n\n");
  printf("\t\t PARAMETERS OF THE RUN \n\n");
  printf("\t Lattice size = %d\n",L);
  printf("\t Time slices = %d\n",Ns);

  printf("/*********************************************/ \n\n");
    printf("\t Iterations per bin:  %d \n",dat->itmax);
    printf("\t Freq. of measurements:  %d \n",dat->mesfr);
    printf("\t Number of blocks:   %d \n",dat->nbin);
    printf("\t Next block to compute : %d \n",dat->itcut);
    printf("\t Starting flag :   %d \n",dat->flag);
    printf("\t Random number seed:   %d \n",dat->seed);
    printf("\n");
    printf("\t Hopping parameter t= %g\n",dat->t);
    printf("\t Hopping parameter in Z direction tz = %g\n",dat->tz);
    printf("\t Coulomb coupling U = %g\n",dat->u);
    printf("\t Inverse Temperature Beta = %g\n",dat->Beta);
    printf("\t Corresponding to delta_tau = %g\n\n",dat->Beta/Ns);
#ifdef HALF_FILLING
    printf("\t Simulation done at half-filling: mu = %g\n\n",dat->u/2);
#else
    printf("\t Chemical potential : mu = %g\n",dat->mu);
#endif
    if(dat->t != dat->tz) printf("\t Simulation with anisot. hopping parameter\n");

  printf("/*********************************************/ \n\n");

#ifdef HISTER
    printf("Interval dBeta %g\n",dat->dBeta);
#endif
    
    return(0);
}         


void Initial(int flag)
{

  if( flag  >= 2)              /*  Restart from a backup Configuration  */
    Read_Conf(0); 
  
  if( flag < 2 )  
    Generate_Ising(flag);      /*  Start from scratch  */ 
  

}

void Write_Hister(int i)
{
  FILE *Fhister;
  char name[100];
  int k;

  sprintf(name,"%s%s_%g.%d",dir,"hister",data.Beta,L);

  Fhister=fopen(name,"a");
  fprintf(Fhister,"%g\t",data.Beta);
  for(k=0;k<NOBS_HISTER;k++)
    fprintf(Fhister,"%g\t",results[k]);
  fprintf(Fhister,"\n");
  fclose(Fhister);

  printf("(Beta = %5.3f, U =%5.3f)\n",data.Beta,data.u);

}



void Jacobi_Generate_Diagonal_K(double **kappa)
{

  int i,j,k,x,y,z,site;
  int sitepx,sitepy,sitemx,sitemy,sitepz,sitemz;
  int n_px,n_py,n_mx,n_my,n_pz,n_mz;

  double **ortog,**temp;
  double d[V+1];
  int nrot;

  ortog = dmatrix(1,V,1,V);
  temp = dmatrix(1,V,1,V);

  for(i=1;i<=V;i++)
    for(j=1;j<=V;j++)
      kappa[i][j] = 0.0;

  site = 1;

  for(z=1;z<=L;z++)
    {
      n_pz = z_p[z];
      n_mz = z_m[z];

      for(y=1;y<=L;y++)             /* Construction of the Hopping matrix */
	{
	  n_py = y_p[y];
	  n_my = y_m[y];
	  
	  for(x=1;x<=L;x++)
	    {
	      n_px = x_p[x];
	      n_mx = x_m[x];
	      
	      sitepx = site + n_px;       /* 6 first neighboors */
	      sitemx = site + n_mx;
	      sitepy = site + n_py;
	      sitemy = site + n_my;
	      sitepz = site + n_pz;
	      sitemz = site + n_mz;
	  

	      /*  Filling the Hopping matrix */
                                         
                                                /* First neighboors coupling */

	      kappa[site][sitepx] = -t;
	      kappa[site][sitemx] = -t;
	      kappa[site][sitepy] = -t;
	      kappa[site][sitemy] = -t;
	      kappa[site][sitepz] = -t;
	      kappa[site][sitemz] = -t;
	      

#ifdef ANISOTROPY

	      kappa[site][sitepz] = -tz;
	      kappa[site][sitemz] = -tz;

#endif
              /*********************************/

	      site++;
	    }
	}
    }

  jacobi(kappa,V,d,ortog,&nrot);     /* Diagonalization via Jacobi rotations */

  /*  for(i=1;i<=V;i++)
      printf("Eigenvalue d[%d] = %g\n",i,d[i]);  */
  
  for(i=1;i<=V;i++)  
    for(j=1;j<=V;j++)  
      temp[i][j] = ortog[j][i];
  
  for(i=1;i<=V;i++)                  /* Similarity transformation */
    for(j=1;j<=V;j++)
      {
	ortog[i][j] = exp(-delta_tau*d[j])*ortog[i][j]; 
      }
  
  matrix_mult(V,ortog,temp,kappa);
  
  free_dmatrix(ortog,1,V,1,V);
  free_dmatrix(temp,1,V,1,V);
  
  
/*
  for(j=1;j<=V;j++)
    for(i=1;i<=V;i++)
      printf("Hopping matrix Kappa(%d %d) = %g\n",i-1,j-1,kappa[i][j]);
*/
  


}


void Generate_K(double **kappa)
{

  int i,j,x,y,z,site;         /* K(i,j) = -t , if i,j are first neighboors */
                              /*           0    otherwise                  */


  int contador,index1,index2,k;
  int sitepx,sitepy,sitemx,sitemy,sitepz,sitemz;
  int n_px,n_py,n_mx,n_my,n_pz,n_mz;
  float argument;
  double **Tpx, **Tpy, **Tpz;
  double **T1, **Tp;


  Tpx = dmatrix(1,V,1,V);
  Tpy = dmatrix(1,V,1,V);
  Tpz = dmatrix(1,V,1,V);
  T1 = dmatrix(1,V,1,V);
  Tp = dmatrix(1,V,1,V);


  for(i=1; i<=V ; i++)         /* kappa is initially set to the identity  */
    {
      kappa[i][i] = 1.0;
      for(j=1; j<i ; j++)
	kappa[i][j] = kappa[j][i] = 0.0;
    }	      
  

  argument = delta_tau * t; 

  site = 1;

  for(z=1;z<=L;z++)
    {
      n_pz = z_p[z];

      for(y=1;y<=L;y++)
	{
	  n_py = y_p[y];
      
	  for(x=1;x<=L;x++)
	    {
	      n_px = x_p[x];  
	      
	      for(i=1; i<=V ; i++)         /*  Initially set to the identity  */
		{
		  Tpx[i][i] = 1.0;
		  Tpy[i][i] = 1.0;
		  Tpz[i][i] = 1.0;
		 	      
		  for(j=1; j<i ; j++)
		    {
		      Tpx[i][j] = Tpx[j][i] = 0.0;
		      Tpy[i][j] = Tpy[j][i] = 0.0;
		      Tpz[i][j] = Tpz[j][i] = 0.0;
		    }
		}

	      /*  Calculation of the first neighboors    */
	      
	      sitepx = site + n_px; 
	      sitepy = site + n_py;
	      sitepz = site + n_pz;
	      
	      /******************  direction +0 *****************************/
	      
	      Tpx[site][site] = cosh(argument);       
	      Tpx[site][sitepx] = sinh(argument);
	      Tpx[sitepx][site] = sinh(argument);
	      Tpx[sitepx][sitepx] = cosh(argument);
	      /******************  direction +1 *****************************/
	      
	      Tpy[site][site] = cosh(argument);      
 	      Tpy[site][sitepy] = sinh(argument);
	      Tpy[sitepy][site] = sinh(argument);
	      Tpy[sitepy][sitepy] = cosh(argument); 	      

	      /******************  direction +2 ****************************/
	      
	      Tpz[site][site] = cosh(argument);      
	      Tpz[site][sitepz] = sinh(argument);
	      Tpz[sitepz][site] = sinh(argument);
	      Tpz[sitepz][sitepz] = cosh(argument);
	      
	      /********************  matrix multiplication *****************/   
	      
	      
	      matrix_mult(V,Tpx,Tpy,T1);

	      matrix_mult(V,T1,Tpz,Tp);

	      matrix_mult(V,kappa,Tp,T1);    /* Single site contribution */

	      for(i=1;i<=V;i++)
		for(j=1;j<=V;j++)
		  kappa[i][j] = T1[i][j];

	      site++;
	  
	    }
	}
      
    }
  

  free_dmatrix(Tpx,1,V,1,V);
  free_dmatrix(Tpy,1,V,1,V);
  free_dmatrix(Tpz,1,V,1,V);
  free_dmatrix(T1,1,V,1,V);
  free_dmatrix(Tp,1,V,1,V);

  

  for(i=1;i<=V;i++)
    for(j=1;j<=V;j++)
      printf("Hopping matrix kappa(%d %d) = %f\n",i-1,j-1,kappa[i][j]);

  

}       


