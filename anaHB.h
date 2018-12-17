# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>


#define ANISOTROPY 

# define L              4
# define vol            (L*L*L)      /* 3D analysis */ 
# define maxit 		5000
# define maxbloque 	 20
# define nom_fich	 "OUT"
# define n_obs_medid     8 /* measured directly in the output */    
# define n_obs_FS	 9 
# define maxfiles	 5000
# define nbetas  	 10    /*  NO hace FS si nbetas = 0 */
# define n_inter	 100
# define n_correl	 200
# define maxmed		 200000
# define n_obs_plot      9


float v_dat[n_obs_FS][maxit],raw_dat[n_obs_medid][maxit];
long int frec[maxbloque][n_inter+1];
double xfrec[n_obs_FS][maxbloque][n_inter+1],SumO[n_obs_FS];

double ymed,x0,h,c1,c2,delta,derlO_v[n_obs_FS][nbetas+1],
        derlO_err[n_obs_FS][nbetas+1];

double coup_maxder[n_obs_FS][maxbloque+1], err_coup[n_obs_FS];
double maxder[n_obs_FS][maxbloque+1], err_maxder[n_obs_FS];
double O_v[n_obs_FS][nbetas+1], O_err[n_obs_FS][nbetas+1];
double derO_v[n_obs_FS][nbetas+1], derO_err[n_obs_FS][nbetas+1];


char   *ficheros[maxfiles],n_coup[10],n_deri[10], nombre[150],nom[100];


int    n_vol,permutar,d,medidas,iteraciones,m;

float Beta,u,t,mu,tz;

int SIGN;

FILE *Foutplt;

int nblo,lblo,n1,n2;


struct s_data
{
  int itmax,         /*  Number of iterations                            */
    mesfr,           /*  frequence of measurements                      */
    nbin,            /*  block number                                  */
    itcut,           /*  next block to compute                        */
    flag,            /*  starting conf: 0(random), 1(cold),2(backup) */
    seed;            /*  random seed                                */
  float t,           /* Hopping parameter                         */
    u,               /* Potential energy                         */
    Beta,             /* Inverse of the Temperature              */
    mu,               /*  Chemical potential                      */
    tz                /* Anisotropic hopping term */
  ;   

 };

struct s_data data;

char cadena[n_obs_FS][100] = {
  "Kinetic Energy",
  "Coulomb Energy",
  "Magnetization Ferro",
  "Magnetization Staggered",
  "Spin-spin 1rst neigh correlation",
  "Spin-spin 2nd neigh correlation",
  "Local magnetization squared",
  "Plane Magnetization",
  "Total Energy"
};









