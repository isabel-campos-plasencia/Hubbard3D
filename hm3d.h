

/***************************************************/
/* hm3d.h : Simulation of the Hubbard model in d=3 */
/*                Header file                      */
/*                                                 */
/* I.Campos (December 00)  version 1.0            */
/***************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include "nr.h"
#include "nrutil.h"

#define ANISOTROPY
#define MAX_MULT  6     /* Stabilization parameter  */


/* #define WRITE_CONF */ 
#define HALF_FILLING 
#define CHECK 
/* #define DETERMINANT  */


/* #define DEBUG  */
/* #define HISTER */


#ifndef CLOCKS_PER_SEC 
#define CLOCKS_PER_SEC 1000000L
#endif 

#define L           4                   /* LATTICE SIZE */
#define Ns          100                  /* TEMPORAL SIZE  */
#define Ns_p1       Ns + 1             /* dim. of arrays in temporal direct. */
#define  V  (L*L*L)

#define EPSILON 1.0e-3              /*  Precision in the inversion */
#define ERROR   1.0e-3               /* %  agreement in Green functions */


#define twopi  6.2831853            /* Some constants */
#define pi 3.1415927


#define maxit 5000      /* maximum number of iterations per bin */

#define n_obs  (9+L/2)       /* number of observables (put the right value) */

#define NOBS_HISTER 5

              /* random number generator (eventually using Martin's) */

#define NormRAN (1.0F/( (float) RAND_MAX+300.0F))
#define  RAN() ( (float) rand() * NormRAN )


#define FNORM   (2.3283063671E-10F)
#define RANDOM  ( (ira[ip++]=ira[ip1++]+ira[ip2++]) ^ira[ip3++] )
#define FRANDOM (FNORM*RANDOM)


#ifdef MAIN
unsigned char ip,ip1,ip2,ip3;
unsigned ira[256];
#else
extern unsigned char ip,ip1,ip2,ip3;
extern unsigned ira[];
#endif

#define randmax NormRAN


struct s_data
{
  int itmax,         /*  Number of iterations                              */
    mesfr,           /*  frequence of measurements                        */
    nbin,            /*  block number                                    */
    itcut,           /*  next block to compute                          */
    flag,            /*  starting conf: 0(random), 1(cold),2(backup)   */
    seed;            /*  random seed                                 */
  float t,           /* Hopping parameter (first neighboors)         */
    u,               /* Potential energy                            */
    Beta,             /* Inverse of the Temperature                 */    
    mu                /* chemical potential                        */
#ifdef ANISOTROPY
    ,tz                /* Anysotropy first neighboors        */
#endif
#ifdef HISTER
        ,dBeta        /* variation of beta in the histeresis        */
#endif
  ;   
};


#define NDATINT   6        /* number of int fields in s_data */
                           /* number of float fields in s_data */
#ifdef HISTER
#define NDATFLOAT 6        /* if nodef ANISOTRPY there is one field less */
#else
#define NDATFLOAT 5
#endif

#define Normener  ( (float) (1.0/((double) Ns*V)))  

#define LPATH 100 

FILE *Finput,*Foutput,*Fconfig;

