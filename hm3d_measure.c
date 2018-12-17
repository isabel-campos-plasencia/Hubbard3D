/**********************************************************/
/* hm3d_measure.c : Simulation of the Hubbard model in d=3  */
/*         functions to measure observables               */
/*                                                        */
/*    I.Campos (December.00)  version 0.0                  */
/*                                                        */
/**********************************************************/
 
#include "hm3d.h"
 
extern double **K;                          /* Hopping matrix */
 
extern double **Bup,**Bdown,**Gup,**Gdown;

extern int x_p[L+1],y_p[L+1],x_m[L+1],y_m[L+1],z_p[L+1],z_m[L+1]; 

extern signed char spin[V][Ns_p1];                        
 
extern int neigh_px,neigh_mx,neigh_py,neigh_my,neigh_mz,neigh_pz;

extern float t,u,Beta,delta_tau,lambda;
 
extern float Ek,Eu,Mstag,S2;
extern int sigmaij,sigmaiijj;
extern float mag_pla[3],mag_stag,mag_ferro;
extern float correlation[];
extern float Mag_Ferro,Mag_Stag,Sigma_1rst,Sigma_2nd;

void Energy()
{

  int x,y,site,sitepx,sitepy,sitemx,sitemy;
  int z,sitepz,sitemz;
 
  site = 1;
 
  for(z=1;z<=L;z++)
    {
      neigh_pz = z_p[z];
      neigh_mz = z_m[z];

      for(y=1;y<=L;y++)
	{
	  neigh_py = y_p[y];
	  neigh_my = y_m[y];
	  
	  for(x=1;x<=L;x++)
	    {
	      neigh_px = x_p[x];
	      neigh_mx = x_m[x];
 
	      Eu += ( (1.0 - Gup[site][site])*(1.0 - Gdown[site][site]) )  ;

	      S2 += (Gup[site][site] + Gdown[site][site] -
		2*Gup[site][site]*Gdown[site][site]);

	      sitepx = site + neigh_px;
	      sitepy = site + neigh_py;
	      sitemx = site + neigh_mx;
	      sitemy = site + neigh_my;
	      sitepz = site + neigh_pz;
	      sitemz = site + neigh_mz;
	      
 
	      Ek += (Gup[sitepx][site] + Gdown[sitepx][site] +
		     Gup[sitepy][site] + Gdown[sitepy][site] + 
		     Gup[sitemx][site] + Gdown[sitemx][site] +
		     Gup[sitemy][site] + Gdown[sitemy][site] +
		     Gup[sitepz][site] + Gdown[sitepz][site] +
		     Gup[sitemz][site] + Gdown[sitemz][site]);
 
 
	      site++;
	    }
	}
    }    
      
}
  

void Magnetization(int t)               /* t indicates time slice  */
{
  int x,y,z,j;
  int k,sig_pla[3],sig_stag,sigx,sigy,sigz;
  long int a0[L+1],a1[L+1],a2[L+1];
  long int c0[L+1],c1[L+1],c2[L+1];

  for(k=0;k<3;k++)
    {
      sig_pla[k] = 1;
      sig_stag = 1;
    }

  for(k=1;k<=L;k++)
    {
      a0[k] = 0;
      a1[k] = 0;
      a2[k] = 0;
    }

  j=1;

  for(z=1;z<=L;z++)
    {
      neigh_pz = z_p[z];
      neigh_mz = z_m[z];

      for(y=1;y<=L;y++)
	{
	  neigh_py = y_p[y];
	  neigh_my = y_m[y];

	  for(x=1;x<=L;x++)
	    {
	      neigh_px = x_p[x];
	      neigh_mx = x_m[x];

	      sig_stag = (x+y+z)&1?1:-1;   /* Equivalent to (-1)^{x+y+z} */

	      sigx = (x)&1?1:-1;
	      sigy = (y)&1?1:-1;
	      sigz = (z)&1?1:-1;

	      sig_pla[0] = sigz;           /* XY plane  */
	      sig_pla[1] = sigy;           /* XZ plane  */
	      sig_pla[2] = sigx;           /* YZ plane  */

	      mag_ferro += spin[j][t];
	      mag_stag += spin[j][t]*sig_stag;

	      for(k=0;k<3;k++)
		mag_pla[k] +=sig_pla[k]*spin[j][t];

	      a0[x] += spin[j][t];
	      a1[y] += spin[j][t];
	      a2[z] += spin[j][t];
	      
	      j++;

	    }          /*  loop in x  */
	}          /*  loop in y  */
    }          /*  loop in z  */


  /* Correlation function between planes at distance K  */

  for(k=1;k<=(L/2 + 1);k++)
    {
      for(j=1;j<=L;j++)
	{
	  c0[k] += a0[j]*a0[(j+k)%L];
	  c1[k] += a1[j]*a1[(j+k)%L];
	  c2[k] += a2[j]*a2[(j+k)%L];
	}

      correlation[k] = (float) ( (c0[k] + c1[k] + c2[k])/3 );

    }     




}


void spin_spin(int n)
{

  int x,y,k,site,sitepx,sitepy,sitemx,sitemy;
  int z,sitepz,sitemz;
  int sitepp[3],sitepm[3],sitemp[3],sitemm[3];

  site = 1;
 
  for(z=1;z<=L;z++)
    {
      neigh_pz = z_p[z];
      neigh_mz = z_m[z];

      for(y=1;y<=L;y++)
	{
	  neigh_py = y_p[y];
	  neigh_my = y_m[y];
	  
	  for(x=1;x<=L;x++)
	    {
	      neigh_px = x_p[x];
	      neigh_mx = x_m[x];
 
	      sitepx = site + neigh_px;        /* 6 first neighboors */
	      sitepy = site + neigh_py;
	      sitemx = site + neigh_mx;
	      sitemy = site + neigh_my;
	      sitepz = site + neigh_pz;
	      sitemz = site + neigh_mz;


	      sitepp[0] = sitepx + neigh_py;    /* 12 second neighboors */
	      sitepp[1] = sitepx + neigh_pz;
	      sitepp[2] = sitepy + neigh_pz;

	      sitepm[0] = sitepx + neigh_my;
	      sitepm[1] = sitepx + neigh_mz;
	      sitepm[2] = sitepy + neigh_mz;

	      sitemp[0] = sitemx + neigh_py;
	      sitemp[1] = sitemx + neigh_pz;
	      sitemp[2] = sitemy + neigh_pz;

	      sitemm[0] = sitemx + neigh_my;
	      sitemm[1] = sitemx + neigh_mz;
	      sitemm[2] = sitemy + neigh_mz;

	      sigmaij +=  spin[site][n]*spin[sitepx][n];
	      sigmaij +=  spin[site][n]*spin[sitepy][n];
	      sigmaij +=  spin[site][n]*spin[sitepz][n];
	      sigmaij +=  spin[site][n]*spin[sitemx][n];
	      sigmaij +=  spin[site][n]*spin[sitemy][n];
	      sigmaij +=  spin[site][n]*spin[sitemz][n];

	      for(k=0;k<3;k++)
		{
		  sigmaiijj +=  spin[site][n]*spin[sitepp[k]][n];
		  sigmaiijj +=  spin[site][n]*spin[sitepm[k]][n];
		  sigmaiijj +=  spin[site][n]*spin[sitemp[k]][n];
		  sigmaiijj +=  spin[site][n]*spin[sitemm[k]][n];

		}


	      site++;

	    }
	}
    }


}


void Measure(int slice)
{

  mag_ferro = 0.;
  mag_stag = 0.;
  sigmaij = sigmaiijj = 0;
  
  Energy();
  
  Magnetization(slice);
  spin_spin(slice);
  
  Mag_Ferro +=  (float) fabs(mag_ferro);
  Mag_Stag +=  (float) fabs(mag_stag);
  Sigma_1rst += (float) fabs(sigmaij);
  Sigma_2nd += (float) fabs(sigmaiijj);


}
