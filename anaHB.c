/***************************************************/
/*                                                 */
/*  Analysis program for the Hubbard Model in D=3  */
/*  I. Campos                                      */
/***************************************************/


#include "anaHB.h"


void describe_data(struct s_data *dat)
{
  printf("/*********************************************/\n\n");
  printf("\t Analysis Program: Hubbard Model in D=3 \n\n");
  if(mu == u/2) printf(" Simulation performed at half-filling\n");
  printf("\n The input corresponds to: \n\n");
  printf("\t\t itmax  %d \n", dat->itmax);
  printf("\t\t mesfr  %d \n", dat->mesfr);
  printf("\t\t nbin  %d \n", dat->nbin);
  printf("\t\t itcut  %d \n", dat->itcut);
  printf("\t\t flag  %d \n", dat->flag);
  printf("\t\t seed  %d \n", dat->seed);
  printf("\t\t Hopping parameter  %g \n", dat->t);
  printf("\t\t Coulomb interaction  %g \n", dat->u);
  printf("\t\t Beta  %g \n", dat->Beta);
  printf("\t\t mu  %g \n", dat->mu);
#ifdef ANISOTROPY
  printf("\t\t tz (anisot.)  %g \n", dat->tz);
#endif
     

  printf("/*********************************************/\n\n");

}

void read_arguments( int argc, char *argv[],
		     struct s_data *data,
		     float *Beta, float *t, float *tz)
{
	int n,nbl,lbl;
	char nombres[100], name[100];
	FILE *Finput, *Fnombres;
	static struct s_data data_old;

	

	switch(argc)
	{
	 case 4: sscanf(argv[3],"%d", &nbl);
	 case 3: sscanf(argv[2],"%d", &n2);
		 sscanf(argv[1],"%d", &n1); break;
	 default: printf("ERROR: to start use\n ");
		  printf(" ana%d (1er fich) (ult. fich)",L);
		  printf(" (n_bloques)\n ");
		  exit(0);
	 }
        nblo=nbl;

        iteraciones = (n2-n1 + 1);
        
	 if(n2-n1+1>maxfiles)
	 {
	    printf("MAXIMUM NUMBER OF FILES = %i \n", maxfiles);
	    exit(0);
	 }

	 printf("\n STARTING \n");
	 printf(" ANALYSIS OF OUT%03d.DAT - OUT%03d.DAT\n",n1,n2);
	 printf(" LATTICE SIZE %dx%d \n",L,L);
	 printf(" JACK-KNIFE BLOCKS %d\n ",nblo);

	 if((nblo)>maxbloque)
	 {
		printf("\n nblo  no puede superar %d \n",maxbloque);
		exit(0);
	 }

        
	 lbl=(n2-n1+1)/nbl;
	 if(lbl==0)
		exit(0);
	 lblo=lbl;
	 n1=n2-lblo*nblo+1;

	 for(n=n1;n<=n2;n++)
	 {
		sprintf(name,"%s%03d.DAT",nom_fich,n);
		if ((Finput=fopen(name,"rb"))==NULL)
		{
			printf("\n THIS FILE DOES NOT EXIST \n",name);
			exit(0);
		}

		fread(data,sizeof(*data),1,Finput);
	
		if(n==n1)
		{
		   memcpy(&data_old,data, sizeof(data_old));
		}
		if(n>n1 && (data->Beta != data_old.Beta   ||
			    data->t != data_old.t   ||
			    data->tz != data_old.tz   ||
			    data->u != data_old.u   ||
			    data->itmax != data_old.itmax ||
			    data->mesfr != data_old.mesfr))
		{
			printf("NOT THE SAME INPUT FROM %d \n",n);
			exit (0);

		}
		fclose(Finput);

		if((ficheros[n-n1]= (char*)malloc(strlen(name)))==NULL)
		{
			printf("Fuera de memoria \n");
			exit(0);
		}
		strcpy(ficheros[n-n1],name);
	
          }

        printf("\n %d \n",data->itmax);
        describe_data(data);
        *Beta = data->Beta;
        *t = data->t;
	*tz = data->tz;

}


void read_data(int n)
{
    int idat,idatnew,imag,j,i,k;
    FILE *Finput;
    struct s_data data;
     
        
    
    Finput=fopen(ficheros[n],"rb");
    fread(&data,sizeof(data),1,Finput);

    fread(&SIGN,sizeof(int),1,Finput);
         
    for(idat=0;idat<n_obs_medid;idat++)	     	     
        fread(&v_dat[idat][0],4,(size_t)data.itmax,Finput);

    for(i=0;i<data.itmax;i++)       /*  Total Energy */
      v_dat[8][i] = v_dat[0][i] + v_dat[1][i];


    /*  for(idat=0;idat<n_obs_medid;idat++)
      for(i=0;i<data.itmax;i++)
	printf("v_dat[%d][%d] =%g\n",idat,i,v_dat[idat][i]);
    */

    fclose(Finput);
}


void Evolution(void)
{       
          FILE *Fevol[n_obs_plot];
          int n,it,obs;
          char nombreevol[n_obs_plot][150];


	  strcpy(nombreevol[0],"KineticE");	  
	  strcpy(nombreevol[1],"CoulombE");	  
	  strcpy(nombreevol[2],"Mferro");	  
	  strcpy(nombreevol[3],"Mstag");	  
	  strcpy(nombreevol[4],"ss1rst");	  
	  strcpy(nombreevol[5],"ss2nd");	  
	  strcpy(nombreevol[6],"Momentsqr");	  
	  strcpy(nombreevol[7],"Mplane");	  


          for(obs=0;obs<n_obs_FS;obs++)                 /* Energies */
          {  
            sprintf(nombreevol[obs],"%s.plt",nombreevol[obs]);
            if((Fevol[obs]=fopen(nombreevol[obs],"w"))==NULL)
              {
                puts("Could not open the result's file\n");
                exit(1);
              }

	

	    fprintf(Fevol[obs],"(set DEVICE postscript file=\n");
	    fprintf(Fevol[obs],"set font duplex\n");
	    fprintf(Fevol[obs],"title top '%s beta=%g U=%g N=%dx%dx%d'\n",
		    cadena[obs],data.Beta,data.u,L,L,L);
	    fprintf(Fevol[obs],"title bottom 'MC Iterations'\n");
	    fprintf(Fevol[obs],"title left '%s'\n",cadena[obs]);
	    fprintf(Fevol[obs],"set limit x %d\t%d\n",0,data.itmax*(n2-n1));
	  }


     
	  medidas = 1;

          for(n=n1;n<=n2;n++)
          {
           read_data(n-n1);

           for(it=0;it<data.itmax;it++)
           {      
              if((it%medidas) == 0)
              {
               for(obs=0;obs<n_obs_medid;obs++)
                {
                  fprintf(Fevol[obs],"%d\t%lf\n",
			  ((n-n1)*data.itmax+it),v_dat[obs][it]);
                }
              }
            }
          }

          for(obs=0;obs< n_obs_medid;obs++)
            {
              fprintf(Fevol[obs],"join dots red\n");
              fclose(Fevol[obs]);
            }
}         


void Correlation (int obs)
{
    int n,iter,itermax,it,i,id,maxdis;
    float ob[maxmed],cor[n_correl],corsum,sum1,sum2,med1[n_correl];
    float normal,med2[n_correl],r;
    FILE *Fcorrel,*FTau;
    char nombcorr[150],nomtau[100];
    double sum,mean,norm,tau;
    
    
    iter=0;
    for(n=n1;n<=n2;n++)
    {        
        read_data(n-n1);
        for(it=0;it<data.itmax;it++)
            ob[iter++]=v_dat[obs][it];
    }
    itermax=iter;

    sprintf(nombcorr,"correl%d_%5.4f_%5.4f_%d.plt",obs,data.Beta,data.u,L);

    Fcorrel=fopen(nombcorr,"w");

    fprintf(Fcorrel,"(set DEVICE postscript file=\n");
    fprintf(Fcorrel,"set font duplex\n");
    fprintf(Fcorrel,"title top 'Corr. in MC time  %s'\n",cadena[obs]);
    fprintf(Fcorrel,"title bottom 'MC Iterations'\n");
    fprintf(Fcorrel, "set symbol 4N\n");
    
    sprintf(nomtau,"tau%d_%5.4f_%5.4f_%d.plt",obs,data.Beta,data.u,L);
    FTau=fopen(nomtau,"w");
    fprintf(FTau,"title left font T '% tau'\n");
    fprintf(FTau,"title top font T 'Tau int, exp %s'\n",cadena[obs]);
    fprintf(FTau,"title bottom font T 'Iterations'\n");
    fprintf(FTau, "set symbol 4N\n");
                                                  

    for(it=0;it<n_correl;it++)
    {
        corsum=0.;
        sum1=0.;
        sum2=0.;
        for(n=0;n<(itermax-it);n++)
        {
            corsum+=ob[n]*ob[n+it];
            sum1+=ob[n];
            sum2+=ob[n+it];
        }
        cor[it]=corsum/(itermax-it);
        med1[it]=sum1/(itermax-it);
        med2[it]=sum2/(itermax-it);
        
        cor[it]-=med1[it]*med2[it];
    }
    normal=cor[0];
    for(n=0;n<n_correl;n++)
    {
        cor[n]=cor[n]/normal;
        fprintf(Fcorrel,"%d\t%f\n",n,cor[n]);
    }
    
    fprintf(Fcorrel,"plot \n");
    fprintf(Fcorrel,"join \n");
    fclose(Fcorrel);


    /* CALCULO AUTOCONSISTENTE DE TAU */
 
    for(i=0;i<itermax;i++)
      {
        sum +=ob[i];
   
      }

    mean = sum / itermax;
    for(id=0;id<n_correl;id++)
      {
        sum =0;
        for(i=0;i<itermax;i++)
	  sum += (ob[i] - mean)*(ob[i+id] -mean);
	

        if(!id) norm = sum;
        cor[id] = sum/norm;
        if( cor[id] <=0)
	  {
            id++;
            break;
	  }
      }
    maxdis = id;
    tau=0.5;
    for(id=1;id < maxdis; id++)
      {
        tau += cor[id];
        if(6*tau<id)
	  break;
      }
    
    for(id=1;id < maxdis; id++)
      {
	r=cor[id-1]/cor[id];
	if(r>1)
	  fprintf(FTau,"%d\t%lf\n",id,1.0/log(r));
      }
    
    fprintf(FTau,"plot\n");
    for(id=1;id < maxdis; id++)
      {
	r=cor[id];
	if(r>0)
	  fprintf(FTau,"%d\t%lf\n",id,(-1.0)*id/log(r));
      }                                                                
    fprintf(FTau,"plot\n");
    fclose(FTau);
     
    printf("Tau integrado de %d = %f\n",obs,tau);      


}  

void Histogram(int acoplo)
{
  int ib,j,sum,sumf2,sumf,i;
  float frecu[n_inter+1],y,mean,disp;
  FILE *Fhistog;
  char nomhist[150];
  double normalization;
  
  sprintf(nomhist,"h%d_%5.4f_%5.4f_%d.plt",acoplo,data.Beta,data.u,L);
  Fhistog=fopen(nomhist,"w");

  fprintf(Fhistog,"(set DEVICE postscript file=\n");
  fprintf(Fhistog,"set font duplex\n");
  fprintf(Fhistog,"title top 'Distribution of  %s'\n",cadena[acoplo]);
  fprintf(Fhistog,"title left 'frequency'\n");
  fprintf(Fhistog,"title bottom 'Spectrum of %s ' \n ",cadena[acoplo]);
  fprintf(Fhistog,"set order x y dummy \n");
  
  for(j=0;j<=n_inter;j++)
    frecu[j]=0.;
  for(j=0;j<=n_inter;j++)
    for(ib=0;ib<nblo;ib++)
      {
	frecu[j]+=frec[ib][j];
      }

  normalization=0;

  for(i=0;i<=n_inter;i++)
    for(j=0;j<nblo;j++)
	normalization+=(double) frec[j][i];

  normalization*=c1;
  
  for (i=0;i<=n_inter;i++)
    {
      sumf=sumf2=0;
      for(ib=0;ib<nblo;ib++)
	{
	  sumf+=frec[ib][i];
	  sumf2+=frec[ib][i]*frec[ib][i];
	}
      
      y=c1*i+c2;
      mean=(double) sumf/ (double) nblo;
      disp=sqrt((sumf2/(double) nblo-mean*mean)/(double) (nblo-1))*
	(double )nblo;
      fprintf(Fhistog,"%10.6f %10.6f %10.6f\n",y,
	      (double) sumf/normalization,disp/normalization);
    }
  
  fprintf(Fhistog,"join solid\n");
  fclose(Fhistog);
}

void FS(void)
{

  int ib,ibb,i,j,k,q,p,volu;
  double SumO[n_obs_FS],SumO2[n_obs_FS],SumderO[n_obs_FS],
    SumderO2[n_obs_FS],SumderlO[n_obs_FS],SumderlO2[n_obs_FS];
  double sumO[n_obs_FS],sumOe[n_obs_FS],sumx[n_obs_FS];
  double expo,expo_frec,O[n_obs_FS],derO[n_obs_FS],derlO[n_obs_FS];
  long int sumf;
  double en, sum,sume,x,y,lO[n_obs_FS];
  char name[100];
  

  for(j=0;j<=nbetas;j++)
    {
      for(k=0;k<n_obs_FS;k++)
	{
	  SumO[k]=0.;
	  SumO2[k]=0.;
	  SumderO[k]=0.;
	  SumderO2[k]=0.;
	  
	}

      for(ib=0;ib<=nblo;ib++)
	{
	  sum=sume=0.;
	  
	  for(k=0;k<n_obs_FS;k++)
	    {
	      sumO[k]=0.;
	      sumOe[k]=0.;
		  
	    }
	  
	  x= x0 - (delta*nbetas) + j*2*delta;
	  
                
	  for(i=0;i<n_inter;i++)
	    {
	      sumf=0.;
	      for(k=0;k<n_obs_FS;k++)
		sumx[k]=0.;
	      y=c1*i+c2;
	      for(ibb=0;ibb<nblo;ibb++)
		if(ib!=ibb)
		  {
		    sumf += frec[ibb][i];

		    for(k=0;k<n_obs_FS;k++)
		      sumx[k] += xfrec[k][ibb][i];
		  }
	      
	      expo=exp((x-x0)*(y-ymed)*vol);
	      expo_frec = expo*sumf;
	      sum += expo_frec;
	      sume += y*expo_frec;
	      
	      for(k=0;k<n_obs_FS;k++)
		{
		  sumO[k] += sumx[k]*expo;
		  sumOe[k] += sumx[k]*y*expo;
		}
	    }
	  
	  en=sume/sum;
	  for(k=0;k<n_obs_FS;k++)
	    {
	      O[k]=sumO[k]/sum;
	      derO[k]=(sumOe[k]/sum - O[k]*en)*vol;
	      
	    }
	  
	  
	  if(ib<nblo)
	    for(k=0;k<n_obs_FS;k++)
	      {
		SumO[k] += O[k];
		SumO2[k] += O[k]*O[k];
		SumderO[k] += derO[k];
		SumderO2[k] += derO[k]*derO[k];
		
	      }
	  
	  for(k=0;k<n_obs_FS;k++)
	    if(fabs(derO[k]) > fabs(maxder[k][ib]))
	      { 
		
		maxder[k][ib]=derO[k];
		coup_maxder[k][ib]=x;
	      }

	  
	} /* nblo */
      for(k=0;k<n_obs_FS;k++)
	{
	  SumO[k]/=nblo;
	  SumderO[k]/=nblo;
	  SumderlO[k]/=nblo;
	  O_v[k][j]=O[k];
	  derO_v[k][j]=derO[k];
	  O_err[k][j]=
	    sqrt(fabs(SumO2[k]/nblo-SumO[k]*SumO[k])*
		 (double)(nblo-1));
	  
	  derO_err[k][j]=
	    sqrt(fabs(SumderO2[k]/nblo -
		      SumderO[k]*SumderO[k])*(double)(nblo-1));

	  
	  
	}


             
    } /* bucle en betas */


  


}


void titulo (void)
{
  int i;
  sprintf(nombre,"title top font T 'L=%d, beta=%f, U=%f, #MC= %dx%d Nbl=%d, Lblo=%d, %s'",L,data.Beta,data.u,data.itmax,data.mesfr,nblo, lblo, nom);


}



void Dibuja_resultados (int v)
{
  int j;
  double x;
  char formato[]= "new plot \n"
    "%s \n"
    "title left font T '%s' \n"
    "title right font T '%s' \n"
    "case               '%s' \n"
    "set order x y dy \n"
    "set labels left font T \n"
    "set labels bottom font T \n"
    "set symbol %s \n";
  
  
  /* OBSERVABLE  */
  
  fprintf(Foutplt,formato,nombre,"O",n_coup,"G- -","9N");
  
  for(j=0;j<=nbetas;j++)
    {
      x=x0- (nbetas*delta) + (2*j*delta);
      fprintf(Foutplt,"%10f  %10f\n",x,O_v[v][j]);
    }
  fprintf(Foutplt,"join \n");
  
  for(j=0;j<=nbetas;j++)
    {
      
      x=x0-(nbetas*delta) + (2*j*delta);
      fprintf(Foutplt,"%10f  %10f  %10f\n",x,O_v[v][j], O_err[v][j]);
    }
  fprintf(Foutplt,"plot\n");
  fprintf(Foutplt,"set symbol 5N \n");
  fprintf(Foutplt, " %10f %10f %10f\n", x0,O_v[v][nbetas/2],
	  O_err[v][nbetas/2]);
  fprintf(Foutplt, "plot \n");
  
  /* DERIVADA  */
  fprintf(Foutplt, formato, nombre,n_deri,n_coup,"G- -","9N");
  
  for(j=0;j<=nbetas;j++)
    {
      x=x0- (nbetas*delta) + (2*j*delta);
      fprintf(Foutplt,"%10f  %10f\n",x,derO_v[v][j]);
    }
  fprintf(Foutplt,"join \n");
	
  for(j=0;j<=nbetas;j++)
    {
      
      x=x0-(nbetas*delta) + (2*j*delta);
      fprintf(Foutplt,"%10f  %10f  %10f\n",
	      x,derO_v[v][j], derO_err[v][j]);
    }
  
  
  fprintf(Foutplt,"plot\n");
  fprintf(Foutplt,"set symbol 3N \n");
  fprintf(Foutplt,"set order x dx y dy \n");
  fprintf(Foutplt,"%10f %10f %10f %10f\n",coup_maxder[v][nblo],
	  err_coup[v],maxder[v][nblo],err_maxder[v]);
  fprintf(Foutplt,"set order x y dy \n");
  fprintf(Foutplt, " %10f %10f %10f \n", x0,derO_v[v][nbetas/2],
	  derO_err[v][nbetas/2]);
  fprintf(Foutplt, "plot \n");
  
  fprintf(Foutplt,"title bottom font T' C=%8.6f+/-%8.6f; D=%8f+/-%8f'\n",
	  coup_maxder[v][nblo],err_coup[v],
	  maxder[v][nblo],err_maxder[v]);
  

}






int main(int argc,char *argv[])
{

  int n,i,ib,il,j,d;
  float e_min,e_max,energy;
  FILE *fevol;

  char nombreresult[150],n_der[3][10],n_cou[3][6];
  
  char nomeigen[50],nombrefichero[150];

  int  n_vo[3],it,p,Num_Ener,maxi,maxj,tt;
  float data_beta[3],sum,sumen,sumen2,ymin,ymax,beta,gamma,coup;
  float en,cv,sigma,y,aut[4],max,s_prom[4],promedio[4],desviacion[4];
  FILE *Foutput;
  int k,mm,kk,itt;
    
    
  strcpy(n_cou[0],"beta");
  strcpy(n_cou[1],"U");
    

  read_arguments(argc,argv,&data,&Beta,&t,&tz);
    
  Evolution();


  printf("Analyzed Observables:\n");
  for(j=0;j< (n_obs_medid);j++)
    printf("          %d\t%s\n",j,cadena[j]);



  for(d=0;d<2;d++)
    {
      sprintf(nombreresult,"%dresult_%5.4f_%5.4f_%d.plt",
	      d,data.Beta,data.u,L);
      if((Foutplt=fopen(nombreresult,"w"))==NULL)
        {
	  puts("Could not open the result file\n");
	  exit(1);
        }


      strcpy(n_coup,n_cou[d]);

      data_beta[0]=data.Beta;
      data_beta[1]=data.u;

      coup=data_beta[d];
      

      sum=sumen=sumen2=0.;
      ymin=2.;               /* Valores maximo y minimo de E */
      ymax=-2.0;
      
      for(n=n1;n<=n2;n++)
        {
	  read_data(n-n1);
	  for(it=0;it<data.itmax;it++)
            {
	      y=v_dat[d][it];
	      sum++;
	      sumen += y;
	      sumen2 += y*y;
	      ymin=(y<ymin)?y:ymin;
	      ymax=(y>ymax)?y:ymax;				
            }
        }
      en=sumen/sum;                  /* Dispersion de la gaussiana */
      cv=fabs(sumen2/sum-en*en);     /* para hacer FS              */

      sigma=sqrt(cv);

      for(n=0;n<n_obs_FS;n++)
	for(ib=0;ib<nblo;ib++)
	  for(j=0;j<=n_inter;j++)
	    {
	      xfrec[n][ib][j]=frec[ib][j]=0.;
	    }


      h=1.0/(ymax-ymin);
      c1=1.0/(float)n_inter/h;
      c2=-0.5/(float)n_inter/h+ymin;
      
      for(n=n1;n<=n2;n++)
	
        {                            /* Clasificacion en intervalos */
	  ib=(n-n1)/lblo;            /* de los Observ. */
	  read_data(n-n1);
	  for(it=0;it<data.itmax;it++)
            {
	      
	      y=v_dat[d][it];

	      i=(int)((y-ymin)*h*n_inter);

	      frec[ib][i]++;
	      
	      for(tt=0;tt<n_obs_FS;tt++)
		xfrec[tt][ib][i]+=v_dat[tt][it];
	      
            }     
        }

      delta = 0.005;

      ymed=en;
      x0=coup;


      Histogram(d);
      FS();

      for(n=0;n<(n_obs_FS );n++)
	{
	  
	  sprintf(nom," O= %s ",cadena[n]);
	  titulo();
	  Dibuja_resultados(n);
	  /* Correlation(n); */

	}

    }
  


  
  for(n=0;n< (n_obs_FS  );n++)
    {
      printf("\n %d. %s= %8.6f +/- %8.6f \n",n,cadena[n],
	     O_v[n][nbetas/2],O_err[n][nbetas/2]);
      sprintf(nombrefichero,"%dvsb_%d.plt",n,L);

      if((Foutput=fopen(nombrefichero,"a"))==NULL)
        {
	  puts("\n No puedo abrir el fichero vsb.plt ");
	  exit(1);
        }
      
      fprintf(Foutput,"%f %f\t%f\t%f\n",data.Beta,data.u,O_v[n][nbetas/2],
	      O_err[n][nbetas/2]);
      fclose(Foutput);
    }
  







  return 0;
  
}          /*  END */
	
