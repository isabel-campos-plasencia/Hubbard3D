/*************************************************************
*   Matrix inversion funcions from Numerical Recipes         *
*                                                            *
*                                                            *
*************************************************************/
#include "hm3d.h"
#define NRANSI
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

#define TINY 1.0e-20;

void ludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=dvector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}

	for (j=1;j<=n;j++) {

		for (i=1;i<j;i++) {

			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;

		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_dvector(vv,1,n);
}

#undef TINY
#undef NRANSI   

void lubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}






void jacobi(double **a, int n, double d[], double **v,int *nrot)
{
  int j,iq,ip,i;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
  
  b=dvector(1,n);
  z=dvector(1,n);
  for (ip=1;ip<=n;ip++) {
    for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=1;ip<=n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;

  for (i=1;i<=500;i++) 
    {
      /* printf("Jacobi iteration number =%d\n",i);        */
      sm=0.0;
      for (ip=1;ip<=n-1;ip++) {
	for (iq=ip+1;iq<=n;iq++)
	  sm += fabs(a[ip][iq]);
      }
      if (sm == 0.0) {
	free_dvector(z,1,n);
	free_dvector(b,1,n);
	return;
      }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++) {
	g=100.0*fabs(a[ip][iq]);
	if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
	    && (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
	  a[ip][iq]=0.0;
	else if (fabs(a[ip][iq]) > tresh) {
	  h=d[iq]-d[ip];
	  if ((float)(fabs(h)+g) == (float)fabs(h))
	    t=(a[ip][iq])/h;
	  else {
	    theta=0.5*h/(a[ip][iq]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq]=0.0;
	  for (j=1;j<=ip-1;j++) {
	    ROTATE(a,j,ip,j,iq)
	      }
	  for (j=ip+1;j<=iq-1;j++) {
	    ROTATE(a,ip,j,j,iq)
	      }
	  for (j=iq+1;j<=n;j++) {
	    ROTATE(a,ip,j,iq,j)
	      }
	  for (j=1;j<=n;j++) {
	    ROTATE(v,j,ip,j,iq)
	      }
	  ++(*nrot);
	}
      }
    }
    for (ip=1;ip<=n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  nrerror("Too many iterations in routine jacobi");
}
#undef ROTATE


void tred2(double **a, int n, double d[], double e[])
{
  int l,k,j,i;
  double scale,hh,h,g,f;
  
  for (i=n;i>=2;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 1) {
      for (k=1;k<=l;k++)
	scale += fabs(a[i][k]);
      if (scale == 0.0)
	e[i]=a[i][l];
      else {
	for (k=1;k<=l;k++) {
	  a[i][k] /= scale;
	  h += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
	e[i]=scale*g;
	h -= f*g;
	a[i][l]=f-g;
	f=0.0;
	for (j=1;j<=l;j++) {
	  a[j][i]=a[i][j]/h;
	  g=0.0;
	  for (k=1;k<=j;k++)
	    g += a[j][k]*a[i][k];
	  for (k=j+1;k<=l;k++)
	    g += a[k][j]*a[i][k];
	  e[j]=g/h;
	  f += e[j]*a[i][j];
	}
	hh=f/(h+h);
	for (j=1;j<=l;j++) {
	  f=a[i][j];
	  e[j]=g=e[j]-hh*f;
	  for (k=1;k<=j;k++)
	    a[j][k] -= (f*e[k]+g*a[i][k]);
				}
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  d[1]=0.0;
	e[1]=0.0;
	/* Contents of this loop can be omitted if eigenvectors not
			wanted except for statement d[i]=a[i][i]; */
	for (i=1;i<=n;i++) {
	  l=i-1;
	  if (d[i]) {
			for (j=1;j<=l;j++) {
			  g=0.0;
			  for (k=1;k<=l;k++)
			    g += a[i][k]*a[k][j];
			  for (k=1;k<=l;k++)
			    a[k][j] -= g*a[k][i];
			}
		}
	  d[i]=a[i][i];
	  a[i][i]=1.0;
	  for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
	}
}


double dpythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+DSQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+DSQR(absa/absb)));
}

void tqli(double d[], double e[], int n, double **z)
{
  double dpythag(double a, double b);
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;
  
  for (i=2;i<=n;i++) e[i-1]=e[i];
  e[n]=0.0;
  for (l=1;l<=n;l++) {
    iter=0;
    do 
      {
	for (m=l;m<=n-1;m++) 
	  {
	    dd=fabs(d[m])+fabs(d[m+1]);
	    if ((float)(fabs(e[m])+dd) == dd) break;
	  }
	if (m != l) 
	  {
	    printf("iter = %d\n",iter);
	    if (iter++ == 30) nrerror("Too many iterations in tqli");
	    g=(d[l+1]-d[l])/(2.0*e[l]);
	    r=dpythag(g,1.0);
	    g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
	    s=c=1.0;
	    p=0.0;
	    for (i=m-1;i>=l;i--) 
	      {
		f=s*e[i];
		b=c*e[i];
		e[i+1]=(r=dpythag(f,g));
		if (r == 0.0) 
		  {
		    d[i+1] -= p;
		    e[m]=0.0;
		    break;
		  }
		s=f/r;
		c=g/r;
		g=d[i+1]-p;
		r=(d[i]-g)*s+2.0*c*b;
		d[i+1]=g+(p=s*r);
		g=c*r-b;
		for (k=1;k<=n;k++) 
		  {
		    f=z[k][i+1];
		    z[k][i+1]=s*z[k][i]+c*f;
		    z[k][i]=c*z[k][i]-s*f;
		  }
	      }
	    if (r == 0.0 && i >= l) continue;
	    d[l] -= p;
	    e[l]=g;
	    e[m]=0.0;
	  }
      } while (m != l);
  }
}





