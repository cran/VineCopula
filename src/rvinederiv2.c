#include "include/vine.h"
#include "include/memoryhandling.h"
#include "include/likelihood.h"
#include "include/deriv.h"
#include "include/deriv2.h"
#include "include/tcopuladeriv.h"
#include "include/tcopuladeriv_new.h"
#include "include/rvine.h"
#include "include/rvinederiv.h"
#include "include/rvinederiv2.h"


#define UMAX  1-1e-10

#define UMIN  1e-10

#define XEPS 1e-4


//////////////////////////////////////////////////////////////
// Function to compute the second derivative of log-likelihood for the pair-copula construction (Rvine)
// (by J.S. and U.S.)
//
// Input:
// T				sample size
// d				dimension (>=2)
// family			copula families: only student //  (1=gaussian, 2=student, 3=clayton, 4=gumbel, 5=frank, 6=joe)
// kk				row number of the parameter wrt to which we want to calculate the derivative
// ii				column number of the parameter wrt to which we want to calculate the derivative
// kkk				row number of the parameter wrt to which we want to calculate the derivative (second parameter)
// iii				column number of the parameter wrt to which we want to calculate the derivative (second parameter)
// maxmat			Maximums-Matrix from RVM object
// matrix			an RVineMatrix in vector form
// condirect, conindirect	Matrizes which tell us where we find the right values 
// par				parameter values (at least d*(d-1)/2 parameters
// par2				second set of parameter values (f.e. for student copulas)
// data				data set for which to compute the derivative of the log-likelihood
// calcupdate, calcupdate2	matrix which tells which terns we need to consider for the calculation of the derivative
// ll_tilde
// vv_tilde, vv2_tilde
// ll_hat
// vv_hat, vv2_hat
// ll
// vv, vv2
// tcop				variable for the t-copula
//
// Output:
// out					second deriv of Loglikelihood
// barvalue_out				array with the contribution to the derivative (for each copula)
// barvdirect_out, barvindirect_out	array for the derivatives of the h-functions  
/////////////////////////////////////////////////////////////


void VineLogLikRvineDeriv2(int* T, int* d, int* family, int* kk, int* ii, int* kkk, int* iii, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data, 
						  double* ll_tilde, double* vv_tilde, double* vv2_tilde, double* ll_hat, double* vv_hat, double* vv2_hat, int* calcupdate, int* calcupdate2,
						  double* out, double* ll, double* vv, double* vv2, double* barvalue_out, double* barvdirect_out, double* barvindirect_out, int* tcop, int* kk_second)
{
	int i, j, k, t, m, **fam, **calc, **calc2, a=1;
	
	double sumloglik=0.0, **theta, **nu, *zr1, *zr2, *tildezr1, *tildezr2, *cop;
	double *handle1, *helpvar, *helpvar2, *helpvar3, *helpvar4, *helpvar5, *helpvar6;
	double *tildezr1_hatk_hati, *tildezr1_tildek_tildei, *barzr1, *tildezr2_hatk_hati, *tildezr2_tildek_tildei, *barzr2;
	double *der, *der_hatk_hati, *der_tildek_tildei, param[2];

	//Rprintf("kk: %d\t ii: %d\t kkk: %d\t iii: %d\t tcop: %d\t kk_second: %d\n", *kk, *ii,*kkk,*iii,*tcop,*kk_second);
	// Zielvariablen
	double ***barvdirect, ***barvindirect;

	
	zr1=(double*) Calloc(*T,double);
	zr2=(double*) Calloc(*T,double);
	handle1=(double*) Calloc(*T,double);
	helpvar=(double*) Calloc(*T,double);
	helpvar2=(double*) Calloc(*T,double);
	helpvar3=(double*) Calloc(*T,double);
	helpvar4=(double*) Calloc(*T,double);
	helpvar5=(double*) Calloc(*T,double);
	helpvar6=(double*) Calloc(*T,double);
	tildezr1=(double*) Calloc(*T,double);
	tildezr2=(double*) Calloc(*T,double);
	tildezr1_hatk_hati=(double*) Calloc(*T,double);
	tildezr1_tildek_tildei=(double*) Calloc(*T,double);
	barzr1=(double*) Calloc(*T,double);
	tildezr2_hatk_hati=(double*) Calloc(*T,double);
	tildezr2_tildek_tildei=(double*) Calloc(*T,double);
	barzr2=(double*) Calloc(*T,double);
	cop=(double*) Calloc(*T,double);
	der=(double*) Calloc(*T,double);
	der_hatk_hati=(double*) Calloc(*T,double);
	der_tildek_tildei=(double*) Calloc(*T,double);
	
	
	
	//Allocate memory
	//value=create_3darray(*d,*d,*T);
	//vdirect = create_3darray(*d,*d,*T);
	//vindirect = create_3darray(*d,*d,*T);
	//tildevalue=create_3darray(*d,*d,*T);
	//tildevdirect = create_3darray(*d,*d,*T);
	//tildevindirect = create_3darray(*d,*d,*T);
	//tildevalue2=create_3darray(*d,*d,*T);
	//tildevdirect2 = create_3darray(*d,*d,*T);
	//tildevindirect2 = create_3darray(*d,*d,*T);
	//barvalue=create_3darray(*d,*d,*T);
	barvdirect = create_3darray(*d,*d,*T);
	barvindirect = create_3darray(*d,*d,*T);
	theta=create_matrix(*d,*d);
	nu=create_matrix(*d,*d);
	fam=create_intmatrix(*d,*d);
	//mmat=create_intmatrix(*d,*d);
	//mat=create_intmatrix(*d,*d);
	calc=create_intmatrix(*d,*d);
	calc2=create_intmatrix(*d,*d);
	
	// Hilfsvariable für die t-copula
	//param=(double*) Calloc(2,double);

	k=0;
	for(i=0;i<(*d);i++)
	{
        for(j=0;j<(*d);j++)
		{
			theta[i][j]=par[(i+1)+(*d)*j-1] ;
			nu[i][j]=par2[(i+1)+(*d)*j-1]    ;
			//mmat[i][j]=maxmat[(i+1)+(*d)*j-1] ;
			//mat[i][j]=matrix[(i+1)+(*d)*j-1] ;
			fam[i][j]=family[(i+1)+(*d)*j-1] ;
			calc[i][j]=calcupdate[(i+1)+(*d)*j-1] ;
			calc2[i][j]=calcupdate2[(i+1)+(*d)*j-1] ;
		}

	}       
	

 
	for(i=0;i<(*d);i++)
	{
        for(j=0;j<(*d);j++){	
			for(t=0;t<(*T);t++ ) {
				// Log-lik
				//vdirect[i][j][t]=vv[(i+1)+(*d)*j+(*d)*(*d)*t-1];
				//vindirect[i][j][t]=vv2[(i+1)+(*d)*j+(*d)*(*d)*t-1];
				//value[i][j][t]=ll[(i+1)+(*d)*j+(*d)*(*d)*t-1];
				// deriv Log-lik
				//tildevdirect[i][j][t]=vv_tilde[(i+1)+(*d)*j+(*d)*(*d)*t-1];
				//tildevindirect[i][j][t]=vv2_tilde[(i+1)+(*d)*j+(*d)*(*d)*t-1];
				//tildevalue[i][j][t]=ll_tilde[(i+1)+(*d)*j+(*d)*(*d)*t-1];
				//tildevdirect2[i][j][t]=vv_hat[(i+1)+(*d)*j+(*d)*(*d)*t-1];
				//tildevindirect2[i][j][t]=vv2_hat[(i+1)+(*d)*j+(*d)*(*d)*t-1];
				//tildevalue2[i][j][t]=ll_hat[(i+1)+(*d)*j+(*d)*(*d)*t-1];
				// Zielvariablen
				barvdirect[i][j][t]=0;
				barvindirect[i][j][t]=0;
				//barvalue[i][j][t]=0;
				barvalue_out[(i+1)+(*d)*j+(*d)*(*d)*t-1]=0;
			}}}
	
	
	if(calc2[(*kk-1)][(*ii-1)]==1)
	{
		//m=mmat[(*kk-1)][*ii-1];
		m=maxmat[(*kk)+(*d)*(*ii-1)-1];
		//k0 = *kk-1;
		for(t=0;t<*T;t++ ) 
		{
			//zr1[t]=vdirect[(*kk-1)][(*ii-1)][t];
			zr1[t]=vv[(*kk)+(*d)*(*ii-1)+(*d)*(*d)*t-1];
			//tildezr1[t]=tildevdirect2[(*kk-1)][(*ii-1)][t];
			tildezr1[t]=vv_hat[(*kk)+(*d)*(*ii-1)+(*d)*(*d)*t-1];
			//cop[t]=exp(value[(*kk-1)][(*ii-1)][t]);
			cop[t]=exp(ll[(*kk)+(*d)*(*ii-1)+(*d)*(*d)*t-1]);
			//der[t]=tildevalue2[(*kk-1)][(*ii-1)][t];
			der[t]=ll_hat[(*kk)+(*d)*(*ii-1)+(*d)*(*d)*t-1];
		}
		
		//if(m == mat[*kk-1][*ii-1])
		if(m == matrix[*kk+(*d)*(*ii-1)-1])
		{	
			for(t=0;t<*T;t++ ) {
				//zr2[t]=vdirect[(*kk-1)][(*d-m)][t];
				zr2[t]=vv[(*kk)+(*d)*(*d-m)+(*d)*(*d)*t-1];
				//tildezr2[t]=tildevdirect2[(*kk-1)][(*d-m)][t];
				tildezr2[t]=vv_hat[(*kk)+(*d)*(*d-m)+(*d)*(*d)*t-1];
			}
			
		}else {
			for(t=0;t<*T;t++)
			{
				//zr2[t]=vindirect[(*kk-1)][(*d-m)][t];
				zr2[t]=vv2[(*kk)+(*d)*(*d-m)+(*d)*(*d)*t-1];
				//tildezr2[t]=tildevindirect2[(*kk-1)][(*d-m)][t];
				tildezr2[t]=vv2_hat[(*kk)+(*d)*(*d-m)+(*d)*(*d)*t-1];
			}	
		}
		if(*kk==*kkk && *ii==*iii)
		{
			param[0]=theta[*kk-1][*ii-1];
			param[1]=nu[*kk-1][*ii-1];
			if(*tcop==1)		//Für die t-copula
			{
				diff2hfunc_rho_tCopula(zr1,zr2,T,param,&fam[*kk-1][*ii-1],barvdirect[*kk-2][(*ii-1)]);
				diff2hfunc_rho_tCopula(zr2,zr1,T,param,&fam[*kk-1][*ii-1],barvindirect[*kk-2][(*ii-1)]);
				diff2PDF_rho_tCopula(zr1,zr2,T,param,&fam[*kk-1][*ii-1],helpvar);
			}
			else if(*tcop==2)
			{
				diff2hfunc_nu_tCopula_new(zr1,zr2,T,param,&fam[*kk-1][*ii-1],barvdirect[*kk-2][(*ii-1)]);
				diff2hfunc_nu_tCopula_new(zr2,zr1,T,param,&fam[*kk-1][*ii-1],barvindirect[*kk-2][(*ii-1)]);
				diff2PDF_nu_tCopula_new(zr1,zr2,T,param,&fam[*kk-1][*ii-1],helpvar);
			}
			else if(*tcop==3)
			{
				diff2hfunc_rho_nu_tCopula_new(zr1,zr2,T,param,&fam[*kk-1][*ii-1],barvdirect[*kk-2][(*ii-1)]);
				diff2hfunc_rho_nu_tCopula_new(zr2,zr1,T,param,&fam[*kk-1][*ii-1],barvindirect[*kk-2][(*ii-1)]);
				diff2PDF_rho_nu_tCopula_new(zr1,zr2,T,param,&fam[*kk-1][*ii-1],helpvar);
			}
			else
			{
				diff2hfunc_mod(zr1,zr2,T,&theta[*kk-1][*ii-1],&fam[*kk-1][*ii-1],barvdirect[*kk-2][(*ii-1)]);
				diff2hfunc_mod(zr2,zr1,T,&theta[*kk-1][*ii-1],&fam[*kk-1][*ii-1],barvindirect[*kk-2][(*ii-1)]);
				diff2PDF_mod(zr1,zr2,T,&theta[*kk-1][*ii-1],&fam[*kk-1][*ii-1],helpvar);
			}
			
			for(t=0;t<*T;t++)
			{
				if(*tcop==3)
				{
					//barvalue[(*kk-1)][(*ii-1)][t] = helpvar[t]/cop[t] - der[t]*tildevalue[(*kk-1)][(*ii-1)][t];
					barvalue_out[(*kk)+(*d)*(*ii-1)+(*d)*(*d)*t-1] = helpvar[t]/cop[t] - der[t]*ll_tilde[(*kk)+(*d)*(*ii-1)+(*d)*(*d)*t-1];
				}
				else
				{
					barvalue_out[(*kk)+(*d)*(*ii-1)+(*d)*(*d)*t-1] = helpvar[t]/cop[t] - der[t]*der[t];
				}
			}
		}
		else
		{
		if(calc2[(*kk)][(*ii-1)]==1)
		{
			param[0]=theta[*kk-1][*ii-1];
			param[1]=nu[*kk-1][*ii-1];
			if(*tcop==1 || (*tcop==3 && *kk_second==0))		//Für die t-copula
			{
				diffPDF_rho_tCopula(zr1,zr2,T,param,&fam[*kk-1][*ii-1],helpvar);
				diff2PDF_rho_u_tCopula_new(zr1,zr2,T,param,&fam[*kk-1][*ii-1],helpvar2);
				diff2hfunc_rho_v_tCopula_new(zr2,zr1,T,param,&fam[*kk-1][*ii-1],barvindirect[*kk-2][(*ii-1)]);
			}
			else if(*tcop==2 || (*tcop==3 && *kk_second==1))
			{
				diffPDF_nu_tCopula_new(zr1,zr2,T,param,&fam[*kk-1][*ii-1],helpvar);
				diff2PDF_nu_u_tCopula_new(zr1,zr2,T,param,&fam[*kk-1][*ii-1],helpvar2);
				diff2hfunc_nu_v_tCopula_new(zr2,zr1,T,param,&fam[*kk-1][*ii-1],barvindirect[*kk-2][(*ii-1)]);
			}
			else
			{
				diffPDF_mod(zr1,zr2,T,&theta[*kk-1][*ii-1],&fam[*kk-1][*ii-1],helpvar);
				diff2PDF_par_u_mod(zr1,zr2,T,&theta[*kk-1][*ii-1],&fam[*kk-1][*ii-1],helpvar2);
				diff2hfunc_par_v_mod2(zr2,zr1,T,&theta[*kk-1][*ii-1],&fam[*kk-1][*ii-1],barvindirect[*kk-2][(*ii-1)]);
			}
			//Rprintf("helpvar: %f",helpvar[0]);
			//Rprintf("helpvar2: %f",helpvar2[0]);
			//Rprintf("tildezr1: %f",tildezr1[0]);
			//Rprintf("cop: %f",cop[0]);
			//Rprintf("der: %f",der[0]);

			for(t=0;t<*T;t++)
			{
				barvalue_out[(*kk)+(*d)*(*ii-1)+(*d)*(*d)*t-1]=der[t]*(-helpvar[t])/cop[t] + helpvar2[t]/cop[t]*tildezr1[t];
				barvdirect[(*kk-2)][(*ii-1)][t]=helpvar[t]*tildezr1[t];
				barvindirect[(*kk-2)][(*ii-1)][t]=barvindirect[(*kk-2)][(*ii-1)][t]*tildezr1[t];
			}
		}
		if(calc2[*kk][*d-m]==1)
		{
			param[0]=theta[*kk-1][*ii-1];
			param[1]=nu[*kk-1][*ii-1];
			if(*tcop==1 || (*tcop==3 && *kk_second==0))		//Für die t-copula
			{
				diff2PDF_rho_u_tCopula_new(zr2,zr1,T,param,&fam[*kk-1][*ii-1],helpvar);
			}
			else if(*tcop==2 || (*tcop==3 && *kk_second==1))
			{
				diff2PDF_nu_u_tCopula_new(zr2,zr1,T,param,&fam[*kk-1][*ii-1],helpvar);
			}
			else
			{
				diff2PDF_par_v_mod(zr1,zr2,T,&theta[*kk-1][*ii-1],&fam[*kk-1][*ii-1],helpvar);
			}
			//Rprintf("helpvar: %f \n",helpvar[0]);
			//Rprintf("tildezr2: %f \n",tildezr2[0]);
			//Rprintf("cop: %f \n",cop[0]);
			for(t=0;t<*T;t++)
			{
				barvalue_out[(*kk)+(*d)*(*ii-1)+(*d)*(*d)*t-1] = barvalue_out[(*kk)+(*d)*(*ii-1)+(*d)*(*d)*t-1] + helpvar[t]/cop[t]*tildezr2[t];
			}
			if(calc2[*kk][*ii-1]==0)
			{
				if(*tcop==1 || (*tcop==3 && *kk_second==0))
				{
					diffPDF_rho_tCopula(zr1,zr2,T,param,&fam[*kk-1][*ii-1],helpvar);
				}
				else if(*tcop==2 || (*tcop==3 && *kk_second==1))
				{
					diffPDF_nu_tCopula_new(zr1,zr2,T,param,&fam[*kk-1][*ii-1],helpvar);
				}
				else
				{
					diffPDF_mod(zr1,zr2,T,&theta[*kk-1][*ii-1],&fam[*kk-1][*ii-1],helpvar);
				}
				//Rprintf("helpvar: %f \n",helpvar[0]);
				//Rprintf("der: %f \n",der[0]);
				for(t=0;t<*T;t++)
				{
					barvalue_out[(*kk)+(*d)*(*ii-1)+(*d)*(*d)*t-1] = barvalue_out[(*kk)+(*d)*(*ii-1)+(*d)*(*d)*t-1] + der[t]*(-helpvar[t])/cop[t];
				}
			}

			if(*tcop==1 || (*tcop==3 && *kk_second==0))
			{
				diff2hfunc_rho_v_tCopula_new(zr1,zr2,T,param,&fam[*kk-1][*ii-1],helpvar);
				diffPDF_rho_tCopula(zr2,zr1,T,param,&fam[*kk-1][*ii-1],helpvar2);
			}
			else if(*tcop==2 || (*tcop==3 && *kk_second==1))
			{
				diff2hfunc_nu_v_tCopula_new(zr1,zr2,T,param,&fam[*kk-1][*ii-1],helpvar);
				diffPDF_nu_tCopula_new(zr2,zr1,T,param,&fam[*kk-1][*ii-1],helpvar2);
			}
			else
			{
				diff2hfunc_par_v_mod(zr1,zr2,T,&theta[*kk-1][*ii-1],&fam[*kk-1][*ii-1],helpvar);
				diffPDF_mod(zr2,zr1,T,&theta[*kk-1][*ii-1],&fam[*kk-1][*ii-1],helpvar2);
			}


			for(t=0;t<*T;t++)
			{
				barvdirect[*kk-2][*ii-1][t] = barvdirect[*kk-2][*ii-1][t] + helpvar[t]*tildezr2[t];
				barvindirect[*kk-2][*ii-1][t] = barvindirect[*kk-2][*ii-1][t] + helpvar2[t]*tildezr2[t];
			}
		}
		}
		
	}
	//else
	//{
	//	k0=*kk-1;
	//}

	for(i=((*ii)-1); i>-1; i--)
	{
		//if(i<(*ii-1))
		//{
		//	k0=*kk-1;
		//}
		for(k=(*kk-2); k>i;k--)
		{
			//m=mmat[k][i];
			m=maxmat[(k+1)+(*d)*i-1];
			for(t=0;t<*T;t++ ) 
			{
				//zr1[t]=vdirect[k][i][t];
				zr1[t]=vv[(k+1)+(*d)*i+(*d)*(*d)*t-1];
				//tildezr1_hatk_hati[t]=tildevdirect2[k][i][t];
				tildezr1_hatk_hati[t]=vv_hat[(k+1)+(*d)*i+(*d)*(*d)*t-1];
				//tildezr1_tildek_tildei[t]=tildevdirect[k][i][t];
				tildezr1_tildek_tildei[t]=vv_tilde[(k+1)+(*d)*i+(*d)*(*d)*t-1];
				barzr1[t]=barvdirect[k][i][t];
				//cop[t]=exp(value[k][i][t]);
				cop[t]=exp(ll[(k+1)+(*d)*i+(*d)*(*d)*t-1]);
				//der_hatk_hati[t]=tildevalue2[k][i][t];
				der_hatk_hati[t]=ll_hat[(k+1)+(*d)*i+(*d)*(*d)*t-1];
				//der_tildek_tildei[t]=tildevalue[k][i][t];
				der_tildek_tildei[t]=ll_tilde[(k+1)+(*d)*i+(*d)*(*d)*t-1];
			}
			
			//if(m==mat[k][i])
			if(m==matrix[(k+1)+(*d)*i-1])
			{
				for(t=0;t<*T;t++ ) 
				{
					//zr2[t]=vdirect[k][(*d-m)][t];
					zr2[t]=vv[(k+1)+(*d)*(*d-m)+(*d)*(*d)*t-1];
					//tildezr2_hatk_hati[t]=tildevdirect2[k][(*d-m)][t];
					tildezr2_hatk_hati[t]=vv_hat[(k+1)+(*d)*(*d-m)+(*d)*(*d)*t-1];
					//tildezr2_tildek_tildei[t]=tildevdirect[k][(*d-m)][t];
					tildezr2_tildek_tildei[t]=vv_tilde[(k+1)+(*d)*(*d-m)+(*d)*(*d)*t-1];
					barzr2[t]=barvdirect[k][(*d-m)][t];
				}
			}
			else
			{
				for(t=0;t<*T;t++ ) 
				{
					//zr2[t]=vindirect[k][(*d-m)][t];
					zr2[t]=vv2[(k+1)+(*d)*(*d-m)+(*d)*(*d)*t-1];
					//tildezr2_hatk_hati[t]=tildevindirect2[k][(*d-m)][t];
					tildezr2_hatk_hati[t]=vv2_hat[(k+1)+(*d)*(*d-m)+(*d)*(*d)*t-1];
					//tildezr2_tildek_tildei[t]=tildevindirect[k][(*d-m)][t];
					tildezr2_tildek_tildei[t]=vv2_tilde[(k+1)+(*d)*(*d-m)+(*d)*(*d)*t-1];
					barzr2[t]=barvindirect[k][(*d-m)][t];
				}
			}			
			for(t=0;t<*T;t++ ) 
			{
				//barvalue[k][i][t]=-der_hatk_hati[t]*der_tildek_tildei[t];
				barvalue_out[(k+1)+(*d)*i+(*d)*(*d)*t-1]=-der_hatk_hati[t]*der_tildek_tildei[t];
				barvdirect[k-1][i][t]=0;
				barvindirect[k-1][i][t]=0;
			}
		
			if(calc2[k+1][i]==1 && calc[k+1][i]==1)
			{
				//Rprintf("Fall 1\n");
				param[0]=theta[k][i];
				param[1]=nu[k][i];
				diff2PDF_u_mod(zr1,zr2,T,param,&fam[k][i],helpvar);
				diffPDF_u_mod(zr1,zr2,T,param,&fam[k][i],helpvar2);

				for(t=0;t<*T;t++)
				{
					LL_mod(&fam[k][i],&a,&zr1[t],&zr2[t],&theta[k][i],&nu[k][i],&helpvar3[t]);
				}
				diffPDF_u_mod(zr1,zr2,T,param,&fam[k][i],helpvar4);
				diffhfunc_v_mod2(zr2,zr1,T,param,&fam[k][i],helpvar5);
				diff2hfunc_v_mod2(zr2,zr1,T,param,&fam[k][i],helpvar6);
				//Rprintf("helpvar %f \n", helpvar[0]);
				//Rprintf("helpvar %f \n", helpvar[1]);
				//Rprintf("helpvar2 %f \n", helpvar2[0]);
				//Rprintf("helpvar2 %f \n", helpvar2[1]);	
				//Rprintf("helpvar3 %f \n", helpvar3[0]);
				//Rprintf("helpvar3 %f \n", helpvar3[1]);
				//Rprintf("helpvar4 %f \n", helpvar4[0]);
				//Rprintf("helpvar4 %f \n", helpvar4[1]);
				//Rprintf("helpvar5 %f \n", helpvar5[0]);
				//Rprintf("helpvar5 %f \n", helpvar5[1]);
				//Rprintf("helpvar6 %f \n", helpvar6[0]);
				//Rprintf("helpvar6 %f \n", helpvar6[1]);
				for(t=0;t<*T;t++ ) 
				{
					//barvalue[k][i][t] = barvalue[k][i][t] + helpvar[t]/cop[t]*tildezr1_hatk_hati[t]*tildezr1_tildek_tildei[t] + helpvar2[t]/cop[t]*barzr1[t];
					barvalue_out[(k+1)+(*d)*i+(*d)*(*d)*t-1] = barvalue_out[(k+1)+(*d)*i+(*d)*(*d)*t-1] + helpvar[t]/cop[t]*tildezr1_hatk_hati[t]*tildezr1_tildek_tildei[t] + helpvar2[t]/cop[t]*barzr1[t];
					barvdirect[k-1][i][t] = barvdirect[k-1][i][t] + exp(helpvar3[t])*barzr1[t] + helpvar4[t]*tildezr1_hatk_hati[t]*tildezr1_tildek_tildei[t];
					barvindirect[k-1][i][t] = barvindirect[k-1][i][t] + helpvar5[t]*barzr1[t] + helpvar6[t]*tildezr1_hatk_hati[t]*tildezr1_tildek_tildei[t];
				}
			}
			
			if(calc2[k+1][*d-m]==1 && calc[k+1][*d-m]==1)
			{
				//Rprintf("Fall2\n");
				//Rprintf("indirect_hat 1 %f \n", tildezr2_hatk_hati[0]);
				//Rprintf("indirect_hat 2 %f \n", tildezr2_hatk_hati[1]);
				//Rprintf("indirect_tilde 1 %f \n", tildezr2_tildek_tildei[0]);
				//Rprintf("indirect_tilde 2 %f \n", tildezr2_tildek_tildei[1]);
				//Rprintf("barzr2 1 %f \n", barzr2[0]);
				//Rprintf("barzr2 2 %f \n", barzr2[1]);
				param[0]=theta[k][i];
				param[1]=nu[k][i];
				diff2PDF_v_mod(zr1,zr2,T,param,&fam[k][i],helpvar);
				diffPDF_v_mod(zr1,zr2,T,param,&fam[k][i],helpvar2);
				diffhfunc_v_mod(zr1,zr2,T,param,&fam[k][i],helpvar3);
				diff2hfunc_v_mod(zr1,zr2,T,param,&fam[k][i],helpvar4);
				for(t=0;t<*T;t++)
				{
					LL_mod2(&fam[k][i],&a,&zr2[t],&zr1[t],&theta[k][i],&nu[k][i],&helpvar5[t]);
				}
				diffPDF_u_mod(zr2,zr1,T,param,&fam[k][i],helpvar6);
				for(t=0;t<*T;t++ ) 
				{
					//barvalue[k][i][t] = barvalue[k][i][t] + helpvar[t]/cop[t]*tildezr2_hatk_hati[t]*tildezr2_tildek_tildei[t] + helpvar2[t]/cop[t]*barzr2[t];
					barvalue_out[(k+1)+(*d)*i+(*d)*(*d)*t-1] = barvalue_out[(k+1)+(*d)*i+(*d)*(*d)*t-1] + helpvar[t]/cop[t]*tildezr2_hatk_hati[t]*tildezr2_tildek_tildei[t] + helpvar2[t]/cop[t]*barzr2[t];
					barvdirect[k-1][i][t] = barvdirect[k-1][i][t] + helpvar3[t]*barzr2[t] + helpvar4[t]*tildezr2_hatk_hati[t]*tildezr2_tildek_tildei[t];
					barvindirect[k-1][i][t] = barvindirect[k-1][i][t] + exp(helpvar5[t])*barzr2[t] + helpvar6[t]*tildezr2_hatk_hati[t]*tildezr2_tildek_tildei[t];
				}
			}
			if(calc2[k+1][i]==1 && calc[k+1][*d-m]==1)
			{
				//Rprintf("Fall3\n");
				//Rprintf("indirect_hat 1 %f \n", tildezr2_hatk_hati[0]);
				//Rprintf("indirect_hat 2 %f \n", tildezr2_hatk_hati[1]);
				//Rprintf("indirect_tilde 1 %f \n", tildezr2_tildek_tildei[0]);
				//Rprintf("indirect_tilde 2 %f \n", tildezr2_tildek_tildei[1]);
				//Rprintf("barzr2 1 %f \n", barzr2[0]);
				//Rprintf("barzr2 2 %f \n", barzr2[1]);
				param[0]=theta[k][i];
				param[1]=nu[k][i];
				diff2PDF_u_v_mod(zr1,zr2,T,param,&fam[k][i],helpvar);
				diffPDF_v_mod(zr1,zr2,T,param,&fam[k][i],helpvar2);
				diffPDF_v_mod(zr2,zr1,T,param,&fam[k][i],helpvar3);				
				for(t=0;t<*T;t++ ) 
				{
					//barvalue[k][i][t] = barvalue[k][i][t] + helpvar[t]/cop[t]*tildezr1_hatk_hati[t]*tildezr2_tildek_tildei[t];
					barvalue_out[(k+1)+(*d)*i+(*d)*(*d)*t-1] = barvalue_out[(k+1)+(*d)*i+(*d)*(*d)*t-1] + helpvar[t]/cop[t]*tildezr1_hatk_hati[t]*tildezr2_tildek_tildei[t];
					barvdirect[k-1][i][t] = barvdirect[k-1][i][t] + helpvar2[t]*tildezr1_hatk_hati[t]*tildezr2_tildek_tildei[t];
					barvindirect[k-1][i][t] = barvindirect[k-1][i][t] + helpvar3[t]*tildezr1_hatk_hati[t]*tildezr2_tildek_tildei[t];
				}
			}
			if(calc2[k+1][*d-m]==1 && calc[k+1][i]==1)
			{
				//Rprintf("Fall4\n");
				//Rprintf("indirect_hat 1 %f \n", tildezr2_hatk_hati[0]);
				//Rprintf("indirect_hat 2 %f \n", tildezr2_hatk_hati[1]);
				//Rprintf("indirect_tilde 1 %f \n", tildezr2_tildek_tildei[0]);
				//Rprintf("indirect_tilde 2 %f \n", tildezr2_tildek_tildei[1]);
				//Rprintf("barzr2 1 %f \n", barzr2[0]);
				//Rprintf("barzr2 2 %f \n", barzr2[1]);
				param[0]=theta[k][i];
				param[1]=nu[k][i];
				diff2PDF_u_v_mod(zr1,zr2,T,param,&fam[k][i],helpvar);
				diffPDF_v_mod(zr1,zr2,T,param,&fam[k][i],helpvar2);
				diffPDF_v_mod(zr2,zr1,T,param,&fam[k][i],helpvar3);				
				for(t=0;t<*T;t++ ) 
				{
					//barvalue[k][i][t] = barvalue[k][i][t] + helpvar[t]/cop[t]*tildezr2_hatk_hati[t]*tildezr1_tildek_tildei[t];
					barvalue_out[(k+1)+(*d)*i+(*d)*(*d)*t-1] = barvalue_out[(k+1)+(*d)*i+(*d)*(*d)*t-1] + helpvar[t]/cop[t]*tildezr2_hatk_hati[t]*tildezr1_tildek_tildei[t];
					barvdirect[k-1][i][t] = barvdirect[k-1][i][t] + helpvar2[t]*tildezr2_hatk_hati[t]*tildezr1_tildek_tildei[t];
					barvindirect[k-1][i][t] = barvindirect[k-1][i][t] + helpvar3[t]*tildezr2_hatk_hati[t]*tildezr1_tildek_tildei[t];
				}
			}
				
		}
	}

	for(k=0;k<(*d);k++)
	{
		for(i=0;i<(*d);i++)
		{
			for(t=0;t<(*T);t++ ) 
			{
				//sumloglik += barvalue[k][i][t];
				sumloglik += barvalue_out[(k+1)+(*d)*i+(*d)*(*d)*t-1];
			}
		}
	}

*out = sumloglik;

for(i=0;i<(*d);i++)
{
	for(j=0;j<(*d);j++){	
		for(t=0;t<*T;t++ ) {
			barvdirect_out[(i+1)+(*d)*j+(*d)*(*d)*t-1]=barvdirect[i][j][t];
			barvindirect_out[(i+1)+(*d)*j+(*d)*(*d)*t-1]=barvindirect[i][j][t];
			//barvalue_out[(i+1)+(*d)*j+(*d)*(*d)*t-1]=barvalue[i][j][t];	
		}}}


//Free memory:
//free_3darray(vdirect,*d,*d); 
//free_3darray(vindirect,*d,*d);
//free_3darray(value,*d,*d);
free_matrix(theta,*d); 
free_matrix(nu,*d); 
free_intmatrix(fam,*d); 
//free_intmatrix(mmat,*d); 
free_intmatrix(calc, *d); 
free_intmatrix(calc2, *d); 
//free_intmatrix(mat, *d);
//free_3darray(tildevindirect,*d,*d); 
//free_3darray(tildevdirect,*d,*d); 
//free_3darray(tildevalue,*d,*d);
//free_3darray(tildevindirect2,*d,*d); 
//free_3darray(tildevdirect2,*d,*d); 
//free_3darray(tildevalue2,*d,*d);
free_3darray(barvindirect,*d,*d); 
free_3darray(barvdirect,*d,*d); 
//free_3darray(barvalue,*d,*d);
Free(zr1); Free(zr2); Free(tildezr1); Free(tildezr2);Free(handle1);Free(cop);Free(der);
Free(helpvar);Free(helpvar2);Free(helpvar3);Free(helpvar4);Free(helpvar5);Free(helpvar6);
Free(tildezr1_hatk_hati); Free(tildezr1_tildek_tildei);Free(barzr1);
Free(tildezr2_hatk_hati); Free(tildezr2_tildek_tildei);Free(barzr2);
Free(der_hatk_hati);Free(der_tildek_tildei);
}



void hesse_step(int* T, int* d, int* family, int* kk, int* ii, int* kkk, int* iii, int* maxmat, int* matrix, int* condirect, int* conindirect, 
				double* par, double* par2, double* data, int* calcupdate, int* calcupdate2, double* out, double* ll, double* vv, double* vv2,
				double* tilde_value, double* tilde_vdirect, double* tilde_vindirect, double* hat_value, double* hat_vdirect, double* hat_vindirect,
				double* barvalue_out, double* barvdirect_out, double* barvindirect_out, int* kk_second, int* kkk_second)
{

	int ii_out, iii_out, kk_out, kkk_out, kk_second_out, kkk_second_out, margin=0;
	int tcop_help=0;
	int tausche=0, t_cop=0;
	int *calcupdate_out, *calcupdate2_out;

	double der1=0, der2=0, *help;

	//double *tilde_vdirect, *tilde_vindirect, *tilde_value, *hat_vdirect, *hat_vindirect, *hat_value;
	//tilde_vdirect = Calloc(((*d)*(*d)*(*T)),double);
	//tilde_vindirect = Calloc(((*d)*(*d)*(*T)),double);
	//tilde_value = Calloc(((*d)*(*d)*(*T)),double);
	//hat_vdirect = Calloc(((*d)*(*d)*(*T)),double);
	//hat_vindirect = Calloc(((*d)*(*d)*(*T)),double);
	//hat_value = Calloc(((*d)*(*d)*(*T)),double);


	ii_out=*ii;
	iii_out=*iii;
	kk_out=*kk;
	kkk_out=*kkk;
	kk_second_out=*kk_second;
	kkk_second_out=*kkk_second;
	calcupdate_out=calcupdate;
	calcupdate2_out=calcupdate2;
	if(*ii > *iii)
	{
		ii_out=*iii;
		iii_out=*ii;
		kk_out=*kkk;
		kkk_out=*kk;
		kk_second_out=*kkk_second;
		kkk_second_out=*kk_second;
		calcupdate_out=calcupdate2;
		calcupdate2_out=calcupdate;
		tausche=1;
	}
	if((*ii==*iii) && (*kk>*kkk))
	{
		kk_out=*kkk;
		kkk_out=*kk;
		kk_second_out=*kkk_second;
		kkk_second_out=*kk_second;
		calcupdate_out=calcupdate2;
		calcupdate2_out=calcupdate;
		tausche=1;
	}



	if((family[(kk_out)+(*d)*(ii_out-1)-1]==2) && kk_second_out==1)
	{
		tcop_help=2;
		VineLogLikRvineDeriv(T, d, family, &kk_out, &ii_out, maxmat, matrix, condirect, conindirect, par, par2, data, 
							  &der1, //Zielvariable
							ll, vv, vv2, calcupdate_out, tilde_vdirect, tilde_vindirect, tilde_value, &tcop_help, &margin);
	}
	else if((family[(kk_out)+(*d)*(ii_out-1)-1]==2) && kk_second_out==0)
	{
		tcop_help=1;
		VineLogLikRvineDeriv(T, d, family, &kk_out, &ii_out, maxmat, matrix, condirect, conindirect, par, par2, data, 
							  &der1, //Zielvariable
							ll, vv, vv2, calcupdate_out, tilde_vdirect, tilde_vindirect, tilde_value, &tcop_help, &margin);
	}
	else
	{
		tcop_help=0;
		VineLogLikRvineDeriv(T, d, family, &kk_out, &ii_out, maxmat, matrix, condirect, conindirect, par, par2, data, 
							  &der1, //Zielvariable
							ll, vv, vv2, calcupdate_out, tilde_vdirect, tilde_vindirect, tilde_value, &tcop_help, &margin);
	}
	if((family[(kkk_out)+(*d)*(iii_out-1)-1]==2) && kkk_second_out==1)
	{
		tcop_help=2;
		VineLogLikRvineDeriv(T, d, family, &kkk_out, &iii_out, maxmat, matrix, condirect, conindirect, par, par2, data, 
							  &der2, //Zielvariable
							ll, vv, vv2, calcupdate2_out, hat_vdirect, hat_vindirect, hat_value, &tcop_help, &margin);
	}
	else if((family[(kkk_out)+(*d)*(iii_out-1)-1]==2) && kkk_second_out==0)
	{
		tcop_help=1;
		VineLogLikRvineDeriv(T, d, family, &kkk_out, &iii_out, maxmat, matrix, condirect, conindirect, par, par2, data, 
							  &der2, //Zielvariable
							ll, vv, vv2, calcupdate2_out, hat_vdirect, hat_vindirect, hat_value, &tcop_help, &margin);
	}
	else
	{
		tcop_help=0;
		VineLogLikRvineDeriv(T, d, family, &kkk_out, &iii_out, maxmat, matrix, condirect, conindirect, par, par2, data, 
							  &der2, //Zielvariable
							ll, vv, vv2, calcupdate2_out, hat_vdirect, hat_vindirect, hat_value, &tcop_help, &margin);
	}

	//subder1=sum(array(tilde_value_array,dim=c(d,d,N))[kk_out,,])

	if(family[(kk_out)+(*d)*(ii_out-1)-1]==2)
	{
		if(kk_second_out==1 && kkk_second_out==1)
		{
			t_cop=2;
		}
		else if((kk_second_out+kkk_second_out)==1)
		{
			t_cop=3;
		}
		else
		{
			t_cop=1;
		}
	}
	else
	{
		t_cop=0;
	}
	
//	Rprintf("kk_out: %d\t ii_out %d \t", kk_out, ii_out);
//	Rprintf("t_cop: %d\n", t_cop);

	// Aufruf der Hauptfunktion
	VineLogLikRvineDeriv2(T, d, family, &kk_out, &ii_out, &kkk_out, &iii_out, maxmat, matrix, condirect, conindirect, par, par2, data, 
						  tilde_value, tilde_vdirect, tilde_vindirect, hat_value, hat_vdirect, hat_vindirect, calcupdate_out, calcupdate2_out,
						  out, ll, vv, vv2, barvalue_out, barvdirect_out, barvindirect_out, &t_cop, &kk_second_out);



	if(tausche==1)
	{
		help=tilde_value;
		tilde_value=hat_value;
		hat_value=help;
	}
	*kk_second=kk_out;
	*kkk_second=kkk_out;


	//Free(tilde_vdirect);Free(tilde_vindirect);Free(tilde_value);
	//Free(hat_vdirect);Free(hat_vindirect);Free(hat_value);
}



void hesse(int* T, int* d, int* family, int* maxmat, int* matrix, int* condirect, int* conindirect, double* par, double* par2, double* data, double* out, double* subhess, double* der, double* subder)
{

	double *loglik, der1=0, der2=0, subder1=0, subder2=0;
	double *ll, *vv, *vv2, *barvalue_out, *barvdirect_out, *barvindirect_out;
	double *tilde_vdirect, *tilde_vindirect, *tilde_value, *hat_vdirect, *hat_vindirect, *hat_value;
	int *calcupdate1, *calcupdate2;
	int i=0, j=0, i1=0, i2=0, k1=0, k2=0, t1=0, t2=0, seperate=1, dd=*d*(*d-1)/2, tt=0, t=0;
	int kk_second=0, kkk_second=0;
	loglik = Calloc(*T,double);
	ll = Calloc(((*d)*(*d)*(*T)),double);
	vv = Calloc(((*d)*(*d)*(*T)),double);
	vv2 = Calloc(((*d)*(*d)*(*T)),double);
	calcupdate1 = Calloc(((*d)*(*d)),int);
	calcupdate2 = Calloc(((*d)*(*d)),int);
	barvalue_out = Calloc(((*d)*(*d)*(*T)),double);
	barvdirect_out = Calloc(((*d)*(*d)*(*T)),double);
	barvindirect_out = Calloc(((*d)*(*d)*(*T)),double);
	tilde_vdirect = Calloc(((*d)*(*d)*(*T)),double);
	tilde_vindirect = Calloc(((*d)*(*d)*(*T)),double);
	tilde_value = Calloc(((*d)*(*d)*(*T)),double);
	hat_vdirect = Calloc(((*d)*(*d)*(*T)),double);
	hat_vindirect = Calloc(((*d)*(*d)*(*T)),double);
	hat_value = Calloc(((*d)*(*d)*(*T)),double);

	for(i=0;i<*d*(*d);i++)
	{
		calcupdate1[i]=1;
	}

	for(i=0;i<*d*(*d)*(*T);i++)
	{
		barvalue_out[i]=0;
		barvdirect_out[i]=0;
		barvindirect_out[i]=0;
	}

	int **fam;

	//Allocate memory
	fam=create_intmatrix(*d,*d);

	for(i=0;i<(*d);i++)
	{
        for(j=0;j<(*d);j++)
		{
			fam[i][j]=family[(i+1)+(*d)*j-1];
			if(family[(i+1)+(*d)*j-1]==2)
			{
				tt++;
			}
		}
	}  

	VineLogLikRvine(T, d, family, maxmat, matrix, condirect, conindirect, par, par2, data, 
		loglik, ll, vv, vv2, calcupdate1, &seperate);


	t1=0;
	for(i1=(*d-1);i1>0;i1--)
	{
		for(k1=*d;k1>i1;k1--)
		{
			t2=0;
			// calcupdate berechen
			calcupdate_func(d, matrix, &k1, &i1, calcupdate1);

			for(i2=(*d-1);i2>0;i2--)
			{
				for(k2=*d;k2>i2;k2--)
				{
					if(t1<=t2)
					{
						if(fam[k1-1][i1-1]==0 || fam[k2-1][i2-1]==0)
						{
							//out[t1][t2]=0;
							out[(t1+1)+((dd+tt)*t2)-1]=0;
							//out[t2][t1]=0;
							out[(t2+1)+((dd+tt)*t1)-1]=0;
						}
						else
						{
							// calcupdate 2 berechen
							calcupdate_func(d, matrix, &k2, &i2, calcupdate2);

							kk_second=0;
							kkk_second=0;
						//	Rprintf("t1: %d\t t2: %d\t",t1,t2);
							hesse_step(T, d, family, &k1, &i1, &k2, &i2, maxmat, matrix, condirect, conindirect, par, par2, data, 
										calcupdate1, calcupdate2, &out[(t1+1)+((dd+tt)*t2)-1], ll, vv, vv2,
										tilde_value, tilde_vdirect, tilde_vindirect, hat_value, hat_vdirect, hat_vindirect,
										barvalue_out, barvdirect_out, barvindirect_out, &kk_second, &kkk_second);
							out[(t2+1)+((dd+tt)*t1)-1]=out[(t1+1)+((dd+tt)*t2)-1];
							
							der1=0;
							der2=0;
							subder1=0;
							subder2=0;
							for(t=0;t<*T;t++)
							{
								der1=0;
								der2=0;
								subder1=0;
								subder2=0;
								for(j=0;j<*d;j++)
								{
									subhess[(t1+1)+((dd+tt)*t2)-1]+=barvalue_out[(k2)+(*d)*j+(*d)*(*d)*t-1];
									subder1+=tilde_value[(k2)+(*d)*j+(*d)*(*d)*t-1];
									subder2+=hat_value[(k1)+(*d)*j+(*d)*(*d)*t-1];
									for(i=0;i<*d;i++)
									{
										der1+=tilde_value[(i+1)+(*d)*j+(*d)*(*d)*t-1];
										der2+=hat_value[(i+1)+(*d)*j+(*d)*(*d)*t-1];
									}
								}
									//Rprintf("%d, %d, %d, %d, der1: %f\t der2: %f\n",i1, k1, i2, k2,der1,der2);
									//Rprintf("t1: %d \t t2: %d \n",t1,t2);
								der[(t1+1)+((dd+tt)*t2)-1]+=der1*der2;
								der[(t2+1)+((dd+tt)*t1)-1]=der[(t1+1)+((dd+tt)*t2)-1];
								subder[(t1+1)+((dd+tt)*t2)-1]+=subder1*subder2;
								subder[(t2+1)+((dd+tt)*t1)-1]=subder[(t1+1)+((dd+tt)*t2)-1];
							}
							subhess[(t2+1)+((dd+tt)*t1)-1]=subhess[(t1+1)+((dd+tt)*t2)-1];
							//der[(t1+1)+((dd+tt)*t2)-1]=der1*der2;
							//der[(t2+1)+((dd+tt)*t1)-1]=der[(t1+1)+((dd+tt)*t2)-1];
							//subder[(t1+1)+((dd+tt)*t2)-1]=subder1*subder2;
							//subder[(t2+1)+((dd+tt)*t1)-1]=subder[(t1+1)+((dd+tt)*t2)-1];
						}
					}
					t2++;
				}
			}
			for(i2=(*d-1);i2>0;i2--)
			{
				for(k2=*d;k2>i2;k2--)
				{
					if(fam[k2-1][i2-1]==2)
					{
						if(t1<=t2)
						{
							if(fam[k1-1][i1-1]==0)
							{
								//out[t1][t2]=0;
								out[(t1+1)+((dd+tt)*t2)-1]=0;
								//out[t2][t1]=0;
								out[(t2+1)+((dd+tt)*t1)-1]=0;
							}
							else
							{
								// calcupdate 2 berechen
								calcupdate_func(d, matrix, &k2, &i2, calcupdate2);

								kk_second=0;
								kkk_second=1;
							//	Rprintf("t1: %d\t t2: %d\t",t1,t2);
								hesse_step(T, d, family, &k1, &i1, &k2, &i2, maxmat, matrix, condirect, conindirect, par, par2, data, 
											calcupdate1, calcupdate2, &out[(t1+1)+((dd+tt)*t2)-1], ll, vv, vv2,
											tilde_value, tilde_vdirect, tilde_vindirect, hat_value, hat_vdirect, hat_vindirect,
											barvalue_out, barvdirect_out, barvindirect_out, &kk_second, &kkk_second);
								out[(t2+1)+((dd+tt)*t1)-1]=out[(t1+1)+((dd+tt)*t2)-1];

								der1=0;
								der2=0;
								subder1=0;
								subder2=0;
								for(t=0;t<*T;t++)
								{
									der1=0;
									der2=0;
									subder1=0;
									subder2=0;
									for(j=0;j<*d;j++)
									{
										subhess[(t1+1)+((dd+tt)*t2)-1]+=barvalue_out[(k2)+(*d)*j+(*d)*(*d)*t-1];
										subder1+=tilde_value[(kk_second)+(*d)*j+(*d)*(*d)*t-1];
										subder2+=hat_value[(kkk_second)+(*d)*j+(*d)*(*d)*t-1];
										for(i=0;i<*d;i++)
										{
											der1+=tilde_value[(i+1)+(*d)*j+(*d)*(*d)*t-1];
											der2+=hat_value[(i+1)+(*d)*j+(*d)*(*d)*t-1];
										}
									}
									der[(t1+1)+((dd+tt)*t2)-1]+=der1*der2;
									der[(t2+1)+((dd+tt)*t1)-1]=der[(t1+1)+((dd+tt)*t2)-1];
									subder[(t1+1)+((dd+tt)*t2)-1]+=subder1*subder2;
									subder[(t2+1)+((dd+tt)*t1)-1]=subder[(t1+1)+((dd+tt)*t2)-1];
								}
								subhess[(t2+1)+((dd+tt)*t1)-1]=subhess[(t1+1)+((dd+tt)*t2)-1];
								//der[(t1+1)+((dd+tt)*t2)-1]=der1*der2;
								//der[(t2+1)+((dd+tt)*t1)-1]=der[(t1+1)+((dd+tt)*t2)-1];
								//subder[(t1+1)+((dd+tt)*t2)-1]=subder1*subder2;
								//subder[(t2+1)+((dd+tt)*t1)-1]=subder[(t1+1)+((dd+tt)*t2)-1];

							}
							
						}
						t2++;
					}
				}
			}
		t1++;
		}
	}
	for(i1=(*d-1);i1>0;i1--)
	{
		for(k1=*d;k1>i1;k1--)
		{
			if(fam[k1-1][i1-1]==2)
			{
				t2=dd;
				// calcupdate berechen
				calcupdate_func(d, matrix, &k1, &i1, calcupdate1);
				/*for(i2=(*d-1);i2>0;i2--)
				{
					for(k2=*d;k2>i2;k2--)
					{
						t2++;
					}
				}*/
				for(i2=(*d-1);i2>0;i2--)
				{
					for(k2=*d;k2>i2;k2--)
					{
						if(fam[k2-1][i2-1]==2)
						{
							if(t1<=t2)
							{
								if(fam[k1-1][i1-1]==0)
								{
									//out[t1][t2]=0;
									out[(t1+1)+((dd+tt)*t2)-1]=0;
									//out[t2][t1]=0;
									out[(t2+1)+((dd+tt)*t1)-1]=0;
								}
								else
								{
									// calcupdate 2 berechen
									calcupdate_func(d, matrix, &k2, &i2, calcupdate2);

									kk_second=1;
									kkk_second=1;
									hesse_step(T, d, family, &k1, &i1, &k2, &i2, maxmat, matrix, condirect, conindirect, par, par2, data, 
											calcupdate1, calcupdate2, &out[(t1+1)+((dd+tt)*t2)-1], ll, vv, vv2,
											tilde_value, tilde_vdirect, tilde_vindirect, hat_value, hat_vdirect, hat_vindirect,
											barvalue_out, barvdirect_out, barvindirect_out, &kk_second, &kkk_second);
									out[(t2+1)+((dd+tt)*t1)-1]=out[(t1+1)+((dd+tt)*t2)-1];

									der1=0;
									der2=0;
									subder1=0;
									subder2=0;
									for(t=0;t<*T;t++)
									{
										der1=0;
										der2=0;
										subder1=0;
										subder2=0;
										for(j=0;j<*d;j++)
										{
											subhess[(t1+1)+((dd+tt)*t2)-1]+=barvalue_out[(k2)+(*d)*j+(*d)*(*d)*t-1];
											subder1+=tilde_value[(kk_second)+(*d)*j+(*d)*(*d)*t-1];
											subder2+=hat_value[(kkk_second)+(*d)*j+(*d)*(*d)*t-1];
											for(i=0;i<*d;i++)
											{
												der1+=tilde_value[(i+1)+(*d)*j+(*d)*(*d)*t-1];
												der2+=hat_value[(i+1)+(*d)*j+(*d)*(*d)*t-1];
											}
										}
										der[(t1+1)+((dd+tt)*t2)-1]+=der1*der2;
										der[(t2+1)+((dd+tt)*t1)-1]=der[(t1+1)+((dd+tt)*t2)-1];
										subder[(t1+1)+((dd+tt)*t2)-1]+=subder1*subder2;
										subder[(t2+1)+((dd+tt)*t1)-1]=subder[(t1+1)+((dd+tt)*t2)-1];
									}
									subhess[(t2+1)+((dd+tt)*t1)-1]=subhess[(t1+1)+((dd+tt)*t2)-1];
									//der[(t1+1)+((dd+tt)*t2)-1]=der1*der2;
									//der[(t2+1)+((dd+tt)*t1)-1]=der[(t1+1)+((dd+tt)*t2)-1];
									//subder[(t1+1)+((dd+tt)*t2)-1]=subder1*subder2;
									//subder[(t2+1)+((dd+tt)*t1)-1]=subder[(t1+1)+((dd+tt)*t2)-1];
								}
							}
							t2++;
						}
					}
				}
			t1++;
			}
		}
	}


	free_intmatrix(fam,*d);
	Free(loglik);
	Free(ll); 
	Free(vv); 
	Free(vv2); 
	Free(calcupdate1); 
	Free(calcupdate2);
	Free(barvalue_out); 
	Free(barvdirect_out); 
	Free(barvindirect_out);
	Free(tilde_vdirect);
	Free(tilde_vindirect);
	Free(tilde_value);
	Free(hat_vdirect);
	Free(hat_vindirect);
	Free(hat_value);
}





////////////////////////////////////

void calcupdate_func(int* d, int* matrix, int* i, int* j, int* calc)
{
	int *g, *c;
	int ii=0, jj=0, kk=0, handle=0, r=0;
	g=Calloc(*d-*i+2,int);
	c=Calloc(*d,int);

	// calc initialisieren
	for(ii=0;ii<(*d)*(*d);ii++)
	{
		calc[ii]=0;
	}


	g[0]=matrix[*j+(*d)*(*j-1)-1];
	for(ii=1;ii<(*d-*i+2);ii++)
	{
		g[ii]=matrix[(*i+ii-1)+(*d)*(*j-1)-1];
	}
	calc[*i+(*d)*(*j-1)-1]=1;

	for(jj=*j;jj>0;jj--)
	{
		for(ii=(*d);ii>(jj);ii--)
		{
			c[0]=matrix[jj+(*d)*(jj-1)-1];
			for(kk=1;kk<(*d-ii+2);kk++)
			{
				c[kk]=matrix[(ii+kk-1)+(*d)*(jj-1)-1];
			}
			handle=0;
			for(r=0;r<(*d-(*i)+2);r++)
			{
				for(kk=0;kk<(*d-ii+2);kk++)
				{
					if(c[kk]==g[r])
					{
						handle++;
					}
				}
			}

			if(handle==(*d-(*i)+2))
			{
				for(kk=jj;kk<ii;kk++)
				{
					calc[(kk+1)+(*d)*(jj-1)-1]=1;
				}
			}
		}
	}

	Free(g);
	Free(c);
}



