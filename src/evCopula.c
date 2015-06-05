#include "include/vine.h"
#include "include/evCopula.h"
#include <math.h>

#define UMAX  1-1e-10

#define UMIN  1e-10

#define XEPS 1e-4

// Some function for the Tawn copula
// (theory based on the extreme value copulas)
// Reference: See help (some master thesis)

// for the calculation of the density as well as for the h-function we need some help functions
// the naming of the functions is due to the notation of the master thesis (and also references therein)

// CDF

///////////////////////////////////
//
// Input:
// t		t-vector
// n		number of observations
// par		first parameter
// par2		second parameter
// par3		third parameter
//
// Output:
// out		ta
//////////////////////////////

void ta(double* t, int* n, double* par, double* par2, double* par3, double* out)	//for CDF
{
	int i=0;
	double t1,t2;
	for(i=0; i<*n;i++)
	{
		t1=pow(*par2*t[i],*par);
		t2=pow(*par3*(1.0-t[i]),*par);
		out[i]=t1+t2;
	}
}

//ta<-function(t,par,par2,par3) {(par2*t)^par+(par3*(1-t))^par}

////////////////////////////////////////////////
// Pickands A for the Tawn copula
// Input:
// t		t-vector
// n		number of observations
// par		first parameter
// par2		second parameter
// par3		third parameter
//
// Output:
// out		Pickands A for the Tawn copula
//////////////////////////////

void Tawn(double* t, int* n, double* par, double* par2, double* par3, double* out)		//für CDF
{
	int i=0, T=1;
	double t1,t2,t3,t4;
	for(i=0; i<*n;i++)
	{
		t1=(1.0-*par3)*(1.0-t[i]);
		t2=(1.0-*par2)*t[i];
		ta(t, &T, par, par2, par3, &t3);
		t4=pow(t3,1.0/(*par));
		out[i]=t1+t2+t4;
	}
}

//Tawn<-function(t,par,par2,par3) {(1-par3)*(1-t)+(1-par2)*t+ta(t,par,par2,par3)^(1/par)}

////////////////////////////////////////////////////
// CDF of Tawn
// Input:
// t		t-vector
// n		number of observations
// par		first parameter
// par2		second parameter
// par3		third parameter
//
// Output:
// out		CDF
/////////////////////////////////////////////////////


void TawnCDF(double* u, double* v, int* n, double* par, double* par2, double* par3, double* out)	// CDF-function
{
	int i=0, T=1;
	double w, A;
	for(i=0; i<*n;i++)
	{
		w=log(v[i])/log(u[i]*v[i]);				// w vector
		Tawn(&w, &T, par, par2, par3, &A);		//Pickands A
		out[i]=pow(u[i]*v[i],A);				// CDF
	}
}


//////////////////////////////////////////////////////////////////
// PDF

// some more help function for the PDF
// see reference for details

void ta2(double* t, int* n, double* par, double* par2, double* par3, double* out)	//für PDF
{
	int i=0;
	double t1,t2;
	for(i=0; i<*n;i++)
	{
		t1=pow(*par3*t[i],*par);
		t2=pow(*par2*(1.0-t[i]),*par);
		out[i]=t1+t2;
	}
}

// something like the first derivative of the ta function

void d1ta(double* t, int* n, double* par, double* par2, double* par3, double* out)	//für PDF
{
	int i=0;
	double t1,t2;
	for(i=0; i<*n;i++)
	{
		t1=*par3 * pow((*par3*t[i]),*par-1.0);
		t2=*par2 * pow(*par2*(1.0-t[i]),*par-1.0);
		out[i]=*par*(t1-t2);
	}
}

//d1ta<-function(t,par,par2,par3) {par*(par3*(par3*t)^(par-1)-par2*(par2*(1-t))^(par-1))}

void d2ta(double* t, int* n, double* par, double* par2, double* par3, double* out)	//für PDF
{
	int i=0;
	double t1,t2;
	for(i=0; i<*n;i++)
	{
		t1=pow(*par3,2) * pow(*par3*t[i],*par-2.0);
		t2=pow(*par2,2) * pow(*par2*(1.0-t[i]),*par-2.0);
		out[i]=*par*(*par-1) * (t1 + t2);
	}
}

//d2ta<-function(t,par,par2,par3) {par*(par-1)*(par3^2*(par3*t)^(par-2)+par2^2*(par2*(1-t))^(par-2))}

// I guess this was some kind of derivative of A (I don't remember, see master thesis)
void Tawn2(double* t, int* n, double* par, double* par2, double* par3, double* out)		//für PDF
{
	int i=0, T=1;
	double t1,t2,t3,t4;
	for(i=0; i<*n;i++)
	{
		t1=(1.0-*par2)*(1.0-t[i]);
		t2=(1.0-*par3)*t[i];
		ta2(&t[i], &T, par, par2, par3, &t3);		//!!!
		t4=pow(t3,1.0/(*par));
		out[i]=t1+t2+t4;
	}
}

//Tawn<-function(t,par,par2,par3) {(1-par2)*(1-t)+(1-par3)*t+ta(t,par,par2,par3)^(1/par)}

void d1Tawn(double* t, int* n, double* par, double* par2, double* par3, double* out)	//für PDF
{
	int i=0, T=1;
	double t2,t1;
	for(i=0; i<*n;i++)
	{
		ta2(&t[i], &T, par, par2, par3, &t1);		//!!
		d1ta(&t[i], &T, par, par2, par3, &t2);
		out[i]=*par2-*par3+1.0/(*par) * pow(t1,(1.0/(*par)-1.0)) * t2;
	}
}

//d1Tawn<-function(t,par,par2,par3) {par2-par3+1/par*ta(t,par,par2,par3)^(1/par-1)*d1ta(t,par,par2,par3)}  Wie in Afunc2Deriv

void d2Tawn(double* t, int* n, double* par, double* par2, double* par3, double* out)	//für PDF
{
	int i=0, T=1;
	double t2,t1,t3;
	for(i=0; i<*n;i++)
	{
		ta2(&t[i], &T, par, par2, par3, &t1);	//!!
		d1ta(&t[i], &T, par, par2, par3, &t2);
		d2ta(&t[i], &T, par, par2, par3, &t3);
		out[i] = 1.0/(*par) * ( (1.0/(*par)-1.0) * pow(t1,(1.0/(*par)-2)) * pow(t2,2) + pow(t1,(1.0/(*par)-1)) * t3   );
	}
}

//d2Tawn<-function(t,par,par2,par3) {1/par*((1/par-1)*ta(t,par,par2,par3)^(1/par-2)*d1ta(t,par,par2,par3)^2+ta(t,par,par2,par3)^(1/par-1)*d2ta(t,par,par2,par3))}

// Ableitung von A nach u
// derivative of A with respect to u (first argument)
// needed for the derivative of c with respect to u
void dA_du(double* u, double* v, int* n, double* par, double* par2, double* par3, double* out)
{
	int i=0, T=1;
	double dA, dw, w;
	for(i=0; i<*n;i++)
	{
		w=log(v[i]) / log(u[i]*v[i]);
		dw=-log(v[i]) / (u[i]*pow(log(u[i]*v[i]),2.0));
		d1Tawn(&w, &T, par, par2, par3, &dA);
		out[i]=dA*dw;
	}
}

//dA_du<-function(u,v,par,par2,par3) {evcBiCopAfuncDeriv(w(u,v),fam,par,par2,par3)*dw_du(u,v)} 

// derivative of A with respect to v

void dA_dv(double* u, double* v, int* n, double* par, double* par2, double* par3, double* out)
{
	int i=0, T=1;
	double dA, dw, w;
	for(i=0; i<*n;i++)
	{
		w=log(v[i])/log(u[i]*v[i]);
		dw=1.0 / (v[i]*log(u[i]*v[i])) - log(v[i]) / (v[i]*pow(log(u[i]*v[i]),2));
		d1Tawn(&w, &T, par, par2, par3, &dA);
		out[i]=dA*dw;
	}
}

// second derivative with respect to u and v

void dA_dudv(double* u, double* v, int* n, double* par, double* par2, double* par3, double* out)
{
	int i=0, T=1;
	double dA, dw_dv, dw_du, w, d2w_dudv, d2A;
	for(i=0; i<*n;i++)
	{
		w=log(v[i])/log(u[i]*v[i]);
		dw_du=-log(v[i])/(u[i]*pow(log(u[i]*v[i]),2));
		dw_dv=1.0 / (v[i]*log(u[i]*v[i])) - log(v[i]) / (v[i]*pow(log(u[i]*v[i]),2));
		d2w_dudv = 2*log(v[i]) / (v[i]*u[i]*pow(log(u[i]*v[i]),3)) - 1.0 / (v[i]*u[i]*pow(log(u[i]*v[i]),2));
		d1Tawn(&w, &T, par, par2, par3, &dA);
		d2Tawn(&w, &T, par, par2, par3, &d2A);
		out[i]=d2A*dw_dv*dw_du + dA*d2w_dudv;
	}
}

//d2A_dudv<-function(u,v,par,par2,par3) {evcBiCopAfunc2Deriv(w(u,v),fam,par,par2,par3)*dw_dv(u,v)*dw_du(u,v)+evcBiCopAfuncDeriv(w(u,v),fam,par,par2,par3)*d2w_dudv(u,v)}


void TawnC(double* u, double* v, int* n, double* par, double* par2, double* par3, double* out)	// für PDF
{
	int i=0, T=1;
	double w, A;
	for(i=0; i<*n;i++)
	{
		w=log(v[i])/log(u[i]*v[i]);
		Tawn2(&w, &T, par, par2, par3, &A);		//!!!
		out[i]=pow(u[i]*v[i],A);
	}
}

// Ableitung von C nach u
// derivative of PDF with respect to u
void dC_du(double* u, double* v, int* n, double* par, double* par2, double* par3, double* out)
{
	int i=0, T=1;
	double w, A, C, dA;
	for(i=0; i<*n;i++)
	{
		w=log(v[i])/log(u[i]*v[i]);
		Tawn2(&w, &T, par, par2, par3, &A);		//!!
		TawnC(&u[i], &v[i], &T, par, par2, par3, &C);	//!!
		dA_du(&u[i], &v[i], &T, par, par2, par3, &dA);
		out[i]=C * (1.0 / u[i] * A + log(u[i]*v[i])*dA);
	}
}

//dC_du<-function(u,v,par,par2,par3) {C(u,v,par,par2,par3) * (1/u*evcBiCopAfunc(w(u,v),fam,par,par2,par3)+log(u*v) * dA_du(u,v,par,par2,par3))}

void TawnPDF(double* u, double* v, int* n, double* par, double* par2, double* par3, double* out)	
{
	int i=0, T=1;
	double w, A, dC, t3, t4, t1, C, t5, t2;
	for(i=0; i<*n;i++)
	{
		w=log(v[i])/log(u[i]*v[i]);
		Tawn2(&w, n, par, par2, par3, &A);
		dC_du(&u[i], &v[i], &T, par, par2, par3, &dC);
		dA_du(&u[i], &v[i], &T, par, par2, par3, &t3);
		dA_dv(&u[i], &v[i], &T, par, par2, par3, &t4);
		t1=dC * (1.0/v[i] * A + log(u[i]*v[i]) * t4);
		TawnC(&u[i], &v[i], &T, par, par2, par3, &C);
		dA_dudv(&u[i], &v[i], &T, par, par2, par3, &t5);
		t2=C * (1.0/u[i]*t4 + 1.0/v[i]*t3 + log(u[i]*v[i])*t5);
		out[i]=t1+t2;
	}
}


// d2C_dvdu<-function(u,v,par,par2,par3) 
// {
// dC_du(u,v,par,par2,par3)*(1/v*evcBiCopAfunc(w(u,v),fam,par,par2,par3) + log(u*v)*dA_dv(u,v,par,par2,par3))+
// C(u,v,par,par2,par3)*
// (1/u*dA_dv(u,v,par,par2,par3) + 1/v*dA_du(u,v,par,par2,par3)+log(u*v)*d2A_dudv(u,v,par,par2,par3))
// }


// Ableitung von C nach v (fuer h-function)
// derivative of PDF with respect to v (for h-func)
void dC_dv(double* u, double* v, int* n, double* par, double* par2, double* par3, double* out)
{
	int i=0, T=1;
	double w, A, C, dA;
	for(i=0; i<*n;i++)
	{
		w=log(v[i])/log(u[i]*v[i]);
		Tawn2(&w, &T, par, par2, par3, &A);		//!!
		TawnC(&u[i], &v[i], &T, par, par2, par3, &C);	//!!
		dA_dv(&u[i], &v[i], &T, par, par2, par3, &dA);
		out[i]=C * (1.0 / v[i] * A + log(u[i]*v[i])*dA);
	}
}
