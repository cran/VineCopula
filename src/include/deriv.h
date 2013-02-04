#if !defined(DERIV_H)
#define DERIV_H

// File deriv.c
void diffPDF_mod(double* u, double* v, int* n, double* param, int* copula, double* out);
void diffPDF(double* u, double* v, int* n, double* param, int* copula, double* out);
void diffPDF_u_mod(double* u, double* v, int* n, double* param, int* copula, double* out);
void diffPDF_u(double* u, double* v, int* n, double* param, int* copula, double* out);
void diffPDF_v_mod(double* u, double* v, int* n, double* param, int* copula, double* out);
void diffhfunc_mod(double* u, double* v, int* n, double* param, int* copula, double* out);
void diffhfunc_mod2(double* v, double* u, int* n, double* param, int* copula, double* out);
void diffhfunc(double* u, double* v, int* n, double* param, int* copula, double* out);
void diffhfunc_v_mod(double* u, double* v, int* n, double* param, int* copula, double* out);
void diffhfunc_v_mod2(double* v, double* u, int* n, double* param, int* copula, double* out);
void diffhfunc_v(double* u, double* v, int* n, double* param, int* copula, double* out);

#endif
