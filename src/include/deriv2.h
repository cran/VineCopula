#if !defined(DERIV2_H)
#define DERIV2_H

// Files für deriv2.c

void diff2PDF_mod(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_u_mod(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_v_mod(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_u(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_u_v_mod(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_u_v(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_par_u_mod(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_par_u(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2PDF_par_v_mod(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_mod(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_mod2(double* v, double* u, int* n, double* param, int* copula, double* out);
void diff2hfunc(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_v_mod(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_v_mod2(double* v, double* u, int* n, double* param, int* copula, double* out);
void diff2hfunc_v(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_par_v_mod(double* u, double* v, int* n, double* param, int* copula, double* out);
void diff2hfunc_par_v_mod2(double* v, double* u, int* n, double* param, int* copula, double* out);
void diff2hfunc_par_v(double* u, double* v, int* n, double* param, int* copula, double* out);

#endif
