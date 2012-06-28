// File hfunc.c
void Hfunc1(int* family,int* n,double* u,double* v,double* theta,double* nu,double* out);
void Hfunc2(int* family,int* n,double* v,double* u,double* theta,double* nu,double* out);
void Hfunc(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out);
void Hinv1(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out);
void Hinv2(int* family, int* n, double* v, double* u, double* theta, double* nu, double* out);
void Hinv(int* family, int* n, double* u, double* v, double* theta, double* nu, double* out);
void HNumInv(int* family, double* u, double* v, double* theta, double* nu, double* out);
void pcondbb1(double* u, double* v, int* n, double* param, double* out);
void pcondbb6(double* u, double* v, int* n, double* param, double* out);
void pcondbb7(double* u, double* v, int* n, double* param, double* out);
void pcondbb8(double* u, double* v, int* n, double* param, double* out);
