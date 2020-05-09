#ifndef BANDDOT_H
#define BANDDOT_H
#include <complex>
double* bandot(int kindex1,int kindex2,int bandnum1,int bandnum2,double volume,int kpoint_total,int bandtotal,std::complex<double>*** kpoint_product,double** occupation,double** bands,double* kweight,double freq);
double smearing(double input,double center,double smearing);
double* sumbands(int kpointstotal,int bandstotal,double volume,std::complex<double>*** kpointsproduct,double** occupation,double** bands,double* weight,double freq);
double searchbandgap(int kpointstotal,int bandstotal,double** occupation,double** bands);
#endif
