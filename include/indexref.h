#ifndef INDEXREF_H
#define INDEXREF_H
#include <complex>
int findindex(int m,int n,int band);
std::complex<double> indexvmatrix(int kpoint,int m,int n,int totalbands,int direction,std::complex<double>*** kpoints_product);
#endif
