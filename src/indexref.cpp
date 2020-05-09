#include "indexref.h"
#include <complex>
#include <stdlib.h>
#include <iostream>
//<m|P|n>=<m|P^{+}|n>=<n|P|m>*
int findindex(int m,int n,int band){
  //enforce m<n
  if(m>n){
    std::cout<<"Call Inner Product function Wrong, should always enforce m<n"<<std::endl;
    std::cout<<" if you want to see the m>=n situation, please use complex conjugate"<<std::endl;
    exit(EXIT_FAILURE);
  }
  else{
    return (2*band-m+1)*m/2+(n-m);
  }
}
std::complex<double> indexvmatrix(int kpoint,int m,int n,int totalbands,int direction,std::complex<double>*** kpoints_product){
 int indextemp;
 if(m<=n){
 indextemp=findindex(m,n,totalbands);
 return kpoints_product[direction][kpoint][indextemp];
 }
 else{
 indextemp=findindex(n,m,totalbands);
 /*return the complex conjugate*/
 return conj(kpoints_product[direction][kpoint][indextemp]);
 }
}
