#include "banddot.h"
#include <iostream>
#include "constant.h"
#include "indexref.h"
#include <cmath>
#include <complex>
#include <mpi.h>
double* bandot(int kindex1,int kindex2,int bandnum1,int bandnum2,double volume,int kpoint_total,int bandtotal,std::complex<double>*** kpoint_product,double** occupation,double** bands,double* kweight,double freq){
  /*Please refer to my OneNote math constant.*/
   int world_size,world_rank;
   MPI_Comm_size(MPI_COMM_WORLD,&world_size);
   MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  double sci=(sci_const::e_q/sci_const::e_mass);
  sci=sci*sci*sci/sci_const::hbar/sci_const::hbar;
  /*deal with the delta function*/
  sci=sci*sci_const::hbar/sci_const::ev2j;
  sci=sci*(-2*sci_const::PI);
  /*velocity matrix units*/
  sci=sci*(sci_const::hbar)*(sci_const::hbar)*(sci_const::hbar)*(2*sci_const::PI/sci_const::alat)*(2*sci_const::PI/sci_const::alat)*(2*sci_const::PI/sci_const::alat);
  /*Integration over dk'*/
  sci=sci;//kweight[kindex1];
  /*Integration over dk''*/
  sci=sci*kweight[kindex2];//kweight[kindex2];
  /*devided by the volume to get the density*/
  sci=sci/(volume*sci_const::alat*sci_const::alat*sci_const::alat);
  /*Changing the light A into E*/
  double prod=sci/(sci_const::ev2j/sci_const::hbar)/(sci_const::ev2j/sci_const::hbar)/4.0;
  /*time the band weight*/
  prod=prod*light::time;
  prod=prod*(occupation[kindex1][bandnum1]-occupation[kindex2][bandnum2]);
  double omega1=bands[kindex1][bandnum1];
  double omega2=bands[kindex2][bandnum2];
  std::complex<double>* P11=new std::complex<double> [3];
  double* result=new double[3];
  for(size_t i=0;i<3;i++){
    result[i]=0.0;
  }
  if(kindex1!=kindex2){
   return result;
  }
  int tempindex;
  for(size_t i=0;i<3;i++){
     P11[i]=indexvmatrix(kindex1,bandnum1,bandnum1,bandtotal,i,kpoint_product);
  
  }
  std::complex<double>* P12=new std::complex<double> [3];
  for(size_t i=0;i<3;i++){
    P12[i]=indexvmatrix(kindex1,bandnum1,bandnum2,bandtotal,i,kpoint_product);
  }
  std::complex<double>* P21=new std::complex<double> [3];
  for(size_t i=0;i<3;i++){
    P21[i]=indexvmatrix(kindex1,bandnum2,bandnum1,bandtotal,i,kpoint_product);
  }
  /*u=+1*/
  double* result1=new double[3];
  std::complex<double> tempcomplex;
  double resultplus=prod*smearing(freq,-(omega2-omega1),gaussian::smearing_ev);
  for(size_t i=0;i<3;i++){
    tempcomplex=(1.0)*sin(light::delta)*(P12[0]*P21[2]-P21[0]*P12[2]);
    result1[i]=(resultplus*tempcomplex*P11[i]).imag();
  }

  /*u=-1*/
  double* result2=new double[3];
  resultplus=prod*smearing(freq,(omega2-omega1),gaussian::smearing_ev);
  for(size_t i=0;i<3;i++){
    tempcomplex=(-1.0)*sin(light::delta)*(P12[0]*P21[2]-P21[0]*P12[2]);
    result2[i]=(resultplus*tempcomplex*P11[i]).imag();
  }
  for(size_t i=0;i<3;i++){
    result[i]=-1*result1[i]-1*result2[i];
  }

  delete [] result1;
  delete [] result2;
  delete [] P12;
  delete [] P21;
  delete [] P11;
  return result;
}
double* sumbands(int kpointstotal,int bandstotal,double volume,std::complex<double>*** kpointsproduct,double** occupation,double** bands,double* kweight,double freq){
   double* totalsum=new double[3];
   double* reducesum=new double[3];
   for(size_t i=0;i<3;i++){
    totalsum[i]=0.0;
    reducesum[i]=0.0;
   }
   double* tempsum;
   int world_rank,world_size;
   MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
   MPI_Comm_size(MPI_COMM_WORLD,&world_size);
   /*sum over kpoints*/
   for(int i=world_rank;i<kpointstotal;i=i+world_size){
    /*sum over n1*/
    for(size_t j=0;j<bandstotal;j++){
    /*sum over n2*/
      for(size_t k=0;k<bandstotal;k++){
          tempsum=bandot(i,i,j,k,volume,kpointstotal,bandstotal,kpointsproduct,occupation,bands,kweight,freq);
          for(size_t m=0;m<3;m++){
            totalsum[m]=totalsum[m]+tempsum[m];
          }
          delete [] tempsum;
      }
    }
   }
   MPI_Reduce(totalsum,reducesum,3,MPI::DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   return reducesum;
}
double searchbandgap(int kpointstotal,int bandstotal,double** occupation,double** bands){
  int world_size,world_rank;
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  double gap,minigap,minigap_root;
  gap=0.0;
  minigap=1000.0;
  for(int i=world_rank;i<kpointstotal;i=i+world_size){
    for(size_t j=0;j<bandstotal-1;j++){
      if((occupation[i][j]-occupation[i][j+1])>0.5){
        gap=bands[i][j+1]-bands[i][j];
      }
    }
    if(gap<minigap){
      minigap=gap;
    }
  }
  MPI_Allreduce(&minigap,&minigap_root,1,MPI::DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  return minigap_root;
}
/*here the smearing represent the sigma, f(x)=1/(sqrt(2*PI)*sigma)*exp(-1/2*((x-center)/sigma)^2)*/
double smearing(double input,double center,double smearing){
  double result=0.0;
  result=-1.0/2.0*((input-center)/smearing)*((input-center)/smearing);
  return 1.0/smearing/sci_const::root2/sci_const::rootpi*exp(result);
}
