#ifndef READQEBANDS_H
#define READQEBANDS_H
#include <string>
#include <iostream>
#include <complex>
void readbands(double** bands,int kpoints,int bandnumber,std::string nscf);
void readoccupation(double** occupationnumber,int kpoints,int bandnumber,std::string nscf);
void readvmatrix(std::complex<double>*** kpoint_product,int kpoints,int bandnumber,std::string pmat);
void readkpoints(double* kweight,double& volume,std::string outnscf);
void readalltogether(std::complex<double>*** kpoint_product,double** occupationnumber,double** bands,int kpoints,int bandnumber,std::string pmat);
void readalltogether_binary(std::complex<double>*** kpoint_product,double** occupationnumber,double** bands,int kpoints,int bandnumber,std::string& pmat);
void readalltogether_binary_cpp(std::complex<double>*** kpoint_product,double** occupationnumber,double** bands,int kpoints,int bandnumber,std::string& pmat);
void readdimension(int& kpoints,int& bandnumber,std::string pmat);
#endif
