#include "readqebands.h"
#include "constant.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "indexref.h"
void readbands(double** bands,int kpoints,int bandnumber,std::string nscf){
  std::fstream fs;
  fs.open(nscf.c_str(),std::fstream::in);
  int count=0;
  std::string temp;
  while(getline(fs,temp)){
    if(temp.find("bands (ev):")!=std::string::npos){
      getline(fs,temp);
      for(size_t i=0;i<bandnumber;i++){
        fs>>bands[count][i];
        bands[count][i]=bands[count][i];
      }
      count=count+1;
    }
  }
  fs.close();
}
void readoccupation(double** occupationnumber,int kpoints,int bandnumber,std::string nscf){
  std::fstream fs;
  fs.open(nscf.c_str(),std::fstream::in);
  int count=0;
  std::string temp;
  while(getline(fs,temp)){
    if(temp.find("occupation numbers")!=std::string::npos){
      for(size_t i=0;i<bandnumber;i++){
        fs>>occupationnumber[count][i];
      }
      count=count+1;
    }
  }
  fs.close();
}
void readvmatrix(std::complex<double>*** kpoint_product,int kpoints,int bandnumber,std::string pmat){
 std::fstream fs;
 std::stringstream ss;
 fs.open(pmat.c_str(),std::fstream::in);
 std::string temp;
 int m,n;
 double temp_double;
 while(getline(fs,temp)){
   ss.clear();
   ss.str(temp);
   ss>>kpoints;
   ss>>m;
   ss>>n;
   if(m<n){
   for(size_t i=0;i<3;i++){
    ss>>temp_double;
    kpoint_product[i][kpoints-1][findindex(m-1,n-1,bandnumber)].real(temp_double);
   }
   for(size_t i=0;i<3;i++){
    ss>>temp_double;
    kpoint_product[i][kpoints-1][findindex(m-1,n-1,bandnumber)].imag(temp_double);
   }
   }
 }
fs.close();
}
void readalltogether(std::complex<double>*** kpoint_product,double** occupationnumber,double** bands,int kpoints,int bandnumber,std::string pmat){
 std::fstream fs;
 std::stringstream ss;
 fs.open(pmat.c_str(),std::fstream::in);
 std::string temp;
 int m,n;
 double temp_double;
 while(getline(fs,temp)){
   ss.clear();
   ss.str(temp);
   ss>>kpoints;
   ss>>m;
   ss>>n;
   if(m<=n){
   for(size_t i=0;i<3;i++){
    ss>>temp_double;
    kpoint_product[i][kpoints-1][findindex(m-1,n-1,bandnumber)].real(temp_double);
   }
   for(size_t i=0;i<3;i++){
    ss>>temp_double;
    kpoint_product[i][kpoints-1][findindex(m-1,n-1,bandnumber)].imag(temp_double);
   }
   /*reading the bands*/
   ss>>temp_double;
   bands[kpoints-1][m-1]=temp_double*sci_const::hatree2ev;
   ss>>temp_double;
   bands[kpoints-1][n-1]=temp_double*sci_const::hatree2ev;
   /*reading the occupations*/
   ss>>temp_double;
   occupationnumber[kpoints-1][m-1]=temp_double/2.0;
   ss>>temp_double;
   occupationnumber[kpoints-1][n-1]=temp_double/2.0;
   }
 }
fs.close();
}

void readalltogether_binary(std::complex<double>*** kpoint_product,double** occupationnumber,double** bands,int kpoints,int bandnumber,std::string& pmat){
  double local[10];
  FILE* fp;
  fp=fopen(pmat.c_str(),"rb");
  for(int i=0;i<kpoints;i++){
  	for(int m=0;m<bandnumber;m++){
		for(int n=0;n<bandnumber;n++){
			fread(local,sizeof(double),10,fp);
			if(m<=n){
			for(int k=0;k<3;k++){
			kpoint_product[k][i][findindex(m,n,bandnumber)].real(local[k]);
			}
			for(int k=0;k<3;k++){
			kpoint_product[k][i][findindex(m,n,bandnumber)].imag(local[k+3]);
			}
			bands[i][m]=local[6]*sci_const::hatree2ev;
			bands[i][n]=local[7]*sci_const::hatree2ev;
			occupationnumber[i][m]=local[8]/2.0;
			occupationnumber[i][n]=local[9]/2.0;
			}
		}
}
  }
  fclose(fp);
}
void readalltogether_binary_cpp(std::complex<double>*** kpoint_product,double** occupationnumber,double** bands,int kpoints,int bandnumber,std::string& pmat){
  double local[10];
  std::fstream fs;
  fs.open(pmat.c_str(),std::ios::binary | std::ios::in);
  for(size_t i=0;i<kpoints;i++){
  	for(size_t m=0;m<bandnumber;m++){
		for(size_t n=0;n<bandnumber;n++){
                        for(size_t j=0;j<10;j++){
                        fs.read((char* )&local[j],sizeof(double));
                        }
			if(m<=n){
			for(size_t k=0;k<3;k++){
			kpoint_product[k][i][findindex(m,n,bandnumber)].real(local[k]);
			}
			for(size_t k=0;k<3;k++){
			kpoint_product[k][i][findindex(m,n,bandnumber)].imag(local[k+3]);
			}
			bands[i][m]=local[6]*sci_const::hatree2ev;
			bands[i][n]=local[7]*sci_const::hatree2ev;
			occupationnumber[i][m]=local[8]/2.0;
			occupationnumber[i][n]=local[9]/2.0;
			}
		}
}
  }
  fs.close();
}
void readkpoints(double* kweight,double& volume,std::string outnscf){
  std::fstream fs;
  fs.open(outnscf.c_str(),std::fstream::in);
  std::string temp;
  std::stringstream ss;
  std::string useless;
  double kx,ky,kz;
  double alat;
  while(getline(fs,temp)){
    if(temp.find("lattice parameter (alat)")!=std::string::npos){
    ss.clear();
      ss.str(temp);
      for(size_t i=0;i<4;i++){
      ss>>useless;
      }
      ss>>alat;
      sci_const::alat=alat*sci_const::rbohr;
    }
    if(temp.find("unit-cell volume")!=std::string::npos){
    ss.clear();
    ss.str(temp);
    ss>>useless;
    ss>>useless;
    ss>>useless;
    ss>>volume;
    volume=volume*sci_const::rbohr*sci_const::rbohr*sci_const::rbohr/sci_const::alat/sci_const::alat/sci_const::alat;
    ss.clear();
    }
  }
  fs.close();
}
void readdimension(int& kpoints,int& bandnumber,std::string pmat){
  std::fstream fs;
  fs.open(pmat.c_str(),std::fstream::in);
  std::string temp;
  std::stringstream ss;
  while(getline(fs,temp)){
    ss.clear();
    ss.str(temp);
    ss>>kpoints;
    ss>>bandnumber;
    ss>>bandnumber;
  }
  fs.close();
}
