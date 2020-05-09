#include <constant.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <complex>
#include "readqebands.h"
#include "indexref.h"
#include "banddot.h"
#include <mpi.h>
#include <sstream>
int main(int argc,char* argv[]){
  MPI_Init(NULL,NULL);
  int world_size,world_rank;
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
  std::fstream fs;
  std::string temp;
  std::string vmatrix;
  std::string bandsfile;
  std::string useless;
  std::stringstream ss;
  int kpointscount,bandnumber;
  if(world_rank==0){
  fs.open(argv[1],std::fstream::in);
  while(getline(fs,temp)){
    if(temp.find("kmesh")!=std::string::npos){
      ss.clear();
      ss.str(temp);
      ss>>useless;
      ss>>sci_const::kx;
      ss>>sci_const::ky;
      ss>>sci_const::kz;
    }
    if(temp.find("datafile")!=std::string::npos){
      ss.clear();
      ss.str(temp);
      ss>>useless;
      ss>>vmatrix;
    }
    if(temp.find("bands")!=std::string::npos){
      ss.clear();
      ss.str(temp);
      ss>>useless;
      ss>>bandsfile;
    }
    if(temp.find("bandnum")!=std::string::npos){
      ss.clear();
      ss.str(temp);
      ss>>useless;
      ss>>sci_const::bandnumber;
    }
  }
  //    readdimension(kpointscount,bandnumber,vmatrix.c_str());
  fs.close();
  }
  else{
  };
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&sci_const::kx,1,MPI::INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&sci_const::ky,1,MPI::INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&sci_const::kz,1,MPI::INT,0,MPI_COMM_WORLD);
  kpointscount=sci_const::kx*sci_const::ky*sci_const::kz;
  bandnumber=sci_const::bandnumber;
  MPI_Bcast(&kpointscount,1,MPI::INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&bandnumber,1,MPI::INT,0,MPI_COMM_WORLD);
  std::complex<double>* kpoint_product=new std::complex<double> [3*kpointscount*bandnumber*(bandnumber+1)/2];
  std::complex<double>*** kpoint_product_array=new std::complex<double>** [3];
  for(size_t i=0;i<3;i++){
    kpoint_product_array[i]=new std::complex<double>* [kpointscount];
    for(size_t j=0;j<kpointscount;j++){
      kpoint_product_array[i][j]=&kpoint_product[(i*kpointscount+j)*bandnumber*(bandnumber+1)/2];
    }
  }
  double* bands=new double [kpointscount*bandnumber];
  double** bands_array=new double* [kpointscount];
  for(size_t i=0;i<kpointscount;i++){
    bands_array[i]=&(bands[i*bandnumber]);
  }
  double* occupation=new double [kpointscount*bandnumber];
  double** occupation_array=new double* [kpointscount];
  for(size_t i=0;i<kpointscount;i++){
    occupation_array[i]=&(occupation[i*bandnumber]);
  }
  double* kpoints=new double [kpointscount*3];
  double** kpoints_array=new double* [kpointscount];
  for(size_t i=0;i<kpointscount;i++){
    kpoints_array[i]=&(kpoints[i*3]);
  }
  double* kweight=new double [kpointscount];
  double volume;
  double kweightsum=0.0;
  if(world_rank==0){
  readkpoints(kweight,volume,bandsfile.c_str());
readalltogether_binary(kpoint_product_array,occupation_array,bands_array,kpointscount,bandnumber,vmatrix);
  for(size_t i=0;i<kpointscount;i++){
    kweight[i]=1.0/kpointscount;
  }
  }
  else{
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&sci_const::alat,1,MPI::DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(bands,kpointscount*bandnumber,MPI::DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(kpoints,kpointscount*3,MPI::DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(kweight,kpointscount,MPI::DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&volume,1,MPI::DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(occupation,kpointscount*bandnumber,MPI::DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(kpoint_product,3*kpointscount*(bandnumber+1)*bandnumber/2,MPI::DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
  double* current_rate;
  double bandgap=searchbandgap(kpointscount,bandnumber,occupation_array,bands_array);
  if(world_rank==0){
    fs.open("spectrum.dat",std::fstream::out);
  }
  for(double photonE=bandgap;photonE<bandgap+1.2;photonE=photonE+0.02){
  current_rate=sumbands(kpointscount,bandnumber,volume,kpoint_product_array,occupation_array,bands_array,kweight,photonE);
  MPI_Barrier(MPI_COMM_WORLD);
  if(world_rank==0){
  fs<<photonE-bandgap<<" "<<current_rate[0]<<" "<<current_rate[1]<<" "<<current_rate[2]<<std::endl;
  }
  }
  if(world_rank==0){
  fs.close();
  }
  delete [] kpoints;
  delete [] occupation;
  delete [] bands;
  delete [] kpoint_product;
 delete [] bands_array;
 delete [] occupation_array;
 for(size_t i=0;i<3;i++){
     delete [] kpoint_product_array[i];
 }
  MPI_Finalize();
}
