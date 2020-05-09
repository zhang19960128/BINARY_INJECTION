#include <math.h>
namespace sci_const{
double e_q=1.60217662*1e-19;
double e_mass=9.10938356*1e-31;
double hbar=1.0545718*1e-34;
double rbohr=5.29177210903*1e-11;
double PI=3.141592653;
double root2=1.41421356237;
double rootpi=1.77245385091;
double ev2j=1.60218*1e-19;
double hztocm=2.9*1e10;
double hatree2ev=27.2114;
double alat;
int kx;
int ky;
int kz;
int bandnumber;
}
namespace light{
double illumination=0.5;/*mw/cm^2*/
double delta=1.5707963265;
/*converting mw/cm^2 to V/m */
/*get the Efield amplitude by E=sqrt(I)*3772.426135064461,refer to the notes in OneNote*/
double mwpercmsq2E=3772.426135064461;
double Eamp=mwpercmsq2E*sqrt(illumination);/*SI units*/
double photonE=3.45;/*ev*/
double freq=photonE*sci_const::ev2j/sci_const::hbar;/*SI units*/
double Ax=Eamp/freq;/*SI units*/
double Ay=0;
double Az=Eamp/freq;/*SI units*/
double time=-1e-14;/*1ps*/
}
namespace gaussian{
double smearing_hertz=0.1*sci_const::ev2j/sci_const::hbar;
double smearing_ev=0.2;/*ev*/
}
