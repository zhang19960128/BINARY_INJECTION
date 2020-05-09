
MPIF90 = ftn
LD = ftn
RM = rm -f

IOTK = /p/home/jiahaoz/q-e-qe-6.1.0/S3DE/iotk
FFT = /opt/cray/pe/fftw/3.3.4.11/broadwell/include

IFLAGS = -I$(IOTK)/src -I$(FFT)
LIBS = -L$(IOTK)/src -L/opt/cray/pe/fftw/3.3.4.11/broadwell/lib -liotk
FFLAGS = -O3 $(IFLAGS) $(LIBS) -fpp -D__MPI #-opt-report=3 -opt-report-phase=vec -xCORE-AVX2 
