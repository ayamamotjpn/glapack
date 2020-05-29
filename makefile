include ./make.inc

ARC = 86
OBJD = obj_$(ARC)
BIND = bin_$(ARC)
CUD = $(shell pwd)

#swp_band2d_g : $(PRG)

FSRC = glapack.f90
#LIB = -L/opt/rh/devtoolset-7/root/lib/gcc/x86_64-redhat-linux/7

# for using CUDAFOR fortran_thunking.c is needed
# but this is not supplied in cuda-10.1
# pgi compiler 19.4 may support a fortran
# wrapper for cubda-10.1
# cuda-10.0 incliude fortran_thunking.c

# Fortran wrapper fortran.c for using cuda
#CSRC = $(CUDADIR)/src/fortran_thunking.c 

ifeq ($(VER),GPU)
	OBJ = $(OBJD)/cublasXt.o $(OBJD)/glapack.o $(OBJD)/fortran.o
else
	OBJ = $(OBJD)/glapack.o
endif


MAGMABLAS = $(MAGMAROOT)/control/magmablas_   # for magmablasf interface
INCC = $(CUDAINC) -I$(CUDADIR)/src -I./       #for using fortran.c fortran interface of cublas

#all : libglapack.so libgblas.so

ifeq ($(VER),GPU)
  all : libglapack_gpu.so libgblas_gpu.so

  libglapack_gpu.so : $(OBJ) 
  #	$(F90) -pg $(OBJ) $(JDIR) $(LIB) -o $(PRG)     # for gpu version without cublas
	$(F90) $(SHARED) $(OBJ) $(LIB) -o libglapack_gpu.so     # for gpu version with magma and cublas

  libgblas_gpu.so : $(OBJD)/fortran.o
	$(CC) $(SHARED) $(OBJD)/fortran.o $(LIB) -o libgblas_gpu.so
	#$(PRG) : $(OBJ)
	#$(F90) -pg $(OBJ) $(JDIR) $(LIB) -o $(PRG)     # for gpu version
else
  all : libglapack_cpu.so libgblas_cpu.so

  libglapack_cpu.so : $(OBJ)	
	$(F90) $(SHARED) $(OBJ) $(LIB) -o libglapack_cpu.so            # for cpu version

  libgblas_cpu.so : $(OBJD)/fortran.o
	$(CC) $(SHARED) $(OBJD)/fortran.o $(LIB) -o libgblas_cpu.so
endif

#libgblas.so : $(OBJD)/fortran.o
#	$(CC) $(SHARED) $(OBJD)/fortran.o $(LIB) -o libgblas.so

#ifeq ($(VER),GPU)
#libfblas.so :   $(OBJD)/fortran.o    # blas fortran wrapper library 
  #	$(F90) -pg $(OBJ) $(JDIR) $(LIB) -o $(PRG)     # for gpu version without cublas
#	$(CC) $(SHARED) $(OBJD)/fortran.o -o libfblas.so     # for gpu version with magma and cublas
#endif

$(OBJD)/cublasXt.o : cublasXt.f90
	$(F90) $(FCOPT) $(INC) cublasXt.f90 $(JDIR) -o $(OBJD)/cublasXt.o
	
$(OBJD)/glapack.o : glapack.f90
	$(F90) $(FCOPT) $(INC) glapack.f90 $(JDIR) -o $(OBJD)/glapack.o

$(OBJD)/fortran.o : $(CUDADIR)/src/fortran.c
	$(CC) $(CCOPT) $(CUDADIR)/src/fortran.c $(CUDAINC) -DCUBLAS_GFORTRAN -o $(OBJD)/fortran.o

#ifneq ($(F90),pgf90)		
#$(OBJD)/fortran.o : $(CSRC)       # for getting fortran interface for cublas
#	$(CC) $(CCOPT) $(CSRC) $(INCC) -DCUBLAS_GFORTRAN -o $(OBJD)/fortran.o 
#endif

#install_intel :
#	cp libglapack.so /usr/local/lib/glapack/intel
#	cp glapack.mod /usr/local/include/glapack/intel

clean :
	rm -f libglapack_*.so libgblas_*.so $(OBJD)/*.mod *.mod $(OBJ)