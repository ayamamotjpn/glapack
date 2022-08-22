include make.inc

ARC = 86
OBJD = obj_$(ARC)
BIND = bin_$(ARC)
CUD = $(shell pwd)

#swp_band2d_g : $(PRG)

FSRC = glapack.f90

# for using CUDAFOR fortran_thunking.c is needed
# but this is not supplied in cuda-10.1
# pgi compiler 19.4 may support a fortran
# wrapper for cubda-10.1
# cuda-10.0 incliude fortran_thunking.c

# Fortran wrapper fortran.c for using cuda
#CSRC = $(CUDADIR)/src/fortran_thunking.c 

OBJ = $(OBJD)/glapack.o
ifeq ($(VER),GPU)
	OBJ += $(OBJD)/cublasXt.o $(OBJD)/fortran.o
	FSRC += cublasXt.f90
endif

MAGMABLAS = $(MAGMAROOT)/control/magmablas_   # for magmablasf interface
INCC = $(CUDAINC) $(CUDADIRINC) -I./          #for using fortran.c fortran interface of cublas

#all : libglapack.so libgblas.so

ifeq ($(VER),GPU)
  all : libglapack_gpu.so libgblas_gpu.so

  libglapack_gpu.so : $(OBJ) 
  #	$(F90) -pg $(OBJ) $(JDIR) $(LIB) -o $(PRG)            # for gpu version without cublas
	$(F90) $(SHARED) $(OBJ) $(LIB) -o libglapack_gpu.so     # for gpu version with magma and cublas

  libgblas_gpu.so : $(OBJD)/fortran.o
	$(CC) $(SHARED) $(OBJD)/fortran.o $(LIB) -o libgblas_gpu.so
	#$(PRG) : $(OBJ)
	#$(F90) -pg $(OBJ) $(JDIR) $(LIB) -o $(PRG)     # for gpu version

else

  all : libglapack_cpu.so          #libgblas_cpu.so

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
	$(F90) $(FCOPT) $(INC) cublasXt.f90 -o $@
	
$(OBJD)/glapack.o : $(MKLBLAS) glapack.f90
	$(F90) $(FCOPT) $(INC) glapack.f90 -o $@

$(OBJD)/mkl_blas.o : $(MKLBLAS)
	$(F90) $(FCOPT) $(INC) $(MKLBLAS) -o $@

$(OBJD)/fortran.o : $(CUDAFRT)/fortran.c
	$(CC) $(CCOPT) $(CUDAFRT)/fortran.c $(CRDAINC) -DCUBLAS_GFORTRAN -o $@

#ifneq ($(F90),pgf90)		
#$(OBJD)/fortran.o : $(CSRC)       # for getting fortran interface for cublas
#	$(CC) $(CCOPT) $(CSRC) $(INCC) -DCUBLAS_GFORTRAN -o $(OBJD)/fortran.o 
#endif

install :
ifeq ($(VER),GPU)
	cp libglapack_gpu.so $(INSTLIBDIR)
	cp $(JDIR)/glapack.mod $(INSTINCDIR)
else
	cp libglapack_cpu.so $(INSTLIBDIR)
	cp $(JDIR)/glapack.mod $(INSTINCDIR)
endif

depend .depend :
ifeq ($(VER),GPU)
	makedepf90 -b $(OBJD) -DMAGMA $(FSRC)  > .depend
else
	makedepf90 -b $(OBJD) $(FSRC)  > .depend
endif

clean :
	rm -f libglapack_*.so libgblas_*.so $(OBJD)/*.mod *.mod $(OBJ)

include .depend
