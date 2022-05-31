VER = GPU
#VER = CPU

##########################################################
########## Select F90 compiler and flags  (Mandatory) 
##########################################################

#BIND = bin_86
OBJD = obj_86
JDIR = obj_86
IDIR = -I$(JDIR)

# for profiler
PG = -pg
PG =

# chose fortran compiler
F90 = ifort
#F90 = gfortran
#F90 = pgf90
#F90 = nvfortran

# chose c compiler
#CC = gcc
CC = icc
#CC = pgcc
#CC = nvcc
#CFLAGS =

# for using magma.moc in MAGMA ifort version user ifort
ifeq ($(VER),GPU)
   F90 = ifort
   CC = icc
endif

ifeq ($(F90),ifort)
  INSD = intel
else ifeq ($(F90),gfortran)
  INSD = gnu
else ifeq ($(F90),nvfortran)
  INSD = pgi
else ifeq ($(F90),pgf90)
  INSD = pgi
endif

CUDAVER = 11.6
CUDADIR = /usr/local/cuda-$(CUDAVER)
CUDALIB = -L$(CUDADIR)/lib64 -lcublas -lcudart
CUDAINC =-I$(CUDADIR)/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/22.3/compilers/include

FGSLINC = -I/usr/local/include/$(INSD)/fgsl
FGSLLIB = -L/usr/local/lib/$(INSD) -lfgsl

MKLINC = -I/opt/intel/mkl/include

GSL = /usr/lib64          # use system libgsl for fedora
#GSL = /usr/lib           # use system libgsl for ubuntu
GSLLIB = -L$(GSL) -lgsl -lgslcblas

MAGMAROOT = /usr/local/magma-2.6.2
MAGMAINC = -I$(MAGMAROOT)/include
MAGMALIB = -L$(MAGMAROOT)/lib -lmagma

# for -DCUDAFOR (only for nvfortran or pgf90)
# CUBLASINC = -I//usr/local/cuda-$(CUDAVER)/targets/x86_64-linux/include

ifeq ($(F90),ifort)
   INSD = intel
   FPIC = -fPIC
   FPP = -fpp
   MJ = -J$(JDIR)
   MKL = -L/opt/intel/mkl/lib/intel64 -L/opt/intel/lib/intel64 -L/usr/lib64 -lblas        # for -lpthread -lm -ldl
endif

ifeq ($(F90),gfortran)
   INSD = gnu
   FPIC = -fPIC
   FPP = -cpp
   MJ = -J$(JDIR)
   MKL = -L/usr/lib64 -lblas
endif

ifeq ($(F90),pgf90)
   INSD = pgi
   FPIC = --compiler-options '-fPIC'
   FPP = -Mpreprocess
   MJ = -module $(JDIR)
   MKL = -L/usr/lib64 -lblas
endif

ifeq ($(F90),nvfortran)
   INSD = pgi
   FPIC = -fPIC
   FPP = -Mpreprocess
   MJ = -module $(JDIR)
   MKL = -L/usr/lib64 -lblas
endif


INSTLIBDIR=/usr/local/lib/$(INSD)
INSTINCDIR=/usr/local/include/$(INCD)
MAGMAINC = -I$(MAGMAROOT)/include
MAGMALIB = -L$(MAGMAROOT)/lib -lmagma

ifeq ($(VER),GPU)
   LIB =  $(MAGMALIB) $(CUDALIB) $(FGSLLIB) $(GSLLIB) $(MKL)   # for magma
   INC = $(FGSLINC) $(FEASTINC_GPU) $(MAGMAINC) $(CUBLASINC) -I$(JDIR)  # for gpu version fgsl.mod  FEAST_GPU and MAGMA
else
   LIB =  $(FGSLLIB) $(GSLLIB) $(MKL)    # for ifort  cpu version  this is better for feast
   INC = $(FGSLINC)  $(MKLINC) -I$(JDIR) # for cpu version fgsl.mod
endif

# *****at the moment, openmp does not work*****

# optimization flag
ifeq ($(F90),ifort)
    #OPT = -g -check all -traceback
    OPT = -O3 -qopenmp
    #OPT = -O3
    SHARED = -shared
endif

ifeq ($(F90),gfortran)
    #OPT = -g -fcheck=all
    OPT = -O3 -fopenmp
    #OPT = -O3
    SHARED = -shared
endif

ifeq ($(F90),pgf90)
	ifeq ($(VER),GPU)
      #CUDAFOR = -DCUDAFOR  -Mcuda   # for using cuda fortran but also can use fortran.c
      #CUDAFOR = -Mcuda   # for using cuda fortran but also can use fortran.c
      CUDAFOR =
      #OPT = -g  -traceback  $(CUDAFOR)    #-dclchk -chkptr  -chkstk
      OPT = -O3 -mp $(CUDAFOR)
    else
      CUDAFOR =
      OPT = -O3 
    endif
    SHARED = -shared    # -Bdynamic
endif

ifeq ($(F90),nvfortran)
   ifeq ($(VER),GPU)
      #CUDAFOR = -DCUDAFOR  -Mcuda   # for using cuda fortran but also can use fortran.c
      #CUDAFOR = -Mcuda   # for using cuda fortran but also can use fortran.c
      CUDAFOR =
      #OPT = -g  -traceback  $(CUDAFOR)    #-dclchk -chkptr  -chkstk
      OPT = -O3 -mp $(CUDAFOR)
    else
      CUDAFOR =  
      OPT = -O3
    endif
    SHARED = -shared    # -Bdynamic
endif


DEVPTR = -Dmagma_devptr_t="integer(kind=8)"
FCOPT = -c $(OPT) $(PG) $(FPP) $(FPIC) $(IDIR) $(MJ) $(DEVPTR) -fPIC
CCOPT = -c $(OPT) $(PG) $(CFLAGS) $(FPIC) 

ifeq ($(VER),GPU)
   FCOPT += -DMAGMA
   CCOPT += -DCUBLAS_USE_THUNKING
endif

#PF90 = $(MPI) 
#PF90FLAGS = $(F90FLAGS) -DMPI

#ifeq ($(F90),ifort)
#   ifeq ($(VER),GPU)
#      F90FLAGS = $(OPT) -DMAGMA $(FPP)  $(DEVPTR) $(FPIC)   # for using MAGMA library
#   else
#      F90FLAGS = $(OPT)  $(FPP) -fPIC
#   endif
#endif
#
#ifeq ($(F90),gfortran)
#   ifeq ($(VER),GPU)
#      F90FLAGS = $(OPT) -DMAGMA $(FPP) $(FOIC)   # for using MAGMA library
#   else
#      F90FLAGS = $(OPT) -m64 -openmp $(FPP) $(FOIC)
#   endif
#endif
#
#ifeq ($(F90),pgf90)
#   ifeq ($(VER),GPU)
#      FMPIFLAGS = $(OPT) -DMAGMA $(FPP) $(FPIC)
#   else
#      F90FLAGS = $(OPT) -mp -openmp $(FPIC)
#   endif
#endif
#
#
#ifeq ($(F90),nvfortran)
#   ifeq ($(VER),GPU)
#      FMPIFLAGS = $(OPT) -DMAGMA $(FPP) $(FPIC)
#   else
#      F90FLAGS = $(OPT) -mp -openmp $(FPIC)
#   endif
#endif
