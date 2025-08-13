##########################################################
########## Select F90 compiler and flags  (Mandatory) 
##########################################################

#BIND = bin_86
OBJD = obj_86
JDIR = obj_86
IDIR = -I$(JDIR)

CUDA_VER = 12.9
MKL_VER = 2025.0
MAGMA_VER = 2.8.0
NV_VER = 25.7

CUDAROOT = $(HOME)/local/cuda-$(CUDA_VER)
CUDAHOME = $(CUDAROOT)/targets/x86_64-linux
CUDALIB = -L$(CUDAHOME)/lib -lcudart
CUDAINC = -I$(CUDAHOME)/include

MAGMAROOT = $(HOME)/local/magma-$(MAGMA_VER)
MAGMAINC = -I$(MAGMAROOT)/include/$(INSD)
MAGMALIB = -L$(MAGMAROOT)/lib -lmagma
MKLROOT = $(HOME)/app/oneapi/mkl/$(MKL_VER)
#NVLIB = -L/opt/nvidia/hpc_sdk/Linux_x86_64/$(NV_VER)/compilers/lib -lgomp

# for profiler
#PG = -pg
PG =

# chose fortran compiler
#F90 = ifort
F90 = ifx
#F90 = gfortran
#F90 = pgf90
#F90 = nvfortran

# chose c compiler
#CC = gcc
#CC = icc
CC = icx
#CC = pgcc
#CC = nvcc

ifeq ($(F90),ifort)
  INSD = intel
else ifeq ($(F90),ifx)
  INSD = intel
else ifeq ($(F90),gfortran)
  INSD = gnu
else ifeq ($(F90),nvfortran)
  INSD = nvidia
else ifeq ($(F90),pgf90)
  INSD = pgi
endif

FGSLINC = -I$(HOME)/local/include/$(INSD)
FGSLLIB = -L$(HOME)/local/lib/$(INSD) -lfgsl

GSL = $(HOME)/local/lib # own libgsl and libgslcblas
GSLLIB = -L$(GSL) -lgsl -lgslcblas

ifeq ($(F90),ifort)
   INSD = intel
   FPIC = -fPIC
   FPP = -fpp
   MJ = -module $(JDIR)
   MKL = -L$(MKLROOT)/lib/intel64 -L$(HOME)/lib
   MKLBLAS = $(MKLROOT)/include/mkl_spblas.f90
   MKLINC = -I$(MKLROOT)/include
else ifeq ($(F90),ifx)
   INSD = intel
   FPIC = -fPIC
   FPP = -fpp
   MJ = -module $(JDIR)
   MKL = -L$(MKLROOT)/lib/intel64 -L$(HOME)/lib
   MKLBLAS = $(MKLROOT)/include/mkl_spblas.f90
   MKLINC = -I$(MKLROOT)/include
else ifeq ($(F90),gfortran)
   INSD = gnu
   FPIC = -fPIC
   FPP = -cpp
   MJ = -J$(JDIR)
else ifeq ($(F90),pgf90)
   INSD = pgi
   FPIC = -fPIC #FPIC = -Xcompiler -fPIC
   FPP = -Mpreprocess
   MJ = -module $(JDIR)
else ifeq ($(F90),nvfortran)
   INSD = nvidia
   FPIC = -fPIC #FPIC = -Xcompiler -fPIC  
   FPP = -Mpreprocess
   MJ = -module $(JDIR)
endif

ifeq ($(CC),nvcc)
   CFLAGS = -Xcompiler -fPIC
   COPT = -O2
else ifeq ($(CC),g++)
   CFLAGS = -fPIC -fopenmp
   COPT = -O2
else ifeq ($(CC),icc)
   CFLAGS = -fPIC -qopenmp
   COPT = -O2 -diag-disable=10448
else ifeq ($(CC),icx)
   CFLAGS = -fPIC -qopenmp
   COPT = -O2
endif

MKL = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

# for issp
#INSTLIBDIR=$(HOME)/local/lib/$(INSD)
#INSTINCDIR=$(HOME)/local/include/$(INSD)

# for PC
INSTLIBDIR=/usr/local/lib/$(INSD)
INSTINCDIR=/usr/local/include/$(INSD)

LIB =  $(FGSLLIB) $(GSLLIB) $(MKL) $(LAPACK) #$(NVLIB)
INC = $(FGSLINC) $(MKLINC) -I$(JDIR) # for cpu version fgsl.mod

# *****at the moment, openmp does not work*****

# optimization flag
ifeq ($(F90),ifort)
    #OPT = -g -qopenmp -check all -traceback -diag-disable=10448
    OPT = -O3 -qopenmp -diag-disable=10448
    SHARED = -shared
endif

ifeq ($(F90),ifx)
    #OPT = -g -qopenmp -check all -traceback
    OPT = -O3 -qopenmp
    SHARED = -shared
endif

ifeq ($(F90),gfortran)
    #OPT = -g -fcheck=all
    OPT = -O3 -fopenmp
    #OPT = -O3
    SHARED = -shared
endif

ifeq ($(F90),pgf90)
   OPT = -O3  -mp
   SHARED = -shared
endif

ifeq ($(F90),nvfortran) 
   OPT = -O3 -mp
   SHARED = -shared
endif

FCOPT = -c $(OPT) $(PG) $(FPP) $(FPIC) $(IDIR) $(MJ)
CCOPT = -c $(COPT) $(PG) $(CFLAGS) 

