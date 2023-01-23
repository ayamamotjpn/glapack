#VER = GPU
VER = CPU

# for profiler
PG = -pg
PG =

# chose fortran compiler
F90 = ifort
#F90 = gfortran
#F90 = pgf90

# chose c compiler
#CC = gcc
CC = icc
#CC = pgcc
#CC = nvcc
#CFLAGS =

INSTDIRROOT = /home/takemori/local/glapack
CUDADIR = /usr/local/cuda-10.1
CUDALIB = -L$(CUDADIR)/lib64 -lcublas -lcudart
CUDAINC =-I$(CUDADIR)/include -I/opt/pgi/linux86-64-llvm/2019/cuda/10.1/include

MKLINC = -I/opt/intel/mkl/include

#GSL = /usr/lib64          # use system libgsl
GSL = /home/takemori/local/gsl-2.4    # use own libgsl

GSLLIB = -L$(GSL) -lgsl -lgslcblas
MAGMAROOT = /usr/local/magma-2.5.1

ifeq ($(F90),ifort)
   FGSLINC = -I/home/takemori/local/fgsl-1.3.0/include/intel
   FGSLLIB = -L/home/takemori/local/fgsl-1.3.0/lib/intel -lfgsl
   MAGMAINC = -I$(MAGMAROOT)/include/intel -I$(MAGMAROOT)/include  # for intel fortran
   MAGMALIB = -L$(MAGMAROOT)/lib/intel -lmagma
   MKL = -L/opt/intel/mkl/lib/intel64 -L/opt/intel/lib/intel64 -L/usr/lib64 -lblas        # for -lpthread -lm -ldl
   #PARDISO =  -L/home/terast4/share/yamamoto/pardiso -lpardiso500-INTEL1301-X86-64
   INSTLIBDIR=$(INSTDIRROOT)/lib/intel
   INSTDIRINC=$(INSTDIRROOT)/include/intel
endif

ifeq ($(F90),gfortran)
   FGSLINC = -I/usr/local/fgsl/include/gnu
   FGSLLIB = -L/usr/local/fgsl/lib/gnu -lfgsl
   MAGMAINC = -I$(MAGMAROOT)/include/gnu -I$(MAGMAROOT)/include
   MAGMALIB = -L$(MAGMAROOT)/lib/gnu  -lmagma
   #MKL = -L/opt/intel/mkl/lib/intel64 -L/opt/intel/lib/intel64 -L/usr/lib64 -lblas
   MKL = -L/usr/lib64 -lblas
   #PARDISO =  -L/home/terast4/share/yamamoto/pardiso -lpardiso500-INTEL1301-X86-64
   INSTLIBDIR=$(INSTDIRROOT)/lib/gnu
   INSTDIRINC=$(INSTDIRROOT)/include/gnu
endif

ifeq ($(F90),pgf90)
   FGSLINC = -I/usr/local/fgsl/include/pgi
   FGSLLIB = -L/usr/local/fgsl/lib/pgi_-lfgsl
   MAGMAINC = -I$(MAGMAROOT)/include/pgi -I$(MAGMAROOT)/include
   MAGMALIB = -L$(MAGMAROOT)/lib/pgi  -lmagma
   #MKL = -L/opt/intel/mkl/lib/intel64 -L/opt/intel/lib/intel64 -L/usr/lib64 -lblas
   MKL = -L/usr/lib64 -lblas
   # for -DCUDAFOR
   CUBLASINC = -I/opt/pgi/linux86-64/18.10/include
   #PARDISO =  -L/home/terast4/share/yamamoto/pardiso -lpardiso500-INTEL1301-X86-64
   INSTLIBDIR=$(INSTDIRROOT)/lib/pgi
   INSTDIRINC=$(INSTDIRROOT)/include/pgi
endif

ifeq ($(VER),GPU)
   LIB =  $(MAGMALIB) $(CUDALIB) $(FGSLLIB) $(GSLLIB) $(MKL)   # for magma
   INC = $(FGSLINC) $(FEASTINC_GPU) $(MAGMAINC) $(CUBLASINC)   # for gpu version fgsl.mod  FEAST_GPU and MAGMA
else
   LIB =  $(FGSLLIB) $(GSLLIB) $(MKL)    # for ifort  cpu version  this is better for feast
   INC = $(FGSLINC)  $(MKLINC)           # for cpu version fgsl.mod
endif

############################################################################
#  Compilers 
############################################################################ 

#####################################################################
########### Two possible options (choose one)
########### Option 1- (any) Fortran Compiler is used alone
########### Option 2- Intel Fortran + (any) C Compiler to 
###########    remove the Fortran runtime library dependencies 
#######################################################################

#OPTION = 1
OPTION = 2

################################################################
########## Select C compiler and flags (mandatory for option #2)
################################################################
ifeq ($(F90),gfortran)
  FFLAGS =
endif

ifeq ($(F90),ifort)
  FFLAGS =
endif

# *****at the moment, openmp does not work*****

# optimization flag
ifeq ($(F90),ifort)
    FPP=-fpp
    #OPT = -g -check all -traceback
    OPT = -O3 -qopenmp #
    #OPT = -O3
    SHARED = -shared
endif

ifeq ($(F90),gfortran)
    FPP=-cpp
    #OPT = -g -fcheck=all
    OPT = -O3 -fopenmp
    #OPT = -O3
    SHARED = -shared
endif

ifeq ($(F90),pgf90)
    FPP=-Mpreprocess
	ifeq ($(VER),GPU)
      #CUDAFOR = -DCUDAFOR  -Mcuda   # for using cuda fortran but also can use fortran.c
      #CUDAFOR = -Mcuda   # for using cuda fortran but also can use fortran.c
      CUDAFOR =
    else
      CUDAFOR =
    endif
    #OPT = -g  -traceback  $(CUDAFOR)    #-dclchk -chkptr  -chkstk
    OPT = -O3 -mp $(CUDAFOR)
    #OPT = -O3                    $(CUDAFOR)
    SHARED = -shared    # -Bdynamic
endif

DEVPTR = -Dmagma_devptr_t="integer(kind=8)"

ifeq ($(VER),GPU)
   FCOPT = -c $(OPT) $(PG) $(FPP) -DMAGMA $(IDIR) $(FFLAGS)  $(DEVPTR) -fPIC  # for gpu version -pg for gprof
   CCOPT = -c  $(PG) -DCUBLAS_USE_THUNKING  $(CFLAGS)   -fPIC               # for gpu version -pg for gprof
else
   FCOPT = -c $(OPT) $(PG) $(FPP) $(IDIR) $(FFLAGS) $(DEVPTR)   -fPIC       # for cpu version
   CCOPT = -c $(OPT) $(PG)  $(CFLAGS)    -fPIC                              # for cpu version
endif

##########################################################
########## Select F90 compiler and flags  (Mandatory) 
##########################################################

#>>>>> Remark: three examples F1-F2-F3 (choose one or create your own)

#F1- Intel Fortran Compiler (work with option 1 above; Mandatory for Option #2 above)

ifeq ($(F90),ifort)
   ifeq ($(VER),GPU)
      F90FLAGS = $(OPT) -DMAGMA $(FPP)  $(DEVPTR) -fPIC   # for using MAGMA library
   else
      F90FLAGS = $(OPT)  $(FPP) -fPIC
   endif
endif

##F2- Gnu Fortran Compiler (work with option 1 above)
ifeq ($(F90),gfortran)
   ifeq ($(VER),GPU)
      F90FLAGS = $(OPT) -DMAGMA $(FPP) -fPIC   # for using MAGMA library
   else
      F90FLAGS = $(OPT) -m64 -openmp $(FPP) -fPIC
   endif
endif

##F3-Portland Group Fortran Compiler (work with option 1 above)
ifeq ($(F90),pgf90)
   ifeq ($(VER),GPU)
      FMPIFLAGS = $(OPT) -openmp -DMAGMA $(FPP) -fPIC
   else
      F90FLAGS = $(OPT) -mp -openmp -fPIC
   endif
endif

##########################################################
########### For FEAST-MPI (if you wish to compile it) 
########### Select your own MPI implementation 
##########################################################

#>>>>> Remark: three example P1-P2-P3 (for example: mpiifort, mpif90.mpich2, mpif90.openmpi, etc.)

##P1- Intel MPI
#MPI= mpiifort -f90=$(F90)

##P2- mpich2
#MPI    = mpif90.mpich2 -f90=$(F90)
#MPI = mpif90
MPI=ifort           # for intel compiler

##P3- openmpi ...requires shell environment variable "OMPI_FC=$(F90)"
#export OMPI_FC=$(F90) # for BASH shell
#MPI	= mpif90.openmpi 


######### MPI is here associated with the F90 compiler above and compiler directives are activated
######### Do not change the following
PF90 = $(MPI) 
PF90FLAGS = $(F90FLAGS) -DMPI

