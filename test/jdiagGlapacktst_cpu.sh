#for gfortran
INSD=gnu
#for ifort
#INSD=intel

MKLVER=2022.1.0
CUDAVER=11.8
#CUDAVER=12.0
CUDALIB=$HOME/local/cuda-$CUDAVER/targets/x86_64-linux/lib
MKLLIB=/opt/intel/oneapi/mkl/$MKLVER/lib/intel64
OMPLIB=/opt/intel/oneapi/compiler/$MKLVER/linux/compiler/lib/intel64_lin
GLAPACK=../lib
export LD_LIBRARY_PATH=$GLAPACK:/usr/local/lib:/usr/local/lib/$INSD:$MKLLIB:$OMPLIB:$LD_LIBRARY_PATH

jdiagGlapacktst_cpu 10 1

