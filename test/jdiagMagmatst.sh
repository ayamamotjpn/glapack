#for gfortran
#INSD=gnu
#for ifort
INSD=intel

MKLVER=2022.1.0
#CUDAVER=11.8
CUDAVER=12.0
CUDALIB=/home/akiji/local/cuda-$CUDAVER/targets/x86_64-linux/lib
MKLLIB=/opt/intel/oneapi/mkl/$MKLVER/lib/intel64
MAGMALIB = /usr/local/magma-2.7.0_intel/lib
GLAPACK=../lib
export LD_LIBRARY_PATH=$GLAPACK:/usr/local/lib:/usr/local/lib/$INSD:$MKLLIB:$CUDALIB:$MAGMALIB:$LD_LIBRARY_PATH

jdiagMagmatst 10 1

jdiagMagmatst 10 2


