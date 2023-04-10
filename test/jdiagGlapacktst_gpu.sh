#for gfortran
#INSD=gnu
#for ifort
INSD=intel

MKLVER=2022.1.0
#CUDAVER=11.8
#MKLVER=2023.0.0
export CUDAVER=12.0
export CUDALIB=$HOME/local/cuda-$CUDAVER/targets/x86_64-linux/lib
export MKLLIB=/opt/intel/oneapi/mkl/$MKLVER/lib/intel64
#GLAPACK=../lib
export GLAPACK=/usr/local/lib/$INSD
export LD_LIBRARY_PATH=$GLAPACK:/usr/local/lib:/usr/local/lib/$INSD:$MKLLIB:$CUDALIB:$LD_LIBRARY_PATH

./jdiagGlapacktst_gpu 10 1

#./jdiagGlapacktst_gpu 10 2


# at the moment intel version does not work
#lrwork invalid
