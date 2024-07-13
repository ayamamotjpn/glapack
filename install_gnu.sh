export INST=/usr
cp ./lib/libglapack_gpu.so $INST/local/lib/gnu
cp ./obj_86/glapack.mod $INST/local/include/gnu
cp ./lib/libgblas_gpu.so $INST/local/lib/gnu
cp ./obj_86/cublasxt.mod $INST/local/include/gnu
