module add cuda/Cuda-6.5
module add compilers/gcc484
export INCLUDE=/usr/local/cuda-6.5/samples/common/inc:$INCLUDE
export LD_LIBRARY_PATH=/usr/local/cuda-6.5/samples/common/lib:$LD_LIBRARY_PATH
export CPATH=$INCLUDE
