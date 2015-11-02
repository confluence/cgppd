module add cuda/Cuda-6.5
module add compilers/gcc484
# Add samples
export INCLUDE=/usr/local/cuda-6.5/samples/common/inc:$INCLUDE
export LD_LIBRARY_PATH=/usr/local/cuda-6.5/samples/common/lib:$LD_LIBRARY_PATH
# I don't know why these aren't set automatically, but we need them
export CPATH=$INCLUDE
export LIBRARY_PATH=$LD_LIBRARY_PATH
