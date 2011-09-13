INCLUDE=-I/$(HOME)/NVIDIA_GPU_Computing_SDK/C/common/inc -I/usr/include/nvidia-current/cuda -I/usr/local/cuda/include
LIBS=-L/$(HOME)/NVIDIA_GPU_Computing_SDK/C/lib -L/usr/lib64/nvidia-current -L/usr/local/cuda/lib64  -L/usr/lib64
LINKS=-lpthread -lgsl -lgslcblas -lcudart -lcutil -lcuda
#TAU_MAKEFILE=$(HOME)/tau2/x86_64/lib/Makefile.tau-pthread-pdt
#COMPILER=tau_cxx.sh -tau_makefile=${TAU_MAKEFILE} -optPDTInst
#COMPILER=g++-4.3
COMPILER=g++
#config for zero.ch.cam.ac.uk
#INCLUDE=-I/usr/local/shared/rhel-4.5/x86_64/cuda/2.3/include -I/usr/local/shared/rhel-4.5/x86_64/cuda/2.3/C/common/inc -I/$(HOME)/include
#LIBS=-L/usr/local/shared/rhel-4.5/x86_64/cuda/2.3/cuda/lib64 -L/usr/local/shared/rhel-4.5/x86_64/cuda/2.3/C/lib -L/$(HOME)/lib

#set CUDA_PROFILE=1

#DEFINE=-DEnableOPENGL -DEnableCUDA -DEnableStreams

CFLAGS=-fno-omit-frame-pointer -O2
# CFLAGS+=-ffast-math
### uncomment to add debug symbols
CFLAGS+=-g

### the following line enables debug output and emulation
#NVCC_COMPILER_FLAGS=-g -deviceemu -D_EMU
#NVCC_COMPILER_FLAGS=-g
NVCC_ARCH=-arch=sm_20

COMMON_OBJS=obj/AminoAcid.o obj/AminoAcids.o obj/Link.o obj/Molecule.o obj/Quaternion.o obj/Replica.o obj/Residue.o obj/TorsionalLookupMatrix.o obj/vector3f.o

1:
	make DEFINE='-DEnableCUDA -DEnableOPENGL' GPUGL test

#if using ubuntu: link to -lgslcblas for gsl rng types
help:
	@echo Make using the following arguments
	@echo "  make DEFINE='-DEnableCUDA' GPU"
	@echo "  make DEFINE='-DEnableCUDA -DEnableStreams' GPU # for async functionality"
	@echo "  make CPU # for CPU-only code"
	@echo "  make DEFINE='-DEnableCUDA -DEnableOPENGL' GPUGL # for OPENGL; don't use it for simulations"

.PHONY: help clean 1

GPU: obj/main.o ${COMMON_OBJS} obj/CudaFunctions.o #CudaFunctions.cubin
	@echo Linking GPU Executable: GPU
	${COMPILER} ${INCLUDE} ${DEFINE} ${CFLAGS} ${LIBS} -o GPU obj/*.o ${LINKS}

# compile this with  DEFINE=-DEnableCUDA -DEnableOPENGL GPU
GPUGL: obj/main.o ${COMMON_OBJS} obj/Camera.o obj/CudaFunctions.o #CudaFunctions.cubin
	@echo Linking GPU Executable: GPUGL
	${COMPILER} ${INCLUDE} ${DEFINE} ${CFLAGS} ${LIBS} -o GPUGL obj/*.o ${LINKS} -lglut -lGL -lGLU


CPU: obj/main.o ${COMMON_OBJS} obj/Camera.o
	@echo Linking CPU Executable: CPU
	${COMPILER} ${INCLUDE} ${DEFINE} ${CFLAGS} ${LIBS} -o CPU obj/*.o ${LINKS}

# TODO: split these out to have separate dependencies?
test: ${COMMON_OBJS} obj/Camera.o obj/CudaFunctions.o tests/all.cpp tests/testReplica.cpp tests/testMolecule.cpp
	${COMPILER} ${INCLUDE} -I/usr/include/unittest++ -Isrc ${DEFINE} ${CFLAGS} ${LIBS} -o test ${COMMON_OBJS} obj/Camera.o obj/CudaFunctions.o tests/*.cpp ${LINKS} -lUnitTest++ -lglut -lGL -lGLU

#DRMS: AminoAcid.o AminoAcids.o Camera.o Link.o Molecule.o Quaternion.o Replica.co Residue.o rngfunctions.o TorsionalLookupMatrix.o vector3f.o drms.o
#	@echo Making DRMS Executable: DRMS
#	g++ ${INCLUDE} ${DEFINE} ${CFLAGS} -o DRMS obj/*.o -lpthread -lgsl -lglut
#	rm -f drms.o

obj/AminoAcid.o: obj/.d src/AminoAcid.cpp src/AminoAcid.h
	${COMPILER} ${CFLAGS} -c src/AminoAcid.cpp -o obj/AminoAcid.o

obj/AminoAcids.o: obj/.d src/AminoAcids.cpp src/AminoAcids.h src/AminoAcid.h
	${COMPILER} ${CFLAGS} -c src/AminoAcids.cpp -o obj/AminoAcids.o

obj/Camera.o: obj/.d src/Camera.cpp src/Camera.h src/vector3f.h
	${COMPILER} ${CFLAGS} -c src/Camera.cpp -o obj/Camera.o

obj/Link.o: obj/.d src/Link.cpp src/Link.h
	${COMPILER} ${CFLAGS} -c src/Link.cpp -o obj/Link.o

obj/main.o: obj/.d src/main.cpp src/definitions.h src/AminoAcid.h src/TorsionalLookupMatrix.h src/AminoAcids.h src/Replica.h src/vector3f.h src/Quaternion.h src/Molecule.h src/Residue.h src/Link.h src/cudaExterns.h

	${COMPILER} ${CFLAGS} ${INCLUDE} ${DEFINE} -c src/main.cpp -o obj/main.o

obj/Molecule.o: obj/.d src/Molecule.cpp src/Molecule.h src/definitions.h src/AminoAcids.h src/AminoAcid.h src/vector3f.h src/Quaternion.h src/TorsionalLookupMatrix.h src/Residue.h src/Link.h
	${COMPILER} ${CFLAGS} ${INCLUDE} ${DEFINE} -c src/Molecule.cpp -o obj/Molecule.o

obj/Quaternion.o: obj/.d src/Quaternion.cpp src/Quaternion.h src/vector3f.h
	${COMPILER} ${CFLAGS} -c src/Quaternion.cpp -o obj/Quaternion.o

obj/Replica.o: obj/.d src/Replica.cpp src/Replica.h src/vector3f.h src/definitions.h src/Quaternion.h src/Molecule.h src/AminoAcids.h src/AminoAcid.h src/TorsionalLookupMatrix.h src/Residue.h src/Link.h
	${COMPILER} ${CFLAGS} ${INCLUDE} ${DEFINE} -c src/Replica.cpp -o obj/Replica.o

obj/CudaFunctions.o: obj/.d src/CudaFunctions.cu src/CudaFunctions.h
	@echo Making CUDA files.
	nvcc ${NVCC_ARCH} ${NVCC_COMPILER_FLAGS} ${INCLUDE} ${DEFINE} -c src/CudaFunctions.cu -o obj/CudaFunctions.o

# CudaFunctions.cubin:
# 	nvcc ${NVCC_COMPILER_FLAGS} ${INCLUDE} ${DEFINE} -cubin src/CudaFunctions.cu

obj/Residue.o: obj/.d src/Residue.cpp src/Residue.h src/vector3f.h src/AminoAcids.h src/AminoAcid.h src/definitions.h
	${COMPILER} ${CFLAGS} -c src/Residue.cpp -o obj/Residue.o

obj/TorsionalLookupMatrix.o: obj/.d src/TorsionalLookupMatrix.cpp src/TorsionalLookupMatrix.h src/AminoAcids.h src/AminoAcid.h src/definitions.h
	${COMPILER} ${CFLAGS} -c src/TorsionalLookupMatrix.cpp -o obj/TorsionalLookupMatrix.o

obj/vector3f.o: obj/.d src/vector3f.cpp src/vector3f.h
	${COMPILER} ${CFLAGS} -c src/vector3f.cpp -o obj/vector3f.o

obj/.d:
	@mkdir -p obj
	@touch obj/.d

clean:
	@echo Cleaning...
	@rm -rf obj
	@rm -f GPUGL GPU CPU test
# 	@rm -f *.cubin
	@echo Done.

# clean.o:
# 	@echo Cleaning object files for safety sake.
# 	@rm -f obj/*
