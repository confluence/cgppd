################################################################################

INCLUDE=-I/usr/include/glog/ -I/usr/include/gflags/
#LIBS=-L/usr/lib64
LIBS=-L/usr/lib
LINKS=-lpthread -lgsl -lgslcblas -lglog
CFLAGS=-fno-omit-frame-pointer -O2 -std=c++11
COMPILER=g++
HGVERSION:= $(shell hg parents --template 'hgid: {node|short}')
DEFINE=
OBJS=AminoAcid AminoAcids Potential Molecule Quaternion Replica Residue Simulation TorsionalLookupMatrix vector3f Geometry

# For automatic dependencies
SUFFIXES+=.d
SOURCES=$(shell find src/ -name "*.cpp")
DEPFILES=$(patsubst %.cpp, %.d, $(SOURCES))
NODEPS=clean help


### the following line enables debug output and emulation
# TODO maybe enable this when debug is on
# NVCC_COMPILER_FLAGS=-g -deviceemu -D_EMU
# NVCC_COMPILER_FLAGS=-g
NVCC_ARCH=-arch=sm_20

################################################################################

# Compilation defaults

CUDA=yes
GL=yes
STREAMS=no
LINKERS=yes
LJ=normal

TEST=yes
DEBUG=no

################################################################################

ifeq ($(CUDA),yes)
#INCLUDE+=-I/$(HOME)/NVIDIA_GPU_Computing_SDK/C/common/inc -I/usr/include/nvidia-current/cuda -I/usr/local/cuda/include
INCLUDE+=-I/$(HOME)/NVIDIA_CUDA_Samples/common/inc -I/usr/local/cuda/include
#LIBS+=-L/$(HOME)/NVIDIA_GPU_Computing_SDK/C/lib -L/usr/lib64/nvidia-current -L/usr/local/cuda/lib64
LIBS+=-L/$(HOME)/NVIDIA_CUDA_Samples/common/lib -L/usr/local/cuda/lib64
#LINKS+=-lcudart -lcutil -lcuda -lglog -lgflags
LINKS+=-lcudart -lcuda -lgflags
OBJS+=CudaFunctions
DEFINE+=-DEnableCUDA
endif

ifeq ($(GL),yes)
LINKS+=-lglut -lGL -lGLU
OBJS+=Camera openglvis
DEFINE+=-DEnableOPENGL
endif

# Always create these, because we might make test with TEST=no (although that would be silly)
TEST_INCLUDE=${INCLUDE} -I/usr/include/cppunit/ -Isrc -Itests
TEST_LINKS=${LINKS} -lcppunit
TEST_SOURCES:=$(shell find tests/ -regex '.*\.\(cpp\|h\)')
# suppress warnings about conversion from string constant to char * when constructing argv in tests
TEST_FLAGS=-Wno-write-strings

ifeq ($(DEBUG),yes)
CFLAGS+=-g
endif

ifeq ($(STREAMS),yes)
DEFINE+=-DEnableStreams
endif

ifeq ($(LINKERS),yes)
DEFINE+=-DEnableFlexibleLinkers
endif

ifeq ($(LJ),off)
DEFINE+=-DEnableLJOff
endif

ifeq ($(LJ),repulsive)
DEFINE+=-DEnableLJRepulsive
endif

OBJFILES=$(patsubst %, obj/%.o, $(OBJS))

################################################################################

ifneq ($(TEST),yes)
cgppd: obj/main.o ${OBJFILES}
else
cgppd: obj/main.o ${OBJFILES} test
endif
	${COMPILER} ${INCLUDE} ${DEFINE} ${CFLAGS} ${LIBS} -o $@ obj/main.o ${OBJFILES} ${LINKS}

test: ${OBJFILES} ${TEST_SOURCES}
	${COMPILER} ${TEST_INCLUDE} ${DEFINE} ${CFLAGS} ${TEST_FLAGS} ${LIBS} -o test ${OBJFILES} ${TEST_SOURCES} ${TEST_LINKS}

obj/CudaFunctions.o: src/CudaFunctions.cu src/CudaFunctions.h
	@echo Making CUDA files.
	nvcc ${NVCC_ARCH} ${NVCC_COMPILER_FLAGS} ${INCLUDE} ${DEFINE} -c src/CudaFunctions.cu -o obj/CudaFunctions.o

ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
-include $(DEPFILES)
endif

src/%.d: src/%.cpp
	$(COMPILER) $(CFLAGS) -MM -MT '$(patsubst src/%.cpp,obj/%.o,$<)' $< -MF $@

obj/main.o: src/main.cpp
	@mkdir -p $(dir $@)
	$(COMPILER) $(CFLAGS) ${INCLUDE} ${DEFINE} -DHGVERSION="\"${HGVERSION}\"" -o $@ -c $<

obj/%.o: src/%.cpp src/%.d src/%.h
	@mkdir -p $(dir $@)
	$(COMPILER) $(CFLAGS) ${INCLUDE} ${DEFINE} -o $@ -c $<

clean:
	@rm -rf obj cgppd test

help:
	@echo "Usage: make [cgppd] [options]"
	@echo "Options:"
	@echo "  CUDA: use CUDA (default: yes)"
	@echo "  GL: use GL (default: yes)"
	@echo "  STREAMS: use CUDA streams (default: no)"
	@echo "  LINKERS: use flexible linkers (default: yes)"
	@echo "  TEST: build unit tests (default: yes)"
	@echo "  DEBUG: compile with debugging symbols (default: no)"

.PHONY: help clean
