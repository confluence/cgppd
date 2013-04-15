################################################################################

INCLUDE=
LIBS=-L/usr/lib64
LINKS=-lpthread -lgsl -lgslcblas
CFLAGS=-fno-omit-frame-pointer -O2
COMPILER=g++
HGVERSION:= $(shell hg parents --template 'hgid: {node|short}')
DEFINE=
OBJS=AminoAcid AminoAcids Link Segment Potential Molecule Quaternion Replica Residue Simulation TorsionalLookupMatrix vector3f

# For automatic dependencies
SUFFIXES+=.d
SOURCES=$(shell find src/ -name "*.cpp")
DEPFILES=$(patsubst %.cpp, %.d, $(SOURCES))
NODEPS=clean help


### the following line enables debug output and emulation
# TODO maybe enable this when debug is on
#NVCC_COMPILER_FLAGS=-g -deviceemu -D_EMU
#NVCC_COMPILER_FLAGS=-g
NVCC_ARCH=-arch=sm_20

################################################################################

# Compilation defaults

CUDA=yes
GL=yes
STREAMS=no
LINKERS=yes

TEST=yes
DEBUG=no

################################################################################

ifeq ($(CUDA),yes)
INCLUDE+=-I/$(HOME)/NVIDIA_GPU_Computing_SDK/C/common/inc -I/usr/include/nvidia-current/cuda -I/usr/local/cuda/include
LIBS+=-L/$(HOME)/NVIDIA_GPU_Computing_SDK/C/lib -L/usr/lib64/nvidia-current -L/usr/local/cuda/lib64
LINKS+=-lcudart -lcutil -lcuda
OBJS+=CudaFunctions
DEFINE+=-DEnableCUDA
endif

ifeq ($(GL),yes)
LINKS+=-lglut -lGL -lGLU
OBJS+=Camera OpenGLController
DEFINE+=-DEnableOPENGL
endif

ifeq ($(TEST),yes)
TEST_INCLUDE=${INCLUDE} -I/usr/include/cppunit/ -Isrc
TEST_LINKS=${LINKS} -lcppunit
TEST_SOURCES:=$(shell find tests/ -name "*.cpp")
endif

ifeq ($(DEBUG),yes)
CFLAGS+=-g
endif

ifeq ($(STREAMS),yes)
DEFINE+=-DEnableStreams
endif

ifeq ($(LINKERS),yes)
DEFINE+=-DEnableFlexibleLinkers
endif

OBJFILES=$(patsubst %, obj/%.o, $(OBJS))

################################################################################

ifneq ($(TEST),yes)
REMC: obj/main.o ${OBJFILES}
	${COMPILER} ${INCLUDE} ${DEFINE} ${CFLAGS} ${LIBS} -o $@ $^ ${LINKS}
else
REMC: obj/main.o ${OBJFILES} ${TEST_SOURCES}
	${COMPILER} ${INCLUDE} ${DEFINE} ${CFLAGS} ${LIBS} -o $@ obj/main.o ${OBJFILES} ${LINKS}
	${COMPILER} ${TEST_INCLUDE} ${DEFINE} ${CFLAGS} ${LIBS} -o test ${OBJFILES} ${TEST_SOURCES} ${TEST_LINKS}
endif

obj/CudaFunctions.o: src/CudaFunctions.cu src/CudaFunctions.h
	@echo Making CUDA files.
	nvcc ${NVCC_ARCH} ${NVCC_COMPILER_FLAGS} ${INCLUDE} ${DEFINE} -c src/CudaFunctions.cu -o obj/CudaFunctions.o

ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
-include $(DEPFILES)
endif

src/%.d: src/%.cpp
	$(COMPILER) $(CFLAGS) -MM -MT '$(patsubst src/%.cpp,obj/%.o,$<)' $< -MF $@

obj/main.o: src/main.cpp src/main.h
	@mkdir -p $(dir $@)
	$(COMPILER) $(CFLAGS) ${INCLUDE} ${DEFINE} -DHGVERSION="\"${HGVERSION}\"" -o $@ -c $<

obj/%.o: src/%.cpp src/%.d src/%.h
	@mkdir -p $(dir $@)
	$(COMPILER) $(CFLAGS) ${INCLUDE} ${DEFINE} -o $@ -c $<

clean:
	@rm -rf obj REMC test

help:
	@echo "Usage: make [REMC] [options]"
	@echo "Options:"
	@echo "  CUDA: use CUDA (default: yes)"
	@echo "  GL: use GL (default: yes)"
	@echo "  STREAMS: use CUDA streams (default: no)"
	@echo "  LINKERS: use flexible linkers (default: yes)"
	@echo "  TEST: build unit tests (default: yes)"
	@echo "  DEBUG: compile with debugging symbols (default: no)"

.PHONY: help clean
