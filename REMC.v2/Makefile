# TODO: set up machine-specific defaults in a separate file
INCLUDE=
LIBS=-L/usr/lib64
LINKS=-lpthread -lgsl -lgslcblas
CFLAGS=-fno-omit-frame-pointer -O2
COMPILER=g++
HGVERSION:= $(shell hg parents --template 'hgid: {node|short}')
DEFINE=

OBJS=obj/AminoAcid.o obj/AminoAcids.o obj/Link.o obj/Segment.o obj/Potential.o obj/Molecule.o obj/Quaternion.o obj/Replica.o obj/Residue.o obj/TorsionalLookupMatrix.o obj/vector3f.o

TEST_SOURCES:=$(shell find tests/ -name "*.cpp")

ifeq ($(CUDA),yes)
INCLUDE+=-I/$(HOME)/NVIDIA_GPU_Computing_SDK/C/common/inc -I/usr/include/nvidia-current/cuda -I/usr/local/cuda/include
LIBS+=-L/$(HOME)/NVIDIA_GPU_Computing_SDK/C/lib -L/usr/lib64/nvidia-current -L/usr/local/cuda/lib64
LINKS+=-lcudart -lcutil -lcuda
OBJS+=obj/CudaFunctions.o
DEFINE+=-DEnableCUDA
endif

ifeq ($(GL),yes)
LINKS+=-lglut -lGL -lGLU
OBJS+=obj/Camera.o
DEFINE+=-DEnableOPENGL
endif

ifeq ($(TEST),yes)
TEST_INCLUDE=${INCLUDE} -I/usr/include/cppunit/ -Isrc
TEST_LINKS=${LINKS} -lcppunit
endif

ifeq ($(DEBUG),yes)
CFLAGS+=-g
endif

ifeq ($(STREAMS),yes)
DEFINE+=-DEnableStreams
endif

### the following line enables debug output and emulation
# TODO maybe enable this when debug is on
#NVCC_COMPILER_FLAGS=-g -deviceemu -D_EMU
#NVCC_COMPILER_FLAGS=-g
NVCC_ARCH=-arch=sm_20

help:
	@echo "Usage: make REMC [CUDA=yes] [GL=yes] [STREAMS=yes] [TEST=yes] [DEBUG=yes]"

.PHONY: help clean

ifneq ($(TEST),yes)
REMC: obj/main.o ${OBJS}
	${COMPILER} ${INCLUDE} ${DEFINE} ${CFLAGS} ${LIBS} -o $@ $^ ${LINKS}
else
REMC: obj/main.o ${OBJS} ${TEST_SOURCES}
	${COMPILER} ${INCLUDE} ${DEFINE} ${CFLAGS} ${LIBS} -o $@ obj/main.o ${OBJS} ${LINKS}
	${COMPILER} ${TEST_INCLUDE} ${DEFINE} ${CFLAGS} ${LIBS} -o test ${OBJS} ${TEST_SOURCES} ${TEST_LINKS}
endif

obj/CudaFunctions.o: src/CudaFunctions.cu src/CudaFunctions.h
	@echo Making CUDA files.
	nvcc ${NVCC_ARCH} ${NVCC_COMPILER_FLAGS} ${INCLUDE} ${DEFINE} -c src/CudaFunctions.cu -o obj/CudaFunctions.o

SUFFIXES+=.d
SOURCES=$(shell find src/ -name "*.cpp")
DEPFILES=$(patsubst %.cpp,%.d,$(SOURCES))
NODEPS=clean help

# TODO: make this also work for blank target
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
	@rm -rf obj src/*.d REMC test
