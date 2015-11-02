################################################################################

INC=-Iinc
LINKS=-lpthread -lgsl -lgslcblas
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
NVCC_ARCH=-arch=sm_20

################################################################################

# Compilation defaults

CUDA=yes
GL=no
LINKERS=yes
LJ=normal

TEST=yes
DEBUG=no

################################################################################

ifeq ($(CUDA),yes)
LINKS+=-lcudart -lcuda
OBJS+=CudaFunctions
DEFINE+=-DEnableCUDA -DEnableStreams
endif

ifeq ($(GL),yes)
LINKS+=-lglut -lGL -lGLU
OBJS+=Camera openglvis
DEFINE+=-DEnableOPENGL
endif

# Always create these, because we might make test with TEST=no (although that would be silly)
TEST_INC=${INC} -Isrc -Itests
TEST_LINKS=${LINKS}
TEST_SOURCES:=$(shell find tests/ -regex '.*\.\(cpp\|h\)')
# suppress warnings about conversion from string constant to char * when constructing argv in tests
TEST_FLAGS=-Wno-write-strings

ifeq ($(DEBUG),yes)
CFLAGS+=-g
endif

ifeq ($(LINKERS),yes)
DEFINE+=-DEnableFlexibleLinkers
endif

ifeq ($(LJ),normal)
APPNAME=cgppd
TESTNAME=test
POLYMERTEST=no
endif

ifeq ($(LJ),off)
DEFINE+=-DEnableLJOff
APPNAME=cgppd_ljoff
TESTNAME=test_ljoff
POLYMERTEST=yes
endif

ifeq ($(LJ),repulsive)
DEFINE+=-DEnableLJRepulsive
APPNAME=cgppd_ljrep
TESTNAME=test_ljrep
POLYMERTEST=yes
endif

ifeq ($(POLYMERTEST),yes)
DEFINE+=-DEnablePolymerTest
endif

OBJFILES=$(patsubst %, obj/%.o, $(OBJS))

################################################################################

ifneq ($(TEST),yes)
${APPNAME}: obj/main.o ${OBJFILES}
else
${APPNAME}: obj/main.o ${OBJFILES} ${TESTNAME}
endif
	${COMPILER} ${INC} ${DEFINE} ${CFLAGS} -o $@ obj/main.o ${OBJFILES} ${LINKS}

${TESTNAME}: ${OBJFILES} ${TEST_SOURCES}
	${COMPILER} ${TEST_INC} ${DEFINE} ${CFLAGS} ${TEST_FLAGS} -DHGVERSION="\"${HGVERSION}\"" -o $@ ${OBJFILES} ${TEST_SOURCES} ${TEST_LINKS}

obj/CudaFunctions.o: src/CudaFunctions.cu src/CudaFunctions.h src/cudaExterns.h src/definitions.h src/constants.h src/AminoAcids.h
	@echo Making CUDA files.
	nvcc ${NVCC_ARCH} ${NVCC_COMPILER_FLAGS} ${INC} ${DEFINE} -c src/CudaFunctions.cu -o obj/CudaFunctions.o

ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
-include $(DEPFILES)
endif

src/%.d: src/%.cpp
	$(COMPILER) $(CFLAGS) ${INC} ${DEFINE} -MM -MT '$(patsubst src/%.cpp,obj/%.o,$<)' $< -MF $@

obj/main.o: src/main.cpp
	@mkdir -p $(dir $@)
	$(COMPILER) $(CFLAGS) ${INC} ${DEFINE} -DHGVERSION="\"${HGVERSION}\"" -o $@ -c $<

obj/%.o: src/%.cpp src/%.d src/%.h
	@mkdir -p $(dir $@)
	$(COMPILER) $(CFLAGS) ${INC} ${DEFINE} -DHGVERSION="\"${HGVERSION}\"" -o $@ -c $<

clean:
	@rm -rf obj $(DEPFILES)

help:
	@echo "Usage: make [cgppd] [options]"
	@echo "Options:"
	@echo "  CUDA: use CUDA, with streams (default: yes)"
	@echo "  GL: use GL (default: no)"
	@echo "  LINKERS: use flexible linkers (default: yes)"
	@echo "  TEST: build unit tests (default: yes)"
	@echo "  DEBUG: compile with debugging symbols (default: no)"
	@echo "  LJ=normal/off/repulsive: compile with modified LJ potentials for polymer tests (default: normal)"

.PHONY: help clean
