# NOTE: needs boost, tclap, and sdsl

CXX=g++ # clang++-3.8 # g+
CPP_FLAGS=-m64 -std=c++0x -Wall 
DEP_PATH=/usr/local
INC_PATH=$(DEP_PATH)/include
INC_PATH_DYN=DYNAMIC/include/ -I hopscotch-map/include
LIB_PATH=$(DEP_PATH)/lib
DEP_FLAGS=-I$(INC_PATH)/ -I$(INC_PATH_DYN)/ -L$(LIB_PATH)/ -lsdsl # -ldivsufsort -ldivsufsort64
DEBUG_FLAGS=-g
NDEBUG_FLAGS=-DNDEBUG
OPT_FLAGS=-O3 -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native
NOPT_FLAGS=-Og
DYN_FLAGS=-DDYNAMIC

# Using Semantic Versioning: http://semver.org/
VERSION=0.1.0
CPP_FLAGS+=-DVERSION=\"$(VERSION)\"

ifeq ($(dynamic),1)
CPP_FLAGS+=$(DYN_FLAGS)
endif

ifeq ($(optimise),0)
CPP_FLAGS+=$(NOPT_FLAGS)
else
CPP_FLAGS+=$(OPT_FLAGS)
endif

ifeq ($(debug),1)
CPP_FLAGS+=$(DEBUG_FLAGS)
else
CPP_FLAGS+=$(NDEBUG_FLAGS)
endif

ifeq ($(verbose),1)
CPP_FLAGS+=-DVERBOSE
endif

ifneq ($(revcomps),0)
CPP_FLAGS+=-DADD_REVCOMPS
endif

ifneq ($(dummies),0)
CPP_FLAGS+=-DALL_DUMMIES
endif

ifeq ($(varord),1)
CPP_FLAGS+=-DVAR_ORDER
endif

BUILD_REQS=debruijn_graph.hpp io.hpp io.o debug.h
DYN_BUILD_REQS=dynBoss.hpp io.hpp io.o debug.h
ASSEM_REQS=debruijn_graph.hpp algorithm.hpp utility.hpp kmer.hpp uint128_t.hpp
PACK_REQS=lut.hpp debug.h io.hpp io.o sort.hpp kmer.hpp dummies.hpp
BINARIES=cosmo-pack cosmo-build-dyn cosmo-add-verify cosmo-delete-verify cosmo-delete-add

default: all

lut.hpp: make_lut.py
		python make_lut.py > lut.hpp

io.o: io.hpp io.cpp debug.h dummies.hpp kmer.hpp
		$(CXX) $(CPP_FLAGS) -c io.cpp

# TODO: Roll these all into one... "cosmo"
cosmo-pack: cosmo-pack.cpp $(PACK_REQS)
		$(CXX) $(CPP_FLAGS) -o $@ $< io.o

cosmo-build: cosmo-build.cpp $(BUILD_REQS)
		$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(DEP_FLAGS) 

cosmo-build-dyn: cosmo-build-dyn.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS) $(DYN_FLAGS) -o cosmo-build-dyn $< io.o $(DEP_FLAGS)

cosmo-verify: cosmo-verify.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS) -o cosmo-verify $< io.o $(DEP_FLAGS)

#cosmo-assemble: cosmo-assemble.cpp $(ASSEM_REQS) wt_algorithm.hpp debruijn_hypergraph.hpp
#		$(CXX) $(CPP_FLAGS) -o $@ $< $(DEP_FLAGS) 

cosmo-benchmark: cosmo-benchmark.cpp $(ASSEM_REQS) wt_algorithm.hpp debruijn_hypergraph.hpp
		$(CXX) $(CPP_FLAGS) -o $@ $< $(DEP_FLAGS)

cosmo-benchmark-dyn: cosmo-benchmark.cpp $(ASSEM_REQS) wt_algorithm.hpp debruijn_hypergraph.hpp
		$(CXX) $(CPP_FLAGS) $(DYN_FLAGS) -o cosmo-benchmark-dyn $< $(DEP_FLAGS)
cosmo-query: cosmo-query.cpp $(ASSEM_REQS) wt_algorithm.hpp debruijn_hypergraph.hpp
		$(CXX) $(CPP_FLAGS) -o $@ $< $(DEP_FLAGS) 

cosmo-add-delete: cosmo-add-delete.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS) -o cosmo-add-delete $< io.o $(DEP_FLAGS)

cosmo-add-verify: cosmo-add-verify.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS) -o cosmo-add-verify $< io.o $(DEP_FLAGS)

cosmo-add-reads: cosmo-add-reads.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS) -o cosmo-add-reads $< io.o $(DEP_FLAGS)

cosmo-delete-verify: cosmo-delete-verify.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS) -o cosmo-delete-verify $< io.o $(DEP_FLAGS)

cosmo-delete-add: cosmo-delete-add.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS) -o cosmo-delete-add $< io.o $(DEP_FLAGS)

cosmo-validate: cosmo-validate.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS) -o cosmo-validate $< io.o $(DEP_FLAGS)	

all: $(BINARIES)

clean:
		rm -rf $(BINARIES) *.o *.dSYM
