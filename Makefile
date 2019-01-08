# NOTE: needs boost, tclap, and sdsl
#LD_LIBRARY_PATH="/ufrc/boucher/baharpan/VARI/3rd_party_inst/boost/lib/"
CXX=g++ # clang++-3.8 # g+
CPP_FLAGS=-m64 -std=c++0x -Wall 
DEP_PATH=../VARI/3rd_party_inst
BOOST_PATH=$(DEP_PATH)/boost
INC_PATH=-isystem $(DEP_PATH)/include/ -isystem $(BOOST_PATH)/include 
INC_PATH_DYN=/ufrc/boucher/baharpan/DYN-BOSS/dyn/include
LIB_PATH=-L$(DEP_PATH)/lib  -L./ -L$(BOOST_PATH)/lib 
BOOST_FLAGS= -lboost_system -lboost_filesystem
DEP_FLAGS=$(INC_PATH) -I$(INC_PATH_DYN) $(LIB_PATH)  $(BOOST_FLAGS)  -lsdsl -fopenmp
DEBUG_FLAGS=-g
NDEBUG_FLAGS=-DNDEBUG
OPT_FLAGS=-O3 -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native
NOPT_FLAGS=-Og
DYN_FLAGS=-DDYNAMIC

# Using Semantic Versioning: http://semver.org/
VERSION=0.5.1
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
BINARIES= cosmo-pack cosmo-build cosmo-build-dyn cosmo-add-delete cosmo-benchmark cosmo-benchmark-dyn # cosmo-assemble

default: all

lut.hpp: make_lut.py
		python make_lut.py > lut.hpp

io.o: io.hpp io.cpp debug.h dummies.hpp kmer.hpp
		$(CXX) $(CPP_FLAGS) -c io.cpp $(DEP_FLAGS)

# TODO: Roll these all into one... "cosmo"
cosmo-pack: cosmo-pack.cpp $(PACK_REQS)
		$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(DEP_FLAGS)

cosmo-build: cosmo-build.cpp $(BUILD_REQS)
	$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(DEP_FLAGS) 

cosmo-build-dyn: cosmo-build.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS) $(DYN_FLAGS) -o cosmo-build-dyn $< io.o $(DEP_FLAGS)

cosmo-add-delete: cosmo-add-delete.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS) -o cosmo-add-delete $< io.o $(DEP_FLAGS)

cosmo-verify: cosmo-verify.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS) -o cosmo-verify $< io.o $(DEP_FLAGS)

cosmo-assemble: cosmo-assemble.cpp $(ASSEM_REQS) wt_algorithm.hpp debruijn_hypergraph.hpp
		$(CXX) $(CPP_FLAGS) -o $@ $< $(DEP_FLAGS) 

cosmo-benchmark: cosmo-benchmark.cpp $(ASSEM_REQS) wt_algorithm.hpp debruijn_hypergraph.hpp
		$(CXX) $(CPP_FLAGS) -o $@ $< $(DEP_FLAGS)

cosmo-benchmark-dyn: cosmo-benchmark.cpp $(ASSEM_REQS) wt_algorithm.hpp debruijn_hypergraph.hpp
		$(CXX) $(CPP_FLAGS) $(DYN_FLAGS) -o cosmo-benchmark-dyn $< $(DEP_FLAGS)
cosmo-query: cosmo-query.cpp $(ASSEM_REQS) wt_algorithm.hpp debruijn_hypergraph.hpp
		$(CXX) $(CPP_FLAGS) -o $@ $< $(DEP_FLAGS) 

all: $(BINARIES)

clean:
		rm -rf $(BINARIES) *.o *.dSYM
