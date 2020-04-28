# NOTE: needs boost, tclap, and sdsl
#LD_LIBRARY_PATH="/ufrc/boucher/baharpan/VARI/3rd_party_inst/boost/lib/"
CXX=g++ # clang++-3.8 # g+
CPP_FLAGS=-m64 -std=c++1z -Wall
DEP_PATH=3rd_party_inst
BOOST_PATH=$(DEP_PATH)/boost
INC_PATH=-isystem $(DEP_PATH)/include/ -isystem $(BOOST_PATH)/include
INC_PATH_DYN=../DYNAMIC/include -I ../hopscotch-map/include
LIB_PATH=-L$(DEP_PATH)/lib  -L./ -L$(BOOST_PATH)/lib
BOOST_FLAGS= -lboost_system -lboost_filesystem
DEP_FLAGS=$(INC_PATH) -I$(INC_PATH_DYN) $(LIB_PATH)  $(BOOST_FLAGS)  -lsdsl -fopenmp
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
DYN_BUILD_REQS=dynBoss.hpp io.hpp io.o debug.h kmer-counter.hpp
ASSEM_REQS=debruijn_graph.hpp algorithm.hpp utility.hpp kmer.hpp uint128_t.hpp
PACK_REQS=lut.hpp debug.h io.hpp io.o sort.hpp kmer.hpp dummies.hpp
BINARIES=cosmo-pack dynamicBOSS load-verify cosmo-delete-verify cosmo-delete-add cosmo-query cosmo-add-verify
BINFile=../bin/

default: all

lut.hpp: make_lut.py
		python make_lut.py > lut.hpp

io.o: io.hpp io.cpp debug.h dummies.hpp kmer.hpp
		$(CXX) $(CPP_FLAGS) -c io.cpp $(DEP_FLAGS)

# TODO: Roll these all into one... "cosmo"
cosmo-pack: cosmo-pack.cpp $(PACK_REQS)
		@ mkdir -p ../bin
		$(CXX) $(CPP_FLAGS) -o $(BINFile)$@ $< io.o $(DEP_FLAGS)

dynamicBOSS: dynamicBOSS.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS)  -o $(BINFile)$@ $< io.o $(DEP_FLAGS)

cosmo-query: cosmo-query.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS)  -o $(BINFile)$@ $< io.o $(DEP_FLAGS)

cosmo-add-verify: cosmo-add-verify.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS) -o $(BINFile)$@ $< io.o $(DEP_FLAGS)

load-verify: load-verify.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS)  -o $(BINFile)$@ $< io.o $(DEP_FLAGS)

cosmo-build-dyn: cosmo-build-dyn.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS) $(DYN_FLAGS) -o $(BINFile)$@ $< io.o $(DEP_FLAGS)

cosmo-add-delete: cosmo-add-delete.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS) -o $(BINFile)$@ $< io.o $(DEP_FLAGS)

cosmo-delete-add: cosmo-delete-add.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS) -o $(BINFile)$@ $< io.o $(DEP_FLAGS)

cosmo-verify: cosmo-verify.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS) -o $(BINFile)$@ $< io.o $(DEP_FLAGS)

cosmo-delete-verify: cosmo-delete-verify.cpp $(DYN_BUILD_REQS)
	$(CXX) $(CPP_FLAGS) -o $(BINFile)$@ $< io.o $(DEP_FLAGS)



all: $(BINARIES)

clean:
		rm -rf $(BINARIES) *.o *.dSYM
