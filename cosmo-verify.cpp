#include <iostream>
#include <fstream>

#include <libgen.h> // basename

#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "io.hpp"
#include "dynBoss.hpp"
#include "debruijn_graph.hpp"
#include "algorithm.hpp"
#include "formatutil.cpp"

using namespace std;
using namespace sdsl;

string extension = ".dbg";

struct parameters_t {
  std::string input_filename = "";
  std::string output_prefix = "";
};

void parse_arguments(int argc, char **argv, parameters_t & params);
void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd("Cosmo Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".packed edge file (output from pack-edges).", true, "", "input_file", cmd);
  string output_short_form = "output_prefix";
  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. Graph will be written to [" + output_short_form + "]" + extension + ". " +
            "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
  params.output_prefix   = output_prefix_arg.getValue();
}

int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);

  ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);

  cerr << "Loading dynamic dBG...\n";
  dyn_boss dbg;
  dbg.load_from_packed_edges( input, "$ACGT" );
  input.close();
  input.clear();
  input.open( p.input_filename, ios::in|ios::binary|ios::ate);
  cerr << "Loading static dBG...\n";
  debruijn_graph<> sdbg = debruijn_graph<>::load_from_packed_edges(input, "$ACGT"/*, &minus_positions*/);

  input.close();

  cerr << "k             : " << dbg.k << ' ' << sdbg.k << endl;
  cerr << "num_nodes()   : " << dbg.num_nodes() << ' ' << sdbg.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg.num_edges() << ' ' << sdbg.num_edges() << endl;
  size_t bs = dbg.bit_size();
  size_t sbs = sdbg.bit_size();
  cerr << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB, "
       << sbs / 8.0 / 1024.0 / 1024.0 << endl;
  cerr << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits, "
       << sbs / static_cast<double>(dbg.num_edges()) << " Bits." << endl;

  cerr << "Verifying edges..." << endl;
  for (size_t i = 0; i < dbg.num_edges(); ++i) {
     //cerr << dbg.edge_label(i) << endl;
     if (dbg.edge_label(i) != sdbg.edge_label(i)) {
	cerr << "Edge verification failed.\n";
	exit(0);
     }
  }

  cerr << "Verifying nodes..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); ++i) {
     //cerr << dbg.edge_label(i) << endl;
     if (dbg.node_label(i) != sdbg.node_label(i)) {
	cerr << "Node verification failed.\n";
	exit(0);
     }
  }

  cerr << "Verifying outdegree..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); ++i) {
     if ( dbg.outdegree( i ) != sdbg.outdegree(i) ) {
	cerr << "outdegree verification failed.\n";
	exit(0);
     }
  }

  cerr << "Verifying indegree..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); ++i) {
     if (dbg.indegree( i ) != sdbg.indegree(i)) {
	cerr << "indegree verification failed.\n";
	exit(0);
     }
  }

  cerr << "Verifying incoming..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); ++i) {
     for (uint64_t j = 0; j < 5; ++j) {
  	//cerr << i << ".incoming(" << j << "): " << dbg.incoming(i, j) << endl;
	if (dbg.incoming( i,j ) != sdbg.incoming(i,j)) {
	   cerr << "incoming verification failed.\n";
	   exit(0);
	}

     }
  }

  cerr << "Verifying outgoing..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); ++i) {
     for (uint64_t j = 0; j < 5; ++j) {
  	//cerr << i << ".incoming(" << j << "): " << dbg.incoming(i, j) << endl;
	if (dbg.outgoing( i,j ) != sdbg.outgoing(i,j)) {
	   cerr << "outgoing verification failed.\n";
	   exit(0);
	}
      

     }
  }

  cerr << "Verifying index...\n";

  unordered_set<kmer_t> kmers;
  unordered_set<kmer_t> edgemers;
  handle_mers( "temp.fasta", dbg.k - 1, kmers, edgemers );
  vector<kmer_t> vKmers( kmers.begin(), kmers.end() );

  std::random_device rd;
  std::mt19937 gen(rd());
  size_t k = dbg.k - 1;
  std::uniform_int_distribution<u_int64_t> sample_dis(0, pow(2, 2*k) - 1);
  size_t nQuery = 5e1;
  kmer_t kk;

  string merToQuery;
  std::uniform_int_distribution<u_int64_t> unif_dist (0, kmers.size() - 1);

  for (unsigned i = 0; i < nQuery; ++i) {
     //generate a random k-mer
     kk = sample_dis( gen ) ;

     merToQuery = get_kmer_str( kk, k );
     //cerr << "Testing " << merToQuery << endl;
     if ( dbg.index( merToQuery.begin() ) != sdbg.index( merToQuery.begin() ) ) {
	cerr << "Index failed.\n";
	exit(0);
     }

     int randnum = unif_dist( gen );
     merToQuery = get_kmer_str( vKmers[ randnum ], k );
     //cerr << "Testing " << merToQuery << endl;
     if ( dbg.index( merToQuery.begin() ) != sdbg.index( merToQuery.begin() ) ) {
	cerr << "Index failed.\n";
	exit(0);
     }

  }

  cerr << "Verification passed!\n";

  // char * base_name = basename(const_cast<char*>(p.input_filename.c_str()));
  // string outfilename = ((p.output_prefix == "")? base_name : p.output_prefix) + extension;
  // //store_to_file(dbg, outfilename);
  // ofstream ofs( outfilename.c_str(), ios::out | ios::binary );
  // dbg.serialize( ofs );
  // ofs.close();
}
