#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <libgen.h> // basename
#include <sys/mman.h> // mlockall

#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_algorithm.hpp>


#include "io.hpp"

//#include "debruijn_graph.hpp"
#include "dynBoss.hpp"

#include "algorithm.hpp"
#include "wt_algorithm.hpp"

#include <stdint.h>
#include "formatutil.cpp"

using namespace std;
using namespace sdsl;

string graph_extension = ".dbg";
string contig_extension = ".fasta";

struct parameters_t {
  std::string input_filename = "";
   std::string fasta_filename = "";
  std::string output_prefix = "";
};

void parse_arguments(int argc, char **argv, parameters_t & params);
void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd("Cosmo Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".dbg file (output from cosmo-build).", true, "", "input_file", cmd);

  string output_short_form = "output_prefix";
  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. Contigs will be written to [" + output_short_form + "]" + contig_extension + ". " +
            "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
  cmd.parse( argc, argv );

  // -d flag for decompression to original kmer biz
  params.input_filename  = input_filename_arg.getValue();
  //params.fasta_filename = fasta_filename_arg.getValue();
  params.output_prefix   = output_prefix_arg.getValue();
}

int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);

  // TO LOAD:
  string fname= p.input_filename + ".solid_kmers_binary.packed";
  cerr << "Constructing dynamic BOSS..." << endl;
  ifstream input(fname, ios::in|ios::binary|ios::ate);
  dyn_boss dbg;
  dbg.load_from_packed_edges( input, "$ACGT" );
  input.close();
  
  cerr << "Constructed graph: " << endl;
  cerr << "k             : " << dbg.k << endl;
  cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg.num_edges() << endl;
  size_t bs = dbg.bit_size();
  cerr << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
  cerr << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;

  // get k-mers and edgemers from file and construct online
  cout << "Constructing graph online..." << endl;
  cout << "Reading k-mer list..." << endl;
  
  unordered_set<kmer_t> kmers;
  unordered_set<kmer_t> edgemers;
  size_t k = dbg.k;
  auto start = std::chrono::system_clock::now();
  string fastaFile= p.input_filename + "_2.fasta";
  //handle_mers( fastaFile, k - 1, kmers, edgemers );
  fname= p.input_filename + "_1.solid_kmers_binary.packed";
  cerr << "Constructing dynamic BOSS..." << endl;
  ifstream input2(fname, ios::in|ios::binary|ios::ate);
  dyn_boss sdbg;
  sdbg.load_from_packed_edges( input2, "$ACGT" );
  input2.close();

  cerr << "Starting online phase...\n";
  size_t nAdded = 0;
  for (size_t i =0; i < dbg.num_edges(); ++i) {
     string mer = dbg.edge_label(i);
     if (mer.find('$') == string::npos) {
	//	cerr << "Adding edge " << mer << endl;
	if (sdbg.index( mer.begin(),0 ) == 0 ){
	   sdbg = Add_Edge(sdbg, dbg.edge_label(i), 1);
	   ++nAdded;
	}
     }
  }

  cout << "Online construction complete: " << nAdded << endl;
  cerr << "k             : " << sdbg.k << endl;
  cerr << "num_nodes()   : " << sdbg.num_nodes() << endl;
  cerr << "num_edges()   : " << sdbg.num_edges() << endl;
  bs = sdbg.bit_size();
  cerr << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
  cerr << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;

  
  cerr << "Verifying edges..." << endl;
  for (size_t i = 0; i < dbg.num_edges(); ++i) {
     //cerr << dbg.edge_label(i) << endl;
     if (dbg.edge_label(i) != sdbg.edge_label(i)) {
	cerr << "Edge verification failed.\n";
	cerr << i << endl;
	cerr << dbg.edge_label(i) << endl;
	cerr << sdbg.edge_label(i) << endl;
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
  
  cerr << "Verification passed!\n";
  
  return 0;

}

