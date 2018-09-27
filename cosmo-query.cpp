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
#include "debruijn_graph.hpp"

#include "debruijn_hypergraph.hpp"
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

  TCLAP::UnlabeledValueArg<std::string> fasta_filename_arg("fasta",
							   ".fasta file (k-mers).",
							   true, "",
							   "fasta_file", cmd);
  string output_short_form = "output_prefix";
  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. Contigs will be written to [" + output_short_form + "]" + contig_extension + ". " +
            "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
  cmd.parse( argc, argv );

  // -d flag for decompression to original kmer biz
  params.input_filename  = input_filename_arg.getValue();
  params.fasta_filename = fasta_filename_arg.getValue();
  params.output_prefix   = output_prefix_arg.getValue();
}

int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);

  // TO LOAD:
  debruijn_graph<> g;
  load_from_file(g, p.input_filename);

  auto dbg_size = size_in_mega_bytes(g);
  cerr << "k             : " << g.k << endl;
  cerr << "num_nodes()   : " << g.num_nodes() << endl;
  cerr << "num_edges()   : " << g.num_edges() << endl;
  cerr << "W size        : " << size_in_mega_bytes(g.m_edges) << " MB" << endl;
  cerr << "L size        : " << size_in_mega_bytes(g.m_node_flags) << " MB" << endl;
  cerr << "DBG size      : " << dbg_size << " MB" << endl;
  cerr << "Bits per edge : " << bits_per_element(g) << " Bits" << endl;

  /*
  cerr << "symbol ends:" << endl;
  for (size_t i = 0; i<5; ++i) {
    cerr << i << ": " << g.m_symbol_ends[i] << endl;
  }
  */

  int num_queries = 2e4;
  size_t min_k = 8; // this could be 0, but it affects longer too much
  size_t max_k = g.k-2; // K is the edge length, and K-1 node lenght.
  // Time will be identical if we measure K-1

  // get k-mers and edgemers from file
  cout << "Reading member k-mer list..." << endl;
  
  unordered_set<kmer_t> kmers;
  unordered_set<kmer_t> edgemers;
  auto start = std::chrono::system_clock::now();
  handle_mers( p.fasta_filename, 27, kmers, edgemers );


  unsigned k = g.k;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<u_int64_t> sample_dis(0, pow(2, 2*k) - 1);
  size_t nQuery = 5e5;
  kmer_t kk;
  double t_elapsed = 0.0;
  clock_t t_start;
  string merToQuery;
  for (unsigned i = 0; i < nQuery; ++i) {
     //generate a random k-mer
     kk = sample_dis( gen ) ;

     // cout << "Testing k-mer: " << get_kmer_str( kk, k ) << ' ' << endl;
     merToQuery = get_kmer_str( kk, k );
     t_start = clock();
     //boopair<size_t,size_t> nn =
     g.index( merToQuery.begin() );

     t_elapsed += double(clock() - t_start);
     //cout << nn.first << " " << nn.second << endl;
   }
   t_elapsed = t_elapsed/ CLOCKS_PER_SEC;

   double avgQueryTime = t_elapsed / nQuery;
   cout << "Average query time: " << avgQueryTime << " s" << endl;

   if (!p.output_prefix.empty()) {
      ofstream ofs( p.output_prefix.c_str(), ofstream::out | ofstream::app );
      ofs << avgQueryTime << ' ';
   }
   
   cout << "Beginning queries of members..." << endl;
   
   std::uniform_int_distribution<u_int64_t> unif_dist (0, kmers.size() - 1);

   vector<kmer_t> vKmers( kmers.begin(), kmers.end() );
   t_elapsed = 0.0;
   bool bMem = false;
   for (unsigned i = 0; i < nQuery; i++) {
      int randnum = unif_dist( gen );

      // Test membership
      merToQuery = get_kmer_str( vKmers[ randnum ], k );
      clock_t t_start = clock();
      bMem = g.index( merToQuery.begin() );
      t_elapsed += double (clock() - t_start);

      assert( bMem );
   }

   avgQueryTime = t_elapsed / nQuery / CLOCKS_PER_SEC;
   cout << "Average query time: " << avgQueryTime << " s" << endl;

   if (!p.output_prefix.empty()) {
      ofstream ofs( p.output_prefix.c_str(), ofstream::out | ofstream::app );
      ofs << avgQueryTime << endl;
   }
}

