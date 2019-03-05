#include <iostream>
#include <fstream>

#include <libgen.h> // basename

#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>


#include "io.hpp"

#include "dynBoss.hpp"

#include "algorithm.hpp"
#include "utility.hpp"
#include "formatutil.cpp"
using namespace std;
using namespace sdsl;

static char base[] = {'$','A','C','G','T'};

string extension = ".dbg";

struct parameters_t {
   std::string input_filename = "";
   std::string kmer_filename = "";
   std::string input2_filename = "";
   std::string output_prefix = "";
};
void parse_arguments(int argc, char **argv, parameters_t & params);
void parse_arguments(int argc, char **argv, parameters_t & params)
{
   TCLAP::CmdLine cmd("Cosmo Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".packed edge file (output from pack-edges).", true, "", "input_file", cmd);
  TCLAP::UnlabeledValueArg<std::string> kmer_filename_arg("kmers",
            "kmers to add in plain text. One per line.", true, "", "input_file", cmd);
  TCLAP::UnlabeledValueArg<std::string> kmer_filename_arg("input2",
            ".packed edge file (output from pack-edges).", true, "", "input_file", cmd);
  string output_short_form = "output_prefix";
  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. Graph will be written to [" + output_short_form + "]" + extension + ". " +
            "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
  params.kmer_filename  = kmer_filename_arg.getValue();
  params.output_prefix   = output_prefix_arg.getValue();
}



int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);

  cerr << "Constructing dynamic BOSS..." << endl;
  ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
  ifstream inputKmers(p.kmer_filename, ios::in|ios::binary|ios::ate);

  dyn_boss dbg;
  dbg.load_from_packed_edges( input, "$ACGT" );
  input.close();
  cerr << "starting graph info: " << endl;
  cerr << "k             : " << dbg.k << endl;
  cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
    cerr << "num_edges()   : " << dbg.num_edges() << endl;
    size_t bs = dbg.bit_size();
    cerr << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
    cerr << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;
    dyn_boss sdbg = dbg;


    static const char alphanum[] = "ACGT";

    vector<string>all_kmers;

    int stringLength = sizeof(alphanum) - 1;

    size_t nOps = 100;
    while (all_kmers.size() < nOps ){
        string kmer;
        while(kmer.size() < dbg.k )
	   kmer += alphanum[rand() % stringLength];
	all_kmers.push_back(kmer);
    }

    cerr<< "Number of kmers to add: " << all_kmers.size() << endl;
    clock_t t_start = clock();
    for (size_t i = 0; i< all_kmers.size();i++){
       dbg = Add_Edge (dbg, all_kmers[i] ,1);
      /*if (dbg.index(all_kmers[i].begin(),0) == 0){
        cerr<<"failed to add kmer "<<i<<endl;
        exit(0);
      }*/
    }
    double t_elapsed = (clock() - t_start) / CLOCKS_PER_SEC;
    cerr << "DONE with addition of all kmers\n";
    cerr << "Time per Op: " << t_elapsed / nOps << endl;

    t_elapsed = (clock() - t_start) / CLOCKS_PER_SEC;

    cerr << "===============================\n";
    cerr << "new greph     : " << endl;
    cerr << "k             : " << dbg.k << endl;
    cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
    cerr << "num_edges()   : " << dbg.num_edges() << endl;
    bs = dbg.bit_size();
    cerr << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
    cerr << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;

  /*   for (size_t i = 0; i<dbg.num_edges(); i++)
         cout<<dbg.edge_label(i)<< " " <<i<<endl;

     for (size_t i=0; i < dbg.num_nodes();i++)
         cout<<dbg.node_label(i)<< " " <<i<<endl;*/
    size_t nVerify = 1000;
    cerr << "Verifying edges..." << endl;
  for (size_t i = 0; i < dbg.num_edges(); i+=dbg.num_edges()/nVerify) {
     if (dbg.edge_label(i) != sdbg.edge_label(i)) {
	cerr << "Edge verification failed.\n";
	exit(0);
     }
  }

  cerr << "Verifying nodes..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); i+=dbg.num_nodes()/nVerify) {
     if (dbg.node_label(i) != sdbg.node_label(i)) {
	cerr << "Node verification failed.\n";
	exit(0);
     }
  }

  cerr << "Verifying outdegree..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); i+=dbg.num_nodes()/nVerify) {
     if ( dbg.outdegree( i ) != sdbg.outdegree(i) ) {
	cerr << "outdegree verification failed.\n";
	exit(0);
     }
  }

  cerr << "Verifying indegree..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); i+=dbg.num_nodes()/nVerify) {
     if (dbg.indegree( i ) != sdbg.indegree(i)) {
	cerr << "indegree verification failed.\n";
	exit(0);
     }
  }

  cerr << "Verifying incoming..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); i+=dbg.num_nodes()/nVerify) {
     for (uint64_t j = 0; j < 5; ++j) {
  	//cerr << i << ".incoming(" << j << "): " << dbg.incoming(i, j) << endl;
	if (dbg.incoming( i,j ) != sdbg.incoming(i,j)) {
	   cerr << "incoming verification failed.\n";
	   exit(0);
	}

     }
  }

  cerr << "Verifying outgoing..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); i+=dbg.num_nodes()/5) {
     for (uint64_t j = 0; j < 5; ++j) {
  	//cerr << i << ".incoming(" << j << "): " << dbg.incoming(i, j) << endl;
	if (dbg.outgoing( i,j ) != sdbg.outgoing(i,j)) {
	   cerr << "outgoing verification failed.\n";
	   exit(0);
	}

     }
  }

cerr << "Verification passed!\n";
cerr << "===============================\n";
cerr << "Deleting some random nodes\n";
for (size_t i = 0; i < dbg.num_nodes(); i += dbg.num_nodes()/3)
  dbg = Delete_node(dbg , dbg.node_label(i));
  cerr << "new greph     : " << endl;
  cerr << "k             : " << dbg.k << endl;
  cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg.num_edges() << endl;
  bs = dbg.bit_size();
  cerr << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
  cerr << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;
}
