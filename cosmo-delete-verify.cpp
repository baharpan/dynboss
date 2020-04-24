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
#include "kmer-counter.hpp"
using namespace std;
using namespace sdsl;

/*This code adds k-mers from a FASTA file into the DynamicBOSS with add-edge
then deletes them with delete-edge and confirms that the starting and
final graphs are the same.
usage: ./cosmo-delete-verify <.dbg file> <.fasta file>*/



struct parameters_t {
   std::string input_filename = "";
   std::string kmer_filename = "";
};

void parse_arguments(int argc, char **argv, parameters_t & params);
void parse_arguments(int argc, char **argv, parameters_t & params)
{
   TCLAP::CmdLine cmd("DynamicBOSS. Copyright (c) Bahar Alipanahi, Alan Kuhnle, Alex Bowe 2019", ' ', VERSION);
   TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".dbg file (output from cosmo-build-dyn).", true, "", "dbg file", cmd);
   TCLAP::UnlabeledValueArg<std::string> kmer_filename_arg("kmers",
            ".fasta file to count kmers to add and delete.", true, "", "fasta file", cmd);

   cmd.parse( argc, argv );

   params.input_filename  = input_filename_arg.getValue();
   params.kmer_filename  = kmer_filename_arg.getValue();
}


int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);

  cerr<<"Loading DynamicBOSS from file:  "<<p.input_filename<<endl;
  dyn_boss sdbg;
  ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
  load_from_file(sdbg, p.input_filename);
  //cout << "Constructing dynamic BOSS from DSK input..." << endl;
  //sdbg.load_from_packed_edges( input, "$ACGT" );
  input.close();
  cout << "graph info: " << endl;
  cout << "k             : " << sdbg.k << endl;
  cout << "num_nodes()   : " << sdbg.num_nodes() << endl;
  cout << "num_edges()   : " << sdbg.num_edges() << endl;
  size_t bs = sdbg.bit_size();
  cout << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
  cout << "Bits per edge : " << bs / static_cast<double>(sdbg.num_edges()) << " Bits" << endl;
  //cout << "matrix:\n";
  //  sdbg.print_boss_matrix( cerr );
  cout << "===============================\n";

  dyn_boss dbg = sdbg;
  cout << "Reading FASTA file " <<p.kmer_filename<< endl;
  size_t nKmers;
  set<string> kmer_2;
  vector<string> kmers_added;

  getKmers( nKmers, dbg.k, kmer_2, p.kmer_filename );
  cout << nKmers << " distinct kmers were counted " << endl;

  cout<<"finding the kmers that are not already in the graph ..."<<endl;
  for (auto it = kmer_2.begin(); it != kmer_2.end(); ++it) {
      string kmer = *it;
      if (!dbg.index_edge_alan( kmer.begin() )) {
         kmers_added.push_back( kmer );
      }
  }

  cout << "===============================\n";
  cout << "Start adding kmers ... "<<endl;
  clock_t t_start = clock();

  for (size_t i = 0; i < kmers_added.size();i++) {
       dbg.add_edge(kmers_added[i]);
   }

  double t_elapsed = (clock() - t_start) / CLOCKS_PER_SEC;
  cout << "Time per Additions (s): " << t_elapsed / kmers_added.size() << endl;
  cout << "DONE with addition of " << kmers_added.size() << " kmers.\n";
  cout << "The remainder (if any) were already in the graph." << endl;


  cout << "===============================\n";
  cout << "Start deleting the added k-mers ..." << endl;

  t_start = clock();

  for (size_t i = 0; i < kmers_added.size();i++) {
     dbg.delete_edge(kmers_added[i]);
  }

  t_elapsed = (clock() - t_start) / CLOCKS_PER_SEC;
  cout << "Time per Deletion (s): " << t_elapsed / kmers_added.size() << endl;


  cout << "===============================\n";
  cout << "final graph   : " << endl;
  cout << "k             : " << dbg.k << endl;
  cout << "num_nodes()   : " << dbg.num_nodes() << endl;
  cout << "num_edges()   : " << dbg.num_edges() << endl;
  bs = dbg.bit_size();
  cout << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
  cout << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;

  cout << "===============================\n";
  cout << "Beginning of graph validation...\n";

  cout << "Verifying edges..." << endl;
  for (size_t i = 0; i < dbg.num_edges(); i+= 1) {
     if (dbg.edge_label(i) != sdbg.edge_label(i)) {
	cout << "Edge verification failed.\n";
	size_t j = i;
	while (dbg.edge_label(j)[0] == '$') {
	   j = dbg._forward( j );
	}

	cout << "The edge at index\n"
	     << j << "\nwith indegree\n"
	     << dbg.indegree( j )
	     << "\n and label\n"
	     << dbg.edge_label(j)
	     << "\n has an incoming dummy string.\n";
	cout << "backward: " << dbg.edge_label(dbg._backward( j )) << endl;

	size_t idxj;
	if (sdbg.index_edge_alan( dbg.edge_label(j).begin(), idxj )) {
	   cout << "This edge is at " << idxj << " in static graph.\n";
	   cout << "It has label\n"
		<< sdbg.edge_label( idxj ) << '\n';
	   cout << "backward: " << sdbg.edge_label(sdbg._backward( idxj )) << endl;
	} else {
	   cout << "This edge is not present in sdbg.\n";
	}
	exit(0);
     }
  }

  cout << "Verifying nodes..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); i+=1) {
     if (dbg.node_label(i) != sdbg.node_label(i)) {
	cout << "Node verification failed.\n";
	exit(0);
     }
  }

  cout << "Verifying outdegree..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); i+= 1) {
     if ( dbg.outdegree( i ) != sdbg.outdegree(i) ) {
	cout << "outdegree verification failed.\n";
	exit(0);
     }
  }

  cout << "Verifying indegree..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); i+= 1) {
     if (dbg.indegree( i ) != sdbg.indegree(i)) {
	cout << "indegree verification failed.\n";
	exit(0);
     }
  }

  cout << "Verifying incoming..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); i+= 1) {
     for (uint64_t j = 0; j < 5; ++j) {
  	//cout << i << ".incoming(" << j << "): " << dbg.incoming(i, j) << endl;
	if (dbg.incoming( i,j ) != sdbg.incoming(i,j)) {
	   cout << "incoming verification failed.\n";
	   exit(0);
	}

     }
  }

  cout << "Verifying outgoing..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); ++i) {
     for (uint64_t j = 0; j < 5; ++j) {
  	//cout << i << ".incoming(" << j << "): " << dbg.incoming(i, j) << endl;
	if (dbg.outgoing( i,j ) != sdbg.outgoing(i,j)) {
	   cout << "outgoing verification failed.\n";
	   exit(0);
	}

     }
  }

  cout << "Verification passed!\n";

  return 0;
}
