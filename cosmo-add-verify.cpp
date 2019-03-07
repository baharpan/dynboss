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
   TCLAP::CmdLine cmd("Dynamic BOSS. Copyright (c) Bahar Alipanahi, Alex Bowe, Alan Kuhnle 2019", ' ', VERSION);
   TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".packed edge file (output from pack-edges).", true, "", "input_file", cmd);
   TCLAP::UnlabeledValueArg<std::string> kmer_filename_arg("kmers",
            "kmers to add in plain text. One per line.", true, "", "input_file", cmd);
   string output_short_form = "output_prefix";
   TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
						  "Output prefix. Graph will be written to [" + output_short_form + "]" + extension + ". " +
						  "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
   cmd.parse( argc, argv );

   params.input_filename  = input_filename_arg.getValue();
   params.kmer_filename  = kmer_filename_arg.getValue();
   params.output_prefix   = output_prefix_arg.getValue();
}

void getKmers( size_t& nKmers,
	       size_t k,
	       vector< string >& kmers,
	       parameters_t& p) {
   cerr << "Getting k-mers from " << p.kmer_filename << "...\n";
   ifstream in(p.kmer_filename );
   string sline;
   vector< string > vline;
   while ( getline( in, sline ) ) {
      vline.push_back( sline );
   }

   size_t pos = 0;
   vector< string > reads;
   string read;
   do {
      if (vline[pos][0] == '>') {
	 //finish current read and start a new one
	 if (!read.empty()) {
	    reads.push_back(read);
	    read.clear();
	 }
      } else {
	 read += vline[pos];
      }
      
      ++pos;
   } while (pos != vline.size());

   if (!read.empty()) //handle the last read
      reads.push_back( read );
   
   for (size_t i = 0; i < reads.size(); ++i) {
      string sline = reads[i];
      size_t read_length = sline.size();
	   
      size_t nMers = read_length - k + 1;
      for (size_t start = 0; start < nMers; ++start) {
	 string kmer = sline.substr( start, k );
	 kmers.push_back( kmer );
      }
   }

   in.close();
   nKmers = kmers.size();
}

int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);

  dyn_boss sdbg;
  
  cout << "Constructing dynamic BOSS from DSK input..." << endl;
  ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);

  sdbg.load_from_packed_edges( input, "$ACGT" );
  input.close();
  cout << "graph info: " << endl;
  cout << "k             : " << sdbg.k << endl;
  cout << "num_nodes()   : " << sdbg.num_nodes() << endl;
  cout << "num_edges()   : " << sdbg.num_edges() << endl;
  // size_t bs = sdbg.bit_size();
  // cout << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
  // cout << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;
  
  cout << "===============================\n";
  cout << "Building dbg online...\n";

  dyn_boss dbg( sdbg.k );
  cout << "Reading k-mers from FASTA file..." << endl;
  size_t nKmers;
  vector<string> kmer_2;
  
  getKmers( nKmers, dbg.k, kmer_2, p );

  cout << nKmers << " read.";

  clock_t t_start = clock();
  for (size_t i = 0; i < kmer_2.size();i++) {
      if ( i % (kmer_2.size() / 100) == 0 ) {
	 cerr << "\r                                            \r";
	 cerr << ((double)i) / kmer_2.size() * 100.0 << "%";
      }
     
      dbg.add_edge(kmer_2[i]);
  }

  double t_elapsed = (clock() - t_start) / CLOCKS_PER_SEC;
  cout << "\nDONE with addition of all kmers\n";
  cout << "Time per Op: " << t_elapsed / kmer_2.size() << endl;

  t_elapsed = (clock() - t_start) / CLOCKS_PER_SEC;

  cout << "===============================\n";
  cout << "new greph     : " << endl;
  cout << "k             : " << dbg.k << endl;
  cout << "num_nodes()   : " << dbg.num_nodes() << endl;
  cout << "num_edges()   : " << dbg.num_edges() << endl;
  // bs = dbg.bit_size();
  // cout << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
  // cout << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;

  cout << "===============================\n";
  cout << "Beginning graph validation...\n";
  
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
