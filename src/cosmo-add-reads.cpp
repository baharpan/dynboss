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
	       vector< vector< string > >& readKmers,
	       parameters_t& p) {
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
      readKmers.push_back( vector< string >() );
      size_t nMers = read_length - k + 1;
      for (size_t start = 0; start < nMers; ++start) {
	 string kmer = sline.substr( start, k );
	 kmers.push_back( kmer );
	 (* (readKmers.end() - 1)).push_back( kmer );
      }
   }

   in.close();
   nKmers = kmers.size();
}

int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);

  {
     dyn_boss dbg;
  
     cout << "Constructing dynamic BOSS from DSK input..." << endl;
     ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);

     dbg.load_from_packed_edges( input, "$ACGT" );
     input.close();
     cout << "graph info: " << endl;
     cout << "k             : " << dbg.k << endl;
     cout << "num_nodes()   : " << dbg.num_nodes() << endl;
     cout << "num_edges()   : " << dbg.num_edges() << endl;
     size_t bs = dbg.bit_size();
     cout << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
     cout << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;
  
     cout << "===============================\n";
     cout << "Reading k-mers from FASTA file..." << endl;
     size_t nKmers;
     vector<string> kmer_2;
     vector< vector<string> > readkmer;
  
     getKmers( nKmers, dbg.k, kmer_2, readkmer, p );

     cerr << nKmers << " read.\n";
     cerr << readkmer.size() << " reads.\n";

     clock_t t_start = clock();
     for (size_t i = 0; i < readkmer.size();i++) {
	// if ( i % (readkmer.size() / 100) == 0 ) {
	// 	 cerr << "\r                                            \r";
	// 	 cerr << ((double)i) / readkmer.size() * 100.0 << "%";
	//  }
     
	dbg.add_read(readkmer[i]);
     }

     double t_elapsed = (clock() - t_start) / CLOCKS_PER_SEC;
     cout << "\nDONE with addition of all reads\n";
     cout << "Time per Read: " << t_elapsed / readkmer.size() << endl;
     cout << "Time per OP: " << t_elapsed / kmer_2.size() << endl;

     t_elapsed = (clock() - t_start) / CLOCKS_PER_SEC;

     cout << "===============================\n";
     cout << "new greph     : " << endl;
     cout << "k             : " << dbg.k << endl;
     cout << "num_nodes()   : " << dbg.num_nodes() << endl;
     cout << "num_edges()   : " << dbg.num_edges() << endl;
     bs = dbg.bit_size();
     cout << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
     cout << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;
  }

  {
     dyn_boss dbg;
  
     cout << "Constructing dynamic BOSS from DSK input..." << endl;
     ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);

     dbg.load_from_packed_edges( input, "$ACGT" );
     input.close();
     cout << "graph info: " << endl;
     cout << "k             : " << dbg.k << endl;
     cout << "num_nodes()   : " << dbg.num_nodes() << endl;
     cout << "num_edges()   : " << dbg.num_edges() << endl;
     size_t bs = dbg.bit_size();
     cout << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
     cout << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;
  
     cout << "===============================\n";
     cout << "Reading k-mers from FASTA file..." << endl;
     size_t nKmers;
     vector<string> kmer_2;
     vector< vector<string> > readkmer;
  
     getKmers( nKmers, dbg.k, kmer_2, readkmer, p );

     cerr << nKmers << " read.\n";
     cerr << readkmer.size() << " reads.\n";

     clock_t t_start = clock();
     for (size_t i = 0; i < kmer_2.size();i++) {
	// if ( i % (readkmer.size() / 100) == 0 ) {
	// 	 cerr << "\r                                            \r";
	// 	 cerr << ((double)i) / readkmer.size() * 100.0 << "%";
	//  }
     
	dbg.add_edge(kmer_2[i], false);
     }

     double t_elapsed = (clock() - t_start) / CLOCKS_PER_SEC;
     cout << "\nDONE with addition of all reads\n";
     cout << "Time per Read: " << t_elapsed / readkmer.size() << endl;
     cout << "Time per OP: " << t_elapsed / kmer_2.size() << endl;

     t_elapsed = (clock() - t_start) / CLOCKS_PER_SEC;

     cout << "===============================\n";
     cout << "new greph     : " << endl;
     cout << "k             : " << dbg.k << endl;
     cout << "num_nodes()   : " << dbg.num_nodes() << endl;
     cout << "num_edges()   : " << dbg.num_edges() << endl;
     bs = dbg.bit_size();
     cout << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
     cout << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;

  }


  return 0;
}
