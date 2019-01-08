#include <iostream>
#include <fstream>
 
#include <libgen.h> // basename

#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "io.hpp"
#ifdef DYNAMIC
#include "dynBoss.hpp"
#else
#include "debruijn_graph.hpp"
#endif
#include "algorithm.hpp"

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
#ifdef DYNAMIC
cout<<"jhe"<<endl;
  dyn_boss dbg;
  dbg.load_from_packed_edges( input, "$ACGT" );
#else
cout<<"aaa"<<endl;
    debruijn_graph<> dbg = debruijn_graph<>::load_from_packed_edges(input, "$ACGT"/*, &minus_positions*/);
    cout<<"aaa"<<endl;
#endif
  input.close();

  cerr << "k             : " << dbg.k << endl;
  cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg.num_edges() << endl;
  size_t bs = dbg.bit_size();
  cerr << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
  cerr << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;

  // cerr << "Edges: " << endl;
  // for (size_t i = 0; i < dbg.num_edges(); ++i) {
  //    cerr << dbg.edge_label(i) << endl;
  // }
  // cerr << "Nodes: " << endl;
  // for (size_t i = 0; i < dbg.num_nodes(); ++i) {
  //    cerr << dbg.node_label( i ) << endl;
  // }

  // for (size_t i = 0; i < dbg.num_nodes(); ++i) {
  //    cerr << "outdegree(" << i << "): ";
  //    cerr << dbg.outdegree( i ) << endl;
  // }

  // for (size_t i = 0; i < dbg.num_nodes(); ++i) {
  //    cerr << "indegree(" << i << "): ";
  //    cerr << dbg.indegree( i ) << endl;
  // }

  // for (size_t i = 0; i < dbg.num_nodes(); ++i) {
  //    for (uint64_t j = 0; j < 5; ++j) {
  // 	cerr << i << ".outgoing(" << j << "): " << dbg.outgoing(i, j) << endl;
  //    }
  // }

  // for (size_t i = 0; i < dbg.num_nodes(); ++i) {
  //    for (uint64_t j = 0; j < 5; ++j) {
  // 	cerr << i << ".incoming(" << j << "): " << dbg.incoming(i, j) << endl;
  //    }
  // }
  
  
  char * base_name = basename(const_cast<char*>(p.input_filename.c_str()));
  string outfilename = ((p.output_prefix == "")? base_name : p.output_prefix) + extension;
  //store_to_file(dbg, outfilename);
  ofstream ofs( outfilename.c_str(), ios::out | ios::binary );
  dbg.serialize( ofs );
  ofs.close();
}
