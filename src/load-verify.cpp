#include <iostream>
#include <fstream>
#include <libgen.h> // basename
#include "tclap/CmdLine.h"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include "io.hpp"
#include "dynBoss.hpp"
#include "algorithm.hpp"
using namespace std;
using namespace sdsl;

/*This code verifies if the load function in "dynBoss.hpp" works correctly.
Once it builds the DynamicBOSS using .packed file, then it loads the DynamicBOSS from
.dbg file and compares them.
usage: ./load-verify <.packed file> <.dbg file>*/

string extension = ".sdbg";

struct parameters_t {
  std::string input_filename = "";
  std::string graph_filename = "";
  std::string output_prefix = "";
};

void parse_arguments(int argc, char **argv, parameters_t & params);
void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd("DynamicBOSS. Copyright (c) Bahar Alipanahi, Alan Kuhnle, Alex Bowe 2019", ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".packed edge file (output from cosmo-pack).", true, "", "input_file", cmd);
  TCLAP::UnlabeledValueArg<std::string> graph_filename_arg("graph",
            ".dbg file (output from cosmo-build-dyn).", true, "", "graph_file", cmd);
  string output_short_form = "output_prefix";
  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. Graph will be written to [" + output_short_form + "]" + extension + ". " +
            "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
  params.graph_filename  = graph_filename_arg.getValue();
  params.output_prefix   = output_prefix_arg.getValue();
}

int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);

  ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
  ifstream input_graph(p.graph_filename, ios::in|ios::binary|ios::ate);

  cout<<"Dynamic library loaded"<<endl;
  dyn_boss dbg;
  dbg.load_from_packed_edges( input, "$ACGT" );
  input.close();

  cerr << "k             : " << dbg.k << endl;
  cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg.num_edges() << endl;
  size_t bs = dbg.bit_size();
  cerr << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
  cerr << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;

  char * base_name = basename(const_cast<char*>(p.input_filename.c_str()));
  string outfilename = ((p.output_prefix == "")? base_name : p.output_prefix) + extension;
  cerr<<"Writing DynamicBOSS in file:  "<<outfilename.c_str()<<endl;
  ofstream ofs( outfilename.c_str(), ios::out | ios::binary );
  dbg.serialize( ofs );
  ofs.close();




  dyn_boss dbg2;
  cerr<<"Loading DynamicBOSS from file:  "<<p.graph_filename<<endl;
  load_from_file(dbg2, p.graph_filename);


  cerr << "k             : " << dbg2.k << endl;
  cerr << "num_nodes()   : " << dbg2.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg2.num_edges() << endl;
  bs = dbg2.bit_size();
  cerr << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
  cerr << "Bits per edge : " << bs / static_cast<double>(dbg2.num_edges()) << " Bits" << endl;

  /*for (size_t i = 0; i<100; i++)
       cout<<i<<" "<<dbg.edge_label(i)<< " " <<dbg2.edge_label(i)<<endl;*/

  size_t nVerify = dbg.num_nodes();
  cerr << "Verifying edges of two dynamic graphs..." << endl;
for (size_t i = 0; i < dbg.num_edges(); i+=dbg.num_edges()/nVerify) {
   if (dbg.edge_label(i) != dbg2.edge_label(i)) {
cerr << "Edge verification failed.\n";
exit(0);
   }
}

cerr << "Verifying nodes of two dynamic graphs..." << endl;
for (size_t i = 0; i < dbg.num_nodes(); i+=dbg.num_nodes()/nVerify) {
   if (dbg.node_label(i) != dbg2.node_label(i)) {
cerr << "Node verification failed.\n";
exit(0);
   }
}

cerr << "Verifying outdegree of two dynamic graphs..." << endl;
for (size_t i = 0; i < dbg.num_nodes(); i+=dbg.num_nodes()/nVerify) {
   if ( dbg.outdegree( i ) != dbg2.outdegree(i) ) {
cerr << "outdegree verification failed.\n";
exit(0);
   }
}

cerr << "Verifying indegree of two dynamic graphs..." << endl;
for (size_t i = 0; i < dbg.num_nodes(); i+=dbg.num_nodes()/nVerify) {
   if (dbg.indegree( i ) != dbg2.indegree(i)) {
cerr << "indegree verification failed.\n";
exit(0);
   }
}

cerr << "Verifying incoming of two dynamic graphs..." << endl;
for (size_t i = 0; i < dbg.num_nodes(); i+=dbg.num_nodes()/nVerify) {
   for (uint64_t j = 0; j < 5; ++j) {
  //cerr << i << ".incoming(" << j << "): " << dbg.incoming(i, j) << endl;
if (dbg.incoming( i,j ) != dbg2.incoming(i,j)) {
   cerr << "incoming verification failed.\n";
   exit(0);
}

   }
}

cerr << "Verifying outgoing of two dynamic graphs..." << endl;
for (size_t i = 0; i < dbg.num_nodes(); i+=dbg.num_nodes()/nVerify) {
   for (uint64_t j = 0; j < 5; ++j) {
  //cerr << i << ".incoming(" << j << "): " << dbg.incoming(i, j) << endl;
if (dbg.outgoing( i,j ) != dbg2.outgoing(i,j)) {
   cerr << "outgoing verification failed.\n";
   exit(0);
}

   }
}

cerr << "Verification passed!\n";
cerr << "===============================\n";

}
