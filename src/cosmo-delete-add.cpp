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

/*This code deletes random edges of the DynamicBOSS with delete_edge, and then add
them back to the graph with add_edge and confirms that the starting and
final graphs are the same.
usage: ./cosmo-delete-add <.dbg file> <number of random edges to delete and add>*/


struct parameters_t {
  std::string input_filename = "";
  std::string number_of_kmers = "";
};
void parse_arguments(int argc, char **argv, parameters_t & params);
void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd("DynamicBOSS. Copyright (c) Bahar Alipanahi, Alan Kuhnle, Alex Bowe 2019", ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".dbg file (output from cosmo-build-dyn).", true, "", "dbg file", cmd);
  TCLAP::UnlabeledValueArg<std::string> number_of_kmers_arg("num_kmers",
            "number of kmers to process dynamically.", true, "", "number_of_kmers", cmd);
  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
  params.number_of_kmers   = number_of_kmers_arg.getValue();
}



int main(int argc, char* argv[]) {
    parameters_t p;
    parse_arguments(argc, argv, p);

    cout<<"Loading DynamicBOSS from file:  "<<p.input_filename<<endl;
    dyn_boss dbg;
    ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
    load_from_file(dbg, p.input_filename);

    //cout << "Constructing dynamic BOSS..." << endl;
    //ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
    //dyn_boss dbg;
    //dbg.load_from_packed_edges( input, "$ACGT" );

    input.close();
    cout << "original graph: " << endl;
    cout << "k             : " << dbg.k << endl;
    cout << "num_nodes()   : " << dbg.num_nodes() << endl;
    cout << "num_edges()   : " << dbg.num_edges() << endl;
    size_t bs = dbg.bit_size();
    cout << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
    cout << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;
    cout << "===============================\n";

    dyn_boss sdbg = dbg;
    vector<string>all_kmers;
    size_t nOps = stoi(p.number_of_kmers) ;
    cout<<"Selecting "<<nOps<<" random edges to delete from the DynamicBOSS ..."<<endl;
    while (all_kmers.size() < nOps ){
       size_t pos;
       string kmer;
       do {
	  pos = rand() % dbg.num_edges();
	  kmer = dbg.edge_label( pos );
       } while ( kmer.find('$') != string::npos );

       all_kmers.push_back(kmer);
    }

    cout << "Beginning of deletion ..." << endl;

    clock_t t_start = clock();
    for (size_t i = 0; i< all_kmers.size();i++){
       dbg.delete_edge( all_kmers[i]);
    }
    double t_elapsed = (clock() - t_start) / CLOCKS_PER_SEC;

    cout << "DONE with deletion of all kmers\n";
    cout << "Time per Op: " << t_elapsed / nOps << endl;
    cout << "===============================\n";

    cout << "Beginning of addition ..." << endl;
    t_start = clock();
    for (size_t i = 0; i< all_kmers.size();i++){
       dbg.add_edge(all_kmers[i]);
    }
    t_elapsed = (clock() - t_start) / CLOCKS_PER_SEC;

    cout << "DONE with addition of all kmers\n";
    cout << "Time per Op: " << t_elapsed / nOps << endl;

    cout << "===============================\n";
    cout << "new greph     : " << endl;
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
