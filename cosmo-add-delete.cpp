#include <iostream>
#include <fstream>

#include <libgen.h> // basename

#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "io.hpp"
//#ifdef DYNAMIC
#include "dynBoss.hpp"
//#else
//#include "debruijn_graph.hpp"
//endif
#include "algorithm.hpp"
#include "utility.hpp"
#include "formatutil.cpp"
using namespace std;
using namespace sdsl;

static char base[] = {'$','A','C','G','T'};

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
dyn_boss Delete_Edge (dyn_boss dbg, std::string kmer);
dyn_boss Delete_node (dyn_boss dbg, std::string kmer);
dyn_boss Add_Edge (dyn_boss dbg,  std::string kmer, bool addition){
    pair<size_t, size_t> start_and_position = dbg.index_finder_for_add_delete(kmer.begin() , 1);
    size_t start = start_and_position.first;
    size_t pos = start_and_position.second;
    vector<size_t> to_add;
    for (size_t i = pos; i < kmer.size(); i++) {
        if (kmer[i]!= '$')
        to_add.push_back (dbg._encode_symbol(kmer[i]));
    }

    size_t i = 0;
    while (i < to_add.size()){
        i++;
        kmer = dbg.node_label(dbg._edge_to_node(start)) += dbg._map_symbol(to_add[i-1]);
        start = dbg.add(start,to_add[i-1],kmer);
        if (dbg._forward(start) == 0)  start++; //Do we need this?
        start = dbg._forward(start);
    }

    //if (dbg.index(kmer.begin(),0) == 1) cout<<"now the node is added"<<endl;
    //else cout<<"still not here"<<endl;

    //If added edge is now at the begining of a read, we should delete the dummies
    //of former first edge of the read (only if their outdegree is 1)
    bool remove_dummy = true;
    if (dbg.outdegree(dbg._edge_to_node(start)) > 1 ) remove_dummy = false;
    string incoming_dummy = '$'+ kmer.substr(1);

    //start is replaced with _forward(start) already at the final step of addition, so to check the indegree "dbg.indegree(dbg._edge_to_node(start))" is correct.
    if (addition && dbg.indegree(dbg._edge_to_node(start)) > 1 && dbg.index(incoming_dummy.begin(),0) == 1 && remove_dummy){
        pair<size_t, size_t> start_and_position = dbg.index_finder_for_add_delete(incoming_dummy.begin() , 0);
        size_t start = start_and_position.first;
        size_t i = 0;
        bool right_place = true; //always start < ref because $A is before AA

        while ( i < dbg.k){
            size_t out = dbg.outdegree(dbg._edge_to_node(start));
            dbg.delete_edge(start,incoming_dummy[dbg.k-1-i]);
            size_t backward_edge;
            start = (right_place) ? start : start-1; //because one edge before start is already deleted
            backward_edge = dbg._backward(start);
            right_place = (backward_edge >= start) ?  false : true;
            start = backward_edge;
            i++;
            if (out > 1 ) break;
          }
    }

    return (dbg);
}

dyn_boss Delete_Edge (dyn_boss dbg, std::string kmer){
    pair<size_t, size_t> start_and_position = dbg.index_finder_for_add_delete(kmer.begin() , 0);
    size_t start = start_and_position.first;
    size_t pos = start_and_position.second;

    //adding dummies of next edge
    if (dbg.edge_label(dbg._forward(start))[dbg.k-1] != '$' && dbg.indegree(dbg._edge_to_node(dbg._forward(start))) == 1) {
        string added_kmer = "$";
        added_kmer += dbg.node_label(dbg._edge_to_node(dbg._forward(start)));
        dbg = Add_Edge(dbg,added_kmer,0);
        start_and_position = dbg.index_finder_for_add_delete(kmer.begin() , 0);
        start = start_and_position.first;
        size_t dif = 0;
        if(dbg._encode_symbol(kmer[dbg.k-1]) - dbg._encode_symbol(dbg.edge_label(start)[dbg.k-1]) > 0) {
            dif = dbg._encode_symbol(kmer[dbg.k-1]) - dbg._encode_symbol(dbg.edge_label(start)[dbg.k-1]);
            start += dif;
          } //careful;
    }

    bool remove_dummy = true;
    if (dbg.outdegree(dbg._edge_to_node(start)) > 1 ) remove_dummy = false;
    string incoming_dummy2 = dbg.edge_label(dbg._backward(start));
    dbg.delete_edge(start,dbg._encode_symbol(kmer[dbg.k-1]));

    //if (dbg.index(kmer.begin(), 0) != 1 ) cout<<"now the edge is deleted"<<endl;

    //if deleting the first edge of a read, all dummies should be deleted
    string incoming_dummy = '$'+ kmer.substr(0,dbg.k-1);
    if (dbg.index(incoming_dummy.begin(),0) == 1 && remove_dummy){
        vector<string> backup;
        size_t ref = start;
        pair<size_t, size_t> start_and_position = dbg.index_finder_for_add_delete(incoming_dummy.begin() , 0);
        size_t start = start_and_position.first;
        size_t i = 0;
        bool right_place;
        right_place = (start >= ref) ?  false : true;

        while (i < dbg.k){
            size_t out = dbg.outdegree(dbg._edge_to_node(start));
             dbg.delete_edge(start,incoming_dummy[dbg.k-1-i]);
             if (out > 1) break;
             size_t backward_edge;
             start = (right_place) ? start : start-1; //because one edge before start is already deleted
             backward_edge = dbg._backward(start);
             right_place = (backward_edge >= start) ?  false : true;
             start = backward_edge;
              i++;
         }
    }
    return (dbg);
}

dyn_boss Delete_node (dyn_boss dbg, std::string kmer){
    for (size_t i = 0; i < 5; i++){
        string for_node = kmer.substr(1) + base[i];
        if (dbg.index(for_node.begin() , 0) == 1 ){
            dbg = Delete_Edge(dbg , kmer+base[i]);
        }
    }
    for (size_t i = 0; i < 5; i++){
        string back_node = base[i] + kmer.substr(0,dbg.k-2);
        if (dbg.index(back_node.begin() , 0) == 1 ){
            dbg = Delete_Edge(dbg , base[i]+kmer);
        }
    }
    return (dbg);
}

int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);

  ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
    dyn_boss dbg;
    dbg.load_from_packed_edges( input, "$ACGT" );
    input.close();
    cerr << "original graph: " << endl;
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
    while (all_kmers.size() < 1000 ){
        string kmer;
        while(kmer.size() <dbg.k )
            kmer += alphanum[rand() % stringLength];
        all_kmers.push_back(kmer);
          }

    cerr<<"number of kmers to process: "<<all_kmers.size()<<endl;
    for (size_t i; i< all_kmers.size();i++){
      if (dbg.index(all_kmers[i].begin(),0) == 0)
        dbg = Add_Edge (dbg, all_kmers[i] ,1);
      if (dbg.index(all_kmers[i].begin(),0) == 0){
        cerr<<"failed to add kmer "<<i<<endl;
        exit(0);
      }
    }
    cerr<<"DONE with addition of all kmers\n";

    for (size_t i; i< all_kmers.size();i++){
     if (dbg.index(all_kmers[i].begin(),0) == 1)
         dbg = Delete_Edge (dbg, all_kmers[i]);
     if (dbg.index(all_kmers[i].begin(),0) == 1){
       cerr<<"failed to delete kmer "<<i<<endl;
       exit(0);
         }
       }
    cerr << "DONE with deletion of all kmers\n";
    cerr << "===============================\n";
    cerr << "new greph     : " << endl;
    cerr << "k             : " << dbg.k << endl;
    cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
    cerr << "num_edges()   : " << dbg.num_edges() << endl;
    bs = dbg.bit_size();
    cerr << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
    cerr << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;

    /* for (size_t i = 0; i<dbg.num_edges(); i++)
         cout<<dbg.edge_label(i)<< " " <<i<<endl;

     for (size_t i=0; i < dbg.num_nodes();i++)
         cout<<dbg.node_label(i)<< " " <<i<<endl;*/

	cerr << "Verifying edges..." << endl;
  for (size_t i = 0; i < dbg.num_edges(); i+=dbg.num_edges()/5) {
     if (dbg.edge_label(i) != sdbg.edge_label(i)) {
	cerr << "Edge verification failed.\n";
	exit(0);
     }
  }

  cerr << "Verifying nodes..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); i+=dbg.num_nodes()/5) {
     if (dbg.node_label(i) != sdbg.node_label(i)) {
	cerr << "Node verification failed.\n";
	exit(0);
     }
  }

  cerr << "Verifying outdegree..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); i+=dbg.num_nodes()/5) {
     if ( dbg.outdegree( i ) != sdbg.outdegree(i) ) {
	cerr << "outdegree verification failed.\n";
	exit(0);
     }
  }

  cerr << "Verifying indegree..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); i+=dbg.num_nodes()/5) {
     if (dbg.indegree( i ) != sdbg.indegree(i)) {
	cerr << "indegree verification failed.\n";
	exit(0);
     }
  }

  cerr << "Verifying incoming..." << endl;
  for (size_t i = 0; i < dbg.num_nodes(); i+=dbg.num_nodes()/5) {
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
cerr << "Deleting some randome nodes\n";
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
