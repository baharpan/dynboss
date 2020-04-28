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

/*This code once query the random kmers from available edges (member_kmers) of
the DynamicBOSS and then query the random kmers counted from a .fasta file
(random_kmers).
usage: ./cosmo-index <.dbg file> <number of member_kmers> <.fasta file> */


struct parameters_t {
  std::string input_filename = "";
  std::string kmer_filename = "";
  std::string number_of_kmers = "";
};


void parse_arguments(int argc, char **argv, parameters_t & params);
void parse_arguments(int argc, char **argv, parameters_t & params)

{

  TCLAP::CmdLine cmd("DynamicBOSS. Copyright (c) Bahar Alipanahi, Alan Kuhnle, Alex Bowe 2019", ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".dbg file (output from cosmo-build-dyn).", true, "", "dbg file", cmd);
  TCLAP::UnlabeledValueArg<std::string> number_of_kmers_arg("num_kmers",
            "number of memebr_kmers to query.", true, "", "number_of_kmers", cmd);
  TCLAP::UnlabeledValueArg<std::string> kmer_filename_arg("kmers",
              ".fasta file to count random_kmers.", true, "", "fasta file", cmd);

  cmd.parse( argc, argv );
  params.input_filename  = input_filename_arg.getValue();
  params.kmer_filename   = kmer_filename_arg.getValue();
  params.number_of_kmers   = number_of_kmers_arg.getValue();
}


int main(int argc, char* argv[]) {
    parameters_t p;
    parse_arguments(argc, argv, p);
    cout<<"Loading DynamicBOSS from file:  "<<p.input_filename<<endl;
    dyn_boss dbg;
    ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
    load_from_file(dbg, p.input_filename);
    input.close();



    cout << "k             : " << dbg.k << endl;
    cout << "num_nodes()   : " << dbg.num_nodes() << endl;
    cout << "num_edges()   : " << dbg.num_edges() << endl;
    size_t bs = dbg.bit_size();
    cout << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
    cout << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;
    cout<<"==============================="<<endl;



    vector<string>all_kmers;
    size_t nOps = stoi(p.number_of_kmers) ;

    cout<<"Selecting "<<nOps<<" random member kmers to query ..."<<endl;
    while (all_kmers.size() < nOps ){
       size_t pos;
       string kmer;
        do {
          pos = rand() % dbg.num_edges();
          kmer = dbg.edge_label( pos );

       } while ( kmer.find('$') != string::npos );

       all_kmers.push_back(kmer);

    }

    cout << "Beginning of querying..." << endl;
    clock_t t_start = clock();

    for (size_t i = 0; i< all_kmers.size();i++){
      dbg.index_edge_alan( all_kmers[i].begin());
    }

    double t_elapsed = (clock() - t_start) / CLOCKS_PER_SEC;

    cout << "DONE with querying "<<nOps<< " of member_kmers\n";
    cout << "Time per Op: " << t_elapsed / nOps << endl;
    cout<<"==============================="<<endl;


    cout << "Reading FASTA file " <<p.kmer_filename<< endl;
    size_t nKmers;
    set<string>kmers;

    getKmers( nKmers, dbg.k, kmers, p.kmer_filename );
    cout << nKmers << " kmers were counted " << endl;

    cout<<"Beginning of querying..."<<endl;
    t_start = clock();

    for (auto it = kmers.begin(); it != kmers.end(); ++it) {
        string kmer = *it;
        dbg.index_edge_alan(kmer.begin());
      }

    t_elapsed = (clock() - t_start) / CLOCKS_PER_SEC;

    cout << "DONE with querying "<<kmers.size()<<" of random_kmers\n";
    cout << "Time per Op: " << t_elapsed /kmers.size() << endl;
}
