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

// usage: ./cosmo-functions <function> <.dbg file> <.fasta file>

string extension = "updated";
struct parameters_t {
   std::string function = "";
   std::string input_filename = "";
   std::string kmer_filename = "";
   std::string output_prefix = "";
   std::string queryfile= "queryResults.tsv";
};

void parse_arguments(int argc, char **argv, parameters_t & params);
void parse_arguments(int argc, char **argv, parameters_t & params)
{
   TCLAP::CmdLine cmd("DynamicBOSS. Copyright (c) Bahar Alipanahi, Alan Kuhnle, Alex Bowe 2019", ' ', VERSION);
   TCLAP::UnlabeledValueArg<std::string> function_arg("functions",
            ".enter add, delete or query.", true, "", "functions", cmd);
   TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".dbg file (output from cosmo-build-dyn).", true, "", "dbg file", cmd);
   TCLAP::UnlabeledValueArg<std::string> kmer_filename_arg("kmers",
            ".fasta file to count kmers to add and delete.", true, "", "fasta file", cmd);
   string output_short_form = "output_prefix";
   TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. The updated graph will be written to [" + output_short_form + "]" + extension + ". " +
            "Default prefix: basename(input_file).", false, "", output_short_form, cmd);

   cmd.parse( argc, argv );
   params.function  = function_arg.getValue();
   params.input_filename  = input_filename_arg.getValue();
   params.kmer_filename  = kmer_filename_arg.getValue();
}


int main(int argc, char* argv[]){
  parameters_t p;
  parse_arguments(argc, argv, p);
  if (p.function != "add" && p.function != "delete" && p.function != "query"){
    cout<<"Please insert one of the functions add, delete, or query\n";
    exit(0);
  }
  cout<<"Function "<<p.function<<" is called."<<endl;

  cerr<<"Loading DynamicBOSS from file: "<<p.input_filename<<endl;
  dyn_boss dbg;
  ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
  load_from_file(dbg, p.input_filename);
  //cout << "Constructing dynamic BOSS from DSK input..." << endl;
  //dbg.load_from_packed_edges( input, "$ACGT" );
  input.close();
  cout << "Original graph: " << endl;
  cout << "k             : " << dbg.k << endl;
  cout << "num_nodes()   : " << dbg.num_nodes() << endl;
  cout << "num_edges()   : " << dbg.num_edges() << endl;
  size_t bs = dbg.bit_size();
  cout << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
  cout << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;
  //cout << "matrix:\n";
  //  dbg.print_boss_matrix( cerr );
  cout << "===============================\n";
  cout << "Reading FASTA file " <<p.kmer_filename<< endl;
  size_t nKmers;
  set<string> kmers;
  getKmers( nKmers, dbg.k, kmers, p.kmer_filename );
  cout << nKmers << " distinct kmers were counted " << endl;
  cout << "===============================\n";

  double time_elapsed = 0;
  size_t counter = 0;
  size_t available_kmers = 0;
  size_t absent_kmers = 0;
  bool warning = 0;


  if (p.function == "add"){
    for (auto it = kmers.begin(); it != kmers.end(); ++it) {
        if (counter % 1000 == 0) {
          cerr<<"\r           \r";
          cerr<<(double)counter*100/kmers.size()<<"% of kmers were processed";
        }
        string kmer = *it;
        if (!dbg.index_edge_alan( kmer.begin() )) {

          clock_t t_start = clock();
          dbg.add_edge( kmer );
          time_elapsed += double (clock() - t_start);

        }
        else
          available_kmers+=1;
        counter+=1;
    }
    cerr<<"\n100% of kmers were processed\n";
    if (available_kmers > 0){
      cout<<available_kmers<<" ("<<available_kmers*100/kmers.size()<<"%) of the kmers were already in the graph\n";
      if (available_kmers == kmers.size()){
        warning = 1;
        cout<<"Warning! all kmers were already in the graph!\n";
      }
    }
    cout<<kmers.size()- available_kmers<< " kmers added to the graph in "<< time_elapsed/CLOCKS_PER_SEC<<" (s)"<<endl;
    cout<<"Time per operaition: "<<(time_elapsed/CLOCKS_PER_SEC)/kmers.size()<<" (s)"<<endl;

  }

  if (p.function == "delete"){
    for (auto it = kmers.begin(); it != kmers.end(); ++it) {
        if (counter % 1000 == 0) {
          cerr<<"\r           \r";
          cerr<<(double)counter*100/kmers.size()<<"% of kmers were processed";
        }
        string kmer = *it;
        if (dbg.index_edge_alan( kmer.begin() )) {

          clock_t t_start = clock();
          dbg.delete_edge( kmer );
          time_elapsed += double (clock() - t_start);

        }
        else
          absent_kmers+=1;
        counter+=1;
    }
    cerr<<"\n100% of kmers were processed\n";
    if (absent_kmers > 0){
      cout<<absent_kmers<<" ("<<absent_kmers*100/kmers.size()<<"%) of the kmers were not in the graph\n";
      if (absent_kmers == kmers.size()){
        warning = 1;
        cout<<"Warning! none of the kmers were in the graph!\n";
      }
    }
    cout<<kmers.size()- absent_kmers<< " kmers deleted from the graph in "<< time_elapsed/CLOCKS_PER_SEC<<" (s)"<<endl;
    cout<<"Time per operaition: "<<(time_elapsed/CLOCKS_PER_SEC)/kmers.size()<<" (s)"<<endl;

}
if (p.function == "query"){
  ofstream ofs(p.queryfile);
  ofs<<"lexo order"<<"\t"<<"kmer"<<"\t"<<"presence"<<endl;
  for (auto it = kmers.begin(); it != kmers.end(); ++it) {
      if (counter % 1000 == 0) {
        cerr<<"\r           \r";
        cerr<<(double)counter*100/kmers.size()<<"% of kmers were processed";
      }
      string kmer = *it;

      clock_t t_start = clock();
      bool present = dbg.index_edge_alan( kmer.begin() );
      time_elapsed += double (clock() - t_start);

      ofs<<counter<<"\t"<<kmer<<"\t"<<present<<endl;
      if (present == 0)
        absent_kmers+=1;
      counter+=1;
  }
  cerr<<"\n100% of kmers were processed\n";
  cout<<kmers.size()-absent_kmers<<" ("<<(kmers.size()-absent_kmers)*100/kmers.size()<<"%) of the kmers were in the graph\n";
  cout<<kmers.size()<< " kmers were queried in "<< time_elapsed/CLOCKS_PER_SEC<<" (s)"<<endl;
  cout<<"Time per operaition: "<<(time_elapsed/CLOCKS_PER_SEC)/kmers.size()<<" (s)"<<endl;
  ofs.close();

}
cout << "===============================\n";
if ((p.function == "add" || p.function == "delete") && !warning){
  cout << "Updated graph: " << endl;
  cout << "k             : " << dbg.k << endl;
  cout << "num_nodes()   : " << dbg.num_nodes() << endl;
  cout << "num_edges()   : " << dbg.num_edges() << endl;
  bs = dbg.bit_size();
  cout << "Total size    : " << bs / 8.0 / 1024.0 / 1024.0 << " MB" << endl;
  cout << "Bits per edge : " << bs / static_cast<double>(dbg.num_edges()) << " Bits" << endl;

  string outfilename = p.input_filename+ "."+extension;
  cout<<"Writing the updated DynamicBOSS in file: "<<outfilename<<endl;
  ofstream ofs( outfilename, ios::out | ios::binary );
  dbg.serialize( ofs );
  ofs.close();
}
if (p.function == "query")
  cout<<"Writing the querying results in file: "<<p.queryfile<<endl;

cout<<"DONE!"<<endl;


}
