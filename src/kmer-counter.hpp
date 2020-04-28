#include <iostream>
#include <fstream>
using namespace std;

void getKmers( size_t& nKmers, size_t k, set< string >& kmers, string& p) {
   ifstream in;
   in.open(p);
   if(in.fail())
    {
        cout<<"Could not read the fasta file called "<<p<<endl;
        exit(0);
    }
   string sline;
   vector< string > vline;
   if (!getline(in, sline) || !(sline[0] == '>')) {
          cout << "Make sure that "<<p<< " is a fatsa file" << endl;
          exit(0);
        }
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
      if (read_length < k){
        cout<<"Warning! There is a read in file "<<p<<" that is shorter than k="<<k<<"!"<<endl;
        exit(0);
      }

      for (size_t i = 0; i < read_length; ++i) {
	       if (sline[i] == 'N'){
	        sline[i] = 'A';
        }
      }
      size_t nMers = read_length - k + 1;
      for (size_t start = 0; start < nMers; ++start) {
	       string kmer = sline.substr( start, k );
	           kmers.insert( kmer );
      }
   }

   in.close();
   nKmers = kmers.size();
}
