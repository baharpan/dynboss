#pragma once
#ifndef _DYNAMIC_BOSS_H
#define _DYNAMIC_BOSS_H

#include <algorithm>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <iostream>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include <dynamic.hpp>
#include "algorithm.hpp"
#include "utility.hpp"
#include "io.hpp"
#include "debug.h"

using namespace std;
using namespace sdsl;
//using namespace dyn;

class dyn_boss {
public:

   const static size_t sigma = 4;
   typedef size_t edge_type;
   typedef pair<edge_type, edge_type> node_type;
   typedef dyn::wt_str::char_type symbol_type;

   size_t           k;
   dyn::gap_bv* p_node_flags;
   //dyn::suc_bv* p_node_flags;
   dyn::wt_str* p_edges;

   array<size_t, 1 + sigma> m_symbol_ends;
   array<size_t, 1 + sigma> m_edge_max_ranks;
   string             m_alphabet;
   size_t             m_num_nodes;
public:
   dyn_boss() {
      p_node_flags = NULL;
      p_edges = NULL;
      m_alphabet = "";
      m_num_nodes = 0;
      k = 0;
   }

   dyn_boss( size_t in_k,
	     dyn::gap_bv* node_flags,
	     dyn::wt_str* edges,
	     array<size_t, 1 + sigma>& symbol_ends,
	     string& alphabet) :
      k( in_k ),
      p_node_flags( node_flags ),
      p_edges( edges ),
      m_symbol_ends(symbol_ends),
      m_alphabet( alphabet )
   {
      _init_max_ranks();
      m_num_nodes = p_node_flags->rank0();
   }

   void _init_max_ranks() {
      //array<size_t, 1+sigma> max_ranks;
      size_t num_edges = p_edges->size();
      for (symbol_type x = 0; x<sigma+1;x++) {
	 this->m_edge_max_ranks[x] = p_edges->rank(num_edges, _with_edge_flag(x, false));
      }
   }

   void load_from_packed_edges(
			       istream & input,
			       string alphabet ) {
      m_alphabet = alphabet;
      // ifstream input(filename, ios::in|ios::binary|ios::ate);
      // check length
      streampos size = input.tellg();
      // should be exceptions...
      assert(size/sizeof(uint64_t) >= (sigma+2)); // space for footer info
      assert(((size_t)size - (sigma+2) * sizeof(uint64_t))%sizeof(uint64_t) == 0); // sequence of uint64_ts

      // read footer
      input.seekg(-(sigma+2) * sizeof(uint64_t), ios::end);
      array< size_t, 1 + sigma >& counts = m_symbol_ends;

      k = 0;
      input.read((char*)&counts[0], (sigma+1) * sizeof(uint64_t));
      input.read((char*)&k, sizeof(uint64_t));

      size_t num_edges = counts[sigma];
#ifndef NDEBUG
      cerr << "counts\n";
      for (size_t i = 0; i <= sigma; ++i) {
	 cerr << i << ' ' << counts[i] << endl;
      }
#endif

      size_t num_blocks = size_t(size)/sizeof(uint64_t) - (sigma+2);
      input.seekg(0, ios::beg); // rewind

      vector<uint64_t> blocks(num_blocks,0);
      input.read((char*)&blocks[0], sizeof(uint64_t) * num_blocks);

      // So we avoid a huge malloc if someone gives us a bad file
      p_node_flags = new dyn::gap_bv();

      // would be nice to fix wavelet trees so the constructor
      // can accept a int_vector<4> instead (which is all we need for DNA)
      p_edges = new dyn::wt_str( 10 );
      //p_edges = new dyn::wt_str( 5 );
      //p_edges = NULL;

      for (size_t i = 0; i < num_edges; i++) {
	 auto x = get_edge(blocks.begin(), i);
	 p_node_flags->push_back( 1-get<1>(x) ); // convert 0s to 1s so we can have a sparse bit vector

	 p_edges->push_back( (get<0>(x) << 1) | !get<2>(x) );
	 //p_edges->push_back( get<0>(x) );
      }

      _init_max_ranks();
      m_num_nodes = p_node_flags->rank0();
   }

   // API
   size_t outdegree(size_t v) {
      auto range = _node_range(v);
      size_t first = get<0>(range);
      size_t last  = get<1>(range);
      size_t count = last - first + 1;
      // This edge access is kinda annoying, but if the ONE edge is $, we don't really have an outgoing edge
      return count - (count == 1 && _strip_edge_flag((*p_edges)[first]) == 0);
   }

   vector<node_type> all_preds(const node_type & v) {

      // node u -> v : edge i -> j
      size_t j = get<0>(v);
      symbol_type y = _symbol_access(j);
      if (y == 0) return vector<node_type>(0);
      size_t i_first = _backward(j);
      size_t i_last  = _next_edge(i_first, y);
      size_t base_rank = p_edges->rank(i_first, _with_edge_flag(y, true));
      size_t last_rank = p_edges->rank(i_last, _with_edge_flag(y, true));
      size_t num_predecessors = last_rank - base_rank + 1;
      // binary search over first -> first + count;
      auto selector = [&](size_t i) -> size_t {
	 return (i == 0)? i_first : p_edges->select(base_rank+i, _with_edge_flag(y,true));
      };

      vector<node_type> result(num_predecessors);
      for (size_t i = 0; i<num_predecessors; i++) {
	 edge_type e_i = selector(i);
	 edge_type e_j = _last_edge_of_node(_edge_to_node(e_i));
	 result.push_back(node_type(e_i, e_j));
      }
      return result;
   }

   size_t indegree(size_t v) {
      // find first predecessor edge i->j
      size_t j = _node_to_edge(v);

      // edge label has to be the last node symbol of v
      symbol_type x = _symbol_access(j);

      if (x == 0) return 0;
      size_t i_first = _backward(j);

      size_t i_last = _next_edge(i_first, x);
      if (i_last = num_edges()) i_last = i_first; //Bahar
      return p_edges->rank(i_last, _with_edge_flag(x, true)) -
      p_edges->rank(i_first, _with_edge_flag(x, true)) + 1;
   }

   // The signed return type is worrying, since it halves the possible answers...
   // but it will only face problems with graphs that have over 2^63 ~= 4^31 edges,
   // which is a saturated debruijn graph of k=31. Which is unlikely.
   // Hopefully an iterator-based API (like BGL) will fix this issue
   ssize_t outgoing(size_t u, symbol_type x) const {
      assert(x < sigma + 1);
      if (x == 0) return -1;
      auto range = _node_range(u);
      size_t first = get<0>(range);
      size_t last  = get<1>(range);


      // Try both with and without a flag
      for (symbol_type c = _with_edge_flag(x,false); c <= _with_edge_flag(x, true); c++) {
	 //size_t most_recent = p_edges->select(p_edges->rank(last+1, c), c);

	 size_t tmpRank = p_edges->rank(last+1,c);
   if (tmpRank == 0) continue;
	 //if (tmpRank > 0) {
	    size_t most_recent = p_edges->select(p_edges->rank(last+1, c) - 1, c);


	    // if within range, follow forward
	    if (first <= most_recent && most_recent <= last) {

	       // Don't have to check fwd for -1 since we checked for $ above
	       return _forward(most_recent);
	    }
	 //}
      }
      return -1;
   }

   // Added for DCC. Will remove the other one later and rename this one.
   ssize_t interval_node_outgoing(const node_type & u, symbol_type x) const {
      assert(x < sigma + 1);
      if (x == 0) return -1;
      //auto range = _node_range(u);
      size_t first = get<0>(u);
      size_t last  = get<1>(u);
      // Try both with and without a flag
      for (symbol_type c = _with_edge_flag(x,false); c <= _with_edge_flag(x, true); c++) {
	 size_t most_recent = p_edges->select(p_edges->rank(last+1, c), c);
	 // if within range, follow forward
	 if (first <= most_recent && most_recent <= last) {
	    // Don't have to check fwd for -1 since we checked for $ above
	    return _edge_to_node(_forward(most_recent));
	 }
      }
      return -1;
   }



   // For DCC
   ssize_t _outgoing_edge_pair(size_t first, size_t last, symbol_type x) const {
      // Try both with and without a flag
      for (symbol_type c = _with_edge_flag(x,false); c <= _with_edge_flag(x, true); c++) {
	 size_t most_recent = p_edges->select(p_edges->rank(last+1, c), c);
	 // if within range, follow forward
	 if (first <= most_recent && most_recent <= last) {
	    // Don't have to check fwd for -1 since we checked for $ above
	    return _forward(most_recent);
	 }
      }
      return -1;
   }

   // incoming
   ssize_t incoming(size_t v, symbol_type x) {
      // This is very similar to indegree, so should maybe be refactored
      assert(x < sigma + 1);
      // node u -> v : edge i -> j
      size_t j = _node_to_edge(v);
      symbol_type y = _symbol_access(j);
      if (y == 0) return -1;
      size_t i_first = _backward(j);
      size_t i_last  = _next_edge(i_first, y);
      size_t base_rank = p_edges->rank(i_first, _with_edge_flag(y, true));
      size_t last_rank = p_edges->rank(i_last, _with_edge_flag(y, true));
      size_t num_predecessors = last_rank - base_rank + 1;
      // binary search over first -> first + count;
      auto selector = [&](size_t i) -> size_t {
	 return (i == 0)? i_first : p_edges->select(base_rank+i-1, _with_edge_flag(y,true));
      };
      auto accessor = [&](size_t i) -> symbol_type { return _first_symbol(selector(i)); };
      ssize_t sub_idx = function_binary_search(0, num_predecessors-1, x, accessor);
      if (sub_idx == -1) return -1;
      return _edge_to_node(selector(sub_idx));
   }




   // string -> node, edge
   // BGL style API

   string node_label(size_t v) {
      size_t i = _node_to_edge(v);
      string label = string(k-1, _map_symbol(symbol_type{}));
      return _node_label_from_edge_given_buffer(i, label);
   }

   string node_label_from_edge(size_t i) {
      string label = string(k-1, _map_symbol(symbol_type{}));
      return _node_label_from_edge_given_buffer(i, label);
   }

   string edge_label(size_t i) {
      string label = string(k, _map_symbol(symbol_type{}));
      _node_label_from_edge_given_buffer(i, label);
      label[k-1] = _map_symbol(_strip_edge_flag((*p_edges)[i]));
      return label;
   }

   size_t num_edges() const { return m_symbol_ends[sigma]; /*_node_flags.size();*/ }
   size_t num_nodes() const { return m_num_nodes; /*m_node_rank(num_edges());*/ }

private:
   size_t _symbol_start(symbol_type x) const {
      assert(x < sigma + 1);
      return (x==0)? 0 : m_symbol_ends[x - 1];
   }

public:
   size_t _node_to_edge(size_t v) {
      //return m_node_flags.select0(v+1);
      return p_node_flags->select0(v);
   }

   size_t _edge_to_node(size_t i) const {
      assert(i < num_edges());
      if (i == num_edges()-1) return num_nodes()-1;
      return p_node_flags->rank0(i+1)-1;
   }

   // This should be moved to a helper file...
   symbol_type _strip_edge_flag(symbol_type x) const {
      return x >> 1;
   }

   // False -> normal edge (but still shifted to fit in the edge alphabet)
   // True -> minus flag
   symbol_type _with_edge_flag(symbol_type x, bool edge_flag) const {
      return (x << 1) | edge_flag;
   }

   symbol_type _symbol_access(size_t i) const {
      assert(i < num_edges());
      // I assume std binary search is optimised for small ranges (if sigma is small, e.g. DNA)
      return upper_bound(m_symbol_ends.begin(), m_symbol_ends.end(), i) - m_symbol_ends.begin();
   }

   node_type get_node(size_t v) const {
      auto r = _node_range(v);
      return node_type(get<0>(r), get<1>(r));
   }

   // provided for DCC paper
   inline symbol_type lastchar(const node_type & v) const {
      return _symbol_access(get<0>(v));
   }

   // more efficient than generating full label for incoming()
   symbol_type _first_symbol(size_t i) {
      symbol_type x = 0;
      for (size_t pos = 1; pos <= k-1; pos++) {
	 // Don't need to map it to the alphabet in this case
	 x = _symbol_access(i);
	 // All are $ before the last $
	 if (x == 0) return x;
	 i = _backward(i);
      }
      return x;
   }

   string _node_label_from_edge_given_buffer(size_t i, string & label) {
      // Calculate backward k times and fill a buffer with symbol_access(edge)
      //string label = string(k-1, _map_symbol(symbol_type{}));
      for (size_t pos = 1; pos <= k-1; pos++) {
	 symbol_type x = _symbol_access(i);
	 label[k-pos-1] = _map_symbol(x);
	 // All are $ before the last $
	 if (x == 0) return label;
	 i = _backward(i);
      }
      return label;
   }

   size_t _next_edge(size_t i, symbol_type x) const {
      if (i >= num_edges() - 1 ) return i;
      // Might not actually occur if out of rank bounds?
      //size_t next_rank = 1 + p_edges->rank(1+i, _with_edge_flag(x, false));
      size_t next_rank =  p_edges->rank(1+i, _with_edge_flag(x, false)); //Bahar
      //we know next_rank >= 1, since i is ostensibly
      //an edge of symbol_type x
      if (next_rank >= m_edge_max_ranks[x]) return num_edges();
      return p_edges->select(next_rank, _with_edge_flag(x, false));
   }

   size_t _rank_distance(size_t a, size_t b) {
      return p_node_flags->rank0(b) - p_node_flags->rank0(a);
   }

   // Return index of first possible edge obtained by following edge i
public:
   ssize_t _forward(size_t i) const {
      symbol_type temp;
      return _forward(i, temp);
   }

   // This is so we can reuse the symbol lookup - save an access during traversal :)
   ssize_t _forward(size_t i, symbol_type & x) const {
      assert(i < num_edges());
      x = _strip_edge_flag((*p_edges)[i]);

      // if x == 0 ($) then we can't follow the edge
      // (should maybe make backward consistent with this, but using the edge 0 loop for node label generation).
      if (x == 0) return 0;
      size_t start = _symbol_start(x);
      size_t nth   = p_edges->rank(i, _with_edge_flag(x, false));
      if ((*p_edges)[i] % 2 == 1) nth-=1;
      // Bahar
      if (x == 4 && p_node_flags->rank0(start+1) + nth - 1 >= num_nodes()) return num_edges();
      size_t next = p_node_flags->select0(p_node_flags->rank0(start+1) + nth - 1);

      return next;
   }


   size_t _backward(size_t i) {
      assert(i < num_edges());
      symbol_type x  = _symbol_access(i);
      // This handles x = $ so that we have the all-$ edge at position 0
      // but probably shouldn't be called this way in the normal case
      size_t x_start = _symbol_start(x);
      // rank is over [0,i) and select is 1-based
      size_t nth = _rank_distance(x_start+1, i+1);
      if (x == 0) return 0;
      // no minus flag because we want the FIRST
      return p_edges->select(nth, _with_edge_flag(x, false));
	    //return p_edges->select(nth, x );
   }

   size_t backward(size_t v) {
      return _edge_to_node(_backward(_node_to_edge(v)));
   }

   symbol_type _map_symbol(symbol_type x) const {
      return (m_alphabet.size() > 0)? m_alphabet[x] : x;
   }

   size_t _first_edge_of_node(size_t v) const {
      //return p_node_flags->select0(v+1);
      return p_node_flags->select0(v);
   }

   size_t _last_edge_of_node(size_t v) const {
      // find the *next* node's first edge and decrement
      // as long as a next node exists!
      assert(v + 1 <= num_nodes());
      if (v+1 == num_nodes()) return num_edges() - 1;
      else return _first_edge_of_node(v+1) - 1;
   }

   pair<size_t, size_t> _node_range(size_t v) const {
      return make_pair(_first_edge_of_node(v), _last_edge_of_node(v));
   }

   // TODO: add first_sibling and edge_range
   // TODO: update to use rank and select for larger alphabets
   size_t _last_sibling(size_t i) const {
      size_t last = i;
      while(last < num_edges() && (*p_node_flags)[last] != 0) {
	 last++;
      }
      // last should be one past end
      return last-1;
   }

   size_t serialize( ostream& out ) {
      size_t w_bytes = 0;
      out.write((char*)&k,sizeof(k));
      w_bytes += sizeof(k);
      w_bytes += p_node_flags->serialize( out );
      w_bytes += p_edges->serialize( out );
      for (size_t i = 0; i < sigma+1; ++i) {
	 out.write( (char*)(&(m_symbol_ends[i])), sizeof(size_t) );
	 w_bytes += sizeof(size_t);
      }
      for (size_t i = 0; i < sigma+1; ++i) {
	 out.write( (char*)(&(m_edge_max_ranks[i])), sizeof(size_t) );
	 w_bytes += sizeof(size_t);
      }
      size_t asize = m_alphabet.size();
      out.write( (char*)&asize, sizeof(size_t));
      w_bytes += sizeof(size_t);
      out.write( m_alphabet.c_str() , m_alphabet.size() );
      w_bytes += m_alphabet.size();
      out.write( (char*)&m_num_nodes, sizeof(size_t));
      w_bytes += sizeof(size_t);
      return w_bytes;
   }

   size_t size() const { return num_edges(); }

   symbol_type _encode_symbol(uint8_t c) const {
      return lower_bound(m_alphabet.begin(), m_alphabet.end(), c) - m_alphabet.begin();
   }

  void delete_edge(size_t start , string kmer){
      size_t in = indegree(_edge_to_node(_forward(start)));
      size_t out = outdegree(_edge_to_node(start));
      size_t forward_index = _forward(start);
      size_t forward_lab = (*p_edges)[_forward(start)];
      size_t forward_to_be_decreased = _encode_symbol(kmer[31]);//(*p_edges)[start]/2;//_encode_symbol(edge_label(_forward(start))[k-2]);

      if (in > 1 && out > 1){

          size_t to_be_decreased = _encode_symbol(kmer[30]);//edge_label(start)[k-2]);
          p_edges->remove(start);
          for (symbol_type i = 0; i < sigma+1; ++i){
              if ( i >= to_be_decreased){
                  m_symbol_ends[i]--;

              }
          }
          size_t next_incoming_edge = start+ outdegree(_edge_to_node(start))-1;
          size_t lab_edge = (*p_edges)[next_incoming_edge];
          if (lab_edge % 2 == 1){

                  p_edges->remove(next_incoming_edge);
                  p_edges->insert(next_incoming_edge,lab_edge-1);

                 }

          if ((*p_node_flags)[start] == 0){
              p_node_flags->remove(start+1);
              }
              else{
                  p_node_flags->remove(start);
              }
        }
      else{
            if (out == 1 && in == 1){


          p_edges->remove(start);   //for removing any edge with outdegree 1
          p_edges->insert(start,0);
      }



          if (in > 1){

              size_t next_incoming_edge = start+ outdegree(_edge_to_node(start));
              p_edges->remove(start);
              p_edges->insert(start,0);

              size_t lab_edge = (*p_edges)[next_incoming_edge];
              if (lab_edge % 2 == 1){

                  p_edges->remove(next_incoming_edge);
                  p_edges->insert(next_incoming_edge,lab_edge-1);
              }
          }


          if (out >  1 ) {

              size_t to_be_decreased = _encode_symbol(kmer[30]);//edge_label(start)[k-2]);

               p_edges->remove(start);
               for (symbol_type i = 0; i < sigma+1; ++i){
               if ( i >= to_be_decreased){
               m_symbol_ends[i]--;
               }
               }

              if ((*p_node_flags)[start] == 0)
                  p_node_flags->remove(start+1);

              else
                  p_node_flags->remove(start);


      }

            if (in == 1 && forward_lab == 0){ //only for removing last edge of a read

                if (out > 1 && start < forward_index){
                  p_edges->remove(forward_index - 1);
                  p_node_flags->remove(forward_index - 1);
                }
                else{
                  p_edges->remove(forward_index);
                  p_node_flags->remove(forward_index);
                }
              for (symbol_type i = 0; i < sigma+1; ++i){
                  if ( i >= forward_to_be_decreased ){
                      m_symbol_ends[i]--;
                  }
              }
              m_num_nodes-=1;
          }

      }

  }


  size_t add(size_t start, symbol_type x ,string  kmer, bool last_node){

     if (start == 0){
	for (size_t m = 0; m <4; ++m){
	   if (edge_label(m) == string(k-1, '$')+kmer[k-1])
	      return start;
	}
     }

     bool last_node_present = false;
     size_t check;

     if (last_node){
	if (index(kmer.substr(1).begin(),0) == 1) {
	   last_node_present = true;
	   size_t i_last  = _next_edge(start, _encode_symbol(kmer[31]));
	   check = (i_last == num_edges() || i_last - start >= 4) ? 0 : i_last - start;
        }
     }

     size_t to_be_increased = _encode_symbol(kmer[k-2]);//_encode_symbol(edge_label(start)[k-2]);
     size_t ins = 2*x;
     int outdeg = -1;
     if (ins > (*p_edges)[start]){
          outdeg = outdegree (_edge_to_node(start));
          if (outdeg == 0) { //outgoing dummy
            p_edges->remove(start);
            p_edges->insert(start,ins);
          }
          else {
              size_t i = 1;
              for (i = 1; i < outdeg; ++i){
                  if (ins < (*p_edges)[start+i])
                      break;
	      }
              start += i;
              p_edges->insert(start,ins);
              p_node_flags->insert(start,1);
          }

     }
     else{
	if ((*p_edges)[start] != ins){
	   p_edges->insert(start,ins);
	   p_node_flags->insert(start+1,1);
	}
     }


     for (symbol_type i = 0; i < sigma+1; ++i){
	if (outdeg != 0 && i >= to_be_increased){
	   m_symbol_ends[i]++;
	}
     }


     int d_insertion = -1;

     if (!last_node_present) d_insertion = _forward(start);



     if (d_insertion != -1){ //destination node is not present
	p_edges->insert(d_insertion , 0);
	p_node_flags->insert(d_insertion , 0);
	m_num_nodes++;

     } //make the indegree > 1
     else {
        p_edges->remove(start+check);
        p_edges->insert(start+check,ins+1);
     }
     if (start >= d_insertion)
	start++;
     for (symbol_type i = 0; i < sigma+1; ++i){

	if (d_insertion != -1 && i >= x){
	   m_symbol_ends[i]++;
	}
     }

     return start;
  }

   template <class InputIterator>
   bool index(InputIterator in , bool edge)  {
         auto c = *in++;
         symbol_type first_symbol = _encode_symbol(c);

         // Range is from first edge of first, to last edge of last
         size_t start = _symbol_start(first_symbol);
         size_t end   = m_symbol_ends[first_symbol]-1;
         size_t first = 0, last = 0;

         // find c labeled pred edge
         // if outside of range, find c- labeled pred edge

         for (size_t i = 0; i < k-2; i++) {
   	       c = *in++;

   	       symbol_type x = _encode_symbol(c);
   	 // update range; Within current range, find first and last occurence of c or c-
   	 // first -> succ(x, first)
   	 for (uint8_t y=x<<1; y<(x<<1)+1; y++) {
   	    size_t tmpRank = p_edges->rank(start, y);
   	    if (tmpRank < p_edges->rank(p_edges->size(), y))
   	       first = p_edges->select(tmpRank, y);

   	    if (start <= first && first <= end) break;
   	 }
      if (!(start <= first && first <= end))
          return false;
   	 // last -> pred(x, last)
   	 if (start == end) {
   	    last = first;
   	 } else {
   	    for (uint8_t y=x<<1; y<(x<<1)+1; y++) {
   	       auto rank_temp = p_edges->rank(end + 1, y);
   	       //last = p_edges->select((rank_temp) - 1, y);
   	       last = p_edges->select((rank_temp) - 1, y);
   	       if (start <= last && last <= end) break;
   	    }
   	 }
   	 if (!(start <= last && last <= end)) {
   	    assert(!"(start <= last && last <= end)");
   	 }

   	 // Follow each edge forward
   	 start = _forward(first, x);

   	 end   = _forward(last, x);
   	 end   = _last_edge_of_node(_edge_to_node(end));
 
         }

         return true;
      }


   template <class InputIterator>

   pair<size_t, size_t> index_finder_for_add_delete(InputIterator in , bool addition)  {

       auto ref = in;
       auto c = 0;
       c = (addition) ?  0 : *in++;
       symbol_type first_symbol = _encode_symbol(c);

          // Range is from first edge of first, to last edge of last
        size_t start = _symbol_start(first_symbol);
        size_t end   = m_symbol_ends[first_symbol]-1;
        size_t first=0, last = 0;

      // find c labeled pred edge
      // if outside of range, find c- labeled pred edge
      for (size_t i = 0; i < k-1 ; i++) { //I changed k (addition) to k-1 for deletion
	       c = *in++;
   	     symbol_type x = _encode_symbol(c);


	 // update range; Within current range, find first and last occurence of c or c-
	 // first -> succ(x, first)
	 for (uint8_t y=x<<1; y<(x<<1)+1; y++) {

	    size_t tmpRank = p_edges->rank(start, y);
	    if (tmpRank < p_edges->rank(p_edges->size(), y))
	       first = p_edges->select(tmpRank, y);
         if (start <= first && first <= end)  break;
        }
          if (!(start <= first && first <= end))
                return (std::make_pair(start,i));


	 if (start == end) {
	    last = first;
	 } else {
	    for (uint8_t y=x<<1; y<(x<<1)+1; y++) {
	       auto rank_temp = p_edges->rank(end + 1, y);
	       //last = p_edges->select((rank_temp) - 1, y);
	       last = p_edges->select((rank_temp) - 1, y);
	       if (start <= last && last <= end) break;
	    }
	 }
	 if (!(start <= last && last <= end)) {
	    assert(!"(start <= last && last <= end)");
	 }
	 start = _forward(first, x);
	 end   = _forward(last, x);
   end   = _last_edge_of_node(_edge_to_node(end));

        }

      return (std::make_pair(first,k)); //edge is present
   }


   size_t bit_size() {
      size_t dSize = 0;
      dSize += p_node_flags->bit_size();
      if (p_edges != NULL) {
	 dSize += p_edges->bit_size();
      }

      //	 dSize += size_in_bytes( m_symbol_ends ) * 8;
      //	 dSize += size_in_bytes( m_edge_max_ranks ) * 8;
      //	 dSize += sizeof( m_num_nodes ) * 8;
      return dSize;
   }
};

#endif
