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

size_t dummyArg;

void rev_complement_string ( std::string& in, std::string& out ) {
   out = in;
   size_t n = in.size();
   for (size_t i = 0 ; i < n; ++i) {
      char c;
      switch( in[ n - 1 - i ] ) {
      case 'A':
	 c = 'T';
	 break;
      case 'C':
	 c = 'G';
	 break;
      case 'G':
	 c = 'C';
	 break;
      case 'T':
	 c = 'A';
	 break;
      }

      out[i] = c;
   }
}

class dyn_boss {
public:

   const static size_t sigma = 4;
   typedef size_t edge_type;
   typedef pair<edge_type, edge_type> node_type;
   typedef dyn::wt_str::char_type symbol_type;

   size_t           k;
   //dyn::gap_bv* p_node_flags;
   dyn::suc_bv* p_node_flags;
   dyn::wt_str* p_edges;
   //dyn::str_check* p_edges;

   array<size_t, 1 + sigma> m_symbol_ends;
   array<size_t, 1 + sigma> m_edge_max_ranks;
   string             m_alphabet;
   size_t             m_num_nodes;

   dyn_boss() {
      p_node_flags = NULL;
      p_edges = NULL;
      m_alphabet = "";
      m_num_nodes = 0;
      k = 0;
   }

   dyn_boss( size_t in_k ) {
      p_node_flags = new dyn::suc_bv();
      p_edges = new dyn::wt_str( 10 );
      m_alphabet = "$ACGT";
      m_num_nodes = 0;
      k = in_k;
      for (size_t i = 0; i < sigma + 1; ++i) {
	 m_symbol_ends[i] = 0;
      }
   }

   dyn_boss( size_t in_k,
	     dyn::suc_bv* node_flags,
	     dyn::wt_str* edges,
	     array<size_t, 1 + sigma>& symbol_ends,
	     string alphabet) :
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

      size_t num_blocks = size_t(size)/sizeof(uint64_t) - (sigma+2);
      input.seekg(0, ios::beg); // rewind

      vector<uint64_t> blocks(num_blocks,0);
      input.read((char*)&blocks[0], sizeof(uint64_t) * num_blocks);

      // So we avoid a huge malloc if someone gives us a bad file
      p_node_flags = new dyn::suc_bv();

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
      //if (i_last == num_edges()) i_last = num_edges() - 1; //Bahar
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

   void print_boss_matrix( std::ostream& os ) {
      for (size_t i = 0; i < num_edges(); ++i) {
	 os << i << ' ' << edge_label(i) << ' ';
	 //os << static_cast<char>(_map_symbol(_strip_edge_flag( p_edges->at(i) ))) << ' ';
	 //os << static_cast<char>(_map_symbol(( p_edges->at(i) ))) << ' ';
	 // switch (p_edges->at(i)) {
	 // default:
	 os << p_edges->at(i) << ' ';

	 // }

	 os << p_node_flags->at(i);

	 if (p_node_flags->at(i) == 0) {
	    size_t j = _edge_to_node( i );
	     os << ' ' << node_label( j ) << ' '
	        << ' ' << indegree( j ) << ' ';
	  }

	 os << '\n';
      }

   }

   size_t _symbol_start(symbol_type x) const {
      assert(x < sigma + 1);
      return (x==0)? 0 : m_symbol_ends[x - 1];
   }


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
      //if (i >= num_edges() - 1 ) return i;
      // Might not actually occur if out of rank bounds?
      //size_t next_rank = 1 + p_edges->rank(1+i, _with_edge_flag(x, false));
      size_t next_rank =  p_edges->rank(1+i, _with_edge_flag(x, false)); //Bahar

      //we know next_rank >= 1, since i is ostensibly
      //an edge of symbol_type x
      //if (next_rank >= m_edge_max_ranks[x]) return num_edges();
      if (next_rank >= p_edges->rank( num_edges(), _with_edge_flag(x, false) ) )
	 return num_edges();
      return p_edges->select(next_rank, _with_edge_flag(x, false));
   }

   size_t _first_edge(size_t i, symbol_type x) const {
      if (i >= num_edges()) {
	 return i;
      }

      if ( p_edges->at( i ) == x ) {
	 return i;
      }

      size_t next_rank =  p_edges->rank(1+i, x );

      //we know next_rank >= 1, since i is ostensibly
      //an edge of symbol_type x
      //if (next_rank >= m_edge_max_ranks[x]) return num_edges();
      size_t rvalue;
      if (next_rank < ( p_edges->rank( num_edges(), x ) ) )
	 rvalue = p_edges->select(next_rank, x );
      else
	 rvalue = num_edges();

      return rvalue;
   }

   size_t _rank_distance(size_t a, size_t b) {
      return p_node_flags->rank0(b) - p_node_flags->rank0(a);
   }

   size_t _rank_distance_alan(size_t a, size_t b) {
      if (a+1 > num_edges())
	 a = num_edges() - 1;
      if (b+1 > num_edges())
	 b = num_edges() - 1;

      return p_node_flags->rank0(b + 1) - p_node_flags->rank0(a + 1);
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


   void load(std::istream& in) {
     read_member(deconst(k), in);
     p_node_flags = new dyn::suc_bv();
     deconst(p_node_flags)->load(in);
     //deconst(m_node_rank).load(in);
     //deconst(m_node_rank).set_vector(&m_node_flags);
     //deconst(m_node_select).load(in);
     //deconst(m_node_select).set_vector(&m_node_flags);
     p_edges = new dyn::wt_str( 10 );
     deconst(p_edges)->load(in);
     read_member(deconst(m_symbol_ends), in);
     read_member(deconst(m_edge_max_ranks), in);
     read_member(deconst(m_alphabet), in);
     read_member(deconst(m_num_nodes), in);

  //size_type size() const { return num_edges(); }
};

   symbol_type _encode_symbol(uint8_t c) const {
      return lower_bound(m_alphabet.begin(), m_alphabet.end(), c) - m_alphabet.begin();
   }

   //This function deletes an edge kmer at position pos
   //
   // It may leave the graph in an inconsistent state,
   // (all nodes may not have incoming edge, etc.),
   // so use carefully and maintain the graph state externally.
   void _delete_edge_at_position( size_t pos, string kmer, bool delete_if_empty = true ) {
      //First, check if this is first edge of multi-edge node
      //If so, the node flag of the next edge must be updated
      size_t node_decrement = 1; //did we remove a node?
      if (pos < num_edges()) {
	 if (p_node_flags->at( pos ) == 0) { //first of node
	    if (pos < num_edges() - 1) {
	       if (p_node_flags->at(pos + 1)) { //second belongs to same node
		  p_node_flags->remove( pos + 1 );
		  p_node_flags->insert( pos + 1, 0 ); //second is now first after deletion at pos
		  node_decrement = 0;
	       } else {
		  //this edge is the only edge of this node
		  //hence, after deletion, this node will be empty
		  if (!delete_if_empty) {
		     //the user doesn't want to delete the node, just the edge
		     node_decrement = 0;
		     _add_edge_to_node( pos, 0 );
		     //this is now a multi-edge node (with an empty outgoing).
		     //the edge to be deleted is at pos + 1
		     pos = pos + 1;
		  }
	       }
	    } else {
	       //this edge is the only edge of this node (and is the last edge in the graph)
	       //hence, after deletion, this node will be empty
	       if (!delete_if_empty) {
		  //the user doesn't want to delete the node, just the edge
		  node_decrement = 0;
		  _add_edge_to_node( pos, 0 );
		  //this is now a multi-edge node (with an empty outgoing).
		  //the edge to be deleted is at pos + 1
		  pos = pos + 1;
	       }
	    }
	 } else {
	    //not first of node
	    node_decrement = 0;
	 }
      } else {
	 //not a valid edge position
	 return;
      }

      //Next, check if there is another edge incoming to the
      //same target of the edge at pos.
      //If this edge has a minus flag, may ignore, as the flag of this edge
      //will not affect others
      bool edge_flag = p_edges->at( pos ) & 1;
      if (edge_flag) {
	 //simply remove edge and node flags, and we're done
	 //removal below

      } else {
	 //need to check if any other edges are incoming to the same target
	 size_t rankPos = p_edges->rank( pos + 1, p_edges->at( pos ) );
	 size_t rankSym = p_edges->rank( num_edges(), p_edges->at( pos ) );
	 symbol_type xminus = p_edges->at( pos ) | 1;
	 if ( rankPos < rankSym ) {
	    //need to check for minuses between this occurrence and the next
	    size_t nextPos = p_edges->select( rankPos, p_edges->at( pos ) );

	    size_t rankMinusPos = p_edges->rank( pos + 1, xminus );
	    size_t rankMinusNextPos = p_edges->rank( nextPos + 1, xminus );
	    if (rankMinusNextPos > rankMinusPos) {
	       //sure enough, we've found one. Let's change the first one
	       //to be not minused.
	       size_t posMinus = p_edges->select( rankMinusPos, xminus );
	       p_edges->remove( posMinus );
	       p_edges->insert( posMinus, p_edges->at( pos ) );
	    }
	 } else {
	    //need to check for minuses after this occurrence
	    size_t rankMinusPos = p_edges->rank( pos + 1, xminus );
	    size_t rankMinusNextPos = p_edges->rank( num_edges(), xminus );
	    if (rankMinusNextPos > rankMinusPos) {
	       //sure enough, we've found one. Let's change the first one
	       //to be not minused.
	       size_t posMinus = p_edges->select( rankMinusPos, xminus );
	       p_edges->remove( posMinus );
	       p_edges->insert( posMinus, p_edges->at( pos ) );
	    }
	 }
	 //Minuses should be handled now, we can safely remove and return

      }

      symbol_type val =  _encode_symbol(kmer[k-2]);//edge_label(start)[k-2]);p_edges->at( pos ) >> 1;

      p_edges->remove( pos );
      p_node_flags->remove( pos );

      //Need to update the number of characters and nodes
      m_num_nodes -= node_decrement;
      for (symbol_type i = 0; i < sigma+1; ++i){
	 if ( i >= val){
	    --m_symbol_ends[i];
	 }
      }
   }

   size_t _delete_edge_from_node( size_t pos, string kmer, bool delete_if_empty = false ) {
      size_t nodeId = _edge_to_node( pos );

      size_t posTarget = _forward( pos );
      size_t inDegTarget = indegree( _edge_to_node( posTarget ) );
      size_t outDegTarget = outdegree( _edge_to_node( posTarget ) );

      if (inDegTarget == 1 && outDegTarget == 0) {
	 //delete target
	 p_node_flags->remove( posTarget );
	 p_edges->remove( posTarget );
	 if (posTarget < pos) {
	    --pos;
	    --nodeId;
	 }
	 symbol_type val = _encode_symbol(kmer[k-1]);
	 for (symbol_type i = 0; i < sigma+1; ++i){
	    if (i >= val)
	       --m_symbol_ends[i];
	 }

	 m_num_nodes -= 1;

      } else {
	 if (inDegTarget == 1 && outDegTarget > 0) {
	    //target needs a dummy chain
	    _add_dummy_chain( kmer.begin() + 1, true );
	    //graph has changed, so reindex the edge to delete
	    index_edge_alan( kmer.begin(), pos );
	    nodeId = _edge_to_node( pos );
	 }
      }

      //next, delete the edge

      //First, check if this is first edge of multi-edge node
      //If so, the node flag of the next edge must be updated
      size_t node_decrement = 1; //did we remove a node?
      if (pos < num_edges()) {
	 if (p_node_flags->at( pos ) == 0) { //first of node
	    if (pos < num_edges() - 1) {
	       if (p_node_flags->at(pos + 1) == 1) { //second belongs to same node
		  p_node_flags->remove( pos + 1 );
		  p_node_flags->insert( pos + 1, 0 ); //second is now first after deletion at pos
		  node_decrement = 0;
	       } else {
		  //this edge is the only edge of this node
		  //hence, after deletion, this node will be empty
		  if (!delete_if_empty) {
		     //the user doesn't want to delete the node, just the edge
		     node_decrement = 0;
		     _add_edge_to_node( pos, 0 );
		     //this is now a multi-edge node (with an empty outgoing).
		     //the edge to be deleted is at pos + 1
		     pos = pos + 1;
		  }
	       }
	    } else {
	       //this edge is the only edge of this node (and is the last edge in the graph)
	       //hence, after deletion, this node will be empty
	       if (!delete_if_empty) {
		  //the user doesn't want to delete the node, just the edge
		  node_decrement = 0;
		  _add_edge_to_node( pos, 0 );
		  //this is now a multi-edge node (with an empty outgoing).
		  //the edge to be deleted is at pos + 1
		  pos = pos + 1;
	       }
	    }
	 } else {
	    //not first of node
	    node_decrement = 0;
	 }
      } else {
	 //not a valid edge position
	 return 0;
      }

      //Next, check if there is another edge incoming to the
      //same target of the edge at pos.
      //If this edge has a minus flag, may ignore, as the flag of this edge
      //will not affect others
      bool edge_flag = p_edges->at( pos ) & 1;
      if (edge_flag) {
	 //simply remove edge and node flags, and we're done
	 //removal below

      } else {
	 //need to check if any other edges are incoming to the same target
	 size_t rankPos = p_edges->rank( pos + 1, p_edges->at( pos ) );
	 size_t rankSym = p_edges->rank( num_edges(), p_edges->at( pos ) );
	 symbol_type xminus = p_edges->at( pos ) | 1;
	 if ( rankPos < rankSym ) {
	    //need to check for minuses between this occurrence and the next
	    size_t nextPos = p_edges->select( rankPos, p_edges->at( pos ) );

	    size_t rankMinusPos = p_edges->rank( pos + 1, xminus );
	    size_t rankMinusNextPos = p_edges->rank( nextPos + 1, xminus );
	    if (rankMinusNextPos > rankMinusPos) {
	       //sure enough, we've found one. Let's change the first one
	       //to be not minused.
	       size_t posMinus = p_edges->select( rankMinusPos, xminus );
	       p_edges->remove( posMinus );
	       p_edges->insert( posMinus, p_edges->at( pos ) );
	    }
	 } else {
	    //need to check for minuses after this occurrence
	    size_t rankMinusPos = p_edges->rank( pos + 1, xminus );
	    size_t rankMinusNextPos = p_edges->rank( num_edges(), xminus );
	    if (rankMinusNextPos > rankMinusPos) {
	       //sure enough, we've found one. Let's change the first one
	       //to be not minused.
	       size_t posMinus = p_edges->select( rankMinusPos, xminus );
	       p_edges->remove( posMinus );
	       p_edges->insert( posMinus, p_edges->at( pos ) );
	    }
	 }
	 //Minuses should be handled now, we can safely remove and return

      }

      symbol_type val =  _encode_symbol(kmer[k-2]);//edge_label(start)[k-2]);p_edges->at( pos ) >> 1;

      p_edges->remove( pos );
      p_node_flags->remove( pos );

      //Need to update the number of characters and nodes
      m_num_nodes -= node_decrement;
      for (symbol_type i = 0; i < sigma+1; ++i){
	 if ( i >= val){
	    --m_symbol_ends[i];
	 }
      }

      return nodeId;
   }

   bool _check_target( size_t pos,
   		       symbol_type nodeSymbol,
   		       symbol_type outSymbol,
   		       size_t& idxIncoming ) {
      //an edge was just inserted at pos, with outgoing symbol outSymbol,
      //into a node with symbol nodeSymbol

      symbol_type& x = outSymbol;
      if ( p_edges->rank( pos, x ) > 0 ) {
   	 //get the preceding edge with x
   	 idxIncoming = p_edges->select( p_edges->rank( pos, x ) - 1, x );
   	 if ( _symbol_access( idxIncoming ) == nodeSymbol ) {


	 }

      }

      if ( p_edges->rank( pos + 1, x ) < p_edges->rank( num_edges(), x )) {
   	 //get the succeeding edge with x
   	 idxIncoming = p_edges->select( p_edges->rank( pos + 1, x ), x );
   	 if ( _symbol_access( idxIncoming ) == nodeSymbol )
   	    return true;
      }

      return false;
   }

   /*
    * @pos -- the first edge of the node to which this edge will be added
    * @x   -- the character of the outgoing edge
    * @targetExists, @idxIncoming
    *      -- if the node corresponding to the target of this edge exists,
    *         targetExists = true, and idxIncoming contains the index of
    *         the first edge incident with target
    */
   size_t _add_edge_to_node( size_t pos,
			     symbol_type x,
			     bool targetExists = false,
			     size_t idxIncoming = 0, size_t* idxTarget = NULL ) {
      bool edge_added = false;

      //we need the character of this node
      symbol_type y = _symbol_access( pos );

      if ( p_edges->at(pos) == 0 ) {
	 p_edges->remove(pos);
	 p_edges->insert(pos, x);
	 //for (symbol_type i = 0; i < sigma+1; ++i){
	 //--m_symbol_ends[i];
	 //}

      } else {

	 //Find where within this node to insert the edge
	 size_t startPos = pos;
	 bool nextNode = false;
	 while (pos < num_edges()) {
	    if (p_edges->at( pos ) < x) {
	       ++pos;
	       if (pos >= num_edges()) {
		  nextNode = true;
		  break;
	       }
	       if (p_node_flags->at(pos) == 0) {
		  nextNode = true;
		  break;
	       }
	    } else {
	       break;
	    }


	 }

	 if (!nextNode) {
	    if ( (p_edges->at(pos) >> 1) == (x >> 1) ) {
	       //this edge already exists
	       return pos;
	    }
	 }

	 edge_added = true;

	 p_edges->insert( pos, x );

	 if (pos > startPos)
	    p_node_flags->insert( pos, 1 );
	 else
	    p_node_flags->insert( pos + 1, 1 );


	 for (symbol_type i = 0; i < sigma+1; ++i){
	    if ( i >= (y)){
	       ++m_symbol_ends[i];
	    }
	 }
      }

      if (x == 0) {
	 //an empty edge was added to this node. let's hope
	 //the user knows what they are doing
	 //no need to worry about targets here.
	 return pos;
      }

      //need to check if target node exists. If not, need to add an
      //empty edge to create it.
      //      size_t idxIncomingTest;
      //bool targetExistsTest = _check_target( pos,
      //y,
      //x,
      //idxIncomingTest );

      if (targetExists) {
	 if (idxTarget != NULL) {
	    if (*idxTarget >= pos) {
	       if (edge_added)
		  ++(*idxTarget);
	    }
	 }

	 //need to determine if the newly added edge should have a minus flag.
	 //if not, the idxIncoming needs to have a minux flag
	 //this is decided simply by their relative order
	 if (idxIncoming >= pos) {
	    if (edge_added)
	       ++idxIncoming; //a new edge has been added before it
	    //and it needs to be replaced with a minus flag
	    p_edges->remove( idxIncoming );
	    p_edges->insert( idxIncoming, x | 1 ); //has the same symbol as the new edge at pos
	 } else {

	    p_edges->remove( pos );
	    p_edges->insert( pos, x | 1 );

	 }
      } else {
	 //target node does not exist
	 //the position it should be at:
	 size_t nodeTarget = p_node_flags->rank0( _symbol_start( x>>1 ) ) + p_edges->rank( pos + 1, x );

	 size_t posTarget = 0; // _symbol_start( x >> 1 ) + p_edges->rank( pos, x );

	 if (nodeTarget > m_num_nodes) {
	    posTarget = p_node_flags->size();
	 } else {
	    posTarget = p_node_flags->select0( nodeTarget - 1);
	 }

	 p_node_flags->insert( posTarget , 0 );
	 m_num_nodes += 1;
	 p_edges->insert( posTarget, 0 ); // dummy edge

	 for (symbol_type i = 0; i < sigma+1; ++i){
	    if ( i >= (x>>1)){
	       ++m_symbol_ends[i];
	    }
	 }

	 if (posTarget <= pos)
	    ++pos;
      }

      return pos;
   }

   /*
    * Adds dummy chain for the node represented by
    * label starting at nodeBegin
    * returns the index of first edge of nodeBegin
    */
   size_t _add_dummy_chain( std::string::iterator nodeBegin, bool targetExists = false ) {

      char c;
      size_t nodeStart = 0;
      size_t outEdge;
      size_t pos = 0;
      dyn_boss::symbol_type x;
      std::string::iterator nodeBeginSaved = nodeBegin;

      //need to make sure the all dummy node exists
      if (m_symbol_ends[ 0 ] == 0) {
	 //$$$$ does not exist.
	 //create it
	 //it is needed
	 p_node_flags->insert( 0, 0 );
	 p_edges->insert( 0, 0 );
	 m_num_nodes += 1;
	 for (dyn_boss::symbol_type i = 0; i < sigma+1; ++i){
	    ++(m_symbol_ends[i]);
	 }
      }

      bool outEdgeExists = true; //optimistic, aren't we?
      do {
	 c = *nodeBegin++;
	 ++pos;
	 x = _encode_symbol( c ) << 1;

	 outEdge = _first_edge( nodeStart, x );

	 size_t dist = _rank_distance_alan( nodeStart, outEdge );

	 if (dist > 0 || outEdge == num_edges() ) {
	    outEdgeExists = false;
	 } else {
	    //the edge exists, so let's follow it
	    nodeStart = _forward( outEdge );

	 }
      } while (outEdgeExists);

      //we have gotten to a dummy node that needs a new outgoing edge
      while (pos < k) {
	 //string dedge( k - pos, '$' );
	 //dedge += kmer.substr( 0, pos );
	 if (pos < k - 1)
	    outEdge = _add_edge_to_node( nodeStart, x );
	 else {
	    if (targetExists) {
	       size_t idxIncoming;
	       index_alan( nodeBeginSaved, idxIncoming );
	       idxIncoming = _node_to_edge( idxIncoming );
	       idxIncoming = _backward( idxIncoming );
	       outEdge = _add_edge_to_node( nodeStart, x, targetExists, idxIncoming );
	    } else {
	       outEdge = _add_edge_to_node( nodeStart, x );
	    }
	 }

	 nodeStart = _forward( outEdge );

	 ++pos;
	 auto c = *nodeBegin++;
	 x = _encode_symbol( c ) << 1;
      }

      return nodeStart;
   }

   // /*
   //  * removes the dummy chain ending at dummyLabel
   //  */
   bool _remove_dummy_chain( std::string dummyLabel, bool removeTarget = false, size_t* targetIdPtr = NULL ) {
      size_t idx;
      bool rvalue = false;
      size_t targetId;
      if (index_edge_alan( dummyLabel.begin(), idx ) ) {
   	 rvalue = true;
   	 vector< size_t > vDelete;
   	 vector< string > vEdgesToDelete;
   	 targetId = _forward( idx );

   	 if (removeTarget) {
   	    string targetLabel = dummyLabel.substr(1) + "$";
   	    vEdgesToDelete.push_back( targetLabel );
   	    vDelete.push_back( targetId );
   	 }

   	 size_t start = idx;
   	 size_t i = 0;
   	 while ( i < k - 1){

   	    if (start <= targetId) {
   	       --targetId;
   	    }

   	    size_t out = outdegree(_edge_to_node(start));

   	    //insert to maintain increasing order
   	    vector<size_t>::iterator itPos = vDelete.begin();
   	    size_t posDelete = 0;
   	    while ( itPos != vDelete.end() ) {
   	       if ( (*itPos) < start ) {
   		  ++itPos;
   		  ++posDelete;
   	       } else {
   		  break;
   	       }
   	    }

   	    vDelete.insert( itPos, start );
   	    vEdgesToDelete.insert( (vEdgesToDelete.begin() + posDelete),
   				   dummyLabel );

   	    if (i == k - 1)
   	       break;

   	    if (out > 1 ) break;
   	    size_t backward_edge;
   	    start = _backward(start);

   	    dummyLabel = '$' + dummyLabel.substr(0,k-1);
   	    i++;

   	 }


   	 for (size_t i = 0; i < vDelete.size(); ++i) {
   	    _delete_edge_at_position( vDelete[i] - i, vEdgesToDelete[i] );

   	 }

   	 if (targetIdPtr != NULL) {
   	    *targetIdPtr = targetId;
   	 }
      }

      return rvalue;
   }

   /*
    * removes the dummy chain ending at start (inclusive)
    */
   void _remove_dummy_chain( size_t start, std::string dummyLabel, bool removeTarget = false ) {

      vector< size_t > vDelete;
      vector< string > vEdgesToDelete;

      if (removeTarget) {
	 size_t targetId = _forward( start );
	 string targetLabel = dummyLabel.substr(1) + "$";
	 vEdgesToDelete.push_back( targetLabel );
	 vDelete.push_back( targetId );
      }

      size_t i = 0;
      while ( i < k - 1){

	 size_t out = outdegree(_edge_to_node(start));

	 //insert to maintain increasing order
	 vector<size_t>::iterator itPos = vDelete.begin();
	 size_t posDelete = 0;
	 while ( itPos != vDelete.end() ) {
	    if ( (*itPos) < start ) {
	       ++itPos;
	       ++posDelete;
	    } else {
	       break;
	    }
	 }

	 vDelete.insert( itPos, start );
	 vEdgesToDelete.insert( (vEdgesToDelete.begin() + posDelete),
				dummyLabel );

	 if (i == k - 1)
	    break;

	 if (out > 1 ) break;
	 size_t backward_edge;
	 start = _backward(start);

	 dummyLabel = '$' + dummyLabel.substr(0,k-1);
	 i++;

      }


      for (size_t i = 0; i < vDelete.size(); ++i) {
	 _delete_edge_at_position( vDelete[i] - i, vEdgesToDelete[i] );

      }


   }

   /*
    * removes the dummy chain ending at start (inclusive)
    */
   void _remove_dummy_chain_alt( size_t start, std::string dummyLabel ) {

      size_t i = 0;
      while ( i < k - 1){

	 size_t nodeStart = _delete_edge_from_node( start, dummyLabel, false );

	 size_t out = outdegree(nodeStart);

	 if (out > 0 ) break;

	 start = _first_edge_of_node( nodeStart );
	 size_t backward_edge;
	 start = _backward(start);

	 dummyLabel = '$' + dummyLabel.substr(0,k-1);
	 i++;

      }
      //check if there is only one dummy
      //node remaining. It is safe to remove it.
      if (p_edges->at(0) == 0) {
	 dummyLabel.assign(k, '$');
	 _delete_edge_at_position( 0, dummyLabel, true );
      }
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

   //input: @in node string to be tested for membership
   //output: @idx = index of node if the returned bool is true
   //        @rvalue = true iff. node exists in dbg
   bool index_alan(std::string::iterator in, size_t& idx = dummyArg ) {
      auto c = *in++;
      symbol_type first_symbol = _encode_symbol(c);

      // Range is from first edge of first, to last edge of last
      size_t start, end;
      size_t rankStart, rankEnd;
      dyn_boss::symbol_type x;
      for (size_t pos = 0; pos < k - 1; ++pos) {
	 if (pos == 0) {
	    start = _symbol_start(first_symbol);
	    end = static_cast<int64_t>( m_symbol_ends[first_symbol]);
	 } else {
	    //move forward
	    size_t first = p_edges->select( rankStart, x );
	    size_t last = p_edges->select( rankEnd - 1, x );

	    start = _forward(first, x);
	    end   = _forward(last, x);
	    end   = _last_edge_of_node(_edge_to_node(end)) + 1;
	 }

	 if (pos == k - 2)
	    break;

	 c = *in++;


	 x = _encode_symbol(c) << 1;
	 rankEnd = p_edges->rank( end, x );
	 rankStart = p_edges->rank( start, x );

	 if (rankEnd == rankStart) {
	    //no outgoing edges with x
	    return false;
	 }
      }


      idx = _edge_to_node( start );

      return true;
   }

   //input: @in edge string to be tested for membership
   //output: @idx = index of edge if the returned bool is true
   //        @rvalue = true iff. edge  exists in dbg
   bool index_edge_alan(std::string::iterator in, size_t& idx = dummyArg) {

      auto c = *in++;
      symbol_type first_symbol = _encode_symbol(c);
      // Range is from first edge of first, to last edge of last
      size_t start, end;
      size_t rankStart, rankEnd;
      dyn_boss::symbol_type x;
      for (size_t pos = 0; pos < k - 1; ++pos) {
   	 if (pos == 0) {
   	    start = _symbol_start(first_symbol);
   	    end = static_cast<int64_t>( m_symbol_ends[first_symbol]);
   	 } else {
   	    //move forward
   	    size_t first = p_edges->select( rankStart, x );
   	    size_t last = p_edges->select( rankEnd - 1, x );

   	    start = _forward(first, x);
   	    end   = _forward(last, x);
   	    end   = _last_edge_of_node(_edge_to_node(end)) + 1;
   	 }

   	 if (pos == k - 2)
   	    break;

	 c = *in++;


   	 x = _encode_symbol(c) << 1;
   	 rankEnd = p_edges->rank( end, x );
   	 rankStart = p_edges->rank( start, x );

   	 if (rankEnd == rankStart) {
   	    //no outgoing edges with x
   	    return false;
   	 }
      }

      //we are at the node that must contain this edge
      c = *in++;
      x = _encode_symbol(c) << 1;
      rankEnd = p_edges->rank( end, x );
      rankStart = p_edges->rank( start, x );
      if (rankEnd == rankStart) {
	 //no outgoing edge with x
	 //try xminus
	 rankEnd = p_edges->rank( end, x | 1 );
	 rankStart = p_edges->rank( start, x | 1 );
	 if (rankEnd == rankStart) {
	    //no outgoing edge with xminus
	    //this edge does not exist
	    return false;
	 } else {
	    idx = p_edges->select( rankStart, x | 1);
	 }
      } else {
	 idx = p_edges->select( rankStart, x );
      }



      return true;
   }


   template <class InputIterator>
   pair<size_t, size_t> index_finder(InputIterator in, string kmer)  {

      auto ref = in;
      auto c = *in++;
      bool new_node = false;
      symbol_type first_symbol = _encode_symbol(c);

      // Range is from first edge of first, to last edge of last
      size_t start = _symbol_start(first_symbol);
      if (m_symbol_ends[first_symbol] == 0) {
	 //we already need to add, this edge is not present
	 return std::make_pair( 0, 0 );
      }

      size_t end   = m_symbol_ends[first_symbol]-1;

      size_t first=0, last = 0;

      // find c labeled pred edge
      // if outside of range, find c- labeled pred edge
      for (size_t i = 0; i < k - 1 ; i++) { //I changed k (addition) to k-1 for deletion
	 c = *in++;
	 symbol_type x = _encode_symbol(c);

	 // update range; Within current range, find first and last occurence of c or c-
	 // first -> succ(x, first)
	 for (uint8_t y=x<<1; y<=(x<<1)+1; y++) {

	    size_t tmpRank = p_edges->rank(start, y);
	    if (tmpRank < p_edges->rank(p_edges->size(), y))
	       first = p_edges->select(tmpRank, y);

	    if (start <= first && first <= end)  break;
	 }
	 if (!(start <= first && first <= end)){
            if (new_node) return (std::make_pair(start,i));
            else return (std::make_pair(start,i+1));
	 }



	 if (start == end) {
	    last = first;
	 } else {
	    for (uint8_t y=x<<1; y <= (x<<1)+1; ++y) {
	       auto rank_temp = p_edges->rank(end + 1, y);
	       //last = p_edges->select((rank_temp) - 1, y);
	       last = p_edges->select((rank_temp) - 1, y);
	       if (start <= last && last <= end) break;
	    }
	 }


	 //assert(!(start <= last && last <= end));

	 start = _forward(first, x);
	 end   = _forward(last, x);
	 end   = _last_edge_of_node(_edge_to_node(end));

      }

      return (std::make_pair(first,k)); //edge is present
   }

   template <class InputIterator>
   pair<size_t, size_t> index_finder_for_add_delete(InputIterator in , bool addition, string kmer)  {
      auto ref = in;
      auto c = 0;
      bool new_node = false;
      if (addition && index(kmer.substr(0,k-1).begin() ,0) != 1)
	 new_node = true;
      else c = *in++;

      symbol_type first_symbol = _encode_symbol(c);

      // Range is from first edge of first, to last edge of last
      size_t start = _symbol_start(first_symbol);
      if (m_symbol_ends[first_symbol] == 0) {
	 if (new_node) return (std::make_pair(start, 0));
	 else return (std::make_pair(start,1));
      }

      size_t end   = m_symbol_ends[first_symbol]-1;

      size_t first=0, last = 0;

      // find c labeled pred edge
      // if outside of range, find c- labeled pred edge
      for (size_t i = 0; i < k - 1 ; i++) { //I changed k (addition) to k-1 for deletion
	 c = *in++;
	 symbol_type x = _encode_symbol(c);

	 // update range; Within current range, find first and last occurence of c or c-
	 // first -> succ(x, first)
	 for (uint8_t y=x<<1; y<=(x<<1)+1; y++) {

	    size_t tmpRank = p_edges->rank(start, y);
	    if (tmpRank < p_edges->rank(p_edges->size(), y))
	       first = p_edges->select(tmpRank, y);

	    if (start <= first && first <= end)  break;
	 }
	 if (!(start <= first && first <= end)){
            if (new_node) return (std::make_pair(start,i));
            else return (std::make_pair(start,i+1));
	 }



	 if (start == end) {
	    last = first;
	 } else {
	    for (uint8_t y=x<<1; y <= (x<<1)+1; ++y) {
	       auto rank_temp = p_edges->rank(end + 1, y);
	       //last = p_edges->select((rank_temp) - 1, y);
	       last = p_edges->select((rank_temp) - 1, y);
	       if (start <= last && last <= end) break;
	    }
	 }


	 //assert(!(start <= last && last <= end));

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


   void add_read( std::vector< std::string >& kmers ) {
      add_edge( kmers[0], false );
      size_t idxSource;
      index_alan( kmers[0].begin() + 1, idxSource );
      for (size_t i = 0; i < kmers.size(); ++i) {
	 idxSource = add_edge_given_position( kmers[i], idxSource );
      }
   }

   //@rvalue: idx of target node of edge after addition
   size_t add_edge_given_position( std::string kmer,
				   size_t idxSource ) {

      size_t idxTarget;
      size_t indegTarget;
      bool targetExists = false;
      string::iterator in = kmer.begin();
      size_t edgeStart = _node_to_edge(idxSource);
      size_t outEdge;
      char c;

      //The source exists and is located at edgeStart
      dyn_boss::symbol_type x = _encode_symbol( kmer.back() ) << 1;
      outEdge = _first_edge( edgeStart, x );
      size_t dist = _rank_distance_alan( edgeStart, outEdge );
      bool edge_added = false;
      if (dist > 0 || outEdge == num_edges() ) {
	 //hmm, let's try again with the minus flag

	 outEdge = _first_edge( edgeStart, x | 1);

	 size_t dist = _rank_distance_alan( edgeStart, outEdge );
	 if (dist > 0 || outEdge == num_edges() ) {

	    //the edge does not exist
	    size_t idxIncoming;
	    targetExists = index_alan( kmer.begin() + 1, idxTarget );

	    if (targetExists) {


	       idxTarget = _node_to_edge( idxTarget );
	       idxIncoming = _backward( idxTarget  );

	    }

	    _add_edge_to_node( edgeStart, x, targetExists, idxIncoming, &idxTarget );
	    edge_added = true;
	 } else {
	    //the edge already exists :/
	    idxTarget = _forward( outEdge );
	 }
      } else {
	 //the edge already exists :/
	 idxTarget = _forward( outEdge ) ;
      }

      // DANGER WILL ROBINSON
      // I have commented this out so idxTarget will be correct.
      // However, dummies are no longer taken care of, which
      // is ill-advised.
      //If target has an incoming dummy, it is no longer necessary
      if (targetExists && edge_added) {
	 indegTarget = indegree( _edge_to_node( idxTarget ) );
      	 if (indegTarget == 1) { //TODO: check to see if indegree is really working.
      	    string incoming_dummy = '$'+ kmer.substr(1); //dummy edge
      	    //remove dummy chain first checks if the dummy exists
      	    _remove_dummy_chain( incoming_dummy, false, &idxTarget );
      	 }
      }

      return _edge_to_node(idxTarget);

   }

   void add_edge (std::string kmer, bool add_reverse_complement = false ) {
      //First, check if the source, target nodes exist
      size_t idxTarget;
      size_t idxSource;
      size_t indegTarget;
      bool sourceExists = index_alan( kmer.begin(), idxSource );
      bool targetExists = false;
      string::iterator in = kmer.begin();
      size_t nodeStart = 0;
      size_t outEdge;
      char c;

      if (!sourceExists) {
	 //add dummies + source
	 nodeStart = _add_dummy_chain( kmer.begin() );
      } else {
	 nodeStart = _node_to_edge(idxSource);
      }

      //The source exists and is located at nodeStart
      dyn_boss::symbol_type x = _encode_symbol( kmer.back() ) << 1;
      outEdge = _first_edge( nodeStart, x );
      size_t dist = _rank_distance_alan( nodeStart, outEdge );
      bool edge_added = false;
      if (dist > 0 || outEdge == num_edges() ) {
	 //hmm, let's try again with the minus flag

	 outEdge = _first_edge( nodeStart, x | 1);

	 size_t dist = _rank_distance_alan( nodeStart, outEdge );
	 if (dist > 0 || outEdge == num_edges() ) {

	    //the edge does not exist
	    size_t idxIncoming;
	    targetExists = index_alan( kmer.begin() + 1, idxTarget );

	    if (targetExists) {

	       indegTarget = indegree( idxTarget );
	       idxTarget = _node_to_edge( idxTarget );
	       idxIncoming = _backward( idxTarget );

	    }

	    _add_edge_to_node( nodeStart, x, targetExists, idxIncoming );
	    edge_added = true;
	 } else {
	    //the edge already exists :/

	 }
      } else {
	 //the edge already exists :/
      }

      //If target has an incoming dummy, it is no longer necessary
      if (targetExists && edge_added) {
	 if (indegTarget == 1) { //TODO: check to see if indegree is really working.
	    string incoming_dummy = '$'+ kmer.substr(1); //dummy edge
	    //remove dummy chain first checks if the dummy exists
	    _remove_dummy_chain( incoming_dummy );
	 }
      }

      if (add_reverse_complement) {
	 string revcomp;
	 rev_complement_string( kmer, revcomp );
	 add_edge( revcomp, false );
      }
   }

   void delete_edge( std::string kmer,
		     bool delete_reverse_complement = false ) {

      size_t indexEdge;
      if ( index_edge_alan( kmer.begin(), indexEdge ) ) {
	 //the edge exists

	 size_t u = _edge_to_node( indexEdge );
	 size_t outdeg = outdegree( u );
	 //size_t indeg = indegree( u );

	 // Delete the edge. At this point, do not want to delete source yet
	 _delete_edge_from_node( indexEdge, kmer, false  );

	 if (outdeg == 1) {
	    //need to delete associated dummies, if any
	    string dummyLabel = "$" + kmer;
	    dummyLabel.pop_back();
	    size_t idx;
	    if (index_edge_alan( dummyLabel.begin(), idx ) ) {
	       //remove dummies as required
	       _remove_dummy_chain_alt( idx, dummyLabel );
	    }

	 }
      }

      if ( delete_reverse_complement ) {
	 string revcomp;
	 rev_complement_string( kmer, revcomp );
	 delete_edge( revcomp, false );
      }
   }
};



#endif
