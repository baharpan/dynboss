// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 * spsi_check.hpp
 *
 *  Created on: Oct 21, 2015
 *      Author: nico
 *
 *  trivial implementation of SPSI structure. To be used to test other SPSIs
 *
 */

#ifndef INTERNAL_SPSI_CHECK_HPP_
#define INTERNAL_SPSI_CHECK_HPP_

#include "includes.hpp"

namespace dyn{

template<uint64_t MAX_SIZE=10000000>
class spsi_check{

public:

	spsi_check(){

		vec = vector<uint64_t>(MAX_SIZE);
		size_=0;

	}

	uint64_t& operator[](ulint i){
		return vec.at(i);
	}

	uint64_t at(uint64_t i) const {
		return vec[i];
	}

	uint64_t psum() const {
		return psum(size_-1);
	}

	uint64_t psum(uint64_t i) const {

		uint64_t s = 0;
		for(uint64_t j=0;j<=i;++j){
			s+=vec[j];
		}

		return s;

	}

	uint64_t search(uint64_t x) const {

		assert(size_>0);
		assert(x<=psum());

		uint64_t i = 0;
		uint64_t s = vec[0];

		while(s<x)	s += vec[++i];

		return i;
	}

	uint64_t search_0(uint64_t x) const {

		assert(size_>0);
		assert(x<=size()-psum());

		uint64_t i = 0;
		uint64_t s = not vec[0];

		while(s<x)	s += (not vec[++i]);

		return i;
	}

	uint64_t search_r(uint64_t x) const {

		assert(size_>0);
		assert(x<=psum()+size());

		uint64_t i = 0;
		uint64_t s = vec[0]+1;

		while(s<x)	s += (vec[++i]+1);

		return i;
	}

	bool contains(uint64_t x ) const {

		if(x==0) return true;

		uint64_t s = 0;

		for(uint64_t i=0;i<size();++i){

			s += vec[i];

			if(s==x) return true;
			if(s>x) return false;

		}

		assert(false);
		return false;

	}

	bool contains_r(uint64_t x ) const {

		if(x==0) return true;

		uint64_t s = 0;

		for(uint64_t i=0;i<size();++i){

			s += (vec[i]+1);

			if(s==x) return true;
			if(s>x) return false;

		}

		assert(false);
		return false;

	}

	void increment(uint64_t i, uint64_t delta, bool subtract = false){

		if(subtract){
			vec[i]-=delta;
		}else{
			vec[i]+=delta;
		}

	}

	void insert(uint64_t i, uint64_t x){

		for(uint64_t j=size_;j>i;j--){

			vec[j] = vec[j-1];

		}

		vec[i] = x;

		size_++;

	}

        void remove(uint64_t i){
		size_--;
		for(uint64_t j=i;j<size_;j++){

			vec[j] = vec[j+1];

		}

        }

	ulint serialize(std::ostream &/*out*/) const {
	    throw std::domain_error("not implemented");
	    return 0;
	}

	void load(std::istream &/*in*/){
	    throw std::domain_error("not implemented");
	}

	uint64_t size() const {
		return size_;
	}

	ulint bit_size() const {return 0;}

private:

	vector<uint64_t> vec;
	uint64_t size_=0;


};

}

#endif /* INTERNAL_SPSI_CHECK_HPP_ */
