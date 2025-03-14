#ifndef LSH_H
#define LSH_H

#include <inttypes.h>
#include <vector>
#include <random>
#include <set>
#include <cmath>
#include <utility>
#include <iostream>
#include <chrono>

#include "common.h"

const unsigned int M_LSH = 4294967291;

using namespace std;

template <typename key_t>
struct hashtable {
	vector< vector< pair<key_t*, int> > > table;
	int capacity;  
	int size;

	// Default constructor
	hashtable() : size(0) {}

	// Constructor that initializes the hashtable with a given capacity
	hashtable(int _capacity) : capacity(_capacity), size(0) {
		table.resize(capacity+1);
	}

	// Insert (p, id) pair in hashtable
	void insert(key_t& p, int id) {
		int pos = id % capacity;
		table[pos].push_back(make_pair(&p, id));
		size++;
	}

	// Return the bucket associated with a given id
	vector< pair<key_t*, int> >& which_bucket(int id) {
		return table[id % capacity];
	}
};


class LSH {
	public:
	
		LSH(int _number_of_tables, int _number_of_h_funcs, vector< point_t >& _points, int _table_capacity, int _W_LSH = 10); 
		
		vector< pair< double, point_t* > > nearest_neighbor_search(int N, point_t& q_point, int dist_type);
		
		// Approximate range search
		vector< pair< double, point_t* > > approximate_range_search(double R, point_t& q_point, int dist_type);
		
		double last_func_time() {
			return timer;
		}
	
	private:
		int W_LSH;
		double timer;
		int number_of_tables;
		int number_of_h_funcs;
		int tables_capacity;

		vector< point_t >& points;

		vector< hashtable< point_t > > tables;

		vector< vector<hash_parameters> > hash_params;

		int h_func(point_t& p, hash_parameters& h_params);

		int id_func(int k, point_t& p, vector< hash_parameters >& g_params);

};


#endif