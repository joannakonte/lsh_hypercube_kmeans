#ifndef CUBE_PROJECTION_H
#define CUBE_PROJECTION_H

#include <inttypes.h>
#include <vector>
#include <random>
#include <set>
#include <cmath>
#include <utility>
#include <iostream>
#include <map>
#include <chrono>

#include "common.h"

using namespace std;

const uint M_cubes = 4294967291;

typedef long long ll;

class hypercube_projection {
	public:
	
		hypercube_projection(int d, vector< point_t >& _points, int _W_cube = 10);
		
		// if max_candidates == -1, it is ignored
		vector< pair< double, point_t* > > nearest_neighbor_search(int N, int probes, int max_candidates, point_t& q_point, int dist_type);
		
		// if max_candidates == -1, it is ignored
		vector< pair< double, point_t* > > approximate_range_search(int R, int probes, int max_candidates, point_t& q_point, int dist_type);
		
		double last_func_time() {
			return timer;
		}

	private:
		int W_cube;
		double timer;
		int d;
		int tables_capacity;

		vector< point_t >& points;

		vector< vector< point_t* > > tables; 		 

		vector< hash_parameters > hash_params;

		vector < map< ll, bool > > f_func;
		
		// Distance function
		int hamming_dist(ll x, ll y);

		int h_func(point_t& p, hash_parameters& h_params);
		ll id_func(point_t& p, int d, vector< hash_parameters >& h_params);

};


#endif