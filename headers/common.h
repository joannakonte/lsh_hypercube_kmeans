#ifndef COMMON_H
#define COMMON_H

#include <inttypes.h>
#include <vector>
#include <random>
#include <cmath>
#include <utility>
#include <iostream>

using namespace std;

struct point_t {
	vector< uint8_t > v;
	uint8_t label;
	int array_pos;
	int cluster;
	double dist_from_cluster;

	bool operator<(const point_t& other) const {
		return label < other.label;
	}
};

struct hash_parameters {
	vector< double > v;		// Normal distribution vector
	double t;				// added noise
	int w;					// window size
	int r;					// random integer

	// Generate random parameters
	hash_parameters() {}
	hash_parameters(int d, int W) {
		v.resize(d);
		default_random_engine generator(random_device{}());
		normal_distribution<double> distribution(0.0, 1.0);
		
		for (int i = 0; i < d; i++) {
			v[i] = abs(distribution(generator));
		}
		
		w = W;
		r = generator();
   		t = uniform_real_distribution<double>(0, w)(generator);
	}
};

// Generate a random number (double)
double getRandomDoubleUniform(double min, double max);


// Distance function
double dist(point_t& x, point_t& y, int dist_type);

#endif