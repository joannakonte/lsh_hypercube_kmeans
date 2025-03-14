#ifndef BRUTE_FORCE_H
#define BRUTE_FORCE_H

#include <inttypes.h>
#include <vector>
#include <random>
#include <set>
#include <cmath>
#include <utility>
#include <iostream>
#include <chrono>

#include "common.h"

using namespace std;

vector< pair< double, point_t* > > brute_force_NN(double& timer_brute_force, int N, vector< point_t >& points, point_t& q_point, int dist_type);

#endif