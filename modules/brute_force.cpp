#include "../headers/brute_force.h"



vector< pair< double, point_t* > > brute_force_NN(double& timer_brute_force, int N, vector< point_t >& points, point_t& q_point, int dist_type) {
	// Record the start time
	auto start = chrono::high_resolution_clock::now();

	set< pair<double, point_t*> > best;
	for (int i=0; i<points.size(); i++) {
		
		best.insert(make_pair(dist(points[i], q_point, 2), &(points[i])));
		
		if (best.size() > N) 
			best.erase(prev(best.end()));	
	}

	vector< pair< double, point_t* > > best_vector(best.begin(), best.end());

	// Record the end time
	auto end = chrono::high_resolution_clock::now();

	chrono::duration<double> duration = chrono::duration_cast<chrono::duration<double>>(end - start);
	timer_brute_force = duration.count();

    return best_vector;
}