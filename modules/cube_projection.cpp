#include "../headers/cube_projection.h"
	
hypercube_projection::hypercube_projection(int _d, vector< point_t >& _points, int _W_cube): 
d(_d), points(_points), tables(1<<d), f_func(d), W_cube(_W_cube) {
	tables_capacity = 1<<d;
	
	hash_params.resize(d);
	
	for (int i = 0; i < d; i++) {
		hash_params[i] = hash_parameters(points[0].v.size(), W_cube);
	}
	
	tables.resize(1<<d);

	for (int i = 0; i < points.size(); i++) {
		tables[id_func(points[i], d, hash_params)].push_back(&(points[i]));
	}

}


vector< pair< double, point_t* > > hypercube_projection::nearest_neighbor_search(int N, int probes, int max_candidates, point_t& q_point, int dist_type=2) {
	auto start = chrono::high_resolution_clock::now();
	
	if (max_candidates == -1) max_candidates = points.size()+1;
	
	// Get hash value of q_point
	ll hash_q = id_func(q_point, d, hash_params);
	
	set< pair< double, point_t* > > best;
	int candidates_checked = 0;
	int corners_checked = 0;
	
	for (int corner_dist = 0; corner_dist <= d; corner_dist++) {
		if (candidates_checked >= max_candidates || corners_checked >= probes) break;
		
		for (int i = 0; i < tables.size(); i++) {
			if (candidates_checked >= max_candidates || corners_checked >= probes) break;

			if (hamming_dist(ll(i), hash_q) != corner_dist) continue; 
			
			corners_checked++;
			
			for (int j = 0; j < tables[i].size(); j++) {
				if (candidates_checked >= max_candidates) break;

				candidates_checked++;
			
				point_t* candidate = tables[i][j];
				double candidate_dist = dist(q_point, (*candidate), dist_type);

				if (best.size() < N) {
					best.insert(make_pair(candidate_dist, candidate));
					continue;
				}

				double cur_largest_dist =  (best.rbegin()->first); // Get the largest distance from the set
				if (candidate_dist < cur_largest_dist) {
					best.erase(prev(best.end()));			       // Remove last element from set
					best.insert(make_pair(candidate_dist, candidate));
				}
			}
		}
	}

	auto end = chrono::high_resolution_clock::now();
	
	chrono::duration<double> duration = chrono::duration_cast<chrono::duration<double>>(end - start);
	timer = duration.count();


	// Convert set to vector
	vector< pair< double, point_t* > > best_vector(best.begin(), best.end());
	return best_vector;
}

vector< pair< double, point_t* > > hypercube_projection::approximate_range_search(int R, int probes, int max_candidates, point_t& q_point, int dist_type) {
	if (max_candidates == -1) max_candidates = points.size()+1;
	
	auto start = chrono::high_resolution_clock::now();
	
	vector< pair<double, point_t*> > points_in_range;

	// Get hash value of q_point
	ll hash_q = id_func(q_point, d, hash_params);
	
	int candidates_checked = 0;
	int corners_checked = 0;

	for (int corner_dist = 0; corner_dist <= d; corner_dist++) {
		for (int i = 0; i < tables.size(); i++) {
			
			if (hamming_dist(ll(i), hash_q) != corner_dist) continue; 
			
			corners_checked++;
			
			for (int j = 0; j < tables[i].size(); j++) {
			
				candidates_checked++;
			
				point_t* candidate = tables[i][j];
				double candidate_dist = dist(q_point, (*candidate), dist_type);

				if (candidate_dist < R) {
					points_in_range.push_back(make_pair(candidate_dist, candidate));
				}
				if (candidates_checked > max_candidates) break;
			}
		
			if (corners_checked >= probes || candidates_checked > max_candidates) break;
		}
		if (corners_checked >= probes || candidates_checked > max_candidates) break;
	}

	auto end = chrono::high_resolution_clock::now();
	
	chrono::duration<double> duration = chrono::duration_cast<chrono::duration<double>>(end - start);
	timer = duration.count();

	return points_in_range;
}


// Hamming distance function
int hypercube_projection::hamming_dist(ll x, ll y) {
	long long xorResult = x ^ y; 
    int distance = __builtin_popcountll(xorResult); // Count the number of bits with value 1
	return distance;
}

int hypercube_projection::h_func(point_t& p, hash_parameters& h_params) {
	double h = 0.0; 
	for (int i = 0; i < p.v.size(); i++) {
		h += p.v[i] * h_params.v[i];
	}
	h += h_params.t;
	h /= h_params.w;
	return static_cast<int>(floor(h));
}

ll hypercube_projection::id_func(point_t& p, int d, vector< hash_parameters >& h_params) {
	ll hash_val = 0;
	
	for (int i = 0; i < d; i++) {
		ll h_value = h_func(p, h_params[i]);

		bool bit = 0;

		auto it = f_func[i].find(h_value);
		
		if (it == f_func[i].end()) {
			bit = rand()%2;
			f_func[i].insert(make_pair(h_value, bit));
		}
		else {
			bit = it->second;
		}

		if (bit) hash_val ^= 1<<i;
	}

	return hash_val;
}