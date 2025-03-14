#include "../headers/lsh.h"
	
LSH::LSH(int _number_of_tables, int _number_of_h_funcs, vector< point_t >& _points, int _table_capacity, int _W_LSH): 
	number_of_tables(_number_of_tables), number_of_h_funcs(_number_of_h_funcs), points(_points), tables_capacity(_table_capacity), W_LSH(_W_LSH) {

	tables.reserve(number_of_tables);;

	for (int i = 0; i < number_of_tables; i++) {
		tables.emplace_back(tables_capacity);
	}

	hash_params.resize(number_of_tables);
	for (int i = 0; i < number_of_tables; i++) {
		hash_params[i].resize(number_of_h_funcs);
		for (int j = 0; j < number_of_h_funcs; j++) {
			hash_params[i][j] = hash_parameters(points[0].v.size(), W_LSH);
		}
	}

	int number_of_points = points.size();
	for (int i = 0; i < number_of_points; i++) {
		for (int j = 0; j < number_of_tables; j++) {
			int id = id_func(number_of_h_funcs, points[i], hash_params[j]);
			tables[j].insert(points[i], id);
		}
	}
}


vector< pair< double, point_t* > > LSH::nearest_neighbor_search(int N, point_t& q_point, int dist_type) {
	// Record the start time
	auto start = chrono::high_resolution_clock::now();
	
	set< pair< double, point_t* > > best;
	int candidates_checked = 0;
	int L = number_of_tables;
	int k = number_of_h_funcs;

	for (int i = 0; i < L; i++){
		// Find the bucket id for the query point q_point using the id_func
		int id_q = id_func(k, q_point, hash_params[i]);
		
		vector< pair< point_t*, int > > bucket = tables[i].which_bucket(id_q);

		for (int j = 0; j < (int)bucket.size(); j++) {
			candidates_checked++;
			point_t* candidate = bucket[j].first;
			int candidate_id = bucket[j].second;



			if ((int)best.size() < N) {
				best.insert(make_pair(dist(q_point, (*candidate), dist_type), candidate));
				continue;
			}
			
			double candidate_dist = dist(q_point, (*candidate), dist_type);

			double cur_largest_dist =  (best.rbegin()->first); // Get the largest distance from the set
			if (candidate_dist < cur_largest_dist) {
				best.erase(prev(best.end()));			       // Remove last element from set
				best.insert(make_pair(candidate_dist, candidate));
			}

			if (candidates_checked > 10*L) break;
		}
		if (candidates_checked > 10*L) break;
	}
	
	// Record the end time
	auto end = chrono::high_resolution_clock::now();
	
	chrono::duration<double> duration = chrono::duration_cast<chrono::duration<double>>(end - start);
	timer = duration.count();
	
	// Convert set to vector
	vector< pair< double, point_t* > > best_vector(best.begin(), best.end());
	return best_vector;
}

vector< pair< double, point_t* > > LSH::approximate_range_search(double R, point_t& q_point, int dist_type) {
	vector< pair<double, point_t*> > points_in_range;
	int candidates_checked = 0;

	int L = number_of_tables;
	int k = number_of_h_funcs;

	for (int i = 0; i < L; i++){
		// Find the bucket id for the query point q_point using the id_func
		int id_q = id_func(k, q_point, hash_params[i]);
		
		vector< pair<point_t*, int> > bucket = tables[i].which_bucket(id_q);

		// For each point in the same bucket as q_point
		for (int j = 0; j < (int)bucket.size(); j++){
			candidates_checked++;
			point_t* candidate = bucket[j].first;
			int candidate_id = bucket[j].second;

			// Calculate the distance between the point and the q_point
			double distance_from_q = dist(q_point, (*candidate), dist_type);
			if (distance_from_q < R){
				points_in_range.push_back(make_pair(distance_from_q, candidate));
			}

			if (candidates_checked > 20*L) break;
		}
		if (candidates_checked > 20*L) break;
	}

	return points_in_range;
}


int LSH::h_func(point_t& p, hash_parameters& h_params) {
	double h = 0.0; 
	for (int i = 0; i < p.v.size(); i++) {
		h += p.v[i] * h_params.v[i];
	}
	h += h_params.t;
	h /= h_params.w;
	return static_cast<int>(floor(h));
}

int LSH::id_func(int k, point_t& p, vector<hash_parameters>& g_params) {
	uint8_t sum = 0;

	//	using modulo properties
	//	(a + b) mod n = [(a mod n) + (b mod n)] mod n.
	//	ab mod n = [(a mod n)(b mod n)] mod n.

	for (int i = 0; i < k; i++) {
		int r = g_params[i].r;
		int M = M_LSH;
		int term = ((r % M) * (h_func(p, g_params[i]) % M)) % M;
		sum = (sum + term) % M;
	} 

	return sum;
}