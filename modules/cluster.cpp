#include "../headers/cluster.h"

clusterer::clusterer(int _clusters, vector<point_t>& _points):
			number_of_clusters(_clusters), points(_points), 
			centroids(_clusters), cluster_size(_clusters) {
}

bool compare(const pair<double, point_t*>& a, const pair<double, point_t*>& b) {
    return a.first < b.first;
}

void clusterer::kmeanspp_init(int dist_type){
	vector< point_t* > c;
	int first_centroid = rand() % points.size();    		// Choose the first centroid randomly
// Initialization++
	c.push_back(&(points[first_centroid]));
	
	int K = number_of_clusters;

	// Initialize distance from cluster
	for (int i = 0; i < points.size(); i++) {
		points[i].cluster = -1;
		points[i].dist_from_cluster = DBL_MAX;
	}
	
	for (int k = 1; k < K; k++) {
		vector< pair<double, point_t*> > distances;			// Vector with all of the distances of the points from the cetnroid
		double max_dist = 0;

		for (int i = 0; i < points.size(); i++) {
			// If point is a centroid continue
			if (find(c.begin(), c.end(), &(points[i])) != c.end())
				continue;

			double min_dist = DBL_MAX;
			// Find the min distance from centroid
			for (int j = 0; j < c.size(); j++) {
				min_dist = min(dist(*(c[j]), points[i], dist_type), min_dist);
			}
			max_dist = max(max_dist, min_dist);
			distances.push_back(make_pair(min_dist, &points[i]));
		}
		
		// Normalize distances
		for (int i = 0; i < distances.size(); i++) {
			distances[i].first /= max_dist;
			distances[i].first *= distances[i].first;
			if (i != 0) 
				// Compute the partial sums
				distances[i].first += distances[i-1].first;
		}
		// Get a random variable from the distances array
		double x = getRandomDoubleUniform(0, distances[distances.size()-1].first);

		pair<double, point_t*> dummy = make_pair(x, &points[0]);
		// Find the next centroid
		int ind = upper_bound(distances.begin(), distances.end(), dummy) - distances.begin();
		c.push_back(distances[ind].second);
	}

	for (int i = 0; i < centroids.size(); i++) {
		centroids[i] = *(c[i]);
	}
}

void clusterer::calc_new_average(point_t& point, int new_cluster, int old_cluster) {
	// fix new cluster
	if (new_cluster == old_cluster) return ;
	int cur_size = cluster_size[new_cluster] + 1;
	for (int j = 0; j < centroids[new_cluster].v.size(); j++) {
		long long dim = (long long)centroids[new_cluster].v[j];
		dim *= cur_size;
		dim += point.v[j];
		dim = ceil((double)dim/(cur_size + 1));
		if (centroids[new_cluster].v[j] != 0 && dim == 0) printf("will change from %d to %lld\n", centroids[new_cluster].v[j], dim);
		centroids[new_cluster].v[j] = (uint8_t)dim;
	}
	cluster_size[new_cluster]++;
	
	// fix previous cluster
	
	if (old_cluster == -1) return; 
	
	cur_size = cluster_size[old_cluster] + 1;
	for (int j = 0; j < centroids[old_cluster].v.size(); j++) {
		long long dim = (long long)centroids[old_cluster].v[j];
		dim *= cur_size;
		dim -= point.v[j];
		dim /= (cur_size - 1);
		centroids[old_cluster].v[j] = dim;
	}
	cluster_size[old_cluster]--;
}

void clusterer::lloyd(int dist_type) {
	// Record the start time
	auto start = chrono::high_resolution_clock::now();

	int how_many_changed = 0;
	do {
		how_many_changed = 0;
		// For every point
		for (int i = 0; i < points.size(); i++) {
			double& min_dist = points[i].dist_from_cluster;
			int prev_cluster = points[i].cluster;
			// Compute distance to every centroid
			for (int j = 0; j < centroids.size(); j++) {
				if (j == prev_cluster) continue;

				double cur_dist = dist(centroids[j], points[i], dist_type);
				
				if (cur_dist < min_dist) {
					points[i].cluster = j;
					// Nearest centroid
					min_dist = cur_dist;
				}

			}

			int current_cluster = points[i].cluster;
			if (current_cluster != prev_cluster) {
				how_many_changed++;
				calc_new_average(points[i], current_cluster, prev_cluster);
			}
		}
	
	} while (how_many_changed > 0);

	// Record the end time
	auto end = chrono::high_resolution_clock::now();
	
	chrono::duration<double> duration = chrono::duration_cast<chrono::duration<double>>(end - start);
	timer = duration.count();
}

void clusterer::lsh_cluster(int number_of_tables, int number_of_h_funcs, int table_capacity, int W_LSH, int dist_type) {
	// Record the start time
	auto start = chrono::high_resolution_clock::now();

	LSH mylsh(number_of_tables, number_of_h_funcs, points, table_capacity, W_LSH);

	int how_many_changed = 0;	
	vector< int > previous_cluster(points.size(), -1);
	int macqueen_loops_cnt = 0;
	do {	// macqueen loop, range starts from the beginning

		// Deassign the points
		for (int i = 0; i < points.size(); i++) {
			previous_cluster[i] = points[i].cluster;
			points[i].cluster = -1;
			points[i].dist_from_cluster = DBL_MAX;
		}


		// Empty the clusters
		for (int i = 0; i < cluster_size.size(); i++) {
			cluster_size[i] = 0;
		}
		
		double R_start = DBL_MAX;	// Range
		// Find distances between centroids
		for (int i = 0; i < centroids.size(); i++) {
			for (int j = i+1; j < centroids.size(); j++) {
				R_start = min(R_start, dist(centroids[i], centroids[j], dist_type));
			}
		}

		how_many_changed = 0;

		double R = R_start/2;

		int loops_cnt = 0;


		do {	// LSH loop, range doubles each time
			// Keep the recently assigned points
			// so that we update their center
			vector< point_t* > who_changed;
			
			// For a specific range, this will be used to resolve conflicts
			// pair<int, double> = (best_cluster_id, dist to that cluster)
			vector< pair<int, double> > unasigned_best(points.size(), pair<int, double>(-1, DBL_MAX)); 

			// Range search on every centroid with the same range
			for (int i = 0; i < centroids.size(); i++) {
				vector< pair< double, point_t* > > neighbours = mylsh.approximate_range_search(R, centroids[i], dist_type);
				for (int j = 0; j < neighbours.size(); j++) {
					point_t* point = neighbours[j].second;
					double new_dist = neighbours[j].first;
					
					if (point->cluster != -1) // Already assigned 
						continue;

					if (unasigned_best[point->array_pos].first == -1) // Cluster will change, so insert in the vector
						who_changed.push_back(point);

					if (new_dist < unasigned_best[point->array_pos].second) {
						unasigned_best[point->array_pos].second = new_dist;
						unasigned_best[point->array_pos].first = i;
					}
				}
			}
			
			for (int i = 0; i < who_changed.size(); i++) {
				point_t* point = who_changed[i];
				point->cluster = unasigned_best[point->array_pos].first;
				point->dist_from_cluster = unasigned_best[point->dist_from_cluster].second;
				calc_new_average(*point, point->cluster, -1);
			}

			R*=2;
			loops_cnt++;
		} while (loops_cnt<5);

		
		// Assign clusters to those with no clusters found in the range search
		for (int i = 0; i < points.size(); i++) {
			if (points[i].cluster != -1) // Already assigned a cluster
				continue;

			// Points that have not been assigned a cluster
			double& min_dist = points[i].dist_from_cluster;

			// Compute distance to every centroid
			for (int j = 0; j < centroids.size(); j++) {

				double cur_dist = dist(centroids[j], points[i], dist_type);
				
				if (cur_dist < min_dist) {
					points[i].cluster = j;
					// Nearest centroid
					min_dist = cur_dist;
				}
			}

			calc_new_average(points[i], points[i].cluster, -1);
		}
		
		
		for (int i = 0; i < points.size(); i++) {
			if (previous_cluster[i] != points[i].cluster) how_many_changed++;
		}
		

		macqueen_loops_cnt++;
	} while (how_many_changed > 0 && macqueen_loops_cnt < 5);

	// Record the end time
	auto end = chrono::high_resolution_clock::now();
	
	chrono::duration<double> duration = chrono::duration_cast<chrono::duration<double>>(end - start);
	timer = duration.count();
}

void clusterer::cube_cluster(int d, int probes, int max_candidates, int W_cube, int dist_type) {
	// Record the start time
	auto start = chrono::high_resolution_clock::now();

	hypercube_projection hypercube(d, points, W_cube);

	int how_many_changed = 0;	
	vector< int > previous_cluster(points.size(), -1);
	int macqueen_loops_cnt = 0;
	do {	// macqueen loop, range starts from the beginning

		// Deassign the points
		for (int i = 0; i < points.size(); i++) {
			previous_cluster[i] = points[i].cluster;
			points[i].cluster = -1;
			points[i].dist_from_cluster = DBL_MAX;
		}


		// Empty the clusters
		for (int i = 0; i < cluster_size.size(); i++) {
			cluster_size[i] = 0;
		}
		
		double R_start = DBL_MAX;	// Range
		// Find distances between centroids
		for (int i = 0; i < centroids.size(); i++) {
			for (int j = i+1; j < centroids.size(); j++) {
				R_start = min(R_start, dist(centroids[i], centroids[j], dist_type));
			}
		}

		how_many_changed = 0;

		double R = R_start/2;

		int loops_cnt = 0;


		do {	// LSH loop, range doubles each time

			// Keep the recently assigned points
			// so that we update their center
			vector< point_t* > who_changed;
			
			// For a specific range, this will be used to resolve conflicts
			// pair<int, double> = (best_cluster_id, dist to that cluster)
			vector< pair<int, double> > unasigned_best(points.size(), pair<int, double>(-1, DBL_MAX)); 

			// Range search on every centroid with the same range
			for (int i = 0; i < centroids.size(); i++) {
				vector< pair< double, point_t* > > neighbours = hypercube.approximate_range_search(R, probes, max_candidates, centroids[i], dist_type);
				for (int j = 0; j < neighbours.size(); j++) {
					point_t* point = neighbours[j].second;
					double new_dist = neighbours[j].first;
					
					if (point->cluster != -1) // Already assigned 
						continue;

					if (unasigned_best[point->array_pos].first == -1) // Cluster will change, so insert in the vector
						who_changed.push_back(point);

					if (new_dist < unasigned_best[point->array_pos].second) {
						unasigned_best[point->array_pos].second = new_dist;
						unasigned_best[point->array_pos].first = i;
					}
				}
			}
			
			for (int i = 0; i < who_changed.size(); i++) {
				point_t* point = who_changed[i];
				point->cluster = unasigned_best[point->array_pos].first;
				point->dist_from_cluster = unasigned_best[point->dist_from_cluster].second;
				calc_new_average(*point, point->cluster, -1);
			}

			R*=2;
			loops_cnt++;
		} while (loops_cnt<5);

		
		// Assign clusters to those with no clusters found in the range search
		for (int i = 0; i < points.size(); i++) {
			if (points[i].cluster != -1) // Already assigned a cluster
				continue;

			// Points that have not been assigned a cluster
			double& min_dist = points[i].dist_from_cluster;

			// Compute distance to every centroid
			for (int j = 0; j < centroids.size(); j++) {

				double cur_dist = dist(centroids[j], points[i], dist_type);
				
				if (cur_dist < min_dist) {
					points[i].cluster = j;
					// Nearest centroid
					min_dist = cur_dist;
				}
			}

			calc_new_average(points[i], points[i].cluster, -1);
		}
		
		
		for (int i = 0; i < points.size(); i++) {
			if (previous_cluster[i] != points[i].cluster) how_many_changed++;
		}
		
		macqueen_loops_cnt++;
	} while (how_many_changed > 0 && macqueen_loops_cnt < 5);

	// Record the end time
	auto end = chrono::high_resolution_clock::now();
	
	chrono::duration<double> duration = chrono::duration_cast<chrono::duration<double>>(end - start);
	timer = duration.count();
}

/* =================== | Silhoutte Score | =================== */

// Put points in clusters 
void clusterer::split_clusters(vector< vector< silhoutte_t > >& cluster_elements) {
	cluster_elements.resize(number_of_clusters);

	for (int i = 0; i < points.size(); i++) {
		cluster_elements[points[i].cluster].push_back((silhoutte_t){&points[i], 0, 0});
	}

}

double clusterer::find_dist(point_t &a, point_t &b, int dist_type) {
	// (x y) -> dist (x y)
	auto search_pair = make_pair(a.array_pos, b.array_pos);
	
	// Make pair with the smallest point as the first element
	if (search_pair.first > search_pair.second) 
		swap(search_pair.first, search_pair.second);
	
	// Get iterator
    auto it = calculated_dists.find(search_pair);	

	// If distance is found, return
    if (it != calculated_dists.end()) {
		return it->second;
	}
	// If distance is not found, calculate it
	else {
		double dist_a_b = dist(a, b, dist_type);
		calculated_dists.insert(make_pair(search_pair, dist_a_b));
		return dist_a_b;
	}

}

vector<double> clusterer::silhoutte(int dist_type) {
	auto start = chrono::high_resolution_clock::now();
	
	vector< vector< silhoutte_t > > cluster_elements;

	split_clusters(cluster_elements);
	// Find the second best cluster for each point
	for (int i = 0; i < cluster_elements.size(); i++) {
		for (auto& element: cluster_elements[i]) {
			
			double min_dist = DBL_MAX;
			double min_pos;
			for (int j = 0; j < centroids.size(); j++) {		// Find the distance to each centroid
				if (j == element.point->cluster) continue;
				
				double cur_dist = dist(*element.point, centroids[j], dist_type);
				if (cur_dist < min_dist) {
					min_dist = cur_dist;
					min_pos = j;
				}
			}
			element.second_best = min_pos;
		}
	}
	vector<double> cluster_silhouette(number_of_clusters, 0);		// Vector with Silhoutte Score for each cluster
	double S_total = 0;

	// Calculate a(i) and b(i)
	for (auto& cluster: cluster_elements) {
		
		for (int i = 0; i < cluster.size(); i++) {
			point_t* point_i = cluster[i].point;
			
			// Calculate distances of points in the same cluster
			for (int j = i + 1; j < cluster.size(); j++) {
				point_t* point_j = cluster[j].point;
				double dist_i_j = dist(*point_i, *point_j, dist_type);
				//dist_i_j /= cluster.size();
				cluster[i].a += dist_i_j;
				cluster[j].a += dist_i_j;
			}
			
			// Calculate distances of points in the second best cluster
			auto& cluster_b = cluster_elements[cluster[i].second_best];			// Get next best cluster for current point i 
			for (int j = 0; j < cluster_b.size(); j++) {						// For all points of second best cluster
				point_t* point_j = cluster_b[j].point;
				double dist_i_j = find_dist(*point_i, *point_j, dist_type);			// Calculate distance of point from all other points in the second best cluster
				cluster[i].b += dist_i_j;
			}

			// Get the average of both
			cluster[i].a /= (cluster.size()-1); // a(i)
			cluster[i].b /= cluster_b.size();	// b(i)
			// Calculate s(i)
			double point_silhouette = (cluster[i].b - cluster[i].a) / max(cluster[i].b, cluster[i].a);
			cluster_silhouette[point_i->cluster] += point_silhouette;
			S_total += point_silhouette;
		}
	}

	// Overall Silhoutte coefficient (averate s(i)/number_of_points)
	for (int i = 0; i < number_of_clusters; i++) {
		cluster_silhouette[i] /= cluster_size[i];
	}

	S_total /= points.size();

	cluster_silhouette.push_back(S_total);

	auto end = chrono::high_resolution_clock::now();
	
	chrono::duration<double> duration = chrono::duration_cast<chrono::duration<double>>(end - start);
	timer = duration.count();

	return cluster_silhouette;
}
