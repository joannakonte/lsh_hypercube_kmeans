#ifndef CLUSTER_H
#define CLUSTER_H

#include <algorithm>
#include <float.h>
#include <chrono>
#include <map>

#include "common.h"
#include "../headers/lsh.h"
#include "../headers/cube_projection.h"

using namespace std;

class clusterer {
	public: 

		clusterer(int _clusters, vector<point_t>& _points);
		
		void kmeanspp_init(int dist_type = 2);
		
		// To be called AFTER kmeans_init
		const vector< point_t >& get_centroids() {
			return centroids;
		}

		// ============= | Clustering Algorithms | ============= //

		// To be called AFTER kmeans_init
		void lloyd(int dist_type = 2);
		
		// To be called AFTER kmeans_init
		void lsh_cluster(int number_of_tables, int number_of_h_funcs, int table_capacity, int W_LSH = 10, int dist_type = 2);
		
		// To be called AFTER kmeans_init
		// max_candidates doesnt really make sense, so preferrably use -1 to ignore it
		void cube_cluster(int d, int probes, int max_candidates = -1, int W_cube=10, int dist_type=2);	

		// ===================================================== //

		// To be called AFTER some clustering algorithm
		
		vector<double> silhoutte(int dist_type);
		
		vector< int > cluster_size;		// Number of points in each cluster

		double last_func_time() {
			return timer;
		}


		
	private:
		double timer;

		// ================== | Silhouette | ================== //
		map< pair< int, int >, double > calculated_dists;

		struct silhoutte_t {
			point_t* point;
			double a;
			double b;
			int second_best;	// Second best cluster
		};

		void split_clusters(vector< vector< silhoutte_t > >& cluster_elements);
		
		// This function is an optimisation to the classic algorithm,
		// it keeps track of previously calculated distances between points
		// and returns it in case we need it again
		double find_dist(point_t &a, point_t &b, int dist_type);
		
		// ==================================================== //

		vector< point_t > centroids;

		vector< point_t >& points;
		
		void calc_new_average(point_t& point, int new_cluster, int old_cluster);

		int number_of_clusters;
		
};

#endif