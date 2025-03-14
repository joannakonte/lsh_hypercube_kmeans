#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <string.h>
#include <iomanip>

#include "../headers/common.h"
#include "../headers/cluster.h"
#include "../headers/dataset_input.h"

using namespace std;

struct parameters {
	string input_file, labels_file, conf_file, output_file;
	int clusters;       // Number of clusters
	int L;      		// Number of vector hash tables
	int	k_LSH;	        // Number of vector hash functions
	int	M;			    // M of Hypercube
	int k_hypercube;    // Number of hypercube dimensions
	int	probes;	        // Probes of hypercube
	bool complete;
	string method;

	parameters() {
		input_file = ""; labels_file = ""; conf_file = ""; 
		output_file = "";
		// Default values
		L = 3;      		// Number of vector hash tables
		k_LSH = 4;	        // Number of vector hash functions
		M = 10;			    // M of Hypercube
		k_hypercube = 3;    // Number of hypercube dimensions
		probes = 2;	        // Probes of hypercube
		complete = false;
	}

	void print() {
		cout << "Input file:" << input_file << endl; 
		cout << "Output file: " << output_file << endl; 
		cout << "Clusters: " << clusters << endl;
		cout << "L: " << L << endl;
		cout << "k_lsh: " << k_LSH << endl;
		cout << "M: " << M << endl;
		cout << "k_hypercube: " << k_hypercube << endl; 
		cout << "probes: " << probes << endl;
		cout << "Complete: " << complete << endl;
		cout << "Method: " << method << endl;
	}
} params;

void loadParametersFromFile(const string &filename) {
	ifstream file(filename);
	string line, key;
	int value;

	if (file.is_open()) {
		while (getline(file, line)) {
			stringstream ss(line);
			getline(ss, key, ':');
			
			if (key == "number_of_clusters") {
				ss >> params.clusters;
			} else if (key == "number_of_vector_hash_tables") {
				ss >> params.L;
			} else if (key == "number_of_vector_hash_functions") {
				ss >> params.k_LSH;
			} else if (key == "max_number_M_hypercube") {
				ss >> params.M;
			} else if (key == "number_of_hypercube_dimensions") {
				ss >> params.k_hypercube;
			} else if (key == "number_of_probes") {
				ss >> params.probes;
			}
		}
		file.close();
	} else {
		cout << "Unable to open file";
	}
}

void parse_flags(int argc, char**& argv) {
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-complete") == 0) 
			params.complete = true;
	}

	if ( (argc-1-params.complete)%2 != 0 ) {
		perror("Wrong flag parameters! Correct usage is:\n$./lsh [ -i <input file> ] [ -c <configuration file> ] [ -o <output file> ] [ -complete <optional> ] [ -m <method> ]\n");
		exit(-100);
	}

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-i") == 0) {
			params.input_file = argv[++i];
		}
		if (strcmp(argv[i], "-labels") == 0) {
			params.labels_file = argv[++i];
		}
		else if (strcmp(argv[i], "-c") == 0) {
			params.conf_file = argv[++i];
		}
		else if (strcmp(argv[i], "-o") == 0) {
			params.output_file = argv[++i];
		}
		else if (strcmp(argv[i], "-complete") == 0) {
			continue;
		}
		else if (strcmp(argv[i], "-m") == 0) {
			params.method = argv[++i];
		}
	}
}

bool compare_by_cluster(point_t& a, point_t& b) {
	return a.cluster < b.cluster;
}

int main(int argc, char* argv[]) {
	parse_flags(argc, argv);
	loadParametersFromFile(params.conf_file);
	params.print();

	vector< point_t > points;
	
	if (params.input_file == "") {
		cout << "Insert input file: ";
		cin >> params.input_file;
	}
	if (params.output_file == "") {
		cout << "Insert output file: ";
		cin >> params.output_file;
	}

	input_file_info file_info;
	store_pointset(points, params.input_file, file_info);

	ofstream outFile(params.output_file);

	clusterer myclusterer(params.clusters, points);

	myclusterer.kmeanspp_init();
	if (params.method == "Classic") {
		myclusterer.lloyd();
		outFile << "Algorithm: Lloyds" << endl; 
	}
	else if (params.method == "LSH") {
		myclusterer.lsh_cluster(params.L, params.k_LSH, points.size()/32);
		outFile << "Algorithm: Range Search LSH" << endl;
	}

	else if (params.method == "Hypercube") {
		myclusterer.cube_cluster(params.k_hypercube, params.probes);
		outFile << "Algorithm: Range Search Hypercube" << endl;
	}

	auto centroids = myclusterer.get_centroids();
	for (int i = 0; i < params.clusters; i++) {
		outFile << "CLUSTER-" << i;
		outFile << " { size: " << myclusterer.cluster_size[i] << " centroid: " << endl;
		// for (int j = 0; j < centroids[i].v.size(); j++) {
		//     outFile << int(centroids[i].v[j]) << " " ;
		// }
		int j = 0;
		for (int r = 0; r < file_info.num_rows; r++) {
			for (int c = 0; c < file_info.num_columns; c++) {
				outFile << setw(3) << int(centroids[i].v[j]) << " ";
				j++;
			}
			outFile << endl;
		}
		outFile << "} " << endl;
	}

	outFile << "clustering_time: " << myclusterer.last_func_time() << endl;
	
	
	if (params.complete) {
		sort(points.begin(), points.end(), compare_by_cluster);
		int idx = 0;
		for (int i = 0; i < params.clusters; i++) {
			outFile << "CLUSTER-" << i << "{";
			while (points[idx].cluster == i) {
				idx++;
				outFile << points[idx].array_pos << ", ";
			}
			outFile << "}";
			break;
		}
	}

	vector<double> silhoutte_scores = myclusterer.silhoutte(2);
	outFile << "Silhouette: [";
	for (int i = 0; i < silhoutte_scores.size()-1; i++)
		outFile << "s" << i <<  " = " << silhoutte_scores[i] << ", ";
	outFile << "stotal = " << silhoutte_scores[silhoutte_scores.size()-1] << "]" << endl;

	outFile << "silhouette time = " << myclusterer.last_func_time() << endl;

	outFile.close();
	return 0;
}