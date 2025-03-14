#include <iostream>
#include <string.h>

#include "../headers/lsh.h"
#include "../headers/common.h"
#include "../headers/brute_force.h"
#include "../headers/dataset_input.h"

#define Q 10

struct parameters {
	string input_file, labels_file, query_file, output_file, query_labels_file;
	int k;	// Νumber of h functions 
	int	L;	// Number of hashtables
	int	N;	// Looking for N Nearest Neighbours
	int	R;	// Range of search
	int labels_print;
	parameters() {
		input_file = ""; labels_file = ""; query_file = ""; 
		output_file = ""; query_labels_file = "";
		k = 4;		// Νumber of h functions 
		L = 5; 		// Number of hashtables
		N = 1; 		// Looking for N Nearest Neighbours
		R = 10000;	
		labels_print = 0;
	}

	void print() {
		cout << "Input file:" << input_file << endl; 
		cout << "Output file: " << output_file << endl;
		cout << query_file << endl;
		cout << k << endl;
		cout << L << endl;
		cout << N << endl; 
		cout << R << endl;
	}
} params;

void parse_flags(int argc, char**& argv) {

	if ( (argc-1)%2 != 0 ) {
		perror("Wrong flag parameters! Correct usage is:\n$./lsh [ -d <input file> ] [ -q <query file> ] [ -k <int> ] [ -L <int> ] [ -labels <labels_file> ] [-qlabels <query_labels_file>] [ -o <output file> ] [ -N <number of nearest> ] [ -R <radius> ] [ -labels_print <labels_print> ]\n");
		exit(-100);
	} 

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-d") == 0) {
			params.input_file = argv[++i];
		}
		else if (strcmp(argv[i], "-q") == 0) {
			params.query_file = argv[++i];
		}
		else if (strcmp(argv[i], "-k") == 0) {
			params.k = atoi(argv[++i]);
		}
		else if (strcmp(argv[i], "-L") == 0) {
			params.L = atoi(argv[++i]);
		}
		else if (strcmp(argv[i], "-labels") == 0) {
			params.labels_file = argv[++i];
		}
		else if (strcmp(argv[i], "-qlabels") == 0) {
			params.query_labels_file = argv[++i];
		}
		else if (strcmp(argv[i], "-o") == 0) {
			params.output_file = argv[++i];
		}
		else if (strcmp(argv[i], "-N") == 0) {
			params.N = atoi(argv[++i]);
		}
		else if (strcmp(argv[i], "-R") == 0) {
			params.R = atoi(argv[++i]);
		}
		else if (strcmp(argv[i], "-labels_print") == 0) {
			params.labels_print = atoi(argv[++i]);
		}
	}
}

int main(int argc, char* argv[]) {
	parse_flags(argc, argv);

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
	if (params.labels_print)
		store_labels(points, params.labels_file, file_info);

	ofstream outFile(params.output_file);
	
	LSH myLSH(params.L, params.k, points, points.size()/16);

	while (1) {
		if (params.query_file == "") {
			cout << "Insert query file: ";
			cin >> params.query_file;
		}

		// Check if the user wants to exit
		if (params.query_file == "exit") {
			break;
		}

		vector< point_t > query_points;

		store_pointset(query_points, params.query_file, file_info);
		if (params.labels_print)
			store_labels(query_points, params.query_labels_file, file_info);
	
		for (int i = 0; i < Q; i++) {
			double timer_brute_force; 
			vector< pair< double, point_t* > > approximate_NN = myLSH.nearest_neighbor_search(params.N, query_points[i], 2);
			vector< pair< double, point_t* > > range_NN = myLSH.approximate_range_search(params.R, query_points[i], 2);
			vector< pair< double, point_t* > > exact_NN = brute_force_NN(timer_brute_force, params.N, points, query_points[i], 2);
			
			outFile << "Query: " << query_points[i].array_pos << endl;
			if (params.labels_print) outFile << " Label: " << int(query_points[i].label) << endl;

			for (int j = 0; j < params.N; j++) {
				double dist_a = approximate_NN[j].first;
				uint8_t label_a = approximate_NN[j].second->label;
				double dist_e = exact_NN[j].first;
				uint8_t label_e = exact_NN[j].second->label;

				outFile << "Nearest neighbor-" << j+1 << ": "<< approximate_NN[j].second->array_pos << endl;
				outFile << "distanceLSH: " << dist_a  << endl;
				if (params.labels_print) outFile << "Label: " << int(label_a) << endl;
				outFile << "distanceTrue: " << dist_e << endl;
				if (params.labels_print) outFile << "Label: " << int(label_e) << endl;
			}
			outFile << "tLSH: " << myLSH.last_func_time() << endl;
			outFile << "tTrue: " << timer_brute_force << endl;
			outFile << "R-near neighbors:" << endl;
			for (int k = 0; k < range_NN.size(); k++) {
				outFile << range_NN[k].second->array_pos << endl;
				if (params.labels_print) outFile << "Label: " << int(range_NN[k].second->label) << endl;
			}
			outFile << endl;
		}

		params.query_file = "";
	}
	outFile.close();
}