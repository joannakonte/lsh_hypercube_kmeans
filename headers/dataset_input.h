#ifndef DATASET_INPUT
#define DATASET_INPUT

#include <fstream>

#include "common.h"

using namespace std;

struct input_file_info {
	uint32_t magic, num_images, num_rows, num_columns;
};

void store_labels(vector< point_t >& points, string& file, input_file_info& file_info);

void store_pointset(vector< point_t >& points, string& file, input_file_info& file_info);


#endif