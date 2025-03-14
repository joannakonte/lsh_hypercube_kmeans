#include "../headers/dataset_input.h"

uint32_t convertBigEndianToNative(uint32_t value) {
	return	((value & 0xFF000000) >> 24)|
			((value & 0x00FF0000) >> 8) |
			((value & 0x0000FF00) << 8) |
			((value & 0x000000FF) << 24);
}

// This function is MNIST database specific!
void store_labels(vector< point_t >& points, string& file, input_file_info& file_info){
 	// Open the label file for reading
    ifstream inFile(file);
	if (!inFile.is_open()) {
        cerr << "Error opening label file!" << endl;
        return;
    }
	
	uint32_t temp;
	
	// Read 2 32numbers from the file 
    inFile.read(reinterpret_cast<char*>(&temp), sizeof(temp));
	inFile.read(reinterpret_cast<char*>(&temp), sizeof(temp));

	// Read labels for each point
	for (int i = 0; i < points.size(); i++) {
		inFile.read(reinterpret_cast<char*>(&points[i].label), sizeof(uint8_t));
	}

	inFile.close();
}

// this function is MNIST database specific!
void store_pointset(vector< point_t >& points, string& file, input_file_info& file_info) {
    // Open the dataset file for reading
    ifstream inFile(file);

    // Read the magic number from the file 
    inFile.read(reinterpret_cast<char*>(&file_info.magic), sizeof(file_info.magic));
    file_info.magic = convertBigEndianToNative(file_info.magic);
    printf("magic number: %d\n", file_info.magic);

    // Read the number of images from the file 
    inFile.read(reinterpret_cast<char*>(&file_info.num_images), sizeof(file_info.num_images));
    file_info.num_images = convertBigEndianToNative(file_info.num_images);
    printf("number of images: %d\n", file_info.num_images);

    // Read the number of rows for each image 
    inFile.read(reinterpret_cast<char*>(&file_info.num_rows), sizeof(file_info.num_rows));
    file_info.num_rows = convertBigEndianToNative(file_info.num_rows);
    printf("number of rows: %d\n", file_info.num_rows);

    // Read the number of columns for each image 
    inFile.read(reinterpret_cast<char*>(&file_info.num_columns), sizeof(file_info.num_columns));
    file_info.num_columns = convertBigEndianToNative(file_info.num_columns);
    printf("number of columns: %d\n", file_info.num_columns);

    // Resize the points vector to hold all images
    points.resize(file_info.num_images);

    // For each image in the dataset 
    for (int i=0; i<file_info.num_images; i++) {
        points[i].v.resize(file_info.num_rows * file_info.num_columns);
		points[i].array_pos = i;
        // For each pixel in the current image
        for (int j = 0; j < file_info.num_rows * file_info.num_columns; j++) {
            inFile.read(reinterpret_cast<char*>(&points[i].v[j]), sizeof(uint8_t));
        }
    }

	inFile.close();
}