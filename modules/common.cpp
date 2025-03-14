#include "../headers/common.h"

// Distance function
double dist(point_t& x, point_t& y, int dist_type) {
	double sum = 0.0;
	for (int i = 0; i < x.v.size(); i++) {
		double diff = x.v[i] - y.v[i];
		sum += pow(fabs(diff), dist_type);
	}
	sum = pow(sum, (1/(double)dist_type));
	return sum;
}

// Generate a random number (double)
double getRandomDoubleUniform(double min, double max) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(min, max);

    return dis(gen);
}
