// NeuralNetwork.hpp 
#include <iostream> 
#include <vector> 
#include "ColVector.hpp"

using namespace std; 

// use typedefs for future ease for changing data types like : float to double 

ColVector::ColVector(std::vector<float> vector) 
{
	this->vector = matrix(vector, vector.size());
};

ColVector::ColVector(uint size) 
{
	std::vector<float> vector (size, 0.0f);
	this->vector = matrix(vector, vector.size());
};

float ColVector::coeffRef(uint pos) 
{
	return(vector[make_pair(0, pos)]);
};

void ColVector::setValue(uint pos, float value) {
	vector[make_pair(pos, 0)] = value;
};
