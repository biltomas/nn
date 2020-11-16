// NeuralNetwork.hpp 
#include <iostream> 
#include <vector> 
#include "RowVector.hpp"

using namespace std; 

// use typedefs for future ease for changing data types like : float to double 

RowVector::RowVector(std::vector<float> vector) 
{
	this->vector = matrix(vector, 1);
};

RowVector::RowVector(uint size) 
{
	std::vector<float> vector (size, 0.0f);
	this->vector = matrix(vector, 1);
};

float RowVector::coeffRef(uint pos) 
{
	return(vector[make_pair(0, pos)]);
};

void RowVector::setValue(uint pos, float value) {
	vector[make_pair(0, pos)] = value;
};

RowVector operator*(RowVector m1, const Matrix& m2) {
	return RowVector((m1.vector * m2.data).to_vector());
};
