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
	// cout << vector.size().first << ", " << vector.size().second << endl;
	vector[make_pair(0, pos)] = value;
};
float RowVector::dot(RowVector vector2) {
	float result = 0;
	for (int i = 0; i < vector.size().second; i++) {
		result += coeffRef(i) * vector2.coeffRef(i);
	}
	return result;
};

RowVector operator*(RowVector m1, const Matrix& m2) {
	return RowVector((m1.vector * m2.data).to_vector());
};

RowVector operator-(RowVector m1, const RowVector& m2) {
	return RowVector((m1.vector - m2.vector).to_vector());
};
