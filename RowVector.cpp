// NeuralNetwork.hpp 
#include <iostream> 
#include <vector> 
#include "RowVector.hpp"

using namespace std; 

// use typedefs for future ease for changing data types like : float to double 

RowVector::RowVector(std::vector<float> vector) 
{
	this->vector = vector;
};

RowVector::RowVector(uint size) 
{
	std::vector<float> vector (size, 0.0f);
	this->vector = vector;
};

float RowVector::coeffRef(uint pos) 
{
	return(vector.at(pos));
};
