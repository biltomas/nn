#pragma once
// NeuralNetwork.hpp 
#include <iostream> 
#include <vector> 
#include "matrix/matrix.hpp"

using namespace std; 

// use typedefs for future ease for changing data types like : float to double 

// neural network implementation class! 
class ColVector { 
public: 
	// constructor 
	ColVector(vector<float> vector);
    ColVector(uint size);
	void setValue(uint pos, float value);
    float coeffRef(uint pos);
    
	matrix vector; // stores the different layers of out network 
}; 
