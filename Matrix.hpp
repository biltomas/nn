// NeuralNetwork.hpp 
#include <iostream> 
#include <vector> 
#include "matrix/matrix.hpp"
using namespace std; 


// use typedefs for future ease for changing data types like : float to double 

// neural network implementation class! 
class Matrix { 
public: 
	// constructor 
	Matrix(vector<vector<float> > data); 
    Matrix(uint x, uint y);
	Matrix(std::vector<float> data, uint rows);
	float coeffRef(uint pos1, uint pos2);
	void setValue(uint pos1, uint pos2, float value);

    int setRandom();

	matrix data; // stores the different layers of out network 
	Matrix transpose();
}; 
