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
	Matrix(vector<vector<float>> data); 
    Matrix(uint x, uint y);

    int setRandom();

	matrix data; // stores the different layers of out network 
}; 
