// NeuralNetwork.hpp 
#include <iostream> 
#include <vector> 
using namespace std; 


// use typedefs for future ease for changing data types like : float to double 

// neural network implementation class! 
class Matrix { 
public: 
	// constructor 
	Matrix(vector<vector<float>>); 
    Matrix(uint x, uint y);

    int setRandom();

	vector<vector<float>> matrix; // stores the different layers of out network 
}; 
