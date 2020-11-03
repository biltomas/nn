// NeuralNetwork.hpp 
#include <iostream> 
#include <vector> 

// use typedefs for future ease for changing data types like : float to double 
typedef float Scalar; 

// neural network implementation class! 
class Matrix { 
public: 
	// constructor 
	Matrix(std::vector<std::vector<Scalar>>); 

	std::vector<std::vector<Scalar>> matrix; // stores the different layers of out network 
}; 
