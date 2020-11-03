// NeuralNetwork.hpp 
#include <iostream> 
#include <vector> 

// use typedefs for future ease for changing data types like : float to double 
typedef float Scalar; 

// neural network implementation class! 
class RowVector { 
public: 
	// constructor 
	RowVector(std::vector<Scalar>); 
    
	std::vector<Scalar> vector; // stores the different layers of out network 
}; 
