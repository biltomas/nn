// NeuralNetwork.hpp 
#include <iostream> 
#include <vector> 
using namespace std; 

// use typedefs for future ease for changing data types like : float to double 

// neural network implementation class! 
class ColVector { 
public: 
	// constructor 
	ColVector(vector<float>); 
    
	vector<float> vector; // stores the different layers of out network 
}; 
