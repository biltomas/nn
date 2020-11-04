// NeuralNetwork.hpp 
#include <iostream> 
#include <vector> 
using namespace std; 

// use typedefs for future ease for changing data types like : float to double 

// neural network implementation class! 
class RowVector { 
public: 
	// constructor 
	RowVector(vector<float> vector);
    RowVector(uint size);
    float coeffRef(uint pos);

	vector<float> vector;
};
