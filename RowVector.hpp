// NeuralNetwork.hpp 
#include <iostream> 
#include <vector> 
#include "matrix/matrix.hpp"
using namespace std; 

// use typedefs for future ease for changing data types like : float to double 

// neural network implementation class! 
class RowVector { 
public: 
	// constructor 
	RowVector(vector<float> vector);
    RowVector(uint size);
	void setValue(uint pos, float value);
    float coeffRef(uint pos);

	matrix vector;
};
RowVector operator*(RowVector m1, const Matrix& m2);
