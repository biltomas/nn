// NeuralNetwork.hpp 
#include <iostream> 
#include <vector> 
#include "Matrix.hpp"
#include <stdlib.h>     /* srand, rand */
#include "matrix/matrix.cpp"

using namespace std; 

// use typedefs for future ease for changing data types like : float to double 

Matrix::Matrix(std::vector<vector<float>> data) 
{
	this->data = matrix(data.size(), data[0].size());
};

Matrix::Matrix(uint x, uint y) 
{
    std::vector<float> vector(x*y, 0.0f);
    
	this->data = matrix(vector, x);
};

int Matrix::setRandom() {
    std::pair size = data.size();

    vector<float> v(size.first * size.second);
    std::generate(v.begin(), v.end(), std::rand);
    data = matrix(v, size.first);
    
};
