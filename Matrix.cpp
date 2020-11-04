// NeuralNetwork.hpp 
#include <iostream> 
#include <vector> 
#include "Matrix.hpp"
#include <stdlib.h>     /* srand, rand */

using namespace std; 

// use typedefs for future ease for changing data types like : float to double 

Matrix::Matrix(std::vector<vector<float>> matrix) 
{
	this->matrix = matrix;
};

Matrix::Matrix(uint x, uint y) 
{
    std::vector<vector<float>> vectorX;
    for (size_t i = 0; i < x; i++)
    {
       std::vector<float> vectorY (y, 0.0f);
       vectorX.push_back(vectorY);
    }
    
	this->matrix = vectorX;
};

int Matrix::setRandom() {
    for (size_t i = 0; i < matrix.size(); i++)
    {
        /* code */
        vector v = matrix.at(i);
        std::generate(v.begin(), v.end(), std::rand);
    }
    
};
