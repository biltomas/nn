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

Matrix::Matrix(std::vector<float> data, uint rows) 
{
	this->data = matrix(data, rows);
};

Matrix::Matrix(uint x, uint y) 
{
    std::vector<float> vector(x*y, 0.0f);
    
	this->data = matrix(vector, x);
};

void Matrix::setRandom() {
    std::pair<size_t, size_t> size = data.size();

    vector<float> v(size.first * size.second);
    // std::generate(v.begin(), v.end(), (float) std::rand/RAND_MAX);
    for (int i = 0; i < v.size(); i++) {
        v[i] = (float) std::rand()/RAND_MAX;
        // cout << v[i] << " ";
    }
    data = matrix(v, size.first);
    
};

Matrix Matrix::transpose() {
    return Matrix(data.transpose().to_vector(), data.size().second);
};

float Matrix::coeffRef(uint pos1, uint pos2) 
{
	return(data[make_pair(pos1, pos2)]);
};

void Matrix::setValue(uint pos1, uint pos2, float value) {
	data[make_pair(pos1, pos2)] = value;
};