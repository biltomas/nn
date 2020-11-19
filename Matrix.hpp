#pragma once
#include <iostream> 
#include <vector> 
#include "matrix/matrix.hpp"
using namespace std; 


// matrix implementation class! 
class Matrix { 
public: 
	// constructor 
	Matrix(vector<vector<float> > data); 
    Matrix(uint x, uint y);
	Matrix(std::vector<float> data, uint rows);
	float coeffRef(uint pos1, uint pos2);
	void setValue(uint pos1, uint pos2, float value);

    int setRandom();

	matrix data; // stores matrix data
	Matrix transpose();
}; 
