#pragma once
// NeuralNetwork.hpp 
#include <iostream> 
#include <vector> 
#include "Matrix.hpp"


// neural network implementation class! 
template <typename T>
class ColVector { 
    std::vector<T>& _vector; // stores the different layers of out network 
public: 
	// constructor 
    // use typedefs for future ease for changing data types like : float to double 

    ColVector(const std::vector<float> vector) : _vector(vector, vector.size()) {}

    ColVector(uint size) : _vector(std::vector(size, T()), size) {}

    float coeffRef(const uint pos) const 
    {
	    return _vector[pos];
    }

    void setValue(const uint pos, const T value) {
        _vector[pos] = value;
    }
    
}; 
