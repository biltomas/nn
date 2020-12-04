#pragma once
#include "Matrix.hpp"
#include "RowVector.hpp"
#include <iostream> 
#include <vector> 
#include <memory>
#include <math.h>

// use typedefs for future ease for changing data types like : float to double 
// typedef Matrix Matrix; 
// typedef RowVector RowVector; 
// typedef ColVector ColVector; 

using std::unique_ptr;
using std::make_unique;

// neural network implementation class! 
class NeuralNetwork { 
public: 
	// constructor 
	NeuralNetwork(vector<size_t> topology, float learningRate = float(0.005)); 

	// function for forward propagation of data 
	void propagateForward(RowVector<float>& input); 

	// function for backward propagation of errors made by neurons 
	void propagateBackward(RowVector<float>& output); 

	// function to calculate errors made by neurons in each layer 
	void calcErrors(RowVector<float>& output); 

	// function to update the weights of connections 
	void updateWeights(); 

	// function to train the neural network give an array of data points 
	void train(std::vector<RowVector<float>*> data, std::vector<RowVector<float>*> output_data); 

	// storage objects for working of neural network 
	/* 
		use pointers when using std::vector<Class> as std::vector<Class> calls destructor of 
		Class as soon as it is pushed back! when we use pointers it can't do that, besides 
		it also makes our neural network class less heavy!! It would be nice if you can use 
		smart pointers instead of usual ones like this 
		*/
	std::vector<unique_ptr<RowVector<float>>> neuronLayers; // stores the different layers of out network 
	std::vector<unique_ptr<RowVector<float>>> cacheLayers; // stores the unactivated (activation fn not yet applied) values of layers 
	std::vector<unique_ptr<RowVector<float>>> deltas; // stores the error contribution of each neurons 
	std::vector<unique_ptr<Matrix<float>>> weights; // the connection weights itself 
	std::vector<size_t> topology;
	float learningRate; 
}; 
