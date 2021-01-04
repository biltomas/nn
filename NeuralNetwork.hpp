#pragma once
#include "Matrix.hpp"
#include "RowVector.hpp"
#include <iostream> 
#include <vector> 
#include <memory>
#include <math.h>
#include <fstream>

using std::unique_ptr;
using std::make_unique;

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
	void train(std::vector<RowVector<float>*> data, std::vector<RowVector<float>>& output_data, float learningRate = float(0.005)); 

	// function to predict output given input and output into file
	void predict(std::vector<RowVector<float>*> data, string outputFile = "./outputFile.csv");

	// function of validation
    void validate(std::vector<RowVector<float>*>& data, std::vector<RowVector<float>>& labels);

	// stores value of every neuron
	std::vector<unique_ptr<RowVector<float>>> neuronLayers;
	// stores inner potentials of neurons
	std::vector<unique_ptr<RowVector<float>>> cacheLayers;
	// stores errors for every neuron
	std::vector<unique_ptr<RowVector<float>>> deltas;
	// input weights of neurons
	std::vector<unique_ptr<Matrix<float>>> weights; 
	// stores number of neurons in every layer
	std::vector<size_t> topology;
	// activation gradients for crossentropy backpropagation
	std::vector<float> activation_gradients_;
	// stores momentum of neuron gradients
    std::vector<unique_ptr<Matrix<float>>> momentum;
	float learningRate; 
}; 
