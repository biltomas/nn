#include "NeuralNetwork.hpp"
#include <iostream> 
#include <vector> 
#include "Matrix.hpp"
#include "ColVector.hpp"
// #include "RowVector.hpp"
#include "RowVector.cpp"
using namespace std; 

// constructor of neural network class 
NeuralNetwork::NeuralNetwork(vector<uint> topology, float learningRate) 
{ 
	this->topology = topology; 
	this->learningRate = learningRate; 
	for (uint i = 0; i < topology.size(); i++) { 
		// initialze neuron layers 
        vector<float> vector;

		if (i == topology.size() - 1) 
			neuronLayers.push_back(new RowVector(topology[i])); 
		else
			neuronLayers.push_back(new RowVector(topology[i] + 1.0f)); 

		// initialize cache and delta vectors 
		cacheLayers.push_back(new RowVector(neuronLayers.size())); 
		deltas.push_back(new RowVector(neuronLayers.size())); 

		// vector.back() gives the handle to recently added element 
		// coeffRef gives the reference of value at that place 
		// (using this as we are using pointers here) 
		if (i != topology.size() - 1) { 
			neuronLayers.back()->vector.at(topology[i]) = 1.0; 
			cacheLayers.back()->vector.at(topology[i]) = 1.0; 
		} 

		// initialze weights matrix 
		if (i > 0) { 
			if (i != topology.size() - 1) { 
				weights.push_back(new Matrix(topology[i - 1] + 1, topology[i] + 1)); 
				weights.back()->setRandom(); 
				weights.back()->matrix.at(topology[i]).assign(weights.back()->matrix.at(topology[i]).size(), 0.0f); 
				weights.back()->matrix.at(topology[i - 1]).at(topology[i]) = 1.0f; 
			} 
			else { 
				weights.push_back(new Matrix(topology[i - 1] + 1, topology[i])); 
				weights.back()->setRandom(); 
			} 
		} 
	} 
}; 
