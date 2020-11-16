#include "NeuralNetwork.hpp"
#include <iostream> 
#include <vector> 
#include "Matrix.hpp"
#include "matrix/matrix.cpp"
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
			neuronLayers.back()->vector[topology[i]] = 1.0; 
			cacheLayers.back()->vector[topology[i]] = 1.0; 
		} 

		// initialze weights matrix 
		if (i > 0) { 
			if (i != topology.size() - 1) { 
				std::pair size = weights.back()->data.size();
				std::vector<float> vector(size.first * size.second, 0.0f);
				for (int x = 0; x< size.first; x++) {
					vector[x + size.second - 1] = 1.0f;
				}
				weights.push_back(new Matrix(topology[i - 1] + 1, topology[i] + 1)); 
				weights.back()->data = matrix(vector, size.first);
			} 
			else { 
				weights.push_back(new Matrix(topology[i - 1] + 1, topology[i])); 
				weights.back()->setRandom(); 
			} 
		} 
	} 
}; 

void NeuralNetwork::propagateForward(RowVector& input) 
{ 
    // set the input to input layer 
    // block returns a part of the given vector or matrix 
    // block takes 4 arguments : startRow, startCol, blockRows, blockCols 
    *neuronLayers.front() = input; 
  
    // propagate the data forawrd 
    for (uint layerIdx = 1; layerIdx < topology.size(); layerIdx++) { 
        // already explained above
		for (uint currPerc = 0; currPerc < neuronLayers[layerIdx]->vector.size(); currPerc++){
			neuronLayers[layerIdx][currPerc] = 0;
			for (uint prevPerc = 0; prevPerc < neuronLayers[layerIdx - 1]->vector.size(); prevPerc++){
				float weight = weights[layerIdx-1][currPerc];
				neuronLayers[layerIdx][currPerc] += neuronLayers[layerIdx - 1][prevPerc] * weights[layerIdx - 1][currPerc];
			}
			(neuronLayers[layerIdx][currPerc]) = (*neuronLayers[layerIdx - 1]) * (*weights[layerIdx - 1]); 
		}
    } 
  
    // apply the activation function to your network 
    // unaryExpr applies the given function to all elements of CURRENT_LAYER 
    for (uint i = 1; i < topology.size() - 1; i++) { 
        neuronLayers[i]->block(0, 0, 1, topology[i]).unaryExpr(std::ptr_fun(activationFunction)); 
    } 
}
