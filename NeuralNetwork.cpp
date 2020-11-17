#include "NeuralNetwork.hpp"
#include <iostream> 
#include <vector> 
// #include "Matrix.cpp"
// #include "matrix/matrix.hpp"
// #include "ColVector.hpp"
// #include "RowVector.hpp"
// #include "RowVector.cpp"
#include <math.h>
#include <cmath>

using namespace std; 

float activationFunction(float x) 
{ 
    return tanhf(x); 
} 
  
float activationFunctionDerivative(float x) 
{ 
    return 1 - tanhf(x) * tanhf(x); 
} 


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
			neuronLayers.back()->setValue(topology[i], 1.0); 
			cacheLayers.back()->setValue(topology[i], 1.0); 
		} 

		// initialze weights matrix 
		if (i > 0) { 
			if (i != topology.size() - 1) { 
				std::pair<std::size_t, std::size_t> size = weights.back()->data.size();
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
	for (uint i = 0; i < input.vector.size().second; i++) {
		neuronLayers.front()->setValue(i, input.coeffRef(i));
	};
  
    // propagate the data forawrd 
    for (uint i = 1; i < topology.size(); i++) { 
        // already explained above 
        (*neuronLayers[i]) = (*neuronLayers[i - 1]) * (*weights[i - 1]); 
    } ;
  
    // apply the activation function to your network 
    // unaryExpr applies the given function to all elements of CURRENT_LAYER 
    for (uint i = 1; i < topology.size() - 1; i++) { 
		for (uint ii = 0; ii < topology[i]; ii++){
			neuronLayers[i]->setValue(ii, activationFunction(neuronLayers[i]->coeffRef(ii)));
		}
    } 
} 

void NeuralNetwork::calcErrors(RowVector& output) 
{ 
    // calculate the errors made by neurons of last layer 
    (*deltas.back()) = output - (*neuronLayers.back()); 
  
    // error calculation of hidden layers is different 
    // we will begin by the last hidden layer 
    // and we will continue till the first hidden layer 
    for (uint i = topology.size() - 2; i > 0; i--) { 
        (*deltas[i]) = (*deltas[i + 1]) * (weights[i]->transpose()); 
    } 
} 

void NeuralNetwork::updateWeights() 
{ 
    // topology.size()-1 = weights.size() 
    for (uint i = 0; i < topology.size() - 1; i++) { 
        // in this loop we are iterating over the different layers (from first hidden to output layer) 
        // if this layer is the output layer, there is no bias neuron there, number of neurons specified = number of cols 
        // if this layer not the output layer, there is a bias neuron and number of neurons specified = number of cols -1 
        if (i != topology.size() - 2) { 
            for (uint c = 0; c < weights[i]->data.cols() - 1; c++) { 
                for (uint r = 0; r < weights[i]->data.rows(); r++) { 
					float num = weights[i]->coeffRef(r, c) + learningRate * deltas[i + 1]->coeffRef(c) * activationFunctionDerivative(cacheLayers[i + 1]->coeffRef(c)) * neuronLayers[i]->coeffRef(r);
                    weights[i]->setValue(r, c, num); 
                } 
            } 
        } 
        else { 
            for (uint c = 0; c < weights[i]->data.cols(); c++) { 
                for (uint r = 0; r < weights[i]->data.rows(); r++) { 
					float num = weights[i]->coeffRef(r, c) + learningRate * deltas[i + 1]->coeffRef(c) * activationFunctionDerivative(cacheLayers[i + 1]->coeffRef(c)) * neuronLayers[i]->coeffRef(r);
                    weights[i]->setValue(r, c, num);
                } 
            } 
        } 
    } 
} 

void NeuralNetwork::propagateBackward(RowVector& output) 
{ 
    calcErrors(output); 
    updateWeights(); 
} 

void NeuralNetwork::train(std::vector<RowVector*> input_data, std::vector<RowVector*> output_data) 
{ 
    for (uint i = 0; i < input_data.size(); i++) { 
        // std::cout << "Input to neural network is : " << *input_data[i] << std::endl; 
        propagateForward(*input_data[i]); 
        // std::cout << "Expected output is : " << *output_data[i] << std::endl; 
        // std::cout << "Output produced is : " << *neuronLayers.back() << std::endl; 
        propagateBackward(*output_data[i]); 
        std::cout << "MSE : " << std::sqrt((*deltas.back()).dot((*deltas.back())) / deltas.back()->vector.size().second) << std::endl; 
    } 
}
