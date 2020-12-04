#include "NeuralNetwork.hpp"
// #include "Matrix.cpp"
// #include "matrix/matrix.hpp"
// #include "ColVector.hpp"
// #include "RowVector.hpp"
// #include "RowVector.cpp"



float activationFunction(float x) 
{ 
    return tanhf(x); 
} 
  
float activationFunctionDerivative(float x) 
{ 
    return 1 - tanhf(x) * tanhf(x); 
} 


// constructor of neural network class 
NeuralNetwork::NeuralNetwork(vector<size_t> topology, float learningRate) 
{ 
	// printf("NeuralNetwork start\n");
	this->topology = topology; 
	this->learningRate = learningRate; 
	for (uint i = 0; i < topology.size(); i++) { 
		// initialze neuron layers 
        std::vector<float> vector;

		if (i == topology.size() - 1) 
			neuronLayers.push_back(make_unique<RowVector<float>>(topology[i])); 
		else
			neuronLayers.push_back(make_unique<RowVector<float>>(topology[i] + 1)); 

		// initialize cache and delta vectors 
		cacheLayers.push_back(make_unique<RowVector<float>>(neuronLayers.back()->length())); 
		deltas.push_back(make_unique<RowVector<float>>(neuronLayers.back()->length())); 
		// printf("NeuralNetwork 1\n");

		// vector.back() gives the handle to recently added element 
		// coeffRef gives the reference of value at that place 
		// (using this as we are using pointers here) 
		if (i != topology.size() - 1) { 
			// printf("NeuralNetwork 1.2\n");
			neuronLayers.back()->setValue(topology[i], 1.0); 
			// printf("NeuralNetwork 1.5\n");
			cacheLayers.back()->setValue(topology[i], 1.0); 
		} 
		// printf("NeuralNetwork 2\n");

		// initialze weights matrix 
		if (i > 0) { 
			if (i != topology.size() - 1) { 
				// printf("NeuralNetwork 2.1\n");
				cout << weights.size() << endl;
				std::size_t size = (topology[i - 1] + 1) * (topology[i] + 1);
				std::vector<float> vector(size, 0.0f);
				
				// printf("NeuralNetwork 2.2\n");
				weights.push_back(make_unique<Matrix<float>>(vector, topology[i - 1] + 1)); 
				// printf("NeuralNetwork 2.3\n");
				weights.back()->setRandom();
				// printf("NeuralNetwork 2.4\n");
                auto fst = topology[i - 1];
                auto snd = topology[i];
				weights.back()->operator[]({fst, snd}) = 1.0f;
				// printf("NeuralNetwork 2.5\n");
				// cout << weights.back()->data.size().first << " , " << weights.back()->data.size().second << endl;
				for (unsigned x = 0; x < (topology[i - 1] + 1); x++) {
					// cout << topology[i - 1] + 1 << ", " << x << endl;
                    Matrix<float>& current_weigths = *weights.back();
                    current_weigths[{x, topology[i]}] = 0.0f;
                    //weights.back()->setValue(x, topology[i], 0.0f);
					// printf("NeuralNetwork 2.6\n");
				}
			} 
			else { 
				// printf("NeuralNetwork 2.3\n");
				weights.push_back(make_unique<Matrix<float>>(topology[i - 1] + 1, topology[i])); 
				weights.back()->setRandom(); 
			} 
		} 
	}
	// printf("NeuralNetwork end\n");
}; 

void NeuralNetwork::propagateForward(RowVector<float>& input) 
{ 
	// printf("propagateForward start\n");
    // set the input to input layer 
    // block returns a part of the given vector or matrix 
    // block takes 4 arguments : startRow, startCol, blockRows, blockCols 
	for (uint i = 0; i < input.length(); i++) {
		cacheLayers.front()->setValue(i, input.coeffRef(i));
	};
	// printf("front set\n");
  
    // propagate the data forawrd 
    for (uint i = 1; i < topology.size(); i++) { 
        // already explained above 
		// cout << "layer: " << i << endl;
		// cout << "updating cache" << endl;
        (*cacheLayers[i]) = (*cacheLayers[i - 1]) * (*weights[i - 1]);
		// cout << "cache " << i << " size " << cacheLayers[i]->vector.cols() << endl;
		// cout << "cache " << i << ": ";
		// for (uint ii = 0; ii < cacheLayers[i]->vector.cols(); ii++){
		// 	cout << cacheLayers[i]->coeffRef(ii) << " ";
		// } 
		// cout << endl;
		// cout << "cacheLayers[i] len: " << (*cacheLayers[i]).vector.size().second << endl;
    };

	// printf("propagated\n");
  
    // apply the activation function to your network 
    // unaryExpr applies the given function to all elements of CURRENT_LAYER 
	// cout << "topology size: " << topology.size() << endl;
    for (uint i = 1; i < topology.size() - 1; i++) { 
		for (uint ii = 0; ii < topology[i]; ii++){
			// cout << i << " " << ii << " inner potential:" << activationFunction(cacheLayers[i]->coeffRef(ii)) << endl;
			neuronLayers[i]->setValue(ii, activationFunction(cacheLayers[i]->coeffRef(ii)));
		}
    } 
	for (uint ii = 0; ii < topology[topology.size() - 1]; ii++) {
		neuronLayers[topology.size() - 1]->setValue(ii, cacheLayers[topology.size() - 1]->coeffRef(ii));
	}
	// printf("propagateForward end\n");
} 

void NeuralNetwork::calcErrors(RowVector<float>& output) 
{ 
	// printf("calcErrors start\n");
    // calculate the errors made by neurons of last layer 
	// cout << "delats size" << deltas.size() << endl;
	// cout << "output size: " << output.vector.size().second << endl;
	// cout << "pred size: " << (*neuronLayers.back()).vector.size().second << endl;
    (*deltas.back()) = output - (*neuronLayers.back()); 
	// printf("first\n");
  
    // error calculation of hidden layers is different 
    // we will begin by the last hidden layer 
    // and we will continue till the first hidden layer 
    for (uint i = topology.size() - 2; i > 0; i--) { 
        (*deltas[i]) = (*deltas[i + 1]) * (weights[i]->transpose()); 
    } 
	// printf("calcErrors end");
} 

void NeuralNetwork::updateWeights() 
{ 
	// printf("updateWeights start\n");
    // topology.size()-1 = weights.size() 
    for (uint i = 0; i < topology.size() - 1; i++) { 
        // in this loop we are iterating over the different layers (from first hidden to output layer) 
        // if this layer is the output layer, there is no bias neuron there, number of neurons specified = number of cols 
        // if this layer not the output layer, there is a bias neuron and number of neurons specified = number of cols -1 
        if (i != topology.size() - 2) { 
            for (unsigned c = 0; c < weights[i]->cols() - 1; c++) { 
                for (uint r = 0; r < weights[i]->rows(); r++) { 
					float num = weights[i]->operator[]({r, c}) + learningRate * 
                        deltas[i + 1]->coeffRef(c) * 
                        activationFunctionDerivative(cacheLayers[i + 1]->coeffRef(c)) * 
                        neuronLayers[i]->coeffRef(r);
					// cout << "weight: " << i << " pos: " << r << ", " << c << " prev value " << weights[i]->coeffRef(r, c) << " set to " << num << endl;
                    weights[i]->operator[]({r, c}) = num; 
                } 
            } 
        } 
        else { 
            for (uint c = 0; c < weights[i]->cols(); c++) { 
                for (uint r = 0; r < weights[i]->rows(); r++) { 
					float num = weights[i]->operator[]({r, c}) + learningRate * deltas[i + 1]->coeffRef(c) * activationFunctionDerivative(cacheLayers[i + 1]->coeffRef(c)) * neuronLayers[i]->coeffRef(r);
					// cout << "weight: " << i << " pos: " << r << ", " << c << " prev value " << weights[i]->coeffRef(r, c) << " set to " << num << endl;
                    weights[i]->operator[]({r, c}) = num;
                } 
            } 
        } 
    } 
	// printf("updateWeights end\n");
} 

void NeuralNetwork::propagateBackward(RowVector<float>& output) 
{ 
    calcErrors(output); 
    updateWeights(); 
} 

void NeuralNetwork::train(std::vector<RowVector<float>*> input_data, std::vector<RowVector<float>*> output_data) 
{ 
	// printf("train start\n");
    for (uint i = 0; i < input_data.size(); i++) { 
        //std::cout << "Input to neural network is : " << input_data[i] << std::endl; 
        propagateForward(*input_data[i]); 
        //std::cout << "Expected output is : " << output_data[i]->coeffRef(0) << std::endl; 
        //std::cout << "Output produced is : ";
		// for(int i=0; i< neuronLayers.back()->vector.rows(); ++i)
  		// 	std::cout << neuronLayers.back()->coeffRef(i) << ' '; 
        propagateBackward(*output_data[i]); 
        std::cout << "\nMSE : " << std::sqrt((*deltas.back()).dot((*deltas.back())) / deltas.back()->length()) << std::endl; 
    } 
	// printf("train end\n");
}
