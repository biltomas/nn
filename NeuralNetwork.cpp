#include "NeuralNetwork.hpp"
// #include "Matrix.cpp"
// #include "matrix/matrix.hpp"
// #include "ColVector.hpp"
// #include "RowVector.hpp"
// #include "RowVector.cpp"


float crossEntropyLoss(const RowVector<float>& probabilities, const RowVector<float>& labels) {
    float total = 0;
    for (size_t i = 0; i < probabilities.length(); i++) {
        total += labels.coeffRef(i) * log(probabilities.coeffRef(i)) +
                (1 - labels.coeffRef(i)) * log(1 - probabilities.coeffRef(i));
    }
    return -1 * total;
}

void softMax (const RowVector<float>& input, RowVector<float>& output) {
	size_t i;
	float m, sum, constant;

	m = -10;
	for (i = 0; i < input.length(); ++i) {
		if (m < input.coeffRef(i)) {
			m = input.coeffRef(i);
		}
	}

	sum = 0.0;
	for (i = 0; i < input.length(); ++i) {
		sum += exp(input.coeffRef(i) - m);
	}

	constant = m + log(sum);
	for (i = 0; i < input.length(); ++i) {
		output.setValue(i, exp(input.coeffRef(i) - constant));
	}
}

void softMaxDerivative (const RowVector<float>& input, RowVector<float>& output) {
	vector<float> exps;
	size_t i;
	float sum = 0.0;
	float temp;
	for (i = 0; i < input.length(); ++i) {
		temp = exp(input.coeffRef(i));
		sum += temp;
		exps.push_back(temp);
	}
	sum *= sum;
	for (i = 0; i <input.length(); i++) {
		float tempSum = 0.0f;
		for (size_t ii = 0; ii < input.length(); ii++) {
			if (ii != i) {
				tempSum += input.coeffRef(ii);
			}
		}
		temp = exps[i] * tempSum / sum;
		output.setValue(i, temp);
	}
}

float outputActivationFunction(float x) 
{ 
    // return 0.5*tanhf(x) + 0.5; 
	return 1 / (1 + exp(-x));
	// return x < 0 ? 0.01*x : x;
} 
  
float outputActivationFunctionDerivative(float x) 
{ 
    // return 0.5*(1 - tanhf(x) * tanhf(x)); 
	return outputActivationFunction(x) * (1 - outputActivationFunction(x));
	// return x < 0 ? 0.01 : 1;
}


float activationFunction(float x) 
{ 
    // return 0.5*tanhf(x) + 0.5; 
	// return 1 / ( 1+exp(-x));
	// return x < 0 ? 0.01*x : x;
	return x < 0 ? 0.24 * (exp(x) - 1) : x;
} 
  
float activationFunctionDerivative(float x) 
{ 
    // return 0.5*(1 - tanhf(x) * tanhf(x)); 
	// return activationFunction(x)*(1-activationFunction(x));
	// return x < 0 ? 0.01 : 1;
	return x < 0 ? activationFunction(x)+0.24 : 1;
}


// constructor of neural network class 
NeuralNetwork::NeuralNetwork(vector<size_t> topology, float learningRate) 
{
	// printf("NeuralNetwork start\n");
	this->topology = topology; 
	this->learningRate = learningRate; 
	for (size_t i = 0; i < topology[topology.size()-1]; i++) {
		activation_gradients_.push_back(0.0f);
	}
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
				// cout << weights.size() << endl;
				std::size_t size = (topology[i - 1] + 1) * (topology[i] + 1);
				std::vector<float> vector(size, 0.0f);
				
				// printf("NeuralNetwork 2.2\n");
				weights.push_back(make_unique<Matrix<float>>(vector, topology[i - 1] + 1));
				momentum.push_back(make_unique<Matrix<float>>(vector, topology[i - 1] + 1)); 
				// printf("NeuralNetwork 2.3\n");
				// weights.back()->setRandom();
				// weights.back()->setNumber(sqrt(3/284));
				weights.back()->setNumber(sqrt(1/2840));
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
				std::size_t size = (topology[i - 1] + 1) * (topology[i]);
				std::vector<float> vector(size, 0.0f);
				momentum.push_back(make_unique<Matrix<float>>(vector, topology[i - 1] + 1)); 
				// weights.back()->setRandom(); 
				// weights.back()->setNumber(sqrt(3/1024));
				weights.back()->setNumber(sqrt(1/1024));
			} 
		} 
	}
	// printf("NeuralNetwork end\n");
}

void NeuralNetwork::propagateForward(RowVector<float>& input) 
{ 
	// printf("propagateForward start\n");
    // set the input to input layer 
    // block returns a part of the given vector or matrix 
    // block takes 4 arguments : startRow, startCol, blockRows, blockCols 
	for (uint i = 0; i < input.length(); i++) {
		neuronLayers.front()->setValue(i, input.coeffRef(i));
		cacheLayers.front()->setValue(i, input.coeffRef(i));
	};
	// printf("front set\n");
  
    // propagate the data forawrd 
    for (uint i = 1; i < topology.size(); i++) { 
        // already explained above 
		// cout << "layer: " << i << endl;
		// cout << "updating cache" << endl;
        matmul(*neuronLayers[i - 1], *weights[i - 1], *neuronLayers[i]);
		*cacheLayers[i] = *neuronLayers[i];
		// cout << "neuron " << i-1 << ": " << (*neuronLayers[i - 1]) << endl; 
		// cout << "cache " << i << ": " << *cacheLayers[i] << endl;
		if (i != topology.size() - 1) {
			for (uint ii = 0; ii < neuronLayers[i]->length(); ii++) {
				neuronLayers[i]->setValue(ii, activationFunction(neuronLayers[i]->coeffRef(ii)));
			}
		} else {
			softMax(*neuronLayers[i], *neuronLayers[i]);
			// for (uint ii = 0; ii < neuronLayers[i]->length(); ii++) {
			// 	neuronLayers[i]->setValue(ii, outputActivationFunction(neuronLayers[i]->coeffRef(ii)));
			// }
		}
		
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
    // for (uint i = 1; i < topology.size(); i++) { 
	// 	for (uint ii = 0; ii < topology[i]; ii++){
	// 		// cout << i << " " << ii << " inner potential:" << activationFunction(cacheLayers[i]->coeffRef(ii)) << endl;
	// 		neuronLayers[i]->setValue(ii, activationFunction(neuronLayers[i]->coeffRef(ii)));
	// 	}
    // } 
	// for (uint ii = 0; ii < topology[topology.size() - 1]; ii++) {
	// 	neuronLayers[topology.size() - 1]->setValue(ii, cacheLayers[topology.size() - 1]->coeffRef(ii));
	// }
	// printf("propagateForward end\n");
} 

void NeuralNetwork::calcErrors(RowVector<float>& output) 
{ 
	// printf("calcErrors start\n");
    // calculate the errors made by neurons of last layer 
	// cout << "delats size" << deltas.size() << endl;
	// cout << "output size: " << output.length() << endl;
	// cout << "pred size: " << (*neuronLayers.back()).length() << endl;
	
	for (size_t i = 0; i < output.length(); i++) {
		deltas.back()->setValue(i, output.coeffRef(i) - neuronLayers.back()->coeffRef(i));
		// cout << output.coeffRef(i) << " - " << neuronLayers.back()->coeffRef(i) << " = " << deltas.back()->coeffRef(i) << endl;
		// cout << deltas.back()->coeffRef(i) << " ";
	}
	for (size_t i = 0; i < output.length(); i++) {
		float activation_grad = 0.0f;
		for (size_t j = 0; j < output.length(); j++) {
			// cout << i << j << endl;
			if (i == j) {
                    activation_grad += neuronLayers.back()->coeffRef(i) * (1.0f - neuronLayers.back()->coeffRef(i))* deltas.back()->coeffRef(j);
                } else {
                    activation_grad += -neuronLayers.back()->coeffRef(i) * neuronLayers.back()->coeffRef(j) * deltas.back()->coeffRef(j);
                }
		}
		// cout << "activation_grad " << activation_grad << endl;
		activation_gradients_[i] = activation_grad;
	}

	// cout << endl;
    // (*deltas.back()) = output - (*neuronLayers.back()); 
	// cout << "deltas.back size: " << (*deltas.back()).length() << endl;
	// printf("first\n");
  
    // error calculation of hidden layers is different 
    // we will begin by the last hidden layer 
    // and we will continue till the first hidden layer 
    for (uint i = topology.size() - 2; i > 0; i--) { 
        matmul(*deltas[i + 1], weights[i]->transpose(), *deltas[i]); 
        weights[i]->transpose(); //now transpose it back
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
		// cout << endl << endl; 
        if (i != topology.size() - 2) { 
            #pragma omp parallel for num_threads(8)
            for (unsigned c = 0; c < weights[i]->cols() - 1; c++) { 
                for (uint r = 0; r < weights[i]->rows(); r++) { 
					float num = weights[i]->operator[]({r, c}) + learningRate * 
                        deltas[i + 1]->coeffRef(c) * 
                        activationFunctionDerivative(cacheLayers[i + 1]->coeffRef(c)) * 
                        neuronLayers[i]->coeffRef(r);
					// cout << "derivative: " << activationFunctionDerivative(cacheLayers[i + 1]->coeffRef(c)) << " of " << cacheLayers[i + 1]->coeffRef(c) << endl;
					// cout << "weight: " << i << " pos: " << r << ", " << c << " prev value " << weights[i]->coeffRef(r, c) << " set to " << num << endl;
                    weights[i]->operator[]({r, c}) = (num + momentum[i]->operator[]({r, c}))/2;
                    momentum[i]->operator[]({r, c}) = (num + momentum[i]->operator[]({r, c}))/2; 
                } 
            } 
        } 
        else { 
            #pragma omp parallel for num_threads(8)
            for (uint c = 0; c < weights[i]->cols(); c++) { 
                for (uint r = 0; r < weights[i]->rows(); r++) { 
					// float num = weights[i]->operator[]({r, c}) + learningRate * deltas[i + 1]->coeffRef(c) * (cacheLayers[i + 1]->coeffRef(c)) * neuronLayers[i]->coeffRef(r);
					float num = weights[i]->operator[]({r, c}) + learningRate * neuronLayers[topology.size()-2]->coeffRef(r) * activation_gradients_[c];
					// cout << "derivative: " << activationFunctionDerivative(cacheLayers[i + 1]->coeffRef(c)) << " of " << cacheLayers[i + 1]->coeffRef(c) << endl;
					// cout << "weight: " << i << " pos: " << r << ", " << c << " prev value " << weights[i]->operator[]({r, c}) << " set to " << num << endl;
                    weights[i]->operator[]({r, c}) = (num + momentum[i]->operator[]({r, c}))/2;
                    momentum[i]->operator[]({r, c}) = (num + momentum[i]->operator[]({r, c}))/2; 
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

void NeuralNetwork::train(std::vector<RowVector<float>*> input_data, std::vector<RowVector<float>>& output_data, float learningRate) 
{ 
	// printf("train start\n");
	this->learningRate = learningRate; 
    for (uint i = 0; i < input_data.size(); i++) { 
        if (i % 1000 == 0) {
            std::cout << "Checkpoint: " << i + 1 << " out of 59,000" << std::endl;
        }
        // std::cout << "Input to neural network is : " << input_data[i] << std::endl; 
        propagateForward(*input_data[i]); 
        // std::cout << "Expected output is : "; 
		// for(unsigned ii=0; ii< neuronLayers.back()->length(); ++ii)
  		// 	std::cout << output_data[i].coeffRef(ii) << ' '; 
        // std::cout << endl << "Output produced is : ";
		// for(unsigned i=0; i< neuronLayers.back()->length(); ++i)
  		// 	std::cout << neuronLayers.back()->coeffRef(i) << ' '; 
        propagateBackward(output_data[i]); 
        // std::cout << "\nMSE : " << std::sqrt((*deltas.back()).dot((*deltas.back())) / deltas.back()->length()) << std::endl; 
    } 
	// printf("train end\n");
}

void NeuralNetwork::predict(std::vector<RowVector<float>*> data, string outputFile)
{
	ofstream myfile;
	myfile.open(outputFile);
	
	for (uint i = 0; i < data.size(); i++) { 
        //std::cout << "Input to neural network is : " << input_data[i] << std::endl; 
        propagateForward(*data[i]); 
		float max = 0;
		int result = 0;
		for (uint ii = 0; ii < neuronLayers.back()->length(); ii++) {
			if (neuronLayers.back()->coeffRef(ii) > max) {
				max = neuronLayers.back()->coeffRef(ii);
				result = ii;
			}
		};
		std::cout << endl << "Output produced is : ";
		for(unsigned i=0; i< neuronLayers.back()->length(); ++i)
  			std::cout << neuronLayers.back()->coeffRef(i) << ' '; 
		cout << endl << "predicted label is: " << result << endl;
		myfile << result << "\n";
    } 
	myfile.close();
}

void NeuralNetwork::validate(std::vector<RowVector<float>*>& data, std::vector<RowVector<float>>& labels) {
    float loss = 0;
    float correct_predictions = 0;
    for (size_t current = 0; current < data.size(); current++) {
        auto& image = *data[current];
        auto& label = labels[current];
        propagateForward(image);
        const std::vector<float>& probabilities = neuronLayers.back()->data().to_vector();
        const std::vector<float>& label_vector = label.data().to_vector();
        //compute the loss
        loss += crossEntropyLoss(*neuronLayers.back(), label);
        unsigned top_class = std::distance(probabilities.begin(), std::max_element(probabilities.begin(),
                                            probabilities.end()));
        unsigned expected_class = std::distance(label_vector.begin(), std::max_element(label_vector.begin(),
                                            label_vector.end()));
        if (top_class == expected_class) {
            correct_predictions++;
        }
    }
    std::cout << "Validation accuracy: " << correct_predictions / data.size() << std::endl;
    std::cout << "Average validation loss: " << loss / data.size() << std::endl;
}
