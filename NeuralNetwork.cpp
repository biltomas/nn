#include "NeuralNetwork.hpp"


float crossEntropyLoss(const RowVector<float>& probabilities, const RowVector<float>& labels) {
    float total = 0;
    for (size_t i = 0; i < probabilities.length(); i++) {
        total += labels.coeffRef(i) * log(probabilities.coeffRef(i)) +
                (1 - labels.coeffRef(i)) * log(1 - probabilities.coeffRef(i));
    }
    return -1 * total;
}

// softMax activation
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

// implementation of elu activation function
float activationFunction(float x) 
{ 
	return x < 0 ? 0.24 * (exp(x) - 1) : x;
} 
  
// elu derivative
float activationFunctionDerivative(float x) 
{ 
	return x < 0 ? activationFunction(x)+0.24 : 1;
}


// constructor of neural network class 
NeuralNetwork::NeuralNetwork(vector<size_t> topology, float learningRate) 
{
	this->topology = topology; 
	this->learningRate = learningRate; 
	// activation gradients for crossentropy backpropagation
	for (size_t i = 0; i < topology[topology.size()-1]; i++) {
		activation_gradients_.push_back(0.0f);
	}
	for (uint i = 0; i < topology.size(); i++) { 
        std::vector<float> vector;

		if (i == topology.size() - 1) 
			neuronLayers.push_back(make_unique<RowVector<float>>(topology[i])); 
		else
			neuronLayers.push_back(make_unique<RowVector<float>>(topology[i] + 1)); 

		// initialize cache and delta vectors 
		cacheLayers.push_back(make_unique<RowVector<float>>(neuronLayers.back()->length())); 
		deltas.push_back(make_unique<RowVector<float>>(neuronLayers.back()->length())); 

		if (i != topology.size() - 1) { 
			neuronLayers.back()->setValue(topology[i], 1.0); 
			cacheLayers.back()->setValue(topology[i], 1.0); 
		} 

		// initialze weights matrix 
		if (i > 0) { 
			if (i != topology.size() - 1) { 
				std::size_t size = (topology[i - 1] + 1) * (topology[i] + 1);
				std::vector<float> vector(size, 0.0f);
				
				weights.push_back(make_unique<Matrix<float>>(vector, topology[i - 1] + 1));
				momentum.push_back(make_unique<Matrix<float>>(vector, topology[i - 1] + 1)); 
				
				weights.back()->setNumber(sqrt(1/2840));
                auto fst = topology[i - 1];
                auto snd = topology[i];
				weights.back()->operator[]({fst, snd}) = 1.0f;
				
				for (unsigned x = 0; x < (topology[i - 1] + 1); x++) {
                    Matrix<float>& current_weigths = *weights.back();
                    current_weigths[{x, topology[i]}] = 0.0f;
				}
			} 
			else { 
				weights.push_back(make_unique<Matrix<float>>(topology[i - 1] + 1, topology[i])); 
				std::size_t size = (topology[i - 1] + 1) * (topology[i]);
				std::vector<float> vector(size, 0.0f);

				momentum.push_back(make_unique<Matrix<float>>(vector, topology[i - 1] + 1)); 
				weights.back()->setNumber(sqrt(1/1024));
			} 
		} 
	}
}

void NeuralNetwork::propagateForward(RowVector<float>& input) 
{ 
    // set first layer
	for (uint i = 0; i < input.length(); i++) {
		neuronLayers.front()->setValue(i, input.coeffRef(i));
		cacheLayers.front()->setValue(i, input.coeffRef(i));
	};
  
    // propagate the data forward 
    for (uint i = 1; i < topology.size(); i++) { 
        *neuronLayers[i] = *neuronLayers[i - 1] * *weights[i - 1];
		*cacheLayers[i] = *neuronLayers[i];

		// apply activation function
		if (i != topology.size() - 1) {
			for (uint ii = 0; ii < neuronLayers[i]->length(); ii++) {
				neuronLayers[i]->setValue(ii, activationFunction(neuronLayers[i]->coeffRef(ii)));
			}
		} else {
			softMax(*neuronLayers[i], *neuronLayers[i]);
		}
    };
} 

void NeuralNetwork::calcErrors(RowVector<float>& output) 
{ 
    // calculate the errors made by neurons of last layer 
	for (size_t i = 0; i < output.length(); i++) {
		deltas.back()->setValue(i, output.coeffRef(i) - neuronLayers.back()->coeffRef(i));
	}
	for (size_t i = 0; i < output.length(); i++) {
		float activation_grad = 0.0f;
		for (size_t j = 0; j < output.length(); j++) {
			if (i == j) {
                    activation_grad += neuronLayers.back()->coeffRef(i) * (1.0f - neuronLayers.back()->coeffRef(i))* deltas.back()->coeffRef(j);
                } else {
                    activation_grad += -neuronLayers.back()->coeffRef(i) * neuronLayers.back()->coeffRef(j) * deltas.back()->coeffRef(j);
                }
		}
		activation_gradients_[i] = activation_grad;
	}
  
    // error calculation of hidden layers 
    for (uint i = topology.size() - 2; i > 0; i--) { 
        *deltas[i] = *deltas[i + 1] * weights[i]->transpose(); 
        weights[i]->transpose(); //now transpose it back
    } 
} 

void NeuralNetwork::updateWeights() 
{ 
    for (uint i = 0; i < topology.size() - 1; i++) { 
        // compute weights in hidden layers
        if (i != topology.size() - 2) { 
            #pragma omp parallel for num_threads(8)
            for (unsigned c = 0; c < weights[i]->cols() - 1; c++) { 
                for (uint r = 0; r < weights[i]->rows(); r++) { 
					float num = weights[i]->operator[]({r, c}) + learningRate * 
                        deltas[i + 1]->coeffRef(c) * 
                        activationFunctionDerivative(cacheLayers[i + 1]->coeffRef(c)) * 
                        neuronLayers[i]->coeffRef(r);
                    weights[i]->operator[]({r, c}) = (num + momentum[i]->operator[]({r, c}))/2;
                    momentum[i]->operator[]({r, c}) = (num + momentum[i]->operator[]({r, c}))/2; 
                } 
            } 
        } 
		// compute weights in output layer
        else { 
            #pragma omp parallel for num_threads(8)
            for (uint c = 0; c < weights[i]->cols(); c++) { 
                for (uint r = 0; r < weights[i]->rows(); r++) { 
					float num = weights[i]->operator[]({r, c}) + learningRate * neuronLayers[topology.size()-2]->coeffRef(r) * activation_gradients_[c];
                    weights[i]->operator[]({r, c}) = (num + momentum[i]->operator[]({r, c}))/2;
                    momentum[i]->operator[]({r, c}) = (num + momentum[i]->operator[]({r, c}))/2; 
                } 
            } 
        } 
    } 
} 

void NeuralNetwork::propagateBackward(RowVector<float>& output) 
{ 
    calcErrors(output); 
    updateWeights(); 
} 

void NeuralNetwork::train(std::vector<RowVector<float>*> input_data, std::vector<RowVector<float>>& output_data, float learningRate) 
{ 
	this->learningRate = learningRate; 
    for (uint i = 0; i < input_data.size(); i++) { 
        if (i % 1000 == 0) {
            std::cout << "Checkpoint: " << i + 1 << " out of 59,000" << std::endl;
        }
        propagateForward(*input_data[i]); 
        propagateBackward(output_data[i]); 
    } 
}

void NeuralNetwork::predict(std::vector<RowVector<float>*> data, string outputFile)
{
	ofstream myfile;
	myfile.open(outputFile);
	
	for (uint i = 0; i < data.size(); i++) { 
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
