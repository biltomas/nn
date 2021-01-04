#include "NeuralNetwork.hpp"
#include "Dataset.hpp"


typedef std::vector<RowVector<float>*> data;

// void genData() 
// { 
//     for (uint r = 0; r < 1000; r++) { 
//         float x = rand() / float(RAND_MAX); 
//         float y = rand() / float(RAND_MAX); 
//         float z = 2 * x + 10 + y; 
//         std::vector<float> out;
//         std::vector<float> in;
//         out.push_back(z);
//         in.push_back(x);
//         in.push_back(y);
//         RowVector out_row = RowVector(out);
//         RowVector in_row = RowVector(in);
//         out_dat.push_back(&out_row);
//         cout << out_dat.size() << endl;
        
//         in_dat.push_back(&in_row);
//     }  
// } 
std::vector<float> oneHotEncode(int out_dat, int num_classes) {
    vector<float> encoded;
    for (int i = 0; i < num_classes; i++) {
        if (out_dat == i) {
            encoded.push_back(1.0f);
        } else {
            encoded.push_back(0.0f);
        };
    };
    return encoded;
};

std::pair<std::vector<RowVector<float>*>, std::vector<RowVector<float>>>
encode_dataset(std::vector<item<float>>& dataset) {
    std::vector<RowVector<float>> out_dat;
    std::vector<RowVector<float>*> in_dat;
    for (uint r = 0; r < dataset.size(); r++) { 
        RowVector<float> label = RowVector<float>(oneHotEncode(dataset[r].label, 10));
        
        out_dat.push_back(label);   
        for (size_t i = 0; i < dataset[r].data.length(); i++) {
            dataset[r].data.setValue(i, dataset[r].data.coeffRef(i)/255);
        }
        in_dat.push_back(&dataset[r].data);

    }
    return {in_dat, out_dat};
}

int main() 
{ 
    std::srand(time(NULL));
    std::cout << "Loading the training dataset..." << std::endl;
    DataLoader loader("../data/fashion_mnist_train_vectors.csv", "../data/fashion_mnist_train_labels.csv");
    auto train_set = loader.load();
    std::cout << "Creating a validation set..." << std::endl;
    auto eval_set = loader.split_into_validation(train_set);
    std::cout << "One-hot encoding the labels of the train set..." << std::endl;
    auto encoded = encode_dataset(train_set);
    auto& in_dat = encoded.first;
    auto& out_dat = encoded.second;
    std::cout << "One-hot encoding the labels of the eval set..." << std::endl;
    auto encoded_eval = encode_dataset(eval_set);
    auto& eval_in_dat = encoded_eval.first;
    auto& eval_out_dat = encoded_eval.second;
    std::cout << "Number of entries of the training set: " << train_set.size() << std::endl;
    std::cout << "Number of entries of the validation set: " << eval_set.size() << std::endl;
    float lr = 0.0532;
    // float lr = 0.03;
    NeuralNetwork n({ 784, 512, 10 }, lr); 
    // std::vector<RowVector*> out_dat; 
    // std::vector<RowVector*> in_dat;
    // cout << out_dat.size() << endl;
    // cout << out_dat.back()->coeffRef(0) << endl;
    // cout << "rows " << out_dat.back()->vector.rows() << endl;
    // cout << "cols "<< out_dat.back()->vector.cols() << endl;
    // cout << "value "<< out_dat.back()->coeffRef(0) << endl;
    for (int i = 0; i < 8; i++) {
        std::cout << "Epoch " << i + 1 << " begins" << std::endl;
        n.train(in_dat, out_dat, lr); 
        lr *= 0.50;
        std::cout << "Validation starts" << std::endl;
        n.validate(eval_in_dat, eval_out_dat);
    }
    
    DataLoader loader1("../data/fashion_mnist_test_vectors.csv", "../data/fashion_mnist_test_labels.csv");
    auto y = loader1.load();

    std::vector<RowVector<float>> test_out_dat;
    std::vector<RowVector<float>*> test_in_dat;
    for (uint r = 0; r < y.size(); r++) { 
        RowVector<float> label = RowVector<float>(oneHotEncode(y[r].label, 10));
        
        test_out_dat.push_back(label);   
        for (size_t i = 0; i<y[r].data.length();i++) {
            y[r].data.setValue(i,y[r].data.coeffRef(i)/255);
        }
        test_in_dat.push_back(&y[r].data);

    }
    n.predict(test_in_dat);

    return 0; 
} 
