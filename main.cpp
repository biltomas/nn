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

int main() 
{ 
    DataLoader loader("../data/fashion_mnist_train_vectors.csv", "../data/fashion_mnist_train_labels.csv");
    auto x = loader.load();
    std::cout << "Number of entries of the dataset: " << x.size() << std::endl;

    std::vector<RowVector<float>> out_dat;
    std::vector<RowVector<float>*> in_dat;
    for (uint r = 0; r < x.size(); r++) { 
        RowVector<float> label = RowVector<float>(oneHotEncode(x[r].label, 10));
        
        out_dat.push_back(label);   
        for (int i = 0; i<x[r].data.length();i++) {
            x[r].data.setValue(i,x[r].data.coeffRef(i)/255);
        }
        in_dat.push_back(&x[r].data);

    }
    NeuralNetwork n({ 784, 100, 10 }); 
    std::cout << "Train data: " << in_dat.size() << std::endl;
    std::cout << "Train data[0] len: " << in_dat[0]->length() << std::endl;
    // std::vector<RowVector*> out_dat; 
    // std::vector<RowVector*> in_dat;
    // cout << out_dat.size() << endl;
    // cout << out_dat.back()->coeffRef(0) << endl;
    // cout << "rows " << out_dat.back()->vector.rows() << endl;
    // cout << "cols "<< out_dat.back()->vector.cols() << endl;
    // cout << "value "<< out_dat.back()->coeffRef(0) << endl;
    for (int i = 0; i < 1; i++)
        n.train(in_dat, out_dat); 
    
    DataLoader loader1("../data/fashion_mnist_test_vectors.csv", "../data/fashion_mnist_test_labels.csv");
    auto y = loader1.load();

    std::vector<RowVector<float>> test_out_dat;
    std::vector<RowVector<float>*> test_in_dat;
    for (uint r = 0; r < y.size(); r++) { 
        RowVector<float> label = RowVector<float>(oneHotEncode(y[r].label, 10));
        
        test_out_dat.push_back(label);   
        for (int i = 0; i<y[r].data.length();i++) {
            y[r].data.setValue(i,y[r].data.coeffRef(i)/255);
        }
        test_in_dat.push_back(&y[r].data);

    }
    n.predict(test_in_dat);

    return 0; 
} 
