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

int main() 
{ 
    DataLoader loader("../data/fashion_mnist_train_vectors.csv", "../data/fashion_mnist_train_labels.csv");
    auto x = loader.load();
    std::cout << "Number of entries of the dataset: " << x.size() << std::endl;

    std::vector<RowVector<float>*> out_dat;
    std::vector<RowVector<float>*> in_dat;
    for (uint r = 0; r < x.size(); r++) { 
        out_dat.push_back(x[r].label);
            
        in_dat.push_back(&x[r].data);

    }
    NeuralNetwork n({ 2, 3, 1 }); 
    // std::vector<RowVector*> out_dat; 
    // std::vector<RowVector*> in_dat;
    // cout << out_dat.size() << endl;
    // cout << out_dat.back()->coeffRef(0) << endl;
    // cout << "rows " << out_dat.back()->vector.rows() << endl;
    // cout << "cols "<< out_dat.back()->vector.cols() << endl;
    // cout << "value "<< out_dat.back()->coeffRef(0) << endl;
    // for (int i = 0; i < 1; i++)
    //     n.train(in_dat, out_dat); 
    return 0; 
} 
