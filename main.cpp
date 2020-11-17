#include "NeuralNetwork.hpp"

typedef std::vector<RowVector*> data;

int main() 
{ 
    NeuralNetwork n({ 2, 3, 1 }); 
    RowVector* out1 = new RowVector({3});
    RowVector* out2 = new RowVector({3});
    RowVector* in1 = new RowVector({1, 2});
    RowVector* in2 = new RowVector({2, 1});
    std::vector<RowVector*> out_dat = {out1, out2}; 
    std::vector<RowVector*> in_dat = {in1, in2};
    n.train(in_dat, out_dat); 
    return 0; 
} 