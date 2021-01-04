#include "../Dataset.hpp"
#include <iostream>



int main(void) {
    DataLoader loader("../data/fashion_mnist_train_vectors.csv", "../data/fashion_mnist_train_labels.csv");
    auto x = loader.load();
    std::cout << "Number of entries of the dataset: " << x.size() << std::endl;
}
